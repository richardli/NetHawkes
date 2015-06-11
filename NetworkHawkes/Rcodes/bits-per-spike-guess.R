## script to simulate range of predictive log likelihood
##


# function to calculate next set of rates 
# TODO: let mu and tau be part of input
#       different mu and tau for each k-k' pair
# return: K+1 vector of baseline rate plus impulses
get_next_rate <- function(allrates, A, W, events, mu, tau, dt){
	# initialize further local variables
	tmax = 10
	binmax = 10 / dt
	g <- function(t, receiver){(tmax / (t * (tmax - t))) * sqrt(tau[,receiver] / 2 / pi) * exp(-tau[,receiver] / 2 *(log(t / (tmax - t)) - mu[,receiver])^2)}
	# plot(1:10, g(1:10), type = "l")
	N <- dim(allrates)[1]
	current <- dim(allrates)[2]
	# 1 2 3 4 5 6 7 8 9 10 11 12 13 (simulate 14)
	# need 6 7 8 9 10 11 12 13, 8 previous events (14 - 10 + 2)
 	start <- max(1, current - binmax + 2)
	
	rate <- matrix(0, N, N+1)
	for(t in start : current){
		# add to rate
		for(i in 1:N){
			# note: add dt
			impulse <- g((current + 1 - t ) * dt, i)
			# rate[i] <- rate[i] + impulse * sum(A[, i] * W[, i] * events[, t])
			rate[i, ] <- c(rate[i, 1], 
						   rate[i, -1] + impulse * A[, i] * W[, i] * events[, t])
		}
	}
	return(rate)
}
# function to simulate events from rate
simulate <- function(rates, dt){
	# note: add dt
	return(rpois(length(rates), rates * dt))
}
Nrep <- 100
## simulate as the paper suggests
T = 1000
T.train = 900
KK = 10
realtest1 <- rep(0, Nrep)
realtest2 <- rep(0, Nrep)

for(rep in 1:Nrep){
	baseline <- rgamma(KK, 1, 1)
	baseline <- matrix(rep(baseline, T), KK, T)
	shape <- 2
	rate <- 10
	W <- matrix(rgamma(KK^2, shape = shape, rate = rate), KK, KK)
	diag(W) <- 0
	p <- 0.3
	A <- matrix(rbinom(KK^2, 1, p), KK, KK)
	mu <- tau <- matrix(0, KK, KK)
	for(k in 1:KK){
		for(kp in 1:KK){
			tau[k, kp] <- rgamma(1, 10, 1)
			mu[k, kp] <- rnorm(1, -1, 1/(10 * tau[k, kp]))				
		}
	}

	allrates <- array(0, dim = c(KK, T, KK+1))
	rates <- baseline
	allrates[, , 1] <- baseline
	events <- matrix(0, KK, T)
	# simulate first time point
	events[, 1] <- simulate(rates[, 1], dt)
	for(t in 2:T){
		impulse_to_add <- get_next_rate(rates[, 1:(t-1),drop=FALSE], 
						  				A, W, 
						  				events[, 1:(t-1), drop=FALSE], 
						  				mu, tau, dt)
		allrates[, t, ] <- allrates[, t, ] + impulse_to_add				   
		rates[, t] <- apply(allrates[, t, ], 1, sum)
		events[, t] <- simulate(rates[, t], dt)
		if(length(which(is.na(events[, t]))) > 0) stop(t)
	}
	rate.homo <- matrix(apply(events[, 1:T.train], 1, mean), KK,  T-T.train)
	rate.true <- rates[, (T.train+1):T]

	ll0 <- sum(log(dpois(events[, ((T.train+1):T)], rate.homo)))
	ll1 <- sum(log(dpois(events[, ((T.train+1):T)], rate.true)))
	ll1_best <- sum(log(dpois(events[, ((T.train+1):T)], events[, ((T.train+1):T)])))

	realtest1[rep] <- (ll1 - ll0) / sum(events[, ((T.train+1):T)])
	realtest2[rep] <- (ll1_best - ll0) /sum(events[, ((T.train+1):T)])

}

##
##
##
##
##
##

# simulate rates for 1000 time points
T <- 1000
train <- 1:900
test <- 901:1000
Nrep <- 100


# assume rate is 0.1 in the training, varying true testing rate
test.rate <- seq(1, 50)
test1 <- rep(0, length(test.rate))
test2 <- rep(0, length(test.rate))
test3 <- rep(0, length(test.rate))
test4 <- rep(0, length(test.rate))

for(case in 1:length(test.rate)){
	bps <- rep(0, Nrep)
	bps2 <- rep(0, Nrep)
	bps3 <- rep(0, Nrep)
	bps4 <- rep(0, Nrep)
	for(i in 1:Nrep){
		# rate <- runif(T, 0.1, 10) + (1:T)/(T/5)
		rate <- c(rep(1, length(train)), rep(test.rate[case], length(test)))
		rate.rev <- c(rep(50, length(train)), rep(test.rate[case], length(test)))

		count <- rpois(T, rate)
		count.rev <- rpois(T, rate.rev)
		rate0 <- mean(count[train])
		rate0.rev <- mean(count.rev[train])

		plot(rate, type = "l", ylim = range(c(rate, count)))
		points(count, col = "red", cex = .5)
		abline(h = rate0, col = "blue")	

		ll0 <- sum(log(dpois(count[test], rate0)))
		ll1 <- sum(log(dpois(count[test], rate[test])))
		ll1_best <- sum(log(dpois(count[test], count[test])))


		ll0.rev <- sum(log(dpois(count.rev[test], rate0.rev)))
		ll1.rev <- sum(log(dpois(count.rev[test], rate.rev[test])))
		ll1_best.rev <- sum(log(dpois(count.rev[test], count.rev[test])))

		bps[i] <- (ll1 - ll0) / sum(count[test])
		bps2[i] <- (ll1_best - ll0) /sum(count[test])

		bps3[i] <- (ll1.rev - ll0.rev) / length(test)
		bps4[i] <- (ll1_best.rev - ll0.rev) / length(test)

		print(c(bps[i], bps2[i], bps3[i], bps4[i]))

	}

	# true rate, as described in the paper 
	test1[case] <- mean(bps)
	# best case, as described in the paper
	test2[case] <- mean(bps2)

	# true rate, in the codes 
	test3[case] <- mean(bps3)
	# best case, in the codes
	test4[case] <- mean(bps4)
}

save.image("bits-per-spike-guess.RData")

pdf("bits-per-spike-guess.pdf", height = 5, width = 14)
par(mfrow = c(1, 3))
library(ggplot2) 
realtest <- c(realtest1, realtest2)
case <- c(rep("True rate", length(realtest1)), rep("best case", length(realtest2)))
realtest <- data.frame(bits_per_spike = realtest, case = case)
# qplot(bits_per_spike, data = realtest, binwidth = 0.01, 
# 	fill = case) 
boxplot(bits_per_spike ~ case, data = realtest, main = "simulated as in paper")

plot(test.rate, test1, type = "l", 
	ylab = "Bit/spike improvement from baseline", 
	xlab = "true rate in test set", 
	main = "Bit/spike improvement from baseline \n fixed training rate = 1")
lines(test.rate, test2, type = "l", col = "red")
legend("topleft", c("knowing true rate", "best case"), 
	lty = c(1,1), col = c("black", "red"))

plot(test.rate, test3, type = "l", xlim = c(20, 50), ylim = c(0, 10),
	ylab = "Bit/spike improvement from baseline", 
	xlab = "true rate in test set", 
	main = "Bit/spike improvement from baseline \n fixed training rate = 50")
lines(test.rate, test4, type = "l", col = "red")
abline(h = 2, col = "red", lty = "longdash")
legend("topright", c("knowing true rate", "best case"), 
	lty = c(1,1), col = c("black", "red"))
dev.off()
# plot(test.rate, test3, type = "l", 
# 	ylab = "bits/spike (normalized by number of bins)", 
# 	xlab = "true rate in test set", 
# 	main = "Bit/spike improvement from baseline \n modified calculation (training rate = 1)")
# lines(test.rate, test4, type = "l", col = "red")
# legend("topleft", c("knowing true rate", "best case"), 
# 	lty = c(1,1), col = c("black", "red"))


