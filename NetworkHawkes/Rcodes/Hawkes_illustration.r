# function to calculate next set of rates 
# TODO: let mu and tau be part of input
#       different mu and tau for each k-k' pair
# return: K+1 vector of baseline rate plus impulses
get_next_rate <- function(allrates, A, W, events){
	# initialize further local variables
	mu = -1
	tau = 1
	tmax = 10
	binmax = 10
	g <- function(t){(tmax / (t * (tmax - t))) * sqrt(tau / 2 / pi) * exp(-tau / 2 *(log(t / (tmax - t)) - mu)^2)}
	# plot(1:10, g(1:10), type = "l")
	N <- dim(allrates)[1]
	current <- dim(allrates)[2]
	# 1 2 3 4 5 6 7 8 9 10 11 12 13 (simulate 14)
	# need 6 7 8 9 10 11 12 13, 8 previous events (14 - 10 + 2)
 	start <- max(1, current - binmax + 2)
	
	rate <- matrix(0, N, N+1)
	for(t in start : current){
		# all impulses
		impulse <- g(current + 1 - t)
		# add to rate
		for(i in 1:N){
			# rate[i] <- rate[i] + impulse * sum(A[, i] * W[, i] * events[, t])
			rate[i, ] <- c(rate[i, 1], 
						   rate[i, -1] + impulse * A[, i] * W[, i] * events[, t])
			# cat(impulse)
			# cat(" ")
			# cat(impulse * sum(A[, i] * W[, i] * events[, t]))
			# cat("\n")
		}
	}
	return(rate)
}
# function to simulate events from rate
simulate <- function(rates){
	# return(round(rates))
	return(rpois(length(rates), rates))
}

N_node = 30
T = 1000
T.train = 900
unstable <- TRUE
count <- 0
# baseline <- runif(1, .3, 1)
baseline <- rgamma(N_net, 1, 1)
# baseline <- 1
while(unstable && count < 1000){
	count <- count + 1
	######################################################
	# follow Linderman's initialization codes
	#    kappa = 2.0
	#    c = np.arange(C).repeat((K // C))
    #    p = 0.4 * np.eye(C) + 0.01 * (1-np.eye(C))
    #    v = kappa * (5 * np.eye(C) + 5.0 * (1-np.eye(C)))
    ######################################################
    shape <- 2
    rate <- shape * 5 
	W <- matrix(rgamma(N_node^2, shape = shape, rate = rate), N_node, N_node)
	# W <- matrix(0.05, N_node, N_node)
	# diag(W) <- 0
	#######################################################
	# Initial idea of generating A
	# p <- runif(1, .01, .02)
	p <- .1 # rbeta(1, 1, 1) 
	A <- matrix(rbinom(N_node^2, 1, p), N_node, N_node)
	# A <- matrix(1, N_node, N_node)
	######################################################
	# initialize A to top 50%? 5%? of W
	# sparsity <- 0.05
	# A <- W
	# large <- which(W >= quantile(W, (1 - sparsity)))
	# small <- which(W < quantile(W, (1 - sparsity)))
	# A[large] <- 1
	# A[small] <- 0
	print(max(abs(eigen(A * W)$values)))

	if(max(abs(eigen(A * W)$values)) < 1) unstable <- FALSE
}
if(count == 1000 && unstable){
	stop("Unstable after 1000 itr!")
}
allrates <- array(0, dim = c(N_node, T, N_node+1))
rates <- array(0, dim= c(N_node, T))
rates[, 1] <- baseline
allrates[, , 1] <- baseline
events <- matrix(0, N_node, T)
# simulate first time point
events[, 1] <- simulate(rates[, 1])

for(t in 2:T){
	impulse_to_add <- get_next_rate(rates[, 1:(t-1),drop=FALSE], 
					  				A, W, events[, 1:(t-1), drop=FALSE])
	# print(range(impulse_to_add))
	
	allrates[, t, ] <- allrates[, t, ] + impulse_to_add				   
	
	rates[, t] <- apply(allrates[, t, ], 1, sum)

	events[, t] <- simulate(rates[, t])
}

if(!unstable){
	print(sum(events))

	rate.homo <- matrix(apply(events[, 1:T.train], 1, mean), N_node, (T-T.train))
	rate.true <- rates[, (T.train+1):T]
	
	ll0sep <- apply(log(dpois(events[, (T.train+1):T], rate.homo)), 1, sum)
	ll1sep <- apply(log(dpois(events[, (T.train+1):T], rate.true)), 1, sum)
	countsep <- apply(events[, (T.train+1):T], 1, sum)

	print("seperate")
	print(
		sum( (ll1sep - ll0sep) / countsep)
	)

	ll0 <- sum(log(dpois(events[, (T.train+1):T], rate.homo)))
	ll1 <- sum(log(dpois(events[, (T.train+1):T], rate.true)))
	ll1_best <- sum(log(dpois(events[, (T.train+1):T], events[, (T.train+1):T]))) 

	print("together")
	print(
	(ll1 - ll0) / sum(events[, (T.train+1):T])
	)
	print(
	(ll1 - ll0) / sum(events[, (T.train+1):T]) * (T - T.train)* N_node
	)
	
	print("best case")
	print(
	(ll1_best - ll0) / sum(events[, (T.train+1):T])
	)
	
	print("In the code")
	print(
	(ll1 - ll0) / prod(dim(events[, (T.train+1):T]))
	)

	print("best case (In the code)")
	print(
	(ll1_best - ll0) /  prod(dim(events[, (T.train+1):T]))
	)

	print("In the code (revised)")
	print(
	(ll1 - ll0) / length(which(events[, (T.train+1):T] != 0))
	)
	
	

	########################
	i = 5
	par(mfrow = c(2, 1))

	plot(rate.true[i, ], type = "l", 
		ylim = range(c(rate.true[i, ], events[i, (T.train+1):T])))
	lines(rate.homo[i, ], col = "blue")
	points(events[i, (T.train+1):T], col = "red")

	plot(rates[i, ], type = "l",
	ylim = range(c(rates[i, ], events[i, ])))
	abline(h = rate.homo[i, 1], col = "blue")
	points(events[i, ], col = "red", cex = .2)
	#########################

}
