T = 1000
T.train = 900
N_net = 10
KK = 30

dt = 1

N_mean = 25732
N_mar = 9425

# function to calculate next set of rates 
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


# store the outcomes
data <- list()
# seed 5 does not allow self-excitement
seed <- 5
# seed 7 allows self-excitement
seed <- 7
# seed 9 for some wild initialization
seed <- 9
set.seed(seed)
total_event_profile <- rep(0, N_net)
seed.save <- NULL

for(i in 1:N_net){
	try <- TRUE
	while(try){
			seed.sub <- round(runif(1, 1e4, 1e6))
			not_explode <- FALSE
			n.explode <- 0
			while(!not_explode){
				set.seed(seed.sub)
				unstable <- TRUE
				isEmpty <- TRUE
				print("+")
			   
				while(unstable || isEmpty){

					cat(".")
					if(seed == 9){
						baseline <- runif(KK, 0.1, 1)
						baseline <- matrix(rep(baseline, T), KK, T)
						W <- matrix(rgamma(KK^2, 2, 5), KK, KK)
						# diag(W) <- 0
						p <- 0.05
						A <- matrix(rbinom(KK^2, 1, p), KK, KK)
						isEmpty <- (sum(A * W) < 0.1)
					}
					if(seed == 7){
						baseline <- runif(KK, 0.1, 0.5)
						baseline <- matrix(rep(baseline, T), KK, T)
						W <- matrix(rgamma(KK^2, 2, 5), KK, KK)
						# diag(W) <- 0
						p <- 0.01
						cutoff <- quantile(as.vector(W[W>0]), 1-p)
						A <- matrix(0, KK, KK)
						for(iii in 1:KK){
							for(jjj in 1:KK){
								if(W[iii, jjj] > cutoff) A[iii, jjj] <- 1
							}
						}
						isEmpty <- (sum(A * W) < 0.1)
					}
					if(seed == 5){
						baseline <- runif(KK, 0.1, 0.5)
						baseline <- matrix(rep(baseline, T), KK, T)
						W <- matrix(rgamma(KK^2, 2, 5), KK, KK)
						diag(W) <- 0
						p <- 0.01
						cutoff <- quantile(as.vector(W[W>0]), 1-p)
						A <- matrix(0, KK, KK)
						for(iii in 1:KK){
							for(jjj in 1:KK){
								if(W[iii, jjj] > cutoff) A[iii, jjj] <- 1
							}
						}
						isEmpty <- (sum(A * W) < 0.1)
					}
					# check stability constraint (maybe doesn't matter)
					if(max(abs(eigen(A * W)$values)) < 10) unstable <- FALSE
				}
				
				allrates <- array(0, dim = c(KK, T, KK+1))
				rates <- array(0, dim= c(KK, T))
				rates <- baseline
				allrates[, , 1] <- baseline
				events <- matrix(0, KK, T)
				# simulate first time point
				events[, 1] <- simulate(rates[, 1], dt)

				# simulate impulse function parameters
				mu <- tau <- matrix(0, KK, KK)
				for(k in 1:KK){
					for(kp in 1:KK){
						tau[k, kp] <- rgamma(1, 10, 1)
						mu[k, kp] <- rnorm(1, -1, 1/(10 * tau[k, kp]))				
					}
				}
				for(t in 2:T){
					impulse_to_add <- get_next_rate(rates[, 1:(t-1),drop=FALSE], 
									  				A, W, 
									  				events[, 1:(t-1), drop=FALSE], 
									  				mu, tau, dt)
					# print(range(impulse_to_add))
					
					allrates[, t, ] <- allrates[, t, ] + impulse_to_add				   
					
					rates[, t] <- apply(allrates[, t, ], 1, sum)
					if(max(rates[, t]) > 10000) {
						break
					}

					events[, t] <- simulate(rates[, t], dt)
				}
				## if break out early
				if(t != T){
					print("too high rate skip")
					seed.sub <- round(runif(1, 1e4, 1e6))
					next
				}

				not_explode <- (sum(events) < KK * 1e4)
				if(sum(events) > KK * 1e4) print(sum(events))
				if(is.na(not_explode)) not_explode <- FALSE
				if(!not_explode){
					n.explode <- n.explode + 1
					seed.sub <- round(runif(1, 1e4, 1e6))
				}
				if(n.explode > 20) stop("too many failure times")
			}	
 
			cat("\n")
			print(paste("baseline rate  : ", apply(baseline, 1, mean)))
			print(paste("max eigenvalue : ", max(abs(eigen(A * W)$values))))
			print(paste("total events   : ", sum(events)))
			print(paste("range in paper : ", 25732 - 9425, "to", 25732 + 9425))
			cat("\n")
			#######################################################################
			#######################################################################

			rates.dense <- matrix(0, dim(rates)[1], dim(rates)[2] * dt)
			events.dense <- matrix(0, dim(rates)[1], dim(rates)[2] * dt)
			for(ii in 1:T){
				i.dense <- trunc((ii-1)*dt)+1
				rates.dense[, i.dense] <- rates.dense[, i.dense] + rates[,ii] * dt
				events.dense[, i.dense] <- events.dense[, i.dense] + events[,ii]
			}

			rate.homo <- matrix(apply(events.dense[, 1:(T.train*dt)], 1, mean), KK, ((T-T.train)*dt))
			rate.true <- rates.dense[, ((T.train+1):T)*dt]
			
			ll0 <- sum(log(dpois(events.dense[, ((T.train+1):T)*dt], rate.homo)))
			ll1 <- sum(log(dpois(events.dense[, ((T.train+1):T)*dt], rate.true)))

			
			print("Bits per spike")
			bps0 = (ll1 - ll0) / sum(events.dense[, ((T.train+1):T)*dt])
			print(bps0)

			print("Bits per second")
			bps3 = (ll1 - ll0) / (T - T.train)
			print(bps3)


			# save only networks with potential ok improvements
			# and to match Linderman's simulated data roughly
			if((bps3 > 1) &&
				(sum(events) < 45000) &&
				(sum(events) > 15000)){
				try <- FALSE
				seed.save <- c(seed.save, seed.sub)
				seed.sub <- round(runif(1, 1e4, 1e6))
			}else{
				seed.sub <- round(runif(1, 1e4, 1e6))
			}

	}
	

	#########################################################################
	# write network information to file
	#########################################################################
	print(paste("save network  ", i))
	total_event_profile[i] <- sum(events)

	cvec <- NULL
	svec <- NULL
	for(t in 1:T){
		cvec <- c(cvec, rep(1:KK, events[, t]))
		svec <- c(svec, rep(t, sum(events[, t])))
	}	
	# note: add dt
	svec <- svec * dt

	#########################################################################
	# save all information
	#########################################################################	
	data[[i]] <- list(baseline = baseline, 
					  p = p, 
					  A = A, 
					  W = W, 
					  rates = rates, 
					  allrates = allrates, 
					  events = events,
					  svec = svec, 
					  cvec = cvec,
					  seed = seed.save, 
					  bps = c(bps0, bps3) )
}




#########################################################################
# write to file
#########################################################################
for(i in 1:N_net){	
	filename1 <- paste("../data/rhawkes-sim/seed", seed, "-net", i,"-time.txt", sep = "")
	filename2 <- paste("../data/rhawkes-sim/seed", seed, "-net", i,"-proc.txt", sep = "")
	file1 <- file(filename1)
	writeLines(as.character(data[[i]]$svec*dt^2), file1, sep = ",")
	# writeLines(as.character(data[[i]]$svec), file1, sep = ",")
	close(file1)
	file2 <- file(filename2)
	writeLines(as.character(data[[i]]$cvec), file2, sep = ",")
	close(file2)
}

print(summary(total_event_profile))
# plot(data[[2]]$svec, data[[2]]$cvec, xlab = "time", ylab = "process", cex = .1)
for(i in 1:N_net){
this <- data[[i]]
save(this, 
	file = paste("../data/rhawkes-sim/seed", seed, "-net", i,"-wk.rda", sep = "")
	)
}
