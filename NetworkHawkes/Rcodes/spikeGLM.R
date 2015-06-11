##
## function to perform GLM fit
##
## K: number of processes
## T: Total time points
## T.train: number of time points used for training
## T.train.sub: number of time points used for CV in training set (if not set to 0)
## events: K by T matrix of event counts
## N.rand: number of random basis to sample
## B: number of basis to use
## tmax: max time points for impulse memory
##

spikeGLM <- function(K = 30, T = 1000, T.train = 900,  T.train.sub = 100, 
	events = NULL, noself = TRUE,
	N.rand = 5000, B = 6, tmax = 10, 
	A = NULL, W = NULL, rates = NULL, verbose = TRUE
	){
	## define impulse function 
	g <- function(t, mu, tau){
		tmaxv <- rep(tmax, length(t))
		return(tmaxv / (t * (tmaxv - t)) * sqrt(tau / 2 / pi) * exp(-tau / 2 *(log(t / (tmaxv - t)) - mu )^2))
	}

	## get random impulse basis 
	base.rand <- matrix(0, N.rand, tmax - 1)
	timelag <- (tmax-1) : 1
	for(i in 1:N.rand){
		tau <- rgamma(1, 10, 1)	
		mu <- rnorm(1, -1, 1/10/tau)
		base.rand[i, ] <- g(timelag, mu, tau)
	}

	## get principle components for the random bases
	pca.fit <- princomp(base.rand)
	pcabasis <- pca.fit$loadings[, 1:B]

	## get basis for each event
	# for easier matrix manipulation, add tmax - 1 columns to events
	eventsnow <- cbind(matrix(0, K, tmax - 1), events)
	# covariate T by K by K*B matrix, 
	# (t, k, kp, b) -> at time t, kp has this influence on k via basis b 
	basis <- array(NA, dim = c(T, K, K*B))
	# for(t in 1: (T + tmax - 1)){
	for(t in 1:T){
		for(k in 1:K){
			for(i in 1: (K*B) ){
				kp <- i %% K 
				if(kp == 0) kp <- K
				b <- (i - kp) / K + 1
				hist <- (t) : (t + tmax - 1 - 1)

				if(noself && kp == k){
					basis[t, k, i] <- 0
				}else{
					basis[t, k, i] <- sum( pcabasis[, b] * eventsnow[kp, hist])
				}
				
			}
		}
	}

	# for L1 regularized version
	alpha <- rep(0, K)
	beta <- array(0, dim = c(K, K, B))
	# for non regularized version
	alpha2 <- rep(0, K)
	beta2 <- array(0, dim = c(K, K, B))

	rate.fit <- matrix(0, K, T - T.train)
	# use the last few iterations for CV
	T.train.last <- T.train - T.train.sub 

	for(k in 1:K){
		dataset <- as.matrix(data.frame(outcome = events[k, ], 
							  basis = basis[, k, ]))
		lambda <- 0

		##  if running CV first
		if(T.train.sub > 0){
			fit0 <- cv.glmnet(y = dataset[T.train.last:T.train, 1] , 
					  x = dataset[T.train.last:T.train, -1], 
					  family = "poisson", type.measure="mse")
			lambda <- fit0$lambda.min
		}
		fit <- glmnet(y = dataset[1:T.train, 1], 
					  x = dataset[1:T.train, -1], 
					  family = "poisson")

		# dataset <- data.frame(dataset)
		# fit <- glm(outcome ~., family = "poisson", data = dataset[1:T.train, ])
		# find the smallest lambda
		lambda <- fit$lambda[length(fit$lambda)]
		# use the least regularized lambda to avoid too many zeros
		pred <- predict(fit, newx =  dataset[(T.train + 1):T, -1], 
			type = "response", s = lambda)
		print(summary(pred))
		rate.fit[k, ] <- pred
		tmp <- coef(fit)
		alpha[k] <- tmp[1]
		beta[k, , ] <- matrix(tmp[-1], K, B)
		## only do non regularized version when calculating AUC
		if(!is.null(A)){
			fit_nonreg <- glm(dataset[1:T.train, 1] ~ dataset[1:T.train, -1], family = "poisson")
			tmp2 <- as.numeric(coef(fit_nonreg))
			alpha2[k] <- tmp2[1]
			beta2[k, , ] <- matrix(tmp2[-1], K, B)			
		}
		cat(".")
	}
	# beta[which(is.na(beta))] <- 0

	# calculate likelihood
	rate.homo <- matrix(apply(events[, 1:T.train], 1, mean),
						K, (T-T.train))
	# avoid zero rate
	rate.fit[which(rate.fit == 0)] <- min(rate.fit[rate.fit > 0]) / 2

	ll.homo <- sum(log(dpois(events[, ((T.train+1):T)], rate.homo)))
	ll.fit <- sum(log(dpois(events[, ((T.train+1):T)], rate.fit)))
	
	if(!is.null(rates)){
		rate.true <- rates[, ((T.train+1):T)]
		ll.true <- sum(log(dpois(events[, ((T.train+1):T)], rate.true)))
	}else{
		rate.true <- NULL
		ll.true <- NULL
	}

	bps0 = (ll.fit - ll.homo) / sum(events[, ((T.train+1):T)])
	bps1 = (ll.fit - ll.homo) / (T - T.train)
	
	if(verbose){
		print("Bits per spike")
		print(bps0)

		print("Bits per second")
		print(bps1)		
	}

	# if ground truth is known
	if(!is.null(A)){
		## see link prediction with L1 regularization
		## since weight is receiver by sender, A is sender by receiver
		beta[which(is.na(beta))] <- 0
		weight <- t(apply(beta, c(1, 2), mean))
		pred.L1 <- prediction(as.vector(weight), as.vector(A))

		## see link prediction with no regularization
		## since weight is receiver by sender, A is sender by receiver
		beta2[which(is.na(beta2))] <- 0
		weight2 <- t(apply(beta2, c(1, 2), mean))
		pred <- prediction(as.vector(weight2), as.vector(A))

		auc.L1 <- performance(pred.L1, "auc")@y.values[[1]]
		auc <- performance(pred, "auc")@y.values[[1]]

		perf.L1 <- performance(pred.L1,"tpr","fpr")
		perf <- performance(pred,"tpr","fpr")
		
		if(verbose){
			print("auc")
			print(auc)
			print("auc L1")
			print(auc.L1)
		}
	}else{
		pred.L1 <- NULL
		pred <- NULL
		auc.L1 <- NULL
		auc <- NULL
		perf.L1 <- NULL
		perf <- NULL
	}

	fit <- list(rate.fit = rate.fit, 
				rate.homo = rate.homo[1],
				rate.true = rate.true,

				ll.homo = ll.homo, 
				ll.true = ll.true, 
				ll.fit = ll.fit, 

				bps0 = bps0, 
				bps1 = bps1, 

				pred.L1 = pred.L1, 
				pred = pred, 
				auc.L1 = auc.L1, 
				auc = auc, 
				perf.L1 = perf.L1, 
				perf = perf)

	return(fit)
}


	
	
	
	

	