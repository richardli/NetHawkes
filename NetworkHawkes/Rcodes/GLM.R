##
## Code implementing GLM for neural spike data
##
## note: results seem not very stable, sometimes AUC very low. 
## need to rerun multiple times to find a good AUC fit. 
## BPS is not too bad though. 
##
remove(list = ls())
library(ROCR)
library(glmnet)
source("spikeGLM.r")


###########################################################################
##  Read Synthetic data and run model
###########################################################################
perf.all.glm <- list()

for(net in 1:10){
	runname = paste("rhawkes-sim/seed9-net", net, sep = "")
	load(paste("../data/", runname, "-wk.rda", sep = ""))
	K <- 30
	T <- 1000
	T.train <- 900

	perf.all.glm[[net]] <- spikeGLM(K = K, T = T, T.train = T.train,  
									T.train.sub = 100, events = this$events, 
									noself = TRUE,
									N.rand = 10000, B = 6, tmax = 10, 
									A = this$A, W = this$W, rates = this$rates, 
									verbose = TRUE)
	auc <- perf.all.glm[[net]]$auc
	if(auc < 0.9){
		perf.all.glm[[net]] <- spikeGLM(K = K, T = T, T.train = T.train,  
									T.train.sub = 100, events = this$events, 
									noself = TRUE,
									N.rand = 10000, B = 6, tmax = 10, 
									A = this$A, W = this$W, rates = this$rates, 
									verbose = TRUE)
	}
}
save(perf.all.glm, file = paste("../data/", runname, "-model-perf-GLM.rda", sep = ""))

###########################################################################
##  Read Synthetic data and run model (without self-connect)
###########################################################################
perf.all.glm <- list()

for(net in 1:10){
	runname = paste("rhawkes-sim/seed5-net", net, sep = "")
	load(paste("../data/", runname, "-wk.rda", sep = ""))
	K <- 30
	T <- 1000
	T.train <- 900

	perf.all.glm[[net]] <- spikeGLM(K = K, T = T, T.train = T.train,  
									T.train.sub = 100, events = this$events, 
									noself = FALSE,
									N.rand = 10000, B = 6, tmax = 10, 
									A = this$A, W = this$W, rates = this$rates, 
									verbose = TRUE)
	## usually the GLM fails randomly, re-run again automatically if auc too low
	auc <- perf.all.glm[[net]]$auc
	if(auc < 0.9){
		perf.all.glm[[net]] <- spikeGLM(K = K, T = T, T.train = T.train,  
									T.train.sub = 100, events = this$events, 
									noself = FALSE,
									N.rand = 10000, B = 6, tmax = 10, 
									A = this$A, W = this$W, rates = this$rates, 
									verbose = TRUE)
	}
}

save(perf.all.glm, file = paste("../data/", runname, "-model-perf-GLM.rda", sep = ""))



###########################################################################
##  Read Stock data and run model
###########################################################################

# read data
# Stock price data
load("../data/sp100/sp-events.rda")
K <- 100
T <- 23290
T.train <- T * 0.8

fit <- spikeGLM(K = K, T = T, T.train = T.train,  
				T.train.sub = 0, events = events, 
				N.rand = 1000, B = 6, tmax = 12, 
				verbose = TRUE)







