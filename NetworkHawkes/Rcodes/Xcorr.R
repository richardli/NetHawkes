## 
## function to perform Cross-correlation prediction
## takes event matrix and return weight adj
##
library(ROCR)

Xcorr <- function(event, maxlag = 10, binsize = NULL, A = NULL, verbose = TRUE){
	K <- dim(event)[1]
	W <- matrix(0, K, K)
	if(!is.null(binsize)){
		M <- trunc((dim(event)[2]  - 1) / binsize) + 1
		event.dense <- matrix(0, K, M)
		for(i in 1:dim(event)[2]){
			j <- trunc((i-1) / binsize) + 1
			event.dense[, j] <- event.dense[, j] + event[, i]
		}
	}else{
		event.dense <- event
	}
	for(i in 1:K){
		for(j in 1:K){
			CC <- ccf(event.dense[i, ], event.dense[j, ], plot = FALSE,
				lag.max = maxlag, type = "correlation")
			for(l in 0:maxlag){
				W[i,j] = W[i,j] + CC$acf[which(CC$lag == l)]
			}
		}
	}
	if(!is.null(A)){
		pred <- prediction(as.vector(W), as.vector(A))
		auc <- performance(pred, "auc")@y.values[[1]]
		perf <- performance(pred,"tpr","fpr")
		if(verbose){print(auc)}
	}else{
		pred  <- NULL
		auc <- NULL
		perf <- NULL
	}

	return(list(W = W, pred = pred, auc = auc, perf = perf))
}


## perform on synthetic data
perf.all.xcorr <- list()

for(net in 1:10){
	runname = paste("rhawkes-sim/seed9-net", net, sep = "")
	load(paste("../data/", runname, "-wk.rda", sep = ""))
	perf.all.xcorr[[net]] <- Xcorr(this$events, maxlag = 2, binsize = 5, A = this$A)
}
save(perf.all.xcorr, 
	file = paste("../data/", runname, "-model-perf-Xcorr.rda", sep = ""))



