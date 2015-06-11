########################################################################
##
## script to load results from java and combine with simulated R data
##
########################################################################
library(ROCR)
# read command line input
args <- commandArgs(TRUE)
case <- as.numeric(args[1])
runname_java <- args[2]

# or
# runname_java <- "rhawkes-sim/seed7" 
# case <- 1

model <- c("Net", "Std", "Empty")[case]
seed <- c(32, 34, 36)[case]

perf.all <- list()
for(net in 1:10){
	runname = paste(runname_java, "-net", net, sep = "")

	Nitr <- 1000
	K <- 30

	bps <- readLines(paste("../data/", runname, seed, "bps-test-out.txt", sep = ""))
	bps <- strsplit(bps, ",")[[1]]
	bps <- as.numeric(as.character(bps))

	bps2 <- readLines(paste("../data/", runname, seed, "bps-second-test-out.txt", sep = ""))
	bps2 <- strsplit(bps2, ",")[[1]]
	bps2 <- as.numeric(as.character(bps2))

	w <- read.table(paste("../data/", runname, seed, "W-out.txt", sep = ""), 
		sep = ",", stringsAsFactors = FALSE)
	w3 <- array(0, dim = c(Nitr, K, K))
	for(i in 1:Nitr){
		for(j in 1:K){
			w3[i, j, ] <- as.numeric(w[(i-1)*K + j, ] )
		}
	}

	a <- read.table(paste("../data/", runname, seed,  "A-out.txt", sep = ""), 
		sep = ",", stringsAsFactors = FALSE)
	a3 <- array(0, dim = c(Nitr, K, K))
	for(i in 1:Nitr){
		for(j in 1:K){
			a3[i, j, ] <- as.numeric(a[(i-1)*K + j, ] )
		}
	}


	load(paste("../data/", runname, "-wk.rda", sep = ""))
	weff <- w3 * a3
	weff.m <- apply(weff, c(2, 3), mean)
	a.m <- apply(a3, c(2, 3), mean)
	A_true <- this$A

	# using effective W for prediction
	pred <- prediction(as.vector(weff.m), as.vector(A_true))
	# using A for prediction
	predA <- prediction(as.vector(a.m), as.vector(A_true))
	
	print(paste("AUC:", performance(pred, "auc")@y.values[[1]], 
						performance(predA, "auc")@y.values[[1]]))

	perfmean <- performance(pred,"tpr","fpr")
	perfmeanA <- performance(predA,"tpr","fpr")

	pdf(paste("../figures/", runname, "ROC.pdf", sep = ""))
	plot(perfmean)
	dev.off()
	pdf(paste("../figures/", runname, "AROC.pdf", sep = ""))
	plot(perfmeanA)
	dev.off()

	print(c(mean(bps), mean(bps2)))
	perf.all[[net]] <- list(pred = pred, 
							pred.A = predA, 
							perf = perfmean, 
							perf.A = perfmeanA,
							bps1 = bps, 
							bps2 = bps2, 
							events = this$events, 
							sum = sum(this$events),
							rates = this$rates)
}

if(model == "Net"){
	save(perf.all, file = paste("../data/", runname, "-model-perf.rda", sep = ""))

}else if(model == "Std"){
	perf.all.SH <- perf.all
	save(perf.all.SH, 
		file = paste("../data/", runname, "-model-perf-SH.rda", sep = ""))

}else if(model == "Empty"){
	perf.all.EH <- perf.all
	save(perf.all.EH, 
		file = paste("../data/", runname, "-model-perf-EH.rda", sep = ""))
}




