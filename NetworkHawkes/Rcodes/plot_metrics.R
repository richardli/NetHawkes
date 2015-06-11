## 
## script to load file
##
library(Hmisc)
library(ggplot2)

# load performance data from R
load("../data/rhawkes-sim/seed9-net10-model-perf.rda")
load("../data/rhawkes-sim/seed9-net10-model-perf-SH.rda")
load("../data/rhawkes-sim/seed9-net10-model-perf-GLM.rda")
load("../data/rhawkes-sim/seed9-net10-model-perf-Xcorr.rda")

#plot(perf.all[[1]]$perf)

##
## helper function to impute on a curve to be evaluated on given x vector
## assume x and xbasis in non-decreasing order, and the same end points(0, 1)
## note some x values repeat in data (note x = 0, need to take separately)
##
impute <- function(x, y, xbasis){
	newy <- rep(0, length(xbasis))
	for(i in 1:length(xbasis)){
		xx <- xbasis[i]
		if(xx == 0){
			newy[i] = y[i]
			next
		}
		# loop over x
		for(j in 2:length(x)){
			# find where it lies in x
			if(x[j-1] <= xx && x[j] >= xx){
				newy[i] = (xx - x[j-1])/(x[j] - x[j-1]) * (y[j] - y[j-1]) + y[j - 1] 
				next
			}
		}
	}
	return(list(x = xbasis, y = newy))
}

##
## function to extract multiple ROC curve from the list
## 
getROC <- function(perf.all, perfname){
	index <- which(names(perf.all[[1]]) == perfname)
	# number of network
	nn <- length(perf.all)
	
	# check all curve has the same x
	x.basis <- perf.all[[1]][[index]]@x.values[[1]]
	for(i in 1:10){
		newx <- perf.all[[i]][[index]]@x.values[[1]]
		if(length(x.basis) < length(newx)){x.basis <- newx}
	}

	# now permute all x vector onto xbasis and 
	# number of points on ROC curve
	np <- length(x.basis)
	x <- y <- matrix(0, nn, np)

	for(i in 1:10){
		newx <- perf.all[[i]][[index]]@x.values[[1]]
		newy <- perf.all[[i]][[index]]@y.values[[1]]
		# impute onto x.basis
		imp <- impute(newx, newy, x.basis)
		newx <- imp$x
		newy <- imp$y			
		# after imputation
		x[i, ] <- newx
		y[i, ] <- newy
	}

	# use only one x vector (FPR)
	xx <- x[1, ]
	# calculate mean, min, max of y by column
	yy <- apply(y, 2, mean)
	y0 <- apply(y, 2, quantile, 0.9)
	y1 <- apply(y, 2, quantile, 0.1)
	return(list(xx = xx, 
				yy = yy, 
				y0 = y0, 
				y1 = y1))
}

roc1 <- getROC(perf.all, "perf")
roc2 <- getROC(perf.all.SH, "perf")
roc3 <- getROC(perf.all.glm, "perf")
roc4 <- getROC(perf.all.xcorr, "perf")

xx1 <- roc1$xx
yy1 <- roc1$yy
xx2 <- roc2$xx
yy2 <- roc2$yy
xx3 <- roc3$xx
yy3 <- roc3$yy
xx4 <- roc4$xx
yy4 <- roc4$yy

pdf("../figures/Synthetic-1.pdf", width = 6, height = 6)
plot(xx1, yy1, type = "l", col = "darkgreen", 
	xlab = "FP rate", ylab = "TP rate", 
	main = "Synthetic Link Prediction ROC", lwd = 2, 
	xlim = c(0, 1), ylim = c(0, 1))
lines(xx2, yy2, lwd = 2, col = "darkblue")
lines(xx3, yy3, lwd = 2, col = "darkred")
lines(xx4, yy4, lwd = 2, col = "darkorange")
abline(a = 0, b = 1, lty = "longdash", lwd = 1.5)
legend("bottomright", 
		c("Net Hawkes", "Std. Hawkes", "GLM", "XCorr"), 
		lty = rep(1, 3), lwd = rep(2, 3),
		col = c("darkgreen", "darkblue", "darkred", "darkorange"))
dev.off()
#
#
# for(i in 1:10){plot(perf.all.glm[[1]]$perf)}
 

##
## Bits per spike
##
nn <- length(perf.all)
nitr <- length(perf.all[[1]]$bps1)
##
## Network Hawkes model
##
bps <- bps.glm <- bps.std <- matrix(0, nitr, nn)
for(i in 1:10){
	#  use Bits per second
	bps[, i] <- perf.all[[i]]$bps2
	bps.glm[, i] <- perf.all.glm[[i]]$bps1   # GLM function does not need to
	bps.std[, i] <- perf.all.SH[[i]]$bps2 
}
colnames(bps) <- colnames(bps.glm) <- colnames(bps.std) <- 1:10
nsave <- 1000

## 
## order the plot by decreasing event counts
##
sum_events <- rep(0, 10)
for(i in 1:10){
	sum_events[i] <- perf.all[[i]]$sum
}
network.order <- order(sum_events, decreasing = FALSE)
network.name <- rep(0, 10)
for(i in 1:10) network.name[network.order[i]] <- i 

bpsdata1 <- data.frame(network = network.name, 
				  mean = apply(bps[(nitr-nsave):nitr, ], 2, mean),
				  # high = apply(bps[(nitr-nsave):nitr, ],  2, quantile, 0.975),
				  # low = apply(bps[(nitr-nsave):nitr, ],  2,  quantile, 1-0.975), 
				  high = apply(bps[(nitr-nsave):nitr, ], 2, mean) + apply(bps[(nitr-nsave):nitr, ],  2, sd),
				  low = apply(bps[(nitr-nsave):nitr, ], 2, mean) - apply(bps[(nitr-nsave):nitr, ],  2, sd), 
				  model = rep("Net.Hawkes", 10))
##
## GLM
##
bpsdata2 <- data.frame(network = network.name, 
				  mean = apply(bps.glm[(nitr-nsave):nitr, ], 2, mean),
				  high = NA,
				  low = NA,
				  model = rep("GLM", 10))

##
## Standard Hawkes
##
bpsdata3 <- data.frame(network = network.name, 
				  mean = apply(bps.std[(nitr-nsave):nitr, ], 2, mean),
				  high = apply(bps.std[(nitr-nsave):nitr, ], 2, mean) + apply(bps.std[(nitr-nsave):nitr, ],  2, sd),
				  low = apply(bps.std[(nitr-nsave):nitr, ], 2, mean) - apply(bps.std[(nitr-nsave):nitr, ],  2, sd), 
				  # low = apply(bps.std[(nitr-nsave):nitr, ],  2,  quantile, 1-0.975), 
				  # high = apply(bps.std[(nitr-nsave):nitr, ],  2, quantile, 0.975),
				  model = rep("Std.Hawkes", 10))


bpsdata <- rbind(bpsdata1, bpsdata2, bpsdata3)
bpsdata$network <- as.factor(bpsdata$network)
bpsdata$model <- factor(bpsdata$model, c("Net.Hawkes", "Std.Hawkes", "GLM"))

model.mean <- aggregate(mean~model, FUN=mean, data = bpsdata)
model.mean[2,2] / model.mean[1, 2]
model.mean[3,2] / model.mean[1, 2]

scale <- matrix(0, 10, 2)
for(i in 1:10){
	scale[i, 1] <- bpsdata2$mean[which(bpsdata2$network == i)] / bpsdata1$mean[which(bpsdata1$network == i)]
	scale[i, 2] <- bpsdata3$mean[which(bpsdata3$network == i)] / bpsdata1$mean[which(bpsdata1$network == i)]
}
m1 <- mean(range(scale[, 1]))
m2 <- mean(range(scale[, 2]))
d1 <- max(scale[, 1]) - m1
d2 <- max(scale[, 2]) - m2
print(paste("GLM vs Net Hawkes", "Mean:", model.mean[3,2]/model.mean[1, 2]))
print(paste("Range: ", m1, "plus-minus", d1))
print(paste("Std vs Net Hawkes", "Mean:", model.mean[2,2]/model.mean[1, 2]))
print(paste("Range: ", m2, "plus-minus", d2))

pdf("../figures/Synthetic-2.pdf", width = 10, height = 5)
g <- ggplot(aes(x = network, y = mean, fill = model), data = bpsdata) 
g <- g + geom_bar(position= position_dodge(width=.9), stat = "identity") 
g <- g + geom_errorbar(aes(ymin =low, ymax = high), 
						position =  position_dodge(width=0.9), width=0.25)
g <- g + theme_bw()
g <- g + xlab("Network") + ylab("Bits/second")
g <- g + ggtitle("Synthetic Predictive Log Likelihood")
g
dev.off()

pdf("../figures/Synthetic-2s.pdf", width = 8, height = 6)
g <- ggplot(aes(x = network, y = mean, fill = model), data = bpsdata) 
g <- g + geom_bar(position= position_dodge(width=.9), stat = "identity") 
g <- g + geom_errorbar(aes(ymin =low, ymax = high), 
						position =  position_dodge(width=0.9), width=0.25)
g <- g + theme_bw()
g <- g + xlab("Network") + ylab("Bits/second")
g <- g + ggtitle("Synthetic Predictive Log Likelihood")
g
dev.off()

# order <- order(apply(bps[(nitr-200):nitr, ], 2, mean), decreasing = TRUE)
# barplot(apply(bps[(nitr-200):nitr, order], 2, mean), col = "darkgreen")
# errbar(x = 1:10, y = apply(bps[(nitr-200):nitr, order],  2, mean), 
# 	  yplus = apply(bps[(nitr-200):nitr, order],  2, quantile, 0.975), 
# 	  yminus = apply(bps[(nitr-200):nitr, order],  2,  quantile, 1-0.975), 
# 	  add = TRUE, col = "darkgreen")