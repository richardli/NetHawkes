###########################################################
# Get all BPS
###########################################################
library(ggplot2)
library(grid)
library(reshape)
library(RColorBrewer)
library(scales)

runname <- rep("", 4)
runname[1] <- "sp100/event13" # Net Hawkes
runname[2] <- "sp100/event15"	# Latent Hawkes
runname[3] <- "sp100/event17" # Std Hawkes
runname[4] <- "sp100/event19" # LGCP background

allBPS <- matrix(0, 4, 3)
colnames(allBPS) <- c("Mean", "Lower", "Upper")
rownames(allBPS) <- c("Net. Hawkes (Erdos-Renyi)", 
	"Net. Hawkes (Latent Distance)", 
	"Std. Hawkes", "Indep. LGCP")

for(i in 1:4){
	bps <- readLines(paste("../data/", runname[i], "bps-test-out.txt", sep = ""))
	bps <- strsplit(bps, ",")[[1]]
	bps <- as.numeric(as.character(bps))
	allBPS[i, 1] <- mean(bps)
	allBPS[i, 2] <- quantile(bps, 0.025)
	allBPS[i, 3] <- quantile(bps, 1 - 0.025)
}
allBPS
allBPS <- data.frame(allBPS)
allBPS$model <- rownames(allBPS)
save(allBPS, file = "../data/sp-BPS.rda")

# bpsdata <- melt(allBPS)
# colnames(bpsdata) <- c("model", "")  
pdf("../figures/StockLatent-3.pdf", width = 6.5, height = 3.5)
g <- ggplot(aes(x = model, y = Mean, fill = model), data = allBPS) 
g <- g + geom_bar(position= position_dodge(width=.5), stat = "identity", width = 0.8) 
g <- g + geom_errorbar(aes(ymin =Lower, ymax = Upper), 
						position =  position_dodge(width=0.5), width=0.25)
g <- g + theme_bw()
g <- g + theme(axis.text.x = element_blank())
g <- g + ylab("Bits/spike")
g <- g + ggtitle("Predictive Log Likelihood")
# g <- g + theme(legend.position="bottom")
g
dev.off()

###########################################################
# Illustration for latent space Hawkes Model
###########################################################
runname <- "sp100/event15"

# T.train <- length(bps)
# T <- length(bps) * 2
# pdf(paste("../data/", runname, "bps-time.pdf", sep = ""), width = 10, height = 6)
# par(mfrow = c(1, 2))
# plot((T.train + 1) : T, bps, type = "l", main = "bits per spike", 
# 	xlab = "iteration", ylab = "Bits/spike")
# plot((T.train + 1) : T, bps2, type = "l", main = "bits per second",
# 		xlab = "iteration", ylab = "Bits/second")
# dev.off()
Nitr <- 1000
K <- 100
###########################################################
#  baseline
###########################################################
baseline <- read.table(paste("../data/", runname, "baseline-out.txt", sep = ""), 
	sep = ",", stringsAsFactors = FALSE)
dim(baseline)
M <- 41
lambda3 <- array(0, dim = c(Nitr, K, M))
for(i in 1:Nitr){
	for(j in 1:K){
		lambda3[i, j, ] <- as.numeric(baseline[(i-1)*K + j, ] )
	}
}

###########################################################
#  extract one network sample for illustration
###########################################################
offset <- Nitr - 1
###########################################################
# read W
w <- read.table(paste("../data/", runname, "W-out.txt", sep = ""), 
	sep = ",", stringsAsFactors = FALSE)
oneW <- matrix(0, K, K)
for(j in 1:K){
	oneW[j, ] <- as.numeric( w[j+offset*100, ] )
}

###########################################################
# read A
a <- read.table(paste("../data/", runname, "A-out.txt", sep = ""), 
	sep = ",", stringsAsFactors = FALSE)

oneA <- matrix(0, K, K)
for(j in 1:K){
	oneA[j, ] <- as.numeric(a[j+offset*100, ] )
}

###########################################################
# read P
p <- read.table(paste("../data/", runname, "P-out.txt", sep = ""), 
	sep = ",", stringsAsFactors = FALSE)

# get one P matrix
oneP <- matrix(0, K, K)
for(j in 1:K){
	oneP[j, ] <- as.numeric(p[j+offset*100, ] )
}

###########################################################
# get all A, W, P matrix
p3 <- a3 <- w3 <- array(0, dim = c(Nitr, K, K))
for(i in 1:Nitr){
	for(j in 1:K){
		p3[i, j, ] <- as.numeric(p[(i-1)*K + j, ] )
		a3[i, j, ] <- as.numeric(a[(i-1)*K + j, ] )
		w3[i, j, ] <- as.numeric(w[(i-1)*K + j, ] )	
	}
}
meanP <- apply(p3, c(2, 3), mean)

max(abs(eigen(oneW * oneA)$values))

data = list(W = oneW, A = oneA, P = oneP, meanP = meanP, 
	allP = p3, allA = a3, allW = w3, alllambda = lambda3)
save(data, 
	file = paste("../data/", runname, "oneNet.rda", sep = ""))


###########################################################
##
## Get 2D embedding of P
##
###########################################################
load(paste("../data/", runname, "oneNet.rda", sep = ""))
A <- data$A
W <- data$W
P <- data$P
lambda <- data$alllambda
pdf("../figures/StockLatent-5.pdf", width = 6, height = 6)
rate.meank <- apply(lambda[1:1000,,], c(2, 3), mean)
rate.mean <- apply(rate.meank, 2, mean)
rate.high <- rate.mean + 2 *  apply(rate.meank, 2, sd)
rate.low <- rate.mean - 2 * apply(rate.meank, 2, sd)
# rate.high <- apply(rate.meank, 2, quantile, 0.9)
# rate.low <- apply(rate.meank, 2, quantile, 0.2)
plot(1:41, rate.mean, type = "l", lwd = 2, 
	ylim = c(0, max(c(rate.mean, rate.high, rate.low))), 
	xaxt = 'n', xlab = NULL, ylab = "rate", yaxt = 'n')
lines(1:41, rate.high, col = "gray", lty = 2)
lines(1:41, rate.low, col = "gray", lty = 2)
abline(v = 41 / 4, col = "gray", lty = 2)
abline(v = 41 / 4 * 2, col = "gray", lty = 2)
abline(v = 41 / 4 * 3, col = "gray", lty = 2)
axis(1, at = c(5, 15, 25, 35), c("Mon", "Tue", "Wed", "Thur"))
dev.off()

pdf("../figures/StockLatent-6.pdf", width = 10, height = 6)
rate.meank <- apply(lambda[1:1000,,], c(2, 3), mean)
plot(1:41, rate.meank[1, ], type = "l", lwd = 2, col = "white",
	ylim = c(0, quantile(rate.meank, 0.99)), 
	xaxt = 'n', xlab = "", ylab = "", yaxt = 'n', 
	main = "Inferred background rate")
for(i in 1:100){
	lines(1:41, rate.meank[i, ], col = "gray", lty = 1)	
}
lines(1:41, rate.mean, col = "red", lwd = 2)
abline(v = 41 / 4, col = "blue", lty = 2, lwd = 2)
abline(v = 41 / 4 * 2, col = "blue", lty = 2, lwd = 2)
abline(v = 41 / 4 * 3, col = "blue", lty = 2, lwd = 2)
axis(1, at = c(5, 15, 25, 35), c("Mon", "Tue", "Wed", "Thur"))
dev.off()

###########################################################
##
## Get 2D embedding of P
##
###########################################################

# P <- data$meanP
load("../data/sp100/sp-events-names.rda")
names <- namesdata$names
sector <- namesdata$sector

# P.toplot <- apply(data$allP[1:1000, , ], c(2, 3), mean)
P.toplot <- data$allP[1, , ]

P.toplot [which(P.toplot  == 0)] <- 1e-9
colnames(P.toplot) <- rownames(P.toplot) <- names

rho <- max(P.toplot) * 2
tau <- 1
dist <- -log(P.toplot / rho) * tau
diag(dist) <- 0
range(dist)

toplot <- c("Information Technology",  "Financials", "Energy", 
	"Health Care", "Consumer Discretionary",  "Industrials")
selected <- setdiff(which(sector %in% toplot), which(names %in% c("nodeleted")))
dist  <- dist[selected, selected]

loc <- cmdscale(dist)
# library(MASS)
# loc <- isoMDS(dist, k=2, p = 2)$points
## just checking
dist.fit <- matrix(0, dim(dist)[1], dim(dist)[1])
for(i in 1:dim(dist)[1]){
	for(j in 1:dim(dist)[1]){
		dist.fit[i, j] <- sum((loc[i, ] - loc[j, ])^2)^.5
	}
}
plot(dist, dist.fit)

## just checking
#
# x0 <- runif(100, 0, 5)
# y0 <- runif(100, 4, 5)
# dist.true <- matrix(0, 100, 100)
# for(i in 1:100){
# 	for(j in 1:100){
# 		dist.true[i, j] <- sqrt((x0[i] - x0[j])^2 + (y0[i] - y0[j])^2)
# 	}
# }
# loc <- cmdscale(dist.true)
# dist.fit <- matrix(0, 100, 100)
# for(i in 1:100){
# 	for(j in 1:100){
# 		dist.fit[i, j] <- sum((loc[i, ] - loc[j, ])^2)^.5
# 	}
# }
# plot(dist.true, dist.fit)

x <- loc[,1]
y <- -loc[,2]
plot(x, y, type="n", xlab="", ylab="", main="cmdscale(eurodist)")
text(x, y, rownames(loc), cex=0.8)

distdata <- data.frame(first = x, second = y, 
	textlabel = y + 0.25,
	names = names[selected], 
	sector = sector[selected])

###########################################################
# first plot
###########################################################
pdf("../figures/StockLatent-1.pdf", width = 9, height = 6)
g <- ggplot(aes(x = first, y = second), data = distdata)
g <- g + geom_point(aes(shape = sector, col = sector, fill = sector), 
					size = 3)
# g <- g + geom_text(aes(label = names, col = sector), size = 3.5)
g <- g + scale_colour_discrete(l=60)
g <- g + scale_shape_manual(values=1:nlevels(distdata$sector))
g <- g + theme_bw() 
g <- g + geom_text(aes(y = textlabel, label = names, col = sector), size = 3.5, 
		data = data.frame(subset(distdata, abs(x) > 2.1)), show_guide  = F )
g <- g + geom_text(aes(y = textlabel, label = names, col = sector), size = 3.5, 
		data = data.frame(subset(distdata, abs(y) > 2.1)), show_guide  = F )
g <- g + xlab("Latent Dimension 1")
g <- g + ylab("Latent Dimension 2")
g <- g + ggtitle("Inferred Embedding of Stocks (six largest sectors)")
g
grid.edit("geom_point.points", grep = TRUE, gp = gpar(lwd = 3))
dev.off()

pdf("../figures/StockLatent-1s.pdf", width = 9, height = 6)
g <- ggplot(aes(x = first, y = second), data = distdata)
g <- g + geom_point(aes(shape = sector, col = sector, fill = sector), 
					size = 3)
# g <- g + geom_text(aes(label = names, col = sector), size = 3.5)
g <- g + scale_colour_discrete(l=60)
g <- g + scale_shape_manual(values=1:nlevels(distdata$sector))
g <- g + theme_bw() 
g <- g + xlim(c(-2, 2)) + ylim(c(-2, 2))
g <- g + xlab("Latent Dimension 1")
g <- g + ylab("Latent Dimension 2")
g <- g + ggtitle("Inferred Embedding of Stocks (six largest sectors)")
g
grid.edit("geom_point.points", grep = TRUE, gp = gpar(lwd = 3))
dev.off()

###########################################################
# second plot, separate for each sector
###########################################################
pdf("../figures/StockLatent-2.pdf", width = 10, height = 8)
g <- ggplot(aes(x = first, y = second), data = distdata )
# g <- g + geom_point(aes(shape = sector, col = sector, fill = sector), 
# 					size = 3)
g <- g + geom_text(aes(label = names, col = sector), size = 3.5)
g <- g + scale_shape_manual(values=1:nlevels(distdata$sector))
g <- g + theme_bw()
g <- g + scale_colour_discrete(guide = FALSE)
g <- g + facet_wrap(~ sector)
g <- g + xlim(c(-2, 2)) + ylim(c(-2, 2))
g <- g + xlab("Latent Dimension 1")
g <- g + ylab("Latent Dimension 2")
g <- g + ggtitle("Inferred Embedding of Stocks (six largest sectors)")
g
dev.off()

###########################################################
# third plot, effective W
###########################################################
##
weff <- data$allA * data$allW
weff.m <- apply(weff, c(2, 3), mean)
mat <- weff.m
mat <- as.matrix(mat)
colnames(mat) <- rownames(mat) <- names
Stock1 <- Stock2 <- rep("", length(mat))
value <- rep(0, length(mat))
K <- dim(mat)[1]
for(i in 1:K){
	for(j in 1:K){
		Stock1[(i-1) * K + j] <- names[i]
		Stock2[(i-1) * K + j] <-  names[j]
		value[(i-1) * K + j] <- mat[i, j]
	}
}
Stock1 <- factor(Stock1, names)
Stock2 <- factor(Stock2, names)
matdata <- data.frame(stock1 = Stock1, stock2 = Stock2, effect.W = value)

# myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
pdf("../figures/StockLatent-4.pdf", width = 12, height = 11)
g <- ggplot(matdata, aes(x=stock1,y=stock2, fill=effect.W)) 
g <- g + geom_raster() 
g <- g + scale_fill_gradient2( low = "white", high = "blue")
# g <- g + scale_fill_gradientn(colours = myPalette(100)) 
g <- g + xlab("") + ylab("")
g
dev.off()


#######################################################
eig <- eigen(A * W)
print(abs(eig$values)[1:6])
