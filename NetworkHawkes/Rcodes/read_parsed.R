##
## This is the initial set of script used for generating data ans illustration
##
##
## File to read parsed stock price and make plots
## 

# load stock list
sp2 <- read.csv("../data/sp100/SP100-2.csv", header = FALSE)[,1]
sp2 <- as.character(sp2)
sp2 <- sp2[-which(sp2 == "GOOG")]

# load sector list
sps <- read.csv("../data/sp100/sp100-sector.csv", header = FALSE)
where <- match(sp2, sps[, 1])
sector <- as.character(sps[where, 2])
# table(sector)
# sp2[which(is.na(where))]

# read time file
times.raw <- scan("../data/sp100/time_points.txt", what = character(0), sep = ",")
times0 <- strptime(times.raw, format = "%Y-%m-%d|%H:%M:%S")
# remove 15 minutes time lag from the vector, add back 3 hrs to make EST
times0 <- times0 - 15 * 60 + 3 * 3600
nt <- length(times0)

prices0 <- matrix(0, 100, nt)
for(i in 1:length(sp2)){
	stock <- sp2[i]
	tmp <- scan(paste("../data/sp100/", stock, ".txt", sep = ""), what = double(4), sep = ",")
	prices0[i, ] <- tmp[1:nt]
}
rownames(prices0) <- sp2

# trading hour starts from 06:30
# trading hour ends at 
hourmin <- as.numeric(format(times0, format = "%H%M"))
period <- intersect(which(hourmin >= 930), which(hourmin < 1600))
# update prices and times
prices <- prices0[, period]
times <- times0[period]
# order stock by sector
prices <- prices[order(sector), ]
sector <- sector[order(sector)]
sector.num <- as.numeric(factor(sector))
# get hour position
hour <- as.numeric(format(times, format = "%H"))
hour.position <- c(1, which(diff(hour) != 0) + 1)
# format to EST
hour.label <- format(times[hour.position], format = "%H:%M")
day.position <- hour.position[which(hour.label == "09:30")]

changes <- t(apply(prices, 1, function(x){return(c(0, diff(x))/x)}))
changes.zero <- prices - prices[, 1]

events <- changes * 0
events[which(abs(changes) > 0.001)] <- 1
save(events, file = "../data/sp100/sp-events.rda")
namesdata <- list(names = rownames(events), sector = sector)
save(namesdata, file = "../data/sp100/sp-events-names.rda")
########################################################################
## save to file
tlist <- NULL
clist <- NULL
for(t in 1:dim(events)[2]){
	for(i in 1:dim(events)[1]){
		if(events[i, t] == 1){
			tlist <- c(tlist, t * 5)
			clist <- c(clist, i)			
		}
	}
}
filename1 <- "../data/sp100/event-time.txt"
filename2 <- "../data/sp100/event-proc.txt" 
file1 <- file(filename1)
writeLines(as.character(tlist), file1, sep = ",")
close(file1)
file2 <- file(filename2)
writeLines(as.character(clist), file2, sep = ",")
close(file2)
########################################################################


index <- which(events == 1)
# y index: which stock
y.index <- index %% 100
y.index[which(y.index == 0)] <- 100
# x index: which time
x.index <- trunc((index-1)/100)+1

library(reshape2)
library(RColorBrewer)
library(scales)
library(ggplot2)

##
## Correlation
##

events.cor <- cor(t(events))
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
pdf("corr-matrix.pdf", width = 12, height = 11)
ggplot(melt(events.cor), aes(Var1,Var2, fill=value)) + geom_raster() + scale_fill_gradientn(colours = myPalette(100)) + xlab("") + ylab("")
dev.off()

punch <- data.frame(Stock = factor(rownames(prices)[y.index], rownames(prices)), 
					Date =  format(times[x.index], "%D"),
					Time = times[x.index], 
					Sector = sector[y.index])

pdf("IT-events.pdf", width = 15, height = 6)
g <- ggplot(aes(x = Time, y = Stock), data = subset(punch, Sector == "Information Technology"))
g <- g + geom_point(aes(color = Stock), alpha = .8, size = 3)
g <- g + facet_grid(. ~ Date, scales = "free")
g <- g + scale_x_datetime(breaks = pretty_breaks(3))
g
dev.off()

pdf("google.pdf", width = 10, height = 6)
plot(changes[which(rownames(prices)=="GOOGL"), ], 
		type = "l", xaxt = 'n', 
		xlab = "Time", ylab = "Price Change (%)", main = "Google (GOOGL) stock price change (%) every 5 seconds", 
		ylim = c(-0.003, 0.003))
abline(h = c(-0.001, 0.001), lty = 2, col = "red", cex = 2)
abline(v = day.position[-1], lty = 2, col = "grey", cex = 2)
axis(1, at = hour.position, labels = hour.label)	
dev.off()

# plot(x = x.index, xaxt = 'n', 
# 	 y = y.index, yaxt = 'n',
# 	 col = sector.num[y.index],
# 	 cex = .2, 
# 	main="Event of stock price changing greater than 0.1%", 
# 	xlab="Time", ylab="Stock")
# abline(v = day.position[-1], lty = 1, col = "red", cex = 5)
# axis(1, at = hour.position, labels = hour.label)	
# axis(2, at = 1:100, labels = rownames(prices), las = 2, cex.axis = .5)

###
### Plots
###
# pdf("IT-ts.pdf", height = 30, width = 10)
# toplot <- which(sector == "Information Technology")
# par(mfrow = c(length(toplot), 1))
# for(i in 1:length(toplot)){
# 	plot(prices[toplot[i], period], type = "l", xaxt = 'n', 
# 		xlab = "Time", ylab = "Price", main = sp2[toplot[i]])
# 	abline(v = day.position[-1], lty = 2, col = "grey")
# 	axis(1, at = hour.position, labels = hour.label)
# }
# dev.off()


# pdf("IT-ts-change.pdf", height = 8, width = 10)
# toplot <- which(sector == "Information Technology")
# # toplot <- 1:100
# plot(changes[toplot[i], ], 
# 		type = "l", xaxt = 'n', 
# 		xlab = "Time", ylab = "Price", main = "IT sector price change (%) every 5 seconds", 
# 		ylim = c(-0.003, 0.003))
# 		# ylim = range(changes[toplot, ]))
# for(i in 1:length(toplot)){
# 	lines(changes[toplot[i], ], 
# 		lty = i, col = i)
# }
# abline(h = c(-0.001, 0.001), lty = 2, col = "red", cex = 2)
# abline(v = day.position[-1], lty = 2, col = "grey", cex = 2)
# axis(1, at = hour.position, labels = hour.label)	
# dev.off()


pdf("IT-ts-changefrom0.pdf", height = 8, width = 10)
toplot <- which(sector == "Information Technology")
plot(changes.zero[toplot[1], ], 
		type = "l", xaxt = 'n', 
		xlab = "Time", ylab = "Price", main = "IT sector price change (Apr 27 to May 1, 2015)", 
		ylim = range(changes.zero[toplot, ]), col = 'white')
for(i in 1:length(toplot)){
	lines(changes.zero[toplot[i], ], 
		lty = i, col = rainbow(length(toplot))[i])
	text(x=length(period), y=changes.zero[toplot[i], length(period)], 
		rownames(prices)[toplot[i]])
}
abline(v = day.position[-1], lty = 2, col = "grey")
axis(1, at = hour.position, labels = hour.label)
# legend("bottomleft", col = 1:length(toplot), sp2[toplot])
dev.off()