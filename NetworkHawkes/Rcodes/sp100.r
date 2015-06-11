###################################################################################
##
## Script to get stock price data using R
##	
##
###################################################################################
#install.packages("quantmod")
library(quantmod)
# two versions
sp <- read.csv("SP100/SP100.csv", header = FALSE)[,1]
sp2 <- read.csv("SP100/SP100-2.csv", header = FALSE)[,1]

sp <- as.character(sp)
sp2 <- as.character(sp2)

sp <- sp[-which(sp == "BRK.B")]
# use only class A google
sp2 <- sp2[-which(sp2 == "GOOG")]


interval <- 1
N <- (6.5 * 60 * 60 +  30 * 60) / interval
dates <- c("2015-04-21", "2015-04-22", "2015-04-23", "2015-04-24")
for(d in dates){
	# wait until open
	file <- paste("SP100/", d, "-", interval, "s.txt", sep = "")
	Sys.sleep( as.numeric(
	as.POSIXct( paste(d, "06:30:00 PDT") )-Sys.time(), units="secs"))

	t <- Sys.time()
	write.table(getQuote(sp2,
	            what=yahooQF(c("Bid","Ask","Last Trade (Price Only)"))),
	            file=file, append=FALSE, col.names = F)
	for (i in 1 : N){
	    left <- t + interval - Sys.time()
		if(left > 0) Sys.sleep(left)
		t <- Sys.time()
	    write.table(getQuote(sp2,
	            what=yahooQF(c("Bid","Ask","Last Trade (Price Only)"))),
	            file=file, append=TRUE, col.names = F)
	    print(t)
	}	
}

# save.image(paste("SP100/", interval, "s-space.RData", sep = ""))

