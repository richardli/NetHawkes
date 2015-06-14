# reinstall this if fortify gives error again
# install.packages('rgdal',repos="http://www.stats.ox.ac.uk/pub/RWin")
# install.packages("gpclib")
library(maptools)
library(rgdal)
library(gpclib)
library(ggmap)
library(foreign)
library(plyr)

gpclibPermit()

## 	download shapefile from:
##	https://data.cityofchicago.org/api/geospatial/cauq-8yn6?method=export&format=Shapefile

overlay <- readOGR("../data/chicago/geo", "CommAreas")
overlay <- spTransform(overlay, CRS("+proj=longlat +datum=WGS84"))
# overlay <- fortify(overlay, region="COMMUNITY") 
overlay.f <- fortify(overlay, region="COMMUNITY") 
overlay2 <- merge(overlay.f, overlay@data, by.x = "id", by.y = "COMMUNITY")


## download data from:
## http://www.icpsr.umich.edu/icpsrweb/NACJD/studies/6399?q=chicago&amp;archive=NACJD
##
## read data 
## for cleaner version, use STATA file
data = read.dta("../data/chicago/ICPSR_06399/DS0001/06399-0001-Data.dta")
data <- data[which(data$GANG == "Yes"), ]
clusterdata <- read.csv("../data/chicago/comcluster.csv")
## clean data
##
## extract year, month, time, community
## output: 1. community
##         2. community identifier
##	       3. which month starting from 0 (Jan 1980 = 0)
##         4. random date starting from 0 (Jan 1, 1980 = 0)
##	       5. cluster
##	       6. cluster identifier
##
# clean community
data$COMAREA <- mapvalues(data$COMAREA, from = "Grand Blvd", to = "Grand Boulevard")
data$COMAREA <- mapvalues(data$COMAREA, from = "West Elston", to = "West Elsdon")
data$COMAREA <- mapvalues(data$COMAREA, from = "Ashburn Gresham", to = "Auburn Gresham")
com_name <- as.character(data$COMAREA)
# to be removed 
toremove <- NULL
# clean numerical community code
com_list <- clusterdata$Community
com <- match(com_name, com_list)
cluster_name <- clusterdata$Cluster[com]
cluster_list <- c("red", "green", "yellow", "blue")
cluster <- match(cluster_name, cluster_list)
toremove <- union(toremove, which(com_name == "Missing"))

# clean date
year <- 1900 + data$INJYEAR
toremove <- union(toremove, which(year < 1980 ))

month <- as.character(data$INJMONTH)
toremove <- union(toremove, which(month == "Missing"))

# real month diff
month.names <- format(ISOdate(2000,1:12,1),"%B")
month.num <- match(month, month.names)
month.count <- (month.num - 1) + (year - 1980) * 12

# random date
date <- sample(size = length(month), x = 1:28, replace = TRUE)
dym <- paste(date, month, year, sep = "-")
date <- as.Date(dym, "%d-%B-%Y")
date.count <- date - as.Date("1-January-1980", "%d-%B-%Y")
date.count <- as.numeric(date.count)

# combine
dataset <- data.frame(com = com, 
					  com_name = com_name, 
					  cluster = cluster, 
					  cluster_name = cluster_name,
					  year = year, 
					  month = month, 
					  month.count = month.count, 
					  date = date, 
					  date.count = date.count)
dataset <- dataset[-toremove, ]
dataset <- dataset[order(dataset$date.count), ]
save(dataset, file = "../data/chicago/clean.rda")
save(com_list, file = "../data/chicago/comlist.rda")


# save for java 
filename0 <- "../data/chicago/event0-proc.txt"  ##  month + community
filename1 <- "../data/chicago/event0-time.txt"  ##  
filename2 <- "../data/chicago/event1-proc.txt"  ##  random day + community
filename3 <- "../data/chicago/event1-time.txt"  ##  
filename4 <- "../data/chicago/event2-proc.txt"  ##  month + cluster
filename5 <- "../data/chicago/event2-time.txt"  ##  
filename6 <- "../data/chicago/event3-proc.txt"  ## random day + cluster
filename7 <- "../data/chicago/event3-time.txt"  ##  

# write time
file1 <- file(filename1)
writeLines(as.character(dataset$month.count), file1, sep = ",")
close(file1)
file1 <- file(filename5)
writeLines(as.character(dataset$month.count), file1, sep = ",")
close(file1)

file1 <- file(filename3)
writeLines(as.character(dataset$date.count), file1, sep = ",")
close(file1)
file1 <- file(filename7)
writeLines(as.character(dataset$date.count), file1, sep = ",")
close(file1)

# write events
file2 <- file(filename0)
writeLines(as.character(dataset$com), file2, sep = ",")
close(file2)
file2 <- file(filename2)
writeLines(as.character(dataset$com), file2, sep = ",")
close(file2)

file2 <- file(filename4)
writeLines(as.character(dataset$cluster), file2, sep = ",")
close(file2)
file2 <- file(filename6)
writeLines(as.character(dataset$cluster), file2, sep = ",")
close(file2)

# counting
count <- table(data$COMAREA)
index <- match(tolower(overlay2$id), tolower(names(count)))
overlay2$counts <- as.numeric(count[index]) 
overlay2$cluster <- clusterdata$Cluster[match(tolower(overlay2$id), tolower(clusterdata$Community))]

# get map					
location <- unlist(geocode('4135 S Morgan St, Chicago, IL 60609'))+c(0,.02)
gmap <- get_map(location=location, zoom = 10, maptype = "terrain", source = "google", col="bw")

# plot
gg <- ggmap(gmap)
gg <- gg + geom_polygon(data=overlay2, aes(x=long, y=lat, group=group, fill =cluster), color="black")
gg <- gg + coord_map()
gg <- gg + theme_bw()
gg <- gg + scale_fill_manual(values = c("blue2","green3", "red2", "orange2"))
gg <- gg + guides(fill=FALSE)
# gg <- gg + scale_fill_gradient(low = "white", high = "black")
gg
