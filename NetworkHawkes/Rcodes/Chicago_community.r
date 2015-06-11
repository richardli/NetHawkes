#install.packages('rgdal',repos="http://www.stats.ox.ac.uk/pub/RWin")
#install.packages("gpclib")
library(maptools)
library(rgdal)
library(ggmap)

# download shapefile from:
#   https://data.cityofchicago.org/api/geospatial/cauq-8yn6?method=export&format=Shapefile

# setwd accordingly

overlay <- readOGR("../data/Chicago", "CommAreas")
overlay <- spTransform(overlay, CRS("+proj=longlat +datum=WGS84"))
overlay <- fortify(overlay, region="COMMUNITY") # it works w/o this, but I figure you eventually want community names

data$COMAREA <- mapvalues(data$COMAREA, from = "Grand Blvd", to = "Grand Boulevard")
data$COMAREA <- mapvalues(data$COMAREA, from = "West Elston", to = "West Elsdon")
data$COMAREA <- mapvalues(data$COMAREA, from = "Ashburn Gresham", to = "Auburn Gresham")

overlay.f <- fortify(overlay, region="COMMUNITY") 
overlay2 <- merge(overlay.f, overlay@data, by.x = "id", by.y = "COMMUNITY")

count <- table(data$COMAREA)
index <- match(tolower(overlay2$id), tolower(names(count)))
overlay2$counts <- as.numeric(count[index]) 
					


overlay2 = SpatialPolygonsDataFrame(p, overlay2)
location <- unlist(geocode('4135 S Morgan St, Chicago, IL 60609'))+c(0,.02)

gmap <- get_map(location=location, zoom = 10, maptype = "terrain", source = "google", col="bw")


Map <- ggplot(sport.f, aes(long, lat, group = group, fill = Partic_Per)) + geom_polygon() + 
    coord_equal() + labs(x = "Easting (m)", y = "Northing (m)", fill = "% Sport Partic.") + 
    ggtitle("London Sports Participation")

gg <- ggmap(gmap)
gg <- gg + geom_polygon(data=overlay2, aes(x=long, y=lat, group=group, fill =counts), color="red")
gg <- gg + coord_map()
gg <- gg + theme_bw()
gg <- gg + scale_fill_gradient(low = "white", high = "black")
gg
