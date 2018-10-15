library(maptools)
library(ggplot2)
library(broom)

Lands <- readOGR(dsn = "C:/ZSwork/glbm/Experiment2/ne_110m_land", layer = "ne_110m_land")
lLands <- as(Lands, "SpatialPolygons")
llLands <- slot(lLands, "polygons")
lpoly <- lapply(llLands, slot, "Polygons") #At this point we have a list of SpatialPolygons
## Change the coords form (-180, 180) to (0 ,360)
coords <- matrix(nrow=0, ncol=2)
for (i in seq_along(lpoly)){
  for (j in seq_along(lpoly[[i]])) {
    crds <- rbind(slot(lpoly[[i]][[j]], "coords"), c(NA, NA)) #the NAs are used to separate the lines
    coords <- rbind(coords, crds)
  }
}
coords[,1] <- ifelse(coords[,1]<0, coords[,1] + 360, coords[,1]) # Because here your levelplot will be ranging from 0 to 360°

lands_data <- tidy(lands) 

## The ogriginal data is in [-180, 180] and we want to change it to [0, 360]
## First split the polygons cross the zero degree
## identify which polygons
lands@polygons
spLand <- SpatialPolygons(lands@polygons)
plot(spLand, col = "grey")
summary(spLand)

library(maptools)
b <- readShapeSpatial("110-m_land.shp") #I used here a world map from Natural Earth.



