#### This file load the GIA prior and GPS data and prepare the data for the BHM model

#### 1 Load GIA prior
if(Sys.info()["sysname"] == "Windows"){
  ice6g <- read.table("Z:/WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5.txt", header = T)
}else if(Sys.info()["sysname"] == "Ubuntu"){
  ice6g <- read.table("~/globalmass/WP2-SolidEarth/GIAforwardModels/textfiles/GIA_Pel-6-VM5.txt", header = T)
}else{
  ice6g <- read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GIA/GIA_Pel-6-VM5.txt", header = T)
}

polycoords <- ice6g[,c(6:13, 6,7)] 
plist <- lapply(ice6g$ID, 
                function(x) Polygons(list(Polygon(cbind(lon = as.numeric(polycoords[x, c(1,3,5,7,9)]), 
                                                        lat = as.numeric(polycoords[x, c(2,4,6,8,10)])))), ID = x))
Plist <- SpatialPolygons(plist, proj4string = CRS("+proj=longlat"))


#### 2 Load GPS data
if(Sys.info()["sysname"] == "Windows"){
  GPSV4b <- read.table("Z:/WP2-SolidEarth/BHMinputs/GPS/GPS_v04b.txt", header = T)
}else if(Sys.info()["sysname"] == "Ubuntu"){
  GPSV4b <- read.table("~/globalmass/WP2-SolidEarth/GPS/NGL/BHMinputFiles/GPS_v04b.txt", header = T)
}else{
  GPSV4b <- read.table("/./projects/GlobalMass/WP2-SolidEarth/BHMinputs/GPS/GPS_v04b.txt", header = T)
}
