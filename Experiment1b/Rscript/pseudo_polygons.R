#### Pseudo observation polygons

#### Use the 14 GIA solutions to find polygons for pseudo observations.


### Assume Gaussian, find means that are within than 2 times sd
GIA_temp$zero <- abs(GIA_temp$trend) < GIA_temp$std * 2
map_temp <- ggplot(data=GIA_temp) + geom_raster(aes(x = x_center, y = y_center, fill = zero)) + 
  coord_fixed() + xlab("Longitude") + ylab("Latitude") + 
  scale_x_continuous(limits=c(0,359),  expand = c(0, 0)) + scale_y_continuous(limits=c(-89.5,89.5),  expand = c(0, 0)) 
map_temp <- map_temp + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                      colour="grey", fill = NA, alpha = 0.5) + ggtitle(GIA_name[i])
map_zero1[[i]] <- map_temp
}

grid.arrange(map_zero1[[1]], map_zero1[[2]],map_zero1[[3]],map_zero1[[4]], ncol = 2)
grid.arrange(map_zero1[[5]], map_zero1[[6]],map_zero1[[7]],map_zero1[[8]], ncol = 2)
grid.arrange(map_zero1[[9]], map_zero1[[10]],map_zero1[[11]],map_zero1[[12]], ncol = 2)
grid.arrange(map_zero1[[13]], map_zero1[[14]], ncol = 2)

## Comments: NOt a good selection criteria for "zero certain areas"; better to use the absolute mean value only.

### 2 deviate from zero in absolute values, say 0.5
map_zero2 <- list()
for (i in 1:14){
  GIA_temp <- GIA_priors[[i]]
  ### Assume Gaussian, find means that are within than 2 times sd
  GIA_temp$zero <- abs(GIA_temp$trend) < 0.5
  map_temp <- ggplot(data=GIA_temp) + geom_raster(aes(x = x_center, y = y_center, fill = zero)) + 
    coord_fixed() + xlab("Longitude") + ylab("Latitude") + 
    scale_x_continuous(limits=c(0,359),  expand = c(0, 0)) + scale_y_continuous(limits=c(-89.5,89.5),  expand = c(0, 0)) 
  map_temp <- map_temp + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                      colour="grey", fill = NA, alpha = 0.5) + ggtitle(GIA_name[i])
  map_zero2[[i]] <- map_temp
}

grid.arrange(map_zero2[[1]], map_zero2[[2]],map_zero2[[3]],map_zero2[[4]], ncol = 2)
grid.arrange(map_zero2[[5]], map_zero2[[6]],map_zero2[[7]],map_zero2[[8]], ncol = 2)
grid.arrange(map_zero2[[9]], map_zero2[[10]],map_zero2[[11]],map_zero2[[12]], ncol = 2)
grid.arrange(map_zero2[[13]], map_zero2[[14]], ncol = 2)

### 3a use the mean of the 14 data sets -- absolute different
GIA_prior_ensemble <- GIA_priors[[1]]
GIA_prior_ensemble$trend <- rowMeans(sapply(GIA_priors, "[[", "trend"))

GIA_prior_ensemble$zero <- abs(GIA_prior_ensemble$trend) < 0.5
map_ensemble <- ggplot(data=GIA_prior_ensemble) + geom_raster(aes(x = x_center, y = y_center, fill = zero)) + 
  coord_fixed() + xlab("Longitude") + ylab("Latitude") + 
  scale_x_continuous(limits=c(0,359),  expand = c(0, 0)) + scale_y_continuous(limits=c(-89.5,89.5),  expand = c(0, 0)) 
map_ensemble <- map_ensemble + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                    colour="grey", fill = NA, alpha = 0.5) + ggtitle("14 GIA ensemble mean")
map_ensemble


### 3b use the mean of the 14 data sets -- statistical difference
GIA_prior_ensemble$sd <- sqrt((rowSums(sapply(GIA_priors, function(x) (x$std)^2)))/14)

GIA_prior_ensemble$zero2 <- abs(GIA_prior_ensemble$trend) < GIA_prior_ensemble$std
map_ensemble2 <- ggplot(data=GIA_prior_ensemble) + geom_raster(aes(x = x_center, y = y_center, fill = zero2)) + 
  coord_fixed() + xlab("Longitude") + ylab("Latitude") + 
  scale_x_continuous(limits=c(0,359),  expand = c(0, 0)) + scale_y_continuous(limits=c(-89.5,89.5),  expand = c(0, 0)) 
map_ensemble2 <- map_ensemble2 + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                            colour="grey", fill = NA, alpha = 0.5) + ggtitle("14 GIA ensemble mean v2")
map_ensemble2
