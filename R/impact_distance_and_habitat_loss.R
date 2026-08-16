# Calculate impact distance from joint likelihood approach and diver displacement

setwd("C:/Users/RPerez/Documents/BC/manuscrito 2")

main.crs = "+init=epsg:32632 +proj=utm +zone=32 +datum=WGS84 +units=km +no_defs +ellps=WGS84 +towgs84=0,0,0"
library(raster)
library(stars)
library(fasterize)
library(distanceto)
library(cubelyr)
library(sf)
library(ggplot2)


# load files
# farms
farms = readRDS("farms18_spring_byyear_updt.rds")

# before
bef <- raster('pred_0108_mean_nosc.tif')
# after
aft <- raster('pred_1721_mean_nosc.tif')
# sign_incr
#decr <- raster('pred_change_975.tif')

decr_tif=read_stars("pred_change_975.tif")
decr=st_as_sf(decr_tif)
plot(decr)

bbox_decr <- st_bbox(decr)

farms_north <- farms[farms$cluster == "West of Sylt" | farms$cluster == "Butendiek", ]
farms_south <- farms[farms$cluster == "UMBO", ]

plot(farms_north)

farms_north_sf <- st_as_sf(farms_north)
dist_farms_north <- distance_raster(farms_north_sf, cellsize= 1,extent= bbox_decr)
plot(dist_farms_north)

farms_south_sf <- st_as_sf(farms_south)
dist_farms_south <- distance_raster(farms_south_sf, cellsize= 1,extent= bbox_decr)
plot(dist_farms_south)

decr0 <- subset(decr,pred_change_975.tif <= 0)
plot(decr0)
decr_line <- st_as_sf(decr0)

decr_union <- st_union(decr_line)
plot(decr_union)

# Convert the multipolygon to a polygon
p <- st_cast(decr_union, "POLYGON")

# Extract the  south and north polygons from the polygon object
p_north <- p[6]
plot(p_north)

plot(p)

p_south <- p[3]
plot(p_south)
plot(farms_south, add=TRUE)

#north
p_spdf_north <- as(p_north, "Spatial") 
bound_north <- as(p_spdf_north, 'SpatialLines')
plot(bound_north)
cells_north <- raster::extract(dist_farms_north, bound_north)
cells_north.df <- data.frame(cells_north[[1]])

#south
p_spdf_south <- as(p_south, "Spatial") 
bound_south <- as(p_spdf_south, 'SpatialLines')
plot(bound_south)
cells_south <- raster::extract(dist_farms_south, bound_south)
cells_south.df <- data.frame(cells_south[[1]])

library(dplyr)
df1 <- cells_north.df %>% rename(distance = cells_north..1..) %>% mutate(area = "north")
df2 <- cells_south.df %>% rename(distance = cells_south..1..) %>% mutate(area = "south")

# combine the two data frames
combined <- bind_rows(df1, df2)
saveRDS(combined, "diver_displacement_area.rds")


ggplot(combined, aes(y= distance, x= area)) + 
  geom_boxplot(fill = "lightblue") +
  theme_bw(base_size = 14)

combined <-readRDS("diver_displacement_area.rds")
library(ggplot2)
library(ggstatsplot)
library(tidyverse)

plt <- ggbetweenstats(
  data = combined,
  x = area,
  y = distance
)+
  labs(x = "Area", y = "Distance") +
  theme_minimal (base_size = 16)
plt
t.test(north$distance, south$distance)


north <- subset(combined, area=="north")
south <- subset(combined, area=="south")
ci <- t.test(north$distance)$conf.int
cat("95% CI: [", round(ci[1], 2), ", ", round(ci[2], 2), "]\n")

ci <- confint(lm(south$distance ~ 1), level = 0.95)

# Print the confidence interval
cat("95% CI: [", round(ci[1], 2), ", ", round(ci[2], 2), "]\n")
# calculate habitat loss

#library(terra)

plot(decr0)
dens_north <- st_intersection(decr0,p_north)
plot(dens_north)
sum(st_area(dens_north))


dens_north$loss <- exp(dens_north$pred_change_975.tif)
summary(dens_north$loss)

sum(st_area(dens_north)) # 1170.577 sq km
sum(dens_north$loss) #220 individuals
sum(dens_north$loss) / sum(st_area(dens_north)) #decrease of 0.188 divers/sqkm



dens_south <- st_intersection(decr0,p_south)
plot(dens_south)
sum(st_area(dens_south)) #879.1082 sq km

dens_south$loss <- exp(dens_south$pred_change_975.tif)
summary(dens_south$loss)

sum(dens_south$loss) #170 indiv
sum(dens_south$loss) / sum(st_area(dens_south)) #decrease of 0.194 divers/sq km


plot(aft)
aft_sqkm <- aft/4.7
plot(aft_sqkm)

plot(p,add=TRUE)

plot(bef)
plot(p,add=TRUE)


diff <- aft-bef
plot(diff)
plot(p,add=TRUE)
diff_north_num <- mask(diff,p_spdf_north)
diff_south_num <- mask(diff,p_spdf_south)

cellStats(diff_north_num, stat="sum") 

cellStats(diff_south_num, stat="sum") 


cell_size<-area(diff_north_num, na.rm=TRUE, weights=FALSE)
cell_size<-cell_size[!is.na(cell_size)]
raster_area<-length(cell_size)*median(cell_size)

  plot(diff_north_num)

diff_south_num <- mask(diff,p_spdf_south)


plot(diff_south_num)
1387/4.7
