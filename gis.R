#map is downloaded from http://www.gadm.org/country
#http://biogeo.ucdavis.edu/data/gadm2.8/rds/BGD_adm1.rds  #up to division
#http://biogeo.ucdavis.edu/data/gadm2.8/rds/BGD_adm2.rds  #up to district
#http://biogeo.ucdavis.edu/data/gadm2.8/rds/BGD_adm3.rds  #up to thana

listOfPackages <- c("sp", "maptools", "dplyr","ggplot2")
if(length(which(!listOfPackages %in% installed.packages()))){
  install.packages(listOfPackages[!listOfPackages %in% installed.packages()])
}

library(sp) # need for this to work with SpatialPolygons
library(maptools) #converting SpatialPolygons into drawable files
library(dplyr)  # merging data
library(ggplot2)

# read map information
# choose file yourself
# BD<- readRDS(choose.files())
# Load file from specific path 
# BD <- readRDS("E:/University/2017.05.Summer/CUSRA/R/CUSRA Qi/BGD_adm1.rds") # load division
# BD <- readRDS("E:/University/2017.05.Summer/CUSRA/R/CUSRA Qi/BGD_adm2.rds") # load district
 BD <- readRDS("E:/University/2017.05.Summer/CUSRA/R/CUSRA Qi/BGD_adm3.rds") # load thana
 BD2 <- fortify(BD) #convert to plotable format
 
# for division
#df <- data.frame(id = rownames(BD@data), regionName = BD@data$NAME_1)
# for district
#df <- data.frame(id = rownames(BD@data), regionName = BD@data$NAME_2)
# for thana
#df <- data.frame(id = rownames(BD@data), regionName = BD@data$NAME_3)
# NOTE: some district Name mismatch

# For CSI
# warning comes from duplicates, no worries
# warning comes from NA value due to mis-print mis-match of original data
df <- left_join(df, csiDF, by = "regionName")
# optional, without this line, region with NA will be grey
#df$regionCsi <- replace(df$regionCsi, is.na(df$regionCsi), 0)
# CANNOT join by regionNmae because NO name for the BD2 data frame
BGD <- left_join(BD2,df, by = "id")

# For Smooth CSI
# warning comes from duplicates, no worries
# warning comes from NA value due to mis-print mis-match of original data
#sdf <- left_join(df, scsiDF, by = "regionName")
# optional, without this line, region with NA will be grey
#sdf$regionCsi <- replace(sdf$regionCsi, is.na(sdf$regionCsi), 0)
# CANNOT join by regionName because NO name for the BD2 data frame
#SBGD <- left_join(BD2,sdf, by = "id")

#create informative map for csi
ggplot(BGD, aes(x = long, y = lat, group = group, fill = regionCsi)) +
  theme_minimal()+ 
  geom_path() +
  geom_polygon(color = "white") +
  scale_fill_distiller(palette = "Spectral") +
  coord_map() + 
  ggtitle("Bangladesh Thana CSI") +
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(fill = "CSI")

#create informative map for smooth csi
ggplot(SBGD, aes(x = long, y = lat, group = group, fill = regionCsi)) +
  theme_minimal()+ 
  geom_path() +
  geom_polygon(color = "white") +
  scale_fill_distiller(palette = "Spectral") +
  coord_map() + 
  ggtitle("Bangladesh District Smooth CSI") +
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(fill = "Smooth CSI")

#####################################
# For single division
# Only a single division/district
DV <- subset(BD, BD@data$NAME_1 %in% c("Dhaka"))
DV2 <- fortify(DV)
df <- data.frame(id = rownames(DV@data), regionName = DV@data$NAME_3)
df <- left_join(df, csiDF, by = "regionName")
#df$regionCsi <- replace(df$regionCsi, is.na(df$regionCsi), 0)
DVM <- left_join(DV2,df, by = "id")

#create informative map for csi for Dhaka ONLY
ggplot(DVM, aes(x = long, y = lat, group = group, fill = regionCsi)) +
  theme_minimal()+ 
  geom_path() +
  geom_text(aes(label = round(regionCsi, 3))) +
  geom_polygon(color = "white") +
  scale_fill_distiller(palette = "Spectral") +
  coord_map() + 
  ggtitle("Dhaka CSI") +
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(fill = "CSI")
