rm(list=ls())

#GOOD 
#https://tpwd.texas.gov/gis/
  
# BAD
# https://www.epa.gov/eco-research/level-iii-and-iv-ecoregions-continental-united-states 
  # edwards plateau is a different shape than Hawkes et al 2017.......


#install.packages("rgdal")
#install.packages("sf")
library(rgdal)
#library(sf)
library(tidyverse)
library(broom)
library(terra)
library(paletteer)
library(svglite)
shape <- readOGR(dsn = "~/Downloads/GouldEcoRegions", layer = "GouldEcoRegions") # this takes a long time to load in
setwd("~/Documents/Texas gradient")
class(shape)

# explore this file: projection?
proj4string(shape) # "+proj=lcc +lat_0=31.1666666666667 +lon_0=-100 +lat_1=27.4166666666667 +lat_2=34.9166666666667 +x_0=1000000 +y_0=1000000 +datum=NAD83 +units=m +no_defs"

shape@data # this lists out names of fields to ssubset by = "Name"
shape$Name
shape@proj4string
make_EPSG(shape)
 #"USA_Contiguous_Albers_Equal_Area_Conic_USGS_version"
# "North American Datum 1983"
  # still dont really understand this.
shape@data

# bring in sampling locations to go on texas map
points <- read.csv("TXGrad_RAPID_sitedata_latlong.csv", stringsAsFactors = FALSE)
points
str(points)
names(points)
head(points)




data.frame(shape.tx) # NA_L3NAME


# subset TX shapefile to only include edwards region.
shape.edwards <- shape[shape$Name=="Edwards Plateau",]

# change projection of edwards region from NAD83 LCC or whatever to longlat!!!!!
shape.edwards.longlat <- spTransform(shape.edwards, CRS("+proj=longlat") )

# change projection of texas from NAD83 LCC or whatever to longlat!
shape.longlat <- spTransform(shape, CRS("+proj=longlat") )


# THIS IS TO MAP EDWARDS ONTO TX
plot(shape , border="grey")
plot(shape.edwards, add = TRUE,  border="black")  # add=true put it on the same plot????


# THIS IS TO MAP POINTS ONTO EDWARDS
plot(shape.edwards.longlat , border ="black")
points(x = points$Lon, y= points$Lat, pch=19)
text(points$Long, points$Lat-0.1, labels = points$SiteCode, cex = 0.3)


# ALL 3
# want all of texass grey
# want edwards region to be black
# want sampling locations as dots 

#png("TXgraph.png", units="in", width=5, height=5, res=300)
svglite("Figure 1a.svg",  width=5, height=5)
plot(shape.longlat , border="grey80")
plot(shape.edwards.longlat , border = "black", add = TRUE)
points(x = points$Lon, y= points$Lat, pch=19, cex = 0.3)
dev.off()
  # add=true put it on the same plot????



################
# drought data?

setwd("~/Documents/Texas gradient")


drought.init <- read.csv("Re TX paper/TX_droughtsevereindex_2023_10_04.csv")
trt <- read.csv("Texas_Metadata.csv")

head(drought.init)
head(trt)
sites <- trt %>% select (c (Site, MAP))
sitemap <- sites[!duplicated(t(apply(sites,1,sort))),]
sitemap 

unique(drought.init$Year)
drought.recent <- drought.init %>% filter(Year<2018)

drought <- merge (drought.recent, sitemap, by = "Site")
names(drought )
head(drought )

drought$label = paste (drought$MAP, drought$Site, sep = ": ")

# merge up with 




MAPcolors <- c("#3B99B1FF",  "#35A4ABFF", "#3BA8A7FF", "#48AEA2FF", "#81BB95FF",
               "#A2C194FF","#B3C58FFF", "#CBC988FF",  "#EAC728FF", "#E9B41FFF" ,
              "#E9AB1CFF", "#E79C15FF",  "#E78D08FF" , "#E78200FF", "#E87600FF" , 
               "#EA6800FF" , "#EC5C00FF" , "#F04105FF", "#F5191CFF")
length(MAPcolors )
  



PDSI = ggplot(drought, aes(x=Year, y=PDSI)) +
  annotate("rect", fill = "blue4", alpha = 0.5, 
           xmin = -Inf, xmax = Inf, ymin = 4, ymax = Inf)+
  annotate("rect", fill = "blue", alpha = 0.5, 
           xmin = -Inf, xmax = Inf, ymin = 3, ymax = 4)+
  annotate("rect", fill = "royalblue3", alpha = 0.5, 
           xmin = -Inf, xmax = Inf, ymin = 2, ymax = 3)+
  annotate("rect", fill = "royalblue2", alpha = 0.5, 
           xmin = -Inf, xmax = Inf, ymin = 2, ymax = 3)+
  annotate("rect", fill = "royalblue1", alpha = 0.5, 
           xmin = -Inf, xmax = Inf, ymin = 1, ymax = 2)+
  annotate("rect", fill = "lightblue", alpha = 0.5, 
           xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 1)+
  annotate("rect", fill = "coral", alpha = 0.5, 
           xmin = -Inf, xmax = Inf, ymin = -1, ymax = -0.5)+
  annotate("rect", fill = "firebrick1", alpha = 0.5, 
           xmin = -Inf, xmax = Inf, ymin = -2, ymax = -1)+
  annotate("rect", fill = "firebrick2", alpha = 0.5, 
           xmin = -Inf, xmax = Inf, ymin = -3, ymax = -2)+
  annotate("rect", fill = "red3", alpha = 0.5, 
           xmin = -Inf, xmax = Inf, ymin = -4, ymax = -3)+
  annotate("rect", fill = "darkred", alpha = 0.5, 
           
           xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -4) + 
  geom_line(aes(color = as.character(MAP)) , linewidth = 1) + 
  scale_color_manual(values = rev(MAPcolors)) + 
  stat_summary(geom = "line", fun = "mean", linewidth =1.4) + 
  geom_vline(xintercept = 2015.9, color = "black", linetype = "dashed")  +   # fall
  geom_vline(xintercept = 2016.3, color = "black", linetype = "dashed")  +   # spring
  geom_vline(xintercept = 2016.6, color = "black", linetype = "dashed")  +   #summer
  geom_hline(yintercept = 0, color = "black")  +   #summer
  scale_x_continuous(breaks = seq(2010, 2017, by = 1)) + 
  theme(legend.key.height = unit(2.5, 'in')) +
  theme_classic(base_size = 15)

PDSI 

#ggsave("TX_PDSI_11624.png", plot = PDSI, device = "png",
 #   width = 7, height = 4.8, dpi = 300)

svglite("Figure 1b.svg",  width=7, height=4.8)
PDSI
dev.off()



