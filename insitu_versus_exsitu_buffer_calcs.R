### Author: Emily Beckman
### Date: April 28, 2020; updated May 6, 2021
### Funding: Institute of Museum and Library Services [award #MA-30-18-0273-18]

################################################################################

### DESCRIPTION:
  ## This script creates circular buffers around in situ points (wild occurrence
	#		records) and ex situ points (wild locations where seeds were collected for
	#		cultivation in botanic gardens or arboreta) to calculate the geographic
	#		(area within buffers) and ecological (number of ecoregions within buffers)
	#		diversity conserved in ex situ living collections. A map showing 50km
	#		buffers around in situ and ex situ points is also created.
	## These data are presented in Fredlock, Zumwalde, et al. (2021)

### INPUTS:
?????	# In situ points: curated by Sean Hoban lab.
	# Ex situ points: source localities of the seedlings used in the genetic
	#		study, from collection trips in 2017 and 2019.
	## All points are in one file (BeckHob_QHOccur_Vetted_ExSituMarked.csv),
	#		with a column marking if the point is represented ex situ.
	## Shapefiles of U.S. EPA Level III and Level IV Ecoregions
	#		Level III (us_eco_l3) -> https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/us_eco_l3.zip
	#		Level IV (us_eco_l4)-> https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/us/us_eco_l4_state_boundaries.zip
	## Shapefile of state boundaries, for mapping
	#		https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_state_20m.zip

### OUTPUTS:
  ## Table of geographic and ecological coverage based on three buffer sizes
	#		(10km, 50km, and 100km) and both ecoregion levels (III and IV)
	## Leaflet map with 50 km buffers around in situ and ex situ points
	## Leaflet map with ecoregions and in situ buffers and points

################################################################################

#################
### LIBRARIES ###
#################

## [code chunk from Shannon M. Still]
my.packages <- c("leaflet","raster","sp","rgeos","dplyr","rgdal","Polychrome",
	"RColorBrewer")
#install.packages (my.packages) # turn on to install current versions
lapply(my.packages, require, character.only=TRUE)
rm(my.packages)

#################
### FUNCTIONS ###
#################

# format text in cell for output table
format.cell <- function(ex_result,in_result,final_result){
	cell <- paste0(round(final_result,2),"%","\n",
								 "(",format(round(ex_result,0),format="d",big.mark=",")," / ",
								 format(round(in_result,0),format="d",big.mark=","),")")
	return(cell)
}

# create buffers around points, using specified projection
create.buffers <- function(df,radius,pt_proj,buff_proj){
	# select coordinate columns
	latlong <- df %>% select(longitude,latitude)
	# turn occurrence point data into a SpatialPointsDataFrame
	sp_df <- SpatialPointsDataFrame(latlong, df, proj4string = pt_proj)
	# reproject SpatialPointsDataFrame to specified projection
	proj_df <- spTransform(sp_df,buff_proj)
	# place buffer around each point
	buffers <- buffer(proj_df,width=radius,dissolve=T)
	# return buffer polygons
	return(buffers)
}

# create buffers around in situ and ex situ spatial points, calculate areas,
#		then compare to calculate percent coverage
compare.buff.area <- function(insitu,exsitu,radius,pt_proj,buff_proj){
	# create buffers
	buffer_insitu <- create.buffers(insitu,radius,pt_proj,buff_proj)
	buffer_exsitu <- create.buffers(exsitu,radius,pt_proj,buff_proj)
	# calculate buffer area
	print(paste("Based on ",radius/1000," km radius..."))
	area_exsitu <- buffer_exsitu@polygons[[1]]@area/1000000
	print(paste("Area covered by ex situ buffers:", round(area_exsitu,0),"km²"))
	area_insitu <- buffer_insitu@polygons[[1]]@area/1000000
	print(paste("Area covered by in situ buffers:", round(area_insitu,0),"km²"))
	# calculate difference between in situ and ex situ buffer areas (% coverage)
	area_diff_percent <- (area_exsitu/area_insitu)*100
	print(paste0("Percent geographic coverage: ", round(area_diff_percent,2), "%"))
	txt <- format.cell(area_exsitu,area_insitu,area_diff_percent)
	return(txt)
}

# create data frame with ecoregion data extracted for area covered by buffers
intersect.eco.buff <- function(df,radius,pt_proj,buff_proj,eco){
	# create buffers
	buffers <- create.buffers(df,radius,pt_proj,buff_proj)
	# make sure ecoregions are in same projection as buffers
	eco_proj <- spTransform(eco,buff_proj)
	# intersect buffers with ecoregions
	buff_join_eco <- raster::intersect(buffers,eco_proj)
	return(buff_join_eco)
}

# create data frame with Level III ecoregion data extracted for area covered by
#		buffers, for both in situ and ex situ points, then compare ecoregion counts
compare.ecol3.count <- function(insitu,exsitu,radius,pt_proj,buff_proj,eco){
	# create data frame of ecoregion-buffer intersection
	eco_insitu <- intersect.eco.buff(insitu,radius,pt_proj,buff_proj,eco)
	eco_exsitu <- intersect.eco.buff(exsitu,radius,pt_proj,buff_proj,eco)
	# count number of ecoregions under buffers
	print(paste("Based on ",radius/1000," km radius..."))
	count_exsitu <- nrow(eco_exsitu@data %>% distinct(US_L3CODE))
	print(paste0("Number of ecoregions under ex situ buffers: ",count_exsitu))
	count_insitu <- nrow(eco_insitu@data %>% distinct(US_L3CODE))
	print(paste0("Number of ecoregions under in situ buffers: ",count_insitu))
	# calculate difference in number of ecoregions
	eco_diff_percent <- (count_exsitu/count_insitu)*100
	print(paste0("Percent ecological coverage: ", round(eco_diff_percent,2), "%"))
	txt <- format.cell(count_exsitu,count_insitu,eco_diff_percent)
	return(txt)
}

# create data frame with Level IV ecoregion data extracted for area covered by
#		buffers, for both in situ and ex situ points, then compare ecoregion counts
compare.ecol4.count <- function(insitu,exsitu,radius,pt_proj,buff_proj,eco){
	# create data frame of ecoregion-buffer intersection
	eco_insitu <- intersect.eco.buff(insitu,radius,pt_proj,buff_proj,eco)
	eco_exsitu <- intersect.eco.buff(exsitu,radius,pt_proj,buff_proj,eco)
	# count number of ecoregions under buffers
	print(paste("Based on ",radius/1000," km radius..."))
	count_exsitu <- nrow(eco_exsitu@data %>% distinct(US_L4CODE))
	print(paste0("Number of ecoregions under ex situ buffers: ",count_exsitu))
	count_insitu <- nrow(eco_insitu@data %>% distinct(US_L4CODE))
	print(paste0("Number of ecoregions under in situ buffers: ",count_insitu))
	# calculate difference in number of ecoregions
	eco_diff_percent <- (count_exsitu/count_insitu)*100
	print(paste0("Percent ecological coverage: ", round(eco_diff_percent,2), "%"))
	txt <- format.cell(count_exsitu,count_insitu,eco_diff_percent)
	return(txt)
}

################################################################################
# A) Set up workspace
################################################################################

# working directory
local_dir <- "/Users/emily/Documents/GitHub/Quercus_havardii_GeoEco_exsitu_conservation"

# define projection of points
#		WGS84
wgs.proj <- sp::CRS(SRS_string="EPSG:4326")
# define projection for calculations (meters must be the unit)
#		U.S. National Atlas Equal Area Projection
#		Warning shouldn't cause any errors; metadata says they're "overly cautious"
aea.proj <- sp::CRS(SRS_string="EPSG:2163")

################################################################################
# B) Read in data
################################################################################

### POINT DATA

# in situ and ex situ points
pts <- read.csv(file.path(local_dir,"BeckHob_QHOccur_Vetted_ExSituMarked.csv"),
	as.is=T, na.strings=c("","NA"))
# in situ subset
insitu <- pts
# ex situ subset
exsitu <- pts[which(!is.na(pts$ex_situ)),]

# create subsets for eastern population and western population
insitu_e <- insitu %>% filter(east_west == "E")
insitu_w <- insitu %>% filter(east_west == "W")
exsitu_e <- exsitu %>% filter(east_west == "E")
exsitu_w <- exsitu %>% filter(east_west == "W")

### POLYGONS

# read in shapefile of Level III ecoregions
ecol3 <- readOGR(file.path(local_dir,
	"us_eco_l3/us_eco_l3.shp"))

# read in shapefile of Level IV ecoregions
ecol4 <- readOGR(file.path(local_dir,
	"us_eco_l4_state_boundaries/us_eco_l4.shp"))

# read in shapefile of state boundaries
state_bound <- readOGR(file.path(local_dir,
	"cb_2018_us_state_20m/cb_2018_us_state_20m.shp"))

################################################################################
# C) Calculate geographic coverage (buffer areas)
################################################################################

### OVERALL

# calculate area based on 100 kilometer buffers
geo_coverage_100 <- compare.buff.area(insitu,exsitu,100000,wgs.proj,aea.proj)
# calculate area based on 50 kilometer buffers
geo_coverage_50 <- compare.buff.area(insitu,exsitu,50000,wgs.proj,aea.proj)
# calculate area based on 10 kilometer buffers
geo_coverage_10 <- compare.buff.area(insitu,exsitu,10000,wgs.proj,aea.proj)

### EAST ONLY

# calculate area based on 100 kilometer buffers
geo_coverage_100e <- compare.buff.area(insitu_e,exsitu_e,100000,wgs.proj,aea.proj)
# calculate area based on 50 kilometer buffers
geo_coverage_50e <- compare.buff.area(insitu_e,exsitu_e,50000,wgs.proj,aea.proj)
# calculate area based on 10 kilometer buffers
geo_coverage_10e <- compare.buff.area(insitu_e,exsitu_e,10000,wgs.proj,aea.proj)

### WEST ONLY

# calculate area based on 100 kilometer buffers
geo_coverage_100w <- compare.buff.area(insitu_w,exsitu_w,100000,wgs.proj,aea.proj)
# calculate area based on 50 kilometer buffers
geo_coverage_50w <- compare.buff.area(insitu_w,exsitu_w,50000,wgs.proj,aea.proj)
# calculate area based on 10 kilometer buffers
geo_coverage_10w <- compare.buff.area(insitu_w,exsitu_w,10000,wgs.proj,aea.proj)

################################################################################
# D) Calculate ecological coverage using buffers (ecoregion counts)
################################################################################

##
#### LEVEL III ECOREGIONS
##

### OVERALL

# count ecoregions under 100 km buffers
eco3_coverage_100 <- compare.ecol3.count(insitu,exsitu,100000,wgs.proj,aea.proj,ecol3)
# count ecoregions under 50 km buffers
eco3_coverage_50 <- compare.ecol3.count(insitu,exsitu,50000,wgs.proj,aea.proj,ecol3)
# count ecoregions under 10 km buffers
eco3_coverage_10 <- compare.ecol3.count(insitu,exsitu,10000,wgs.proj,aea.proj,ecol3)

### EAST ONLY

# count ecoregions under 100 km buffers
eco3_coverage_100e <- compare.ecol3.count(insitu_e,exsitu_e,100000,wgs.proj,aea.proj,ecol3)
# count ecoregions under 50 km buffers
eco3_coverage_50e <- compare.ecol3.count(insitu_e,exsitu_e,50000,wgs.proj,aea.proj,ecol3)
# count ecoregions under 10 km buffers
eco3_coverage_10e <- compare.ecol3.count(insitu_e,exsitu_e,10000,wgs.proj,aea.proj,ecol3)

### WEST ONLY

# count ecoregions under 100 km buffers
eco3_coverage_100w <- compare.ecol3.count(insitu_w,exsitu_w,100000,wgs.proj,aea.proj,ecol3)
# count ecoregions under 50 km buffers
eco3_coverage_50w <- compare.ecol3.count(insitu_w,exsitu_w,50000,wgs.proj,aea.proj,ecol3)
# count ecoregions under 10 km buffers
eco3_coverage_10w <- compare.ecol3.count(insitu_w,exsitu_w,10000,wgs.proj,aea.proj,ecol3)

##
#### LEVEL IV ECOREGIONS
##

### OVERALL

# count ecoregions under 100 km buffers
eco4_coverage_100 <- compare.ecol4.count(insitu,exsitu,100000,wgs.proj,aea.proj,ecol4)
# count ecoregions under 50 km buffers
eco4_coverage_50 <- compare.ecol4.count(insitu,exsitu,50000,wgs.proj,aea.proj,ecol4)
# count ecoregions under 10 km buffers
eco4_coverage_10 <- compare.ecol4.count(insitu,exsitu,10000,wgs.proj,aea.proj,ecol4)

### EAST ONLY

# count ecoregions under 100 km buffers
eco4_coverage_100e <- compare.ecol4.count(insitu_e,exsitu_e,100000,wgs.proj,aea.proj,ecol4)
# count ecoregions under 50 km buffers
eco4_coverage_50e <- compare.ecol4.count(insitu_e,exsitu_e,50000,wgs.proj,aea.proj,ecol4)
# count ecoregions under 10 km buffers
eco4_coverage_10e <- compare.ecol4.count(insitu_e,exsitu_e,10000,wgs.proj,aea.proj,ecol4)

### WEST ONLY

# count ecoregions under 100 km buffers
eco4_coverage_100w <- compare.ecol4.count(insitu_w,exsitu_w,100000,wgs.proj,aea.proj,ecol4)
# count ecoregions under 50 km buffers
eco4_coverage_50w <- compare.ecol4.count(insitu_w,exsitu_w,50000,wgs.proj,aea.proj,ecol4)
# count ecoregions under 10 km buffers
eco4_coverage_10w <- compare.ecol4.count(insitu_w,exsitu_w,10000,wgs.proj,aea.proj,ecol4)

################################################################################
# E) Create summary table of results
################################################################################

summary_tbl <- data.frame(
	Method = c("Geographic","Geographic","Geographic",
						 "Ecological, Level III","Ecological, Level III","Ecological, Level III",
						 "Ecological, Level IV","Ecological, Level IV","Ecological, Level IV"),
	Buffer_km = c("100","50","10",
								"100","50","10",
								"100","50","10"),
	Overall = c(geo_coverage_100,geo_coverage_50,geo_coverage_10,
							eco3_coverage_100,eco3_coverage_50,eco3_coverage_10,
							eco4_coverage_100,eco4_coverage_50,eco4_coverage_10),
	East = c(geo_coverage_100e,geo_coverage_50e,geo_coverage_10e,
					 eco3_coverage_100e,eco3_coverage_50e,eco3_coverage_10e,
					 eco4_coverage_100e,eco4_coverage_50e,eco4_coverage_10e),
	West = c(geo_coverage_100w,geo_coverage_50w,geo_coverage_10w,
					 eco3_coverage_100w,eco3_coverage_50w,eco3_coverage_10w,
					 eco4_coverage_100w,eco4_coverage_50w,eco4_coverage_10w),
	stringsAsFactors=F)
summary_tbl

# write table to csv
write.csv(summary_tbl, file.path(local_dir,
	"Qhavardii_GeoEco_exsitu_SummaryTbl.csv"),row.names = F)

################################################################################
# F) Map points and buffers
################################################################################

# create map to visualize buffer and point data
geo_map <- leaflet() %>%
	## background
	addProviderTiles("Esri.WorldGrayCanvas",
		options = providerTileOptions(maxZoom = 10)) %>%
	## state boundaries
	addPolygons(data = state_bound,
		fillOpacity = 0, color = "#757575", weight = 1.5, opacity = 1) %>%
	## in situ buffers
	addPolygons(data = create.buffers(insitu,50000,wgs.proj,wgs.proj),
		smoothFactor = 0.5,	weight = 2, color = "#de5f5f", fillOpacity = 0.35) %>%
	## in situ points
	addCircleMarkers(data = insitu, lng = ~longitude, lat = ~latitude,
		popup = ~paste("In situ:", Pop),
		radius = 4, fillOpacity = 0.9, stroke = F, color = "#de5f5f") %>%
	## ex situ buffers
	addPolygons(data = create.buffers(exsitu,50000,wgs.proj,wgs.proj),
		smoothFactor = 0.5, weight = 2, color = "#2426bd", fillOpacity = 0.35) %>%
	## ex situ points
	addCircleMarkers(data = exsitu, lng = ~longitude, lat = ~latitude,
		popup = ~paste("Ex situ:", Pop),
		radius = 4, fillOpacity = 0.9, stroke = F, color = "#2426bd") %>%
	## title, legend, scalebar, and zoom
	addControl("Quercus havardii in situ distribution and wild collection sites of ex situ accessions",
		position = "topright") %>%
	addLegend(labels =
		c(paste0("Geographic range (occurrence points and 50 km buffers)"),
			paste0("Populations sampled for ex situ (collection locations and 50 km buffers)")),
		colors = c("#de5f5f","#2426bd"), title = "Legend",
		position = "bottomright", opacity = 0.8) %>%
	addScaleBar(position = "bottomleft",
		options = scaleBarOptions(maxWidth = 150)) %>%
	setView(-106, 36, zoom = 6)

# view map
geo_map

# save map
htmlwidgets::saveWidget(geo_map, file = "Quercus_havardii_buffer_map_50km.html")

################################################################################
# G) Map ecoregions
################################################################################

# transform ecoregion polygons to WGS84, for mapping
ecol3_wgs <- spTransform(ecol3,wgs.proj)
ecol4_wgs <- spTransform(ecol4,wgs.proj)

# select only ecoregions that are within the buffers; otherwise map takes a
#		long time to load in browser and colors are not distinct enough
inter <- intersect.eco.buff(insitu,50000,wgs.proj,wgs.proj,ecol3)
codes <- unique(inter@data$US_L3CODE)
ecol3_sel <- ecol3_wgs[ecol3_wgs@data$US_L3CODE %in% codes,]
inter <- intersect.eco.buff(insitu,50000,wgs.proj,wgs.proj,ecol4)
codes <- unique(inter@data$US_L4CODE)
ecol4_sel <- ecol4_wgs[ecol4_wgs@data$US_L4CODE %in% codes,]

# create ecoregion color palette; can run display.brewer.all() to see options
ecol3_pal <- colorFactor(palette = "Set2",domain = ecol3_sel@data$US_L3NAME,
	reverse = F, na.color = "white")
ecol4_pal <- colorFactor(palette = "Set2",domain = ecol4_sel@data$US_L4NAME,
	reverse = F, na.color = "white")

# comment out/uncomment based on which level you're mapping
level <- "III"
level <- "IV"

# create map to visualze species distribution and ecoregions
# 	uncomment the section for Level III or Level IV ecoregions, depending
#			on which you're wanting to map
eco_map <- leaflet() %>%
	## background
	addProviderTiles("Esri.WorldGrayCanvas",#"Esri.WorldShadedRelief",#"Esri.WorldTerrain",
		options = providerTileOptions(maxZoom = 10)) %>%
	addPolygons(
	## EPA Level III ecoregions
		#data = ecol3_sel,
		#fillColor = ~ecol3_pal(ecol3_sel@data$US_L3NAME), fillOpacity = 0.6,
		#color = ~ecol3_pal(ecol3_sel@data$US_L3NAME), weight = 1, opacity = 1) %>%
	## EPA Level IV ecoregions
		data = ecol4_sel,
		fillColor = ~ecol4_pal(ecol4_sel@data$US_L4NAME), fillOpacity = 0.6,
		color = ~ecol4_pal(ecol4_sel@data$US_L4NAME), weight = 0.5, opacity = 1) %>%
	## state boundaries
	addPolygons(data = state_bound,
		fillOpacity = 0, color = "#969696", weight = 1.5, opacity = 1) %>%
	## in situ buffers
	addPolygons(data = create.buffers(insitu,50000,wgs.proj,wgs.proj),
		smoothFactor = 0.5,	weight = 2, color = "#363636", fillOpacity = 0.2) %>%
	## in situ points
	addCircleMarkers(data = insitu, lng = ~longitude, lat = ~latitude,
		popup = ~paste("In situ:", Pop),
		radius = 4, fillOpacity = 0.9, stroke = F, color = "#363636") %>%
	## title, legend, scalebar, and zoom
	addControl(paste0("Quercus havardii in situ distribution and U.S. EPA Level ",level," Ecoregions"),
		position = "topright") %>%
	addLegend(labels =
		c("Geographic range of Quercus havardii (occurrence points and 50 km buffers)",
			paste0("Background colors show U.S. EPA Level ",level," Ecoregions that fall within buffers")),
		colors = c("#363636","#d6ccb4"), title = "Legend", position = "bottomright",
		opacity = 0.8) %>%
	addScaleBar(position = "bottomleft",
		options = scaleBarOptions(maxWidth = 150)) %>%
	setView(-104, 36, zoom = 6)

# view map
eco_map

# save map
#htmlwidgets::saveWidget(eco_map, file = "Quercus_havardii_ecoregions_L3.html")
#htmlwidgets::saveWidget(eco_map, file = "Quercus_havardii_ecoregions_L4.html")
