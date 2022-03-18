# Adapted from https://github.com/trhermes/Ust-Biyke-R-Map

library(dplyr)
library(ggmap)          # ggmap() get_stamenmap()
library(ggrepel)        # geom_label_repel()
library(rgdal)          # readOGR()
library(broom)          # tidy()
library(ggsn)           # scalebar()
library(grDevices)      # cairo_pdf()
library(rgeos)          # gCentroid()
library(sf)             # st_coordinates()
library(ggspatial)      # annotation_north_arrow()

# Read in site location data in csv format: site name, latitude, longitude - (decimal degrees)
arch_sites <- read.csv("data/input/site_coords.csv", sep=",")

dir.create("tmp")

# Download shapefiles of Natural Earth Data admin0 borders, and Diva GIS KG admin2
download.file(
  "https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip", 
  "tmp/admin_0.zip", "auto")
download.file(
  "https://biogeo.ucdavis.edu/data/diva/adm/KGZ_adm.zip", 
  "tmp/admin_2.zip", "auto")

# Unzip files
unzip("tmp/admin_0.zip", exdir = "tmp")
unzip("tmp/admin_2.zip", exdir = "tmp")

# Read in shapefiles
borders <- readOGR(dsn="tmp/ne_10m_admin_0_countries.shp", stringsAsFactors = FALSE)
districts <- readOGR(dsn="tmp/KGZ_adm2.shp", stringsAsFactors = FALSE) 
j_o_border <- districts[districts$NAME_2 %in% c("Djety-Oguz"), ]

# add  centroid of Jeti Oguz area to site list
all_sites <- gCentroid(j_o_border,byid=TRUE) %>% 
  st_as_sf() %>% 
  st_coordinates() %>%
  as_tibble() %>% 
  mutate(site = "Jeti Oguz") %>% 
  rename(long = X,
         lat = Y) %>% 
  bind_rows(arch_sites)

# Convert into tabular format and group features for mapping
borders2 <- tidy(borders, group=group) 
j_o_border2 <- tidy(j_o_border, group=group) 

# Manually specify bounding box for map
# One could also generate the bounds based on min and max lat/long in sites
map_borders <- c(bottom  = min(all_sites$lat) - 4, 
                 top     = max(all_sites$lat) + 1.5,
                 left    = min(all_sites$long) - 12,
                 right   = max(all_sites$long) + 13)

# Download map tiles
map <- get_stamenmap(map_borders, zoom = 9, maptype = "terrain-background", force=T)

# Text labels for plotting map features
mountains <- tibble::tribble(
  ~name,       ~lat, ~long, ~rotate,
  "Pamir Mnts.",     38.1,  73.1,   0,
  "Tian Shan Mnts.", 40.75,  74.5, 20,
  "Tian Shan Mnts.", 43,    84,   20,
  "Dzhungar Mnts.",  44.55,  82.5, -15,
  "Altai Mnts.",     49.5,  84.8, -39.8
)
countries <- tibble::tribble(
  ~name,        ~lat, ~long,
  "Kyrgyzstan", 42.1, 73.3,
  "Uzbekistan", 40.4,  66.2,
  "Tajikistan", 38.9,  70,
  "Kazakhstan", 46,    70, 
  "China",      45.5,  85.5,
  "Mongolia",   49,    90.4, 
  "Russia",     50.2,  86.8 
)

# Map it
figure1 <- ggmap(map) +
  geom_path(data=borders2, aes(x=long, y=lat, group = group), size=0.5, alpha = 0.5) +
  geom_polygon(data=j_o_border2, aes(x=long, y=lat, group = group), size=2, alpha = 0.5) +
  geom_point(data=all_sites, stroke=1, size = 5, aes(x=long, y=lat), shape=21) +
  xlab(expression(paste("Longitude (", degree,"E)"))) + 
  ylab(expression(paste("Latitude (", degree,"N)"))) +
  scalebar(x.min=81, x.max=90, y.min=38.6, y.max=81, dist = 250, height = 0.007, 
           st.dist = 0.008, st.size=4, dist_unit = "km",
           transform = TRUE, model = "WGS84", location = "bottomleft") +
  geom_text(data = countries, 
            aes(x = long, y = lat, label = name),
            fontface="italic",
            size = 4) +
  geom_text(data = mountains,
            aes(x = long, y = lat, label = name, angle = rotate),
            fontface="italic",
            alpha = 0.5,
            size = 3.2) +
  geom_label_repel(data=all_sites, aes(x=long, y=lat, label=site), size=4.5, 
                   box.padding = 0.5, 
                   #nudge_x = 1,
                   #nudge_y = .6,
                   label.padding = 0.2) +
  annotation_north_arrow(location = "br", width = unit(1.2, "cm")) +
  theme(axis.text = element_text(size = 12), 
        axis.title = element_text(size = 12),
        plot.margin=grid::unit(c(0,0,0,0), "mm"))  

ggsave("plots/Figure_1.pdf", figure1, 
       width = 21,
       height = 13,
       units = "cm",
       scale = 1.1, device = grDevices::cairo_pdf)

#print(figure1)
