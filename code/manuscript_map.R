library(sf)
library(tidyverse)
library(ggspatial)
library(patchwork)

first_bbox <- st_bbox(
  c(xmin = -126, xmax = -110,
    ymin = 35, ymax = 53),
  crs = 4326)

osoyoos_lake<- read_sf("D:/Margot/3 STOCKWELL_GIS/OKANAGAN/SHAPEFILES_NEW/US_Osoyoos.shp") |> 
  st_transform(3005)
border <- read_sf("D:/Margot/3 STOCKWELL_GIS/OKANAGAN/SHAPEFILES_NEW/Border.shp") |> 
  st_transform(3005)
OK_lakes <- read_sf("D:/Margot/3 STOCKWELL_GIS/OKANAGAN/SHAPEFILES_NEW/Okan_Lakes_Basic.shp") |> 
  st_transform(3005)
transects <- read_sf("D:/Margot/3 STOCKWELL_GIS/OKANAGAN/SHAPEFILES OKN/VDS_McIntyre.shp") |> 
  st_transform(3005)
streams <- read_sf("D:/Margot/3 STOCKWELL_GIS/OKANAGAN/SHAPEFILES OKN/Spwn_Grnds_Strms.shp") |> 
  st_transform(3005)

rivers  <- st_read("../../data/spatial_data/HydroRIVERS_v10_na.gdb", quiet = FALSE) 
columbia_all <- rivers |> filter(MAIN_RIV == 70386275) |> 
  st_transform(st_crs(osoyoos_lake))

Okanagan_dam <- st_sf(name = "Okanagan Dam",   geometry = st_sfc(st_point(c(-119.5891, 49.5286)),crs = 4326)) |> 
  st_transform(st_crs(VDS13))


i <- st_nearest_feature(Okanagan_dam, columbia_all)

start_id <- columbia_all$HYRIV_ID[i]

#get downstream reaches from Okanagan####
ids <- start_id
current <- start_id
while(TRUE) {

  next_id <- columbia_all$NEXT_DOWN[
    match(current, columbia_all$HYRIV_ID)
  ]
  
  if(length(next_id) == 0 || is.na(next_id) || next_id == 0)
    break
  
  ids <- c(ids, next_id)
  current <- next_id
}

okanagan_downstream <- columbia_all[columbia_all$HYRIV_ID %in% ids, ]

ggplot(okanagan_downstream)+
  geom_sf()

columbia_bbox <- st_bbox(st_union(OK_lakes, okanagan_downstream)) + c(-200000, -10000, 10000, -60000)

OK_lakes <- OK_lakes |> 
  st_crop(columbia_bbox)

countries <- ne_countries(country = c("Canada", "United States of America"), scale = "large") |>
  st_transform(st_crs(columbia_bbox)) |> 
  st_crop(columbia_bbox)

columbia_crop <- columbia_all |> st_crop(columbia_bbox)

ggplot(columbia_all |> filter(ORD_FLOW<6), aes(color = ORD_FLOW))+
  geom_sf()+
  scale_color_viridis_c()

lake_labels <- data.frame(
  label = c("Okanagan L.",
            "Skaha L.",
            "Vaseux L.",
            "Osoyoos L."),
  x = c(1472500, 1473000, 1477500, 1483000)-2000,
  y = c(523000, 510000, 496000, 473000)
)


index_limits <- subset(transects, NAME %in% c("McIntyre Dam", "VDS12"))
y_top <- 492708.5
y_bottom <- 484498.2

VDS13 <- transects |> filter(NAME == "VDS13") |> st_line_sample(sample = 0.5)

ECCC_station <- st_sf(name = "Station 08NM085",   geometry = st_sfc(st_point(c(-119.566389, 49.114444)),crs = 4326)) |> 
  st_transform(st_crs(VDS13))

can_us_label <- st_as_sf(data.frame(lon = -119.7, lat = c(49.015, 48.99), label = c("CANADA", "USA")),
  coords = c("lon", "lat"),
  crs = 4326) |> 
  st_transform(st_crs(VDS13))

town_labels <- data.frame(
  name = c("Malott", "Oroville", "Omak", "Oliver", "Okanagan Falls", "Penticton", "Brewster", "Nighthawk", "Princeton", "Osooyos"),
  lon  = c(-119.73, -119.43, -119.53, -119.55, -119.58, -119.60, -119.78, -119.62,  -120.50, -119.47),
  lat  = c(  48.28,   48.94,   48.41,   49.18,   49.34,   49.50,   48.10,   48.99,    49.48, 49.03),
  nudge_x = c(1000,   1000,   1000,    1000,    1000,    1000,    1000,    1000,     1000, -1000),
  nudge_y = c(   0,       0,       0,       -500,       0,       0,    0.03,    0.04,       0, 0),
  hjust   = c(   1,       1,       1,       0,       0,       0,       1,     0.5,        0, 1)) |> 
  filter(name != "Nighthawk")

towns_sf <- st_as_sf(town_labels, coords = c("lon", "lat"), crs = 4326) |> 
  st_transform(st_crs(can_us_label))

map1bb <- st_bbox(c( xmin = 1455000, xmax = 1486000, ymin = 460000, ymax = 525000), crs = st_crs(can_us_label))

map1bb_object <- st_as_sfc(map1bb, crs = st_crs(can_us_label))


study_area <- ggplot() +
  geom_sf(data = streams, colour = "grey90", linewidth = 0.2) +
  geom_sf(data = columbia_crop |> filter(ORD_FLOW %in% c(6)), color = "grey75")+
  geom_sf(data = columbia_crop |> filter(ORD_FLOW %in% c(5)) |> filter(!(HYRIV_ID %in% okanagan_downstream$HYRIV_ID)), color = "grey75")+
  geom_sf(data = OK_lakes, fill = "grey85", linewidth = 0.3) +
  geom_sf(data = osoyoos_lake, fill = "grey85", linewidth = 0.3) +
  geom_sf(data = border, color = "grey50", linewidth = 0.3)+
  geom_text(data = lake_labels, aes(x, y, label = label), fontface = "italic", size = 3)+
  geom_sf(data = VDS13, pch = 19, size = 2)+
  geom_sf(data = ECCC_station, pch = 19, size = 2)+
  geom_sf(data = towns_sf, pch = 19, size = 2)+
  geom_sf_text(data = towns_sf, aes(label = name, hjust = hjust), size = 2.2, nudge_x = towns_sf$nudge_x, nudge_y = towns_sf$nudge_y) +
  annotate("segment", x = 1469500,  xend = 1480000, y = y_top,  yend = y_top, linetype = 2,   linewidth = 0.3) +
  annotate("segment", x = 1469500,  xend = 1480000, y = y_bottom, yend = y_bottom, linetype = 2, linewidth = 0.3)+
  annotate("text", x = 1476000, y = 489000, label = "Index\nspawning\nsection",size = 3)+
  annotate("text",  x = 1468000, y = 478000, label = "Okanagan River", angle = -75, fontface = "italic",  size = 3)+
  annotate("text", label = "VDS 13", x = st_coordinates(VDS13)[1]-3500, y = st_coordinates(VDS13)[2], size = 2.5)+
  annotate("text", label = "Station\n08NM085", x = st_coordinates(ECCC_station)[1]+3500, y = st_coordinates(ECCC_station)[2]+2000, size = 2.5)+
  geom_sf_text(data = can_us_label, aes(label = label), size = 2.4, color = "grey30", hjust = 0) +
  
  scale_x_continuous(breaks = seq(-119.8, -119.2, by = 0.2))+
  scale_y_continuous(breaks = seq(49, 49.5, by = 0.1))+
  coord_sf(xlim = map1bb[c("xmin", "xmax")], ylim = map1bb[c("ymin", "ymax")],, expand = FALSE)+
  
  annotation_north_arrow(location = "tl", which_north = "true", style = north_arrow_minimal, pad_x = unit(0.1, "cm"), pad_y = unit(0.25, "cm"), height = unit(1, "cm"), width = unit(1, "cm"),)+
  annotation_scale(location = "tr",  width_hint = 0.18, pad_y = unit(1, "cm")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title = element_blank())

study_area

#larger Columbia map####

can_us_label2 <- st_as_sf(data.frame(lon = -122, lat = c(49.15, 48.85), label = c("CANADA", "USA")),
                         coords = c("lon", "lat"),
                         crs = 4326) |> 
  st_transform(st_crs(VDS13))

Malott <- towns_sf |> filter(name == "Malott")
Lake_pateros <- st_as_sf(data.frame(lon = -119.743314, lat = 48.092655, label = "L. Pateros"),
                         coords = c("lon", "lat"),
                         crs = 4326) |> 
  st_transform(st_crs(VDS13))

regional_map <- ggplot(okanagan_downstream)+
  geom_sf(data = countries, fill = "grey95")+
  geom_sf(data = columbia_crop |> filter(ORD_FLOW<=5), color = "grey75")+
  geom_sf(data = columbia_crop |> filter(ORD_FLOW<=3), aes(linewidth = DIS_AV_CMS))+
  geom_sf()+
  scale_linewidth(range = c(0.1, 0.8), guide = "none") +
  geom_sf(data = OK_lakes)+
  geom_sf(data = map1bb_object,fill = NA, linetype = 2, size = 0.8)+
  geom_sf(data = Malott)+
  geom_sf_text(data = Malott, aes(label = name), hjust = 1.15, size = 3)+
  geom_sf(data = Lake_pateros)+
  geom_sf_text(data = Lake_pateros, aes(label = label), hjust = 1.15, size = 3)+
  coord_sf(expand = FALSE)+
  theme_bw()+
  geom_sf_text(data = can_us_label2, aes(label = label), size = 3.5, color = "grey30", hjust = 0)+
  annotate("text",  x = 1495000, y = 415000, label = "Okanogan R.", angle = -95, fontface = "italic",  size = 3)+
  annotate("text",  x = 1350000, y = 110000, label = "Columbia R.", angle = 4.5, fontface = "italic",  size = 4)+
  annotate("text", x = 1420000, y = 520000, label = "Panel B",size = 3.5)+
  labs(x = NULL, y = NULL)

#combine maps####
regional_map + study_area+
  plot_annotation(tag_levels = "A")
ggsave("./figures/manuscript_map.pdf", width = 9, height = 6)

