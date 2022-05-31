# library(ggplot2)
# library(sf)
# library(sp)
# library(mapdata)
# library(maps)
# library(maptools)
# library(rnaturalearth)
# library(rnaturalearthdata)
# library(rnaturalearthhires)
# library(ggthemes)
# library(scales)
# 
# lon.range <- c((-130), (-78))
# lat.range <- c(7, 55)
# 
# world <- ne_countries(scale = "medium", returnclass = "sf")
# world.sp <- ne_countries(returnclass = "sp")
# country_points <- data.frame(world.sp$name, getSpPPolygonsLabptSlots(world.sp))
# names(country_points) <- c("name", "X", "Y")
# country_points <- country_points[which(country_points$name %in% c("Canada", "Costa Rica", "Guatemala", "Honduras",
#                                                                   "Nicaragua",  "Mexico", "Panama", "El Salvador", "United States")),]
# country_points$name <- c("Canada", "Costa Rica", "Guat.", "Hond.", "Mexico", "Nic.", "Panama", "El Sal.", "U.S.A.")
# points <- list()
# points[[1]] <- country_points[c(3,4,5,6),]
# points[[2]] <- country_points[c(8),]
# points[[3]] <- country_points[c(7),]
# points[[4]] <- country_points[c(2),]
# points[[5]] <- country_points[c(1),]
# points[[8]] <- country_points[c(9),]
# points[[7]] <- data.frame(name = "Islas Revillagigedo", X= -123, Y = 22)
# class(world)
# 
# states.sp <- ne_states(country = c("united states of america", "canada","mexico"), returnclass = "sp")
# states.sf <- ne_states(country = c("united states of america", "canada","mexico"), returnclass = "sf")
# states.sf <- states.sf[which(states.sf$name %in% c("Oregon", "California", "Washington", "British Columbia")),]
#                                                    # "Chiapas","Oaxaca","Guerrero","Michoacan","Colima","Jalisco","Nayarit","MichoacÃ¡n")),]
# state_points <- data.frame(states.sp$name, getSpPPolygonsLabptSlots(states.sp))
# names(state_points) <- c("name", "X", "Y")
# # Colima_point <- state_points[which(state_points$name %in% c("Colima")),]
# # Colima_point$name <- "Colima"
# state_points <- state_points[which(state_points$name %in% c("Oregon", "California", "Washington", "British Columbia")),]
#                                                             # "Chiapas","Oaxaca","Guerrero","Michoacan","Jalisco","Nayarit","MichoacÃ¡n")),]
# state_points$name <- c("Washington", "British Columbia", "California","Oregon")
# 
# coords.REV <- data.frame(x = c(-116, -112), y = c(21.5, 18.9))
# # coords.Colima <- data.frame(x = c(-106, -103.9), y= c(18.2, 19.14))
# 
# g1 <- ggplot(data = world) +
#   geom_sf() +
#   geom_sf(data = states.sf, color = "gray50") +
#   geom_sf(data = world, fill=NA) +
#   geom_text(data = points[[1]], mapping=aes(x=X, y=Y, label=name), color="black", fontface = "italic", size = 4, check_overlap=FALSE) +
#   geom_text(data = points[[2]], mapping=aes(x=X, y=Y, label=name), hjust = 1, vjust = 1, nudge_y = -0.5, color="black", fontface = "italic",  size = 4, check_overlap=TRUE) +
#   geom_text(data = points[[3]], mapping=aes(x=X, y=Y, label=name), nudge_x = 1, nudge_y = 1.5, color="black", fontface = "italic",  size = 4, check_overlap=TRUE) +
#   geom_text(data = points[[4]], mapping=aes(x=X, y=Y, label=name), nudge_x = -2, nudge_y = -1, color="black", fontface = "italic",  size = 4, check_overlap=TRUE) +
#   geom_text(data = points[[5]], mapping=aes(x=X, y=Y, label=name), nudge_x = -5, nudge_y = -2, color="black", fontface = "italic",  size = 4, check_overlap=TRUE) +
#   geom_text(data = points[[8]], mapping=aes(x=X, y=Y, label=name), nudge_x = -8, color="black", fontface = "italic",  size = 4, check_overlap=TRUE) +
# #  geom_text(data = points[[6]], mapping=aes(x=X, y=Y, label=name), hjust = 1, size = 5, fontface = "bold", check_overlap=TRUE) +
#   geom_text(data = points[[7]], mapping=aes(x=X, y=Y, label=name), nudge_x = 4, nudge_y = -0.5, fontface = "italic",  size = 3, check_overlap=TRUE) +
#   geom_text(data = state_points, mapping=aes(x=X, y=Y, label=name), color="black", fontface = "italic", size = 3, check_overlap=FALSE) +
#   # geom_text(data = Colima_point, mapping=aes(x=X, y=Y, label=name), nudge_x = -3.3, nudge_y = -1.3, color="black", fontface = "italic", size = 3, check_overlap=FALSE) +
# #  geom_text(data = state_points, mapping=aes(x=X, y=Y, label=name), color="black", fontface = "italic", size = 3, check_overlap=FALSE) +
#   geom_segment(data = coords.REV, aes(x=coords.REV$x[1], y=coords.REV$y[1], xend = coords.REV$x[2], yend = coords.REV$y[2])) +
#   # geom_segment(data = coords.Colima, aes(x=coords.Colima$x[1], y=coords.Colima$y[1], xend = coords.Colima$x[2], yend = coords.Colima$y[2])) +
#   coord_sf(xlim = lon.range, ylim = lat.range) +
#   theme(panel.background = element_rect(fill = 'white'), 
#         panel.border = element_rect(linetype = 1, fill = NA), 
#         axis.text=element_text(size=8), plot.margin = margin(0, 0, 0, 0, "cm"),
#         axis.title.x=element_blank(), axis.title.y=element_blank())
# 
# ggsave(plot=g1, filename="bigmap.tiff", 
#        width=170, height = 180, units="mm", dpi = 1200, compression = "lzw")
# 
# 
###########################################################

# Figure 1 map extended
library(ggplot2)
library(sf)
library(sp)
library(mapdata)
library(maps)
library(maptools)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(ggthemes)
library(scales)

lon.range <- c((-130), (-77))
lat.range <- c(7, 52)

world <- ne_countries(scale = "medium", returnclass = "sf")
world.sp <- ne_countries(returnclass = "sp")
country_points <- data.frame(world.sp$name, getSpPPolygonsLabptSlots(world.sp))
names(country_points) <- c("name", "X", "Y")
country_points <- country_points[which(country_points$name %in% c("Canada", "Costa Rica", "Guatemala", "Honduras",
                                                                  "Nicaragua",  "Mexico", "Panama", "El Salvador", "United States")),]
# country_points$name <- c("Canada", "Costa Rica", "Guat.", "Hond.", "Mexico", "Nic.", "Panama", "El Sal.", "U.S.A.")
country_points$name <- c("Canada", "C.R.", "Guat.", "Hon.", "Mexico", "Nic.", "Pan.", "E.S.", "U.S.A.")

states.sp <- ne_states(country = c("united states of america", "canada","mexico"), returnclass = "sp")
states.sf <- ne_states(country = c("united states of america", "canada","mexico"), returnclass = "sf")
states.sf <- states.sf[which(states.sf$name %in% c("Oregon", "California", "Washington", "British Columbia",
                                                   "Chiapas","Oaxaca","Guerrero","Michoacan","Colima","Jalisco","Nayarit","Michoacán")),]
state_points <- data.frame(states.sp$name, getSpPPolygonsLabptSlots(states.sp))
names(state_points) <- c("name", "X", "Y")
state_points <- state_points[which(state_points$name %in% c("Oregon", "California", "Washington", "British Columbia",
                                                            "Chiapas","Oaxaca","Guerrero","Michoacan","Colima","Jalisco","Nayarit","Michoacán")),]
state_points$name <- c("WA", "BC", "CA","CS","","NA","JA","MI","GR","OA","OR")

g1 <- 
  ggplot(data = world, size=0.4) +
  geom_sf() +
  geom_sf(data = states.sf, color = "gray40", size = 0.2, linetype = "solid") +    #  fontface = "italic"
  geom_sf(data = world, fill=NA) +
  geom_text(data = country_points[1,], mapping=aes(x=X, y=Y, label=name), color="black", size = 3, nudge_y = -6) +
  geom_text(data = country_points[2,], mapping=aes(x=X, y=Y, label=name), nudge_x = 0.3, color="black",  size = 2.3) +
  geom_text(data = country_points[3,], mapping=aes(x=X, y=Y, label=name), nudge_y = -0.3, color="black",  size = 2.5) +
  geom_text(data = country_points[c(4,6),], mapping=aes(x=X, y=Y, label=name), color="black",  size = 2.5) +
  geom_text(data = country_points[5,], mapping=aes(x=X, y=Y, label=name), nudge_y = 1, color="black", size = 3) +
  geom_text(data = country_points[7,], mapping=aes(x=X, y=Y, label=name), nudge_x = -1.2, nudge_y = 0.02, color="black", size = 2.3) +
  geom_text(data = country_points[8,], mapping=aes(x=X, y=Y, label=name), nudge_x = 0.4, nudge_y = -0.1, color="black", size = 2) +
  geom_text(data = country_points[9,], mapping=aes(x=X, y=Y, label=name), color="black", size = 3) +
  geom_text(data = state_points, mapping=aes(x=X, y=Y, label=name), color="gray40", size = 2.3) +
  geom_text(data = state_points[2,], mapping=aes(x=X, y=Y, label="BC"), nudge_x = 3, nudge_y = -3, color="gray40", size = 2.3) +
  geom_text(aes(x=-118, y=13, label="Pacific Ocean"), color="gray65", size = 4, fontface="italic") +
  geom_segment(aes(x=state_points[5,2], xend=-100.5, y=state_points[5,3], yend=21.15), size = 0.2, color="gray40", linetype = "solid") +
  annotate("text", x=-99.8, y=21.3, size = 2.3, color="gray40", label="CL") +
  coord_sf(xlim = lon.range, ylim = lat.range) + 
  theme(panel.background = element_rect(fill = 'white'), 
        panel.border = element_rect(linetype = 1, fill = NA), 
        axis.text=element_text(size=8), plot.margin = margin(0, 0, 0, 0, "cm"),
        axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_segment(aes(x = -79.8, y = 7.25, xend = -78.3, yend = 7.25),  
               linejoin = "mitre", size = 1, color="black",
               arrow = arrow(length = unit(0.13, "cm"), type = "closed")) +
  geom_segment(aes(x = -94.1, y = 14.5, xend = -92.6, yend = 14.5),  
               linejoin = "mitre", size = 1, color="gray60",   
               arrow = arrow(length = unit(0.13, "cm"), type = "closed")) +
  geom_segment(aes(x = -105, y = 17, xend = -102.8, yend = 18),  
               linejoin = "mitre", size = 1, color="black", 
               arrow = arrow(length = unit(0.13, "cm"), type = "closed")) +
  geom_segment(aes(x = -105, y = 17, xend = -105, yend = 19.2),  
               linejoin = "mitre", size = 1, color="black",
               arrow = arrow(length = unit(0.13, "cm"), type = "closed")) +
  geom_segment(data = data.frame(x = rep(-105,113), y = rep(17,113), xend = -105 + 1.1 * cos(seq(pi/7,pi/2,0.01)), yend = 17 + 1.1 * sin(seq(pi/7,pi/2,0.01))), 
               aes(x=x, y=y, xend=xend, yend=yend),
               linejoin = "mitre", size = 1, color="black") +
  geom_segment(aes(x=-140, xend=-78 ,y = 7.25, yend=7.25),  size = 0.4, color="black", linetype = "dashed") +
  geom_segment(aes(x=-140,xend=-105.02, y=19.2, yend=19.2), size = 0.4, color="black", linetype = "dashed") +
  geom_segment(aes(x=-140,xend=-103.92,y=18.6,yend=18.6), size = 0.4, color="gray60", linetype = "dashed") +
  geom_segment(aes(x=-140,xend=-105.82,y=20.4,yend=20.4), size = 0.4, color="gray60", linetype = "dashed")

ggsave(plot=g1, filename="Fig1_Map_ext.tiff", 
       width=165, height = 160, units="mm", dpi = 1200, compression = "lzw")
