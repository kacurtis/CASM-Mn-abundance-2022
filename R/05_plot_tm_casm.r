# plots for Tech Memo on CAm abundance

# create data sets
source("01_prepdata_camsm_gh.r")

# recaptures
table(rowSums(array2binom(apply(y3d.camsm.1921, sum, MARGIN=c(1,3))))) 
table(rowSums(array2binom(apply(y3d.camsm.1921, sum, MARGIN=c(1,2))))) 

# load libraries
require(ggplot2)
require(gridExtra)
require(grid)


# Figure 1: map
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

lon.range <- c((-108), (-77))
lat.range <- c(6.5, 21.5)

world <- ne_countries(scale = "medium", returnclass = "sf")
world.sp <- ne_countries(returnclass = "sp")
country_points <- data.frame(world.sp$name, getSpPPolygonsLabptSlots(world.sp))
names(country_points) <- c("name", "X", "Y")
country_points <- country_points[which(country_points$name %in% c("Costa Rica", "Guatemala", "Honduras",
                                                                  "Nicaragua",  "Mexico", "Panama", "El Salvador")),]
country_points$name <- c("CR", "GT", "HN", "MX", "NI", "PA", "SV")

states.sp <- ne_states(country = "mexico", returnclass = "sp")
states.sf <- ne_states(country = "mexico", returnclass = "sf")
states.sf <- states.sf[which(states.sf$name %in% c("Chiapas","Oaxaca","Guerrero","Michoacan","Colima","Jalisco","Nayarit","Michoacán")),]
state_points <- data.frame(states.sp$name, getSpPPolygonsLabptSlots(states.sp))
names(state_points) <- c("name", "X", "Y")
state_points <- state_points[which(state_points$name %in% c("Chiapas","Oaxaca","Guerrero","Michoacan","Colima","Jalisco","Nayarit","Michoacán")),]
state_points$name <- c("CS","CL","NA","JA","MI","GR","OA")

g1 <- 
  ggplot(data = world, size=0.4) +
  geom_sf() +
  geom_sf(data = states.sf, color = "gray50", size = 0.2, linetype = "solid") +
  geom_sf(data = world, fill=NA) +
  geom_text(data = country_points[1,], mapping=aes(x=X, y=Y, label=name), color="black", fontface = "italic", size = 3.3, , nudge_x = 0.2) +
  geom_text(data = country_points[2,], mapping=aes(x=X, y=Y, label=name), hjust = 1, vjust = 1, nudge_y = -0.3, color="black", fontface = "italic",  size = 3.3) +
  geom_text(data = country_points[c(3,5),], mapping=aes(x=X, y=Y, label=name), color="black", fontface = "italic",  size = 3.3) +
  geom_text(data = country_points[4,], mapping=aes(x=X, y=Y, label=name), nudge_y = -4, nudge_x = 4, color="black", fontface = "italic",  size = 3.3) +
  geom_text(data = country_points[6,], mapping=aes(x=X, y=Y, label=name), nudge_x = -1, nudge_y = -0.2, color="black", fontface = "italic",  size = 3.3) +
  geom_text(data = country_points[7,], mapping=aes(x=X, y=Y, label=name), nudge_x = -0.45, nudge_y = 0.12, color="black", fontface = "italic",  size = 3.3) +
  geom_text(data = state_points, mapping=aes(x=X, y=Y, label=name), color="gray50", fontface = "italic", size = 2.5, check_overlap=FALSE) +
  coord_sf(xlim = lon.range, ylim = lat.range) + 
  theme(panel.background = element_rect(fill = 'white'), 
        panel.border = element_rect(linetype = 1, fill = NA), 
        axis.text=element_text(size=8), plot.margin = margin(0, 0, 0, 0, "cm"),
        axis.title.x=element_blank(), axis.title.y=element_blank()) +
  geom_segment(aes(x = -79.5, y = 7.25, xend = -78.15, yend = 7.25),  
               linejoin = "mitre", size = 1, color="black",
               arrow = arrow(length = unit(0.13, "cm"), type = "closed")) +
  geom_point(aes(-78.17, 7.25),  size = 3, color="black", shape = 8, stroke=1) +
  geom_segment(aes(x = -93.8, y = 14.5, xend = -92.4, yend = 14.5),  
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
  geom_point(aes(-105.02, 19.2),  size = 3, color="black", shape = 8, stroke=1) +
  geom_point(aes(-103.92, 18.6),  size = 3, color="gray60", shape = 8, stroke=1) +
  geom_point(aes(-105.82, 20.4),  size = 3, color="gray60", shape = 8, stroke=1)

ggsave(plot=g1, filename="Fig1_Map.tiff", 
       width=170, height = 90, units="mm", dpi = 1200, compression = "lzw")


# Figure 2: spatiotemporal distribution of captures
# devtools::install_github("NicolasH2/ggbrace")
library(ggbrace)

sight.camsm.1921.ann <- 
  sight.camsm %>% 
  filter(occ %in% 2019:2021) %>% 
  group_by(id, occ, trapid, traplat) %>% 
  summarize(lat = if_else(any(!is.na(lat)), mean(lat, na.rm=TRUE), mean(traplat))) %>% ungroup()

g1 <- 
  ggplot(data=sight.camsm.1921.ann, aes(x=factor(occ), y=lat)) +
  # add horizontal lines at trap breaks
  geom_hline(yintercept = c(8.25, 10.25, 11.05, 12, 13.25, 14.5, 16.5), linetype = 2, color = "gray") +
  geom_dotplot(binaxis = "y", binwidth = 0.05, method = "histodot", stackdir = "center", dotsize = 0.45, stackratio = 2.25, origin = 7.25) +
  labs(x="Year",y="Latitude") +
  scale_y_continuous(breaks=seq(8,19,1), limits = c(7.25, 19.25), expand = c(0,0)) +
  theme_bw() + theme(panel.grid.major= element_line(colour="white"), panel.grid.minor = element_line(colour="white"),
                     text = element_text(size = 10)) +
  geom_brace(aes(x=c(3.65, 3.8), y=c(7.25,8.25), label = "PA"), inherit.data=F, rotate=90, labelsize=3, color = "gray60") + 
  geom_brace(aes(x=c(3.65, 3.8), y=c(8.25,11.08), label = "CR"), inherit.data=F, rotate=90, labelsize=3, color = "gray60") + 
  geom_brace(aes(x=c(3.65, 3.8), y=c(11.08, 13), label = "NI"), inherit.data=F, rotate=90, labelsize=3, color = "gray60") + 
  geom_brace(aes(x=c(3.65, 3.8), y=c(13, 13.75), label = "SV"), inherit.data=F, rotate=90, labelsize=3, color = "gray60") + 
  geom_brace(aes(x=c(3.65, 3.8), y=c(13.75, 14.5), label = "GT"), inherit.data=F, rotate=90, labelsize=3, color = "gray60") + 
  geom_brace(aes(x=c(3.65, 3.8), y=c(14.5, 19.25), label = "S.\nMX"), inherit.data=F, rotate=90, labelsize=3, color = "gray60") + 
  coord_cartesian(x = c(1,3), clip = "off") +
  theme(plot.margin = unit(c(0.01, 0.11, 0, 0), units="npc"))
  
ggsave(plot=g1, filename="Fig2_DistSightings1921_alt.tiff", 
       width=110, height=140, units="mm", dpi = 1200, compression = "lzw")


# Figure 3: see 02_permute_space.r


# Figure 4: conceptual diagram of population limits vs model domain boundaries
opar <- par()

ac <- runif(7,0,100)
b <- seq(0.5,2,length.out=7)
s <- 20
plim <- 80
pmem <- ac <= plim

tiff("Fig4_poplimconcept.tif", width = 2, height = 5, units = "in", compression="lzw", res=1200,
     pointsize = 10)
par(mai = c(0,0,0,0))
plot(0,0, xlim=c(0,10),ylim=c(0,105), bty="n", xlab = "", ylab = "", main = "", 
     type="n", xaxt = "n", yaxt = "n")
lines(x=c(0,10), y=c(0,0), lwd = 1.5)
lines(x=c(0,10), y=c(100,100), lwd = 1.5)
lines(x=c(0,10), y=c(0,0), lwd = 1.5, col=2, lty=2)
lines(x=c(0,10), y=c(plim, plim), lwd = 1.5, col=2)
#abline(h=c(10,35,55,75), lty="dashed", lwd = 1.5, col="gray")
segments(b, ac-s, b, ac+s, lwd = 1.5, col=1)
segments(b[pmem], ac[pmem]-s, b[pmem], ac[pmem]+s, lwd = 1.5, col=2)
points(b, ac, pch = 21, bg = 1)
points(b[pmem], ac[pmem], pch=21, bg=2)
text(10,102,"Northern MDB", pos=2)
text(10,plim + 2,"Northern PL", pos=2)
text(10,2,"Southern PL and MDB", pos=2)
text(10, c(10,33,55,70) + 2, paste0("Trap ", 1:4), pos=2)
segments(rep(10,5), c(0.75,20.75,45.75,65.75), rep(10,5), c(19.25,44.25,64.25,79.25), lwd = 3, col="gray")
dev.off()

par(opar)


# Figure 5: posterior distributions for SCR model
library(coda)
attach("nim.SCR0pjk1D.camsm.1921.QAB.19.2.rdata")
attach("Nprocessed.SCR0pjk.1921.rdata")
pdn <- data.frame(
  N = as.vector(cf * nlive.casm.19))
pds <- data.frame(
  sigma = as.matrix(out.nim.1921.QAB.19$fit1$samples)[,"sigma"])
detach(2)

p1 <- ggplot(data=pdn, aes(N)) + geom_density(fill=I("gray70"), col=I("black"), size = 0.2) +
  theme_bw() + theme(panel.grid.major= element_line(colour="white"), panel.grid.minor = element_line(colour="white"),
                     text = element_text(size = 10)) + 
  xlab(expression('Adjusted '~italic(N))) +
  scale_x_continuous(expand = c(0, 0), limits = c(500,3100), breaks=c(500,1000,1500,2000,2500,3000)) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  geom_vline(xintercept=quantile(pdn$N, 0.2), linetype=2) +
  annotate("text", label=expression(N[min]), x=1220, y=0.0005, angle=90, size = 3)

p2 <- ggplot(data=pds, aes(sigma)) + geom_density(fill=I("gray70"), col=I("black"), size = 0.2) +
  theme_bw() + theme(panel.grid.major= element_line(colour="white"), panel.grid.minor = element_line(colour="white"),
                     text = element_text(size = 10)) + 
  xlab(expression(sigma~'('*degree~'Latitude)')) + 
  scale_x_continuous(limits = c(2, 7)) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  
ggsave(plot=marrangeGrob(grobs = list(p1, p2), nrow = 2, ncol=1, 
                         left = textGrob("Probability density", gp=gpar(fontsize=10), rot = 90), top = NULL), 
       filename="Fig5_posteriors.tiff", 
       width=120, height=120, units="mm", compression = "lzw")
