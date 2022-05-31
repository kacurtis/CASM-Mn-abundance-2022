# spatial randomization test

# load libraries
library(ggplot2)

# compare all distances among within-ind interannual sightings to all distances over all indiv
## use randomization of lat/trapid within year after sampling single point per ind-year-loc, prior to selection of my animals
### remove missing locations, years
sightloc <- sight.camsm %>% 
  # filter(!duplicated(paste(id, trapid, occ))) %>%   # keep only first sighting per occ-trapid of each individual? results are the same. 
  select(id, occ, lat, trapid, traplat) %>% 
  group_by(id) %>% mutate(ny=n_distinct(occ)) %>% 
  ungroup() %>% 
  mutate(lat = if_else(is.na(lat), traplat, lat))

# define function for distance metric
dist.ind.mean <- function(x) {
  # filter multiyear individuals
  my.iyl <- x %>% filter(ny>1)                                    
  # calculate mean distance among within-individual sightings 
  d.ind <- my.iyl %>% 
    group_by(id) %>% tidyr::nest() %>% 
    # mutate(dists = purrr::map(data, ~ geosphere::distm(data.frame(lon=.x$lon,lat=.x$lat))/1000),
    #        xdist = purrr::map(dists, function(x) mean(x[upper.tri(x)])))   # mean per individual
    # mutate(xdist = purrr::map(data, ~mean(dist(as.numeric(.x$trapid)))))   # mean per indiv
    mutate(xdist = purrr::map(data, ~mean(dist(as.numeric(.x$lat)))))   # mean per indiv
  return(mean(unlist(d.ind$xdist)))   
}

# permute and calculate metric
nrep <- 1000
truth <- as.numeric(rep(NA, nrep))
test <- as.numeric(rep(NA, nrep))
set.seed(1020)
for (i in 1:nrep) {
  all.iyl <- sightloc %>% 
    group_by(id, occ) %>% sample_n(1) %>% ungroup() %>% arrange(occ, id) 
  truth[i] <- dist.ind.mean(all.iyl)
  perm <- all.iyl %>% 
    group_by(occ) %>% 
    mutate(r=sample(1:n()), id=id[r], ny=ny[r]) %>%
    ungroup() %>% select(-r)
  test[i] <- dist.ind.mean(perm)
}
rm(all.iyl, perm, i)

# plot results
## lat
df <- data.frame(d = c(truth, test), Type = factor(c(rep("true", nrep), rep("permuted", nrep))))
p1 <- 
  ggplot(df, aes(d)) +
  geom_histogram(aes(fill = Type), color = "black") +
  scale_fill_manual(values = c("true" = "gray50", "permuted" = "gray80")) +
  xlab("Degrees Latitude") + ylab("Frequency") +
  scale_x_continuous(limits = c(2,4), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,850), expand = c(0,0)) +
  theme_bw() + theme(panel.grid.major= element_line(colour="white"), panel.grid.minor = element_line(colour="white"),
                     text = element_text(size = 10),
                     legend.position = c(0.95, 0.95), legend.justification = c(1, 1)) +
  theme(legend.key.size = unit(0.4,"cm"))

ggsave(plot=p1, filename="Fig3_PermutationTest.tiff", 
       width=100, height=80, units="mm", dpi = 1200, compression = "lzw")

sum(truth<test)/nrep  # 1.000

# ## km - results look the same
# hist(truth, col=2, xlim=c(600,1000), breaks=10*(65:100),
#      main="Mean Within-Individual Interannual Distance", xlab="km")
# hist(test, breaks=10*(65:120), add = TRUE)
# sum(truth<test)  # 1000

# ## trapid
# hist(truth, col=2, xlim=c(1,3), breaks=seq(1,3,0.1),
#      main="Mean within-individual interannual distance", xlab="geographic sampling units")
# hist(test, breaks=seq(1,3,0.1), add=T)
# legend("topright", col=2:1, pch=c(15,22), legend=c("true","permuted"))
