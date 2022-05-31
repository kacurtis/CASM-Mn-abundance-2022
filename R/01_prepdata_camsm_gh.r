############ CREATE CENTRAL AMERICA / SOUTHERN MEXICO DATA SETS ############

# import libraries
# devtools::install_github("kacurtis/CRCutils")
library(CRCutils)
library(dplyr)

# define utility functions
source("dim.drop0.r")

# import data
mnid <- import_mnid(datdir = "data/", mnidfile = "MN ID 20210831.mdb")

# initial data filtering, editing, and selection
source("selectfields_crcdata.r")

# filter and edit sight for CAm/SMex 
sight.camsm <- sight %>%
  # filter out Q C, D, 3, 4, PQ
  filter(is.na(Quality) | !(Quality %in% c("C","D","PQ","3","4"))) %>%
  # drop Dead animals
  filter(!Dead) %>%
  # filter Central America and S Mexico sightings
  filter(loc %in% c(10:17,19)) %>%  
  # create occ to shift Nov/Dec into next year
  mutate(occ = if_else(mm %in% 11:12, year+1, as.numeric(year))) %>%
  # filter out is.na(occ) (ids for records without occ all have other CASM records)
  filter(!is.na(occ)) %>% 
  # drop record with incorrect lat (other records for same id-date correct)
  filter(is.na(Position.Type) | is.na(Locality) | !(Position.Type=="GENERAL" & Locality=="Padre Ramos" & year %in% 2016:2017)) %>%
  # drop records with wrong loc
  filter(!(id==11141 & year==2005 & loc==11) & !(id==10750 & year==2015 & loc==17)) %>% 
  # group by id and date for filtering and create variables with which to identify problems
  group_by(id, date) %>% 
  mutate(nidd=n()) %>% 
  # drop records with NA positions on days with known-position record
  filter(nidd==1 | all(is.na(lat)) | !is.na(lat)) %>%
  # recalculate records per day, calculate max lat diff per day
  mutate(nidd=n()) %>% 
  # drop incorrect duplicated WGRP records 
  filter(is.na(Research.Group) | is.na(lat) | !(Research.Group=="WGRP" & (lat<17.35 | lat>17.8) & nidd>1)) %>% 
  # recalculate records per day, calculate max lat difference per day
  mutate(nidd=n(), maxdd = if(sum(!is.na(lat))>1) max(dist(lat),na.rm=TRUE) else NA_real_) %>% 
  # keep only "PRECISE" and "Precise" when "GENERAL" or "General" in same day (20210831 update added "Estimated")
  mutate(postype = any(Position.Type %in% c("GENERAL","General","Estimated")) & any(Position.Type %in% c("PRECISE","Precise"))) %>%
  filter((postype & Position.Type %in% c("PRECISE","Precise")) | !postype) %>% 
  ungroup() %>% 
  # assign trapids from explore.cam.r
  mutate(trapid = case_when(is.na(lat) ~ NA_character_,
                            lat <= 8.25 ~ "1",
                            lat>8.25 & lat<=10.25 ~ "2",
                            lat>10.25 & lat<=11.05 ~ "3",
                            lat>11.05 & lat<=12 ~ "4",
                            lat>12 & lat<=13.25 ~ "5",
                            lat>13.25 & lat<=14.5 ~ "6",
                            lat>14.5 & lat<=16.5 ~ "7",
                            lat>16.5 ~ "8")) %>% 
  # assign is.na(lat) records to trapids
  mutate(trapid = case_when(is.na(lat) & Locality=="Zih. Area" ~ "8",
                            is.na(lat) & loc==11 ~ "2",
                            is.na(lat) & loc==16 ~ "6",
                            TRUE ~ trapid)) %>% 
  # collapse sightings to daily to eliminate remaining duplicates and avoid overweighting particular geographic coordinates
  group_by(id, occ, loc, trapid, date) %>%     # add Research.Group here and split by lat=13.75 to split Guat/El Sal
  summarize(lat = if_else(any(!is.na(lat)), mean(lat, na.rm=T), NA_real_)) %>% 
  ungroup() %>% 
  # calculate trap coordinates
  group_by(trapid) %>% mutate(traplat = mean(lat, na.rm=T)) %>% ungroup()
  
# subset without SMex
sight.cam <- sight.camsm %>% filter(loc %in% 10:17) 

# create non-spatial annualized histories
ch.camsm <- array2binom(unclass(with(sight.camsm, table(id, occ))))
ch.cam <- array2binom(unclass(with(sight.cam, table(id, occ))))
## drop first three occasions
ch.camsm.96 <- array2binom(unclass(with(sight.camsm %>% filter(occ >= 1996), table(id, occ))))
ch.cam.96 <- array2binom(unclass(with(sight.cam %>% filter(occ >= 1996), table(id, occ))))

# create spatial annualized capture histories
## trap information
traplocs <- sort(unique(sight.camsm$traplat))
## spatial capture histories
y3d.camsm <- array2binom(unclass(with(sight.camsm, table(id, trapid, occ))))
y3d.cam <- array2binom(unclass(with(sight.cam, table(id, trapid, occ))))

# closed population data structures
## trap information
yy <- 2018:2021
traps.1821 <- sight.camsm %>% filter(occ %in% yy) %>% select(trapid, traplat) %>% distinct() %>% arrange(trapid)
traplocs.1821 <- array(traps.1821$traplat, dim=length(traps.1821$traplat), dimnames=list(traps.1821$trapid))
traps.cam.1821 <- sight.cam %>% filter(occ %in% yy) %>% select(trapid, traplat) %>% distinct() %>% arrange(trapid)
traplocs.cam.1821 <- array(traps.cam.1821$traplat, dim=length(traps.cam.1821$traplat), dimnames=list(traps.cam.1821$trapid))
yy <- 2019:2021
traps.1921 <- sight.camsm %>% filter(occ %in% yy) %>% select(trapid, traplat) %>% distinct() %>% arrange(trapid)
traplocs.1921 <- array(traps.1921$traplat, dim=length(traps.1921$traplat), dimnames=list(traps.1921$trapid))
traps.cam.1921 <- sight.cam %>% filter(occ %in% yy) %>% select(trapid, traplat) %>% distinct() %>% arrange(trapid)
traplocs.cam.1921 <- array(traps.cam.1921$traplat, dim=length(traps.cam.1921$traplat), dimnames=list(traps.cam.1921$trapid))
yy <- 2004:2006
traps.cam.0406 <- sight.cam %>% filter(occ %in% yy) %>% select(trapid, traplat) %>% distinct() %>% arrange(trapid)
traplocs.cam.0406 <- array(traps.cam.0406$traplat, dim=length(traps.cam.0406$traplat), dimnames=list(traps.cam.0406$trapid))
rm(yy)

## spatial capture history
### CAm with SMex
y3d.camsm.1821 <- dim.drop0(dim.drop0(y3d.camsm[,,dimnames(y3d.camsm)$occ %in% 2018:2021], 1), 2)
y3d.camsm.1921 <- dim.drop0(dim.drop0(y3d.camsm[,,dimnames(y3d.camsm)$occ %in% 2019:2021], 1), 2)
y3d.cam.1821 <- dim.drop0(dim.drop0(y3d.cam[,,dimnames(y3d.cam)$occ %in% 2018:2021], 1), 2)
y3d.cam.1921 <- dim.drop0(dim.drop0(y3d.cam[,,dimnames(y3d.cam)$occ %in% 2019:2021], 1), 2)
y3d.cam.0406 <- dim.drop0(dim.drop0(y3d.cam[,,dimnames(y3d.cam)$occ %in% 2004:2006], 1), 2)
