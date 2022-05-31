# assess sex heterogeneity in sighting prob off CAm/SMex

# define functions
source("sim.CRclosed.r")

# assess sex het
## filter and edit sight for CAm/SMex 
sight.camsm.all <- sight %>%
  # drop Dead animals
  filter(!Dead) %>%
  # filter Central America and S Mexico sightings
  filter(loc %in% c(10:17,19)) %>%  
  # create occ to shift Nov/Dec into next year
  mutate(occ = if_else(mm %in% 11:12, year+1, as.numeric(year))) %>%
  # filter out is.na(occ) (ids for records without occ all have other CASM records)
  filter(!is.na(occ)) %>% 
  # filter out Q C, D, 3, 4, PQ
  filter(is.na(Quality) | !(Quality %in% c("C","D","PQ","3","4"))) %>%
  # drop record with incorrect lat (other records for same id-date correct)
  filter(is.na(Position.Type) | is.na(Locality) | !(Position.Type=="GENERAL" & Locality=="Padre Ramos" & year %in% 2016:2017)) %>%
  # drop single record with wrong loc
  filter(!(id==11141 & year==2005 & loc==11)) %>% 
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
  select(-c(postype, nidd, maxdd))

## edit and merge genetic sex data fom sample table
### any disagreements among SWFSC sex assignments of same animal?
misid <- biop %>% filter(!is.na(id) & SWFSC.Sex %in% c("F","M")) %>% 
  group_by(id) %>% mutate(n=n()) %>% filter(n>1) %>%
  select(id, SWFSC.Sex) %>% distinct() %>% group_by(id) %>% mutate(n=n()) %>% filter(n>1) %>% 
  arrange(id)
### non-CAm biopsies from WC and Mex
id.sex.caor <- biop %>% 
  filter(!is.na(id) & ((OSU.Sex %in% c("M","F")) | (SWFSC.Sex %in% c("M","F")))) %>% 
  filter(as.numeric(loc)>=30 & as.numeric(loc)<=73) %>%    # do not gain biopsies if use loc>=20
  filter(Certainty.of.Sample.ID %in% c("Positive","Probable") | is.na(Certainty.of.Sample.ID)) %>% 
  filter(Sample.Number != "JKJ-20181112-1") %>% 
  # where sex assignments disagree, choose M
  group_by(id) %>% summarize(sex = case_when(any(SWFSC.Sex %in% c("F","M")) ~ SWFSC.Sex[SWFSC.Sex %in% c("F","M")][1],
                                             any(OSU.Sex=="M") ~ "M",
                                             TRUE ~ "F"),
                             sexorg = case_when(all(SWFSC.Sex %in% c("F","M")) & all(OSU.Sex %in% c("F","M")) & all(SWFSC.GENDETER==10) ~ "osu",
                                                any(SWFSC.Sex %in% c("F","M")) & any(OSU.Sex %in% c("F","M")) ~ "sw.osu",
                                                any(SWFSC.Sex %in% c("F","M")) ~ "sw",
                                                any(OSU.Sex %in% c("F","M")) ~ "osu"))
sight.camsm.sex.caor <- left_join(sight.camsm.all, id.sex.caor) %>% filter(!is.na(sex))
## (1) estimate heterogeneity from capture histories
sight.camsm.sex.caor %>% group_by(id,sex) %>% summarize(nocc=n_distinct(occ)) %>% 
  group_by(sex) %>% summarize(n=n(), ns=mean(nocc), ns.se=sd(nocc)/sqrt(n), nalt=sum(nocc>1))
# including Certainty.of.Sample.ID=NA
# sex       n    ns   ns.se   nalt
# F        34  1.24   0.104   6
# M        68  1.81   0.133   32    #### high M:F ratio IN 2021
# approximate ratio as lognormal, with median=0.81/0.24=3.38, analytic cv=sqrt((0.104/0.24)^2+(0.133/0.81)^2)=0.463
# then sigma=sqrt(log(.463^2+1))=0.441, mu=log(3.38)=1.22, mean=exp(1.22+0.5*0.441^2)=3.73
# qlnorm(p=c(0.16, 0.84),meanlog=1.22,sdlog=0.441)  # 2.2, 5.3


# estimate potential bias in population estimate ignoring sex heterogeneity
# must depend on p, K (diff in overall p among classes)
# (RMark Closed Mt performed same as Mt in Rcapture, RMark Closedhet performed abysmally)
library(Rcapture)
## simulate 
N = 1100
nsim = 5000
n.chaoneg = rep(NA, nsim)
n.mt = rep(NA, nsim)
aic.mt = rep(NA, nsim)
bic.mt = rep(NA, nsim)
n.chao = rep(NA, nsim)
aic.chao = rep(NA, nsim)
bic.chao = rep(NA, nsim)
n.pois = rep(NA, nsim)
n.darr = rep(NA, nsim)
for (i in 1:nsim) {
  test <- sim.closed(N=N, K=3, p=c(0.06,.06,.23), hetfac = 1/3.4, array2d=TRUE)
  mrneg <- closedp(test$Y, neg=FALSE)
  mr <- closedp(test$Y)
  n.chaoneg[i]= mrneg$results[7,1]
  n.mt[i] = mr$results[2,1]
  aic.mt[i] = mr$results[2,5]
  bic.mt[i] = mr$results[2,6]
  n.chao[i]= mr$results[7,1]
  aic.chao[i] = mr$results[7,5]
  bic.chao[i] = mr$results[7,6]
  n.pois[i]= mr$results[8,1]
  n.darr[i]= mr$results[9,1]
  rm(test)
}
#n.aic <- if_else(aic.mt<aic.chao, n.mt, n.chao)
#n.bic <- if_else(bic.mt<bic.chao, n.mt, n.chao)
mean(n.mt)/N   # 0.917(-8.3%), 0.810(-19.0%), 0.719(-28.1%),
#   for hetfac =     1/2.2,         1/3.4,         1/5.3,    respectively; 
# since bias estimated for +/-SE of estimated heterogeneity (2.2-5.3) is approximately symmetrical about the bias estimated for the median (3.4),
# uncertainty in bias from heterogeneity can be parameterized as -19.0(+/-9.9), so bias cf estim as 1/(1-.105+.052-rnorm(-.190,.099))
mean(n.chao)/N  
mean(n.chaoneg)/N
mean(n.pois)/N
mean(n.darr)/N
median(n.mt)/N
median(n.chao)/N
median(n.aic)/N
median(n.bic)/N
median(n.chaoneg)/N
median(n.pois)/N
median(n.darr)/N

# collated mean and median results in Nclosed.sims.xlsx
# changing N, p made little difference
