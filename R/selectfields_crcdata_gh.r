# load libraries
library(magrittr)
library(lubridate)

# unlist data
sight <- mnid$sight
loccode <- mnid$loccode
biop <- mnid$biop

# edit field types, rename, select
sight %<>% rename(vessel = Vessel, lon = Dec.Long, lat = Dec.Lat, 
                  # note that there are some records with is.na(date) but with a year
                  loc = Loc.Code, year = Year, date = Date, id = CRC.ID) %>%
  mutate(loc = as.character(loc), id = as.character(id),
         mm = month(date), mmlab = month(date, label = TRUE)) %>% 
  select(-c(AKA, Region, Sub.area, Lat, Long, Photographer:Frame, Gen.Sex, Associated.IDs, 
            Other.ID:InjuryDescription, Temp.ID, HW_enc_id:HW_enc_url, Nickname:HW_media_lic, 
            HW_enc_created_on:HW_enc_user_groups, HW_enc_public:HW_ind_alternate_id)) %>% 
  arrange(date)
loccode %<>% rename(loc = Loc.Code, rgn = Region, subarea = Sub.area) %>% 
  mutate(loc = as.character(loc))
biop %<>% rename(id = CRC.ID, loc = Loc.Code, lat = Dec.Lat, lon = Dec.Long) %>% 
  mutate(id = as.character(id), loc = as.character(loc)) %>% 
  select(-c(Preservation:Sampled.by, Sample.Photo, SourceDB, DateTimeAdded)) %>% 
  arrange(Date)

# merge regions from loccode to sight based on loc
sight %<>% left_join(loccode %>% select(loc, rgn, subarea)) %>%
  select(1:6, rgn, subarea, everything())
