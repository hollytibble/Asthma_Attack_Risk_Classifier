# pre-amble
library(dplyr)
library(readr)
library(lubridate)
library(zoo)
library(tidyverse)
library(mltools)
library(data.table)
rm(list=ls())
set.seed(12345)

######################################################################################
### Loading in pre-cleaned files from other locations
######################################################################################

###  Source the previous Dictionary and Setup files necessary to run this script
#source("1. Read Code Dictionary.R")
#source("2. Read Code Setup.R")
#source("3. Prescribing Record Setup")

### Import Read Code Setup and Prescribing Record Setup Files
load("Data/Read Code Records.RData")
load("Data/BTS Step Records.RData")
load("Data/steroid_bursts")
load("Data/neb_saba")
load("Data/reliever_use")
load("Data/controller_use")
load("Data/prescribing_nasal")
load("Data/CSA_measures.RData")
load("Data/CMA7.RData")

population<-intersect(GP_names,unique(prescribing_asthma$ID))
rm(GP_names)

######################################################################################
### Cleaning other events files 
######################################################################################

### A&E 
A_E <- read_csv("../../linked_data/AandEWhole.csv")
A_E$date <-as.Date(A_E$ArrivalDate,format="%m/%d/%Y")
names(A_E)[1]<-"ID"   
# min(A_E$date) #2007-06-01
# max(A_E$date) # 2017-09-30
A_E<- A_E %>% 
  filter((format(date,"%Y")<2017) | 
               (format(date,"%Y")==2017 & 
                  as.numeric(format(date,"%m"))<4)) %>%
  filter(!is.na(PresentingComplaintText) | 
          !is.na(Disease1Code))
nrow(A_E) 
length(unique(A_E$ID)) 
A_E$asthma_flag<-grepl("ASTHMA",A_E$PresentingComplaintText,fixed=TRUE) +
  grepl("asthma",A_E$PresentingComplaintText,fixed=TRUE) +
  grepl("J45",A_E$Disease1Code,fixed=TRUE)+
  grepl("J46",A_E$Disease1Code,fixed=TRUE)
A_E<-A_E[which(A_E$asthma_flag!=0),] 
nrow(A_E) 
length(unique(A_E$ID)) 
A_E<-A_E %>% filter(ID %in% population)
nrow(A_E) 
length(unique(A_E$ID)) 
A_E<-A_E %>%
  mutate(Event="A&E") %>%
  dplyr::select(ID, date, Event)

### Inpatient hospital admissions
inpatient<-read_csv("../../linked_data/1516-0489_SMR01.csv")
names(inpatient)[1]<-"ID" 
inpatient$date<-as.Date(paste(substr(inpatient$ADMISSION_DATE,1,4),
                              substr(inpatient$ADMISSION_DATE,5,6),
                              substr(inpatient$ADMISSION_DATE,7,8),sep="/"),format="%Y/%m/%d")
# min(inpatient$date) # 2000-01-01
# max(inpatient$date) # 2017-03-31
# count(is.na(inpatient$date)) # no missing dates
# table(format(inpatient$date,"%Y")) # no impossible years
inpatient$asthma_flag<- grepl("J45",inpatient$condition_1,fixed=TRUE)+
  grepl("J46",inpatient$condition_1,fixed=TRUE)
inpatient<-inpatient[which(inpatient$asthma_flag!=0),] 
inpatient<- inpatient[which(inpatient$ID %in% population),] # 109469
inpatient<-inpatient %>%
  mutate(Event="Inpatient") %>%
  dplyr::select(ID, date, Event)

### Death records
mortality <- read_csv("../../linked_data/1516-0489_NRS_Deaths_Stillbirths.csv")
mortality<-mortality[which(mortality$NRSDeaths==1),] 
names(mortality)[1]<-"ID" 
mortality$date<-as.Date(paste(substr(mortality$DEATH_DATE,1,4),substr(mortality$DEATH_DATE,5,6),substr(mortality$DEATH_DATE,7,8),sep="/"),format="%Y/%m/%d")
min(mortality$date) # 2000-01-01
max(mortality$date) # 2017-03-31
# sum(is.na(mortality$date)) # no missing dates
# table(format(mortality$date,"%Y")) # no impossible years
right_censor_dates<-mortality %>%
  dplyr::rename(DOD=date) %>%
  dplyr::select(ID, DOD)
mortality$asthma_flag<- grepl("J45",mortality$PRIMARY_CAUSE_OF_DEATH,fixed=TRUE)+
  grepl("J46",mortality$PRIMARY_CAUSE_OF_DEATH,fixed=TRUE)
mortality_asthma<-mortality[which(mortality$asthma_flag!=0),] 
mortality_asthma<- mortality_asthma[which(mortality_asthma$ID %in% population),]
nrow(mortality %>% filter(asthma_flag==0 & ID %in% population))
mortality_asthma<-mortality_asthma %>%
  mutate(Event="Asthma.Death") %>%
  dplyr::select(ID, date, Event)

adherence<-adherence %>% filter(ID %in% population) 
CSA_measures<-CSA_measures %>% filter(ID %in% population) 
neb_saba<-neb_saba %>% filter(ID %in% population)
controller_use<-controller_use %>% select(-PrescDate) %>% filter(ID %in% population) 
GP<-GP %>% filter(ID %in% population)
prescribing_asthma<-prescribing_asthma %>% filter(ID %in% population)
steroid<-steroid %>% filter(ID %in% population) 
prescribing_nasal<-prescribing_nasal %>% filter(ID %in% population)

reliever_use<-reliever_use %>% 
  filter(ID %in% population &
           is.finite(as.numeric(Event_value))) %>%
  mutate(Event_value=as.character(Event_value))

######################################################################################
### Censoring and cleaning workspace - including adding missing records to my shell 
######################################################################################

### Convert these events, and those from the GP records, into list format
events<-bind_rows(A_E,adherence,controller_use,
                  CSA_measures,GP, inpatient,
                  mortality_asthma, neb_saba,
                 prescribing_asthma,prescribing_nasal,
                 reliever_use, steroid) %>%
  select(-dt1_copd)
length(unique(events$ID))

# Flag first prescription of step 0 or above for censoring (plus six months)
left_censor_dates_PIS<-prescribing_asthma %>%
  group_by(ID) %>%
  mutate(left_censor_PIS = min(date)+months(6)) %>%
  select(ID, left_censor_PIS) %>%
  dplyr::slice(1)

# Flag last prescription of step 0 or above for censoring (plus six months)
right_censor_dates_PIS<-prescribing_asthma %>%
  group_by(ID) %>%
  mutate(right_censor_PIS = max(date)+months(6)) %>%
  select(ID, right_censor_PIS) %>%
  dplyr::slice(1)

###left & right censoring dates from data
censors<-left_join(as.data.frame(population) %>% 
                     mutate(ID = as.character(population)) %>% select(ID),
                   right_censor_dates_diag)
censors<-left_join(censors,left_censor_dates_PIS)
censors<-left_join(censors,right_censor_dates_PIS)
censors<-left_join(censors,right_censor_dates)
censors<-left_join(censors,GP %>% ungroup %>% distinct(ID,dt1_copd))

# COPD Reference Category
censors<-censors %>%
  mutate(temp = as.Date(dt1_copd)-left_censor_PIS,
         COPD_ref = case_when(is.na(temp) ~ "No COPD",
                              temp<0 ~ "COPD First",
                              temp<=(5*365.25) ~ "Within 5 years",
                              temp>(5*365.25) ~ "More than 5 years later")) %>%
  select(-temp)

### clean up
censors <- censors %>%
  mutate(right_censor=as.Date(pmin.int(DOD,
                                       as.Date("31-03-2017","%d-%m-%Y"), # uniform end of data
                                       na.rm=T),origin='1970-01-01'),
         soft_right_censor=as.Date(pmin.int(as.Date("31-03-2016","%d-%m-%Y"), # one year prior to uniform end of data
                                            right_censor_diag,
                                            right_censor_PIS,
                                            na.rm=T),origin='1970-01-01'))  

events<-left_join(events,censors) %>% 
  filter(date>=left_censor_PIS) %>% 
  filter(date<=right_censor) 
length(unique(events$ID))

### Work out follow-up for each person
followup_1<-events %>% # time without COPD
  group_by(ID,COPD_ref,dt1_copd) %>%
  filter(COPD_ref!="COPD First") %>%
  dplyr::slice(1) %>%
  mutate(right_censor = pmin(right_censor,dt1_copd,na.rm = T),
         followup = right_censor-left_censor_PIS+1) %>%
  filter(followup>0) %>%
  ungroup %>%
  select(ID,followup)

followup_2<-events %>% # time with COPD-asthma overlap
  distinct(ID,COPD_ref,dt1_copd,left_censor_PIS,right_censor) %>%
  filter(!is.na(dt1_copd)) %>%
  mutate(left_censor_PIS = pmax(left_censor_PIS,dt1_copd,na.rm=T),
         followup = right_censor-left_censor_PIS+1) %>%
  filter(followup>0) %>%
  ungroup %>%
  select(ID,followup,COPD_ref)

events<-events %>% select(-left_censor_PIS,-right_censor_diag,-DOD,-right_censor,-right_censor_PIS)

## clear out workspace to save memory
rm(list=setdiff(ls(),c("events","followup_1","followup_2")))

######################################################################################
### Covariate Coding
######################################################################################

# All dates of Primary care events
dates<-events %>% 
  filter(!Event %in% c("Inpatient","A&E","Asthma.Death","steroid",
                       "CMA7_2","controllers")) %>%
  distinct(ID,date,COPD_ref,dt1_copd,soft_right_censor) %>% 
  mutate(Event = "zzz.GP Visit",Event_value = NA)

### Change to wide format
events.wide<-rbind(events,dates) %>%
  distinct() %>%
  mutate(temp=1) %>%
  mutate(Event_desc = Event) %>%
  pivot_wider(names_from=Event,values_from=temp) %>%
  dplyr::select(-c(BMI,`BMI Category`,Height,Weight,`zzz.GP Visit`)) %>%
  mutate(Event_desc = ifelse(Event_desc=='BMI',
                             'X_BMI',
                             ifelse(Event_desc=='BMI Category',
                                    'X_BMI_Category',
                                    Event_desc))) %>%
  arrange(ID,date,Event_desc)
rm(dates,events)

### CCI Comorbidity 
for(var in names(events.wide)[which(substr(names(events.wide),1,4)=="Com.")]) {
  events.wide<-events.wide %>%
    mutate(temp_var = ifelse(Event_desc==var,1,NA)) %>%
    group_by(ID) %>%
    mutate(temp_var = na.locf(temp_var, na.rm = FALSE))
  events.wide[,c(var)]<-as.numeric(ifelse(is.na(events.wide$temp_var),
                                          0,
                                          events.wide$temp_var))
  events.wide$temp_var<-NULL
}

### Diagnosis Variables & Nasal Sprays
for(var in c("Eczema","GERD","Rhinitis","Anxiety_Depression",
             "Nasal.Polyps","Nasal.Spray")) {
  events.wide<-events.wide %>%
    mutate(temp_var = ifelse(Event_desc==var,date,NA)) %>%
    group_by(ID) %>%
    mutate(temp_var = as.Date(na.locf(temp_var, na.rm = FALSE)))
  events.wide[,c(var)]<-ifelse(is.na(events.wide$temp_var),
                               "Never",
                               ifelse(is.na(events.wide$date -years(1)), # leap years mess this up!
                                      ifelse(events.wide$temp_var>=(events.wide$date +days(1)-years(1)),
                                             "In_Last_Year",
                                             ifelse(events.wide$temp_var>=(events.wide$date +days(1)-years(5)),
                                                    "In_Last_5_Years",
                                                    "Longer_than_5_Years_Ago")),
                                      ifelse(events.wide$temp_var>=(events.wide$date -years(1)),
                                      "In_Last_Year",
                                      ifelse(events.wide$temp_var>=(events.wide$date -years(5)),
                                             "In_Last_5_Years",
                                             "Longer_than_5_Years_Ago"))))
  events.wide$temp_var<-NULL
}
rm(var)


### Smoking Status
events.wide<-events.wide %>%
    mutate(temp_var = ifelse(Event_desc=="Smoking",Event_value,NA)) %>%
    group_by(ID) %>%
    mutate(temp_var = na.locf(temp_var, na.rm = FALSE),
           Smoking = ifelse(is.na(temp_var),"Never",temp_var)) %>%
  select(-temp_var)


### BMI / Height + Weight
events.wide<-events.wide %>%
  # Raw BMI
  mutate(temp_var = ifelse(Event_desc=="X_BMI",Event_value,NA),
         Obese = ifelse(!is.na(temp_var),
                        ifelse(temp_var>=30,
                               "Obese",
                               "Not"),
                        NA)) %>%
  # BMI Category
  mutate(temp_var=ifelse(Event_desc=="X_BMI_Category",Event_value,NA),
         Obese = ifelse(!is.na(temp_var),
                        ifelse(temp_var=="Obese",
                               "Obese",
                               "Not"),
                        Obese)) %>%
  # Height and Weight
  group_by(ID) %>%
  mutate(temp_var = na.locf(ifelse(Event_desc=="Height",Event_value,NA), na.rm = FALSE),
         temp_var2 = na.locf(ifelse(Event_desc=="Weight",Event_value,NA), na.rm = FALSE),
         Obese = ifelse(is.na(Obese) & !is.na(temp_var) & !is.na(temp_var2),
                        ifelse((as.numeric(temp_var2)^2/as.numeric(temp_var))>=30,
                               "Obese",
                               "Not"),
                        Obese)) %>%
  select(-temp_var,-temp_var2) %>%
  mutate(Obese = ifelse(is.na(Obese),"Not", Obese))


### Peak Flow
events.wide<-events.wide %>%
  mutate(temp_var = as.numeric(ifelse(Event_desc=="Peak Flow",Event_value,NA))) %>%
  group_by(ID) %>%
  mutate(temp_var = na.locf(temp_var, na.rm = FALSE),
         temp_varp = ifelse(is.na(temp_var),0,temp_var), #cummax thinks NA is higher than 0
         maxPF = cummax(temp_varp), #highest PF recorded so far
         last_PF = as.numeric(date-na.locf(ifelse(Event_desc=="Peak Flow",
                                                  date,NA),na.rm=F)), # time since last peak flow
         temp_var2 = ifelse(last_PF<=7 & maxPF!=0 & temp_varp!=0,
                            temp_varp/maxPF,
                            999),
         Peak_Flow = ifelse(temp_var2==999,
                            "missing",
                            ifelse(temp_var2>0.9,
                                   "gt_0.9",
                                   ifelse(temp_var2>0.8,
                                          "0.8_0.9",
                                          ifelse(temp_var2>0.7,
                                                 "0.7_0.8",
                                                 "lt_0.7")))))  %>%
    dplyr:: select(-c(temp_var,temp_varp,maxPF,last_PF,temp_var2,`Peak Flow`))


### Blood Eosinophil Counts
events.wide<-events.wide %>%
  mutate(temp_var = as.numeric(ifelse(Event_desc=="Blood.Eosinophil.Counts",Event_value,NA))) %>%
  group_by(ID) %>%
  mutate(temp_var = na.locf(temp_var, na.rm = FALSE)) %>%
  mutate(Blood_Eosinophil_Counts=ifelse(is.na(temp_var),"Missing",
                                          ifelse(as.numeric(temp_var)<0.4,
                                                 "le_0.4",
                                                 "ge_0.4"))) %>%
  dplyr::select(-temp_var,-`Blood.Eosinophil.Counts`)

### Reliever Use
median_reliever_use<-median(as.numeric(unlist(events.wide %>% ungroup %>% 
                                            filter(Event_desc=="reliever_use") %>% 
                                            select(Event_value))))
events.wide<-events.wide %>%
  mutate(temp_var = as.numeric(ifelse(Event_desc=="reliever_use",Event_value,NA)),
         temp_var2 = as.Date(na.locf(ifelse(Event_desc=="NEB_SABA",date,NA),na.rm=F))) %>%
  group_by(ID) %>%
  mutate(reliever_use = replace_na(na.locf(temp_var,na.rm=F),median_reliever_use),
         time_since_last = date-temp_var2,
         NEB_SABA = ifelse(!is.na(time_since_last) & time_since_last<=90,1,0)) %>%
  select(-temp_var,-temp_var2,-time_since_last)
rm(median_reliever_use)

### Controller Use
events.wide<-events.wide %>%
  group_by(ID) %>% 
  mutate(temp_var = as.numeric(ifelse(Event_desc=="controllers",Event_value,NA)),
         controllers = replace_na(na.locf(temp_var, na.rm=F),0)) %>%
  select(-temp_var)

### Adherence
events.wide <-events.wide %>%
  mutate(temp_var = ifelse(Event_desc=="CSA_3",as.numeric(Event_value),NA),
         temp_var2 = ifelse(Event_desc=="CMA7_2",as.numeric(Event_value),NA),
         CSA_3 = replace_na(na.locf(temp_var,na.rm=F),0),
         CMA7_2 = replace_na(na.locf(temp_var2,na.rm=F),0)) %>%
  select(-temp_var,-temp_var2)

###  Respiratory Infections
temp <- events.wide %>%
  filter(Event_desc=="RI") %>%
  mutate(year = year(date)+1) %>% # want to link to the previous year
  count(ID,year,name="prev_year_ARI") %>%
  select(ID,year,prev_year_ARI)
temp2 <- events.wide %>%
  filter(Event_desc=="RI") %>%
  group_by(ID,year(date)) %>%
  mutate(this_year_ARI =row_number()) %>%
  ungroup %>% 
  select(ID,date,this_year_ARI)
events.wide<-left_join(events.wide %>% mutate(year=year(date)),temp)
events.wide<-left_join(events.wide,temp2) %>%
  group_by(ID,year) %>%
  mutate(this_year_ARI = na.locf(this_year_ARI,na.rm=F)) %>%
  group_by(ID) %>%
  mutate(recent_ARI = 1*(pmax(replace_na(prev_year_ARI,0),replace_na(this_year_ARI,0))>1),
         temp_var = as.numeric(date-lag(na.locf(ifelse(replace_na(RI,0)==1,date,NA),na.rm=F))),
         last_ARI=ifelse(is.na(temp_var),
                         "6.gt_2yrs_unknown",
                         ifelse(temp_var<14,
                                "1.lt_2wks",
                                ifelse(temp_var<60,
                                       "2.2wks_2mon",
                                       ifelse(temp_var<180,
                                              "3.2_6mon",
                                              ifelse(temp_var<365,
                                                     "4.6_12mon",
                                                     ifelse(temp_var<730,
                                                            "5.1_2yrs",
                                                            "6.gt_2yrs_unknown"))))))) %>%
  ungroup %>%
  select(-c(RI,prev_year_ARI,this_year_ARI,temp_var))
rm(temp,temp2)

### Asthma Encounters
temp<-events.wide %>%
  filter(Event_desc=="Asthma.Encounter") %>%
  count(ID,year, name="prev_year_encounters") %>%
  select(ID,year,prev_year_encounters) %>%
  mutate(year=year+1)  # want this to match to the year after, so add one to years
events.wide<-left_join(events.wide,temp)
temp2 <- events.wide %>%
  filter(Event_desc=="Asthma.Encounter") %>%
  group_by(ID,year) %>%
  mutate(this_year_encounters =row_number()) %>%
  ungroup %>% 
  select(ID,date,this_year_encounters)
events.wide<-left_join(events.wide,temp2) %>%
  group_by(ID,year) %>%
  mutate(this_year_encounters = na.locf(this_year_encounters,na.rm=F)) %>%
  group_by(ID) %>%
  mutate(recent_asthma_encounters = 1*(pmax(replace_na(prev_year_encounters,0),
                                         replace_na(this_year_encounters,0))>1))  %>%
  ungroup %>%
  select(-c(`Asthma.Encounter`,prev_year_encounters,this_year_encounters))
rm(temp,temp2)

### BTS Step
events.wide <-events.wide %>%
  group_by(ID) %>%
  mutate(temp_var = as.numeric(ifelse(Event_desc=="BTS Step",Event_value,NA)),
         BTS_Step = replace_na(na.locf(temp_var,na.rm=F),0)) %>%
  select(-temp_var,-`BTS Step`)

### Steroid prescriptions in the last two years
temp<-events.wide %>%
  filter(Event_desc=="steroid") %>%
  count(ID,year, name="prev_year_steroids") %>%
  select(ID,year,prev_year_steroids) %>%
  mutate(year=year+1)  # want this to match to the year after, so add one to years
events.wide<-left_join(events.wide,temp)
temp2 <- events.wide %>%
  filter(Event_desc=="steroid") %>%
  group_by(ID,year) %>%
  mutate(this_year_steroids =row_number()) %>%
  ungroup %>% 
  select(ID,date,this_year_steroids)
events.wide<-left_join(events.wide,temp2) %>%
  group_by(ID,year) %>%
  mutate(this_year_steroids = na.locf(this_year_steroids,na.rm=F)) %>%
  group_by(ID) %>%
  mutate(recent_steroids = 1*(pmax(replace_na(prev_year_steroids,0),
                                         replace_na(this_year_steroids,0))>1))  %>%
  ungroup %>%
  select(-c(prev_year_steroids,this_year_steroids))
rm(temp,temp2)

### Previous Asthma Attacks - Determined from Read Codes and PC prescriptions ONLY
events.wide<-events.wide %>%
  mutate(temp_var = ifelse(Event_desc %in% c("Asthma.attack","steroid"), 1,NA)) %>%
  group_by(ID) %>%
  # have they had asthma attacks before?
    mutate(temptemp = as.numeric(date-lag(na.locf(ifelse(temp_var==1,date,NA),na.rm=F))),
           last_PC_attack=ifelse(is.na(temptemp),
                                 "gt_2yrs_unknown",
                                 ifelse(temptemp<30,
                                        "lt_1mon",
                                        ifelse(temptemp<90,
                                               "1_3mon",
                                               ifelse(temptemp<180,
                                                      "3_6mon",
                                                      ifelse(temptemp<365,
                                                             "6_12mon",
                                                             ifelse(temptemp<730,
                                                                    "1_2yrs",
                                                                    "gt_2yrs_unknown"))))))) %>%
   select(-c(temptemp,temp_var,steroid,Asthma.attack))

######################################################################################
### Outcomes
######################################################################################

events.wide <- events.wide %>% filter(date<soft_right_censor) # 1 year prior to uniform end of data))
events.wide$soft_right_censor<-NULL
length(unique(events.wide$ID))

### Asthma Outcomes
events.wide<-events.wide  %>%
  group_by(ID) %>%
  mutate(temp=ifelse(Event_desc %in% c("A&E","Inpatient", "Asthma.Death","steroid"),1, NA), # ALL attacks + deaths
         last_attack = as.numeric(date-lag(na.locf(ifelse(temp==1,date,NA),na.rm=F))), # time since last attack of each kind (PC and ANY)
         ttemp = ifelse(!is.na(temp) & (is.na(last_attack) | last_attack>=14),1, 0), 
         next_attack = as.numeric(as.Date(lead(na.locf(ifelse(ttemp==1,date,NA),na.rm=F,fromLast = T)))-date)) %>%
  select(-c(temp,ttemp,last_attack,`A&E`,Inpatient,Asthma.Death,year))

events.wide<-as.data.frame(events.wide) %>% 
  group_by(ID,date) %>%
  mutate(FLAG = sum(Event_desc %in% c("steroid","A&E","Inpatient"))) %>% 
  filter(Event_desc=="zzz.GP Visit" & FLAG==0) %>%
  select(-FLAG) %>%
  ungroup()

events.wide$next_attack<-replace_na(events.wide$next_attack,9999999)
events.wide$attack_4wks<-as.numeric(events.wide$next_attack<=7*4)
events.wide$attack_12wks<-as.numeric(events.wide$next_attack<=7*12)
events.wide$attack_26wks<-as.numeric(events.wide$next_attack<=7*26)
events.wide$attack_52wks<-as.numeric(events.wide$next_attack<=7*52)

events.wide<-events.wide %>% select(-Event_value, -next_attack, -Event_desc) 
  
length(unique(events.wide$ID))

### check if any variables have missing data
apply(events.wide,2,function(x) sum(is.na(x)))

######################################################################################
### Demographics dataset 
######################################################################################

demographics <- read_csv("../../linked_data/1516-0489_PrimaryCareDemographics.csv")
names(demographics)[names(demographics)=="Index1"]<-"ID" 
# data cleaning
sum(is.na(demographics$Sex))==0
# Add a new variable indicating the end of their data
demographics <- demographics %>%
  mutate(temp = if_else(is.na(DeductionDate), as_date("2018-03-31"), as_date(DeductionDate))) %>%
  # calculate fake DOB from age and end of study per row
  mutate(DoB = temp - years(Age)) %>%
  # Most recent SIMD
  mutate(SIMD = ifelse(!is.na(SIMD2012quintile),
                       as.character(SIMD2012quintile),
                       ifelse(!is.na(SIMD2009quintile),
                              as.character(SIMD2009quintile),
                              "missing"))) %>%
  # urban rurality
  mutate(UR6 = ifelse(is.na(UR6_Code),
                      "missing",
                      as.character(UR6_Code))) %>%
  # choose baseline record - first registration
  mutate(rand = runif(nrow(demographics),0,1)) %>%
  arrange(ID,RegDate,rand) %>%
  group_by(ID) %>%
  dplyr::slice(n()) %>%
  dplyr::select(ID, Sex, DoB, DataZone, SIMD, UR6)
# Datazone to NUTS3 conversion
conversion <- read_csv("../../linked_data/Reference/datazone to nuts3.csv")
conversion <- conversion %>%
  dplyr::rename(DataZone = ZONECODE) %>%
  dplyr::rename(NUTS3 = NUTS3_CODE) %>%
  select(DataZone, NUTS3)
demographics<-left_join(demographics,conversion)
sum((!is.na(demographics$DataZone) & is.na(demographics$NUTS3)))==0
demographics$DataZone<-NULL
rm(conversion)

### Add in Demographics dataset
### Dropping people based on age
analysis_data<-left_join(events.wide,demographics) 
length(unique(analysis_data$ID)) 

analysis_data<-analysis_data %>% 
  mutate(age = floor(as.numeric(date-DoB)/365.25),
         NUTS3 = ifelse(is.na(NUTS3),"missing",NUTS3), 
         month = as.factor(month(date))) %>%
  filter(!is.na(DoB) & !is.na(Sex))
length(unique(analysis_data$ID)) 

analysis_data<-analysis_data %>%  
  filter(age>=18) %>%
  select(-DoB)  
length(unique(analysis_data$ID)) 

analysis_data$Obese<-(analysis_data$Obese=="Obese")*1
analysis_data$CMA7_2<-as.numeric(analysis_data$CMA7_2)
analysis_data$CSA_3<-as.numeric(analysis_data$CSA_3)
analysis_data$controllers<-as.numeric(analysis_data$controllers)

##################################################################
#      Split holdout and sensitivity
##################################################################

rm(list=setdiff(ls(),c("analysis_data","followup_1","followup_2")))

analysis_data<-analysis_data %>% 
  group_by(ID) %>% 
  mutate(first_sample = min(date),
         no_data_pre_COPD = (!is.na(dt1_copd) & first_sample>=dt1_copd),
         sens_flag = (COPD_ref=="COPD First" | no_data_pre_COPD)) %>%
  ungroup %>%
  select(-first_sample,-no_data_pre_COPD)

sensitivity_COPD_IDs<-unlist(unique(analysis_data %>% filter(sens_flag) %>% select(ID)))
followup_IDs<-unlist(unique(analysis_data %>% filter(!sens_flag) %>% select(ID)))
### 10% data for hold-out analysis
IDs_holdout<-followup_IDs[sample(1:length(followup_IDs),floor(length(followup_IDs)*0.1))]
IDs_train<-setdiff(followup_IDs,IDs_holdout)

# If you don't have any asthma consultations after the COPD diagnosis date, you should be removed from followup2
followup_2<-followup_2 %>% 
  inner_join(analysis_data %>% 
               filter(date>=dt1_copd) %>%
               distinct(ID))

followup_test<-followup_1  %>% filter(ID %in% IDs_holdout)
followup_train<-followup_1  %>% filter(ID %in% IDs_train)
followup_COPD<-followup_2 %>% filter(ID %in% sensitivity_COPD_IDs | ID %in% IDs_holdout)
save(followup_test,file="Data/followup_test.RData")
save(followup_train,file="Data/followup_train.RData")
save(followup_COPD,file="Data/followup_COPD.RData")
rm(followup_1,followup_2,followup_IDs,followup_test,followup_train,followup_COPD)

analysis_data_test<-analysis_data %>% filter(ID %in% IDs_holdout & (is.na(dt1_copd) | date<dt1_copd))
analysis_data_train<-analysis_data %>% filter(ID %in% IDs_train & (is.na(dt1_copd) | date<dt1_copd))
analysis_data_COPD<-analysis_data %>% filter(ID %in% sensitivity_COPD_IDs |
                                               (ID %in% IDs_holdout & !is.na(dt1_copd) & date>=dt1_copd))
nrow(analysis_data)==nrow(analysis_data_COPD)+nrow(analysis_data_test)+nrow(analysis_data_train)+
  nrow(analysis_data %>% filter(ID %in% IDs_train & !(is.na(dt1_copd) | date<dt1_copd)))
save(analysis_data_test,file="Data/analysis_data_test.RData")
save(analysis_data_train,file="Data/analysis_data_train.RData")
save(analysis_data_COPD,file="Data/analysis_data_COPD.RData")

##################################################################
#      Enrichment setup
##################################################################

## set primary outcome, and remove the others for this section
outcome<-as.factor(analysis_data$attack_52wks)
other_outcomes<-analysis_data %>% 
  select(c(names(analysis_data)[which(substr(names(analysis_data),1,6)=="attack")],"ID","dt1_copd","date"))
other_outcomes_train<-other_outcomes %>% filter(ID %in% IDs_train & (is.na(dt1_copd) | date<dt1_copd)) %>% 
  select(-date,-dt1_copd)
other_outcomes_test<-other_outcomes %>% filter(ID %in% IDs_holdout  & (is.na(dt1_copd) | date<dt1_copd)) %>% 
  select(-date,-dt1_copd)
save(other_outcomes_train,file="Data/other_outcomes_train.RData")
save(other_outcomes_test,file="Data/other_outcomes_test.RData")
rm(other_outcomes)


### one hot encoding for factors
analysis_data<-analysis_data %>% 
  select(-names(analysis_data)[which(substr(names(analysis_data),1,6)=="attack")])

x<-one_hot(as.data.table(lapply(analysis_data %>%
                                  select(Smoking, Peak_Flow, Blood_Eosinophil_Counts, last_ARI, 
                                         Eczema, Nasal.Polyps, Anxiety_Depression,
                                         GERD, Nasal.Spray, Rhinitis, Sex, SIMD, UR6, NUTS3, 
                                         month, last_PC_attack),
                                factor)))

analysis_data<-analysis_data %>% select(-c(Smoking,Peak_Flow,Blood_Eosinophil_Counts,last_ARI,
                                           Eczema,Nasal.Polyps, Anxiety_Depression,
                                           GERD, Nasal.Spray, Rhinitis, Sex,SIMD,UR6,NUTS3,
                                           month,last_PC_attack))
data<-as.data.frame(cbind(analysis_data,x,outcome))

test_data<-data %>% filter(ID %in% IDs_holdout  & (is.na(dt1_copd) | date<dt1_copd)) %>% 
  select(-ID,-date,-dt1_copd)
train_data<-data %>% filter(ID %in% IDs_train & (is.na(dt1_copd) | date<dt1_copd)) %>% 
 select(-date,-dt1_copd)   # keep ID in because we need it for the enrichment pre-analysis
COPD_data<-data %>% filter(ID %in% sensitivity_COPD_IDs |
                             (ID %in% IDs_holdout & !is.na(dt1_copd) & date>=dt1_copd)) %>%
  select(-ID,-date,-dt1_copd)
save(test_data,file="Data/test_data.RData")
save(train_data,file="Data/train_data.RData")
save(COPD_data,file="Data/COPD_data.RData")


