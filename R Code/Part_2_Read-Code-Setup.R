### pre-amble
library(stringr)
library(lubridate)
library(readr)
library(dplyr)
rm(list=ls())
set.seed(12345)

### Load Dictionary Workspace
load("../Data/Dictionary.RData")

### Read Codes Dataset Loading
GP <- read_csv("../../../linked_data/1516-0489_Readcodes.csv")
names(GP)[names(GP)=="Index1"]<-"ID"

#####################################################################
###  Data Cleaning
#####################################################################

### Dates
GP$date<-as.Date(GP$EventDate,format="%m/%d/%Y")
min(GP$date) 
max(GP$date)  
GP<-GP[which((format(GP$date,"%Y")<2017) | 
               (format(GP$date,"%Y")==2017 & as.numeric(format(GP$date,"%m"))<4)),] 

### Missing Read Codes
GP<-GP[-which(is.na(GP$ReadCode)),] 

### Remove duplicates 
GP <- distinct(GP, ID, EventDate, ReadCode, Data1, Data2, Data3, .keep_all = TRUE) 

### Check Read Codes formatted correctly
GP<-GP[which(str_length(GP$ReadCode)==5),]

length(unique(GP$ID)) #48975

#####################################################################
###  Potential Exclusions: COPD diagnosis
#####################################################################

### COPD
COPD<-GP %>% 
  filter(ReadCode %in% codes_COPD) %>%
  group_by(ID) %>%
  summarise(dt1_copd = min(EventDate))

#####################################################################
###  Asthma Resolution
#####################################################################

GP$asthma_resolution_flag <- ifelse(GP$ReadCode=="212G.", 1, 0)

### find dates of asthma resolution 
right_censor_dates_diag <- GP %>%
  filter(asthma_resolution_flag==1) %>%
  group_by(ID) %>%
  mutate(right_censor_diag = min(date)) %>%
  select(ID, right_censor_diag) %>%
  dplyr::slice(1)
GP$asthma_resolution_flag <- NULL

#####################################################################
###  Asthma Encounters
#####################################################################

### Asthma Encounters
### have to do seperately and merge in as codes in this list will also be used later
asthma_encounter_codes<-unique(c(codes_asthma_diagnosis,codes_asthma_related))
asthma_encounters<-GP[which(GP$ReadCode %in% asthma_encounter_codes),]
asthma_encounters$Event<-"Asthma.Encounter"
asthma_encounters <- asthma_encounters %>%
  select(c("ID","date","Event")) %>%
  group_by(ID, date) %>%
  dplyr::slice(1)
length(unique(asthma_encounters$ID))

#####################################################################
###  Read Coding 
#####################################################################

### Initialise my new data columns
GP$Event<-NA
GP$Event_value<-NA

### Asthma Attack
GP$Event<-ifelse(GP$ReadCode %in% codes_asthma_attack, "Asthma.attack", GP$Event)

### Other Chronic Pulmonary Diseases
GP$Event<-ifelse(GP$ReadCode %in% codes_charlson_pulmonary, "Com.Chronic.pulmonary.disease", GP$Event)


### Nasal Polpys
GP$Event<-ifelse(GP$ReadCode %in% codes_nasal_polyps, "Nasal.Polyps", GP$Event)


### Peak Flow Recording
GP$Event<-ifelse(GP$ReadCode %in% codes_peak_flow & 
                   (GP$Data1>20 & GP$Data1<2000),
                 "Peak Flow", GP$Event)
GP$Event_value<-ifelse(GP$ReadCode %in% codes_peak_flow & 
                         (GP$Data1>20 & GP$Data1<2000),
                       as.character(GP$Data1), GP$Event_value)

# Respiratory Infections
GP$Event<-ifelse(GP$ReadCode %in% codes_RI, "RI", GP$Event)

### Smoking
GP$Event<-ifelse(GP$ReadCode %in% codes_smoking_all, "Smoking", GP$Event)
GP$Event_value<-ifelse(GP$ReadCode %in% codes_smoking_never, "Never", 
                       ifelse(GP$ReadCode %in% codes_smoking_former, "Former",
                              ifelse(GP$ReadCode %in% codes_smoking_current, "Current",
                                     ifelse(GP$ReadCode %in% codes_smoking_current_conditional & 
                                              (GP$Data1>0 | (!is.na(GP$Data2) & GP$Data2!="0")),
                                            "Current", 
                                            GP$Event_value))))

### BMI
### raw value - will introduce some NAs as am converting some characters to numbers
GP$Event<-ifelse(GP$ReadCode %in% codes_BMI_val, "BMI", GP$Event)
GP$Event_value<-ifelse(GP$ReadCode %in% codes_BMI_val & !is.na(GP$Data1) & GP$Data1>10 & GP$Data1<100,
                       as.character(GP$Data1),
                       ifelse(GP$ReadCode %in% codes_BMI_val & !is.na(GP$Data3) & GP$Data3!="kg/m2" & GP$Data3!="Unknown" &
                                as.numeric(GP$Data3)>10 & as.numeric(GP$Data3)<100,
                              as.character(as.numeric(GP$Data3)),
                              GP$Event_value))
### categorical
GP$Event<-ifelse(GP$ReadCode %in% codes_BMI_cat, "BMI Category", GP$Event)
GP$Event_value<-ifelse(GP$ReadCode %in% codes_BMI_under, "Underweight / Low BMI", GP$Event_value)
GP$Event_value<-ifelse(GP$ReadCode %in% codes_BMI_normal, "Normal BMI", GP$Event_value)
GP$Event_value<-ifelse(GP$ReadCode %in% codes_BMI_overweight, "Overweight / High BMI", GP$Event_value)
GP$Event_value<-ifelse(GP$ReadCode %in% codes_BMI_obese, "Obese", GP$Event_value)
### height and weight
GP$Event<-ifelse(GP$ReadCode %in% codes_height, "Height", GP$Event)
GP$Event_value<-ifelse(GP$ReadCode %in% codes_height & is.na(GP$Data1),GP$Event_value,
                       ifelse(GP$ReadCode %in% codes_height & GP$Data1>22000, GP$Event_value,
                              ifelse(GP$ReadCode %in% codes_height & GP$Data1>3000, as.character(GP$Data1/100),
                                     ifelse(GP$ReadCode %in% codes_height & !is.na(GP$Data3) & GP$Data3=="cm" & GP$Data1>30 & GP$Data1<272, as.character(GP$Data1),
                                            ifelse(GP$ReadCode %in% codes_height & GP$Data1<2.2,as.character(GP$Data1*100), 
                                                   ifelse(GP$ReadCode %in% codes_height & GP$Data1>30 & GP$Data1<272,as.character(GP$Data1),
                                                          GP$Event_value))))))
GP$Event<-ifelse(GP$ReadCode %in% codes_weight, "Weight", GP$Event)
GP$Event_value<-ifelse(GP$ReadCode %in% codes_weight & is.na(GP$Data1) & is.na(GP$Data3),GP$Event_value,
                       ifelse(GP$ReadCode %in% codes_weight & !is.na(GP$Data3) & (GP$Data3=="Kg" | GP$Data3=="kg") & GP$Data1>20 & GP$Data1<300, as.character(GP$Data1), 
                              ifelse(GP$ReadCode %in% codes_weight & GP$Data1>20 & GP$Data1<300, as.character(GP$Data1),
                                     GP$Event_value)))

### Rhinitis
GP$Event<-ifelse(GP$ReadCode %in% codes_rhinitis, "Rhinitis", GP$Event)

### Anxiety/Depression
GP$Event<-ifelse(GP$ReadCode %in% codes_anx_dep, "Anxiety_Depression", GP$Event)

### Eczema / Dermitits
GP$Event<-ifelse(GP$ReadCode %in% codes_eczema, "Eczema", GP$Event)

### GERD
GP$Event<-ifelse(GP$ReadCode %in% codes_GERD, "GERD", GP$Event)


### Blood Eosinophil Counts
GP$Event<-ifelse(GP$ReadCode %in% codes_eosinophils & 
                   (is.na(GP$Data3) | GP$Data3!="MEA037") & 
                   !is.na(GP$Data1),
                 "Blood.Eosinophil.Counts", GP$Event)
GP$Event_value<-ifelse(GP$ReadCode %in% codes_eosinophils & !is.na(GP$Data1),
                       as.character(GP$Data1),
                       GP$Event_value)

#####################################################################
###  Format data and save 
#####################################################################

### get rid of Read codes for non-specified events and superfluous variables
GP<-GP[which(!is.na(GP$Event)),c("ID","date","Event","Event_value")]

# Get rid of ones which are missing a value, when required
GP<-GP %>% filter(!Event %in% c("BMI","Height","Weight",
                                "Smoking","Blood.Eosinophil.Counts") | 
                    !is.na(Event_value)) 

### add in asthma encounters
GP<- full_join(GP,asthma_encounters)

GP<-left_join(GP,COPD)
rm(COPD)

GP_names<-unlist(GP %>% 
  filter(Event %in% c("Asthma.Encounter","Asthma.attack","Blood.Eosinophil.Counts","Peak Flow","RI")) %>% 
  distinct(ID))
GP<-GP %>% filter(ID %in% GP_names)

### sort by person and date
rm(list=setdiff(ls(), c("GP_names","GP","right_censor_dates_diag")))

save.image("../Data/Read Code Records.RData")
