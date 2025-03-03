### pre-amble
rm(list=ls())
library(tidyverse)
library(lubridate)
library(zoo)
library(pracma)
library(readr)
library(doSNOW)
set.seed(12345)

### Load Dictionary Workspace
load("../Data/Dictionary.RData")

##########################################################################
###  Record Cleaning 
##########################################################################

prescribing <- read_csv("../../../linked_data/20190919 amended_pis/20190919 amended_pis.csv")
names(prescribing)[names(prescribing)=="Index10"]<-"ID"
nrow(prescribing) 
length(unique(prescribing$ID))

prescribing<-prescribing %>% filter(PIBNFChapterCode %in% c("03","06","12"))

### data cleaning
prescribing$PrescDate<-as.Date(prescribing$PrescDate)
prescribing$DispDate<-as.Date(prescribing$DispDate)
prescribing<-prescribing[which((format(prescribing$PrescDate,"%Y")!="2017" |
                                 format(prescribing$PrescDate,"%m")!="04") &
                                 as.numeric(format(prescribing$PrescDate,"%Y"))>=2009),]


### Load in the dataset which contains all dose information
dose <- read_csv("../../../linked_data/PIS/1516-0489_Final_dose_instructions_and_variables2.csv")
prescribing<-left_join(prescribing,dose)
rm(dose)

prescribing$dose<-toupper(prescribing$ePRNativeDoseInstructions)
prescribing$drugname<-ifelse(prescribing$PIApprovedName==prescribing$PIPrescribableItemName,
                                    prescribing$PIPrescribableItemName,
                                    paste(prescribing$PIPrescribableItemName,prescribing$PIApprovedName))
prescribing<-prescribing %>%
  select(-ePRNativeDoseInstructions,-PIApprovedName,-PIPrescribableItemName)

# Some have been deleted
prescribing<-prescribing[which(substr(prescribing$dose,1,7)!="DELETED"),]

nrow(prescribing) 
length(unique(prescribing$ID)) 

##########################################################################
###  Extracting steroid bursts
##########################################################################

steroid<-prescribing %>%
  filter(PIBNFParagraphCode=="0603020" & 
           !is.na(prescribed_quantity) &
           prescribed_quantity>0) %>%
  rename(date = PrescDate)

filter<-vapply(steroid_ind,
               function(x) str_detect(steroid$dose, x), 
               logical(nrow(steroid)))
steroid$flag<-apply(filter,1,sum)
steroid<-steroid %>%  filter(flag==0) %>% select(-flag)

load("../Data/Read Code Records.RData")
dates<-GP %>% group_by(ID) %>%  distinct(ID,date)
steroid<-inner_join(steroid,dates) 
rm(dates,filter,GP,right_censor_dates_diag)

steroid$strength<-NA
for (strength in c(2.5,500,125,120,100,80,40,30, 25,20,16,15,10,5,4,2,1)) {
  steroid$strength<-ifelse(is.na(steroid$strength) &
                             (str_detect(steroid$dose,paste0(strength," MG"))==T |
                                str_detect(steroid$dose,paste0(strength,"MG"))==T |
                                str_detect(steroid$dose,paste0(strength," MILLIGRAM"))==T |
                                str_detect(steroid$dose,paste0(strength,"MILLIGRAM"))==T |
                                str_detect(steroid$ePRNDName,paste0(strength," MG"))==T |
                                str_detect(steroid$ePRNDName,paste0(strength,"MG"))==T |
                                str_detect(steroid$ePRNDName,paste0(strength," MILLIGRAM"))==T |
                                str_detect(steroid$ePRNDName,paste0(strength,"MILLIGRAM"))==T),
                           strength,steroid$strength)
}
steroid$strength<-ifelse(is.na(steroid$strength) &
                           (str_detect(steroid$ePRNDName,"1 G") | str_detect(steroid$ePRNDName,"1G")),
                         1000,steroid$strength)
sum(is.na(steroid$strength))
for (strength in c(2.5,500,125,120,100,80,40,30, 25,20,16,15,10,5,4,2,1)) {
  steroid$strength<-ifelse(is.na(steroid$strength) &
                             (str_detect(steroid$`PIItemStrength/UOM`,paste0(strength," MG"))==T |
                                str_detect(steroid$`PIItemStrength/UOM`,paste0(strength,"MG"))==T |
                                str_detect(steroid$`PIItemStrength/UOM`,paste0(strength," MILLIGRAM"))==T |
                                str_detect(steroid$`PIItemStrength/UOM`,paste0(strength,"MILLIGRAM"))==T),
                           strength,steroid$strength)
}
sum(is.na(steroid$strength))
steroid$strength<-ifelse(steroid$`PIItemStrength/UOM`=="500 MICROGRAMS",0.5,steroid$strength)
sum(is.na(steroid$strength))

steroid<-steroid %>% 
  mutate(Event="steroid") %>%
  select(ID,date,Event)

save(steroid,file="../Data/steroid_bursts")
prescribing<-prescribing %>% filter(PIBNFChapterCode!="06")

##########################################################################
###   Drug Classification
##########################################################################

presc_unique<-prescribing %>% 
  select(drugname,dose,ePRNDName,PIDrugFormulation,
         PIBNFChapterCode,PIBNFSectionCode) %>% 
  add_count(drugname,dose,ePRNDName,PIDrugFormulation,
            PIBNFChapterCode,PIBNFSectionCode) %>%
  distinct() %>%
  arrange(-n)

# keep sprays & corticosteroids for rhinitis drugs
presc_unique_nasal<-presc_unique %>% filter(PIDrugFormulation=="SPRAY" | PIBNFSectionCode=="1202")
sum(presc_unique_nasal$n)
filter<-vapply(med_ICS,
               function(x) str_detect(presc_unique_nasal$drugname, x), 
               logical(nrow(presc_unique_nasal)))
presc_unique_nasal$flag<-apply(filter,1,sum)
presc_unique_nasal<-presc_unique_nasal %>% filter(flag==1)
prescribing_nasal<-inner_join(prescribing,presc_unique_nasal) %>%
  dplyr::select(ID,PrescDate) %>%
  dplyr::rename(date = PrescDate) %>%
  mutate(Event = "Nasal.Spray")
save(prescribing_nasal,file="../Data/prescribing_nasal")
rm(presc_unique_nasal,prescribing_nasal)

presc_unique <- presc_unique %>% filter(PIBNFChapterCode!="12")

# Identifying and assigning drug inclusion keywords
presc_unique$key<-""
for (med in med_list) {
  filter<-vapply(get(paste0("med_",med)),
                 function(x) str_detect(presc_unique$drugname, x), 
                 logical(nrow(presc_unique)))
  presc_unique[,paste0("flag_",med)]<-apply(filter,1,sum)
  for (medx in get(paste0("med_",med))) {
    filter<-vapply(get(paste0("brand_",medx)),
                   function(x) str_detect(presc_unique$drugname, x), 
                   logical(nrow(presc_unique)))
    presc_unique[,paste0("flag_",med)]<-(presc_unique[,paste0("flag_",med)]+apply(filter,1,sum))>0
    temp<-str_detect(presc_unique$drugname,medx) | as.logical(apply(filter,1,sum))
    presc_unique[temp,"key"]<-paste(presc_unique$key[temp],medx,sep="_")
  }
}
presc_unique$key<-ifelse(str_detect(presc_unique$drugname,"RELVAR ELLIPTA") |
                           str_detect(presc_unique$drugname,"VILANTEROL"),
                         paste(presc_unique$key,"VILANTEROL",sep="_"),
                         presc_unique$key)
presc_unique$key<-ifelse(str_detect(presc_unique$drugname,"MOXISLYTE"),
                         "_",
                         presc_unique$key)
presc_unique$key<-str_replace(substr(presc_unique$key,2,30)," ","_")

## How many flags does each prescription have?
presc_unique$flag_sum<-rowSums(presc_unique[,c("flag_LAMA","flag_ICS","flag_LABA", 
                                               "flag_LTRA", "flag_SABA", 
                                               "flag_Steroid", "flag_Theophylline","flag_MAb")])

### Drop if they aren't in any of these categories
presc_unique<-presc_unique[which(presc_unique$flag_sum!=0),]

### code the drug class in cases with only one flag
presc_unique$drug_class<-NA
for (med in med_list) {
  presc_unique$drug_class<-ifelse(presc_unique[,paste0("flag_",med)]==T & presc_unique$flag_sum==1,
                                  med,
                                  presc_unique$drug_class) 
}

# Code the two-flag cases
presc_unique$drug_class<-ifelse(presc_unique$flag_ICS==T & 
                                  presc_unique$flag_LABA==T,
                                "ICS+LABA",presc_unique$drug_class)
presc_unique$drug_class<-ifelse(presc_unique$flag_LAMA==T & 
                                  presc_unique$flag_SABA==T,"LAMA",presc_unique$drug_class)

sum(presc_unique %>% filter(drug_class=="MAb") %>% select(n))
presc_unique<-presc_unique %>% filter(drug_class!="MAb")
sum(presc_unique$n)
sum(presc_unique %>% filter(drug_class %in% c("ICS","ICS+LABA")) %>% select(n))

## Getting rid of records matching ICS exclusion keywords
for (keyword in exclusion_keywords) {
  #print(presc_unique %>% filter(str_detect(presc_unique$dose,keyword)) %>% select(dose))
  #print(presc_unique %>% filter(str_detect(presc_unique$drugname,keyword)) %>% select(drugname))
  #print(presc_unique %>% filter(str_detect(presc_unique$ePRNDName,keyword)) %>% select(ePRNDName))
  presc_unique<-presc_unique %>%
    filter((str_detect(dose,keyword)==F & 
             str_detect(drugname,keyword)==F &
             str_detect(ePRNDName,keyword)==F))
}
sum(presc_unique$n)
sum(presc_unique %>% filter(drug_class %in% c("ICS","ICS+LABA")) %>% select(n))

# Removing records matching exclusion brand names 
for (brand in exclusion_brands) {
  #print(presc_unique %>% filter(str_detect(presc_unique$drugname,brand)) %>% select(drugname))
  presc_unique<-presc_unique[which(str_detect(presc_unique$drugname,brand)==F),]
}
sum(presc_unique$n)
sum(presc_unique %>% filter(drug_class %in% c("ICS","ICS+LABA")) %>% select(n))

## Flagging brand names
presc_unique$brandname<-""
for (med in brand_list) {
  for (brand in get(paste0("brand_",med))) {
    presc_unique$brandname <- ifelse(str_detect(presc_unique$drugname, brand)==T,
                                     brand,
                                     presc_unique$brandname)
  }
}
presc_unique$brandname <-ifelse(presc_unique$brandname %in% c("","BECLOMETHASONE","CROMOGLICATE"),
                                "",
                                presc_unique$brandname)

presc_unique$disp_sub<-ifelse(presc_unique$brandname=="",
                                 0,
                                 ifelse(str_detect(presc_unique$ePRNDName,presc_unique$brandname),
                                        0,
                                        1))

### Changing the drug class of solutions
filters_solution<-vapply(formulation_keywords_sol,
                         function(x) str_detect(paste0(presc_unique$dose,
                                                       presc_unique$drugname,
                                                       presc_unique$ePRNDName),x), 
                         logical(nrow(presc_unique)))
sum(presc_unique %>% filter(drug_class=="ICS") %>% select(n))
presc_unique<-presc_unique %>%
  mutate(solution=apply(filters_solution,1,sum),
         drug_class=ifelse(drug_class=="ICS" & (is.na(solution) | solution>0 | PIDrugFormulation=="SOL"),
                           "ICS_SOL",
                           drug_class)) %>%
  select(-solution)
sum(presc_unique %>% filter(drug_class=="ICS") %>% select(n))

prescribing_asthma<-inner_join(prescribing,presc_unique) %>%
  dplyr::select(ID,key,PrescDate,drugname,ePRNDName,`PIItemStrength/UOM`,prescribed_quantity,
                dispensed_quantity,drug_class,ndx,brandname,disp_sub,dose,PIDrugFormulation) %>%
  ungroup
rm(presc_unique)

##########################################################################
###  Reliever medication use
##########################################################################

neb_saba<-prescribing_asthma %>%
  filter(str_detect(ePRNDName,"NEB|ORAL|TAB|SYRUP|SOL|INJ") |
           PIDrugFormulation %in%  c("TABS","SOL","SOLN","SYRUP") |
           (PIDrugFormulation=="CAPS" & !str_detect(ePRNDName,"INH")) |
           str_detect(`PIItemStrength/UOM`,"ML")) %>%
  dplyr::rename(date = PrescDate) %>%
  select(ID,date) %>%
  mutate(Event = "NEB_SABA",Event_value="1")
save(neb_saba,file="../Data/neb_saba")
rm(neb_saba)

reliever_use<-prescribing_asthma %>%
  filter(!str_detect(ePRNDName,"NEB|ORAL|TAB|SYRUP|SOL|INJ") &
           (PIDrugFormulation!="CAPS" | str_detect(ePRNDName,"INH")) &
           !PIDrugFormulation %in% c("TABS","SOL","SOLN","SYRUP") &
           str_detect(`PIItemStrength/UOM`,"ML")==F) %>%
  mutate(strength = ifelse(str_detect(ePRNDName,
                                      "400MICRO|400MCG|400 MICRO|400 MCG")>0,400,
                           ifelse(str_detect(ePRNDName,
                                             "200MICRO|200MCG|200 MICRO|200 MCG")>0,200,
                                  ifelse(str_detect(ePRNDName,
                                                    "100MICRO|100MCG|100 MICRO|100 MCG")>0,100,
                                         ifelse(str_detect(ePRNDName,
                                                           "95MICRO|95MCG|95 MICRO|95 MCG")>0,95,NA)))))
sum(is.na(reliever_use$strength))
reliever_use$strength = ifelse(is.na(reliever_use$strength),
                               ifelse(str_detect(reliever_use$`PIItemStrength/UOM`,
                                                 "400MICRO|400MCG|400 MICRO|400 MCG")>0,400,
                                      ifelse(str_detect(reliever_use$`PIItemStrength/UOM`,
                                                        "200MICRO|200MCG|200 MICRO|200 MCG")>0,200,
                                             ifelse(str_detect(reliever_use$`PIItemStrength/UOM`,
                                                               "100MICRO|100MCG|100 MICRO|100 MCG")>0,100,
                                                    ifelse(str_detect(reliever_use$`PIItemStrength/UOM`,
                                                                      "95MICRO|95MCG|95 MICRO|95 MCG")>0,95,NA)))),
                               reliever_use$strength)
reliever_use<-reliever_use %>% 
  mutate(strength = ifelse(brandname=="AIROMIR", 100,
                           ifelse(brandname=="PULVINAL SALBUTAMOL", 200,
                                  ifelse(brandname=="SALAMOL", 100,
                                         ifelse(brandname=="VENTOLIN" & strength==95, 100,
                                                ifelse(brandname=="VENTOLIN" & strength==400, 200,
                                                       ifelse(brandname=="SALBULIN", 100,
                                                              ifelse(brandname=="ASMASAL", 95,
                                                                     strength))))))))
# sum(is.na(reliever_use$strength))
# table(reliever_use$strength,useNA = "ifany")*100/nrow(reliever_use)


reliever_use$doses<-ifelse(str_detect(reliever_use$ePRNDName,
                                      "100 DOSE|100DOSE|100-DOSE")>0,100,
                           ifelse(str_detect(reliever_use$ePRNDName,
                                             "200 DOSE|200DOSE|200-DOSE")>0,200,
                                  ifelse(str_detect(reliever_use$ePRNDName,
                                                    "60 DOSE|60DOSE|60-DOSE")>0,60,
                                         ifelse(str_detect(reliever_use$ePRNDName,
                                                           "120 DOSE|120DOSE|120-DOSE")>0,120,NA))))
# table(reliever_use$doses,useNA = "ifany")*100/nrow(reliever_use)
# View(reliever_use %>% filter(!is.na(doses)) %>% 
#        count(brandname,strength,doses) %>% 
#        arrange(brandname,strength,n))

reliever_use$doses<-ifelse(is.na(reliever_use$doses),
                           ifelse(reliever_use$strength==200,
                                  ifelse(reliever_use$brandname=="PULVINAL SALBUTAMOL",100,
                                         ifelse(reliever_use$brandname=="VENTOLIN" |
                                                  reliever_use$brandname=="",60,200)),
                                  200),
                           reliever_use$doses)

reliever_use$qty<-ifelse(!is.na(reliever_use$dispensed_quantity) & 
                           reliever_use$dispensed_quantity>=1 & 
                           reliever_use$dispensed_quantity<=4,
                         reliever_use$dispensed_quantity,1)

reliever_use<-reliever_use %>%
  group_by(ID) %>%
  arrange(ID, PrescDate) %>%
  mutate(reliever_gap =as.numeric(lead(PrescDate)-PrescDate),
         obtained_doses = qty*doses,
         obtained_mcg = qty*doses*strength) %>%
  select(ID, PrescDate, reliever_gap, obtained_mcg,obtained_doses,qty) %>%
  group_by(ID, PrescDate) %>%
  summarise(qty = sum(qty),
            obtained_doses = sum(obtained_doses),
            obtained_mcg = sum(obtained_mcg),
            reliever_gap = sum(reliever_gap)) %>% 
  mutate(Event_value=obtained_mcg/as.numeric(reliever_gap),
         Event = "reliever_use") %>%
  dplyr::rename(date = PrescDate) %>%
  select(ID,date,Event,Event_value)
save(reliever_use,file="../Data/reliever_use")
rm(reliever_use)
prescribing_asthma<-prescribing_asthma %>% filter(drug_class!="SABA")

##########################################################################
###  ICS daily medicine amount used calculation
##########################################################################

prescribing_asthma_ICS<-prescribing_asthma %>%
  filter(drug_class %in% c("ICS","ICS+LABA"))
rownames(prescribing_asthma_ICS)<-NULL
nrow(prescribing_asthma_ICS)

prescribing_asthma<-prescribing_asthma %>%
  filter(!drug_class %in% c("ICS","ICS+LABA"))

# Dose frequency - as required/needed coded as minimum
filters_freq_one<-vapply(dose_freq_one,function(x) str_detect(prescribing_asthma_ICS$dose,x),
                         logical(nrow(prescribing_asthma_ICS)))
filters_freq_two<-vapply(dose_freq_two,function(x) str_detect(prescribing_asthma_ICS$dose,x),
                         logical(nrow(prescribing_asthma_ICS)))
filters_freq_four<-vapply(dose_freq_four,function(x) str_detect(prescribing_asthma_ICS$dose,x),
                          logical(nrow(prescribing_asthma_ICS)))
filters_freq_daily<-vapply(dose_freq_daily,function(x) str_detect(prescribing_asthma_ICS$dose,x),
                           logical(nrow(prescribing_asthma_ICS)))
prescribing_asthma_ICS$freq_dose<-ifelse(apply(filters_freq_one,1,sum)>=1,1,NA)
prescribing_asthma_ICS$freq_dose<-ifelse(is.na(prescribing_asthma_ICS$freq_dose) &
                                           (apply(filters_freq_two,1,sum)>=1 |
                                              ((str_detect(prescribing_asthma_ICS$dose, "MORN") |
                                                  str_detect(prescribing_asthma_ICS$dose, "AM") |
                                                  str_detect(prescribing_asthma_ICS$dose, "A.M") |
                                                  str_detect(prescribing_asthma_ICS$dose, "MANE")) &
                                                 (str_detect(prescribing_asthma_ICS$dose, "EVE") |
                                                    str_detect(prescribing_asthma_ICS$dose, "NIGHT") |
                                                    str_detect(prescribing_asthma_ICS$dose, "BEDTIME") |
                                                    str_detect(prescribing_asthma_ICS$dose, "PM") |
                                                    str_detect(prescribing_asthma_ICS$dose, "P.M") |
                                                    str_detect(prescribing_asthma_ICS$dose, "NOCTE")))),
                                         2,prescribing_asthma_ICS$freq_dose)
prescribing_asthma_ICS$freq_dose<-ifelse(apply(filters_freq_four,1,sum)>=1,4,prescribing_asthma_ICS$freq_dose)
prescribing_asthma_ICS$freq_dose<-ifelse(is.na(prescribing_asthma_ICS$freq_dose),
                                         ifelse(apply(filters_freq_daily,1,sum)>=1,
                                                1,
                                                prescribing_asthma_ICS$freq_dose),
                                         prescribing_asthma_ICS$freq_dose)
round(table(prescribing_asthma_ICS$freq_dose, useNA = "ifany")*100/nrow(prescribing_asthma_ICS),1)
# table(prescribing_asthma_ICS$freq_dose, prescribing_asthma_ICS$key, useNA = "ifany")
sum(is.na(prescribing_asthma_ICS$freq_dose) & prescribing_asthma_ICS$key %in% c("CICLESONIDE","FLUTICASONE_VILANTEROL"))
sum(is.na(prescribing_asthma_ICS$freq_dose) & !prescribing_asthma_ICS$key %in% c("CICLESONIDE","FLUTICASONE_VILANTEROL"))
prescribing_asthma_ICS$freq_dose<-ifelse(is.na(prescribing_asthma_ICS$freq_dose),
                                         ifelse(prescribing_asthma_ICS$key %in% c("CICLESONIDE","FLUTICASONE_VILANTEROL"),
                                                1,2),
                                         prescribing_asthma_ICS$freq_dose)
table(prescribing_asthma_ICS$freq_dose, useNA = "ifany")*100/nrow(prescribing_asthma_ICS)


## Dose Quantity
prescribing_asthma_ICS$daily_dose<-NA
dose_quantities<-c("1","2","3","4","ONE","TWO","THREE","FOUR")
dose_quantities_num<-rep(c(1,2,3,4),2)
for (k in 1:8) {
  prescribing_asthma_ICS$daily_dose<-ifelse(str_detect(prescribing_asthma_ICS$dose,
                                                       paste0("TAKE ",dose_quantities[k]))==TRUE |
                                              str_detect(prescribing_asthma_ICS$dose,
                                                         paste0("INHALE ",dose_quantities[k]))==TRUE |
                                              str_detect(prescribing_asthma_ICS$dose,
                                                         paste0(dose_quantities[k]," AT "))==TRUE |
                                              str_detect(prescribing_asthma_ICS$dose,
                                                         paste0(dose_quantities[k]," TO BE TAKEN "))==TRUE |
                                              str_detect(prescribing_asthma_ICS$dose,
                                                         paste0(dose_quantities[k],"PUF"))==TRUE |
                                              str_detect(prescribing_asthma_ICS$dose,
                                                         paste0(dose_quantities[k]," PUF"))==TRUE |
                                              str_detect(prescribing_asthma_ICS$dose,
                                                         paste0(dose_quantities[k]," P "))==TRUE |
                                              str_detect(prescribing_asthma_ICS$dose,
                                                         paste0(dose_quantities[k],"P "))==TRUE |
                                              str_detect(prescribing_asthma_ICS$dose,
                                                         paste0(dose_quantities[k]," DAILY"))==TRUE,
                                            dose_quantities_num[k],
                                            prescribing_asthma_ICS$daily_dose)
  for (keyword in dose_quant) {
    prescribing_asthma_ICS$daily_dose<-ifelse(str_detect(prescribing_asthma_ICS$dose,
                                                         paste0(dose_quantities[k]," ",keyword))==TRUE,
                                              dose_quantities_num[k],
                                              prescribing_asthma_ICS$daily_dose)
  }
}
round(table(prescribing_asthma_ICS$daily_dose, useNA = "ifany")*100/nrow(prescribing_asthma_ICS),1)
# table(prescribing_asthma_ICS$daily_dose, prescribing_asthma_ICS$key, useNA = "ifany")
list<-c("BUDESONIDE","CICLESONIDE","FLUTICASONE_SALMETEROL","FLUTICASONE_VILANTEROL","MOMETASONE")
sum(is.na(prescribing_asthma_ICS$daily_dose) & prescribing_asthma_ICS$key %in% list)
sum(is.na(prescribing_asthma_ICS$daily_dose) & !prescribing_asthma_ICS$key %in% list)
prescribing_asthma_ICS$daily_dose<-ifelse(is.na(prescribing_asthma_ICS$daily_dose),
                                          ifelse(prescribing_asthma_ICS$key %in% list,
                                                 1,
                                                 2),
                                          prescribing_asthma_ICS$daily_dose)
table(prescribing_asthma_ICS$daily_dose, useNA = "ifany")*100/nrow(prescribing_asthma_ICS)

## How much do they take per day
prescribing_asthma_ICS$daily_dose_units<-prescribing_asthma_ICS$freq_dose*prescribing_asthma_ICS$daily_dose

##########################################################################
###  Adherence: Part 1
##########################################################################

#sum(is.na(prescribing_asthma_ICS$dispensed_quantity))
#sum(is.na(prescribing_asthma_ICS$dispensed_quantity))*100/nrow(prescribing_asthma_ICS)
#sum(prescribing_asthma_ICS$prescribed_quantity!=prescribing_asthma_ICS$dispensed_quantity,na.rm=T)
#sum(prescribing_asthma_ICS$prescribed_quantity!=prescribing_asthma_ICS$dispensed_quantity,na.rm=T)*100/nrow(prescribing_asthma_ICS)
#sum(prescribing_asthma_ICS$prescribed_quantity!=prescribing_asthma_ICS$dispensed_quantity &prescribing_asthma_ICS$prescribed_quantity==0,na.rm=T)

prescribing_asthma_ICS$qty<-ifelse(!is.na(prescribing_asthma_ICS$dispensed_quantity),
                                   ifelse(prescribing_asthma_ICS$dispensed_quantity<1,
                                          1,
                                          prescribing_asthma_ICS$dispensed_quantity),
                                   ifelse(prescribing_asthma_ICS$prescribed_quantity<1,
                                          1,
                                          prescribing_asthma_ICS$prescribed_quantity))

round(table(prescribing_asthma_ICS$qty)*100/nrow(prescribing_asthma_ICS),1)
#table(prescribing_asthma_ICS$qty,prescribing_asthma_ICS$drug_class)

# doses per unit extraction
prescribing_asthma_ICS$doses_per_pack<-ifelse(prescribing_asthma_ICS$qty>=14,prescribing_asthma_ICS$qty,NA)
for (dose in pack_size_doses) {
  prescribing_asthma_ICS$doses_per_pack<-ifelse(is.na(prescribing_asthma_ICS$doses_per_pack) &
                                                  (str_detect(prescribing_asthma_ICS$ePRNDName, paste0(dose," DOSE"))==T |
                                                     str_detect(prescribing_asthma_ICS$ePRNDName, paste0(dose,"DOSE"))==T |
                                                     str_detect(prescribing_asthma_ICS$ePRNDName, paste0(dose,"-DOSE"))==T |
                                                     str_detect(prescribing_asthma_ICS$ePRNDName, paste0(dose," X "))==T) , 
                                                dose, 
                                                prescribing_asthma_ICS$doses_per_pack)
}

#View(prescribing_asthma_ICS %>% filter(is.na(doses_per_pack) & 
#                                         str_detect(ePRNDName,"DOSE") &
#                                         str_detect(ePRNDName,"/DOSE")==F &
#                                         str_detect(ePRNDName,"UNIT DOSE")==F &
#                                         str_detect(ePRNDName,"METERED DOSE")==F) %>% 
#       count(ePRNDName) %>% arrange(-n))

sum(!is.na(prescribing_asthma_ICS$doses_per_pack))*100/nrow(prescribing_asthma_ICS) # 15.2

##########################################################################
###  Medication Strength
##########################################################################

## Look for the values followed by keywords
prescribing_asthma_ICS$strength<-NA
for (k in doses_mcg) {
  prescribing_asthma_ICS$strength<-ifelse(is.na(prescribing_asthma_ICS$strength) &
                                            ((str_detect(prescribing_asthma_ICS$ePRNDName, paste0(k, "/"))==T &
                                                prescribing_asthma_ICS$drug_class=="ICS+LABA") |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste0(k, "MCG"))==T |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste(k, "MCG", sep=" "))==T |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste(k, "MICROGRAM", sep=" "))==T |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste0(k, "MICROGRAM"))==T),
                                          as.numeric(k),
                                          prescribing_asthma_ICS$strength)
}
for (k in doses_mg) {
  prescribing_asthma_ICS$strength<-ifelse(is.na(prescribing_asthma_ICS$strength) &
                                            (str_detect(prescribing_asthma_ICS$ePRNDName, paste0(k, "MG"))==T |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste(k, "MG", sep=" "))==T |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste0(k, "MILLIGRAM"))==T |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste(k, "MILLIGRAM", sep=" "))==T),
                                          as.numeric(k)*1000,
                                          prescribing_asthma_ICS$strength)
}
sum(is.na(prescribing_asthma_ICS$strength))*100/nrow(prescribing_asthma_ICS)
#View(prescribing_asthma_ICS %>% filter(is.na(strength)) %>% count(ePRNDName) %>% arrange(-n))

# from manual comparison against our appendix dose strength table
# a lot more need removing than in SIVEII because of the dose matching problem
nrow(prescribing_asthma_ICS)
prescribing_asthma_ICS<-prescribing_asthma_ICS %>%
  filter(!(key=="BECLOMETASONE" & !is.na(strength) & !strength %in% c(50,100,200,250,400))) %>%
  filter(!(key=="BECLOMETASONE_FORMOTEROL" & !is.na(strength) & !strength %in% c(100,200))) %>%
  filter(!(key=="BUDESONIDE" & !is.na(strength) & !strength %in% c(100,200,250,400))) %>%
  filter(!(key=="BUDESONIDE_FORMOTEROL" & !is.na(strength) & !strength %in% c(50,100,160,200,320,400))) %>%
  filter(!(key=="FLUTICASONE" & !is.na(strength) & !strength %in% c(50,100,125,250,500))) %>%
  filter(!(key=="FLUTICASONE_FORMOTEROL" & !is.na(strength) & !strength %in% c(50,125,250))) %>%
  filter(!(key=="FLUTICASONE_SALMETEROL" & !is.na(strength) & !strength %in% c(50,100,125,250,500))) %>%
  filter(!(key=="FLUTICASONE_VILANTEROL" & !is.na(strength) & !strength %in% c(92,184))) %>%
  filter(!(key=="MOMETASONE" & !is.na(strength) & !strength %in% c(200,400))) 
nrow(prescribing_asthma_ICS)
length(unique(prescribing_asthma_ICS$ID))

# mostly accuhalers and evohalers etc
for (k in doses_mcg) {
  prescribing_asthma_ICS$strength<-ifelse(is.na(prescribing_asthma_ICS$strength) &
                                            (str_detect(prescribing_asthma_ICS$ePRNDName, paste0(k," ACCUHALER")) |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste0(k," CLICKHALER")) |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste0(k," EVOHALER")) |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste0(k," TURBOHALER")) |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste0("QVAR ",k)) |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste0("SERETIDE ",k)) |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste0("INHAL ",k)) |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste0("EVOHALER ",k)) |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste0("SERETIDE MDI ",k)) |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste0("SERETIDE",k)) |
                                               str_detect(prescribing_asthma_ICS$ePRNDName, paste0("ALVESCO ",k))),
                                          as.numeric(k),
                                          prescribing_asthma_ICS$strength)
}
sum(is.na(prescribing_asthma_ICS$strength))*100/nrow(prescribing_asthma_ICS)

mode_strength_key<-prescribing_asthma_ICS %>%
  filter(!is.na(strength)) %>%
  dplyr::count(key,strength) %>%
  arrange(key,n) %>%
  group_by(key) %>%
  dplyr::slice(n()) %>%
  dplyr::rename(mode_strength = strength) %>%
  select(-n)
prescribing_asthma_ICS<-left_join(prescribing_asthma_ICS,mode_strength_key) %>%
  mutate(strength = ifelse(is.na(strength),mode_strength,strength)) %>%
  select(-mode_strength)
sum(is.na(prescribing_asthma_ICS$strength))
rm(mode_strength_key)

#Table of strengths
#table(prescribing_asthma_ICS$strength,prescribing_asthma_ICS$key, useNA = "ifany")

#round(prop.table(table(prescribing_asthma_ICS$strength,
#                       prescribing_asthma_ICS$PIApprovedName, useNA = "ifany"),
#                 margin=2)*100,1)


##########################################################################
###  Adherence: Part 2
##########################################################################

# doses per unit: modal values by brandname, key, and strength
temp_mode_doses_per_pack<-prescribing_asthma_ICS %>%
  filter(!is.na(doses_per_pack)) %>%
  add_count(key, strength, brandname, name="total") %>%
  add_count(key, strength, brandname, doses_per_pack, name = "n") %>%
  group_by(key, strength, brandname) %>%
  filter(n==max(n)) %>% 
  rename(mode_doses_per_pack = doses_per_pack) %>%
  mutate(percentage = round(100*n/total,3)) %>%
  select(key, mode_doses_per_pack, strength,brandname,percentage) %>%
  arrange(-percentage) %>%
  dplyr::slice(1) 

prescribing_asthma_ICS<-left_join(prescribing_asthma_ICS,temp_mode_doses_per_pack) %>%
  mutate(flag = ifelse(is.na(doses_per_pack) & !is.na(mode_doses_per_pack),1,0),
         doses_per_pack=ifelse(is.na(doses_per_pack),mode_doses_per_pack,doses_per_pack)) %>%
  select(-mode_doses_per_pack)

# imputed<-prescribing_asthma_ICS %>% filter(flag==1)
# summary(imputed$percentage)
# sum(imputed$percentage<60)
# sum(imputed$percentage<60)*100/nrow(imputed)
# imputed %>% filter(percentage<60) %>% count(key,brandname,strength)
# rm(imputed)
# sum(!is.na(prescribing_asthma_ICS$doses_per_pack))*100/nrow(prescribing_asthma_ICS) # 98.8

# impute missing values from medicines.org.uk
prescribing_asthma_ICS$doses_per_pack<-ifelse(is.na(prescribing_asthma_ICS$doses_per_pack) & 
                                                prescribing_asthma_ICS$key=="FLUTICASONE_FORMOTEROL",
                                              120,prescribing_asthma_ICS$doses_per_pack)
prescribing_asthma_ICS$doses_per_pack<-ifelse(is.na(prescribing_asthma_ICS$doses_per_pack) & 
                                                prescribing_asthma_ICS$key=="FLUTICASONE_SALMETEROL",
                                              ifelse(prescribing_asthma_ICS$strength==500,
                                                     60,120),
                                              prescribing_asthma_ICS$doses_per_pack)
prescribing_asthma_ICS$doses_per_pack<-ifelse(is.na(prescribing_asthma_ICS$doses_per_pack) & 
                                                prescribing_asthma_ICS$key=="BECLOMETASONE_FORMOTEROL",
                                              120,prescribing_asthma_ICS$doses_per_pack)
prescribing_asthma_ICS$doses_per_pack<-ifelse(is.na(prescribing_asthma_ICS$doses_per_pack) & 
                                                prescribing_asthma_ICS$key=="BECLOMETASONE",
                                              100,prescribing_asthma_ICS$doses_per_pack)
prescribing_asthma_ICS$doses_per_pack<-ifelse(is.na(prescribing_asthma_ICS$doses_per_pack) & 
                                                prescribing_asthma_ICS$key=="BUDESONIDE_FORMOTEROL",
                                              ifelse(prescribing_asthma_ICS$strength==400,
                                                     60,120),
                                              prescribing_asthma_ICS$doses_per_pack)
sum(is.na(prescribing_asthma_ICS$doses_per_pack))==0

prescribing_asthma_ICS$quantity<-ifelse(prescribing_asthma_ICS$qty<14,
                                        prescribing_asthma_ICS$qty*prescribing_asthma_ICS$doses_per_pack,
                                        prescribing_asthma_ICS$doses_per_pack)

## Save File to Use for Quantity coding and EHR prescriptions analyses
#save(prescribing_asthma_ICS,file="prescriptions_Secondary_Analysis.RData")

rm(list=setdiff(ls(),c("prescribing_asthma_ICS","prescribing_asthma")))

##########################################################################
###  Adherence
##########################################################################

# Calculating some prescription-level information
# and Combining the prescriptions on the same day
df<-prescribing_asthma_ICS %>%
  mutate(supply_duration = quantity/daily_dose_units) %>%
  group_by(ID) %>%
  arrange(ID, PrescDate) %>%
  mutate(interval_duration = as.numeric(lead(PrescDate)-PrescDate)) %>%
  select(ID,PrescDate,supply_duration,interval_duration) %>%
  mutate(supply_duration=ifelse(row_number()!=1 & !is.na(lag(interval_duration)) & lag(interval_duration==0),
                                supply_duration+lag(supply_duration),
                                supply_duration)) %>%
  filter(is.na(interval_duration) | interval_duration!=0)

CSA_measures<-df %>%
  arrange(ID,PrescDate) %>%
  mutate(CSA=supply_duration/interval_duration,
         CSA_cum=cumsum(CSA),
         CSA_3 = (CSA_cum-lag(CSA_cum,3))/3) %>%
  filter(!is.na(interval_duration)) %>%
  group_by(ID) %>%
  mutate(Event="CSA_3",
         row = row_number(),
         Event_value = as.character(ifelse(row==1,
                                           NA,
                                           ifelse(row==2,lag(CSA),
                                                  ifelse(row==3,(lag(CSA)+lag(CSA,n=2))/2,
                                                         ifelse(row==4,(lag(CSA)+lag(CSA,n=2)+lag(CSA,n=3))/3,
                                                                CSA_3)))))) %>%
  rename(date = PrescDate) %>%
  filter(!is.na(Event_value)) %>%
  select(ID,date,Event,Event_value)
save(CSA_measures,file="../Data/CSA_measures.RData")
rm(CSA_measures)

df<-df %>%
  mutate(refill=row_number(),
         excess = pmax(0,supply_duration-interval_duration),
         gap=pmax(0,interval_duration-supply_duration),
         cum_excess=cumsum(excess),
         cum_gap=cumsum(gap*(cum_excess>0 & gap>0)), 
         leftovers = ifelse(cum_excess>0 & 
                              gap>0 &
                              cum_excess>=cum_gap,
                            -1*gap,
                            excess),
         last=PrescDate) %>%
  select(ID,refill,PrescDate,supply_duration,leftovers,interval_duration,last)

# Here I extract the information needed for assembling the daily level data frame on which to affix the records
framex<-df %>% 
  group_by(ID) %>%
  arrange(ID,PrescDate) %>%
  dplyr::slice(1) %>% 
  mutate(startdate=PrescDate,
         end=as.numeric(as.Date("2017-03-31")-PrescDate)) %>% 
  select(ID,end,startdate)

cl<-makeCluster(16)   # detectCores()
registerDoSNOW(cl)
names<-unique(df$ID)

adherence<-foreach(i=1:length(names), 
                   .packages=c("dplyr","tidyr","data.table","zoo"),
                   .combine=rbind) %dopar% {
                     frame<-framex %>% filter(ID==names[i])
                     
                     # Here I am making the frame, and merging in the data.
                     CMA7<-setDT(frame)[,list(ID=ID,startdate=startdate,PrescDate=seq(1,end+1)), by=1:nrow(frame)] %>%
                       mutate(PrescDate=startdate+PrescDate-1, 
                              ID= as.character(ID)) %>%
                       select(-startdate,-nrow) %>%
                       left_join(df) %>%
                       arrange(PrescDate) %>%
                       # I am now passing down information to the days following a prescription event 
                       mutate(refill=na.locf(refill),
                              last=na.locf(last),
                              days_since_last = PrescDate-last,
                              leftovers = ifelse(is.na(leftovers),0,leftovers),
                              duration_locf=na.locf(supply_duration),
                              year= year(PrescDate))  %>%
                       group_by(year) %>%
                       filter(sum(days_since_last==0)>=1) %>%
                       ungroup  %>%
                       group_by(ID,year) %>%
                       mutate(supply_duration = ifelse(row_number()==1 & duration_locf>days_since_last,
                                                       duration_locf-days_since_last,supply_duration),
                              supply_all_over = ifelse(refill==min(refill),
                                                       supply_duration,
                                                       supply_duration+lag(cumsum(leftovers))),
                              duration_locf = na.locf(supply_duration,na.rm=F),
                              all_locf=na.locf(supply_all_over,na.rm=F),
                              supply2=ifelse(!is.na(supply_duration),1,pmax(0,(pmin(duration_locf*2,all_locf)-days_since_last)>=1)), # capped carryover
                              window_length = as.numeric(max(PrescDate)-min(PrescDate)+1),
                              CMA7_2=sum(supply2,na.rm = T)/window_length) %>%
                       distinct(ID,year,CMA7_2) 
                     
                   }  

stopCluster(cl)

adherence<-adherence %>% 
  mutate(date = as.Date(paste0("01-01-",as.character(year-1)),"%d-%m-%Y"),
         Event = "CMA7_2",
         Event_value = as.character(CMA7_2)) %>%
  ungroup %>% 
  select(ID,date,Event,Event_value)
save(adherence,file="../Data/CMA7.RData")
rm(adherence,df,cl,framex,names)

##########################################################################
###  ICS Dose Category
##########################################################################

## How much do they take per day
prescribing_asthma_ICS$daily_dose_max<-prescribing_asthma_ICS$daily_dose_units*prescribing_asthma_ICS$strength

prescribing_asthma_ICS$dose_category<-ifelse((prescribing_asthma_ICS$key %in% 
                                                c("BUDESONIDE","BECLOMETASONE") & 
                                                prescribing_asthma_ICS$brandname=="") |
                                               prescribing_asthma_ICS$brandname %in%
                                               c("PULMICORT","SYMBICORT","CLENIL","BECODISKS",
                                                 "BECLAZONE","AEROBEC", "FILAIR"),
                                             ifelse(prescribing_asthma_ICS$daily_dose_max<=400,"low",
                                                    ifelse(prescribing_asthma_ICS$daily_dose_max<=800,"medium",
                                                           ifelse(prescribing_asthma_ICS$daily_dose_max<=3200,"high",
                                                                  "unknown"))),
                                             NA)

prescribing_asthma_ICS$dose_category<-ifelse((prescribing_asthma_ICS$brandname=="" & 
                                                prescribing_asthma_ICS$key %in% 
                                                c("BECLOMETASONE_FORMOTEROL")) |
                                               prescribing_asthma_ICS$brandname %in%
                                               c("QVAR","FOSTAIR","PULVINAL BECLOMETASONE"),
                                             ifelse(prescribing_asthma_ICS$daily_dose_max<=200,"low",
                                                    ifelse(prescribing_asthma_ICS$daily_dose_max<=400,"medium",
                                                           ifelse(prescribing_asthma_ICS$daily_dose_max<=1600,"high",
                                                                  "unknown"))),
                                             prescribing_asthma_ICS$dose_category)

prescribing_asthma_ICS$dose_category<-ifelse(prescribing_asthma_ICS$key=="CICLESONIDE",
                                             ifelse(prescribing_asthma_ICS$daily_dose_max<=160,"low",
                                                    ifelse(prescribing_asthma_ICS$daily_dose_max<=320,"medium",
                                                           "unknown")),
                                             prescribing_asthma_ICS$dose_category)

prescribing_asthma_ICS$dose_category<-ifelse((prescribing_asthma_ICS$brandname=="" & 
                                                prescribing_asthma_ICS$key %in% 
                                                c("FLUTICASONE","FLUTICASONE_SALMETEROL",
                                                  "FLUTICASONE_FORMOTEROL")) |
                                               prescribing_asthma_ICS$brandname %in%
                                               c("FLIXOTIDE","FLUTIFORM","SERETIDE"),
                                             ifelse(prescribing_asthma_ICS$daily_dose_max<=200,"low",
                                                    ifelse(prescribing_asthma_ICS$daily_dose_max<=500,"medium",
                                                           ifelse(prescribing_asthma_ICS$daily_dose_max<=2000,"high",
                                                                  "unknown"))),
                                             prescribing_asthma_ICS$dose_category)

prescribing_asthma_ICS$dose_category<-ifelse(prescribing_asthma_ICS$key=="FLUTICASONE_VILANTEROL",
                                             ifelse(prescribing_asthma_ICS$daily_dose_max<46 |
                                                      prescribing_asthma_ICS$daily_dose_max>368,
                                                    "unknown",
                                                    ifelse(prescribing_asthma_ICS$daily_dose_max<=92,
                                                           "medium",
                                                           "high")),
                                             prescribing_asthma_ICS$dose_category)

prescribing_asthma_ICS$dose_category<-ifelse(prescribing_asthma_ICS$brandname=="ASMABEC",
                                             ifelse(prescribing_asthma_ICS$daily_dose_max<=200,"low",
                                                    ifelse(prescribing_asthma_ICS$daily_dose_max<=400,"medium",
                                                           "unknown")),
                                             prescribing_asthma_ICS$dose_category)

prescribing_asthma_ICS$dose_category<-ifelse(prescribing_asthma_ICS$brandname=="BUDELIN",
                                             ifelse(prescribing_asthma_ICS$daily_dose_max<400 |
                                                      prescribing_asthma_ICS$daily_dose_max>3200,
                                                    "unknown",
                                                    ifelse(prescribing_asthma_ICS$daily_dose_max<=800,
                                                           "medium",
                                                           "high")),
                                             prescribing_asthma_ICS$dose_category)

prescribing_asthma_ICS$dose_category<-ifelse(prescribing_asthma_ICS$key=="MOMETASONE",
                                             ifelse(prescribing_asthma_ICS$daily_dose_max<=400,"low",
                                                    ifelse(prescribing_asthma_ICS$daily_dose_max<=800,"medium",
                                                           "unknown")),
                                             prescribing_asthma_ICS$dose_category)

prescribing_asthma_ICS$dose_category<-ifelse(prescribing_asthma_ICS$brandname=="SIRDUPLA",
                                             ifelse(prescribing_asthma_ICS$daily_dose_max<250 |
                                                      prescribing_asthma_ICS$daily_dose_max>2000,
                                                    "unknown",
                                                    ifelse(prescribing_asthma_ICS$daily_dose_max<=500,
                                                           "medium",
                                                           "high")),
                                             prescribing_asthma_ICS$dose_category)

prescribing_asthma_ICS$dose_category<-ifelse(prescribing_asthma_ICS$brandname=="AIRFLUSAL",
                                             ifelse(prescribing_asthma_ICS$daily_dose_max<500 |
                                                      prescribing_asthma_ICS$daily_dose_max>2000,
                                                    "unknown",
                                                    "high"),
                                             prescribing_asthma_ICS$dose_category)

prescribing_asthma_ICS$dose_category<-ifelse(prescribing_asthma_ICS$key=="BUDESONIDE_FORMOTEROL",
                                             ifelse(prescribing_asthma_ICS$daily_dose_max<=320,
                                                    "low",
                                                    ifelse(prescribing_asthma_ICS$daily_dose_max<=640,
                                                           "medium",
                                                           ifelse(prescribing_asthma_ICS$daily_dose_max<=2560,
                                                                  "high",
                                                                  "unknown"))),
                                             prescribing_asthma_ICS$dose_category)

sum(is.na(prescribing_asthma_ICS$dose_category))

##########################################################################
###  Number of controller medication prescriptions in year
##########################################################################

controller_use <- prescribing_asthma_ICS %>%
  ungroup %>%
  mutate(date = as.Date(paste0("01-01-",(year(PrescDate)+1)),"%d-%m-%Y"),
         Event = "controllers") %>%
  arrange(ID, date) %>%
  group_by(ID, date) %>%
  add_count() %>%
  dplyr::slice(1) %>%
  mutate(Event_value = as.character(n)) %>%
  ungroup() %>%
  select(ID,PrescDate,date,Event,Event_value)
save(controller_use,file="../Data/controller_use")
rm(controller_use)

##########################################################################
###  Inclusion criteria
##########################################################################

### add the ICS/ICS+LABA into the not ICS records, rather than merging
### as some ICS/ICS+LABA records have been removed here
prescribing_asthma<-bind_rows(prescribing_asthma,prescribing_asthma_ICS)
rm(prescribing_asthma_ICS)
prescribing_asthma<-prescribing_asthma %>% select(ID,PrescDate,ePRNDName,key,
                                                  brandname,dose,dose_category,qty,drug_class)

inclusion <- prescribing_asthma %>%
  filter(drug_class %in% c("ICS","ICS+LABA")) %>%
  select(ID) %>%
  group_by(ID) %>%
  slice(1)
prescribing_asthma <- prescribing_asthma %>% filter(ID %in% inclusion$ID)

# now that we have made the secondary analysis file we can further restrict the data for speed
# keep only people with at least one GP record, or they will never have a query sample
load("../Data/Read Code Records.RData")
prescribing_asthma <- prescribing_asthma %>% filter(ID %in% GP_names)
nrow(prescribing_asthma) 
rm(GP,GP_names)

##########################################################################
###  BTS Steps
##########################################################################

### grace period of previous prescriptions
grace<-120

### how long has it been since your last prescription of each stage?
prescribing_asthma<-prescribing_asthma %>%
  group_by(ID) %>%
  arrange(ID,PrescDate) %>%
  mutate(last_ICS=replace_na(as.numeric(PrescDate-na.locf(as.Date(ifelse(drug_class=="ICS",
                                                                         PrescDate,NA)),na.rm=F)),9999),
         last_combo=replace_na(as.numeric(PrescDate-na.locf(as.Date(ifelse(drug_class=="ICS+LABA",
                                                                           PrescDate,NA)),na.rm=F)),9999),
         last_ICS_low=replace_na(as.numeric(PrescDate-na.locf(as.Date(ifelse(drug_class %in% c("ICS","ICS+LABA") &
                                                                               dose_category=="low",
                                                                             PrescDate,NA)),na.rm=F)),9999),
         last_ICS_medium=replace_na(as.numeric(PrescDate-na.locf(as.Date(ifelse(drug_class %in% c("ICS","ICS+LABA") &
                                                                                  dose_category=="medium",
                                                                                PrescDate,NA)),na.rm=F)),9999),
         last_ICS_high=replace_na(as.numeric(PrescDate-na.locf(as.Date(ifelse(drug_class %in% c("ICS","ICS+LABA") &
                                                                                dose_category=="high",
                                                                              PrescDate,NA)),na.rm=F)),9999),
         last_ICS_dose=ifelse(last_ICS_low<last_ICS_high & last_ICS_low<last_ICS_medium,"low",
                              ifelse(last_ICS_medium<last_ICS_high,"medium","high")),
         last_LTRA=replace_na(as.numeric(PrescDate-na.locf(as.Date(ifelse(drug_class=="LTRA",
                                                                          PrescDate,NA)),na.rm=F)),9999),
         last_LABA=replace_na(as.numeric(PrescDate-na.locf(as.Date(ifelse(drug_class=="LABA",
                                                                          PrescDate,NA)),na.rm=F)),9999),
         last_LAMA=replace_na(as.numeric(PrescDate-na.locf(as.Date(ifelse(drug_class=="LAMA",
                                                                          PrescDate,NA)),na.rm=F)),9999),
         last_Theo=replace_na(as.numeric(PrescDate-na.locf(as.Date(ifelse(drug_class=="Theophylline",
                                                                          PrescDate,NA)),na.rm=F)),9999),
         last_ICS_SOL=replace_na(as.numeric(PrescDate-na.locf(as.Date(ifelse(drug_class=="ICS_SOL",
                                                                             PrescDate,NA)),na.rm=F)),9999))


### Number of add-on therapies in last grace period
prescribing_asthma<-prescribing_asthma %>%
  mutate(number_add_on=(as.numeric(last_LAMA<=grace)+ 
                          as.numeric(last_LTRA<=grace)+
                          as.numeric(last_ICS_SOL<=grace)+
                          as.numeric(last_Theo<=grace)))

### BTS Step Coding 
### STEP 0 
prescribing_asthma<-prescribing_asthma %>%  
  mutate(BTS_Step=ifelse(last_ICS>grace &
                           last_combo>grace,
                         0,NA)) 

### STEP 1 
prescribing_asthma<-prescribing_asthma %>% 
  mutate(BTS_Step=ifelse(last_ICS<=grace & 
                           last_ICS_dose=="low" &
                           last_LABA>grace &
                           last_combo>grace & 
                           number_add_on==0,  
                         1,BTS_Step))


### STEP 2 - low ICS + LABA (inc. as combo)
prescribing_asthma<-prescribing_asthma %>% 
  mutate(BTS_Step=ifelse(((last_ICS<=grace & last_LABA<=grace) | 
                            last_combo<=grace) &
                           last_ICS_dose=="low",
                         2,BTS_Step)) 


### STEP 3 
prescribing_asthma<-prescribing_asthma %>% 
  mutate(BTS_Step=ifelse((last_ICS<=grace & 
                            last_ICS_dose=="low" &
                            last_LABA>grace &
                            last_combo>grace &
                            number_add_on>0) |
                           ((last_combo<=grace | last_ICS<=grace) &
                              last_ICS_dose=="medium" & 
                              number_add_on==0), 3,BTS_Step))



### STEP 4
prescribing_asthma<-prescribing_asthma %>% 
  mutate(BTS_Step=ifelse((last_ICS<=grace | last_combo<=grace) &
                           (last_ICS_dose=="high" |
                              (last_ICS_dose=="medium" &
                                 number_add_on>0)),
                         4,BTS_Step))

### If the BTS_Step is missing, and the next step is the same as the previous step,
### carry the last step through
prescribing_asthma<-prescribing_asthma %>%
  group_by(ID) %>% 
  mutate(BTS_Step = ifelse(is.na(BTS_Step),
                           ifelse(!is.na(na.locf(BTS_Step,na.rm=F)),
                                  ifelse(!is.na(na.locf(BTS_Step,na.rm=F, fromLast = T)) & 
                                           na.locf(BTS_Step,na.rm=F)==na.locf(BTS_Step,na.rm=F, fromLast = T),
                                         na.locf(BTS_Step,na.rm=F),
                                         0),
                                  0),
                           BTS_Step))


##########################################################################
###  Final Cleaning
##########################################################################

prescribing_asthma <- prescribing_asthma %>%
  select(ID, PrescDate, BTS_Step) %>%
  # put it into the same format as the other datasets
  rename(date = PrescDate) %>%
  rename(Event_value = BTS_Step)
prescribing_asthma$Event<-"BTS Step"
prescribing_asthma$Event_value<-as.character(prescribing_asthma$Event_value)

# don't keep multiple records per day
prescribing_asthma<-prescribing_asthma %>%
  group_by(ID,date) %>%
  filter(Event_value==max(Event_value)) %>%
  ungroup

rm(grace)
save.image("../Data/BTS Step Records.RData")
