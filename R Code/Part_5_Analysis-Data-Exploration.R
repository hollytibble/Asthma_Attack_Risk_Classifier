##############################################################
#  Analysis Data Exploration
##############################################################

rm(list=ls())
library(epitools)

##############################################################
# People and records
##############################################################

# this is the number of records
# recall that train_data and test_data have removed the records AFTER COPD diagnosis (unlike analysis_data)
# and COPD_data is all the records after COPD diagnosis
load("Data/train_data.RData")
nrow(train_data)
load("Data/test_data.RData")
nrow(test_data)
load("Data/COPD_data.RData")
nrow(COPD_data)

# this is the number of people
load("Data/analysis_data_train.RData")
nrow(analysis_data_train %>% filter(is.na(dt1_copd) | date<dt1_copd) %>% distinct(ID))
load("Data/analysis_data_test.RData")
nrow(analysis_data_test %>% filter(is.na(dt1_copd) | date<dt1_copd) %>% distinct(ID))
load("Data/analysis_data_COPD.RData")
nrow(analysis_data_COPD %>% distinct(ID))

##############################################################
# Followup
##############################################################

load("Data/followup_train.RData")
sum(as.numeric(followup_train$followup))/365.25
summary(as.numeric(followup_train$followup))/365.25

load("Data/followup_test.RData")
sum(as.numeric(followup_test$followup))/365.25
summary(as.numeric(followup_test$followup))/365.25

load("Data/followup_COPD.RData")
sum(as.numeric(followup_COPD$followup))/365.25
summary(as.numeric(followup_COPD$followup))/365.25

rm(followup_COPD,followup_test,followup_train)

##############################################################
#  Demography
##############################################################

attack_ever_train<-analysis_data_train %>%
  group_by(ID) %>%
  summarise(value = sum(attack_52wks)>0) %>%
  ungroup %>%
  mutate(key = "Attack ever") %>%
  add_count(key,name="total") %>% 
  group_by(key,value) %>%
  summarise(train = paste0(n()," (",round(n()*100/total,2),"%)")) %>%
  distinct(key,value,train)

demographics_train<-analysis_data_train %>%
  # this is where we weed out those who were dropped from data_train as they had no eligible consultations before COPD
  filter((is.na(dt1_copd) | date<dt1_copd)) %>% 
  group_by(ID) %>%
  arrange(date) %>%
  dplyr::slice(1) %>%
  ungroup %>% 
  mutate(age_cat = case_when(age<=35 ~ "18-35",
                         age<=45 ~ "36-45",
                         age<=60 ~ "46-60",
                         age<=75 ~ "61-75",
                         age>75 ~ "76+")) %>%
  select(c(attack_52wks,age_cat,Sex, SIMD, UR6,
           names(analysis_data_train)[which(substr(names(analysis_data_train),1,4)=="Com.")],
           Anxiety_Depression, Eczema, GERD, Nasal.Polyps, Rhinitis,
           Obese, BTS_Step)) %>%
  gather("key","value") %>%
  rbind(analysis_data_train %>%
          group_by(ID) %>%
          summarise(value = sum(attack_52wks)>0) %>%
          ungroup %>%
          mutate(key = "Attack ever") %>%
          select(-ID)) %>%
  add_count(key,name="total") %>%
  group_by(key,value) %>%
  summarise(train = paste0(n()," (",round(n()*100/total,2),"%)")) %>%
  distinct(key,value,train)

demographics_test<-analysis_data_test %>%
  # this is where we weed out those who were dropped from data_train as they had no eligible consultations before COPD
  filter((is.na(dt1_copd) | date<dt1_copd)) %>% 
  group_by(ID) %>%
  arrange(date) %>%
  dplyr::slice(1) %>%
  ungroup %>% 
  mutate(age_cat = case_when(age<=35 ~ "18-35",
                             age<=45 ~ "36-45",
                             age<=60 ~ "46-60",
                             age<=75 ~ "61-75",
                             age>75 ~ "76+")) %>%
  select(c(attack_52wks,age_cat,Sex, SIMD, UR6,
           names(analysis_data_train)[which(substr(names(analysis_data_train),1,4)=="Com.")],
           Anxiety_Depression, Eczema, GERD, Nasal.Polyps, Rhinitis,
           Obese, BTS_Step)) %>%
  gather("key","value") %>%
  rbind(analysis_data_test %>%
          group_by(ID) %>%
          summarise(value = sum(attack_52wks)>0) %>%
          ungroup %>%
          mutate(key = "Attack ever") %>%
          select(-ID)) %>%
  add_count(key,name="total") %>%
  group_by(key,value) %>%
  summarise(test = paste0(n()," (",round(n()*100/total,2),"%)")) %>%
  distinct(key,value,test)

demographics_COPD<-rbind(analysis_data_test,analysis_data_COPD) %>%
  filter(!is.na(dt1_copd) & date>=dt1_copd) %>% 
  group_by(ID) %>%
  arrange(date) %>%
  dplyr::slice(1) %>%
  ungroup %>% 
  mutate(age_cat = case_when(age<=35 ~ "18-35",
                             age<=45 ~ "36-45",
                             age<=60 ~ "46-60",
                             age<=75 ~ "61-75",
                             age>75 ~ "76+")) %>%
  select(c(attack_52wks,age_cat,Sex, SIMD, UR6,
           names(analysis_data_train)[which(substr(names(analysis_data_train),1,4)=="Com.")],
          Anxiety_Depression, Eczema, GERD, Nasal.Polyps, Rhinitis,
           Obese, BTS_Step)) %>%
  gather("key","value") %>%
  rbind(rbind(analysis_data_test,analysis_data_COPD) %>%
          filter(!is.na(dt1_copd) & date>=dt1_copd) %>%
          group_by(ID) %>%
          summarise(value = sum(attack_52wks)>0) %>%
          ungroup %>%
          mutate(key = "Attack ever") %>%
          select(-ID))%>%
  add_count(key,name="total") %>%
  group_by(key,value)%>%
  summarise(COPD = paste0(n()," (",round(n()*100/total,2),"%)")) %>%
  distinct(key,value,COPD)


demographics_table<-demographics_train %>%
  left_join(demographics_test) %>%
  left_join(demographics_COPD)
write.csv(demographics_table,file="Data_to_extract/demographics_table.csv")

##############################################################
#  Demography by sample
##############################################################

demographics_train_samples<-analysis_data_train %>%
  # this is where we weed out those who were dropped from data_train as they had no eligible consultations before COPD
  filter((is.na(dt1_copd) | date<dt1_copd)) %>% 
  mutate(age_cat = case_when(age<=35 ~ "18-35",
                             age<=45 ~ "36-45",
                             age<=60 ~ "46-60",
                             age<=75 ~ "61-75",
                             age>75 ~ "76+")) %>%
  select(c(attack_52wks,age_cat,Sex, SIMD, UR6,
           names(analysis_data_train)[which(substr(names(analysis_data_train),1,4)=="Com.")],
           Anxiety_Depression, Eczema, GERD, Nasal.Polyps, Rhinitis,
           Obese, BTS_Step)) %>%
  gather("key","value") %>%
  add_count(key,name="total") %>%
  group_by(key,value) %>%
  summarise(train = paste0(n()," (",round(n()*100/total,2),"%)")) %>%
  distinct(key,value,train)

demographics_test_samples<-analysis_data_test %>%
  # this is where we weed out those who were dropped from data_train as they had no eligible consultations before COPD
  filter((is.na(dt1_copd) | date<dt1_copd)) %>% 
  mutate(age_cat = case_when(age<=35 ~ "18-35",
                             age<=45 ~ "36-45",
                             age<=60 ~ "46-60",
                             age<=75 ~ "61-75",
                             age>75 ~ "76+")) %>%
  select(c(attack_52wks,age_cat,Sex, SIMD, UR6,
           names(analysis_data_train)[which(substr(names(analysis_data_train),1,4)=="Com.")],
           Anxiety_Depression, Eczema, GERD, Nasal.Polyps, Rhinitis,
           Obese, BTS_Step)) %>%
  gather("key","value") %>%
  add_count(key,name="total") %>%
  group_by(key,value) %>%
  summarise(test = paste0(n()," (",round(n()*100/total,2),"%)")) %>%
  distinct(key,value,test)

demographics_COPD_samples<-rbind(analysis_data_test,analysis_data_COPD) %>%
  filter(!is.na(dt1_copd) & date>=dt1_copd) %>% 
  mutate(age_cat = case_when(age<=35 ~ "18-35",
                             age<=45 ~ "36-45",
                             age<=60 ~ "46-60",
                             age<=75 ~ "61-75",
                             age>75 ~ "76+")) %>%
  select(c(attack_52wks,age_cat,Sex, SIMD, UR6,
           names(analysis_data_train)[which(substr(names(analysis_data_train),1,4)=="Com.")],
           Anxiety_Depression, Eczema, GERD, Nasal.Polyps, Rhinitis,
           Obese, BTS_Step)) %>%
  gather("key","value") %>%
  add_count(key,name="total") %>%
  group_by(key,value)%>%
  summarise(COPD = paste0(n()," (",round(n()*100/total,2),"%)")) %>%
  distinct(key,value,COPD)


demographics_table_samples<-demographics_train_samples %>%
  left_join(demographics_test_samples) %>%
  left_join(demographics_COPD_samples)
write.csv(demographics_table_samples,file="Data_to_extract/demographics_table_samples.csv")

##############################################################
# Primary Care Contacts
##############################################################

PCC<-analysis_data_train %>%
  filter((is.na(dt1_copd) | date<dt1_copd)) %>% 
  mutate(year=year(date))  %>%
  group_by(ID,year) %>%
  summarise(Sex =first(Sex),
            SIMD =first(SIMD),
            UR6 = first(UR6),
            age = first(age),
            Smoking = first(Smoking),
            BTS = median(BTS_Step),
            n =  n())
summary(PCC$n)
model<-glm(n~as.numeric(SIMD)+as.numeric(UR6)+age+Sex+Smoking+BTS,
           PCC %>% filter(SIMD!="missing" & UR6!='missing' & year>2009 & year<2016),
           family=poisson(link="log"))
summary(model)
round(model$coefficients,3)
rm(PCC,model)

##############################################################
#  Time Between Primary Care Contacts
##############################################################

TBPC<-analysis_data_train %>%
  filter((is.na(dt1_copd) | date<dt1_copd)) %>% 
  group_by(ID) %>%
  mutate(gap = as.numeric(date-lag(date))) %>%
  filter(year(date)>2009 & !is.na(gap)) 
summary(TBPC$gap)
rm(TBPC)

##############################################################
#  Continuous Feature Distribution
##############################################################

summary(train_data$age)
summary(train_data$reliever_use)
summary(train_data$controllers)
summary(train_data$CSA_3)
summary(train_data$CMA7_2)

sum(train_data$controllers>16)*100/nrow(train_data)
sum(train_data$CSA_3>(0.91*4))*100/nrow(train_data)
sum(train_data$reliever_use>(690*4))*100/nrow(train_data)

##############################################################
#  Categorical Feature Distribution
##############################################################

categorical_train<-analysis_data_train %>% 
  filter((is.na(dt1_copd) | date<dt1_copd)) %>% 
  select(-c(names(analysis_data_train)[which(substr(names(analysis_data_train),1,6)=="attack")],
            ID,controllers,CMA7_2,CSA_3,age,reliever_use,date,dt1_copd,COPD_ref)) %>%
  gather("key","value") %>%
  add_count(key,name="total") %>%
  group_by(key,value)%>%
  summarise(train = paste0(n()," (",round(n()*100/total,2),"%)")) %>%
  distinct(key,value,train)

categorical_test<-analysis_data_test %>% 
  filter((is.na(dt1_copd) | date<dt1_copd)) %>% 
  select(-c(names(analysis_data_train)[which(substr(names(analysis_data_train),1,6)=="attack")],
            ID,controllers,CMA7_2,CSA_3,age,reliever_use,date,dt1_copd,COPD_ref)) %>%
  gather("key","value") %>%
  add_count(key,name="total") %>%
  group_by(key,value)%>%
  summarise(test = paste0(n()," (",round(n()*100/total,2),"%)")) %>%
  distinct(key,value,test)

categorical_COPD<-analysis_data_COPD %>% 
  filter(!is.na(dt1_copd) & date>=dt1_copd) %>% 
  select(-c(names(analysis_data_train)[which(substr(names(analysis_data_train),1,6)=="attack")],
            ID,controllers,CMA7_2,CSA_3,age,reliever_use,date,dt1_copd,COPD_ref)) %>%
  gather("key","value") %>%
  add_count(key,name="total") %>%
  group_by(key,value)%>%
  summarise(COPD = paste0(n()," (",round(n()*100/total,2),"%)")) %>%
  distinct(key,value,COPD)

categorical_table<-categorical_train %>%
  left_join(categorical_test) %>%
  left_join(categorical_COPD)
write.csv(categorical_table,file="Data_to_extract/categorical_table.csv")

