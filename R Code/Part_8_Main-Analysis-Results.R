########################################################################################
#      Training Data Enrichment Working Script
#      Created by Holly Tibble
########################################################################################

rm(list=ls()) # remove all variables from workspace

# library boots
library(tidyverse)
library(scales)
library(ranger)
library(ROSE)
library(pROC)
library(tibble)
library(predtools)

##################################################################
#      Main Results: threshold
##################################################################

results <- read_csv("../Data_to_extract/main_results v6.csv")
results[,1]<-NULL

AUC<-results$AUC
write.csv(AUC,file="../Data_to_extract/main_results AUC.csv")
round(summary(AUC*100),1)

Prevalence<-results$Prevalence
write.csv(Prevalence,file="../Data_to_extract/main_results Prevalence")
round(summary(Prevalence),1)

outcomes<-results %>% 
  select(TP:FN) %>%
  mutate(temp=ifelse(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))==0,
                     1,sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))),
         MCC = (TP*TN - FP*FN)*100/temp,
         accuracy = (TP+TN)*100/(TN+FP+TP+FN),
         sensitivity = TP*100/(TP+FN),
         specificity = TN*100/(TN+FP),
         PPV = TP*100/(FP+TP),
         NPV = TN*100/(FN+TN),
         balanced_acc = (sensitivity+specificity)/2) %>%
  select(-temp,-TP,-TN,-FP,-FN)
write.csv(outcomes,file="../Data_to_extract/Main_results_performance.csv")

round(summary(outcomes$MCC),1)
round(summary(outcomes$accuracy),1)
round(summary(outcomes$sensitivity),1)
round(summary(outcomes$specificity),1)
round(summary(outcomes$PPV),1)
round(summary(outcomes$NPV),1)
round(summary(outcomes$balanced_acc),1)

### feature importance
Importance<-results[,3:174] %>%
  mutate(iteration = row_number()) %>%
  gather("variable","value",-iteration) %>%
  mutate(key = ifelse(substr(variable,(nchar(variable)-1),nchar(variable))=="_p","p_value","coefficient"),
         variable = ifelse(substr(variable,(nchar(variable)-1),nchar(variable))=="_p",
                       substr(variable,1,(nchar(variable)-2)),variable)) %>%
  spread(key="key", value = "value")
write.csv(Importance,file="../Data_to_extract/Main_results_importance.csv")

Average_Importance <-Importance %>%
  group_by(variable) %>%
  summarise(mean_c = round(mean(coefficient),5),
            min_c = round(min(coefficient),5),
            max_c = round(max(coefficient),5),
            perc_sig = sum(p_value<0.05)*100/n()) 
View(Average_Importance %>% mutate(abs = abs(mean_c)) %>% arrange(-abs))

##################################################################
#      Full retrain
##################################################################

load("../Data/train_data.RData")
train_data<-train_data %>%
  select(-COPD_ref,
         -Anxiety_Depression_Never,
         -Blood_Eosinophil_Counts_Missing,
         -Eczema_Never,
         -GERD_Never,
         -Nasal.Polyps_Never,
         -Nasal.Spray_Never,
         -Peak_Flow_missing,
         -Rhinitis_Never,
         -last_PC_attack_gt_2yrs_unknown,
         -Sex_M,
         -month_12,
         -Smoking_Never,
         -last_ARI_6.gt_2yrs_unknown,
         -SIMD_3,
         -SIMD_missing,
         -UR6_3,
         -NUTS3_UKM25,
         -NUTS3_missing,
         -sens_flag,
         -ID)

train_data$CSA_3<-ifelse(train_data$CSA_3>(0.91*4),(0.91*4),train_data$CSA_3)
train_data$controllers<-ifelse(train_data$controllers>16,16,train_data$controllers)
train_data$reliever_use<-ifelse(train_data$reliever_use>(690*4),(690*4),train_data$reliever_use)

train_data$CSA_3<-scale(train_data$CSA_3, center = 0, scale = (0.91*4))
train_data$controllers<-scale(train_data$controllers, center = 0, scale = 16)
train_data$reliever_use<-scale(train_data$reliever_use, center = 0, scale = (690*4))
train_data$age<-scale(train_data$age, center = 18, scale = 102)

model_full<-glm(outcome ~ ., data=train_data, family=binomial(link="logit"))
save(model_full,file="../Data/Full model v5.RData")

model_summary<-as.data.frame(coef(summary(model_full))) %>%
  rownames_to_column()
write.csv(model_summary,file="../Data_to_extract/full_model_summary.csv")
View(model_summary %>% select(-`Std. Error`, -`z value`) %>% 
       mutate(Estimate = round(Estimate,3),
              `Pr(>|z|)` = round(`Pr(>|z|)`,3),
              Abs = abs(Estimate),
              OR = round(exp(Estimate),3)) %>%
       arrange(-OR))

########################################################################################

train_data2<-train_data %>%
  select(-Anxiety_Depression_Longer_than_5_Years_Ago,
         -Eczema_Longer_than_5_Years_Ago,
         -GERD_Longer_than_5_Years_Ago,
         -Nasal.Polyps_Longer_than_5_Years_Ago,
         -Nasal.Polyps_Longer_than_5_Years_Ago,
         -Rhinitis_Longer_than_5_Years_Ago)

View(as.data.frame(coef(summary(glm(outcome ~ ., data=train_data2, family=binomial(link="logit"))))) %>%
       rownames_to_column() %>% 
       select(-`Std. Error`, -`z value`) %>% 
       mutate(Estimate = round(Estimate,3),
              `Pr(>|z|)` = round(`Pr(>|z|)`,3),
              Abs = abs(Estimate),
              OR = round(exp(Estimate),3)))

########################################################################################

rm(list=setdiff(ls(),c("model_summary","model_full")))
load("../Data/test_data.RData")

test_data$CSA_3<-ifelse(test_data$CSA_3>(0.91*4),(0.91*4),test_data$CSA_3)
test_data$controllers<-ifelse(test_data$controllers>16,16,test_data$controllers)
test_data$reliever_use<-ifelse(test_data$reliever_use>(690*4),(690*4),test_data$reliever_use)

test_data$CSA_3<-scale(test_data$CSA_3, center = 0, scale = (0.91*4))
test_data$controllers<-scale(test_data$controllers, center = 0, scale = 16)
test_data$reliever_use<-scale(test_data$reliever_use, center = 0, scale = (690*4))
test_data$age<-scale(test_data$age, center = 18, scale = 102)

predictions<-predict(model_full,test_data,type="response")

ROC<-as.data.frame(cbind(outcome = test_data$outcome,
                         prob = predictions))
write.csv(ROC,file="../Data_to_extract/holdout_ROC.csv")
roc.curve((ROC$outcome-1),ROC$prob)

results<-read_csv("../Data_to_extract/enrichment_results v5.csv", col_types = cols(...1 = col_skip()))
threshold<-median(results$thr_GLMb)
rm(results)

stats<-as.data.frame(t(c(auc(test_data$outcome,predictions),
                         sum(test_data$outcome==1)*100/nrow(test_data),
                         as.numeric(sum(test_data$outcome==1 & predictions>threshold)),
                         as.numeric(sum(test_data$outcome==0 & predictions<=threshold)),
                         as.numeric(sum(test_data$outcome==0 & predictions>threshold)),
                         as.numeric(sum(test_data$outcome==1 & predictions<=threshold)))))
names(stats)<-c("AUC","Prev","TP","TN","FP","FN")
stats<-stats %>%
  mutate(temp=ifelse(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))==0,
                     1,sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))),
         MCC = (TP*TN - FP*FN)/temp,
         accuracy = (TP+TN)*100/(TN+FP+TP+FN),
         sensitivity = TP*100/(TP+FN),
         specificity = TN*100/(TN+FP),
         PPV = TP*100/(FP+TP),
         NPV = TN*100/(FN+TN),
         balanced_acc = (sensitivity+specificity)/2) %>%
  select(-temp)
write.csv(stats,file="../Data_to_extract/holdout_results v5.csv")

discrimination<-test_data %>%
  select(Com.Chronic.pulmonary.disease,recent_ARI,
         recent_asthma_encounters,BTS_Step,recent_steroids,
         Smoking_Current,Smoking_Former,Smoking_Never,
         Peak_Flow_missing,Blood_Eosinophil_Counts_Missing,
         last_PC_attack_gt_2yrs_unknown,COPD_ref,outcome) %>%
  cbind(predictions) %>%
  mutate(outcome=as.numeric(outcome)-1,
         group = case_when(outcome==1 & predictions>threshold ~ "TP",
                           outcome==1 & predictions<=threshold ~ "FN",
                           outcome==0 & predictions>threshold ~ "FP",
                           outcome==0 & predictions<=threshold ~ "TN"))


summary(lm(outcome~predictions,data=discrimination))

ggplot(discrimination) + 
  geom_smooth(aes(x=predictions,y=outcome)) +
  geom_line(aes(x=predictions,y=predictions)) 

########################################################################################

density_plot<-discrimination %>% select(predictions,outcome)
write.csv(density_plot,file="../Data_to_extract/density_plot.csv")
ggplot(discrimination) + geom_density(aes(x=predictions,fill=as.factor(outcome)),alpha=0.5)

subgroup_analysis<-rbind(discrimination %>% count(Com.Chronic.pulmonary.disease,group) %>%
                           mutate(var = "Com.Chronic.pulmonary.disease") %>% rename(level = Com.Chronic.pulmonary.disease),
                         discrimination %>% count(COPD_ref,group) %>%
                           mutate(var = "COPD_ref") %>% rename(level = COPD_ref),
                         discrimination %>% count(recent_ARI,group) %>%
                           mutate(var = "recent_ARI") %>% rename(level = recent_ARI),
                         discrimination %>% count(recent_asthma_encounters,group) %>%
                           mutate(var = "recent_asthma_encounters") %>% rename(level = recent_asthma_encounters),
                         discrimination %>% count(BTS_Step,group) %>%
                           mutate(var = "BTS_Step") %>% rename(level = BTS_Step),
                         discrimination %>% count(recent_steroids,group) %>%
                           mutate(var = "recent_steroids") %>% rename(level = recent_steroids),
                         discrimination %>% count(Peak_Flow_missing,group) %>%
                           mutate(var = "Peak_Flow_missing") %>% rename(level = Peak_Flow_missing),
                         discrimination %>% count(Blood_Eosinophil_Counts_Missing,group) %>%
                           mutate(var = "Blood_Eosinophil_Counts_Missing") %>% rename(level = Blood_Eosinophil_Counts_Missing),
                         discrimination %>% count(last_PC_attack_gt_2yrs_unknown,group) %>%
                           mutate(var = "last_PC_attack_gt_2yrs_unknown") %>% rename(level = last_PC_attack_gt_2yrs_unknown),
                         discrimination %>% count(Smoking_Current,Smoking_Former,Smoking_Never,group) %>%
                           mutate(var = "Smoking",
                                  level = case_when(Smoking_Current==1 ~ "Current",
                                                    Smoking_Former==1 ~ "Former",
                                                    Smoking_Never==1 ~ "Never")) %>%
                           select(var,level,group,n)) %>%
  pivot_wider(names_from = group,values_from = n,values_fill = 0) %>%
  mutate(prevalence = round((TP + TN + FP + FN)*100/nrow(test_data),1),
         attack_rate = round((TP + FP)*100/(TP + TN + FP + FN),1),
         sensitivity = round(TP*100/(TP+FN),1),
         specificity = round(TN*100/(TN+FP),1),
         PPV = round(TP*100/(FP+TP),1),
         NPV = round(TN*100/(FN+TN),1)) %>%
  select(-TP,-TN,-FP,-FN)
write.csv(subgroup_analysis,file="../Data_to_extract/subgroup_analysis.csv")

########################################################################################

calib1<-discrimination %>% select(outcome,predictions,Com.Chronic.pulmonary.disease)
write.csv(calib1,file="../Data_to_extract/calib1.csv")
calibration_plot(data=calib1,obs="outcome",pred="predictions",group="Com.Chronic.pulmonary.disease")

calib2<-discrimination %>% select(outcome,predictions,BTS_Step)
write.csv(calib2,file="../Data_to_extract/calib2.csv")
calibration_plot(data=calib2,obs="outcome",pred="predictions",group="BTS_Step")

calib3<-discrimination %>% select(outcome,predictions,Peak_Flow_missing)
write.csv(calib3,file="../Data_to_extract/calib3.csv")
calibration_plot(data=calib3,obs="outcome",pred="predictions",group="Peak_Flow_missing")

calib4<-discrimination %>% select(outcome,predictions,Blood_Eosinophil_Counts_Missing)
write.csv(calib4,file="../Data_to_extract/calib4.csv")
calibration_plot(data=calib4,obs="outcome",pred="predictions",group="Blood_Eosinophil_Counts_Missing")

calib5<-discrimination %>% select(outcome,predictions,recent_ARI)
write.csv(calib5,file="../Data_to_extract/calib5.csv")
calibration_plot(data=calib5,obs="outcome",pred="predictions",group="recent_ARI")

calib6<-discrimination %>% select(outcome,predictions,recent_asthma_encounters)
write.csv(calib6,file="../Data_to_extract/calib6.csv")
calibration_plot(data=calib6,obs="outcome",pred="predictions",group="recent_asthma_encounters")

calib7<-discrimination %>% select(outcome,predictions,recent_steroids)
write.csv(calib7,file="../Data_to_extract/calib7.csv")
calibration_plot(data=calib7,obs="outcome",pred="predictions",group="recent_steroids")

calib8<-discrimination %>% select(outcome,predictions,Smoking_Never)
write.csv(calib8,file="../Data_to_extract/calib8.csv")
calibration_plot(data=calib8,obs="outcome",pred="predictions",group="Smoking_Never")

calib9<-discrimination %>% select(outcome,predictions,last_PC_attack_gt_2yrs_unknown)
write.csv(calib9,file="../Data_to_extract/calib9.csv")
calibration_plot(data=calib9,obs="outcome",pred="predictions",group="last_PC_attack_gt_2yrs_unknown")

calib10<-discrimination %>% select(outcome,predictions,COPD_ref)
write.csv(calib10,file="../Data_to_extract/calib10.csv")
calibration_plot(data=calib10,obs="outcome",pred="predictions",group="COPD_ref")
