########################################################################################
#      Training Data Enrichment Working Script
#      Created by Holly Tibble
########################################################################################

rm(list=ls()) # remove all variables from workspace

# library boots
library(e1071)
library(ggplot2)
library(pROC)
library(dplyr)
library(gtools)
library(ROSE)
library(foreach)
library(doParallel)
library(tidyr)
library(stringr)
library(doSNOW)
library(pryr)

##################################################################
#      Dataset Specification
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
         -sens_flag)

names<-names(train_data %>% select(-ID,-outcome))

round(summary(train_data$CSA_3),2)
round(summary(train_data$controllers),2)
round(summary(train_data$reliever_use),2)

train_data$CSA_3<-ifelse(train_data$CSA_3>(0.91*4),(0.91*4),train_data$CSA_3)
train_data$controllers<-ifelse(train_data$controllers>16,16,train_data$controllers)
train_data$reliever_use<-ifelse(train_data$reliever_use>(690*4),(690*4),train_data$reliever_use)

train_data$CSA_3<-scale(train_data$CSA_3, center = 0, scale = (0.91*4))
train_data$controllers<-scale(train_data$controllers, center = 0, scale = 16)
train_data$reliever_use<-scale(train_data$reliever_use, center = 0, scale = (690*4))
train_data$age<-scale(train_data$age, center = 18, scale = 102)

# create list of samples to include in each iteration
ID_list<-unique(train_data$ID)
indices<-replicate(100,ID_list[sample(1:length(ID_list),
                                     floor(length(ID_list)*0.9))], 
                   simplify = FALSE)

results<-read_csv("../Data_to_extract/enrichment_results v5.csv", col_types = cols(...1 = col_skip()))
summary(results$thr_GLMb)
threshold<-median(results$thr_GLMb)
rm(results)

##################################################################
#     Classification Models
##################################################################

cl<-makeCluster(4)   
# detectCores()
registerDoSNOW(cl)
start=Sys.time()
results<-foreach(i=1:100, .combine=rbind, 
                 .packages=c("pROC","e1071","dplyr","tibble","tidyr"),
                 .errorhandling = "pass") %dopar% {

                   # split the data into two partitions: training and validation
                   train<-train_data %>% 
                     filter(ID %in% as.character(indices[[i]])) %>%  
                     select(-ID)
                   train<-train[vapply(train, function(x) length(unique(x))>1, logical(1L))]
                   
                   yvalidate <- unlist(train_data %>% 
                                         filter(!ID %in% as.character(indices[[i]])) %>%  
                                         select(outcome))
                   xvalidate<-train_data %>% 
                     filter(!ID %in% as.character(indices[[i]])) %>%  
                     select(-ID,-outcome)
                   
                   ### Train Model
                   model<-glm(outcome ~ ., data=train, family=binomial(link="logit"))
                   # save the assessment of the predictions
                   predictions_model<-predict(model,xvalidate,type="response")
                   coefficients<-unlist(data.frame(model$coefficients) %>% 
                                        select(model.coefficients))
                   p_values<-coef(summary(model))[,4]
                   
                   df<-c(TP=as.numeric(sum(yvalidate==1 & predictions_model>threshold)),
                         TN=as.numeric(sum(yvalidate==0 & predictions_model<=threshold)),
                         FP=as.numeric(sum(yvalidate==0 & predictions_model>threshold)),
                         FN=as.numeric(sum(yvalidate==1 & predictions_model<=threshold)))
                   
                   c(auc(as.numeric(yvalidate),predictions_model)[1],
                     sum(train$outcome==1)*100/nrow(train),
                     coefficients,
                     p_values,
                     df)
                 }

stopCluster(cl)
Sys.time()-start

results<-as.data.frame(results)
rownames(results)<-NULL
names(results)<-c("AUC","Prevalence",
                  "Intercept",names,
                  "Intercept_p",paste0(names,"_p"),
                  "TP","TN","FP","FN")
rm(list=setdiff(ls(),"results"))
write.csv(results,file="../Data_to_extract/main_results v6.csv")