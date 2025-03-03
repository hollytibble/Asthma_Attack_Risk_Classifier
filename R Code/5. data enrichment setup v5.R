########################################################################################
#      Training Data Enrichment Working Script
#      Created by Holly Tibble
########################################################################################

rm(list=ls()) # remove all variables from workspace

# library boots
library(e1071)
library(pROC)
library(ranger)
library(doParallel)
library(doSNOW)
library(gtools)
library(ROSE)
library(foreach)
library(tidyr)
library(stringr)
library(DescTools)
library(pryr)
library(xgboost)
library(performanceEstimation)
library(tidyverse)

##################################################################
#      Dataset Specification
##################################################################

load("../Data/train_data.RData")

train_data$COPD_ref<-NULL

train_data$CSA_3<-ifelse(train_data$CSA_3>(0.91*4),(0.91*4),train_data$CSA_3)
train_data$controllers<-ifelse(train_data$controllers>16,16,train_data$controllers)
train_data$reliever_use<-ifelse(train_data$reliever_use>(690*4),(690*4),train_data$reliever_use)

train_data$CSA_3<-as.numeric(scale(train_data$CSA_3, center = 0, scale = (0.91*4)))
train_data$controllers<-as.numeric(scale(train_data$controllers, center = 0, scale = 16))
train_data$reliever_use<-as.numeric(scale(train_data$reliever_use, center = 0, scale = (690*4)))
train_data$age<-as.numeric(scale(train_data$age, center = 18, scale = 102))

train_data<-train_data %>%
  select(-Nasal.Polyps_Longer_than_5_Years_Ago, 
         -Eczema_Longer_than_5_Years_Ago)

# create list of samples to include in each iteration
ID_list<-unique(train_data$ID)
indices<-replicate(10,ID_list[sample(1:length(ID_list),
                             floor(length(ID_list)*0.9))], 
                   simplify = FALSE)

prevalence<-round(sum(train_data$outcome==1)/nrow(train_data),3)

##################################################################
#     Classification Models
##################################################################

cl<-makeCluster(4)   
# detectCores()
registerDoSNOW(cl)
start=Sys.time()
results<-foreach(k=1:40, .combine=rbind, 
                 .packages=c("pROC","e1071","smotefamily","randomForest","performanceEstimation",
                             "ranger","xgboost","dplyr","DescTools","tidyverse"),
                 .errorhandling = "pass") %dopar% {
                   
                   i<-((k-1) %% 10)+1
                   enrichment<-floor((k-1)/10)+1
                   
                   ### predictions
                   predictions_func<-function(temp,data){
                     if(substr(temp$call,1,10)[1] =="naiveBayes") {
                       predictions<-predict(temp,data,type="raw")[,2]
                     } else if (substr(temp$call,1,6)[1] %in% c("ranger","glm(fo")) {
                       predictions<-predict(temp,data,type="response")$predictions[,2]
                     } else if (substr(temp$call,1,3)[1]=="glm") {
                       predictions<-predict(temp,data,type="response")
                     } else if (substr(temp$call,1,3)[1]=="xgb") {
                       predictions<-predict(temp,
                                            xgb.DMatrix(data=as.matrix(data),
                                                        label=rep(0,nrow(data))))
                     } else {
                       predictions<-attr(predict(temp,newdata=data,probability=TRUE), "probabilities")[,2]
                     }
                     return(predictions)
                   }
                   
                   ### function for finding optimum threshold
                   thresholder<-function(predictions) { 
                     ### locate the minimum of the function using a Golden Section Line Search
                     result <- optimize(
                       MCC_threshold, # the function to be minimized
                       c(0,1), # the bounds on th function parameter
                       predictions=predictions,
                       maximum=TRUE, # we are concerned with the function maxima
                       tol=1e-8) # the size of the final bracketing
                     ### re-assign to the half of the validation set not used in the optimisation
                     return(result$maximum)
                     rm(result)
                   }
                   
                   ### thresholding version of the MCC measure
                   MCC_threshold<-function(threshold,predictions){
                     TP<-as.numeric(sum(train_x$outcome==1 & predictions>threshold))
                     TN<-as.numeric(sum(train_x$outcome==0 & predictions<=threshold))
                     FP<-as.numeric(sum(train_x$outcome==0 & predictions>threshold))
                     FN<-as.numeric(sum(train_x$outcome==1 & predictions<=threshold))
                     temp<-ifelse(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))==0,
                                  1,sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
                     MCC <- (TP*TN - FP*FN)/temp
                     return(MCC)
                   }
                   
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
                     select(-ID,-outcome,-sens_flag)
                   
                   # specify the appropriate enriched data for training
                   invisible(ifelse(enrichment==1, train_x<-train, # version 1
                                    ifelse(enrichment==2, 
                                           train_x<-SMOTE(train[,-102],train[,102],K=5, dup_size = 2)$data,
                                           ifelse(enrichment==3, 
                                                  train_x<-SMOTE(train[,-102],train[,102],K=5, dup_size = 3)$data,  
                                                  train_x<-SMOTE(train[,-102],train[,102],K=5, dup_size = 4)$data))))
                   rm(train)
                   invisible(ifelse(enrichment==1, 
                          train_x<-train_x, 
                          train_x<-train_x %>% mutate(outcome = as.factor(class)) %>% select(-class)))
                   
                   
                   ### Naive Bayes Classifier
                   # Train the model - no hyperparameters to tune
                   nbc<-naiveBayes(outcome ~ ., data=train_x)
                   # save the assessment of the predictions
                   predictions_nbc_train<-predictions_func(nbc,train_x[,-ncol(train_x)])
                   threshold_nbc<-invisible(thresholder(predictions_nbc_train))
                   balanced_threshold_nbc<-mean(c(threshold_nbc,prevalence))
                   predictions_nbc<-predictions_func(nbc,xvalidate)
                   nbc_stats<-c(BrierScore(as.numeric(yvalidate)-1,predictions_nbc),
                                auc(as.numeric(yvalidate),predictions_nbc)[1],
                                as.numeric(sum(yvalidate==1 & predictions_nbc>0.5)),
                                as.numeric(sum(yvalidate==0 & predictions_nbc<=0.5)),
                                as.numeric(sum(yvalidate==0 & predictions_nbc>0.5)),
                                as.numeric(sum(yvalidate==1 & predictions_nbc<=0.5)),
                                as.numeric(sum(yvalidate==1 & predictions_nbc>threshold_nbc)),
                                as.numeric(sum(yvalidate==0 & predictions_nbc<=threshold_nbc)),
                                as.numeric(sum(yvalidate==0 & predictions_nbc>threshold_nbc)),
                                as.numeric(sum(yvalidate==1 & predictions_nbc<=threshold_nbc)),
                                as.numeric(sum(yvalidate==1 & predictions_nbc>prevalence)),
                                as.numeric(sum(yvalidate==0 & predictions_nbc<=prevalence)),
                                as.numeric(sum(yvalidate==0 & predictions_nbc>prevalence)),
                                as.numeric(sum(yvalidate==1 & predictions_nbc<=prevalence)),
                                as.numeric(sum(yvalidate==1 & predictions_nbc>balanced_threshold_nbc)),
                                as.numeric(sum(yvalidate==0 & predictions_nbc<=balanced_threshold_nbc)),
                                as.numeric(sum(yvalidate==0 & predictions_nbc>balanced_threshold_nbc)),
                                as.numeric(sum(yvalidate==1 & predictions_nbc<=balanced_threshold_nbc)),
                                threshold_nbc,balanced_threshold_nbc)
                   rm(nbc,predictions_nbc,threshold_nbc,predictions_nbc_train,balanced_threshold_nbc)
                   
                   ### Logistic Regression
                   LR<-glm(outcome ~ ., data=train_x, family=binomial(link="logit"))
                   # save the assessment of the predictions
                   predictions_LR_train<-predictions_func(LR,train_x[,-ncol(train_x)])
                   threshold_LR<-invisible(thresholder(predictions_LR_train))
                   balanced_threshold_LR<-mean(c(threshold_LR,prevalence))
                   predictions_LR<-predictions_func(LR,xvalidate)
                   LR_stats<-c(BrierScore(as.numeric(yvalidate)-1,predictions_LR),
                               auc(as.numeric(yvalidate),predictions_LR)[1],
                               as.numeric(sum(yvalidate==1 & predictions_LR>0.5)),
                               as.numeric(sum(yvalidate==0 & predictions_LR<=0.5)),
                               as.numeric(sum(yvalidate==0 & predictions_LR>0.5)),
                               as.numeric(sum(yvalidate==1 & predictions_LR<=0.5)),
                               as.numeric(sum(yvalidate==1 & predictions_LR>threshold_LR)),
                               as.numeric(sum(yvalidate==0 & predictions_LR<=threshold_LR)),
                               as.numeric(sum(yvalidate==0 & predictions_LR>threshold_LR)),
                               as.numeric(sum(yvalidate==1 & predictions_LR<=threshold_LR)),
                               as.numeric(sum(yvalidate==1 & predictions_LR>prevalence)),
                               as.numeric(sum(yvalidate==0 & predictions_LR<=prevalence)),
                               as.numeric(sum(yvalidate==0 & predictions_LR>prevalence)),
                               as.numeric(sum(yvalidate==1 & predictions_LR<=prevalence)),
                               as.numeric(sum(yvalidate==1 & predictions_LR>balanced_threshold_LR)),
                               as.numeric(sum(yvalidate==0 & predictions_LR<=balanced_threshold_LR)),
                               as.numeric(sum(yvalidate==0 & predictions_LR>balanced_threshold_LR)),
                               as.numeric(sum(yvalidate==1 & predictions_LR<=balanced_threshold_LR)),
                               threshold_LR,balanced_threshold_LR)
                   rm(RF1,predictions_LR,threshold_LR,predictions_LR_train,balanced_threshold_LR)
                   
                   ### Random Forest 1
                   RF1<-ranger(outcome ~ ., train_x, mtry = floor(sqrt(ncol(train_x)-1)),probability = TRUE)
                   # save the assessment of the predictions
                   predictions_RF1_train<-predictions_func(RF1,train_x[,-ncol(train_x)])
                   threshold_RF1<-invisible(thresholder(predictions_RF1_train))
                   balanced_threshold_RF1<-mean(c(threshold_RF1,prevalence))
                   predictions_RF1<-predictions_func(RF1,xvalidate)
                   RF1_stats<-c(BrierScore(as.numeric(yvalidate)-1,predictions_RF1),
                                auc(as.numeric(yvalidate),predictions_RF1)[1],
                                as.numeric(sum(yvalidate==1 & predictions_RF1>0.5)),
                                as.numeric(sum(yvalidate==0 & predictions_RF1<=0.5)),
                                as.numeric(sum(yvalidate==0 & predictions_RF1>0.5)),
                                as.numeric(sum(yvalidate==1 & predictions_RF1<=0.5)),
                                as.numeric(sum(yvalidate==1 & predictions_RF1>threshold_RF1)),
                                as.numeric(sum(yvalidate==0 & predictions_RF1<=threshold_RF1)),
                                as.numeric(sum(yvalidate==0 & predictions_RF1>threshold_RF1)),
                                as.numeric(sum(yvalidate==1 & predictions_RF1<=threshold_RF1)),
                                as.numeric(sum(yvalidate==1 & predictions_RF1>prevalence)),
                                as.numeric(sum(yvalidate==0 & predictions_RF1<=prevalence)),
                                as.numeric(sum(yvalidate==0 & predictions_RF1>prevalence)),
                                as.numeric(sum(yvalidate==1 & predictions_RF1<=prevalence)),
                                as.numeric(sum(yvalidate==1 & predictions_RF1>balanced_threshold_RF1)),
                                as.numeric(sum(yvalidate==0 & predictions_RF1<=balanced_threshold_RF1)),
                                as.numeric(sum(yvalidate==0 & predictions_RF1>balanced_threshold_RF1)),
                                as.numeric(sum(yvalidate==1 & predictions_RF1<=balanced_threshold_RF1)),
                                threshold_RF1,balanced_threshold_RF1)
                   rm(RF1,predictions_RF1,threshold_RF1,predictions_RF1_train,balanced_threshold_RF1)
                   
                   ### Random Forest 2
                   RF2<-ranger(outcome ~ ., train_x, mtry = floor(2*sqrt(ncol(train_x)-1)),probability = TRUE)
                   # save the assessment of the predictions
                   predictions_RF2_train<-predictions_func(RF2,train_x[,-ncol(train_x)])
                   threshold_RF2<-invisible(thresholder(predictions_RF2_train))
                   balanced_threshold_RF2<-mean(c(threshold_RF2,prevalence))
                   predictions_RF2<-predictions_func(RF2,xvalidate)
                   RF2_stats<-c(BrierScore(as.numeric(yvalidate)-1,predictions_RF2),
                                auc(as.numeric(yvalidate),predictions_RF2)[1],
                                as.numeric(sum(yvalidate==1 & predictions_RF2>0.5)),
                                as.numeric(sum(yvalidate==0 & predictions_RF2<=0.5)),
                                as.numeric(sum(yvalidate==0 & predictions_RF2>0.5)),
                                as.numeric(sum(yvalidate==1 & predictions_RF2<=0.5)),
                                as.numeric(sum(yvalidate==1 & predictions_RF2>threshold_RF2)),
                                as.numeric(sum(yvalidate==0 & predictions_RF2<=threshold_RF2)),
                                as.numeric(sum(yvalidate==0 & predictions_RF2>threshold_RF2)),
                                as.numeric(sum(yvalidate==1 & predictions_RF2<=threshold_RF2)),
                                as.numeric(sum(yvalidate==1 & predictions_RF2>prevalence)),
                                as.numeric(sum(yvalidate==0 & predictions_RF2<=prevalence)),
                                as.numeric(sum(yvalidate==0 & predictions_RF2>prevalence)),
                                as.numeric(sum(yvalidate==1 & predictions_RF2<=prevalence)),
                                as.numeric(sum(yvalidate==1 & predictions_RF2>balanced_threshold_RF2)),
                                as.numeric(sum(yvalidate==0 & predictions_RF2<=balanced_threshold_RF2)),
                                as.numeric(sum(yvalidate==0 & predictions_RF2>balanced_threshold_RF2)),
                                as.numeric(sum(yvalidate==1 & predictions_RF2<=balanced_threshold_RF2)),
                                threshold_RF2,balanced_threshold_RF2)
                   rm(RF2,predictions_RF2,threshold_RF2,predictions_RF2_train,balanced_threshold_RF2)
                   
                   ### Random Forest 3
                   RF3<-ranger(outcome ~ ., train_x, mtry = floor(4*sqrt(ncol(train_x)-1)),probability = TRUE)
                   # save the assessment of the predictions
                   predictions_RF3_train<-predictions_func(RF3,train_x[,-ncol(train_x)])
                   threshold_RF3<-invisible(thresholder(predictions_RF3_train))
                   balanced_threshold_RF3<-mean(c(threshold_RF3,prevalence))
                   predictions_RF3<-predictions_func(RF3,xvalidate)
                   RF3_stats<-c(BrierScore(as.numeric(yvalidate)-1,predictions_RF3),
                                auc(as.numeric(yvalidate),predictions_RF3)[1],
                                as.numeric(sum(yvalidate==1 & predictions_RF3>0.5)),
                                as.numeric(sum(yvalidate==0 & predictions_RF3<=0.5)),
                                as.numeric(sum(yvalidate==0 & predictions_RF3>0.5)),
                                as.numeric(sum(yvalidate==1 & predictions_RF3<=0.5)),
                                as.numeric(sum(yvalidate==1 & predictions_RF3>threshold_RF3)),
                                as.numeric(sum(yvalidate==0 & predictions_RF3<=threshold_RF3)),
                                as.numeric(sum(yvalidate==0 & predictions_RF3>threshold_RF3)),
                                as.numeric(sum(yvalidate==1 & predictions_RF3<=threshold_RF3)),
                                as.numeric(sum(yvalidate==1 & predictions_RF3>prevalence)),
                                as.numeric(sum(yvalidate==0 & predictions_RF3<=prevalence)),
                                as.numeric(sum(yvalidate==0 & predictions_RF3>prevalence)),
                                as.numeric(sum(yvalidate==1 & predictions_RF3<=prevalence)),
                                as.numeric(sum(yvalidate==1 & predictions_RF3>balanced_threshold_RF3)),
                                as.numeric(sum(yvalidate==0 & predictions_RF3<=balanced_threshold_RF3)),
                                as.numeric(sum(yvalidate==0 & predictions_RF3>balanced_threshold_RF3)),
                                as.numeric(sum(yvalidate==1 & predictions_RF3<=balanced_threshold_RF3)),
                                threshold_RF3,balanced_threshold_RF3)
                   rm(RF3,predictions_RF3,threshold_RF3,predictions_RF3_train,balanced_threshold_RF3)
                   
                   ### Random Forest 4
                   RF4<-ranger(outcome ~ ., train_x, mtry = floor(8*sqrt(ncol(train_x)-1)),probability = TRUE)
                   # save the assessment of the predictions
                   predictions_RF4_train<-predictions_func(RF4,train_x[,-ncol(train_x)])
                   threshold_RF4<-invisible(thresholder(predictions_RF4_train))
                   balanced_threshold_RF4<-mean(c(threshold_RF4,prevalence))
                   predictions_RF4<-predictions_func(RF4,xvalidate)
                   RF4_stats<-c(BrierScore(as.numeric(yvalidate)-1,predictions_RF4),
                                auc(as.numeric(yvalidate),predictions_RF4)[1],
                                as.numeric(sum(yvalidate==1 & predictions_RF4>0.5)),
                                as.numeric(sum(yvalidate==0 & predictions_RF4<=0.5)),
                                as.numeric(sum(yvalidate==0 & predictions_RF4>0.5)),
                                as.numeric(sum(yvalidate==1 & predictions_RF4<=0.5)),
                                as.numeric(sum(yvalidate==1 & predictions_RF4>threshold_RF4)),
                                as.numeric(sum(yvalidate==0 & predictions_RF4<=threshold_RF4)),
                                as.numeric(sum(yvalidate==0 & predictions_RF4>threshold_RF4)),
                                as.numeric(sum(yvalidate==1 & predictions_RF4<=threshold_RF4)),
                                as.numeric(sum(yvalidate==1 & predictions_RF4>prevalence)),
                                as.numeric(sum(yvalidate==0 & predictions_RF4<=prevalence)),
                                as.numeric(sum(yvalidate==0 & predictions_RF4>prevalence)),
                                as.numeric(sum(yvalidate==1 & predictions_RF4<=prevalence)),
                                as.numeric(sum(yvalidate==1 & predictions_RF4>balanced_threshold_RF4)),
                                as.numeric(sum(yvalidate==0 & predictions_RF4<=balanced_threshold_RF4)),
                                as.numeric(sum(yvalidate==0 & predictions_RF4>balanced_threshold_RF4)),
                                as.numeric(sum(yvalidate==1 & predictions_RF4<=balanced_threshold_RF4)),
                                threshold_RF4,balanced_threshold_RF4)
                   rm(RF4,predictions_RF4,threshold_RF4,predictions_RF4_train,balanced_threshold_RF4)
                   
                   ### XGBoost1
                   XGB1<-xgboost(xgb.DMatrix(data=as.matrix(train_x %>% select(-outcome)),
                                             label=as.integer(train_x$outcome)-1), 
                                 nrounds=100,
                                 eta=0.1,
                                 objective="binary:logistic",
                                 verbose=0)
                   # save the assessment of the predictions
                   predictions_XGB1_train<-predictions_func(XGB1,train_x[,-ncol(train_x)])
                   threshold_XGB1<-invisible(thresholder(predictions_XGB1_train))
                   balanced_threshold_XGB1<-mean(c(threshold_XGB1,prevalence))
                   predictions_XGB1<-predictions_func(XGB1,xvalidate)
                   XGB1_stats<-c(BrierScore(as.numeric(yvalidate)-1,predictions_XGB1),
                                 auc(as.numeric(yvalidate),predictions_XGB1)[1],
                                as.numeric(sum(yvalidate==1 & predictions_XGB1>0.5)),
                                as.numeric(sum(yvalidate==0 & predictions_XGB1<=0.5)),
                                as.numeric(sum(yvalidate==0 & predictions_XGB1>0.5)),
                                as.numeric(sum(yvalidate==1 & predictions_XGB1<=0.5)),
                                as.numeric(sum(yvalidate==1 & predictions_XGB1>threshold_XGB1)),
                                as.numeric(sum(yvalidate==0 & predictions_XGB1<=threshold_XGB1)),
                                as.numeric(sum(yvalidate==0 & predictions_XGB1>threshold_XGB1)),
                                as.numeric(sum(yvalidate==1 & predictions_XGB1<=threshold_XGB1)),
                                as.numeric(sum(yvalidate==1 & predictions_XGB1>prevalence)),
                                as.numeric(sum(yvalidate==0 & predictions_XGB1<=prevalence)),
                                as.numeric(sum(yvalidate==0 & predictions_XGB1>prevalence)),
                                as.numeric(sum(yvalidate==1 & predictions_XGB1<=prevalence)),
                                as.numeric(sum(yvalidate==1 & predictions_XGB1>balanced_threshold_XGB1)),
                                as.numeric(sum(yvalidate==0 & predictions_XGB1<=balanced_threshold_XGB1)),
                                as.numeric(sum(yvalidate==0 & predictions_XGB1>balanced_threshold_XGB1)),
                                as.numeric(sum(yvalidate==1 & predictions_XGB1<=balanced_threshold_XGB1)),
                                threshold_XGB1,balanced_threshold_XGB1)
                   rm(XGB1,predictions_XGB1,threshold_XGB1,predictions_XGB1_train,balanced_threshold_XGB1)
                   
                   ### XGBoost2
                   XGB2<-xgboost(xgb.DMatrix(data=as.matrix(train_x %>% select(-outcome)),
                                             label=as.integer(train_x$outcome)-1), 
                                 nrounds=100,
                                 eta=0.25,
                                 objective="binary:logistic",
                                 verbose=0)
                   # save the assessment of the predictions
                   predictions_XGB2_train<-predictions_func(XGB2,train_x[,-ncol(train_x)])
                   threshold_XGB2<-invisible(thresholder(predictions_XGB2_train))
                   balanced_threshold_XGB2<-mean(c(threshold_XGB2,prevalence))
                   predictions_XGB2<-predictions_func(XGB2,xvalidate)
                   XGB2_stats<-c(BrierScore(as.numeric(yvalidate)-1,predictions_XGB2),
                                 auc(as.numeric(yvalidate),predictions_XGB2)[1],
                                 as.numeric(sum(yvalidate==1 & predictions_XGB2>0.5)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB2<=0.5)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB2>0.5)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB2<=0.5)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB2>threshold_XGB2)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB2<=threshold_XGB2)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB2>threshold_XGB2)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB2<=threshold_XGB2)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB2>prevalence)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB2<=prevalence)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB2>prevalence)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB2<=prevalence)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB2>balanced_threshold_XGB2)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB2<=balanced_threshold_XGB2)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB2>balanced_threshold_XGB2)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB2<=balanced_threshold_XGB2)),
                                 threshold_XGB2,balanced_threshold_XGB2)
                   rm(XGB2,predictions_XGB2,threshold_XGB2,predictions_XGB2_train,balanced_threshold_XGB2)
                   
                   ### XGBoost3
                   XGB3<-xgboost(xgb.DMatrix(data=as.matrix(train_x %>% select(-outcome)),
                                             label=as.integer(train_x$outcome)-1), 
                                 nrounds=100,
                                 eta=0.5,
                                 objective="binary:logistic",
                                 verbose=0)
                   # save the assessment of the predictions
                   predictions_XGB3_train<-predictions_func(XGB3,train_x[,-ncol(train_x)])
                   threshold_XGB3<-invisible(thresholder(predictions_XGB3_train))
                   balanced_threshold_XGB3<-mean(c(threshold_XGB3,prevalence))
                   predictions_XGB3<-predictions_func(XGB3,xvalidate)
                   XGB3_stats<-c(BrierScore(as.numeric(yvalidate)-1,predictions_XGB3),
                                 auc(as.numeric(yvalidate),predictions_XGB3)[1],
                                 as.numeric(sum(yvalidate==1 & predictions_XGB3>0.5)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB3<=0.5)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB3>0.5)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB3<=0.5)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB3>threshold_XGB3)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB3<=threshold_XGB3)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB3>threshold_XGB3)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB3<=threshold_XGB3)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB3>prevalence)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB3<=prevalence)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB3>prevalence)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB3<=prevalence)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB3>balanced_threshold_XGB3)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB3<=balanced_threshold_XGB3)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB3>balanced_threshold_XGB3)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB3<=balanced_threshold_XGB3)),
                                 threshold_XGB3,balanced_threshold_XGB3)
                   rm(XGB3,predictions_XGB3,threshold_XGB3,predictions_XGB3_train,balanced_threshold_XGB3)
                   
                   ### XGBoost4
                   XGB4<-xgboost(xgb.DMatrix(data=as.matrix(train_x %>% select(-outcome)),
                                             label=as.integer(train_x$outcome)-1), 
                                 nrounds=200,
                                 eta=0.1,
                                 objective="binary:logistic",
                                 verbose=0)
                   # save the assessment of the predictions
                   predictions_XGB4_train<-predictions_func(XGB4,train_x[,-ncol(train_x)])
                   threshold_XGB4<-invisible(thresholder(predictions_XGB4_train))
                   balanced_threshold_XGB4<-mean(c(threshold_XGB4,prevalence))
                   predictions_XGB4<-predictions_func(XGB4,xvalidate)
                   XGB4_stats<-c(BrierScore(as.numeric(yvalidate)-1,predictions_XGB4),
                                 auc(as.numeric(yvalidate),predictions_XGB4)[1],
                                 as.numeric(sum(yvalidate==1 & predictions_XGB4>0.5)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB4<=0.5)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB4>0.5)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB4<=0.5)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB4>threshold_XGB4)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB4<=threshold_XGB4)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB4>threshold_XGB4)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB4<=threshold_XGB4)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB4>prevalence)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB4<=prevalence)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB4>prevalence)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB4<=prevalence)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB4>balanced_threshold_XGB4)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB4<=balanced_threshold_XGB4)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB4>balanced_threshold_XGB4)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB4<=balanced_threshold_XGB4)),
                                 threshold_XGB4,balanced_threshold_XGB4)
                   rm(XGB4,predictions_XGB4,threshold_XGB4,predictions_XGB4_train,balanced_threshold_XGB4)
                   
                   ### XGBoost5
                   XGB5<-xgboost(xgb.DMatrix(data=as.matrix(train_x %>% select(-outcome)),
                                             label=as.integer(train_x$outcome)-1), 
                                 nrounds=200,
                                 eta=0.25,
                                 objective="binary:logistic",
                                 verbose=0)
                   # save the assessment of the predictions
                   predictions_XGB5_train<-predictions_func(XGB5,train_x[,-ncol(train_x)])
                   threshold_XGB5<-invisible(thresholder(predictions_XGB5_train))
                   balanced_threshold_XGB5<-mean(c(threshold_XGB5,prevalence))
                   predictions_XGB5<-predictions_func(XGB5,xvalidate)
                   XGB5_stats<-c(BrierScore(as.numeric(yvalidate)-1,predictions_XGB5),
                                 auc(as.numeric(yvalidate),predictions_XGB5)[1],
                                 as.numeric(sum(yvalidate==1 & predictions_XGB5>0.5)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB5<=0.5)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB5>0.5)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB5<=0.5)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB5>threshold_XGB5)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB5<=threshold_XGB5)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB5>threshold_XGB5)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB5<=threshold_XGB5)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB5>prevalence)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB5<=prevalence)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB5>prevalence)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB5<=prevalence)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB5>balanced_threshold_XGB5)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB5<=balanced_threshold_XGB5)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB5>balanced_threshold_XGB5)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB5<=balanced_threshold_XGB5)),
                                 threshold_XGB5,balanced_threshold_XGB5)
                   rm(XGB5,predictions_XGB5,threshold_XGB5,predictions_XGB5_train,balanced_threshold_XGB5)
                   
                   ### XGBoost6
                   XGB6<-xgboost(xgb.DMatrix(data=as.matrix(train_x %>% select(-outcome)),
                                             label=as.integer(train_x$outcome)-1), 
                                 nrounds=200,
                                 eta=0.5,
                                 objective="binary:logistic",
                                 verbose=0)
                   # save the assessment of the predictions
                   predictions_XGB6_train<-predictions_func(XGB6,train_x[,-ncol(train_x)])
                   threshold_XGB6<-invisible(thresholder(predictions_XGB6_train))
                   balanced_threshold_XGB6<-mean(c(threshold_XGB6,prevalence))
                   predictions_XGB6<-predictions_func(XGB6,xvalidate)
                   XGB6_stats<-c(BrierScore(as.numeric(yvalidate)-1,predictions_XGB6),
                                 auc(as.numeric(yvalidate),predictions_XGB6)[1],
                                 as.numeric(sum(yvalidate==1 & predictions_XGB6>0.5)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB6<=0.5)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB6>0.5)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB6<=0.5)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB6>threshold_XGB6)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB6<=threshold_XGB6)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB6>threshold_XGB6)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB6<=threshold_XGB6)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB6>prevalence)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB6<=prevalence)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB6>prevalence)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB6<=prevalence)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB6>balanced_threshold_XGB6)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB6<=balanced_threshold_XGB6)),
                                 as.numeric(sum(yvalidate==0 & predictions_XGB6>balanced_threshold_XGB6)),
                                 as.numeric(sum(yvalidate==1 & predictions_XGB6<=balanced_threshold_XGB6)),
                                 threshold_XGB6,balanced_threshold_XGB6)
                   rm(XGB6,predictions_XGB6,threshold_XGB6,predictions_XGB6_train,balanced_threshold_XGB6)
                   
                   # Record sample sizes
                   c(enrichment,
                     i,
                     nrow(train_x[which(train_x$outcome==0),]),
                     nrow(train_x[which(train_x$outcome==1),]),
                     nbc_stats,
                     LR_stats,
                     RF1_stats,
                     RF2_stats,
                     RF3_stats,
                     RF4_stats,
                     XGB1_stats,
                     XGB2_stats,
                     XGB3_stats,
                     XGB4_stats,
                     XGB5_stats,
                     XGB6_stats)
                 }

stopCluster(cl)
Sys.time()-start
rm(train_data,indices,cl,start)

results<-data.frame(results)
names(results)<-c("enrichment","iteration","class_0","class_1",
                  "BS_nbc","AUC_nbc","TP_nbc","TN_nbc","FP_nbc","FN_nbc",
                  "TP_nbc5","TN_nbc5","FP_nbc5","FN_nbc5",
                  "TP_nbcp","TN_nbcp","FP_nbcp","FN_nbcp",
                  "TP_nbcb","TN_nbcb","FP_nbcb","FN_nbcb",
                  "thr_nbc5","thr_nbcb",
                  "BS_GLM","AUC_GLM","TP_GLM","TN_GLM","FP_GLM","FN_GLM",
                  "TP_GLM5","TN_GLM5","FP_GLM5","FN_GLM5",
                  "TP_GLMp","TN_GLMp","FP_GLMp","FN_GLMp",
                  "TP_GLMb","TN_GLMb","FP_GLMb","FN_GLMb",
                  "thr_GLM5","thr_GLMb",
                  "BS_RF1","AUC_RF1","TP_RF1","TN_RF1","FP_RF1","FN_RF1",
                  "TP_RF15","TN_RF15","FP_RF15","FN_RF15",
                  "TP_RF1p","TN_RF1p","FP_RF1p","FN_RF1p",
                  "TP_RF1b","TN_RF1b","FP_RF1b","FN_RF1b",
                  "thr_RF15","thr_RF1b",
                  "BS_RF2","AUC_RF2","TP_RF2","TN_RF2","FP_RF2","FN_RF2",
                  "TP_RF25","TN_RF25","FP_RF25","FN_RF25",
                  "TP_RF2p","TN_RF2p","FP_RF2p","FN_RF2p",
                  "TP_RF2b","TN_RF2b","FP_RF2b","FN_RF2b",
                  "thr_RF25","thr_RF2b",
                  "BS_RF3","AUC_RF3","TP_RF3","TN_RF3","FP_RF3","FN_RF3",
                  "TP_RF35","TN_RF35","FP_RF35","FN_RF35",
                  "TP_RF3p","TN_RF3p","FP_RF3p","FN_RF3p",
                  "TP_RF3b","TN_RF3b","FP_RF3b","FN_RF3b",
                  "thr_RF35","thr_RF3b",
                  "BS_RF4","AUC_RF4","TP_RF4","TN_RF4","FP_RF4","FN_RF4",
                  "TP_RF45","TN_RF45","FP_RF45","FN_RF45",
                  "TP_RF4p","TN_RF4p","FP_RF4p","FN_RF4p",
                  "TP_RF4b","TN_RF4b","FP_RF4b","FN_RF4b",
                  "thr_RF45","thr_RF4b",
                  "BS_XB1","AUC_XB1","TP_XB1","TN_XB1","FP_XB1","FN_XB1",
                  "TP_XB15","TN_XB15","FP_XB15","FN_XB15",
                  "TP_XB1p","TN_XB1p","FP_XB1p","FN_XB1p",
                  "TP_XB1b","TN_XB1b","FP_XB1b","FN_XB1b",
                  "thr_XB15","thr_XB1b",
                  "BS_XB2","AUC_XB2","TP_XB2","TN_XB2","FP_XB2","FN_XB2",
                  "TP_XB25","TN_XB25","FP_XB25","FN_XB25",
                  "TP_XB2p","TN_XB2p","FP_XB2p","FN_XB2p",
                  "TP_XB2b","TN_XB2b","FP_XB2b","FN_XB2b",
                  "thr_XB25","thr_XB2b",
                  "BS_XB3","AUC_XB3","TP_XB3","TN_XB3","FP_XB3","FN_XB3",
                  "TP_XB35","TN_XB35","FP_XB35","FN_XB35",
                  "TP_XB3p","TN_XB3p","FP_XB3p","FN_XB3p",
                  "TP_XB3b","TN_XB3b","FP_XB3b","FN_XB3b",
                  "thr_XB35","thr_XB3b",
                  "BS_XB4","AUC_XB4","TP_XB4","TN_XB4","FP_XB4","FN_XB4",
                  "TP_XB45","TN_XB45","FP_XB45","FN_XB45",
                  "TP_XB4p","TN_XB4p","FP_XB4p","FN_XB4p",
                  "TP_XB4b","TN_XB4b","FP_XB4b","FN_XB4b",
                  "thr_XB45","thr_XB4b",
                  "BS_XB5","AUC_XB5","TP_XB5","TN_XB5","FP_XB5","FN_XB5",
                  "TP_XB55","TN_XB55","FP_XB55","FN_XB55",
                  "TP_XB5p","TN_XB5p","FP_XB5p","FN_XB5p",
                  "TP_XB5b","TN_XB5b","FP_XB5b","FN_XB5b",
                  "thr_XB55","thr_XB5b",
                  "BS_XB6","AUC_XB6","TP_XB6","TN_XB6","FP_XB6","FN_XB6",
                  "TP_XB65","TN_XB65","FP_XB65","FN_XB65",
                  "TP_XB6p","TN_XB6p","FP_XB6p","FN_XB6p",
                  "TP_XB6b","TN_XB6b","FP_XB6b","FN_XB6b",
                  "thr_XB65","thr_XB6b")
rownames(results)<-NULL
write.csv(results,file="../Data_to_extract/enrichment_results v5.csv")


#---------------------------------------------------------------------------------------------------#
#
#  Results
#
#---------------------------------------------------------------------------------------------------#

results<-read_csv("../Data_to_extract/enrichment_results v5.csv", col_types = cols(...1 = col_skip()))
library(scales)

outcomes<-results %>%
  select(-class_0,-class_1) %>%
  gather("key","value",-enrichment,-iteration) %>%
  mutate(statistic = str_replace(substr(key,1,3),"_",""),
         modelx = str_replace(key,paste0(statistic,"_"),""),
         Threshold = ifelse(nchar(modelx)==3,"Fixed",
                            ifelse(substr(modelx,4,4)=="5","Variable",
                                   ifelse(substr(modelx,4,4)=="b","Balanced","Prevalence"))),
         model = substr(modelx,1,3),
         Algorithm = ifelse(model %in% c("GLM","nbc"),
                            toupper(model),
                            ifelse(substr(model,1,2)=="RF",
                                   "RF","XGB"))) %>%
  select(-key,-modelx) %>%
  spread(statistic,value) %>%
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

enrichment_class_sizes<-results %>%
  mutate(total = class_0+class_1,
         MCP = class_1*100/total) %>%
  select(enrichment,iteration,class_0,class_1,total,MCP) %>%
  gather(key = "class", value = "size", -iteration,-enrichment,-MCP)  %>%
  group_by(enrichment,class) %>%
  summarize(avg_size = round(mean(size),1),
            min_size = round(min(size),1),
            max_size=round(max(size),1))

# sample size by enrichment method
ggplot(enrichment_class_sizes %>% filter(class!="total")) +
  geom_bar(aes(x=as.factor(enrichment),y=avg_size,fill=class),
           position="stack",stat="identity") +
  scale_y_continuous(labels=comma, limits=c(0,800000),
                     expand=c(0,0),breaks=seq(0,800000,100000)) +
  theme_bw() +
  xlab("Enrichment Method") +
  ylab("Sample Size") +
  labs(fill="Class") +
  scale_fill_hue(labels=c("Positive","Negative"))

# average AUC across iterations for enrichment and models
ggplot(outcomes) +
  geom_boxplot(aes(x=as.factor(enrichment), y=AUC,fill=Algorithm)) +
  theme_bw() +
  xlab("Enrichment Method") +
  ylab("Area Under the Curve") +
  theme(text = element_text(size=15))

ggplot(outcomes) +
  geom_boxplot(aes(x=as.factor(enrichment), y=AUC,fill=as.factor(enrichment))) +
  theme_bw() +
  xlab("Enrichment Method") +
  ylab("Area Under the Curve") +
  facet_wrap(vars(Algorithm),nrow=1) +
  theme(text = element_text(size=15),
        legend.position = "none")

# then look at the other measures by enrichment and thresholds
ggplot(outcomes) +
  geom_boxplot(aes(x=as.factor(enrichment), y=MCC,fill=Threshold)) +
  theme_bw() +
  xlab("Enrichment Method") +
  ylab("MCC") +
  facet_grid(Algorithm~Threshold) +
  theme(text = element_text(size=15))

ggplot(outcomes) +
  geom_boxplot(aes(x=as.factor(enrichment), y=accuracy,fill=Threshold)) +
  theme_bw() +
  xlab("Enrichment Method") +
  ylab("Accuracy") +
  facet_grid(Algorithm~Threshold) +
  theme(text = element_text(size=15))

ggplot(outcomes) +
  geom_boxplot(aes(x=as.factor(enrichment), y=balanced_acc,fill=Threshold)) +
  theme_bw() +
  xlab("Enrichment Method") +
  ylab("Balanced Accuracy") +
  facet_grid(Algorithm~Threshold) +
  theme(text = element_text(size=15))

ggplot(outcomes) +
  geom_boxplot(aes(x=as.factor(enrichment), y=sensitivity,fill=Threshold)) +
  theme_bw() +
  xlab("Enrichment Method") +
  ylab("Sensitivity") +
  facet_grid(Algorithm~Threshold) +
  theme(text = element_text(size=15))

ggplot(outcomes) +
  geom_boxplot(aes(x=as.factor(enrichment), y=specificity,fill=Threshold)) +
  theme_bw() +
  xlab("Enrichment Method") +
  ylab("Specificity") +
  facet_grid(Algorithm~Threshold) +
  theme(text = element_text(size=15))

ggplot(outcomes) +
  geom_boxplot(aes(x=as.factor(enrichment), y=PPV,fill=Threshold)) +
  theme_bw() +
  xlab("Enrichment Method") +
  ylab("PPV") +
  facet_grid(Algorithm~Threshold) +
  theme(text = element_text(size=15))

ggplot(outcomes) +
  geom_boxplot(aes(x=as.factor(enrichment), y=NPV,fill=Threshold)) +
  theme_bw() +
  xlab("Enrichment Method") +
  ylab("NPV") +
  facet_grid(Algorithm~Threshold) +
  theme(text = element_text(size=15))

summary(lm(MCC ~ enrichment + Threshold + Algorithm, data=outcomes))
summary(lm(MCC ~ as.factor(enrichment)*Threshold, data=outcomes %>% filter(Algorithm=="GLM")))

#-----------------------------------------------------------------------------
# Thresholds
#-----------------------------------------------------------------------------

GLM<-outcomes %>% filter(Algorithm=="GLM" & enrichment==1)
  
ggplot(GLM) +
  geom_boxplot(aes(x=as.factor(Threshold), y=MCC)) +
  theme_bw() +
  xlab("Thresholding Approach") +
  ylab("MCC") +
  theme(text = element_text(size=15))

ggplot(GLM) +
  geom_boxplot(aes(x=as.factor(Threshold), y=sensitivity)) +
  theme_bw() +
  xlab("Thresholding Approach") +
  ylab("Sensitivity") +
  theme(text = element_text(size=15))

ggplot(GLM) +
  geom_boxplot(aes(x=as.factor(Threshold), y=PPV)) +
  theme_bw() +
  xlab("Thresholding Approach") +
  ylab("PPV") +
  theme(text = element_text(size=15))



temp<-outcomes %>% filter(model=="GLM" & enrichment==1 & !is.na(thr)) %>% select(Threshold,thr)
tapply(temp$thr,temp$Threshold,summary)
