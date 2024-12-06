library(performanceEstimation)
library(tidyverse)

df<-iris %>% 
  rename(Feature1 = Sepal.Length,
         Feature2 = Sepal.Width)  %>%
  mutate(data = "Original",
         Class = factor(ifelse(Species=="setosa","Minor","Major"))) %>%
  select(Feature1, Feature2,Class)

ggplot(df) + geom_point(aes(x=Feature1, y = Feature2, colour = Class)) +
  theme_bw()


df_smote<-smote(Class ~ .,df,
      perc.over= 1, 
      perc.under = 0.1)
table(df$Class)
table(df_smote$Class)

ggplot(df) + geom_point(aes(x=Feature1, y = Feature2, colour = Class), alpha=0.2) +
  theme_bw()
ggplot(bind_rows(df,df_smote)) + geom_point(aes(x=Feature1, y = Feature2, colour = Class), alpha=0.2) +
  theme_bw()

df_2<-bind_rows(df, 
                df_smote %>% mutate(data = "SMOTEd")) %>%
  mutate(Class = ifelse(Class =="Major","Major",
                        ifelse(is.na(data),
                               "Minor - SMOTEd",
                               "Minor - Original"))) %>%
  pivot_longer(cols=Feature1:Feature2, names_to = "Feature",values_to = "Value")

ggplot(df_2) + 
  geom_density(aes(x=Value,fill=Class),alpha=0.2) +
  theme_bw() +
  facet_wrap(vars(Feature),scales="free_x")




