library(glmnet)


setwd("~/Desktop/dog_aging_clock/")

##################################################################
## Load data
##################################################################
load("data/data_atac_filtered.Rdata") 
load("data/data_dnam_filtered.RData")
load("data/metadata_flow_survey.rda")

age = sample.info$Age..years.
names(age) = sample.info$Study.LID..RRBS.

##################################################################
## Add meta data features
##################################################################
names(sample.info)
tmp = sample.info[,c(17, 19:49, 71)]
row.names(tmp) = sample.info$Study.LID..RRBS.
tmp$sex2 = as.numeric(factor(tmp$sex2))
tmp$weight.cat = as.numeric(factor(tmp$weight.cat, levels = c("small", "medium", "large")))

tmp = data.frame(t(tmp))
tmp = tmp[,sample.info$Study.LID..RRBS.]

sample_ids = intersect(names(data.atac.pre), names(filtered_perc_meth))

data.all.0 = rbind(filtered_perc_meth[,sample_ids], data.atac.pre[,sample_ids], tmp[,sample_ids])

##################################################################
## Prepare for model
##################################################################
## Center data 
data.all.center = scale(data.all.0, center = T, scale = F)
data.all.cscale= scale(data.all.0, center = T, scale = T)



age = sample.info$Age..years.
names(age) = sample.info$Study.LID..RRBS.
age = age[sample_ids]


## Choose which dataset to use
data.all = data.all.cscale


##################################################################
## LOOCV2
##################################################################
set.seed(123)
glm.loo.comb = c()
for(i in 1:length(age)){
  
  print(i)
  
  ## Remove held-out sample from training
  train <- data.all[,-i]
  test <- data.all[,i]
  trainage <- age[-i]
  testage <- age[i]
  testid <- sample_ids[i]
  trainid <- sample_ids[-i]
  
  ## LOOCV
  model <- cv.glmnet(data.matrix(t(train)), trainage, alpha = 0.5, nfolds=length(trainage), family="gaussian")
  
  ## Predict age of the held-out sample
  predicted <- predict(model, newx=t(test), s=model$lambda.min)
  
  ## Extract non-zero beta coefficients and associated features
  betas <- coef(model, s=model$lambda.min)
  index_beta <-  which(!betas[,1] == 0)
  betacoefs <- betas[index_beta]
  features <- row.names(data.all)[index_beta]
  
  glm.loo.comb[[testid]] = list(glm = model, pred.age = predicted, actual.age = testage, coef = features)
  
}

pred = lapply(glm.loo.comb, function(x) x$pred.age)
actual = lapply(glm.loo.comb, function(x) x$actual.age)
plot(actual, pred)


save(glm.loo.comb, file = "data/glm.loo.comb_alpha0.5_20230331.rda")
