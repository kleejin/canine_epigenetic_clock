library(glmnet)
##################################################################
## Load data
##################################################################
sample.info = readRDS("metadata_flow_survey.rds")
filtered_perc_meth = readRDS("prefiltered_dnam_features.rds")
data.atac.pre = readRDS("prefiltered_atac_features.rds")



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
age = sample.info$Age..years.
names(age) = sample.info$Study.LID..RRBS.
age = age[sample_ids]

#Center and scale features
data.all = scale(data.all.0, center = T, scale = T)


##################################################################
## LOOCV glmnet
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
  ## predict age of the held-out sample
  predicted <- predict(model, newx=t(test), s=model$lambda.min)
  ## Extract non-zero beta coefficients and associated features
  betas <- coef(model, s=model$lambda.min)
  index_beta <-  which(!betas[,1] == 0)
  betacoefs <- betas[index_beta]
  features <- row.names(data.all)[index_beta]
  
  glm.loo.comb[[testid]] = list(glm = model, pred.age = predicted, actual.age = testage, coef = features)

}

## Quick view of predicted versus actual ages
pred = lapply(glm.loo.comb, function(x) x$pred.age)
actual = lapply(glm.loo.comb, function(x) x$actual.age)
plot(actual, pred) #quick view of predicted versus actual ages

## Save results for later
save(glm.loo.comb, file = "glm.loo.comb_alpha0.5_20230331.rda")
