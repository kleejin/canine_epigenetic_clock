library(glmnet)

##################################################################
## Load data
##################################################################
sample.info = readRDS("metadata_flow_survey.rds")
filtered_perc_meth = readRDS("prefiltered_dnam_features.rds")



##################################################################
## Add meta data features to feature matrix
##################################################################
names(sample.info)
tmp = sample.info[,c(17, 19:49, 71)]
row.names(tmp) = sample.info$Study.LID..RRBS.
tmp$sex2 = as.numeric(factor(tmp$sex2))
tmp$weight.cat = as.numeric(factor(tmp$weight.cat, levels = c("small", "medium", "large")))

tmp = data.frame(t(tmp))
tmp = tmp[,sample.info$Study.LID..RRBS.]

percMethData = rbind(tmp,filtered_perc_meth[,sample.info$Study.LID..RRBS.])


##################################################################
## Get ready for model
##################################################################
sample_ids = sample.info$Study.LID..RRBS.
percMethData = percMethData[,sample_ids]
age = sample.info$Age..years.
names(age) = sample.info$Study.LID..RRBS.

# Center/scale features
percMethData = t(scale(t(percMethData), center = T, scale = T))


##################################################################
## LOOCV glmnet
##################################################################
set.seed(123)
glm.loo.meth = c()
for(i in 1:length(age)){
  
  print(i)
  
  ## Remove held-out sample from training
  train <- percMethData[,-i]
  test <- percMethData[,i]
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
  features <- row.names(percMethData)[index_beta]
  
  glm.loo.meth[[testid]] = list(glm = model, pred.age = predicted, actual.age = testage, coef = features)
  
}

## Quick view of predicted versus actual ages
pred = unlist(lapply(glm.loo.meth, function(x) x$pred.age))
actual = unlist(lapply(glm.loo.meth, function(x) x$actual.age))
plot(actual, pred)
summary(lm(pred~actual))

## Save results for later
save(glm.loo.meth, file = "glm.loo.meth_alpha0.5_20230331.rda")
