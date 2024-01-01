library(glmnet)

##################################################################
## Load data
##################################################################
sample.info = readRDS("metadata_flow_survey.rds")
data.atac.pre = readRDS("prefiltered_atac_features.rds")


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

data.atac.pre = rbind(tmp,data.atac.pre[,sample.info$Study.LID..RRBS.])


##################################################################
## Get ready for model
##################################################################
sample_ids = names(data.atac.pre)
age = sample.info$Age..years.
names(age) = sample.info$Study.LID..RRBS.
age = age[sample_ids]

## Center/scale features
data.atac.pre = t(scale(t(data.atac.pre), center = T, scale = T))


##################################################################
## LOOCV glmnet
##################################################################
set.seed(123)
glm.loo.atac = c()
for(i in 1:length(age)){
  
  print(i)
  
  ## Remove held-out sample from training
  train <- data.atac.pre[,-i]
  test <- data.atac.pre[,i]
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
  features <- row.names(data.atac.pre)[index_beta]
  
  glm.loo.atac[[testid]] = list(glm = model, pred.age = predicted, actual.age = testage, coef = features)
  
}

## Quick view of predicted versus actual ages
pred = lapply(glm.loo.atac, function(x) x$pred.age)
actual = lapply(glm.loo.atac, function(x) x$actual.age)
plot(actual, pred)

## Save results for later
save(glm.loo.atac, file = "glm.loo.atac_alpha0.5_20230331.rda")
