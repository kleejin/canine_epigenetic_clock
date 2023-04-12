
library(glmnet)


setwd("~/Desktop/dog_aging_clock/")

##################################################################
## Load data
##################################################################
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


percMethData = rbind(tmp,filtered_perc_meth[,sample.info$Study.LID..RRBS.])



##################################################################
## Get ready for model
##################################################################
sample_ids = sample.info$Study.LID..RRBS.
percMethData = percMethData[,sample_ids]

# Center/scale model
percMethData = t(scale(t(percMethData), center = T, scale = T))



##################################################################
## LOOCV1
##################################################################
set.seed(123)

train = percMethData
trainage = age

final.model = cv.glmnet(data.matrix(t(train)), trainage, alpha = 0.5, nfolds=length(trainage), family="gaussian")
predicted = predict(final.model, newx=t(train), s=final.model$lambda.min)

plot(predicted ~ trainage)

betas = coef(final.model, s=final.model$lambda.min)
index_beta =  which(!betas[,1] == 0)
betacoefs = betas[index_beta]
features = row.names(percMethData)[index_beta]

## Store results
glm.alls.meth = list(glm = final.model, pred.age = predicted, actual.age = age, coef = features)

save(glm.alls.meth, file = "data/glm.alls.meth_alpha0.5_20230331.rda")
