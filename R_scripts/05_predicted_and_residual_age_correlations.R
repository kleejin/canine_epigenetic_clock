library(wesanderson)
library(ggplot2)
library(gridExtra)
library(dplyr)

##################################################################
## Load predicted age results
##################################################################
load("glm.loo.atac_alpha0.5_20230331.rda")
load("glm.loo.meth_alpha0.5_20230331.rda")
load("glm.loo.comb_alpha0.5_20230331.rda")
sample.info = readRDS("metadata_flow_survey.rds")


## Extract info
atac.results = data.frame(id = names(unlist(lapply(glm.loo.atac, function(x) x$pred.age))),
                              pred.age = unlist(lapply(glm.loo.atac, function(x) x$pred.age)),
                              actual.age = unlist(lapply(glm.loo.atac, function(x) x$actual.age)))
atac.results$data = "ATAC"
atac.results$model = "LOOCV2"

meth.results = data.frame(id = names(unlist(lapply(glm.loo.meth, function(x) x$pred.age))),
                              pred.age = unlist(lapply(glm.loo.meth, function(x) x$pred.age)),
                              actual.age = unlist(lapply(glm.loo.meth, function(x) x$actual.age)))
meth.results$data = "Methylation"
meth.results$model = "LOOCV2"


comb.results = data.frame(id = names(unlist(lapply(glm.loo.comb, function(x) x$pred.age))),
                              pred.age = unlist(lapply(glm.loo.comb, function(x) x$pred.age)),
                              actual.age = unlist(lapply(glm.loo.comb, function(x) x$actual.age)))
comb.results$data = "Combined"
comb.results$model = "LOOCV2"

## Combine all results together
all.results = rbind(atac.results, meth.results, comb.results)
all.results$data = factor(all.results$data, levels = c("ATAC", "Methylation", "Combined"))
sample.keep = comb.results$id



##################################################################
## Residual age calculation and plotting
##################################################################
library(GGally)
pal <- wes_palette("FantasticFox1", 3, type = "continuous")

## Residual Age correlated between models
tmp1 = residuals(lm(pred.age ~ actual.age, data = atac.results))
tmp2 = residuals(lm(pred.age ~ actual.age, data = meth.results))
tmp3 = residuals(lm(pred.age ~ actual.age, data = comb.results))

residual.summary = data.frame(id = names(tmp1), ATAC = tmp1, Methylation = tmp2, Combined = tmp3)



##################################################################
## Linear models from Fig 3a,c; Fig 5
##################################################################
## Age
atac.results = left_join(atac.results, sample.info, by = c("id" = "Study.LID..RRBS."))
meth.results = left_join(meth.results, sample.info, by = c("id" = "Study.LID..RRBS."))
comb.results = left_join(comb.results, sample.info, by = c("id" = "Study.LID..RRBS."))

## --> All data: LOOCV
tmp = summary(lm(pred.age ~ actual.age, data = loo.atac.results))
tmp
tmp$coefficients
tmp = summary(lm(pred.age ~ actual.age, data = loo.meth.results))
tmp
tmp$coefficients
tmp = summary(lm(pred.age ~ actual.age, data = loo.comb.results))
tmp
tmp$coefficients

loo.atac.results = left_join(loo.atac.results, sample.info, by = c("id" = "Study.LID..RRBS."))
loo.meth.results = left_join(loo.meth.results, sample.info, by = c("id" = "Study.LID..RRBS."))
loo.comb.results = left_join(loo.comb.results, sample.info, by = c("id" = "Study.LID..RRBS."))


### LOOCV
## --> Large
summary(lm(pred.age ~ actual.age, data = loo.atac.results[loo.atac.results$weight.cat == "large",]))
summary(lm(pred.age ~ actual.age, data = loo.meth.results[loo.meth.results$weight.cat == "large",]))
summary(lm(pred.age ~ actual.age, data = loo.comb.results[loo.comb.results$weight.cat == "large",]))

## --> Med
summary(lm(pred.age ~ actual.age, data = loo.atac.results[loo.atac.results$weight.cat == "medium",]))
summary(lm(pred.age ~ actual.age, data = loo.meth.results[loo.meth.results$weight.cat == "medium",]))
summary(lm(pred.age ~ actual.age, data = loo.comb.results[loo.comb.results$weight.cat == "medium",]))

## --> Small
summary(lm(pred.age ~ actual.age, data = loo.atac.results[loo.atac.results$weight.cat == "small",]))
summary(lm(pred.age ~ actual.age, data = loo.meth.results[loo.meth.results$weight.cat == "small",]))
summary(lm(pred.age ~ actual.age, data = loo.comb.results[loo.comb.results$weight.cat == "small",]))


### Allsamples
## --> Large
summary(lm(pred.age ~ actual.age, data = alls.atac.results[alls.atac.results$weight.cat == "large",]))
summary(lm(pred.age ~ actual.age, data = alls.meth.results[alls.meth.results$weight.cat == "large",]))
summary(lm(pred.age ~ actual.age, data = alls.comb.results[alls.comb.results$weight.cat == "large",]))

## --> Med
summary(lm(pred.age ~ actual.age, data = alls.atac.results[alls.atac.results$weight.cat == "medium",]))
summary(lm(pred.age ~ actual.age, data = alls.meth.results[alls.meth.results$weight.cat == "medium",]))
summary(lm(pred.age ~ actual.age, data = alls.comb.results[alls.comb.results$weight.cat == "medium",]))

## --> Small
summary(lm(pred.age ~ actual.age, data = alls.atac.results[alls.atac.results$weight.cat == "small",]))
summary(lm(pred.age ~ actual.age, data = alls.meth.results[alls.meth.results$weight.cat == "small",]))
summary(lm(pred.age ~ actual.age, data = alls.comb.results[alls.comb.results$weight.cat == "small",]))


## Residual age
residual.summary = left_join(residual.summary, sample.info[,c("Study.LID..RRBS.", "final.weight..kg.", "weight.cat", "weight.cat2")],
                             by = c("id" = "Study.LID..RRBS."))
summary(lm(ATAC ~ final.weight..kg., data = residual.summary))
summary(lm(Methylation ~ final.weight..kg., data = residual.summary))
summary(lm(Combined ~ final.weight..kg., data = residual.summary))



##################################################################
## RMSE calculations from Fig3a,c
##################################################################
## LOOCV
sqrt(mean((loo.atac.results$actual.age - loo.atac.results$pred.age)^2))
sqrt(mean((loo.meth.results$actual.age - loo.meth.results$pred.age)^2))
sqrt(mean((loo.comb.results$actual.age - loo.comb.results$pred.age)^2))

## By weight --LOOCV
tmp.atac = loo.atac.results[loo.atac.results$weight.cat == "small",]
tmp.meth = loo.meth.results[loo.meth.results$weight.cat == "small",]
tmp.comb = loo.comb.results[loo.comb.results$weight.cat == "small",]
sqrt(mean((tmp.atac$actual.age - tmp.atac$pred.age)^2))
sqrt(mean((tmp.meth$actual.age - tmp.meth$pred.age)^2))
sqrt(mean((tmp.comb$actual.age - tmp.comb$pred.age)^2))

tmp.atac = loo.atac.results[loo.atac.results$weight.cat == "medium",]
tmp.meth = loo.meth.results[loo.meth.results$weight.cat == "medium",]
tmp.comb = loo.comb.results[loo.comb.results$weight.cat == "medium",]
sqrt(mean((tmp.atac$actual.age - tmp.atac$pred.age)^2))
sqrt(mean((tmp.meth$actual.age - tmp.meth$pred.age)^2))
sqrt(mean((tmp.comb$actual.age - tmp.comb$pred.age)^2))

tmp.atac = loo.atac.results[loo.atac.results$weight.cat == "large",]
tmp.meth = loo.meth.results[loo.meth.results$weight.cat == "large",]
tmp.comb = loo.comb.results[loo.comb.results$weight.cat == "large",]
sqrt(mean((tmp.atac$actual.age - tmp.atac$pred.age)^2))
sqrt(mean((tmp.meth$actual.age - tmp.meth$pred.age)^2))
sqrt(mean((tmp.comb$actual.age - tmp.comb$pred.age)^2))

