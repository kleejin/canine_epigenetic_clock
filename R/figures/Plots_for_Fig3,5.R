## Examine model results

setwd("~/Desktop/dog_aging_clock/")
library(wesanderson)
library(ggplot2)
library(gridExtra)
library(dplyr)

##################################################################
## Predicated Age
##################################################################

## LOOCV1 results
load("data/glm.alls.atac_alpha0.5_20230331.rda")
load("data/glm.alls.meth_alpha0.5_20230331.rda")
load("data/glm.alls.comb_alpha0.5_20230331.rda")

## LOOCV2 results
load("data/glm.loo.atac_alpha0.5_20230331.rda")
load("data/glm.loo.meth_alpha0.5_20230331.rda")
load("data/glm.loo.comb_alpha0.5_20230331.rda")



## Extract info
loo.atac.results = data.frame(id = names(unlist(lapply(glm.loo.atac, function(x) x$pred.age))),
                              pred.age = unlist(lapply(glm.loo.atac, function(x) x$pred.age)),
                              actual.age = unlist(lapply(glm.loo.atac, function(x) x$actual.age)))
loo.atac.results$data = "ATAC"
loo.atac.results$model = "LOOCV2"

loo.meth.results = data.frame(id = names(unlist(lapply(glm.loo.meth, function(x) x$pred.age))),
                              pred.age = unlist(lapply(glm.loo.meth, function(x) x$pred.age)),
                              actual.age = unlist(lapply(glm.loo.meth, function(x) x$actual.age)))
loo.meth.results$data = "Methylation"
loo.meth.results$model = "LOOCV2"


loo.comb.results = data.frame(id = names(unlist(lapply(glm.loo.comb, function(x) x$pred.age))),
                              pred.age = unlist(lapply(glm.loo.comb, function(x) x$pred.age)),
                              actual.age = unlist(lapply(glm.loo.comb, function(x) x$actual.age)))
loo.comb.results$data = "Combined"
loo.comb.results$model = "LOOCV2"



alls.atac.results = data.frame(id = names(glm.alls.atac$pred.age[,1]),
                              pred.age = glm.alls.atac$pred.age[,1],
                              actual.age = glm.alls.atac$actual.age)
alls.atac.results$data = "ATAC"
alls.atac.results$model = "LOOCV1"


alls.meth.results = data.frame(id = names(glm.alls.meth$pred.age[,1]),
                               pred.age = glm.alls.meth$pred.age[,1],
                               actual.age = glm.alls.meth$actual.age)
alls.meth.results$data = "Methylation"
alls.meth.results$model = "LOOCV1"


alls.comb.results = data.frame(id = names(glm.alls.comb$pred.age[,1]),
                               pred.age = glm.alls.comb$pred.age[,1],
                               actual.age = glm.alls.comb$actual.age)
alls.comb.results$data = "Combined"
alls.comb.results$model = "LOOCV1"


atac.results = rbind(loo.atac.results, alls.atac.results)
meth.results = rbind(loo.meth.results, alls.meth.results)
comb.results = rbind(loo.comb.results, alls.comb.results)


all.results = rbind(atac.results, meth.results, comb.results)
all.results$data = factor(all.results$data, levels = c("ATAC", "Methylation", "Combined"))
sample.keep = alls.comb.results$id

## Plot only non-normalized data
pal <- wes_palette("FantasticFox1", 3, type = "continuous")
pal
p = all.results %>%
  ggplot(aes(y = pred.age, x = actual.age, color = data))+
  geom_point(alpha = 0.5, pch = 16)+
  scale_color_manual(values = pal)+
  facet_grid(vars(model), vars(data))+
  geom_abline(intercept = 0, slope = 1, lty = 2, color = "grey")+
  stat_smooth(method="lm", se=T, lwd = 0.5)+
  theme_bw()
ggsave(p, filename = "dogclock_alpha0.5.pdf", width = 6.75, height = 3.2, useDingbats = F)


## Add weight
load("data/metadata_flow_survey.rda")

all.results = left_join(all.results, sample.info[,c("Study.LID..RRBS.", "final.weight..kg.", "weight.cat", "weight.cat2")],
                        by = c("id" = "Study.LID..RRBS."))

p1 = all.results[all.results$model == "LOOCV1",] %>% 
  ggplot(aes(y = pred.age, x = actual.age, color = data))+
  # geom_point(aes(y = age), color = "black")+
  geom_point(alpha = 0.5, pch = 16)+
  scale_color_manual(values = pal)+
  # geom_point(data = age.only, aes(x = id, y = age), color = "black")+
  # facet_wrap(~weight.cat, nrow = 1, scales = "free_x")+
  geom_abline(intercept = 0, slope = 1, lty = 2, color = "grey")+
  facet_grid(vars(data), vars(weight.cat))+
  stat_smooth(method="lm", se=T, lwd = 0.5)+
  theme_bw()+
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))+
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12))+
  theme(legend.position="none")


p2 = all.results[all.results$model == "LOOCV2",] %>% ## <-- Use this one
  ggplot(aes(y = pred.age, x = actual.age, color = data))+
  # geom_point(aes(y = age), color = "black")+
  geom_point(alpha = 0.5, pch = 16)+
  scale_color_manual(values = pal)+
  # geom_point(data = age.only, aes(x = id, y = age), color = "black")+
  # facet_wrap(~weight.cat, nrow = 1, scales = "free_x")+
  geom_abline(intercept = 0, slope = 1, lty = 2, color = "grey")+
  facet_grid(vars(data), vars(weight.cat))+
  stat_smooth(method="lm", se=T, lwd = 0.5)+
  theme_bw()+
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))+
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12))+
  theme(legend.position="none")

grid.arrange(p1, p2, nrow = 1)

ggsave(p1, filename = 'dogclock_by_size_smaller_alls_alpha0.5.pdf', width = 5, height = 5, useDingbats = F)
ggsave(p2, filename = 'dogclock_by_size_smaller_loo_alpha0.5.pdf', width = 5, height = 5, useDingbats = F)






##################################################################
## Residual age
##################################################################
library(GGally)

## Residual Age correlated between models
tmp1 = residuals(lm(pred.age ~ actual.age, data = alls.atac.results))
tmp2 = residuals(lm(pred.age ~ actual.age, data = alls.meth.results))
tmp3 = residuals(lm(pred.age ~ actual.age, data = alls.comb.results))

residual.summary = data.frame(id = names(tmp1), ATAC = tmp1, Methylation = tmp2, Combined = tmp3)

lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "blue") +
    geom_smooth(method = method, color = "red", ...)
  p
}


p = ggpairs(residual.summary[,2:4],lower = list(continuous = wrap(lowerFn, method = "lm")),
        diag = list(continuous = wrap("barDiag", colour = "black")),
        upper = list(continuous = wrap("cor", size = 10)))
ggsave(p, filename = "ggpairs_residualage_alls_alpha0.5.pdf", width = 8, height = 6, useDingbats = F)





tmp1 = residuals(lm(pred.age ~ actual.age, data = loo.atac.results))
tmp2 = residuals(lm(pred.age ~ actual.age, data = loo.meth.results))
tmp3 = residuals(lm(pred.age ~ actual.age, data = loo.comb.results))

residual.summary = data.frame(id = names(tmp1), ATAC = tmp1, Methylation = tmp2, Combined = tmp3)

lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "blue") +
    geom_smooth(method = method, color = "red", ...)
  p
}


p = ggpairs(residual.summary[,2:4],lower = list(continuous = wrap(lowerFn, method = "lm")),
            diag = list(continuous = wrap("barDiag", colour = "black")),
            upper = list(continuous = wrap("cor", size = 10)))
ggsave(p, filename = "ggpairs_residualage_loo_alpha0.5.pdf", width = 8, height = 6, useDingbats = F)







## Residual age versus weight
tmp1 = residuals(lm(pred.age ~ actual.age, data = alls.atac.results))
tmp2 = residuals(lm(pred.age ~ actual.age, data = alls.meth.results))
tmp3 = residuals(lm(pred.age ~ actual.age, data = alls.comb.results))

all.results.loocv1 = all.results[all.results$model == "LOOCV1",]

all.results.loocv1$resid.age = NA

for(i in sample.keep){
  all.results.loocv1[all.results.loocv1$id == i & all.results.loocv1$data == "ATAC", "resid.age"] = tmp1[i]
  all.results.loocv1[all.results.loocv1$id == i & all.results.loocv1$data == "Methylation", "resid.age"] = tmp2[i]
  all.results.loocv1[all.results.loocv1$id == i & all.results.loocv1$data == "Combined", "resid.age"] = tmp3[i]
}

p = all.results.loocv1 %>%
  ggplot(aes(y = resid.age, x = final.weight..kg., color = data))+
  geom_point(alpha = 0.5, pch = 16)+
  facet_wrap(~data, nrow = 1)+
  scale_color_manual(values = pal)+
  stat_smooth(method="lm", se=T, lwd = 0.5)+
  xlab("Estimated Breed Weight (kg)")+
  ylab("Residual Age")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))+
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12))+
  theme(legend.position="none")

ggsave(p, filename = 'redisual_age_v_weight_alls_alpha0.5.pdf', width = 9, height = 3, useDingbats = F)





tmp1 = residuals(lm(pred.age ~ actual.age, data = loo.atac.results))
tmp2 = residuals(lm(pred.age ~ actual.age, data = loo.meth.results))
tmp3 = residuals(lm(pred.age ~ actual.age, data = loo.comb.results))
all.results.loocv2 = all.results[all.results$model == "LOOCV2",]

all.results.loocv2$resid.age = NA

for(i in sample.keep){
  all.results.loocv2[all.results.loocv2$id == i & all.results.loocv2$data == "ATAC", "resid.age"] = tmp1[i]
  all.results.loocv2[all.results.loocv2$id == i & all.results.loocv2$data == "Methylation", "resid.age"] = tmp2[i]
  all.results.loocv2[all.results.loocv2$id == i & all.results.loocv2$data == "Combined", "resid.age"] = tmp3[i]
}


p = all.results.loocv2 %>%
  ggplot(aes(y = resid.age, x = final.weight..kg., color = data))+
  geom_point(alpha = 0.5, pch = 16)+
  facet_wrap(~data, nrow = 1)+
  scale_color_manual(values = pal)+
  stat_smooth(method="lm", se=T, lwd = 0.5)+
  xlab("Estimated Breed Weight (kg)")+
  ylab("Residual Age")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 12), strip.text.y = element_text(size = 12))+
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12))+
  theme(legend.position="none")

ggsave(p, filename = 'redisual_age_v_weight_loo_alpha0.5.pdf', width = 9, height = 3, useDingbats = F)




## Correlated predicted age
predage.summary = data.frame(id = sample.keep)
row.names(predage.summary) = predage.summary$id
predage.summary$ATAC = loo.atac.results[sample.keep, "pred.age"]
predage.summary$Methylation = loo.meth.results[sample.keep, "pred.age"]
predage.summary$Combined = loo.comb.results[sample.keep, "pred.age"]

p1 = ggpairs(predage.summary[,2:4],lower = list(continuous = wrap(lowerFn, method = "lm")),
        diag = list(continuous = wrap("barDiag", colour = "blue")),
        upper = list(continuous = wrap("cor", size = 10)))


predage.summary = data.frame(id = sample.keep)
row.names(predage.summary) = predage.summary$id
predage.summary$ATAC = alls.atac.results[sample.keep, "pred.age"]
predage.summary$Methylation = alls.meth.results[sample.keep, "pred.age"]
predage.summary$Combined = alls.comb.results[sample.keep, "pred.age"]

p2 = ggpairs(predage.summary[,2:4],lower = list(continuous = wrap(lowerFn, method = "lm")),
             diag = list(continuous = wrap("barDiag", colour = "blue")),
             upper = list(continuous = wrap("cor", size = 10)))





ggsave(p1, filename = "ggpairs_predictedage_loo_alpha0.5.pdf", width = 8, height = 6, useDingbats = F)
ggsave(p2, filename = "ggpairs_predictedage_alls_alpha0.5.pdf", width = 8, height = 6, useDingbats = F)




##################################################################
## Linear models
##################################################################
## Age
atac.results = left_join(atac.results, sample.info, by = c("id" = "Study.LID..RRBS."))
meth.results = left_join(meth.results, sample.info, by = c("id" = "Study.LID..RRBS."))
comb.results = left_join(comb.results, sample.info, by = c("id" = "Study.LID..RRBS."))

## --> All data: allsamples
tmp = summary(lm(pred.age ~ actual.age, data = alls.atac.results))
tmp
tmp$coefficients
tmp = summary(lm(pred.age ~ actual.age, data = alls.meth.results))
tmp
tmp$coefficients
tmp = summary(lm(pred.age ~ actual.age, data = alls.comb.results))
tmp
tmp$coefficients

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

alls.atac.results = left_join(alls.atac.results, sample.info, by = c("id" = "Study.LID..RRBS."))
alls.meth.results = left_join(alls.meth.results, sample.info, by = c("id" = "Study.LID..RRBS."))
alls.comb.results = left_join(alls.comb.results, sample.info, by = c("id" = "Study.LID..RRBS."))

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
## RMSE
##################################################################
## LOOCV2
sqrt(mean((loo.atac.results$actual.age - loo.atac.results$pred.age)^2))
sqrt(mean((loo.meth.results$actual.age - loo.meth.results$pred.age)^2))
sqrt(mean((loo.comb.results$actual.age - loo.comb.results$pred.age)^2))

## LOOCV1
sqrt(mean((alls.atac.results$actual.age - alls.atac.results$pred.age)^2))
sqrt(mean((alls.meth.results$actual.age - alls.meth.results$pred.age)^2))
sqrt(mean((alls.comb.results$actual.age - alls.comb.results$pred.age)^2))

## By weight --LOOCV2
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



## By weight --LOOCV1
tmp.atac = alls.atac.results[alls.atac.results$weight.cat == "small",]
tmp.meth = alls.meth.results[alls.meth.results$weight.cat == "small",]
tmp.comb = alls.comb.results[alls.comb.results$weight.cat == "small",]
sqrt(mean((tmp.atac$actual.age - tmp.atac$pred.age)^2))
sqrt(mean((tmp.meth$actual.age - tmp.meth$pred.age)^2))
sqrt(mean((tmp.comb$actual.age - tmp.comb$pred.age)^2))

tmp.atac = alls.atac.results[alls.atac.results$weight.cat == "medium",]
tmp.meth = alls.meth.results[alls.meth.results$weight.cat == "medium",]
tmp.comb = alls.comb.results[alls.comb.results$weight.cat == "medium",]
sqrt(mean((tmp.atac$actual.age - tmp.atac$pred.age)^2))
sqrt(mean((tmp.meth$actual.age - tmp.meth$pred.age)^2))
sqrt(mean((tmp.comb$actual.age - tmp.comb$pred.age)^2))

tmp.atac = alls.atac.results[alls.atac.results$weight.cat == "large",]
tmp.meth = alls.meth.results[alls.meth.results$weight.cat == "large",]
tmp.comb = alls.comb.results[alls.comb.results$weight.cat == "large",]
sqrt(mean((tmp.atac$actual.age - tmp.atac$pred.age)^2))
sqrt(mean((tmp.meth$actual.age - tmp.meth$pred.age)^2))
sqrt(mean((tmp.comb$actual.age - tmp.comb$pred.age)^2))
