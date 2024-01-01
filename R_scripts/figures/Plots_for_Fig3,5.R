library(wesanderson)
library(ggplot2)
library(gridExtra)
library(dplyr)

##################################################################
## Load results
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
## Plot Fig 3a,c
##################################################################
## Fig3a
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


## Fig3c
all.results = left_join(all.results, sample.info[,c("Study.LID..RRBS.", "final.weight..kg.", "weight.cat", "weight.cat2")],
                        by = c("id" = "Study.LID..RRBS."))

p = all.results[all.results$model == "LOOCV2",] %>% ## <-- Use this one
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

ggsave(p, filename = 'dogclock_by_size_smaller_loo_alpha0.5.pdf', width = 5, height = 5, useDingbats = F)




##################################################################
## Plot Fig 3b and SFig 3a
##################################################################
lambda.dist.meth = unlist(lapply(glm.loo.meth, function(x) x$glm$lambda.min))
lambda.dist.atac = unlist(lapply(glm.loo.atac, function(x) x$glm$lambda.min))
lambda.dist.comb = unlist(lapply(glm.loo.comb, function(x) x$glm$lambda.min))

ncoef.dist.meth = unlist(lapply(glm.loo.meth, function(x) length(x$coef)))
ncoef.dist.atac = unlist(lapply(glm.loo.atac, function(x) length(x$coef)))
ncoef.dist.comb = unlist(lapply(glm.loo.comb, function(x) length(x$coef)))

lambda.df = data.frame(lambda = c(lambda.dist.meth, lambda.dist.atac, lambda.dist.comb),
                       ncoef = c(ncoef.dist.meth, ncoef.dist.atac, ncoef.dist.comb),
                       data = c(rep("DNAm", length(lambda.dist.meth)), rep("ATAC", length(lambda.dist.atac)), rep("zComb", length(lambda.dist.comb))))

pal <- wes_palette("FantasticFox1", 3, type = "continuous")
g1 = lambda.df %>%
  ggplot(aes(x = lambda, fill = data))+
  geom_density(alpha = 0.5, adjust = 1, color = NA)+
  geom_vline(xintercept = mean(lambda.dist.atac), color = pal[1])+
  geom_vline(xintercept = mean(lambda.dist.comb), color = pal[3])+
  geom_vline(xintercept = mean(lambda.dist.meth), color = pal[2])+
  scale_fill_manual(values = pal)+
  theme_bw()

g2 = lambda.df %>%
  ggplot(aes(x = ncoef, fill = data))+
  geom_density(alpha = 0.5, adjust = 1, color = NA)+
  geom_vline(xintercept = mean(ncoef.dist.atac), color = pal[1])+
  geom_vline(xintercept = mean(ncoef.dist.comb), color = pal[3])+
  geom_vline(xintercept = mean(ncoef.dist.meth), color = pal[2])+
  scale_fill_manual(values = pal)+
  theme_bw()

ggsave(g1, filename = "lambda_loocv2_distributions.pdf", width = 4.5, height = 3.5, useDingbats = F)
ggsave(g2, filename = "ncoef_loocv2_distributions.pdf", width = 4.5, height = 3.5, useDingbats = F)





##################################################################
## Plot Fig 5; SFig 3b
##################################################################
library(GGally)

## Supplemental Figure 3b: Residual Age correlated between models
tmp1 = residuals(lm(pred.age ~ actual.age, data = atac.results))
tmp2 = residuals(lm(pred.age ~ actual.age, data = meth.results))
tmp3 = residuals(lm(pred.age ~ actual.age, data = comb.results))

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


## Plot Fig 5
tmp1 = residuals(lm(pred.age ~ actual.age, data = atac.results))
tmp2 = residuals(lm(pred.age ~ actual.age, data = meth.results))
tmp3 = residuals(lm(pred.age ~ actual.age, data = comb.results))
all.results.loocv2 = all.results[all.results$model == "LOOCV2",]

all.results.loocv2$resid.age = NA

for(i in sample.keep){
  all.results.loocv2[all.results.loocv2$id == i & all.results.loocv2$data == "ATAC", "resid.age"] = tmp1[i]
  all.results.loocv2[all.results.loocv2$id == i & all.results.loocv2$data == "Methylation", "resid.age"] = tmp2[i]
  all.results.loocv2[all.results.loocv2$id == i & all.results.loocv2$data == "Combined", "resid.age"] = tmp3[i]
}
all.results.loocv2 = left_join(all.results.loocv2, sample.info, by = c("id" = "Study.LID..RRBS."))


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

