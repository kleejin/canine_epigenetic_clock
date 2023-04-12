options(stringsAsFactors = F)
setwd("/Volumes/kellyjin_externalhd/Dog ATAC/")

library(tidyverse)
library(gprofiler2)
library(GenomicRanges)
# library(amap)
library(RColorBrewer)
library(gridExtra)

#########################################################
## Load stuff
#########################################################
load("data/metadata_flow_survey.rda")
names(sample.info)
meta.feat = names(sample.info)[c(17, 19:49, 58:65, 71)]


#########################################################
## Meta data v Age
#########################################################
print(meta.feat)

## Flow features
flow.lm.results = data.frame(meta.feature = meta.feat[2:32], p = NA, effectsize = NA, adj.r2 = NA)
row.names(flow.lm.results) = flow.lm.results$meta.feature
for(i in meta.feat[2:32]){
  age = as.numeric(sample.info$Age..years.)
  feat = as.numeric(sample.info[,i])
  lm.tmp = summary(lm(feat ~ age))
  flow.lm.results[i,"p"] = lm.tmp$coefficients[2,4]
  flow.lm.results[i,"effectsize"] = lm.tmp$coefficients[2,1]
  flow.lm.results[i,"adj.r2"] = lm.tmp$adj.r.squared
}

## Survey features
survey.lm.results = data.frame(meta.feature = meta.feat[-c(2:32)], p = NA, f = NA)
row.names(survey.lm.results) = survey.lm.results$meta.feature
for(i in meta.feat[-c(2:32)]){
  age = as.numeric(sample.info$Age..years.)
  feat = sample.info[,i]
  lm.tmp = aov(age ~ feat)
  survey.lm.results[i,"p"] = summary(lm.tmp)[[1]][["Pr(>F)"]][1]
  survey.lm.results[i,"f"] = summary(lm.tmp)[[1]][["F value"]][1]
}

#########################################################
## Plots
#########################################################
## Flow
flow.lm.results$adj.p = p.adjust(flow.lm.results$p, method = "fdr")
signif.flow = flow.lm.results$meta.feature[flow.lm.results$adj.p < 0.05]
print(signif.flow)

g1 = sample.info %>%
  ggplot(aes(x = Age..years., y = CD62L.CD44Hi.CD8.T.cells))+
  scale_y_log10()+
  geom_point(alpha = 0.5, pch = 16)+
  xlab("Age (years)")+ylab("CD62L CD44Hi CD8 T cells")+
  stat_smooth(method="lm", se=T, lwd = 0.5)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=15))

g2 = sample.info %>%
  ggplot(aes(x = Age..years., y = CD62L.CD44Hi.DN.T.cells))+
  scale_y_log10()+
  geom_point(alpha = 0.5, pch = 16)+
  xlab("Age (years)")+ylab("CD62L CD44Hi DN T cells")+
  stat_smooth(method="lm", se=T, lwd = 0.5)+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=15))


## Survey
sample.info$exercise_type_short = gsub(" (.*)", "", sample.info$exercise_type)
g3 = sample.info %>%
  ggplot(aes(x = exercise_type_short, y = Age..years.))+
  geom_boxplot()+
  geom_jitter(width = 0.2)+
  xlab("Exercise Type")+ylab("Age (years)")+
  theme_bw()+
  theme(axis.text=element_text(size=14), axis.title=element_text(size=15))

pdf(file = "survey_flow_age_corr.pdf", width = 10, height = 8)
grid.arrange(g1, g2, g3, nrow = 2)
dev.off()

