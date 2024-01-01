options(stringsAsFactors = F)
library(tidyverse)
library(gprofiler2)
library(GenomicRanges)
library(RColorBrewer)
library(gridExtra)

#########################################################
## Load objects
#########################################################
sample.info = readRDS("metadata_flow_survey.rds")
names(sample.info)
meta.feat = names(sample.info)[c(17, 19:49, 58:65, 71)] #select meta data features to correlate with age (includes flow features and survey questions)


#########################################################
## Run linear models of age versus meta data features
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


flow.lm.results$p.adj = p.adjust(flow.lm.results$p, method = "BH")
survey.lm.results$p.adj = p.adjust(survey.lm.results$p, method = "BH")

head(flow.lm.results) #includes pval, adj.p, effectsize, and adj.r2 of each lm test
head(survey.lm.results)


## Save results for later
saveRDS(flow.lm.results, file = "linearmodel_results_flow.rds")
saveRDS(survey.lm.results, file = "linearmodel_results_survey.rds")




