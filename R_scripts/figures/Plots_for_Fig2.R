options(stringsAsFactors = F)
library(tidyverse)
library(gprofiler2)
library(RColorBrewer)
library(gridExtra)

#########################################################
## Load 
#########################################################
sample.info = readRDS("metadata_flow_survey.rds")
names(sample.info)


#########################################################
## Plot Fig 2
#########################################################
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

