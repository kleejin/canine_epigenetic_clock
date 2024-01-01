#########################################################
## Plot Fig 1
#########################################################
sample.info = readRDS("metadata_flow_survey.rds")

par(mfrow = c(2,2))
#Age distr
hist(sample.info$Age..years., xlab = "Age (years)", breaks = 10,
     col = "grey", las = 1, main = "", cex.lab = 1.4, cex.axis = 1.25)
#Weight distr
hist(sample.info$final.weight..kg., xlab = "Est. Breed Weight (kg)", breaks = 10,
     col = "grey", las = 1, main = "", cex.lab = 1.4, cex.axis = 1.25)

#Weight by Age
plot(sample.info$final.weight..kg., sample.info$Age..years., xlab = "Est. Breed Weight (kg)", ylab = "Age (years)",
     las = 1, pch = 19, cex.lab = 1.4, cex.axis = 1.25, col = rgb(0,0,0,0.3))
#Top 10 Breeds
graphics.off()

pdf("breed_summary.pdf", width = 6, height = 6)
par(mar = c(8,10,2,2))
tmp = sort(sort(table(sample.info$Breed), decreasing = T)[1:10])
barplot(tmp, horiz = T, las = 1, main = "Top 10 Breeds", cex.lab = 1.4, cex.axis = 1.25, cex.main = 1.5,
        col = "grey", xlab = "Count")
dev.off()


#########################################################
## Plot Fig S1
#########################################################
library(ggplot2)
library(dplyr)
library(gridExtra)

g1 = sample.info %>%
  ggplot(aes(x = Age..years., fill = sex2)) + 
  geom_histogram(bins = 16, alpha = 0.5, position = "identity") +
  # geom_density(alpha = 0.5)+
  scale_fill_manual(values = c("indianred", "skyblue"))+
  theme_bw() #+ theme(legend.position="none")

g2 = sample.info %>%
  ggplot(aes(x = final.weight..kg., fill = sex2)) + 
  geom_histogram(bins = 20, alpha = 0.5, position = "identity") +
  # geom_density(alpha = 0.5)+
  scale_fill_manual(values = c("indianred", "skyblue"))+
  theme_bw()

pdf("cohort_summary_by_sex.pdf", width = 10, height = 4)
par(mar = c(8,10,2,2))
grid.arrange(g1, g2, nrow = 1)
dev.off()

