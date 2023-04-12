load("data/metadata_flow_survey.rda")

## Figure 1
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

pdf("breed_summary_20230331.pdf", width = 6, height = 6)
par(mar = c(8,10,2,2))
tmp = sort(sort(table(sample.info$Breed), decreasing = T)[1:10])
barplot(tmp, horiz = T, las = 1, main = "Top 10 Breeds", cex.lab = 1.4, cex.axis = 1.25, cex.main = 1.5,
        col = "grey", xlab = "Count")
dev.off()

