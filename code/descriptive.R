# Basic descriptive stats

dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

suppressPackageStartupMessages({
  library(pheatmap)
})

covars <- readRDS("../data/covars_formatted.rds")
meth <- readRDS("../data/meth_formatted.rds")
proteins <- readRDS("../data/proteins_formatted.rds")

# basic descriptive
str(covars)
table(covars$QC.Warning)
hist(covars$storage_time)
hist(covars$age.sample)
table(covars$gender)
table(covars$cohort)
table(covars$centre)
table(covars$smoking_status)
hist(covars$tsq_smoking)
hist(covars$smok_duration)
hist(covars$smok_intensity)
hist(covars$packyears)
hist(covars$CSI)
table(covars$subtype)
hist(covars$ttd)
hist(covars$bmi)
table(covars$education_plco)
hist(covars$Monocytes)
hist(covars$B)
hist(covars$CD4T)
hist(covars$CD8T)
hist(covars$NK)
hist(covars$Eosinophils)
hist(covars$Neutrophils)
table(covars$Gel.status)
table(covars$chip)
table(covars$chip.pos)

# PROTEINS

# proteins heatmap
pheatmap(cor(proteins), breaks = seq(-1, 1, length.out = 100),
         show_rownames = FALSE, show_colnames = FALSE,
         filename = "../figures/proteins_heatmap.png")

# run a pca on proteins
pcaX = prcomp(proteins)
ev = with(pcaX, sdev**2/sum(sdev**2))

png("../figures/proteins_pca.png", width = 350, height = 350)
plot(ev, pch = 19, col = "navy", xlab = "# of PCs",
     ylab = "Proportion of EV", ylim = c(0, 1.2), cex = 0.3)
points(cumsum(ev), pch = 19, col = "tomato", cex = 0.3)
legend("top", pch = 19, col = c("navy", "tomato"),
       legend = c("EV", "Cumulative EV"), cex = 0.9, horiz = T)
dev.off()
