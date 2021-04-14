# Denoise the methylation data
# 
dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

suppressPackageStartupMessages({
  library(pheatmap)
  library(lme4)
})

source("./omics/cpma.R")
source("./omics/mlm.R")
source("./omics/mlmer.R")
source("./omics/pqq.R")
source("./omics/private.R")
source("./omics/setutils.R")
source("./omics/utils.R")

#### Get the data

meth <- readRDS("../data/meth_formatted.rds")
meth_tech <- readRDS("../data/meth_tech_formatted.rds")

# check all row names are equal
print(paste0("Check rownames: ", all(row.names(meth) == row.names(meth_tech))))

# Pre-denoising PCA

print("Pre-denoising PCA starting...")
t0=Sys.time()

pcaX = prcomp(meth)
ev = with(pcaX, sdev**2/sum(sdev**2))

png("../figures/meth_pca.png", width = 350, height = 350)
plot(ev, pch = 19, col = "navy", xlab = "# of PCs",
     ylab = "Proportion of EV", ylim = c(0, 1.2), cex = 0.3)
points(cumsum(ev), pch = 19, col = "tomato", cex = 0.3)
legend("top", pch = 19, col = c("navy", "tomato"),
       legend = c("EV", "Cumulative EV"), cex = 0.9, horiz = T)
dev.off()

t1=Sys.time()
print("Pre-denoising PCA complete")
print(t1-t0)

# Denoise

print("Denoising starting...")
t0=Sys.time()

model = suppressWarnings(mlmer(meth ~ (1 | chip) + (1 | chip.pos),
                               data = meth_tech, save.residuals = TRUE,
                               save.ranks = FALSE))
meth_denoised <- model$residuals
saveRDS(meth_denoised, "../data/meth_denoised.rds")

t1=Sys.time()
print("Denoising Complete")
print(t1-t0)

# Post-denoising PCA

print("Post-denoising PCA starting...")
t0=Sys.time()

pcaX = prcomp(meth_denoised)
ev = with(pcaX, sdev**2/sum(sdev**2))

png("../figures/meth_denoised_pca.png", width = 350, height = 350)
plot(ev, pch = 19, col = "navy", xlab = "# of PCs",
     ylab = "Proportion of EV", ylim = c(0, 1.2), cex = 0.3)
points(cumsum(ev), pch = 19, col = "tomato", cex = 0.3)
legend("top", pch = 19, col = c("navy", "tomato"),
       legend = c("EV", "Cumulative EV"), cex = 0.9, horiz = T)
dev.off()

t1=Sys.time()
print("Post-denoising PCA complete")
print(t1-t0)
