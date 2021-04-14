# Multivariate SPLS Calibration

dev.off
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

suppressPackageStartupMessages({
  library(sgPLS)
  library(pheatmap)
  library(parallel)
})

source("spls_calibration.R")

# input
X <- readRDS("../data/meth_denoised.rds")
Y <- readRDS("../data/proteins_formatted.rds")
MaxVarX=100

# 100 repeats takes 12 mins
nrepeat=100

print("Parameters:")
print(paste0("X dimensions: ", dim(X)))
print(paste0("Y dimensions: ", dim(Y)))
print(paste0("Repeats: ", nrepeat))
print(paste0("Cores available: ", detectCores()))

print("Starting Calibration...")
print(paste0("Time: ", Sys.time()))
t0=Sys.time()
Summary <- MultivariatePLS(X, Y, nrepeat=nrepeat, stepX=1, stepY=1, MinVarX=1, MaxVarX=MaxVarX)

t1=Sys.time()
print("Finished calibration.")
print(difftime(t1, t0, units = "mins"))
print(paste0("Time: ", Sys.time()))

print("Saving results")
saveRDS(Summary, "../pls_results/multi_pls_calibration.rds")
