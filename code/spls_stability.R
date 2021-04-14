# SPLS Stability Analysis

dev.off
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

suppressPackageStartupMessages({
  library(sgPLS)
  library(parallel)
})

source("spls_calibration.R")

# data
X <- readRDS("../data/meth_denoised.rds")
Y <- readRDS("../data/proteins_formatted.rds")

# 10000 repeats takes 66 mins with this
nrepeat=10000
MaxVarX=50

print("Parameters:")
print(paste0("X dimensions: ", dim(X)))
print(paste0("Y dimensions: ", dim(Y)))
print(paste0("Repeats: ", nrepeat))
print(paste0("MaxVar: ", MaxVarX))
print(paste0("Cores available: ", detectCores()))

print("Starting stability analysis...")
print(paste0("Time: ", Sys.time()))
t0=Sys.time()
Summary <- StabilityPLS(X, Y, nrepeat, MaxVarX=MaxVarX)
t1=Sys.time()
print("Finished stability analysis.")
print(difftime(t1, t0, units = "mins"))
print(paste0("Time: ", Sys.time()))

print("Saving results")
saveRDS(Summary, "../pls_results/pls_stability.rds")
