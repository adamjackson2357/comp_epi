# Denoise the methylation data

dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

nchunks = 450

suppressPackageStartupMessages({
  library(lme4)
  library(parallel)
})

print("Getting data")

meth <- readRDS("../data/meth_formatted.rds")
meth <- meth[,1:1000]
meth_tech <- readRDS("../data/meth_tech_formatted.rds")
covars <- readRDS("../data/covars_formatted.rds")

print("Check all row names are equal")
print(paste0("Check rownames: ", all(row.names(meth) == row.names(meth_tech)) & all(row.names(meth) == row.names(covars))))

print("Reformatting")
covars <- cbind(covars, meth_tech)

# define the denoise function
denoise <- function(X){
  model = lmer(X ~ case + bmi + age.sample + smoking_status + (1 | chip) + (1 | chip.pos),
               data = covars, REML = FALSE)
  beta = fixef(model)[c("(Intercept)", "case1")]
  X = cbind(rep(1, length(covars$case)), as.numeric(covars$case))
  return(X %*% beta + resid(model))
}

print("Starting denoising...")

ids=as.character(cut(1:ncol(meth), breaks = nchunks, labels = 1:nchunks))

t0=Sys.time()
no_cores=min(detectCores(), nchunks) - 1
print(paste0("Number of cores: ", no_cores))

cl <- makeCluster(no_cores)
clusterExport(cl, c("meth", "covars", "ids", "denoise", "nchunks"))
clusterEvalQ(cl, library(lme4))

denoised=parSapply(cl=cl, 1:nchunks, FUN=function(k){
  X_chunk=meth[,ids==k]
  return(apply(X_chunk, 2, denoise))
})

stopCluster(cl)

print("Finished denoising")
t1=Sys.time()
print(t1-t0)

print("Unlisting denoised data")
meth_denoised = data.frame(row.names = row.names(meth))
for (k in 1:length(denoised)){
  meth_denoised <- cbind(meth_denoised, denoised[[k]])
}

print("Saving data")
saveRDS(meth_denoised, "../data/meth_denoised.rds")
