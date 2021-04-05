# Univariate analysis

# dev.off()
# rm(list=ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

suppressPackageStartupMessages({
  library(parallel)
})

nchunks = 450

print("Getting data")
meth <- readRDS("../data/meth_formatted.rds")
proteins <- readRDS("../data/proteins_formatted.rds")

X <- meth
Y <- proteins[,1:2]

# get the beta and p-value for each linear regression
linear_reg = function(x, y) {
  model = lm(y~x)
  res = c(summary(model)$coefficients["x",1],
          summary(model)$coefficients[,"Pr(>|t|)"][2])
  names(res) = c("coef", "pval")
  return(res)
}

#### loop through proteins and parApply CpGs

print("Starting univariate analysis...")
no_cores=min(detectCores(), nchunks)
print(paste0("Number of cores: ", no_cores))

t0 = Sys.time()

ids=as.character(cut(1:ncol(X), breaks = nchunks, labels = 1:nchunks))
cl <- makeCluster(no_cores)

univ_res = NULL
for (i in 1:ncol(Y)){
  print(i)
  y <- Y[,i]
  print(length(y))

  clusterExport(cl, c("X", "y", "ids", "linear_reg", "nchunks"))

  univL=parLapply(cl=cl, 1:nchunks, fun=function(k){
    X_chunk=X[,ids==k]
    return(t(apply(X_chunk, 2, FUN = linear_reg, y)))
  })

  univ = data.frame(row.names = colnames(X))
  for (j in 1:length(univL)){univ <- rbind(univ, univL[[j]])}
  univ_res <- rbind(univ_res, univ)
}

stopCluster(cl)

print("Finished univariate analysis.")
t1 = Sys.time()
print(t1 - t0)

print("Post-processing")
univ_res <- data.frame(univ_res, protein = rep(colnames(Y), each = ncol(X)))

bonferroni <- 0.05 / nrow(univ_res)
univ_res <- univ_res[univ_res$pval <= bonferroni,]

print("Save")
saveRDS(univ_res, "../data/univ_res.rds")


#### using loop and apply - old version


# get the beta and p-value for each linear regression
# anova version - don't think it's necessary as only one variable
# linear_reg = function(x, y) {
#   null_model = lm(y~1)
#   model = lm(y~x)
#   res = c(summary(model)$coefficients["x",1],
#           anova(null_model, model)[,"Pr(>F)"][2])
#   names(res) = c("coef", "pval")
#   return(res)
# }

# t0 = Sys.time()
# all_results = NULL
# for (i in 1:ncol(Y)){
#   y <- Y[,i]
#   univ = t(apply(X, 2, FUN = linear_reg, y))
#   all_results <- rbind(all_results, univ)
# }
# 
# # add the protein names
# all_results <- data.frame(all_results, protein = rep(colnames(Y), each = ncol(X)))
# 
# t1 = Sys.time()
# print(t1 - t0)
