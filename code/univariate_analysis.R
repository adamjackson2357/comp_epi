# Univariate analysis

dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

suppressPackageStartupMessages({
  library(parallel)
})

print("Getting data")
meth <- readRDS("../data/meth_denoised.rds")
proteins <- readRDS("../data/proteins_formatted.rds")

X <- meth
Y <- proteins

# get the beta and p-value for each linear regression
linear_reg = function(x, y) {
  model = lm(y~x)
  res = c(summary(model)$coefficients["x",1],
          summary(model)$coefficients[,"Pr(>|t|)"][2])
  names(res) = c("coef", "pval")
  return(res)
}

#### loop through proteins and parApply CpGs
nchunks = 24
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

dim(univ_res)

print("Post-processing")
univ_res <- data.frame(univ_res, protein = rep(colnames(Y), each = ncol(X)))
univ_res$bonf <- p.adjust(univ_res$pval, method = "bonf")
univ_res$bh <- p.adjust(univ_res$pval, method = "BH")

print("Save")
saveRDS(univ_res, "../data/univ_res.rds")

#### ANALYSIS

# Volcano plot
Volcano = function(results, title, annot = NULL, bonf) {
  par(mar = c(4.5, 4.5, 1, 1))
  plot(results$coef, -log10(results$pval), pch = 19,
       las = 1, cex = 0.5, xlab = expression(beta),
       ylab = expression(-log[10](p[value])),
       xlim=c(-0.6, 0.6), ylim = c(0, 12), main = title,
       col = ifelse(results$bh < 0.05,
                    yes = "tomato", no = "darkgrey"))
  if (!is.null(annot)) {

    text(results$coef, -log10(results$pval), pos = 3,
         offset = 0.2, cex = 0.5, labels = ifelse(results$bonf < 0.05, yes = annot, no = ""))
  }
  abline(v = 0, lty = 3)
  abline(h = -log10(bonf), lty = 2, col = "darkred")
  legend("topleft", col = c("darkred", "tomato", "darkgrey"),
         lty = c(2, NA, NA), pch = c(NA, 19, 19), cex = 0.7,
         legend = c("Bonferroni threshold at 0.05", "FDR significant hits", "Not significant"))
}

# get the bonferroni threshold
bonf <- 0.05/nrow(univ_res)

# volcano plots for all significant proteins
sig_proteins <- unique(univ_res[univ_res$bonf <= 0.05,][,"protein"])
png(file="../figures/sig_volcano.png", width=800, height=400)
old.par <- par(mfrow=c(1, 2))

for (i in 1:length(sig_proteins)){
  protein <- sig_proteins[i]
  results = univ_res[univ_res$protein == protein,]
  Volcano(results, title=protein, annot=rownames(results), bonf)
}

par(old.par)
dev.off()

# volcano plots for all proteins
proteins <- unique(univ_res$protein)
pdf(file="../figures/volcano_plots.pdf")
old.par <- par(mfrow=c(3, 3))

for (i in 1:length(proteins)){
  protein <- proteins[i]
  results = univ_res[univ_res$protein == protein,]
  Volcano(results, title=protein, annot=rownames(results), bonf)
}

par(old.par)
dev.off()


head(univ_res)



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
