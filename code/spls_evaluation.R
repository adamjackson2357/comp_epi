# SPLS Evaluation

dev.off
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

suppressPackageStartupMessages({
  library(sgPLS)
  library(pheatmap)
  library(ggplot2)
  library(formattable)
})
source("spls_calibration.R")

# data
meth <- readRDS("../data/meth_denoised.rds")
proteins <- readRDS("../data/proteins_formatted.rds")
meth_annot <- readRDS("../data/meth_annot_formatted.rds")
protein_names <- colnames(proteins)
linked <- table(data.frame(protein=meth_annot$alt.name, cpg=row.names(meth_annot)))
non_linked <- linked == 0

# get the summary dataframes
univ_calib = readRDS("../pls_results/univ_pls_calibration.rds")
multi_calib = readRDS("../pls_results/pls_calibration.rds")
multi_stab = readRDS("../pls_results/pls_stability.rds")
multi_stab=multi_stab[order(nrow(multi_stab):1),]

#### UNIVARIATE EVALUATION - best models

keepX_univ = NULL
y = unique(univ_calib$y)
for (i in 1:length(y)){
  y_summary <- univ_calib[univ_calib$y == i,]
  keepx <- as.numeric(y_summary[which.min(y_summary$MSEP),]$NVarX)
  keepX_univ <- rbind(keepX_univ, keepx)
}
keepX_univ = data.frame(proteins=colnames(proteins), keepX_univ)
rownames(keepX_univ) <- NULL

# number of features with lowest MSEP
ggplot(keepX_univ, aes_string(x="keepX_univ")) +
  geom_histogram(bins=100) +
  theme_classic() +
  ylab("Frequency")

#### Get the best model features
set.seed(1)
X_loadings = NULL
X_explained = NULL
for (i in protein_names){
  keepx = as.numeric(keepX_univ[keepX_univ$protein == i,][2])
  model <- spls(meth, proteins[,i], keepX=keepx, ncomp=1, mode='regression')
  x_loadings <- t(model$loadings$X[,1])
  x_explained <- round(model$explained_variance$X, 2)
  X_loadings <- rbind(X_loadings, x_loadings)
  X_explained <- rbind(X_explained, x_explained)
}

# only cpgs which have at least one non-zero loading
non_zero_cpgs_best <- apply(X_loadings != 0, MARGIN=2, sum)
col_names <- paste0(keepX_univ$proteins, " (", keepX_univ$keepX_univ,", ", X_explained, ")")
pheatmap(X_loadings, cluster_rows = FALSE, cluster_cols = FALSE,
         labels_col="",labels_row=col_names, width=10, height=6.6,
         filename = "../figures/univ_best_heatmap.png")

#### UNIVARIATE EVALUATION - median features

# set the number of features as the median number of cpgs per protein
num_features <- median(table(meth_annot$alt.name))
set.seed(1)
X_loadings = NULL
X_explained = NULL
for (i in protein_names){
  keepx = as.numeric(keepX_univ[keepX_univ$protein == i,][2])
  model <- spls(meth, proteins[,i], keepX=num_features, ncomp=1, mode='regression')
  x_loadings <- t(model$loadings$X[,1])
  x_explained <- round(model$explained_variance$X, 2)
  X_loadings <- rbind(X_loadings, x_loadings)
  X_explained <- rbind(X_explained, x_explained)
}


# only cpgs which have at least one non-zero loading
non_zero_cpgs_median <- apply(X_loadings != 0, MARGIN=2, sum)
col_names <- paste0(keepX_univ$proteins, " (", num_features,", ", X_explained, ")")
pheatmap(X_loadings, cluster_rows = FALSE, cluster_cols = FALSE,
         labels_col="", labels_row=col_names, width=10, height=6.6,
         filename = "../figures/univ_median_heatmap.png")

#### UNIVARIATE EVALUATION - fixed features as number of related proteins

# set the number of features as the number of cpgs per protein
num_features <- table(meth_annot$alt.name)
set.seed(1)
X_loadings = NULL
X_explained = NULL
for (i in protein_names){
  keepx = as.numeric(keepX_univ[keepX_univ$protein == i,][2])
  model <- spls(meth, proteins[,i], keepX=num_features[i], ncomp=1, mode='regression')
  x_loadings <- t(model$loadings$X[,1])
  x_explained <- round(model$explained_variance$X, 2)
  X_loadings <- rbind(X_loadings, x_loadings)
  X_explained <- rbind(X_explained, x_explained)
}

# only cpgs which have at least one non-zero loading
non_zero_cpgs_proteins <- apply(X_loadings != 0, MARGIN=2, sum)
col_names <- paste0(keepX_univ$proteins, " (", num_features,", ", X_explained, ")")
pheatmap(X_loadings, cluster_rows = FALSE, cluster_cols = FALSE,
         labels_row=col_names, width=10, height=6.6,
         labels_col = "", filename = "../figures/univ_protein_heatmap.png")
dim(X_loadings)

#### MULTIVARIATE EVALUATION

# plot
png(filename="../figures/multi_pls_calibration.png", width=750, height=500)
PlotCalibsPLS(multi_calib)
dev.off()

# get the number of features with the best MSEP
keepX <- as.numeric(multi_calib[which.min(multi_calib$MSEP),]$NVarX)

# Get the best model features
set.seed(1)
model <- spls(meth, proteins, keepX= keepX, ncomp=1, mode='regression')

# non-zero loadings for X and Y
X_loadings <- model$loadings$X[,1][model$loadings$X[,1] !=0]
Y_loadings <- model$loadings$Y

# explained variance in X and Y
model$explained_variance$X
model$explained_variance$Y

# visualise the cpg loadings
cpg_loadings <- data.frame(cpg=factor(names(X_loadings), levels=names(X_loadings)),
                           loading=X_loadings)
ggplot(cpg_loadings, aes(x=cpg, y=loading)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylim(c(-0.5, 0.5))

# visualise the protein loadings
protein_loadings <- data.frame(protein=factor(rownames(Y_loadings), levels=rownames(Y_loadings)),
                               loading=Y_loadings[,1])
ggplot(protein_loadings, aes(x=protein, y=loading)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylim(c(-0.5, 0.5))

#### STABILITY ANALYSIS MULTIVARIATE

selected_X = 37

# filter for CpGs that were selected at least once
selected_once = multi_stab[multi_stab$NVarX == selected_X,][-1]
selected_once = selected_once > 0

# select only stably selected column names for the given optimum number of variables
selected_stab = multi_stab[multi_stab$NVarX == selected_X,][-1][,selected_once]
selected_stab = data.frame(CpG=colnames(selected_stab), proportion=t(selected_stab))
head(selected_stab)
colnames(selected_stab)[2] <- "Proportion"
selected_stab = selected_stab[order(-selected_stab$Proportion),]
selected_stab$CpG <- factor(selected_stab$CpG, levels=selected_stab$CpG)
tail(selected_stab)
str(selected_stab)
dim(selected_stab)

ggplot(selected_stab, aes_string(x="CpG", y="Proportion")) +
  geom_bar(stat="identity", colour="white") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  geom_hline(yintercept=0.8, linetype="dashed", color = "red", size=0.5) +
  geom_hline(yintercept=0.2, linetype="dashed", color = "red", size=0.5)

# those that are selected at least 20% of the time
# put breaks in at 80% and at the first 37

above_80 <- nrow(selected_stab[selected_stab$Proportion >= 0.8,])
gaps_col <- c(above_80)
gaps_row <- c(nrow(multi_stab)-selected_X, nrow(multi_stab)-selected_X+1)

selected <- selected_stab[selected_stab$Proportion >= 0.2,]$CpG
stable_summary <- multi_stab[,as.character(selected)]
pheatmap(stable_summary, cluster_rows = FALSE, cluster_cols = FALSE,
         gaps_col = gaps_col, gaps_row = gaps_row,
         angle_col = 45)

#### OVERALL Comparison

summary <- data.frame(CpG=names(non_zero_cpgs_best),
                      Best=non_zero_cpgs_best,
                      Median=non_zero_cpgs_median,
                      Protein=non_zero_cpgs_proteins)

summary <- merge(summary, selected_stab,by="CpG")

summary <- summary[order(-summary$Proportion),]
write.csv(summary, "../figures/summary.csv")
# png("../figures/univariate_comparisons.png", width = 350, height = 350)
# formattable(summary, 
#             align =c("l","c","c","c"), 
#             list(`Indicator Name` = formatter(
#               "span", style = ~ style(color = "grey",font.weight = "bold"))))
# dev.off()

head(selected_stab)
