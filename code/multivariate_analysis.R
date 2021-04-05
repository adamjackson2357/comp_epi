# Multivariate analysis

dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# args=commandArgs(trailingOnly=TRUE)
# step=as.numeric(args[1])
# sample=as.numeric(args[2])

suppressPackageStartupMessages({
  library(sgPLS)
  library(pheatmap)
  library(parallel)
})

source("functions_spls.R")

# data
proteins <- readRDS("../data/proteins_formatted.rds")
meth <- readRDS("../data/meth_formatted.rds")

print(all(rownames(proteins) == rownames(meth)))

# select X and Y
X <- meth[,1:25]
Y <- proteins[,1:10]
dim(X)
dim(Y)

# params
set.seed(1)
nComp=2
nrepeat=1
sparsity='X'
stepY=1
stepX=10
MinVarX=1
MaxVarX=20
MinVarY=1
MaxVarY=NULL
narrow_stepX = 1
narrow_stepY = NULL

res <- CalibratesPLS(X, Y, nComp=nComp, nrepeat, sparsity, 
                     stepX, stepY, MinVarX, MaxVarX,
                     MinVarY, MaxVarY, narrow_stepX, narrow_stepY)
MysPLSX <- spls(X,Y,keepX=res$keepX, ncomp=nComp, mode='regression')

# all X loadings
MysPLSX$loadings$X

# selected X loadings
MysPLSX$loadings$X[,1][MysPLSX$loadings$X[,1] !=0]
MysPLSX$loadings$X[,2][MysPLSX$loadings$X[,2] !=0]

# all Y loadings
head(MysPLSX$loadings$Y)

# explained variance
MysPLSX$explained_variance$X
MysPLSX$explained_variance$Y
