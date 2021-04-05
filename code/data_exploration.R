# data exploration

dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
suppressPackageStartupMessages({
  library(yaml)
})

config <- read_yaml("../configs/main.yml")
proteins_fname <- config$proteins_fname
meth_fname <- config$meth_fname
meth_tech_fname <- config$meth_tech_fname
covars_fname <- config$covars_fname

## Get the data

# get the protein data
proteins <- readRDS(proteins_fname)
dim(proteins)

# get the Methylation data
# meth <- readRDS("../data/Methylation.rds")
# saveRDS(meth[,1:10], "../data/meth_sample.rds")
meth <- readRDS(meth_fname)
dim(meth)
head(meth)

# get the chips
meth_tech <- readRDS(meth_tech_fname)
dim(meth_tech)
head(meth_tech)
meth_tech %>%
  group_by(chip) %>%
  summarise(n = n())

# get the covariates
covars <- readRDS(covars_fname)
dim(covars)
head(covars)

## Get data for the same set of ids

ids <- intersect(intersect(intersect(row.names(meth), row.names(meth_tech)), row.names(proteins)), row.name(covars))

meth <- meth[row.names(meth) %in% ids,]
meth <- meth[match(ids, row.names(meth)),]

meth_tech <- meth_tech[row.names(meth_tech) %in% ids,]
meth_tech <- meth_tech[match(ids, row.names(meth_tech)),]

proteins <- proteins[row.names(proteins) %in% ids,]
proteins <- proteins[match(ids, row.names(proteins)),]

covars <- covars[row.names(covars) %in% ids,]
covars <- covars[match(ids, row.names(covars)),]

# qa
nrow(meth) == nrow(meth_tech)
nrow(meth) == nrow(proteins)
nrow(meth) == nrow(covars)
all(row.names(meth) == row.names(meth_tech))
all(row.names(meth) == row.names(proteins))
all(row.names(meth) == row.names(covars))

# add the chip to the covariates
covars$meth_chip <- meth_tech$chip


### Table One


head(covars)

colSums(is.na(covars))
dim(covars)

