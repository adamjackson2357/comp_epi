# Format the data

dev.off()
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

suppressPackageStartupMessages({
  library(dplyr)
})

#### Get the data

meth <- readRDS("../data/Methylation.rds")
meth_tech <- readRDS("../data/Methylation_technical_confounders.rds")
covars <- readRDS("../data/Covariates_no_rep.rds")
proteins <- readRDS("../data/Proteins_selected_denoised_re.rds")

#### QA

# remove those who didn' pass the QC checks
covars <- covars %>% filter(QC.Warning == "Pass")
(table(covars$QC.Warning))

# checking for duplicates
length((row.names(covars))) == length(unique(row.names(covars)))

#### Data formatting

# subset and reorder data based on id
id_reorder <- function(df, ids){
  df <- df[row.names(df) %in% ids,]
  df <- df[match(ids, row.names(df)),]
  return(df)
}

# get ids that are in meth, meth_tech and covars
ids <- intersect(intersect(intersect(row.names(meth), row.names(meth_tech)), row.names(covars)), row.names(proteins))
length(ids)

# subset and re-order
meth <- id_reorder(meth, ids)
meth_tech <- id_reorder(meth_tech, ids)
covars <- id_reorder(covars, ids)
proteins <- id_reorder(proteins, ids)

# check all row names are equal
all(row.names(meth) == row.names(meth_tech)) &
  all(row.names(meth) == row.names(covars)) &
  all(row.names(meth) == row.names(proteins))

# Save
saveRDS(meth, "../data/meth_formatted.rds")
saveRDS(meth_tech, "../data/meth_tech_formatted.rds")
saveRDS(covars, "../data/covars_formatted.rds")
saveRDS(proteins, "../data/proteins_formatted.rds")
