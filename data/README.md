The following data is available:

EPIC N=382 and NOWAC FEM N=266
Replicated measurements (N=2) for 56 EPIC
16 controls (same sample measured 16 times)
720 measurements (704 samples and 16 controls) over 8 plates (N=90 samples each)

# Covariates_no_rep.rds: 
individual characteristics (row names with IDs)

###############

# Gene_expression_annotation.rds: 
annotation of gene expression data

# Gene_expression.rds: 
gene expression data

###############

Metabolite data has been prepared in the following way:
log2-transformation, exclusion of the features missing in more than 50% of the samples, de-noising for technical confounders and background subtraction (see lecture 7.2 in Computational Epidemiology)

# Metabolites_negative_imputed_selected_summarised_pc_centroid.rds: 
summarised clusters of metabolomics features (negative mode)

# Metabolites_negative_imputed_selected.rds: 
metabolomics features (negative mode)

# Metabolites_positive_imputed_selected_summarised_pc_centroid.rds: 
summarised clusters of metabolomics features (positive mode)

# Metabolites_positive_imputed_selected.rds: 
metabolomics features (positive mode)

# Clusters_pc_centroid_metab_neg.rds
groups of metabolomics features identified by clustering (negative mode)

# Clusters_pc_centroid_metab_pos.rds: 
groups of metabolomics features identified by clustering (positive mode)

################

# Methylation_technical_confounders.rds: 
technical covariates for methylation data (chip and position on the chip), not corrected for in the methylation data

# Methylation.rds: 
methylation data
M-value transformation (logit2), exclusion of the features missing in more than 30% of the samples, imputation, scaling

# hm450.rds: 
annotation of methylation data

################

# Proteins.rds: 
inflammatory proteins
raw data (log2-scale) as received from OLINK

# Sample_sheet_proteins.rds: 
technical covariates for inflammatory proteins
