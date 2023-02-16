setwd("/Users/sofiaillescas/Desktop/LAB/LCR/metabolomica")

library(purrr)
library(readxl)
library(imputeTS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbiplot)

#Importing excel file with one sheet per analysis (e.g. lipids I, )
df_list <- map(set_names(excel_sheets("P21145_Results_final.xlsx")),
               read_excel, path = "P21145_Results_final.xlsx", 
)


df_list <- lapply(df_list,t)
df_list <- lapply(df_list,data.frame)


#Excluding plasma results to focus only on cerebrospinal fluid samples
df_csf <- df_list[names(df_list) %in% "Plasma_GCMS" == FALSE]
#Removing plasma samples form remaining results columns
df_csf <- lapply(df_csf,function(x) x[,1:42])

#Naming columns as COS code
coscode <- unlist(c(df_csf$CSF_GCMS["COS Code",]))
df_csf <- lapply(df_csf, set_names, coscode)#change to using colnames?
#Removing extra rows
df_csf <- lapply(df_csf,slice,-c(1:3))

#Making metadata file with genes and patient codes
df_w_genes <- read_excel("P21145_Results.xlsx", col_names=TRUE)
df_w_genes <- na_replace(df_w_genes,0)
metadata <- df_w_genes[,1:3]
metadata$Gene <- factor(metadata$Gene)
metadata$Label <- as.numeric(metadata$Gene)

#Converting character values to numbers and replacing NA with 0s
df_csf$`Tryptophan & metabolites` <- df_csf$`Tryptophan & metabolites`[rowSums(is.na(df_csf$`Tryptophan & metabolites`)) != ncol(df_csf$`Tryptophan & metabolites`), ]
df_csf$`Lip-I` <- df_csf$`Lip-I`[rowSums(is.na(df_csf$`Lip-I`)) != ncol(df_csf$`Lip-I`), ]
df_csf$`Lip-II` <- df_csf$`Lip-II`[rowSums(is.na(df_csf$`Lip-II`)) != ncol(df_csf$`Lip-II`), ]
df_csf$`CSF_GCMS` <- df_csf$`CSF_GCMS`[rowSums(is.na(df_csf$`CSF_GCMS`)) != ncol(df_csf$`CSF_GCMS`), ]

df_csf <- lapply(df_csf, apply, c(1,2), as.numeric)
df_csf <- lapply(df_csf, function(x) na_replace(x,0))
df_csf <- lapply(df_csf,data.frame)

#Normalization by change into log scale
df_csf_log <- lapply(df_csf,function(x) x+1)
df_csf_log <- lapply(df_csf_log,log)

df_csf_log$`Tryptophan & metabolites`[1,] <- metadata$Label
df_csf_log$`Lip-I`[1,] <- metadata$Label
df_csf_log$`Lip-II`[1,] <- metadata$Label
df_csf_log$CSF_GCMS[1,] <- metadata$Label

write.csv(df_csf_log$`Tryptophan & metabolites`, file="tryptophan_and_metabolites.csv")
write.csv(df_csf_log$`Lip-I`, file="lip1s.csv")
write.csv(df_csf_log$`Lip-II`, file="lip2.csv")
write.csv(df_csf_log$CSF_GCMS, file="csfgcms.csv")

# PCA based on normalized counts(normalized between samples)
pca_log1 <- prcomp(t(df_csf_log$`Tryptophan & metabolites`), center=TRUE, scale=TRUE)
ggbiplot(pca_log1,var.axes=FALSE, groups=metadata$Gene)

pca_log2 <- prcomp(t(df_csf_log$`Lip-I`), center=TRUE, scale=TRUE)
ggbiplot(pca_log2,var.axes=FALSE, groups=metadata$Gene)

pca_log3 <- prcomp(t(df_csf_log$`Lip-II`), center=TRUE, scale=TRUE)
ggbiplot(pca_log3,var.axes=FALSE, groups=metadata$Gene)

pca_log4 <- prcomp(t(df_csf_log$CSF_GCMS), center=TRUE, scale=TRUE)
ggbiplot(pca_log4,var.axes=FALSE, groups=metadata$Gene)

summary(pca_log)
