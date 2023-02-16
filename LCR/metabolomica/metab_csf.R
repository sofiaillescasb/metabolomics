setwd("/Users/sofiaillescas/Desktop/LAB/LCR/metabolomica")

library(purrr)
library(readxl)
library(imputeTS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbiplot)
library(geneplotter)
library(limma)
library(caret)
library(factoextra)

#Making metadata file with genes and patient codes
df_w_genes <- read_excel("P21145_Results.xlsx", col_names=TRUE)
metadata <- df_w_genes[,1:3]
metadata$Gene <- factor(metadata$Gene)
metadata$Label <- as.numeric(metadata$Gene)

write.csv(metadata, file="metadata.csv")

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
df_csf <- lapply(df_csf,function(x) x[-c(1:3),])

df_all <- rbind(df_csf$`Lip-I`,df_csf$`Lip-II`)
colnames(df_all) <- metadata$`COS Code`

mec <- unlist(metadata[metadata$Gene=="MeCP2",1])
ctrl <- unlist(metadata[metadata$Gene=="Control",1])

df_lip_mec <- df_all[,colnames(df_all)%in%mec|colnames(df_all)%in%ctrl]
write.table(read.csv("mecp2_lip.csv", sep=",", header=TRUE, row.names = 1), "mecp2_lip.tab", row.names=FALSE, col.names = FALSE)

df_all <- df_lip_mec

#,df_csf$CSF_GCMS
# df_all <- rbind(df_all,df_csf$`Lip-I`)
# df_all <- rbind(df_all,df_csf$`Lip-II`)



stx <- unlist(metadata[metadata$Gene=="STXBP1",1])
grin <- unlist(metadata[metadata$Gene=="GRIN",1])
mec <- unlist(metadata[metadata$Gene=="MeCP2",1])
ctrl <- unlist(metadata[metadata$Gene=="Control",1])
pat <- unlist(metadata[metadata$Gene!="Control",1])


write.csv(df_all, file="ad_all_raw.csv")
#Removing rows with no values and changing NAs to a small number

df_all <- df_all[rowSums(is.na(df_all)) != ncol(df_all), ]
df_all_num <- apply(df_all, c(1,2), function(x) as.numeric(x))
minim <- min(df_all_num,na.rm=TRUE)
# df_all<- na_replace(df_all_num,minim/5)
df_all <- na.omit(df_all_num)
df_all <- data.frame(df_all)

# colnames(df_all) <- metadata$`COS Code`

#Normalization
df_all_log <- data.frame(normalizeMedianAbsValues(df_all))
# colnames(df_all_log) <- metadata$`COS Code`

multidensity(as.list(df_all_log),
             xlab="log 2", legend=NULL, main="After median normalization", las=1)
# multidensity(as.list(df_all_log[,pat]),
#              xlab="log 2", legend=NULL, main="After median normalization", las=1)

df_all_log <- log(df_all_log, base=10)
# colnames(df_all_log) <- metadata$`COS Code`

multidensity(as.list(df_all_log),
             xlab="log 2", legend=NULL, main="After log transformation", las=1)
# multidensity(as.list(df_all_log[,pat]),
#              xlab="log 2", legend=NULL, main="After log transformation", las=1)

df_all_log <- data.frame(scale(df_all_log))
# colnames(df_all_log) <- metadata$`COS Code`

multidensity(as.list(df_all_log),
             xlab="log 2", legend=NULL, main="After scaling", las=1)
# multidensity(as.list(df_all_log[,pat]),
#              xlab="log 2", legend=NULL, main="After scaling", las=1)

colnames(df_all_log) <- gsub("F.","F ",colnames(df_all_log))

write.csv(df_all_log, file="log_ad_all.csv", row.names = TRUE)

##Dimensionality reduction and clustering

# PCA based on normalized counts(normalized between samples)
pca_log <- prcomp(t(df_all_log), center=TRUE, scale=TRUE)
summary(pca_log)
ggbiplot(pca_log,var.axes=FALSE, groups=metadata[metadata$Gene%in%c("MeCP2","Control"),]$Gene, ellipse=FALSE)


#Selecting patients of interest

# df_stx <- df_all_log[colnames(df_all_log)%in%stx|colnames(df_all_log)%in%ctrl]
# df_grin <- df_all_log[colnames(df_all_log)%in%grin|colnames(df_all_log)%in%ctrl]
# df_stx_grin <- df_all_log[colnames(df_all_log)%in%stx|colnames(df_all_log)%in%grin]
# 
# 
# ## PCA based on normalized counts(normalized between samples) for patients of interest
# pca_log_stx <- prcomp(t(df_stx), center=TRUE, scale=TRUE)
# summary(pca_log_stx)
# ggbiplot(pca_log_stx,var.axes=FALSE, groups=metadata$Gene[metadata$`COS Code`%in%colnames(df_stx)], ellipse=FALSE)
# 
# pca_log_grin <- prcomp(t(df_grin), center=TRUE, scale=TRUE)
# summary(pca_log_grin)
# ggbiplot(pca_log_grin,var.axes=FALSE, groups=metadata$Gene[metadata$`COS Code`%in%colnames(df_grin)], ellipse=FALSE)
# 
# pca_log_stxgrin <- prcomp(t(df_stx_grin), center=TRUE, scale=TRUE)
# summary(pca_log_stxgrin)
# ggbiplot(pca_log_stxgrin,var.axes=FALSE, groups=metadata$Gene[metadata$`COS Code`%in%colnames(df_stx_grin)], ellipse=FALSE)
# 
# pca_log_stxgrin_noo <- prcomp(t(df_stx_grin[,-c(4,5)]), center=TRUE, scale=TRUE)
# summary(pca_log_stxgrin_noo)
# ggbiplot(pca_log_stxgrin_noo,var.axes=FALSE, groups=metadata$Gene[metadata$`COS Code`%in%colnames(df_stx_grin[,-c(4,5)])], ellipse=FALSE)
# 
# #Finding relevant metabolites
# loads_stx <- summary(pca_log_stx)$importance
# pc_min_stx <- min(which(loads_stx[3,]>0.7))
# pca_log_stx$rotation <- pca_log_stx$rotation[,1:pc_min_stx]
# stx_mean <- apply(pca_log_stx$rotation, 2, mean)
# 
# loads_grin <- summary(pca_log_grin)$importance
# #how many components are needed to account for over 70% of the variability?
# pc_min_grin <- min(which(loads_grin[3,]>0.7))
# ord_rot_grin <- data.frame(pca_log_grin$rotation)
# in_ord_grin <- apply(ord_rot_grin[, 1:pc_min_grin], 2, function(x) order(abs(x), decreasing=TRUE))
# #getting names of metabolites with highest scores for the first 5 components
# main_stx <- data.frame(apply(in_ord_grin, 2, function(x) rownames(ord_rot_grin)[x][1:150]))
# 
# qplot(c(1:12), loads_grin[2,]) + 
#   geom_line() + 
#   xlab("Principal Component") + 
#   ylab("Variance Explained") +
#   ggtitle("Scree Plot") +
#   ylim(0, 1)

##Univariate analysis
#student t-test of metabolites STXBP1 vs control
# stx.only <- df_stx[,colnames(df_stx)%in%stx]
# ctrl.only <- df_stx[,colnames(df_stx)%in%ctrl]
# 
# df_sc <- cbind(rbind(rep("STXBP1",length(stx.only)),stx.only),rbind(rep("Control",length(ctrl.only)),ctrl.only))
# rownames(df_sc)[1] <- "Label"
# stx_var <- apply(t(df_sc)[,-1], 2, bartlett.test, t(df_sc)[,1])
# stx_var_p <- lapply(stx_var, function(x) x$p.value)
# stx_var_p_adj <- p.adjust(unlist(stx_var_p),method="fdr")
# min(stx_var_p_adj)
# stx_noth <- names(stx_var_p_adj[which(stx_var_p_adj<0.05)])
# 
# stx_t <- list()
# 
# for (i in 1:length(rownames(stx.only))) {
#   stx_t <- append(stx_t,list(t.test(stx.only[i,],ctrl.only[i,])))
# }
# 
# names(stx_t) <- rownames(stx.only)
# stx_p <- lapply(stx_t, function(x) x$p.value)
# stx_p_adj <- p.adjust(unlist(stx_p),method="fdr")
# stx_int_t <- names(stx_p_adj[which(stx_p_adj<0.05)])

#student t-test of metabolites GRIN vs control

#Confirming homogeneity of variance
# grin.only <- df_grin[,colnames(df_grin)%in%grin]
# ctrl.only <- df_grin[,colnames(df_grin)%in%ctrl]
# 
# df_gc <- cbind(rbind(rep("GRIN",length(grin.only)),grin.only),rbind(rep("Control",length(ctrl.only)),ctrl.only))
# rownames(df_gc)[1] <- "Label"
# grin_var <- apply(t(df_gc)[,-1], 2, bartlett.test, t(df_gc)[,1])
# grin_var_p <- lapply(grin_var, function(x) x$p.value)
# grin_var_p_adj <- p.adjust(unlist(grin_var_p),method="fdr")
# min(grin_var_p_adj)
# 
# grin_t <- list()
# 
# for (i in 1:length(rownames(grin.only))) {
#   grin_t <- append(grin_t,list(t.test(grin.only[i,],ctrl.only[i,],var.equal = TRUE)))
# }
# 
# names(grin_t) <- rownames(grin.only)
# 
# grin_p <- lapply(grin_t, function(x) x$p.value)
# grin_p_adj <- p.adjust(unlist(grin_p),method="fdr")
# grin_int_t <- names(grin_p_adj[which(grin_p_adj<0.05)]) #7 significantly different metabolites

mec.only <- df_all_log[,colnames(df_all_log)%in%mec]
ctrl.only <- df_all_log[,colnames(df_all_log)%in%ctrl]

df_m <- cbind(rbind(rep("MeCP2",length(mec.only)),mec.only),rbind(rep("Control",length(ctrl.only)),ctrl.only))
rownames(df_m)[1] <- "Label"
m_var <- apply(t(df_m)[,-1], 2, bartlett.test, t(df_m)[,1])
m_var_p <- lapply(m_var, function(x) x$p.value)
m_var_p_adj <- p.adjust(unlist(m_var_p),method="fdr")
min(m_var_p_adj)
m_noth <- names(m_var_p_adj[which(m_var_p_adj<0.05)])

m_t <- list()

for (i in 1:length(rownames(mec.only))) {
  m_t <- append(m_t,list(t.test(mec.only[i,],ctrl.only[i,])))
}

names(m_t) <- rownames(mec.only)
m_p <- lapply(m_t, function(x) x$p.value)
m_p_adj <- p.adjust(unlist(m_p),method="fdr")
m_int_t <- names(m_p_adj[which(m_p_adj<0.05)])


