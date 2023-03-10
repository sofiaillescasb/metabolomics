#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Please provide the name of the analysis you wish to use (tryp, lip1, lip2, gcms)", call.=FALSE)
}

setwd("/Users/sofiaillescas/Desktop/LAB/LCR/metabolomica")

library(factoextra)
library(purrr)
library(readxl)
library(imputeTS)
library(reshape2)
library(ggpubr)
library(rstatix)


metdat <- read.csv("metadata.csv", header=TRUE, row.names = 1)

#Importing excel file with one sheet per analysis (e.g. lipids I, )
df_list <- map(set_names(excel_sheets("P21145_Results_final.xlsx")),
               read_excel, path = "P21145_Results_final.xlsx",
)


df_list <- lapply(df_list,t)
df_list <- lapply(df_list,data.frame)

#Excluding plasma results to focus only on cerebrospinal fluid samples
df_csf <- df_list[names(df_list) %in% "Plasma_GCMS" == FALSE]
#Removing plasma samples form remaining results columns
df_csf <- lapply(df_csf,function(x) x[,1:43])
df_csf <- lapply(df_csf,function(x) x[-c(1:3),])

arg_lst <- list()
arg_lst <- sapply(names(df_csf), append, arg_lst)
names(arg_lst) <- c("tryp", "lip1", "lip2", "gcms")
args1 <- which(names(arg_lst)==args[1])


df_all <- df_csf[[args1]]
rownames(df_all) <- df_all$X1
df_all <- df_all[,-1]
colnames(df_all) <- metdat$COS.Code

#Removing all features with over 25% missing values

df_all <- df_all[rowSums(is.na(df_all))<2,]
df_all_num <- apply(df_all, c(1,2), function(x) as.numeric(x))

mins <- apply(df_all_num,2,min,na.rm=TRUE)

df_all2 <- na_replace(df_all_num)
df_all2 <- df_all2+mins
df_all <- data.frame(df_all2)
colnames(df_all) <- gsub("F.","F ",colnames(df_all))

#Standardizing lipid species names
if (args1%in%c("Lip-I","Lip-II")) {
  df_all <- cbind(df_all,rownames(df_all))
  colnames(df_all)[length(colnames(df_all))] <- "Species"

  df_all$Species <- gsub("-sn1","",rownames(df_all))
  df_all$Species <- gsub("ChoE","CE",df_all$Species)
  df_all$Species <- gsub("-sn2","",df_all$Species)
  df_all$Species <- gsub("-iso1","",df_all$Species)
  df_all$Species <- gsub("-iso2","",df_all$Species)
  df_all$Species <- gsub("-iso3","",df_all$Species)
  df_all$Species <- gsub("GCDCA","ST 24:1;O4;G",df_all$Species)
  df_all$Species <- gsub("CDCA","ST 24:1;O4",df_all$Species)
  df_all$Species <- gsub("GCA","ST 24:1;O5;G",df_all$Species)
  df_all$Species <- gsub("LCA-S","ST 24:1;O3",df_all$Species)
  df_all$Species <- gsub("TCA","ST 24:1;O5;T",df_all$Species)
  df_all$Species <- gsub("Dehydroisoandrosterone 3-sulfate","DHEA-S",df_all$Species)
  df_all$Species[df_all$Species=="9(s)-HODE"] <- "FA 18:2;O"


  sapply(df_all$Species,isValidLipidName)

  df_all <- aggregate(.~Species,data=df_all,sum)
  colnames(df_all) <- gsub("F.","F ",colnames(df_all))
  df_all <- data.frame(df_all,row.names = 1)
}



##########################################Normalization###############################################

#Log transformation

df_all_log <- log2(df_all)
colnames(df_all_log) <- gsub("F.","F ",colnames(df_all_log))

genes <- unique(metdat$Gene)
sep_pat <- lapply(genes,function(x) unlist(metdat[metdat$Gene==x,1]))
names(sep_pat) <- genes

sep_vals <- lapply(sep_pat, function(x) df_all_log[,x])
names(sep_vals) <- genes

for (i in 1:length(sep_pat)) {
  hist(t(df_all_log[,sep_pat[[i]]]),
       xlab="log 2", legend=NULL, main= paste0(names(sep_pat[i])," after log2 transformation"), las=1)
}



#Adding labels to separate data by group

df_norm <- df_all_log
df_norm_lab <- data.frame(cbind(as.factor(metdat$Gene), t(df_norm)))
colnames(df_norm_lab)[1] <- "Label"
colnames(df_norm_lab)[-1] <- rownames(df_all_log)
log_df <- melt(df_norm_lab,id = "Label")
log_df$Label <- as.factor(log_df$Label)

maxim_val <- apply(df_norm_lab[,-1], 2, max)
log_df$Class <- sapply(log_df$variable, function(x) substring(x,1,3))
log_df$Class <- gsub("-","-O",log_df$Class)

stat.test <- function(g1,g2=10) {
  log_df[which(log_df$Label%in%c(g1,g2)),] %>%
  group_by(Class,variable) %>%
  wilcox_test(value ~ Label) %>%
  adjust_pvalue(method="fdr") %>%
  add_significance("p.adj") %>%
  mutate(y.position = 1)
}

uni_lab <- unique(log_df$Label[log_df$Label!=10])
sta_lst <- lapply(uni_lab, stat.test)
names(sta_lst) <- unique(metdat$Gene[metdat$Gene!="MeCP2"])
m_int_w <- lapply(sta_lst, function(x) x$variable[x$p.adj.signif!="ns"])
m_int_w <- unlist(m_int_w)

write.table(m_int_w, file="int_features.txt", append = TRUE, col.names = FALSE, row.names = TRUE, quote = FALSE)
