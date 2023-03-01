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

#Making metdat file with genes and patient codes
df_w_genes <- read_excel("P21145_Results.xlsx", col_names=TRUE)
metdat <- df_w_genes[,1:3]
metdat$Gene <- factor(metdat$Gene)
metdat$Label <- as.numeric(metdat$Gene)

write.csv(metdat, file="metdat.csv")

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

df_all <- rbind(df_csf$`Lip-I`,df_csf$`Lip-II`)
rownames(df_all) <- df_all$X1
df_all <- df_all[,-1]
colnames(df_all) <- metdat$`COS Code`


stx <- unlist(metdat[metdat$Gene=="STXBP1",1])
grin <- unlist(metdat[metdat$Gene=="GRIN",1])
mec <- unlist(metdat[metdat$Gene=="MeCP2",1])
ctrl <- unlist(metdat[metdat$Gene=="Control",1])
pat <- unlist(metdat[metdat$Gene!="Control",1])


#Removing all features with NAs

df_all <- df_all[rowSums(is.na(df_all)) != ncol(df_all), ]
df_all_num <- apply(df_all, c(1,2), function(x) as.numeric(x))
df_all <- na.omit(df_all_num)
df_all <- data.frame(df_all)
df_all <- cbind(df_all,rownames(df_all))
colnames(df_all)[length(colnames(df_all))] <- "Species"
df_all$Species <- gsub("-sn1","",rownames(df_all))
df_all$Species <- gsub("-sn2","",df_all$Species)
df_all$Species <- gsub("-iso1","",df_all$Species)
df_all$Species <- gsub("-iso2","",df_all$Species)
df_all$Species <- gsub("-iso3","",df_all$Species)
df_all$Species <- gsub("GCDCA","ST 24:1;O4;G",df_all$Species)
df_all$Species <- gsub("CDCA","ST 24:1;O4",df_all$Species)
df_all$Species <- gsub("GCA","ST 24:1;O5;G",df_all$Species)
df_all$Species <- gsub("LCA-S","ST 24:1;O3",df_all$Species)
df_all$Species <- gsub("TCA","ST 24:1;O5;T",df_all$Species)
df_all$Species[df_all$Species=="9(s)-HODE"] <- "FA 18:2;O"
df_all$Species <- gsub("ChoE","CE",df_all$Species)

sapply(df_all$Species,isValidLipidName)

df_all <- aggregate(.~Species,data=df_all,sum)
colnames(df_all) <- gsub("F.","F ",colnames(df_all))
df_all <- data.frame(df_all,row.names = 1)


##Normalization

#Median normalization
#Log transformation
df_all_log <- log2(df_all)
colnames(df_all_log) <- gsub("F.","F ",colnames(df_all_log))

hist(t(df_all_log[,mec]),
     xlab="log 2", legend=NULL, main="Rett after log2 transformation", las=1)


df_stx <- df_all_log[colnames(df_all_log)%in%stx|colnames(df_all_log)%in%ctrl]
df_grin <- df_all_log[colnames(df_all_log)%in%grin|colnames(df_all_log)%in%ctrl]


## PCA based on normalized counts(normalized between samples) for patients of interest
pca_log_stx <- prcomp(t(df_stx), center=TRUE, scale=TRUE)
summary(pca_log_stx)
ggbiplot(pca_log_stx,var.axes=FALSE, groups=metdat$Gene[metdat$`COS Code`%in%colnames(df_stx)], ellipse=FALSE)

pca_log_grin <- prcomp(t(df_grin), center=TRUE, scale=TRUE)
summary(pca_log_grin)
ggbiplot(pca_log_grin,var.axes=FALSE, groups=metdat$Gene[metdat$`COS Code`%in%colnames(df_grin)], ellipse=FALSE)

ctrl.only <- df_all_log[,colnames(df_all_log)%in%ctrl]
stx.only <- df_all_log[,colnames(df_all_log)%in%stx]


s_w <- list()

for (i in 1:length(rownames(stx.only))) {
  s_w <- append(s_w,list(wilcox.test(unlist(stx.only[i,]),unlist(ctrl.only[i,]))))
}

names(s_w) <- rownames(stx.only)
sw_p <- lapply(s_w, function(x) x$p.value)
sw_p_adj <- p.adjust(unlist(sw_p),method="fdr")
s_int_w <- names(sw_p_adj[which(sw_p_adj<0.05)])

df_norm <- df_stx
df_norm_lab <- data.frame(cbind(as.factor(metdat$Gene[metdat$`COS Code` %in% c(stx,ctrl)]), t(df_norm)))
colnames(df_norm_lab)[1] <- "Label"
colnames(df_norm_lab)[-1] <- rownames(df_all_log)
log_df <- melt(df_norm_lab,id = "Label")
log_df$Label <- as.factor(log_df$Label)

maxim_val <- apply(df_norm_lab[,-1], 2, max)

stat.test <- log_df %>%
  group_by(variable) %>%
  wilcox_test(value ~ Label) %>%
  adjust_pvalue(method="fdr") %>%
  add_significance("p.adj") %>%
  mutate(y.position = maxim_val)
stat.test


p <- ggplot(log_df, aes(x = variable, y = value, color = as.factor(Label))) +  # ggplot function
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),legend.position = c(10, 0.7)) +
  stat_pvalue_manual(stat.test, x="variable",size = 3)

p 

######################GRIN####################################
ctrl.only <- df_all_log[,colnames(df_all_log)%in%ctrl]
grin.only <- df_all_log[,colnames(df_all_log)%in%grin]


g_w <- list()

for (i in 1:length(rownames(grin.only))) {
  g_w <- append(g_w,list(wilcox.test(unlist(grin.only[i,]),unlist(ctrl.only[i,]))))
}

names(g_w) <- rownames(grin.only)
gw_p <- lapply(g_w, function(x) x$p.value)
gw_p_adj <- p.adjust(unlist(gw_p),method="fdr")
g_int_w <- names(gw_p_adj[which(gw_p_adj<0.05)])

df_norm <- df_grin
df_norm_lab <- data.frame(cbind(as.factor(metdat$Gene[metdat$`COS Code` %in% c(stx,ctrl)]), t(df_norm)))
colnames(df_norm_lab)[1] <- "Label"
colnames(df_norm_lab)[-1] <- rownames(df_all_log)
log_df <- melt(df_norm_lab,id = "Label")
log_df$Label <- as.factor(log_df$Label)

maxim_val <- apply(df_norm_lab[,-1], 2, max)

stat.test <- log_df %>%
  group_by(variable) %>%
  wilcox_test(value ~ Label) %>%
  adjust_pvalue(method="fdr") %>%
  add_significance("p.adj") %>%
  mutate(y.position = maxim_val)
stat.test


p <- ggplot(log_df, aes(x = variable, y = value, color = as.factor(Label))) +  # ggplot function
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),legend.position = c(10, 0.7)) +
  stat_pvalue_manual(stat.test, x="variable",size = 3)

p 

