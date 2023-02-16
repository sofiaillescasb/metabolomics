setwd("/Users/sofiaillescas/Desktop/LAB/LCR/metabolomica")

library(factoextra)
library(purrr)
library(readxl)
library(imputeTS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbiplot)
library(geneplotter)
library(limma)
library(factoextra)

#pdf("clustering_MeCP2_lip.pdf")

metadata <- read.csv("metadata.csv", header=TRUE, row.names = 1)

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
colnames(df_all) <- metadata$COS.Code
rownames(df_all) <- gsub("ChoE","CE",rownames(df_all))

mec <- unlist(metadata[metadata$Gene=="MeCP2",1])
ctrl <- unlist(metadata[metadata$Gene=="Control",1])

df_all <- df_all[colnames(df_all)%in%mec|colnames(df_all)%in%ctrl]

#Removing all features with NAs

df_all <- df_all[rowSums(is.na(df_all)) != ncol(df_all), ]
df_all_num <- apply(df_all, c(1,2), function(x) as.numeric(x))
df_all <- na.omit(df_all_num)
df_all <- data.frame(df_all)

##Normalization

#Median normalization
df_all_log <- data.frame(normalizeMedianValues(df_all))
colnames(df_all_log) <- gsub("F.","F ",colnames(df_all_log))

hist(t(df_all_log[,mec]),
             xlab="log 2", legend=NULL, main="Rett after median normalization", las=1)
hist(t(df_all_log[,ctrl]),
     xlab="log 2", legend=NULL, main="Control after scaling and centering", las=1)
#Log transformation
df_all_log <- log(df_all_log)

hist(t(df_all_log[,mec]),
     xlab="log 2", legend=NULL, main="Rett after log2 transformation", las=1)

hist(t(df_all_log[,ctrl]),
     xlab="log 2", legend=NULL, main="Control after log2 transformation", las=1)

#Scaling and centering
df_all_log <- data.frame(scale(df_all_log))
colnames(df_all_log) <- gsub("F.","F ",colnames(df_all_log))

hist(t(df_all_log[,mec]),
     xlab="log 2", legend=NULL, main="Rett after scaling and centering", las=1)

hist(t(df_all_log[,ctrl]),
     xlab="log 2", legend=NULL, main="Control after scaling and centering", las=1)

df_norm <- df_all_log
df_norm_lab <- data.frame(cbind(as.factor(metadata$Gene[metadata$COS.Code %in% c(mec,ctrl)]), t(df_norm)))
colnames(df_norm_lab)[1] <- "Label"

##Univariate Analysis

#t-test
ctrl.only <- df_all_log[,colnames(df_all_log)%in%ctrl]
mec.only <- df_all_log[,colnames(df_all_log)%in%mec]
m_t <- list()

for (i in 1:length(rownames(mec.only))) {
  m_t <- append(m_t,list(t.test(mec.only[i,],ctrl.only[i,])))
}

names(m_t) <- rownames(mec.only)
m_p <- lapply(m_t, function(x) x$p.value)
m_p_adj <- p.adjust(unlist(m_p),method="fdr")
m_int_t <- names(m_p_adj[which(m_p_adj<0.05)])



c_t <- list()

for (i in 1:length(rownames(ctrl.only))) {
  c_t <- append(c_t,list(t.test(ctrl.only[i,])))
}

names(c_t) <- rownames(ctrl.only)
c_p <- lapply(c_t, function(x) x$p.value)
c_p_adj <- p.adjust(unlist(c_p),method="fdr")
c_int_t <- names(c_p_adj[which(c_p_adj<0.05)])
cat(c_int_t, sep='\n')
cat(gsub("-...","",c_int_t, perl=TRUE), sep='\n')
#PLS-DA to find most important features that separate patients from controls
library(mixOmics)

X <- df_norm_lab[,-1]
Y <- df_norm_lab$Label
pd_m <- plsda(X,Y)
plotIndiv(pd_m, centroid=TRUE, ellipse = TRUE)
plotVar(pd_m,comp.select = 1)
plotLoadings(pd_m)
vip_m <- data.frame(vip(pd_m))
vip_ms <- vip_m[vip_m$comp1>=1,]
var_ms <- rownames(vip_ms)
##Classification and clustering using selected variables
#Random Forest

library(randomForest)
set.seed(2023)
trainm <- sample(nrow(df_norm_lab), 0.7*nrow(df_norm_lab), replace = FALSE)#number of columns unimportant
TrainSetm <- df_norm_lab[trainm,c("Label",var_ms)]
ValidSetm <- df_norm_lab[-trainm,c("Label",var_ms)]
modelmec <- randomForest(TrainSetm[,-1],as.factor(TrainSetm[,1]), data=TrainSetm, ntree = 1000, importance = TRUE, mtry = 10, proximity = TRUE)
modelmec
predict(modelmec,ValidSetm,proximity = TRUE)


#k-means clustering
#All features
fviz_nbclust(df_norm_lab, k.max = length(rownames(df_norm_lab))-1,kmeans, method = "wss")
fviz_nbclust(df_norm_lab, k.max = length(rownames(df_norm_lab))-1,kmeans, method = "silhouette")
fviz_nbclust(df_norm_lab, k.max = length(rownames(df_norm_lab))-1,kmeans, method = "gap_stat")

#k-means algorithm and plotting for all features
k <- kmeans(df_norm_lab, centers=2, nstart=30)
fviz_cluster(k, data=df_norm_lab)


#Only selected features
#Optimal number of clusters:
fviz_nbclust(df_norm_lab[,var_ms], k.max = length(rownames(df_norm_lab[,var_ms]))-1,kmeans, method = "wss")
fviz_nbclust(df_norm_lab[,var_ms], k.max = length(rownames(df_norm_lab[,var_ms]))-1,kmeans, method = "silhouette")
fviz_nbclust(df_norm_lab[,var_ms], k.max = length(rownames(df_norm_lab[,var_ms]))-1,kmeans, method = "gap_stat")

#k-means algorithm and plotting for selected variables
k2 <- kmeans(df_norm_lab[,var_ms], centers=2, nstart=30)
fviz_cluster(k2, data=df_norm_lab[,var_ms])

#dev.off()

m_int_t
