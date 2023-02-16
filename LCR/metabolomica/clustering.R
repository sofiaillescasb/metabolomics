setwd("/Users/sofiaillescas/Desktop/LAB/LCR/metabolomica")

library(randomForest)
library(caret)
library(mixOmics)
library(factoextra)


metdat <- read.csv("metadata.csv", header=TRUE, row.names = 1)
df_norm <- read.csv("log_ad_all.csv", header=TRUE, row.names=1)
#colnames(df_norm) <- metdat$COS.Code
df_norm_lab <- data.frame(cbind(as.factor(metdat$Gene[metdat$Gene=="MeCP2"|metdat$Gene=="Control"]), t(df_norm)))
colnames(df_norm_lab)[1] <- "Label"

# lst <- list()
# mt <- unique(metdat$Gene[metdat$Gene!="Control"])
# for (i in mt) {
#   lst <- append(lst, list(df_norm_lab[metdat$Gene==i|metdat$Gene=="Control",]))
# }
# names(lst) <- mt

## sPLS-DA
#pdf("clustering_MeCP2_lip.pdf")
#Still need to separate into blocks metabolomics vs lipidomics

# X <- df_norm_lab[,-1]
# Y <- df_norm_lab$Label
# splsda_all <- splsda(X,Y)
# plotIndiv(splsda_all, centroid=TRUE, ellipse = TRUE)


# for (df in 1:length(lst)) {
#   X <- lst[[df]][,-1]
#   Y <- lst[[df]]$Label
#   spd <- splsda(X,Y)
#   plotIndiv(spd, centroid=TRUE, ellipse = TRUE)
#   plotVar(spd,comp.select = 1)
#   assign(paste0("features", df), selectVar(spd)$name)
#   plotLoadings(spd)
# }

# df_red <- df_norm_lab[,]

X <- df_norm_lab[,-1]
Y <- df_norm_lab$Label
spd_m <- splsda(X,Y)
plotIndiv(spd_m, centroid=TRUE, ellipse = TRUE)
plotVar(spd_m,comp.select = 1)
var_m <- selectVar(spd_m)$name[1:70]
plotLoadings(spd_m)
var_m[order(var_m)]

##Random Forest

set.seed(224)
# train <- sample(nrow(df_norm_lab), 0.7*nrow(df_norm_lab), replace = FALSE)
# TrainSet <- df_norm_lab[train,]
# ValidSet <- df_norm_lab[-train,]
# model <- randomForest(TrainSet[,-1],as.factor(TrainSet[,1]), data=TrainSet, ntree = 1000, importance = TRUE)
# model
# 

# for (j in 1:length(lst)) {
#   train <- sample(nrow(lst[[1]]), 0.7*nrow(lst[[1]]), replace = FALSE)
#   TrainSet <- df_norm_lab[train,]
#   ValidSet <- df_norm_lab[-train,]
#   m <- randomForest(TrainSet[,-1],as.factor(TrainSet[,1]), data=TrainSet, ntree = 1000, importance = TRUE)
#   print(m)
# }

# #STXBP1 vs control
# set.seed(224)
# train <- sample(nrow(df_stx), 0.7*nrow(df_stx), replace = FALSE)
# TrainSetst <- df_stx[train,]
# ValidSetst <- df_stx[-train,]
# modelstx <- randomForest(TrainSetst[,-1],as.factor(TrainSetst[,1]), data=TrainSetst, ntree = 1000, importance=TRUE)
# modelstx
# 
#MeCP2 vs control
set.seed(224)
trainm <- sample(nrow(df_norm_lab), 0.7*nrow(df_norm_lab), replace = FALSE)#number of columns unimportant
TrainSetm <- df_norm_lab[trainm,c("Label",var_m)]
ValidSetm <- df_norm_lab[-trainm,c("Label",var_m)]
modelmec <- randomForest(TrainSetm[,-1],as.factor(TrainSetm[,1]), data=TrainSetm, ntree = 1000, importance = TRUE)
modelmec
# 
# #GRIN vs control
# set.seed(224)
# traing <- sample(nrow(df_grin), 0.7*nrow(df_grin), replace = FALSE)
# TrainSetg <- df_grin[traing,]
# ValidSetg <- df_grin[-traing,]
# modelgrin <- randomForest(TrainSetg[,-1],as.factor(TrainSetg[,1]), data=TrainSetg, ntree = 10000, importance = TRUE)
# modelgrin



##K-means

# fviz_nbclust(lst$GRIN[,-1], k.max = length(rownames(lst$GRIN[,-1]))-1,kmeans, method = "wss")
# fviz_nbclust(lst$GRIN[,-1], k.max = length(rownames(lst$GRIN[,-1]))-1,kmeans, method = "silhouette")
# fviz_nbclust(lst$GRIN[,-1], k.max = length(rownames(lst$GRIN[,-1]))-1,kmeans, method = "gap_stat")
# k <- kmeans(lst$GRIN[,-1], centers=3, nstart=30)
# fviz_cluster(k, data=lst$GRIN[,-1])
# 
# fviz_nbclust(lst$STXBP1[,-1], k.max = length(rownames(lst$STXBP1[,-1]))-1,kmeans, method = "wss")
# fviz_nbclust(lst$STXBP1[,-1], k.max = length(rownames(lst$STXBP1[,-1]))-1,kmeans, method = "silhouette")
# fviz_nbclust(lst$STXBP1[,-1], k.max = length(rownames(lst$STXBP1[,-1]))-1,kmeans, method = "gap_stat")
# k <- kmeans(lst$STXBP1[,-1], centers=2, nstart=30)
# fviz_cluster(k, data=lst$STXBP1[,-1])
# 
# fviz_nbclust(lst$KearnSyre[,-1], k.max = length(rownames(lst$KearnSyre[,-1]))-1,kmeans, method = "wss")
# fviz_nbclust(lst$KearnSyre[,-1], k.max = length(rownames(lst$KearnSyre[,-1]))-1,kmeans, method = "silhouette")
# fviz_nbclust(lst$KearnSyre[,-1], k.max = length(rownames(lst$KearnSyre[,-1]))-1,kmeans, method = "gap_stat")
# k <- kmeans(lst$KearnSyre[,-1], centers=2, nstart=30)
# fviz_cluster(k, data=lst$KearnSyre[,-1])
# 
# fviz_nbclust(lst$FOLR1[,-1], k.max = length(rownames(lst$FOLR1[,-1]))-1,kmeans, method = "wss")
# fviz_nbclust(lst$FOLR1[,-1], k.max = length(rownames(lst$FOLR1[,-1]))-1,kmeans, method = "silhouette")
# fviz_nbclust(lst$FOLR1[,-1], k.max = length(rownames(lst$FOLR1[,-1]))-1,kmeans, method = "gap_stat")
# k <- kmeans(lst$FOLR1[,-1], centers=4, nstart=30)
# fviz_cluster(k, data=lst$FOLR1[,-1])
# 
# fviz_nbclust(lst$CDKL5[,-1], k.max = length(rownames(lst$CDKL5[,-1]))-1,kmeans, method = "wss")
# fviz_nbclust(lst$CDKL5[,-1], k.max = length(rownames(lst$CDKL5[,-1]))-1,kmeans, method = "silhouette")
# fviz_nbclust(lst$CDKL5[,-1], k.max = length(rownames(lst$CDKL5[,-1]))-1,kmeans, method = "gap_stat")
# k <- kmeans(lst$CDKL5[,-1], centers=3, nstart=30)
# fviz_cluster(k, data=lst$CDKL5[,-1])

# fviz_nbclust(lst$FOXG1[,-1], k.max = length(rownames(lst$FOXG1[,-1]))-1,kmeans, method = "wss")
# fviz_nbclust(lst$FOXG1[,-1], k.max = length(rownames(lst$FOXG1[,-1]))-1,kmeans, method = "silhouette")
# fviz_nbclust(lst$FOXG1[,-1], k.max = length(rownames(lst$FOXG1[,-1]))-1,kmeans, method = "gap_stat")
# k <- kmeans(lst$FOXG1[,-1], centers=3, nstart=30)
# fviz_cluster(k, data=lst$FOXG1[,-1])

fviz_nbclust(df_norm_lab[,var_m], k.max = length(rownames(df_norm_lab[,var_m]))-1,kmeans, method = "wss")
fviz_nbclust(df_norm_lab[,var_m], k.max = length(rownames(df_norm_lab[,var_m]))-1,kmeans, method = "silhouette")
fviz_nbclust(df_norm_lab[,var_m], k.max = length(rownames(df_norm_lab[,var_m]))-1,kmeans, method = "gap_stat")
k <- kmeans(df_norm_lab[,var_m], centers=2, nstart=30)
fviz_cluster(k, data=df_norm_lab[,var_m])

# rm_ctrl <- c("CSF.15", "CSF.16","CSF.31")
# keep <- rownames(df_norm_lab)[!rownames(df_norm_lab)%in%rm_ctrl]

# fviz_nbclust(lst$BCKDK[,-1], k.max = length(rownames(lst$BCKDK[,-1]))-1,kmeans, method = "wss")
# fviz_nbclust(lst$BCKDK[,-1], k.max = length(rownames(lst$BCKDK[,-1]))-1,kmeans, method = "silhouette")
# fviz_nbclust(lst$BCKDK[,-1], k.max = length(rownames(lst$BCKDK[,-1]))-1,kmeans, method = "gap_stat")
# k <- kmeans(lst$BCKDK[,-1], centers=3, nstart=30)
# fviz_cluster(k, data=lst$BCKDK[,-1])
# 
# fviz_nbclust(lst$EEM[,-1], k.max = length(rownames(lst$EEM[,-1]))-1,kmeans, method = "wss")
# fviz_nbclust(lst$EEM[,-1], k.max = length(rownames(lst$EEM[,-1]))-1,kmeans, method = "silhouette")
# fviz_nbclust(lst$EEM[,-1], k.max = length(rownames(lst$EEM[,-1]))-1,kmeans, method = "gap_stat")
# k <- kmeans(lst$EEM[,-1], centers=2, nstart=30)
# fviz_cluster(k, data=lst$EEM[,-1])
# 
# fviz_nbclust(lst$NKH[,-1], k.max = length(rownames(lst$NKH[,-1]))-1,kmeans, method = "wss")
# fviz_nbclust(lst$NKH[,-1], k.max = length(rownames(lst$NKH[,-1]))-1,kmeans, method = "silhouette")
# fviz_nbclust(lst$NKH[,-1], k.max = length(rownames(lst$NKH[,-1]))-1,kmeans, method = "gap_stat")
# k <- kmeans(lst$NKH[,-1], centers=3, nstart=30)
# fviz_cluster(k, data=lst$NKH[,-1])

# fviz_nbclust(lst$EnoilCoA[,-1], k.max = length(rownames(lst$EnoilCoA[,-1]))-1,kmeans, method = "wss")
# fviz_nbclust(lst$EnoilCoA[,-1], k.max = length(rownames(lst$EnoilCoA[,-1]))-1,kmeans, method = "silhouette")
# fviz_nbclust(lst$EnoilCoA[,-1], k.max = length(rownames(lst$EnoilCoA[,-1]))-1,kmeans, method = "gap_stat")
# k <- kmeans(lst$EnoilCoA[,-1], centers=3, nstart=30)
# fviz_cluster(k, data=lst$EnoilCoA[,-1])



#dev.off()
