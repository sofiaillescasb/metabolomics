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
library(rgoslin)
library(reshape2)
library(ggpubr)
library(rstatix)

# pdf("clustering_MeCP2_lip.pdf")

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

df_all <- rbind(df_csf$`Lip-I`,df_csf$`Lip-II`)
rownames(df_all) <- df_all$X1
df_all <- df_all[,-1]
colnames(df_all) <- metdat$COS.Code


rownames(df_all) <- gsub("ChoE","CE",rownames(df_all))


mec <- unlist(metdat[metdat$Gene=="MeCP2",1])
ctrl <- unlist(metdat[metdat$Gene=="Control",1])

df_all <- df_all[colnames(df_all)%in%mec|colnames(df_all)%in%ctrl]

#Removing all features with over 25% missing values

na_m <- rowSums(is.na(df_all[,mec]))
na_c <- rowSums(is.na(df_all[,ctrl]))
df_all <- df_all[na_m < 0.25*(length(mec)-na_m)|na_c < 0.25*(length(ctrl)-na_c),]
df_all_num <- apply(df_all, c(1,2), function(x) as.numeric(x))

mins <- apply(df_all_num,2,min,na.rm=TRUE)

df_all2 <- na_replace(df_all_num)
df_all2 <- df_all2+mins
df_all <- data.frame(df_all2)
colnames(df_all) <- gsub("F.","F ",colnames(df_all))

#Standardizing lipid species names
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
df_all$Species <- gsub("Dehydroisoandrosterone 3-sulfate","DHEA-S",df_all$Species)
df_all$Species[df_all$Species=="9(s)-HODE"] <- "FA 18:2;O"


sapply(df_all$Species,isValidLipidName)

df_all <- aggregate(.~Species,data=df_all,sum)
colnames(df_all) <- gsub("F.","F ",colnames(df_all))
df_all <- data.frame(df_all,row.names = 1)


##########################################Normalization###############################################

#Log transformation
df_all_log <- log2(df_all)
colnames(df_all_log) <- gsub("F.","F ",colnames(df_all_log))

hist(t(df_all_log[,mec]),
     xlab="log 2", legend=NULL, main="Rett after log2 transformation", las=1)


hist(t(df_all_log[,ctrl]),
     xlab="log 2", legend=NULL, main="Control after log2 transformation", las=1)


########################################Univariate Analysis############################################

#Wilcoxon rank sum test

ctrl.only <- df_all_log[,colnames(df_all_log)%in%ctrl]
mec.only <- df_all_log[,colnames(df_all_log)%in%mec]

#Adding labels to separate data by group

df_norm <- df_all_log
df_norm_lab <- data.frame(cbind(as.factor(metdat$Gene[metdat$COS.Code %in% c(mec,ctrl)]), t(df_norm)))
colnames(df_norm_lab)[1] <- "Label"
colnames(df_norm_lab)[-1] <- rownames(df_all_log)
log_df <- melt(df_norm_lab,id = "Label")
log_df$Label <- as.factor(log_df$Label)

maxim_val <- apply(df_norm_lab[,-1], 2, max)
log_df$Class <- sapply(log_df$variable, function(x) substring(x,1,3))
log_df$Class <- gsub("-","-O",log_df$Class)

stat.test <- log_df %>%
  group_by(Class,variable) %>%
  wilcox_test(value ~ Label) %>%
  adjust_pvalue(method="fdr") %>%
  add_significance("p.adj") %>%
  mutate(y.position = maxim_val)

stat.test
m_int_w <- stat.test$variable[stat.test$p.adj.signif!="ns"]

log_df$Label <- sub(1,"Control",log_df$Label)
log_df$Label <- sub(2,"Rett",log_df$Label)
Group <- as.factor(log_df$Label)

#Plotting results as boxplots separated by lipid class
p <- ggplot(log_df, aes(x = variable, y = value, color = Group)) +  # ggplot function
  geom_boxplot() + 
  facet_wrap(~Class,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 7),legend.position = "none") +
  stat_pvalue_manual(stat.test, x="variable",size = 3)

p 

#########################################Multivariate Analysis###########################################

#PCA

pca_m <- prcomp(df_norm_lab[,-1], center = TRUE, scale = TRUE)
ggbiplot(pca_m,var.axes=FALSE, groups=metdat[metdat$Gene%in%c("MeCP2","Control"),]$Gene, ellipse=FALSE)
pca_m <- opls(df_norm_lab[,-1], predI = 2)
plot(pca_m,typeVc ="x-score",parAsColFcVn =df_norm_lab[,1])
plot(pca_m,typeVc ="x-loading",parAsColFcVn =df_norm_lab[,1])
l_df <- abs(getLoadingMN(pca_m))

#Random Forest 
library(vivid)
library(randomForest)
library(pROC)

set.seed(223)
train <- sample(nrow(df_norm_lab), 0.7*nrow(df_norm_lab), replace = FALSE)
TrainSet <- df_norm_lab[train,]
ValidSet <- df_norm_lab[-train,]
modelmec2 <- randomForest(TrainSet[,-1],as.factor(TrainSet[,1]), data=TrainSet, ntree = 5000, mtry=5, importance = )
modelmec2
pred_val <- predict(modelmec2,ValidSet[,-1])
table(pred_val,ValidSet[,1])
varImpPlot(modelmec2)
imprt  <- data.frame(round(importance(modelmec2), 2))
roc(ValidSet[,1],as.numeric(pred_val),plot=TRUE,auc.polygon=TRUE,grid=TRUE,print.auc=TRUE)

#Making heatmap
df_norm_lab$Label <- as.numeric(df_norm_lab$Label)
rf_fit <- vivi(fit = modelmec2, 
               data = df_norm_lab, 
               response = "Label",
               gridSize = 10,
               importanceType = "MeanDecreaseGini",
               nmax = 100,
               reorder = TRUE,
               class = 1,
               predictFun = NULL,
               )

#Selection of important lipids with orthogonal PLS-DA

library(ropls)

pd_m <- opls(df_norm_lab[,-1],df_norm_lab$Label,predI = 1,orthoI = 2)
plot(pd_m,typeVc ="x-score",)
vip_m <- data.frame(getVipVn(pd_m))
colnames(vip_m) <- "VIP"
vip_m$"Species" <- rownames(vip_m)
vip_m_o <- vip_m[order(vip_m$VIP,decreasing = TRUE),]

#lollipop plot with VIP scores
vip_m_o %>%
  filter(!is.na(VIP)) %>%
  arrange(VIP) %>%
  tail(20) %>%
  mutate(Species=factor(Species, Species)) %>%
  ggplot( aes(x=Species, y=VIP) ) +
  geom_segment( aes(x= Species,xend=Species, y=0, yend=VIP), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("VIP scores lipids")

var_m <- vip_m$Species[which(vip_m$VIP>1)]

interst <- rownames(vip_m_o)[rownames(vip_m_o)%in%m_int_w]


#tsne

library(M3C)

perpo <- round(length(rownames(df_norm_lab))^(1/2)) 
tsne(t(df_norm_lab[,-1]),perplex = perpo,labels=as.factor(df_norm_lab$Label))
tsne(t(df_norm_lab[,-1]),perplex = perpo-1,labels=as.factor(df_norm_lab$Label))
tsne(t(df_norm_lab[,-1]),perplex = perpo-2,labels=as.factor(df_norm_lab$Label))
tsne(t(df_norm_lab[,-1]),perplex = perpo+1,labels=as.factor(df_norm_lab$Label))

#umap

umap(t(df_norm_lab[,-1]),labels=as.factor(df_norm_lab$Label))


# dev.off()


