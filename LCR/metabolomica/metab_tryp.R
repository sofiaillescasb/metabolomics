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

pdf("tryp_met_rett.pdf")

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
df_csf <- lapply(df_csf,function(x) x[-c(1:2),])

df_all <- df_csf$`Tryptophan & metabolites`
rownames(df_all) <- df_all$X1
df_all <- df_all[,-1]
colnames(df_all) <- metdat$COS.Code

mec <- unlist(metdat[metdat$Gene=="MeCP2",1])
ctrl <- unlist(metdat[metdat$Gene=="Control",1])

df_all <- df_all[colnames(df_all)%in%mec|colnames(df_all)%in%ctrl]

#Removing all features with NAs

df_all <- df_all[rowSums(is.na(df_all)) != ncol(df_all), ]
df_all_num <- apply(df_all, c(1,2), function(x) as.numeric(x))
df_all <- na.omit(df_all_num)
df_all <- data.frame(df_all)
df_all <- cbind(df_all,rownames(df_all))
colnames(df_all)[length(colnames(df_all))] <- "Metabolite"


df_all <- aggregate(.~Metabolite,data=df_all,sum)
colnames(df_all) <- gsub("F.","F ",colnames(df_all))
df_all <- data.frame(df_all,row.names = 1)


##Normalization

#Median normalization
#Log transformation
df_all_log <- log2(df_all)
colnames(df_all_log) <- gsub("F.","F ",colnames(df_all_log))

hist(t(df_all_log[,mec]),
     xlab="log 2", legend=NULL, main="Rett after log2 transformation", las=1)

hist(t(df_all_log[,ctrl]),
     xlab="log 2", legend=NULL, main="Control after log2 transformation", las=1)
##Univariate Analysis

#Wilcoxon rank sum test

ctrl.only <- df_all_log[,colnames(df_all_log)%in%ctrl]
mec.only <- df_all_log[,colnames(df_all_log)%in%mec]

m_w <- list()

for (i in 1:length(rownames(mec.only))) {
  m_w <- append(m_w,list(wilcox.test(unlist(mec.only[i,]),unlist(ctrl.only[i,]))))
}

names(m_w) <- rownames(mec.only)
mw_p <- lapply(m_w, function(x) x$p.value)
mw_p_adj <- p.adjust(unlist(mw_p),method="fdr")
m_int_w <- names(mw_p_adj[which(mw_p_adj<0.05)])

df_norm <- df_all_log
df_norm_lab <- data.frame(cbind(as.factor(metdat$Gene[metdat$COS.Code %in% c(mec,ctrl)]), t(df_norm)))
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
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  stat_pvalue_manual(stat.test, x="variable")

p 

#oPLS-DA to find most important features that separate patients from controls

library(ropls)


pd_m <- opls(df_norm_lab[,-1],df_norm_lab$Label,predI = 1,orthoI = 1)
plot(pd_m,typeVc ="x-score",)
vip_m <- data.frame(getVipVn(pd_m))
colnames(vip_m) <- "VIP"
vip_m$"Metabolite" <- rownames(vip_m)
vip_m_o <- vip_m[order(vip_m$VIP,decreasing = TRUE),]


vip_m_o %>%
  filter(!is.na(VIP)) %>%
  arrange(VIP) %>%
  tail(20) %>%
  mutate(Metabolite=factor(Metabolite, Metabolite)) %>%
  ggplot( aes(x=Metabolite, y=VIP) ) +
  geom_segment( aes(x= Metabolite,xend=Metabolite, y=0, yend=VIP), color="grey") +
  geom_point(size=3, color="#69b3a2") +
  coord_flip() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none"
  ) +
  xlab("") +
  ylab("VIP scores GCMS")



var_m <- vip_m$Metabolite[which(vip_m$VIP>1)]

interst <- rownames(vip_m_o)[rownames(vip_m_o)%in%m_int_w]

##Classification and clustering using selected variables
#Random Forest #Very bad, no good results

library(randomForest)

train <- sample(nrow(df_norm_lab), 0.7*nrow(df_norm_lab), replace = FALSE )
TrainSet <- df_norm_lab[train,]
ValidSet <- df_norm_lab[-train,]
modelmec2 <- randomForest(TrainSet[,-1],as.factor(TrainSet[,1]), data=TrainSet, ntree = 5000, mtry=8,importance = TRUE)
modelmec2
pred_val <- predict(modelmec2,ValidSet[,-1])
table(pred_val,ValidSet[,1])
varImpPlot(modelmec2)


#PCA
pca_m <- prcomp(df_norm_lab[,-1])
ggbiplot(pca_m,var.axes=FALSE, groups=metdat[metdat$Gene%in%c("MeCP2","Control"),]$Gene, ellipse=FALSE)
pca_m <- opls(df_norm_lab[,-1], predI = 2)
plot(pca_m,typeVc ="x-score",parAsColFcVn =df_norm_lab[,1])
plot(pca_m,typeVc ="x-loading",parAsColFcVn =df_norm_lab[,1])
l_df <- abs(getLoadingMN(pca_m))

#tsne

library(M3C)

perpo <- round(length(rownames(df_norm_lab))^(1/2)) 
tsne(t(df_norm_lab[,-1]),perplex = perpo,labels=as.factor(df_norm_lab$Label))
tsne(t(df_norm_lab[,-1]),perplex = perpo-1,labels=as.factor(df_norm_lab$Label))
tsne(t(df_norm_lab[,-1]),perplex = perpo-2,labels=as.factor(df_norm_lab$Label))
tsne(t(df_norm_lab[,-1]),perplex = perpo+1,labels=as.factor(df_norm_lab$Label))

#umap

umap(t(df_norm_lab[,-1]),labels=as.factor(df_norm_lab$Label))



dev.off()


