library(purrr)
library(tidyverse)
library(reshape2)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(knitr)
library(readxl)

metdat <- read_excel("//hsjdbcn.es/dfsroot/Recursos/metabolismosinaptico/SOFIA/metabolomics_rett_grin/Supplementary_1.xlsx",col_names = TRUE)

df_list <- map(purrr::set_names(excel_sheets("//hsjdbcn.es/dfsroot/Recursos/metabolismosinaptico/SOFIA/metabolomics_rett_grin/metabolomica/P21145_Results_final.xlsx")),
               read_excel, path = "//hsjdbcn.es/dfsroot/Recursos/metabolismosinaptico/SOFIA/metabolomics_rett_grin/metabolomica/P21145_Results_final.xlsx",
)
library(writexl)


aas <- c("Phenylalanine",
         "Proline",
         "Serine",
         "Threonine",
         "Tryptophan",
         "Valine",
         "Alanine", 
         "Aspartic acid",
         "Glutamine",
         "Glutamic acid",
         "Glycine",
         "Isoleucine",
         "Leucine")

df_list <- lapply(df_list,t)
df_list <- lapply(df_list,data.frame)

df_csf <- df_list[names(df_list) %in% "Plasma_GCMS" == FALSE]
df_csf <- lapply(df_csf,function(x) x[,1:42])

df_all <- rbind(df_csf$"Tryptophan & metabolites",df_csf$"CSF_GCMS")
df_all <- df_all %>%
  row_to_names(row_number = 1)

m <- read.csv("C:/Users/killescas/OneDrive - San Juan de Dios/Desktop/TESIS/LCR/metabolomica/metdat.csv",sep=";")
df_all <- df_all[,sort(colnames(df_all))]

m <- m[order(m$COS.Code),]

colnames(df_all) <- m$ID
metnames <- rownames(df_all)
df_all <- data.frame(t(df_all))
colnames(df_all) <- metnames
df_all <- df_all[order(rownames(df_all)),]
df_all <- subset(df_all,rownames(df_all)%in%metdat$ID)
df_all <- df_all[,aas]
metdat <- metdat[1:24,]
metdat <- metdat[order(metdat$ID),]
tab_aa <- cbind(metdat,df_all)

write_xlsx(tab_aa, path="C:/Users/killescas/OneDrive - San Juan de Dios/Desktop/TESIS/LCR/metabolomica/tab_aa.xlsx")

library(LaCroixColoR)
metdat <- metdat[order(metdat$Disease),]
lcp <- lacroix_palette("Tangerine", 5, type = "discrete")
ggboxplot(metdat, x = "Disease", y = "Age", fill = "Disease", palette = lcp[-2]) + geom_point()