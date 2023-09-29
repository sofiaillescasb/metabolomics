
library(purrr)
library(readxl)
library(janitor)
library(imputeTS)
library(reshape2)
library(writexl)
library(qwraps2)
library(tidyverse)
library(here)

metdat <- read.csv(here("metadata/metdat.csv"), header=TRUE,row.names = 1,sep=";")

# metdat <- read.csv("/Users/sofiaillescas/Desktop/LAB/LCR/metabolomica/metdat.csv", header=TRUE, row.names = 1)

metdat_char <- read.csv(here("metadata/patient_characteristics.csv"),sep=";")

df_list <- map(purrr::set_names(excel_sheets(here("data_raw/P21145_Results_final.xlsx"))),
               read_excel, path = here("data_raw/P21145_Results_final.xlsx"),
)

# df_list <- map(purrr::set_names(excel_sheets("/Users/sofiaillescas/Desktop/LAB/LCR/metabolomica/P21145_Results_final.xlsx")),
#   read_excel, path = "/Users/sofiaillescas/Desktop/LAB/LCR/metabolomica/P21145_Results_final.xlsx",)

df_list <- lapply(df_list,t)
df_list <- lapply(df_list,data.frame)


df_csf <- df_list[names(df_list) %in% "Plasma_GCMS" == FALSE]
df_csf <- lapply(df_csf,function(x) x[,1:41])

df_all <- rbind(df_csf$"Tryptophan & metabolites",df_csf$"CSF_GCMS")
df_all <- df_all %>%
  row_to_names(row_number = 1)
df_all <- subset(df_all, rownames(df_all)!="COS Code1")
# df_all <- subset(df_all, rownames(df_all)!="Sample ID1")
# df_all <- subset(df_all, rownames(df_all)!="Sample ID")
df_all <- subset(df_all,rowSums(is.na(df_all))!=ncol(df_all))
df_all_num <- apply(df_all, c(1,2), as.numeric)
rows_to_remove <- grep("Carn", rownames(df_all_num))
df_all_num <- df_all_num[-rows_to_remove,]
mins <- apply(df_all_num,1,min,na.rm=TRUE)
mins <- mins/2
df_all <- data.frame(na_replace(df_all_num))
df_all <- df_all+mins
colnames(df_all) <- gsub("\\."," ",colnames(df_all))
metdat <- subset(metdat, COS.Code != "CSF.18")

df_all_log <- log2(df_all)
metdat <- metdat[order(metdat$ID),]
df_all_log <- df_all_log[,order(colnames(df_all_log))]
df_norm_lab <- data.frame(t(df_all_log))
df_norm_lab <- cbind(metdat$Label,df_norm_lab)
colnames(df_norm_lab)[1] <- "Label" 
colnames(df_norm_lab) <- c("Label",rownames(df_all_log))
df_norm_lab <- subset(df_norm_lab,Label%in%c(3,8,10,12,2))
df_norm_lab$Label <- gsub(3,1,df_norm_lab$Label)
df_norm_lab1 <- rownames_to_column(df_norm_lab)

log_df <- melt(df_norm_lab1, id.vars = c("rowname","Label"))
log_df$Label <- as.factor(log_df$Label)
all_mets <- data.frame(colnames(df_norm_lab)[-1])
colnames(all_mets) <- "Detected metabolites"



write_xlsx(all_mets, path=here("data_processed/supplementary_1.xlsx"))
saveRDS(log_df, here("data_processed/long_format_df.rds"))
saveRDS(df_norm_lab, here("data_processed/wide_format_df.rds"))

df_all_2 <- data.frame(cbind(metdat$Label,t(df_all)),check.names = FALSE)
colnames(df_all_2)[1] <- "Label" 
# colnames(df_all_2) <- c("Label",rownames(df_all_2))
df_all_2 <- subset(df_all_2,Label%in%c(3,8,10,12,2))
df_all_2$Label <- gsub(12,8, df_all_2$Label)
df_all_2$Label <- gsub(2,10, df_all_2$Label)
df_all_2 <- rownames_to_column(df_all_2 ,var = "ID")



library(qwraps2)

# For continuous variables raw

vars <- colnames(df_all_2)[-c(1,2)]
data <- df_all_2

summary_list <- list()

# Iterate over the variable names
for (l in vars) {
  summary_list[[l]] <- 
    list("min" = as.formula(paste0("~ min(`", l, "`)")),
         "median (iqr)" = as.formula(paste0("~ median_iqr(`", l, "`)")),
         "max" = as.formula(paste0("~ max(`", l, "`)")),
         "mean (CI)" = as.formula(paste0("~ qwraps2::frmtci(mean_ci(`", l, "`,show_level=TRUE))")),
         "SD" = as.formula(paste0("~ sd(`", l, "`)")))
  
}

st_lst <- summary_table(dplyr::group_by(data,Label),summary_list)
rnms <- rownames(st_lst)
st_lst <- data.frame(cbind(rep(vars,each=5),st_lst), check.names = FALSE)
colnames(st_lst)[1] <- "Metabolite"
st_lst <- rownames_to_column(st_lst, var="Statistic")
st_lst$Statistic <- rnms
pivoted_data <- pivot_wider(st_lst, names_from = Metabolite, values_from = colnames(st_lst)[-c(1,2)])
