setwd("/Users/sofiaillescas/Desktop/LAB/zandl-lang2022/lda-proce")

library(lipidr)
library(purrr)
library(readxl)

filenames <- list.files(".", pattern="*.xlsx", full.names=TRUE)

# df_list <- lapply(filenames, function(x) map(set_names(excel_sheets(x)),
#                read_excel, path = x))


df_list <- readRDS("df_list")
names(df_list) <- gsub("./","",filenames)
names(df_list) <- gsub("shortlipids","",names(df_list))
names(df_list) <- gsub("shortLipids","",names(df_list))
names(df_list) <- gsub(".xlsx","",names(df_list))
saveRDS(df_list, file="df_list")

ctrl <- as.character(c(42,43,44,46,47,49,50,51,52,53,54,55,56))
rett <- as.character(c(35,36,37,38,39,40,41))
df_list <- df_list[names(df_list)%in%paste0(ctrl,"_pos")|names(df_list)%in%paste0(ctrl,"_neg")|names(df_list)%in%paste0(rett,"_pos")|names(df_list)%in%paste0(rett,"_neg")]
#For some reason 46 has its pos and neg names inverted, though that shouldn't really matter because I will add both
df_over <- lapply(names(df_list), function(x) df_list[[x]][grep("Overview",names(df_list[[x]]))])
names(df_over) <- names(df_list)

for (i in 1:length(df_over)) {
  names(df_over[[i]]) <- gsub("- Overview","",names(df_over[[i]]))
  for (j in 1:length(df_over[[i]])) {
    df_over[[i]][[j]] <- df_over[[i]][[j]][,c(1:2)]
    df_over[[i]][[j]]$Species<- substring(df_over[[i]][[j]]$Species,1,4)
    colnames(df_over[[i]][[j]]) <- c("Species","Area")
  }
}

#Adding up the areas of each species           
for (k in 1:length(df_over)) {
  for (l in 1:length(df_over[[k]])) {
    if (length(rownames(df_over[[k]][[l]]))!=0) {
      df_over[[k]][[l]] <- aggregate(Area~Species, data=df_over[[k]][[l]], FUN=sum)
    }
  }
}

#Making dataframe with all the species of each sample together

for (m in 1:length(df_over)) {
  for (n in 1:length(df_over[[m]])) {
    if (length(rownames(df_over[[m]][[n]]))!=0) {
      df_over[[m]][[n]]$Species <- paste0(rep(names(df_over[[m]][n]),length(rownames(df_over[[m]][[n]]))),df_over[[m]][[n]]$Species)
    }
    
  }
}

list_rbind(df_over$`35_neg`)
df_over <- lapply(df_over, function(x) list_rbind(x))
