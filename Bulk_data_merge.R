library(tidyverse)
library(readxl)
library(tidyverse)
library(GEOquery)
library(tidyverse)
library(GEOquery)
library(limma) 
library(affy)
library(stringr)
# gse752 <- load(file = "GSE138518_pro/gse138.rds")
gse752 <- data_processed
# rm(gse138)
rm(data_processed)
gse874 <- data_processed
rm(data_processed)
gse752$all_data$exp %>% dim()
gse874$all_data$exp %>% dim()

same_gene <- intersect(gse752$all_data$exp %>% rownames(),gse874$all_data$exp %>% rownames())
length(same_gene)
exp752 <- gse752$all_data$exp
exp874 <- gse874$all_data$exp
# exp752  <- normalizeBetweenArrays(exp138)
# exp193 <- gse193$all_data$exp %>% normalizeBetweenArrays()
boxplot(exp752 )
boxplot(exp874)

exp752  <- exp752[same_gene,] %>% as.data.frame()
exp874 <- exp874[same_gene,] %>% as.data.frame()

# exp_all <-  merge(exp45, exp255, by = "Row.names")

identical(rownames(exp752),rownames(exp874))

# 合并数据
exp_all <- cbind(exp752, exp874)




pd752<- gse752$all_data$pd %>% select(title,`disease:ch1`,group)
pd752$batch <- "gse752"
pd874 <- gse874$all_data$pd %>% select(title,`disease:ch1`,group)
pd874$batch <- "gse874"
pd_all  <- rbind(pd752,pd874)





data_all = list(exp_all = exp_all,
                pd_data = pd_all)

save(data_all,file = "UC_data.rda")
getwd()
setwd("bulk_data/")
####
load("UC_data.rda")
exp_all = data_all$exp_all
pd_all = data_all$pd_data
boxplot(exp_all,outline=T, notch=T,col=pd_all, las=2)
library(FactoMineR)
library(factoextra)
#绘制去批次前PCA图
dat.pca <- PCA(as.data.frame(t(exp_all)), graph = FALSE)

pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "points",#仅显示"点(points)"（但不是“文本(text)”）（show "points" only (but not "text")）
                         col.ind = pd_all$batch,
                         palette = c("#00AFBB", "#E7B800","navy"),
                         addEllipses = TRUE, 
                         legend.title = "Groups")
pca_plot

#### 使用sva包计算批次效应
library(sva)
library(tidyverse)
pd_all$batch <- pd_all$batch %>% as.factor()
exp_all_combat <- ComBat(exp_all, batch = pd_all$batch) # batch为批次信息
# 查看去除批次效应的结果如何
library(tinyarray)
group <- pd_all$batch
draw_pca(exp = exp_all_combat, group_list = pd_all$batch)
draw_pca(exp = exp_all_combat, group_list = pd_all$group)
boxplot(exp_all_combat)
?draw_pca
#做校正的PCA分析
pdf(file = "9_PCA_after.pdf",width = 7,height = 6)
pre.pca <- PCA(t(exp_all),graph = FALSE)
fviz_pca_ind(pre.pca,
             geom= "point",
             col.ind = pd_all$Batch,
             addEllipses = TRUE,
             legend.title="Group")
dat.pca <- PCA(t(exp_all_combat),graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "points",#仅显示"点(points)"（但不是“文本(text)”）（show "points" only (but not "text")）
             col.ind = pd_all$batch,
             palette = c("#00AFBB", "#E7B800","navy"),
             addEllipses = TRUE, 
             legend.title = "Groups")
boxplot(exp_all_combat)
dev.off()


data_merge <- list(exp_all = exp_all_combat,
                   pd_all = pd_all)
save(data_merge,file = "UC_merge.rda")
boxplot(data_all)
