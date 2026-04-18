library(stringr)
?dir
getwd()
setwd("part4_sc/")
library(stringr)
?dir
getwd()
setwd("part4_sc/")
fs = paste0("GSE231993_RAW/",dir("GSE231993_RAW/"))
fs
# list.files("GSE231993_RAW/")
# samples <- fs %>%
#   str_extract("(young|middle)-[0-9]+") %>%
#   unique();samples
samples = dir("GSE231993_RAW/") %>% str_split_i(pattern = "_",i = 1) %>% unique();samples
#多个样本，每个样本单独一个文件夹，每个文件夹里3个固定名称的文件
ctr = function(s){
  ns = paste0("01_data/",s)
  if(!file.exists(ns))dir.create(ns,recursive = T)
}
lapply(samples, ctr)

lapply(fs, function(s){
  #s = fs[1]
  for(i in 1:length(samples)){
    #i = 1
    if(str_detect(s,samples[[i]])){
      file.copy(s,paste0("01_data/",samples[[i]]))
    }
  }
})

#所有文件改名，去掉前缀（工作目录之下开始写）
on = paste0("01_data/",dir("01_data/",recursive = T));on
# nn = str_remove(on,"GSM\\d+_sample\\d_");nn
# 使用正则表达式提取第三个短横线之后的内容
# nn= extracted <- str_extract(on, "(?<=-)[^-]+$")
# nn = str_remove(on,"GSM\\d+_\\d+[A-Z]-[A-Z]-\\d+-\\d+-");nn

nn <- str_remove(on, "GSM\\d+_\\d+[A-Z]-[A-Z]-.*-");nn
# nn <- str_remove(on, "GSM/([^/]+)\\.(barcode|matrix|feature)");nn
file.rename(on,nn)
?file.rename()




#所有文件改名，去掉前缀（工作目录之下开始写）
on = paste0("01_data/",dir("01_data/",recursive = T));on
# nn = str_remove(on,"GSM\\d+_sample\\d_");nn
# 使用正则表达式提取第三个短横线之后的内容
# nn= extracted <- str_extract(on, "(?<=-)[^-]+$")
# nn = str_remove(on,"GSM\\d+_\\d+[A-Z]-[A-Z]-\\d+-\\d+-");nn

nn <- str_remove(on, "GSM\\d+_\\d+[A-Z]-[A-Z]-.*-");nn
# nn <- str_remove(on, "GSM/([^/]+)\\.(barcode|matrix|feature)");nn
file.rename(on,nn)
?file.rename()
#2.批量读取
rm(list = ls())
library(Seurat)
rdaf = "sce.all.Rdata"
if(!file.exists(rdaf)){
  f = dir("01_data/")
  scelist = list() #创建空的列表，下面的for循环每执行一次，scelist里面就会多一个元素。
  for(i in 1:length(f)){
    pda <- Read10X(paste0("01_data/",f[[i]]))
    scelist[[i]] <- CreateSeuratObject(counts = pda, 
                                       project = f[[i]],
                                       min.cells = 3,
                                       min.features = 200)
    print(dim(scelist[[i]]))#输出每个文件的基因数和细胞数
  }
  sce.all = merge(scelist[[1]],scelist[-1]) #合并多个对象
  sce.all = JoinLayers(sce.all) 
  #merge后，每个样本的表达矩阵是一个独立的的layer，JoinLayers是合并为一个表达矩阵
  set.seed(335)
  # sce.all = subset(sce.all,downsample=700)#每个样本抽700个细胞
  save(sce.all,file = rdaf)
}
table(sce.all$orig.ident)
sum(table(Idents(sce.all)))



#3.质控指标
#过滤线粒体、核糖体和红细胞基因
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")



#根据小提琴图指定指标去掉离群值（并没有用核糖体和红细胞基因过滤）
sce.all = subset(sce.all,percent.mt < 25&
                   nCount_RNA < 40000 &
                   nFeature_RNA < 6000)
table(sce.all@meta.data$orig.ident)

SaveSeuratRds(sce.all,file = "UC_used.rds")
getwd()
ncol(sce.all)


#4.整合降维聚类分群
#多样本的整合，使用harmony
f = "obj.Rdata"
library(harmony)
if(!file.exists(f)){
  sce.all = sce.all %>% 
    NormalizeData() %>%  
    FindVariableFeatures() %>%  
    ScaleData(features = rownames(.)) %>%  
    RunPCA(pc.genes = VariableFeatures(.))  %>%
    RunHarmony("orig.ident") %>%
    FindNeighbors(dims = 1:15, reduction = "harmony") %>% 
    FindClusters(resolution = 0.5) %>% 
    RunUMAP(dims = 1:15,reduction = "harmony") %>% 
    RunTSNE(dims = 1:15,reduction = "harmony")
  save(sce.all,file = f)
}
load(f)
ElbowPlot(sce.all)


# 将样本、细胞类型、细胞个数制成表格
sample_table = as.data.frame(table(scobj@meta.data$group,
                                   scobj@meta.data$celltype))
# 定义列名
names(sample_table ) = c("Samples","Celltype","CellNumber")
### 横轴为细胞类型，纵轴为不同样本来源细胞比例
# 将样本、细胞类型、细胞个数制成表格
sample_table2 = as.data.frame(table(scobj@meta.data$celltype,
                                    scobj@meta.data$group))
# 定义列名
names(sample_table2) = c("celltype","Samples","CellNumber")

# 指定颜色
sample_colors = c("#FF6F61", "#6B5B95", "#70C1B3", "#FFBE0B", "#FF85A1", "#A05195")
# 横轴为样本,纵轴为细胞比例
p1 = ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=Celltype)) +  
  geom_bar(position="fill",width = 0.7, linewidth = 0.5,colour = '#222222') + 
  scale_fill_manual(values = colors) +  
  theme(panel.grid = element_blank(),        
        panel.background = element_rect(fill = "transparent",colour = NA),        
        axis.line.x = element_line(colour = "black"),       
        axis.line.y = element_line(colour = "black"),      
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size = 15)) + 
  labs(y="Percentage") + 
  RotatedAxis()
p1

## 横轴为细胞比例，纵轴为样本
p2 = ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=Celltype)) + 
  geom_bar(position="fill",width = 0.7, inewidth = 0.5,colour = '#222222') + 
  scale_fill_manual(values = colors) +   
  theme(panel.grid = element_blank(),    
        panel.background = element_rect(fill = "transparent",colour = NA),       
        axis.line.x = element_line(colour = "black"),       
        axis.line.y = element_line(colour = "black"),       
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size = 15)) +
  labs(y="Percentage") + 
  RotatedAxis()+coord_flip() 
p2

library(dplyr)
library(stringr)
if( length(ncell1)>2 && length(ncell2)>2 ) {
  sce.markers <- FindMarkers(sce, 
                             logfc.threshold=0.20, 
                             min.pct=0.05, 
                             test.use = "wilcox",
                             fc.name = "avg_log2FC",
                             ident.1=paste0(celltype[i], "_UC"), 
                             ident.2=paste0(celltype[i], "_control") 
  )
  sce.markers$gene <- rownames(sce.markers)
  sce.markers$Group <- celltype[i]
  sce.markers$Regulated <- if_else(sce.markers$avg_log2FC > 0, "Up", "Down")
} else{
  print(paste0("Only One group or group cells num < 3: ",celltype[i]))
  cat("\n")
  next
}
