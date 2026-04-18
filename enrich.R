####
library(gground)
library(ggprism)
library(tidyverse)
library(org.Hs.eg.db)
library(clusterProfiler)
library(package)
hvgs <- DEG_limma_voom %>% filter(change !="NOT")
nrow(hvgs)
toxic_targets <- fread("toxic_targets_merge/toxic_targets.csv",header = T)

toxic_targets <- toxic_targets %>% select(x)
toxic_pcos_genes <- intersect(rownames(hvgs),toxic_targets$x)
toxi
rownames(hvgs)
save(toxic_pcos_genes,hvgs,file="deg_intersect.rda")


# data_fdr<- data.table::fread('gens_fdr.csv',
#                           
# colnames(data_fdr)[1]<- "name"
# deg <- data_fdr
# symbol 转entrzid
setwd("enrich/")

df <- bitr(unique(library(gground)
                  library(ggprism)
                  library(tidyverse)
                  library(org.Hs.eg.db)
                  library(clusterProfiler)
                  library(package)
                  # data_fdr<- data.table::fread('gens_fdr.csv',
                  #                           
                  # colnames(data_fdr)[1]<- "name"
                  # deg <- data_fdr
                  # symbol 转entrzid
                  
                  gene_enrich <- intersect(deg_sig$symbol,all__toxic_targets)
                  length(gene_enrich)
                  df <- bitr(unique(gene_enrich), fromType = "SYMBOL",
                             toType = c( "ENTREZID"),
                             OrgDb = org.Hs.eg.db)
                  
                  
GO_all<-enrichGO(gene=df$ENTREZID,#基因列表(转换的ID)
                 keyType="ENTREZID",#指定的基因ID类型，默认为ENTREZID
                 OrgDb=org.Hs.eg.db,#物种对应的org包
                 ont="ALL",#CC细胞组件，MF分子功能，BF生物学过程，ALL以上三个
                 pvalueCutoff=0.05,#p值阈值,一般设置为0.05如果报错提示此阈值下的基因过少可以将此处的阈值调高（报错信息见下）  
                 pAdjustMethod="fdr",#多重假设检验校正方式
                 minGSSize=10,#注释的最小基因集，默认为10
                 maxGSSize=500,#注释的最大基因集，默认为500
                 qvalueCutoff=0.05,#q值阈值，一般设置为0.05如果报错提示此阈值下的基因过少如果报错可以将此处的阈值调高（报错信息见下）
                 readable=TRUE)#基因ID转换为基因名

GO_result<-data.frame(GO_all)



kk <- enrichKEGG(gene = df$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

kk <- DOSE::setReadable(kk, OrgDb='org.Hs.eg.db', keyType='ENTREZID')#按需替换
barplot(kk)

kegg_result <- kk@result
View(kegg_result)
head(kk)[,1:6]


write.csv(kegg_result,file = "kegg_result.csv")
# write.csv(kegg_result,file = "kegg_result.csv")
# 加载必要的库
library(ggplot2)
library(ggrepel)  # 用于更好的文本标签布局
library(ggforce)
# 查看数据
head(dt)

# 定义颜色
colors <- c(
  "BP" = "#FFB6C1",
  "CC" = "#ADD8E6",   
  "MF" = "#E6E6FA")
# )
colnames(GO_result)
# 创建基础绘图
# dt<- dt %>% arrange(ONTOLOGY)
p <- ggplot(dt) +
  geom_link(aes(x = 0, y = Description,
                xend = -log10(pvalue), yend = Description,
                # alpha = after_stat(),
                color = ONTOLOGY,
                size = after_stat(index)),
            n = 500, show.legend = TRUE) 
p
colnames(dt)
# 添加点和文本标签
p1 <- p +
  geom_point(aes(x = -log10(pvalue),y = Description), color = "black", fill = "white",size = 6,shape = 21) +
  geom_text(aes(x = -log10(pvalue), y = Description), label=dt$Count, size=3, nudge_x=0.05) + 
  theme_classic() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        #axis.text = element_text(color = "black"),
        axis.line = element_line(color = "black", size = 0.6), # 加粗x轴和y轴的线条
        axis.text = element_text(face = "bold"), # 加粗x轴和y轴的标签
        axis.title = element_text( size = 13)    # 加粗x轴和y轴的标题
  ) +
  xlab("-Log10 Pvalue") + ylab("") + 
  scale_color_manual(values = colors)
p1

p2 <- p1 + 
  facet_wrap(~ONTOLOGY)
# ,scales = "free")
# ,ncol = 2) 
p2

