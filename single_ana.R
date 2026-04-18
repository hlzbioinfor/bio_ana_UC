
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
