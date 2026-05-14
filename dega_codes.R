library(limma)
library(tidyverse)
library(dplyr)
library(stringr)
library(data.table)
library(tibble)
library(ggpubr)
library(patchwork)
group <- data_merge$pd_all$group
group %>% table
library(limma)
exp <- data_merge$exp_all
exprSet <- exp%>% as.data.frame()

# group <-Group
group_list <- group
group %>% table
group_list <- factor(group_list,levels = c("Normal","Disease"))
group_list
design <- model.matrix(~0 + group_list)
design

# design <- design %>% as.data.frame()
colnames(design) <- c("Normal","Disease")
rownames(design) <- colnames(exprSet)
design


# 明确 uc vs normal
contrast.matrix <- makeContrasts(Disease - Normal, levels = design)
contrast.matrix

is.na(exprSet) %>% sum()



fit <- lmFit(exprSet, design)
fit <- contrasts.fit(fit = fit, contrast.matrix)
fit <- eBayes(fit)
allDiff <- topTable(fit, number = Inf)

head(allDiff, 12)

write.csv(allDiff,file = "allDiff_lihc.csv")
DEG_limma_voom <- allDiff
logFC = 0.5
k1 <- (DEG_limma_voom$adj.P.Val < 0.05) & (DEG_limma_voom$logFC < -logFC)
k2 <- (DEG_limma_voom$adj.P.Val < 0.05) & (DEG_limma_voom$logFC > logFC)
DEG_limma_voom <- mutate(DEG_limma_voom, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))



DEG_limma_voom$symbol <- rownames(DEG_limma_voom)
head(DEG_limma_voom)
table(DEG_limma_voom$change)



colnames(DEG_limma_voom)
# 普通火山图
p <- ggplot(data = DEG_limma_voom, 
            aes(x = logFC, 
                y = -log10(P.Value))) +  # 设置x轴为logFC，y轴为-P.Value的对数值
  geom_point(alpha = 0.4, size = 1, 
             aes(color = change)) +      # 添加散点，根据change列着色
  ylab("-log10(Pvalue)") +               # 设置y轴标签
  scale_color_manual(values = mypal) +  # 设置颜色映射，蓝色表示下调，灰色表示稳定，红色表示上调
  geom_vline(xintercept = c(-logFC, logFC), lty = 4, col = "black", lwd = 0.8) +  # 添加垂直参考线，用于标记logFC阈值
  geom_hline(yintercept = -log10(P.Value), lty = 4, col = "black", lwd = 0.8) +   # 添加水平参考线，用于标记-P.Value阈值
  theme_classic2()  # 使用网格白底主题
p
