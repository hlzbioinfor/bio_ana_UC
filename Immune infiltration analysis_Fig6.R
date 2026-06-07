library(tidyverse)
library(dplyr)
library(stringr)
library(data.table)
library(tibble)
library(ggpubr)
library(patchwork)
library(paletteer)
library(IOBR)
packageVersion("IOBR")
packageVersion("IOBR")


exp <- data_merge$exp_all
range(exp)
pd<-data_merge$pd_all
identical(colnames(exp),rownames(pd))
# 我们使用tpm数据进行演示，首先进行log2转换，然后去掉一些低质量的基因。
# exp_immu <- exp[apply(exp,1,sd)>0.5,]
# dim(exp_immu)
exp_immu <- exp
identical(rownames(pd),colnames(exp))

tme_deconvolution_methods
```
```{r}
# CIBERSORT
im_cibersort <- deconvo_tme(eset = exp_immu,
                            method = "cibersort",
                            arrays = T,
                            perm = 1000
)

epic_res <- EPIC(
  bulk      = expr_mat,
  reference = "TRef"
)
head(epic_res)

EPIC_re<-epic_res$mRNAProportions

