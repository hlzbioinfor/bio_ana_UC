library(pROC)
library(ggplot2)
library(caret)
library(ggprism)
library(data.table)
library(dplyr) #数据处理使用
library(data.table) #数据读取使用
library(caTools) #LR模型使用
library(ggpubr) #绘图使用
library(glmnet)  # 加载LASSO回归核心包
# load("merge_phe_exp.rda")
data <- data_merge$exp_all[aim_gens,] %>% as.data.frame()
group <- data_merge$pd_all$group
identical(data %>% colnames(),rownames(data_merge$pd_all))
# data <- data%>% mutate(group = factor(ifelse(group == "Disease",1,0)))
data$group %>% table
data <- t(data) %>% as.data.frame()
data$group <- group
data <- data[,c(ncol(data),1:(ncol(data)-1))]
ncol(data)

split <- sample.split(data$group, SplitRatio = 0.7)  # 将数据按照指定比例分割
train_data <- subset(data, split == TRUE)  #
test_data <- subset(data, split == FALSE)  #
x <- model_data[,-1] %>% as.matrix()
y <- model_data[,1] %>% as.matrix()
# model <- glm(cluster ~ ., data = train_data, family = binomial)
# 使用交叉验证确定最优lambda值
cv_model <- cv.glmnet(x,y,
                      alpha =1,                  # alpha=1表示Lasso回归
                      family ="binomial",        # 二分类问题
                      type.measure ="deviance",# 使用偏差作为评估指标
                      nfolds =10)                # 10折交叉验证



# 保存交叉验证结果图
# 可视化交叉验证结果（选择最优Lambda）
plot(cv_model, main = "交叉验证误差曲线")  
# setwd("")
plot(cv_model)
# pdf(file.path("cv_lambda.pdf"), width =10, height =8)
# plot(cv_model)
# dev.off()
# dev.new()
# 输出最优lambda值
cat("最优lambda (lambda.min):", cv_model$lambda.min,"\n")
cat("1个标准误的lambda (lambda.1se):", cv_model$lambda.1se,"\n")
####
# 训练最终模型（使用全数据与最优Lambda）
lasso <- glmnet(
  x = x,
  y = y,
  family = "binomial",
  alpha = 1,
  lambda = cv_model$lambda # 使用交叉验证中的Lambda序列
)
# 绘制系数路径图（观察特征随Lambda变化）
plot(lasso, xvar = "lambda", label = TRUE, main = "LASSO系数收缩路径")
plot(lasso, xvar = "lambda", label = F, main = "LASSO系数收缩路径")
# 使用lambda.1se构建最终模型（更简约的模型）
lasso_model <- glmnet(x, y,
                      alpha =1,
                      family ="binomial",
                      lambda = cv_model$lambda.1se)

# lasso_model_1se <- glmnet(x, y,
#                       alpha =1,
#                       family ="binomial",
#                       lambda = cv_model$lambda.1se)




# 查看系数
coef_lasso <- coef(lasso_model)
print(coef_lasso)
# coef_lasso_1se <- coef(lasso_model_1se)
# selected_vars_1se<- rownames(coef_lasso_1se)[which(coef_lasso_1se !=0)]
# 找出被选中的变量（系数不为零）
selected_vars <- rownames(coef_lasso)[which(coef_lasso !=0)]
selected_vars <- selected_vars[selected_vars !="(Intercept)"]# 移除截距项
cat("Lasso选中的变量:", paste(selected_vars, collapse =", "),"\n")

# 保存选中的变量及其系数
selected_coef <- data.frame(
  Variable = rownames(coef_lasso)[which(coef_lasso !=0)],
  Coefficient =as.numeric(coef_lasso[which(coef_lasso !=0)])
)
length(selected_vars)
write.csv(selected_coef,file = "lasso_model.csv")
getwd()
# setwd("part3_bulk_ana/lasso/")

# 导入必要的包,没有安装的可以先安装一下
library(dplyr) #数据处理使用
library(data.table) #数据读取使用
library(randomForest) #RF模型使用
packageVersion("randomForest")
library(caret) # 调参和计算模型评价参数使用
library(pROC) #绘图使用
library(ggplot2) #绘图使用
library(ggpubr) #绘图使用
library(ggprism) #绘图使用
# 导入必要的包,没有安装的可以先安装一下
library(dplyr) #数据处理使用
library(data.table) #数据读取使用
library(randomForest) #RF模型使用
library(caret) # 调参和计算模型评价参数使用
library(pROC) #绘图使用
library(ggplot2) #绘图使用
library(ggpubr) #绘图使用
library(ggprism) #绘图使用

# 读取数据
load("merge_phe_exp.rda")
merge_phe_exp <- merge_phe_exp%>% mutate(group = factor(ifelse(group == "Disease",1,0)))
data <- merge_phe_exp  # 替换为你的数据文件名或路径# 分割数据为训练集和测试集
set.seed(123)  # 设置随机种子，保证结果可复现
# split <- sample.split(data$type, SplitRatio = 0.8)  # 将数据按照指定比例分割
# train_data <- subset(data, split == TRUE)  # 训练集
# test_data <- subset(data, split == FALSE)  # 测试集
train_data <- data# 训练集
test_data <-data  # 测试集
# 定义训练集特征和目标变量
X_train <- train_data[, -1]
y_train <- as.factor(train_data[, 1])

# 创建随机森林分类模型
fit <- randomForest(x = X_train, y = y_train,ntree = 300)
plot(fit)
print(fit)
# 输出默认参数下的模型性能
# print(model)


# 进行参数调优
# 创建训练控制对象
ctrl <- trainControl(method = "cv", number = 5)#使用五折交叉验证，也可以选择10折交叉验证。
# 定义参数网格
grid <- expand.grid(mtry = c(2:8))  # 每棵树中用于分裂的特征数量，这里只是随便给的测试，主要为了介绍如何调参，并非最优选择。

# 使用caret包进行调参
set.seed(12345) 
rf_model <- train(x = X_train, y = y_train,
                  method = "rf",
                  trControl = ctrl,
                  tuneGrid = grid)

# 输出最佳模型和参数
print(rf_model)
plot(rf_model)
####
# setwd("../rf/")
# 调整Caret没有提供的参数
# 如果我们想调整的参数Caret没有提供，可以用下面的方式自己手动调参。
# 用刚刚调参的最佳mtry值固定mtry
grid <- expand.grid(mtry = c(4))  # 每棵树中用于分裂的特征数量

# 定义模型列表，存储每一个模型评估结果
modellist <- list()

# 调整的参数是决策树的数量
for (ntree in c(100,200, 300,400,500)) {
  set.seed(123)
  fit <- train(x = X_train, y = y_train, method="rf", 
               metric="Accuracy", tuneGrid=grid, 
               trControl=ctrl, ntree=ntree)
  key <- toString(ntree)
  modellist[[key]] <- fit
}

# compare results
results <- resamples(modellist)
# 输出最佳模型和参数
summary(results)

# 使用最佳参数训练最终模型
final_model <- randomForest(x = X_train, y = y_train,mtry = 4,ntree = 100)
# 输出最终模型
print(final_model)
plot(final_model)
final_model$importance

varImpPlot(final_model)



# set.seed(647)
res <- rfcv(trainx = dat.cox[,-1],trainy = dat.cox[,1],
            cv.fold = 5,
            recursive = T
)
res$n.var #变量个数












# 在测试集上进行预测
X_test <- test_data[, -1]
y_test <- as.factor(test_data[, 1])
test_predictions <- predict(final_model, newdata = test_data)

# 计算模型指标
confusion_matrix <- confusionMatrix(test_predictions, y_test)
accuracy <- confusion_matrix$overall["Accuracy"]
precision <- confusion_matrix$byClass["Pos Pred Value"]
recall <- confusion_matrix$byClass["Sensitivity"]
f1_score <- confusion_matrix$byClass["F1"]

# 输出模型指标
print(confusion_matrix)
print(paste("Accuracy:", accuracy))
print(paste("Precision:", precision))
print(paste("Recall:", recall))
print(paste("F1 Score:", f1_score))


load("../lasso/lasso_ana.rds")
# 绘制混淆矩阵热图
# 将混淆矩阵转换为数据框

confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Normal","Disease")
rownames(confusion_matrix_df) <- c("Normal","Disease")
draw_data <- round(confusion_matrix_df / rowSums(confusion_matrix_df),2)
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)

ggplot(draw_data, aes(real,variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(value))) +
  scale_fill_gradient(low = "#F0F0F0", high = "#3575b5") +
  labs(x = "True", y = "Guess", title = "Confusion matrix") +
  theme_prism(border = T)+
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none")

final_model <- randomForest(x = X_train, y = y_train,mtry = 4,ntree = 100)
write.csv(mach_gene$rf_genes %>% as.data.frame(),file = "revise/补充材料/step4machine_learning/rf_genes.csv",)
###ROC
library(tidyverse)
library(dplyr)
library(pROC)
colnames(t_roc_data)[ncol(t_roc_data)] = "Group"
colnames(t_roc_data)
genes = colnames(t_roc_data)[1:8]
genes
roc_formula <- paste("Group ~", paste(genes, collapse = "+")) %>% as.formula
roc_formula
data = t_roc_data
roc <- roc(roc_formula, data = data,
           percent = F, 
           ci= TRUE, 
           smoooth = F, 
           levels=c("Normal", "Disease")
)

roc
### 绘制多条ROC曲线
auc_labels <- sapply(names(roc), function(gene) {
  sprintf("%s AUC=%.3f", gene, roc[[gene]]$auc)
})
# 拼接
subtitle_text <- paste(auc_labels, collapse = "\n")

ggroc(roc, legacy.axes = TRUE )+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=1, linewidth = 0.4)+
  ggtitle('')+
  ggsci::scale_color_lancet()+
  theme_test()+
  ggprism::theme_prism(border = T)+
  annotate("text",x=0.75,y=0.05,label= subtitle_text )
ggroc(roc, legacy.axes = TRUE )+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="darkgrey", linetype=1, linewidth = 0.4)+
  ggtitle('')+
  ggsci::scale_color_lancet()+
  theme_test()+
  ggprism::theme_prism(border = T)+
  annotate("text",x=0.75,y=0.22,label= subtitle_text )
ls()
roc_data <- list(
  data= data,
  roc= roc,
  t_roc_data= t_roc_data
)
save(roc_data,file = "roc_new_vali.rda")
