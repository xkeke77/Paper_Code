install.packages('glmnet')
library(glmnet)
phe<-read.csv('phe.csv',quote='',header=T,row.names=1,check.names=F)  #获取表型数据
data<-read.csv('cox_exprdata.csv',quote='',header=T,row.names=1,check.names=F)  #获取表达数
set.seed(1)  #设置种子，为了可以复现结果
fit_cv <- cv.glmnet(t(data), as.matrix(phe$event), alpha=1, family = 'binomial', type.measure='auc') #交叉验证，多次建模，找到最合适的参数lambda
plot(fit_cv)  #绘制交叉验证结果
get_coe <- function(the_fit,the_lamb){  #提取结果的函数
Coefficients <- coef(the_fit, s = the_lamb)
Active.Index <- which(Coefficients != 0)  #至少coef(回归系数)不能为0啊
Active.Coefficients <- Coefficients[Active.Index]
re <- data.frame(rownames(Coefficients)[Active.Index],Active.Coefficients)
re <- data.table('var_names'=rownames(Coefficients)[Active.Index], 'coef'=Active.Coefficients)
re$expcoef <- exp(re$coef)  
return(re[order(expcoef)])
}
get_coe(fit_cv,fit_cv$lambda.min)  #提取模型系数
model<-as.data.frame(get_coe(fit_cv,fit_cv$lambda.min)) #存为数据框
model<-model[-1,]  #去掉截距，你要看好截距在哪一行啊
model<-model[,-3]  #去掉expcoef
names(model)<-c('id','coef') #重新命名
data$id <- rownames(data)
model_expr<-merge(model,data,by='id')  #整合表达
rownames(model_expr)<-model$id  #设置行名
model_expr<-model_expr[,-1]
sample<-model$coef*model_expr  #对每一个基因乘以对应的系数
sample<-apply(sample,2,sum)  #对列求和，2是列，1是行
summary(sample)
write.csv(sample,file = 'riskscore.csv')