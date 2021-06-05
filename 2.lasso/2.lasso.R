install.packages('glmnet')
library(glmnet)
phe<-read.csv('phe.csv',quote='',header=T,row.names=1,check.names=F)  #Obtain phenotypic data
data<-read.csv('lasso_input.csv',quote='',header=T,row.names=1,check.names=F)  #Get the number of expressions
train_sub = sample(nrow(data),7/10*nrow(data))
train_data = data[train_sub,]
test_data = data[-train_sub,]
set.seed(1)  
fit_cv <- cv.glmnet(t(train_data), as.matrix(phe$status), alpha=1, family = 'binomial', type.measure='auc') 
plot(fit_cv)  
get_coe <- function(the_fit,the_lamb){  #Function to extract results
Coefficients <- coef(the_fit, s = the_lamb)
Active.Index <- which(Coefficients != 0)  
Active.Coefficients <- Coefficients[Active.Index]
re <- data.frame(rownames(Coefficients)[Active.Index],Active.Coefficients)
re <- data.table('var_names'=rownames(Coefficients)[Active.Index], 'coef'=Active.Coefficients)
re$expcoef <- exp(re$coef)  
return(re[order(expcoef)])
}
get_coe(fit_cv,fit_cv$lambda.min)  #Extract model coefficients
model<-as.data.frame(get_coe(fit_cv,fit_cv$lambda.min)) #Save as data frame
model<-model[-1,]  #Cut off the intercept
model<-model[,-3]  #Remove expcoef
names(model)<-c('id','coef')
data$id <- rownames(data)
model_expr<-merge(model,data,by='id')  
rownames(model_expr)<-model$id  
model_expr<-model_expr[,-1]
sample<-model$coef*model_expr  #Multiply each gene by the corresponding coefficient
sample<-apply(sample,2,sum)  
summary(sample)
write.csv(sample,file = 'trainriskscore.csv')