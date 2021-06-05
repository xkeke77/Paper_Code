## 加载survival包
library(survival)

## 输入数据与单因素Cox分析类似
data = read.table("示例数据.txt",sep = "\t", header = T, row.names = 1,check.names = F)

cox = coxph(Surv(time, status) ~ ., data)## 对全部数据分析与预后的相关性

## 输出分析结果，与单因素Cox相似
coxsummary = summary(cox)
genes = rownames(coxsummary$coefficients)
outTab = cbind(Gene = genes, HR = round(coxsummary$coefficients[,"exp(coef)"],2),
               z=round(coxsummary$coefficients[,"z"],2),
               "95% CI"=paste(round(coxsummary$conf.int[,3], 2), round(coxsummary$conf.int[,4],2), sep = "-"),
               pvalue=round(coxsummary$coefficients[,"Pr(>|z|)"],2))
write.table(outTab,"multivariateCox.result.txt", row.names = F, sep = "\t", quote = F)




cox = step(cox, direction = "both") ## 筛选最合适的模型

## 输出分析结果，与单因素Cox相似
coxsummary = summary(cox)

genes = rownames(coxsummary$coefficients)
outTab = cbind(Gene = genes, HR = round(coxsummary$coefficients[,"exp(coef)"],2),
               z=round(coxsummary$coefficients[,"z"],2),
               "95% CI"=paste(round(coxsummary$conf.int[,3], 2), round(coxsummary$conf.int[,4],2), sep = "-"),
               pvalue=round(coxsummary$coefficients[,"Pr(>|z|)"],2))
write.table(outTab,"risk.model.txt", row.names = F, sep = "\t", quote = F)

## 计算风险评分，并将样品按风险评分分成高低风险组
riskscore = predict(cox, type = "risk", newdata = data)


risk = ifelse(riskscore>median(riskscore),"high","low")
riskresult = cbind(patient = rownames(data),data[,c("time","status",genes)], riskscore, risk)
write.table(riskresult,"riskscore.txt", sep = "\t", quote = F, row.names = F)
