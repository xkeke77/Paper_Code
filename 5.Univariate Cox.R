# 加载survival包
library(survival)

# 读取数据，前2列为生存状态和时间，之后为需要分析的各基因表达量
data = read.table("示例数据.txt",sep = "\t", header = T, row.names = 1,check.names = F)
genes = colnames(data)[3:ncol(data)]

# 初始化用于输出结果的表格
outTab = data.frame()
for (i in genes) {
  expr = data[,i]
  cox = coxph(Surv(time, status) ~ expr, data)
  coxsummary = summary(cox)
  outTab=rbind(outTab,cbind(gene=i,HR=round(coxsummary$coefficients[,"exp(coef)"],2), ## HR表示风险比
                            z=round(coxsummary$coefficients[,"z"],2), ## z值
                          "95% CI"=paste(round(coxsummary$conf.int[,3], 2), round(coxsummary$conf.int[,4],2), sep = "-"), #95%置信区间                         
                          pvalue=round(coxsummary$coefficients[,"Pr(>|z|)"],2))) ## p值
}
# 保存输出结果
write.table(outTab,file="univariateCox.result.txt",sep="\t",row.names=F,quote=F,col.names = c("Gene","HR","z","95% CI","pValue"))
