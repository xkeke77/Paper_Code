install.packages("survival")
library(survival)

# Read the data, the first 2 columns are the survival state and time, used to build the cox model
data = read.table("train_roc.txt",sep = "\t", header = T, row.names = 1,check.names = F)
genes = colnames(data)[3:ncol(data)]


outTab = data.frame()
for (i in genes) {
  expr = data[,i]
  cox = coxph(Surv(time, status) ~ expr, data)
  coxsummary = summary(cox)
  outTab=rbind(outTab,cbind(gene=i,HR=round(coxsummary$coefficients[,"exp(coef)"],2), 
                            z=round(coxsummary$coefficients[,"z"],2), 
                          "95% CI"=paste(round(coxsummary$conf.int[,3], 2), round(coxsummary$conf.int[,4],2), sep = "-"),                        
                          pvalue=round(coxsummary$coefficients[,"Pr(>|z|)"],2))) 
}
write.table(outTab,file="univariateCox.result.txt",sep="\t",row.names=F,quote=F,col.names = c("Gene","HR","z","95% CI","pValue"))
