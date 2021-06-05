library(survival)

data = read.table("train_roc.txt",sep = "\t", header = T, row.names = 1,check.names = F)

cox = coxph(Surv(time, status) ~ ., data)## Analyze the correlation between all indicators and prognosis

coxsummary = summary(cox)
genes = rownames(coxsummary$coefficients)
outTab = cbind(Gene = genes, HR = round(coxsummary$coefficients[,"exp(coef)"],2),
               z=round(coxsummary$coefficients[,"z"],2),
               "95% CI"=paste(round(coxsummary$conf.int[,3], 2), round(coxsummary$conf.int[,4],2), sep = "-"),
               pvalue=round(coxsummary$coefficients[,"Pr(>|z|)"],2))
write.table(outTab,"multivariateCox.result.txt", row.names = F, sep = "\t", quote = F)
