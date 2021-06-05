#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")
#biocLite("impute")

logFoldChange=1
adjustP=0.05

library(limma)
library("impute")
setwd("C:\\Users\\lenovo\\Desktop\\1")
rt=read.table("input.txt",sep="\t",header=T)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
#Making expression gene matrix
dimnames=list(rownames(exp),colnames(exp))
exp=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

#impute missing expression data
mat=impute.knn(exp)
rt=mat$data

rt=avereps(rt)     #Gene corresponds to multiple probes to take the mean
#normalize
pdf(file="rawBox.pdf")
boxplot(rt,col = "blue",xaxt = "n",outline = F)
dev.off()
rt=normalizeBetweenArrays(as.matrix(rt))
pdf(file="normalBox.pdf")
boxplot(rt,col = "red",xaxt = "n",outline = F)
dev.off()


#differential
#class <- c("con","con","con","treat","treat","treat")
class <- c(rep("con",16),rep("treat",6))    #Modify it according to the grouping
design <- model.matrix(~0+factor(class))
colnames(design) <- c("con","treat")
#Linear model fitting
fit <- lmFit(rt,design)
#Build a difference comparison matrix
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
#Bayes test
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)

#write table
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & P.Value < adjustP )), ]
write.table(diffSig,file="diff.xls",sep="\t",quote=F)
diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & P.Value < adjustP )), ]
write.table(diffUp,file="up.xls",sep="\t",quote=F)
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & P.Value < adjustP )), ]
write.table(diffDown,file="down.xls",sep="\t",quote=F)

#write expression level of diff gene
hmExp=rt[rownames(diffSig),]
diffExp=rbind(id=colnames(hmExp),hmExp)
write.table(diffExp,file="diffExp.txt",sep="\t",quote=F,col.names=F)

#volcano
pdf(file="vol.pdf")
xMax=max(-log10(allDiff$P.Value))
yMax=max(abs(allDiff$logFC))
plot(-log10(allDiff$P.Value), allDiff$logFC, xlab="-log10(P.Value)",ylab="logFC",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.4)
diffSub=subset(allDiff, P.Value<adjustP & logFC>logFoldChange)
points(-log10(diffSub$P.Value), diffSub$logFC, pch=20, col="red",cex=0.4)
diffSub=subset(allDiff, P.Value<adjustP & logFC<(-logFoldChange))
points(-log10(diffSub$P.Value), diffSub$logFC, pch=20, col="green",cex=0.4)
abline(h=0,lty=2,lwd=3)
dev.off()
