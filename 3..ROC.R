#1.常见ROC曲线
#--------------------------------------ROCR包---------------------------------------------------
library(ROCR)
mydata=read.table("data.txt",sep="\t",header=T)
View(mydata)
#predictions为风险评分，labels为状态
pred<-prediction(mydata$predictions,mydata$labels) 
perf<-performance(pred,measure = "tpr",x.measure = "fpr")
plot(perf,col=rainbow(10))
auc <- performance(pred,'auc')
auc=unlist(slot(auc,"y.values"))
rm(list = ls())
options(stringsAsFactors = F)
#2.时间依赖的ROC曲线
#载入数据
KM.input<-read.csv(file = "KM.input.csv",header = T)
head(KM.input)
# X event time_year RiskScore
# 1 TCGA-A1-A0SE-01A-11R-A085-13     0 3.6191781 0.3686996
# 2 TCGA-A1-A0SH-01A-11R-A085-13     0 3.9369863 1.8872376
# 3 TCGA-A1-A0SJ-01A-11R-A085-13     0 1.1397260 0.9330429
# 4 TCGA-A1-A0SK-01A-12R-A085-13     1 2.6493151 0.5337272
# 5 TCGA-A1-A0SM-01A-11R-A085-13     0 0.6630137 3.0229704
# 6 TCGA-A1-A0SO-01A-22R-A085-13     0 2.3342466 0.8083318
#数据包括样本名称、事件(生：0,死：1),生存时间，风险值。该数据源自cox风险比例模型

#--------------------------------------survivalROC包---------------------------------------------------
#1.载入R包
#install.packages("survivalROC")
#install.packages("timeROC") #2个包都可以绘制生存时间依赖的ROC曲线
#1.使用SurvivalROC （不能显示置信区间和SD）
library(survival)
library(survivalROC)
?survivalROC #查看这个函数的格式
#输入数据
Survival_ROC_input<-KM.input
#这里我预测五年的生存率
survival_ROC<-survivalROC(Stime=Survival_ROC_input$time_year, #生存时间，Event time or censoring time for subjects
                          status=Survival_ROC_input$event, #生存状态,dead or alive
                          marker=Survival_ROC_input$RiskScore, #风险得分，Predictor or marker value
                          predict.time=5, #预测5年的生存时间
                          method="KM" #使用KM法进行拟合，默认的方法是method="NNE"
)
survival_ROC_AUC<-round(survival_ROC$AUC,3)#ROC曲线的AUC保留3位小数（文章保留了3位）
#画图
#x轴为False positive rate，y轴为True positive rate
plot(survival_ROC$FP,survival_ROC$TP,type="l",xlim=c(0,1),ylim=c(0,1),
     xlab="False positive rate",  
     ylab="True positive rate",
     main=paste0("ROC Curve", " (", "AUC = ",survival_ROC_AUC," )"),  #标题
     cex.main=1.5,#放大标题
     cex.lab=1.3,#坐标轴标签（名称）的缩放倍数
     cex.axis=1.3, font=1.2, #放大轴上的数字
     lwd=1.5, #指定线条宽度
     col="red"  #红色
)
abline(a=0,b=1) #添加一条y=x的线
#计算最佳cutoff
cutoff<-survival_ROC$cut.values[which.max(survival_ROC$TP-survival_ROC$FP)]
cutoff

#---------------------------------------timeROC包------------------------------------------------------------------------
#2.使用timeROC(可以计算置信区间和SD）
#Time ROC可以时计算多个时间的AUC
#文章没有说明计算的几年的，我这里计算3，5，10年
library(timeROC)
library(survival)
?timeROC #看一下说明书
#输入数据
time_ROC_input<-KM.input
time_ROC<-timeROC(T=time_ROC_input$time_year, #生存时间(dead和alive的生存时间).
                  delta=time_ROC_input$event, #生存结局，Censored的样本必须用0表示
                  marker=time_ROC_input$RiskScore, #预测的变量，这里是风险评分，在没有特殊情况下，值越大，越容易发生危险事件
                  cause=1, #阳性结局的赋值（必须是1或2），也就是dead的赋值，这里dead是1表示的
                  weighting = "marginal", #权重计算方法，这是默认方法，选择KM计算删失分布，weighting="aalen" [选用COX]，weighting="aalen" [选用Aalen]
                  times = c(3,5,10), #计算3、5、10年的ROC曲线
                  ROC=TRUE,
                  iid=TRUE #计算AUC
)
time_ROC #查看结果，可以看到，还包括了SE
#绘制ROC曲线啦
summary(time_ROC) #返回12个参数
time_ROC$AUC
#绘制3年的ROC
plot(time_ROC,time=3,col="red",title=FALSE,lwd=2) #将生成一条两倍于 默认宽度的线条
#在此基础上添加5年的ROC
plot(time_ROC,time=5,col="blue",add=TRUE,title=FALSE,lwd=2)
#add 10年的ROC
plot(time_ROC,time=10,col="black",add=TRUE,title=FALSE,lwd=2)
#添加图例
?legend
legend("bottomright", #图例设置在右下角
       c(paste0("AUC at 3 years = ", round(time_ROC$AUC[1],3)),
         paste0("AUC at 5 years = ", round(time_ROC$AUC[2],3)),
         paste0("AUC at 10 years = ", round(time_ROC$AUC[3],3))),
       col=c("red","blue","black"),lwd=2,bty="n")