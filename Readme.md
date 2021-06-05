## Project Name
Identification of key mRNAs as prediction models for early metastasis of pancreatic cancer based on LASSO.
## Summary
We summarize the method of this study, including method code and input and output files.
## Environmental dependence
R version 3.6.2 
## Method Code
### Difference Analysis
In order to identify EMT differential genes, we screened 994 EMT-related genes from the GSEA gene set. We then want to explore the relationship between EMT-related genes and pancreatic cancer metastasis. We included the public database TCGA human pancreatic cancer transcriptome sequencing data, our own cell line BXPC-3 and BXPC-3-M8 transcriptome sequencing data, and two sets of microarray data from the public database GEO (GSE23952, GSE21654). We matched the expression data of 994 EMT-related genes in human samples and cell lines. Next, we analyzed the difference between the free group and the metastatic group of each set of data. The "[Limma](/limma/)" package of [R Software 3.6.2](https://www.r-project.org/)  was used to process chip data, and the "[edgeR](/edgeR/)" package of R software was used to screen the mRNAs differentially expressed between different groups using the obtained sequencing data. The criteria ![1](https://latex.codecogs.com/svg.image?%5Cleft%7C%20log%5C%20Fold%5C%20Change%5Cright%7C%3E1) and ![2](https://latex.codecogs.com/svg.image?P%5C%20value%3C0.05)  were set as the threshold criteria. We extracted the overlapping differentially expressed genes (DEGs) from the GEO, BXPC-3, BXPC-3-M8, and TCGA datasets for the following analyses.
### [LASSO Regression Model Construction](/2.lasso.R)
To construct a risk score model for pancreatic cancer metastasis prediction, we developed risk scores using the LASSO regression algorithm. The “glmnet” package in R software was used to build the lasso feature screening model. Screened genes and their coefficients were defined by minimum binominal deviance. The formula for the risk score was as follows:

![3](https://latex.codecogs.com/svg.image?Risk%5C%20Score=%5Csum_%7Bi=1%7D%5E%7Bn%7DCoef_%7Bi%7D*x_%7Bi%7D)

where ![4](https://latex.codecogs.com/svg.image?Coef_%7Bi%7D) is coefficient, and ![5](https://latex.codecogs.com/svg.image?x_%7Bi%7D)  is the expression level of the corresponding gene in the sample.
### [ROC curve](/3..ROC.R)
The area under the ROC curve (AUC) is an accurate indicator in diagnostic tests, which is used to evaluate the quality of the model. The "ROCR" package in R software was used to generate the ROC curve.
### [Multi-group Heat Map Analysis](/4.multi-group%20heat%20map%20analysis.R)
Through the analysis of multi-group heat map, we explored the relationship between metastasis status and patient characteristics (tumor grade, TNM stage, and survival status, age, gender). The "pheatmap" package in R software was used to  generate heat maps.
### Univariate and Multivariate Cox Analyses
To further verify the effect of risk score and patient characteristic features on the prognosis of pancreatic cancer, [univariate](5.Univariate%20Cox.R) and [multivariate](5.Multivariate%20cox.R) Cox analyses were performed for risk score, age, sex, histological grade, and TNM stage. The “survival” package in R software was used to construct the Cox model.
## Contact
If there is any problem with the code, please contact us: 221909252059@zust.edu.cn

