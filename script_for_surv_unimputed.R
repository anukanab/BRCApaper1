#   AnuSA.R
#  /DS/Nathan/AnuSurvival
#  January 11, 2021
in1DataPath="/Volumes/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana-New-Analysis/TF/SurvivalAnalysis/unimputed/"
outDataPath=in1DataPath
setwd(in1DataPath)

library(survival)
#setwd("/Volumes/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana-New-Analysis/TF/SurvivalAnalysis/unimputed/")
#pull merged file with clinical data and cluster IDs
f=paste0(in1DataPath,"eventannotation_TF_T.txt")
#BRCA_data <- read.table("/Volumes/salomonis2/NCI-R01/TCGA-BREAST-CANCER/Anukana-New-Analysis/TF/SurvivalAnalysis/unimputed/eventannotation_TF_T.txt", sep="\t",header = TRUE)
BRCA_data=read.table(file=f,sep="\t",
            header=T,stringsAsFactors = F)
length(BRCA_data)
#remove - and replace with _
x=colnames(BRCA_data)
colnames(BRCA_data)=gsub("-","_",x)
x=colnames(BRCA_data)
sample(colnames(BRCA_data)[9983:9999]) # to check if the colnames written below matches with how R reads it
colnames(BRCA_data)[9982]
#Survival (Cancer and normal samples)
meta=Surv(BRCA_data$days_to_death.Clinical.Significance.BREAST.CANCER,
          BRCA_data$vital_status.Clinical.Significance.BREAST.CANCER)
#  for assessing viability of SA
meta.df=data.frame(
   daysTo=BRCA_data$days_to_death.Clinical.Significance.BREAST.CANCER,
   vitalStat=BRCA_data$vital_status.Clinical.Significance.BREAST.CANCER)

#create data frame to save results
results=as.data.frame(matrix(ncol=5, nrow=9983))
colnames(results)=c("column", "LRT", "Wald", "SlrT","Zscore")
results[1:9982,1]=colnames(BRCA_data)[2:9983]

#iterate through columns of interest and do coxph, save pvalues for 3 tests from summary in results
coxph.minN=3
for(i in 2:9983){
# i=2
#  values of predictor var
  xv=BRCA_data[,x[i]]
#  join with survival info
  assess.df=cbind(meta.df,xv)
#  remove rows w/ NA value
  assess.df=assess.df[!is.na(assess.df$xv),]
#  how many vitalStat==1 ?
   n1s=sum(assess.df$vitalStat) 
#  only run coxph if enough vitalStat==1
  if (n1s >= coxph.minN) {
  form=as.formula(paste0("meta~", paste0("BRCA_data$", x[i])))
  res.cox=coxph(form)
  summ=summary(res.cox)
#  print(summ)
  results[(i-1), 2]=summ$logtest[3]
  results[(i-1), 3]=summ$waldtest[3]
  results[(i-1), 4]=summ$sctest[3]
  results[(i-1), 5]=summ$coefficients[1,4] }
# fill results row w/ -1s if was not run
  if (n1s < coxph.minN) { results[(i-1),2:5]=rep(-1,4)}
}

sum(results$LRT==-1,na.rm = T) # 184

#write results
#write.table(results, "results_combined.txt", sep="\t", row.names=F)

meta1=Surv(BRCA_data$years_to_birth.Clinical.Significance.BREAST.CANCER)

#iterate through columns of interest and do coxph, save pvalues for 3 tests from summary in results
coxph.minN=3
for(i in 2:9983){
  # i=2
  #  values of predictor var
  xv=BRCA_data[,x[i]]
  #  join with survival info
  assess.df=cbind(meta.df,xv)
  #  remove rows w/ NA value
  assess.df=assess.df[!is.na(assess.df$xv),]
  #  how many vitalStat==1 ?
  n1s=sum(assess.df$vitalStat) 
  #  only run coxph if enough vitalStat==1
  if (n1s >= coxph.minN) {
    form=as.formula(paste0("meta1~", paste0("BRCA_data$", x[i])))
    res.cox=coxph(form)
    summ=summary(res.cox)
    #  print(summ)
    results[(i-1), 2]=summ$logtest[3]
    results[(i-1), 3]=summ$waldtest[3]
    results[(i-1), 4]=summ$sctest[3]
    results[(i-1), 5]=summ$coefficients[1,4] }
  # fill results row w/ -1s if was not run
  if (n1s < coxph.minN) { results[(i-1),2:5]=rep(-1,4)}
}

sum(results$LRT==-1,na.rm = T) # 184

#write results
write.table(results, "results_combined_ageC.txt", sep="\t", row.names=F)

