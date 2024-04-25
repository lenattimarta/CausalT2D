### global causal model: causal discovery (n=100 repetitions of the learning process, each time on a random 90% portion of the training set)
library("readxl")
library("bnlearn")
library("Rgraphviz")
install.packages("caret")
library(caret)
#set working directory
setwd("...") 
# import auxiliary functions
source(".utils/add_arcs.R")
source(".utils/completeCPT.R")

#load data
data=read.table("./discretizedDataset.csv", sep=",", header=TRUE) 
data=data.frame(data)
summary(data) # training set size: 86618 rows, 20 columns {age,sex,Pressure,BMI,LDL,HDL,TG,FPG,Total_Cholesterol,Depression,HTN,OA,COPD,antidepressant,cholesterol_lowering_meds,antihypertensive_meds,corticosteroids,smoking_meds,smoking, T2D}

reduce_var=FALSE
if(reduce_var){
  reduced_var=c("age","sex", "BMI", "Pressure","FPG", "LDL", "HDL", "TG", "T2D", "HTN")
  data=data[,reduced_var]
}
n_bootstrap_rep=100
for(idx in 1:n_bootstrap_rep)
{ 
  #factorise data
  for(i in 1:(dim(data)[2])){
    data[,i]=factor(data[,i])
  }
  
  #blacklist definition based on prior knowledge and experts opinion (iterative process)
  name_vec=names(data)[-c(1:2)]
  medications=c("antidepressant", "cholesterol_lowering_meds", "antihypertensive_meds","corticosteroids","smoking_meds")
  comorb_and_smok=c("HTN","COPD","OA","T2D","Depression","smoking")
  biomarkers=c("LDL","HDL","TG","Total_Cholesterol","FPG","BMI","Pressure")
  bl_age=data.frame("from"=c(name_vec,names(data)[2]),"to"=rep(names(data)[1],length(name_vec)+1)) #BL of nodes entering to age
  bl_sex=data.frame("from"=c(name_vec,names(data)[1]),"to"=rep(names(data)[2],length(name_vec)+1)) #BL of nodes entering to sex
  bl_T2D=data.frame("from"="T2D","to"=names(data)[-c(1,2)]) #BL of nodes starting from T2D
  bl_medication=expand.grid("from"=medications,"to"=medications)
  bl_BMI=data.frame("from"=c("Pressure","LDL","HDL","TG","FPG"),"to"=c("BMI","BMI","BMI","BMI","BMI")) #extra BLs
  bl_extra=data.frame("from"=c("BMI","BMI"),"to"=c("LDL","HDL"))
  bl_extra2=expand.grid("from"=medications,"to"=comorb_and_smok)
  bl_extra3=expand.grid("from"=comorb_and_smok[1:(length(comorb_and_smok)-1)],"to"=comorb_and_smok[1:(length(comorb_and_smok)-1)])
  bl_extra4=expand.grid("from"=biomarkers,"to"="smoking")
  bl_extra5=data.frame("from"=c("HDL","HTN"),"to"=c("antihypertensive_meds","cholesterol_lowering_meds"))
  bl_extra6=data.frame("from"="HDL","to"="LDL")
  bl=rbind(bl_age,bl_sex, bl_T2D,bl_medication, bl_BMI,bl_extra,bl_extra2,bl_extra3,bl_extra4,bl_extra5,bl_extra6)
  
  #Extraction of a random portion of the dataset (90%), with the same frequency of T2D patients as in the whole dataset
  table(data$T2D) #initial frequency of subjects with future T2D onset 
  
  indices <- createDataPartition(data$T2D, p = 0.90, list = FALSE)
  data_train <- data[indices, ] # Sampled dataset stratified by output class
  table(data_train$T2D)#check frequency after stratification
  summary(data_train)
  
  #structure learning (Hill-Climbing) on the selected portion of the dataset
  bn_hc=hc(data_train,bl=bl)
  
  graph_file_name<- file.path("./Bootstrap_runs_global/graphs", paste0('graph', "_", idx, ".RData"))
  save(bn_hc, file = graph_file_name)
  graphviz.plot(bn_hc,  shape="ellipse",layout="dot", main="Global graph learnt on a random 90% portion of the training set") 
  graph_noLV_file_name <- file.path("./Bootstrap_runs_global/plots", paste0('graph', "_", idx, ".pdf"))
  #plot the global graph learnt on the selected 90\% of the training set
  pdf(file=graph_noLV_file_name, height=10, width=10)
  par(mfrow=c(1, 1))
  graphviz.plot(bn_hc,  shape="ellipse",layout="dot", main="Global graph learnt on a random 90% portion of the training set") 
  dev.off()
  #save arcs
  archs_noLV_file_name<- file.path("./Bootstrap_runs_global/tables", paste0('arcs', "_", idx, ".csv"))
  write.table(bn_hc$arcs,file=archs_noLV_file_name,sep=",",row.names=FALSE)
}