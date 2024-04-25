### global causal model: check stability among different portions of the training set (90% of the training data at each run)
library("readxl")
library("bnlearn")
library("Rgraphviz")
#set working directory
setwd("...") 
# import auxiliary functions
source(".utils/add_arcs.R")
source(".utils/completeCPT.R")

#load data
data=read.table(".\\discretizedDataset.csv", sep=",", header=TRUE) 
data=data.frame(data)
summary(data)

reduce_var=FALSE
if(reduce_var){
  reduced_var=c("age","sex", "BMI", "Pressure","FBS", "LDL", "HDL", "TG", "Diabetes", "HTN")
  data=data[,reduced_var]
}

#factorise data
for(i in 1:(dim(data)[2])){
  data[,i]=factor(data[,i])
}


#read bootstrapped networks file names as a list
filenames <- list.files("./Bootstrap_runs_global/graphs", pattern="*.RData", full.names=TRUE)
ldf <- lapply(filenames, load)
res <- lapply(ldf, summary)

list_of_networks <- list() #create empty list

for (i in filenames) { #loop through the files
  print(i)
  list_of_networks[[i]] <- get(load( i)) #add files to list position
}
#compute strength of each arc (i.e., frequency of appearance in the bootstrapped networks)
strength_after_bootstrapping=custom.strength(list_of_networks, names(data), weights = NULL, cpdag = TRUE, debug = FALSE)
# Compute the averaged global network structure by considering only arcs with frequency above a certain threshold
threshold=0.9 # strength threshold for inclusion in the averaged network structure.
final_BN=averaged.network(strength_after_bootstrapping, threshold)
inclusion.threshold(strength_after_bootstrapping)
#plot the average global network
graph_final_file_name<- file.path("./Bootstrap_runs_global", "averaged_global_graph_thr0.9_data_stratified90.pdf")
pdf(file=graph_final_file_name, height=10, width=10)
par(mfrow=c(1, 1))
graphviz.plot(final_BN,  shape="ellipse",layout="dot", main="averaged global network") 
dev.off()

#save archs and corresponding strengths
archs_strength_file_name<- file.path("./Bootstrap_runs_global", 'arcs_strength_bootstrap_100_stratified.csv')
write.table(strength_after_bootstrapping,file=archs_strength_file_name,sep=",",row.names=FALSE)

#compare the average global network after bootstrapping with the global network learnt on the whole training dataset
name_vec=names(data)[-c(1:2)]
medications=c("antidepressant", "cholesterol_lowering_meds", "antihypertensive_meds","corticosteroids","smoking_meds")
comorb_and_smok=c("HTN","COPD","OA","Diabetes","Depression","smoking")
biomarkers=c("LDL","HDL","TG","Total_Cholesterol","FBS","BMI","Pressure")
#blacklist definition based on prior knowledge and experts opinion (iterative process)
bl_age=data.frame("from"=c(name_vec,names(data)[2]),"to"=rep(names(data)[1],length(name_vec)+1)) #BL of nodes entering to age
bl_sex=data.frame("from"=c(name_vec,names(data)[1]),"to"=rep(names(data)[2],length(name_vec)+1)) #BL of nodes entering to sex
bl_diabetes=data.frame("from"="Diabetes","to"=names(data)[-c(1,2)]) #BL of nodes starting from diabetes
bl_medication=expand.grid("from"=medications,"to"=medications)
bl_BMI=data.frame("from"=c("Pressure","LDL","HDL","TG","FBS"),"to"=c("BMI","BMI","BMI","BMI","BMI")) #extra BLs
bl_extra=data.frame("from"=c("BMI","BMI"),"to"=c("LDL","HDL"))
bl_extra2=expand.grid("from"=medications,"to"=comorb_and_smok)
bl_extra3=expand.grid("from"=comorb_and_smok[1:(length(comorb_and_smok)-1)],"to"=comorb_and_smok[1:(length(comorb_and_smok)-1)])
bl_extra4=expand.grid("from"=biomarkers,"to"="smoking")
bl_extra5=data.frame("from"=c("HDL","HTN"),"to"=c("antihypertensive_meds","cholesterol_lowering_meds"))
bl_extra6=data.frame("from"="HDL","to"="LDL")
bl=rbind(bl_age,bl_sex, bl_diabetes,bl_medication, bl_BMI,bl_extra,bl_extra2,bl_extra3,bl_extra4,bl_extra5,bl_extra6)

#structure learning (Hill-Climbing)
bn_hc_ALL_DATA=hc(data,bl=bl)
graph_file_name<- file.path("./Bootstrap_runs_global", "global_graph_ALL.RData")
save(bn_hc_ALL_DATA, file = graph_file_name)
graphviz.plot(bn_hc_ALL_DATA,  shape="ellipse",layout="dot", main="global graph (whole training set)") 
graph_file_name <- file.path("./Bootstrap_runs_global", "global_graph_ALL.pdf")
#plot the structure learned on the whole dataset
pdf(file=graph_file_name, height=10, width=10)
par(mfrow=c(1, 1))
graphviz.plot(bn_hc_ALL_DATA,  shape="ellipse",layout="dot", main="global graph (whole training set)") 
dev.off()

#Compare the global network learned on the whole dataset with the averaged global network obtained after bootstrapping on 90% of the dataset
graphs_comparison=compare(bn_hc_ALL_DATA, final_BN, arcs=TRUE)