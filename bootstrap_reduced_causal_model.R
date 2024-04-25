### reduced causal model: causal discovery with latent variable 'Lifestyle' (n=100 repetitions of the learning process)

library("readxl")
library("bnlearn")
library("Rgraphviz")
#set working directory
setwd("...") 
# import auxiliary functions
source(".utils/add_arcs.R")
source(".utils/completeCPT.R")
set.seed(27272)
n_bootstrap_rep=100
#load data
data=read.table("./discretizedDataset.csv", sep=",", header=TRUE) #load data
data=data.frame(data)# training set size: 86618 rows, 20 columns {age,sex,Pressure,BMI,LDL,HDL,TG,FPG,Total_Cholesterol,Depression,HTN,OA,COPD,antidepressant,cholesterol_lowering_meds,antihypertensive_meds,corticosteroids,smoking_meds,smoking, T2D}
summary(data)
#select subset of features to be considered in the reduced model
reduce_var=TRUE
if(reduce_var){
  reduced_var=c("age","sex", "BMI", "Pressure","FPG", "LDL", "HDL", "TG", "T2D")
  data=data[,reduced_var]
}
for(idx in 1:n_bootstrap_rep)
{ 
  print(idx)
  #add latent variable
  num_LV=1
  data[ "Lifestyle"] = factor(rep(NA, nrow(data)), levels=c("1", "2", "3", "4"))
  #factorise data
  for(i in 1:(dim(data)[2])){
    data[,i]=factor(data[,i])
}
  data$Lifestyle=factor(data$Lifestyle,levels=c("1", "2", "3", "4"))
  
  summary(data)

  #creation of the initial blacklist based on prior knowledge and experts opinion (iterative process)
  which_LV=(dim(data)[2]-(num_LV-1)):(dim(data)[2])
  name_vec=names(data)[-c(1:2,which_LV)]
  bl_age=data.frame("from"=c(name_vec,names(data)[2]),"to"=rep(names(data)[1],length(name_vec)+1)) #BL of nodes entering to age (root node)
  bl_sex=data.frame("from"=c(name_vec,names(data)[1]),"to"=rep(names(data)[2],length(name_vec)+1)) #BL of nodes entering to sex (root node)
  bl_T2D=data.frame("from"="T2D","to"=names(data)[-c(1,2,which_LV)]) #BL of nodes starting from T2D
  bl_BMI=data.frame("from"=c("Pressure","LDL","HDL","TG","FPG"),"to"=c("BMI","BMI","BMI","BMI","BMI")) #extra BLs
  bl_extra=data.frame("from"=c("BMI","BMI"),"to"=c("LDL","HDL"))
  
  bl_start_noLV=rbind(bl_age,bl_sex, bl_T2D,bl_BMI,bl_extra)
  bl_file_name <- file.path("./Bootstrap_runs/tables", paste0('black_list_bootstrap', "_", idx, ".csv"))
  write.table(bl_start_noLV,file=bl_file_name, sep=",",row.names=FALSE)
  
  #creation of the whitelist for LV based on prior knowledge and experts opinion (iterative process), one per LV
  wl_Lifestyle=data.frame("from"=rep("Lifestyle",7),"to"=c("BMI","Pressure","FPG","TG","LDL","HDL","T2D")) 
  wl_em=rbind(wl_Lifestyle)# row binding of whitelists(if more than one LV is present)
  
  #creation of a blacklist including archs pointing from latent variables to non-modifiable (root) nodes 
  bl_ageLV=expand.grid(from=names(data[,which_LV]),to="age")
  bl_sexLV=expand.grid(from=names(data[,which_LV]),to="sex")
  #creation of a blacklist including archs pointing to LVs (LVs must be modeled as root nodes)
  bl_Lifestyle_in=expand.grid(from=names(data[,-which(names(data)=="Lifestyle")]), to="Lifestyle") 
  bl_Lifestyle=expand.grid(from="Lifestyle",to=setdiff(c("age","sex", name_vec),wl_Lifestyle$to)) #block all archs that are not in the white list
  #set blacklist (no LV)
  bl=rbind(bl_start_noLV)
  #set blacklist for structural EM
  bl_em=rbind(bl,bl_ageLV,bl_sexLV,bl_Lifestyle, bl_Lifestyle_in)
  data_train=data
  
  #1) train starting graph (Hill-Climbing) excluding LVs
  bn_hc=hc(data_train[,-which_LV],bl=bl)
  graph_noLV_file_name<- file.path("./Bootstrap_runs/graphs", paste0('graph_noLV', "_", idx, ".RData"))
  save(bn_hc, file = graph_noLV_file_name)
  graphviz.plot(bn_hc,  shape="ellipse",layout="dot", main="graph without LV") 
  graph_noLV_file_name <- file.path("./Bootstrap_runs/plots", paste0('graph_noLV', "_", idx, ".pdf"))
  
  pdf(file=graph_noLV_file_name, height=10, width=10)
  par(mfrow=c(1, 1))
  graphviz.plot(bn_hc,  shape="ellipse",layout="dot", main="graph without LV") #no latent variables
  dev.off()
  
  #2) EM to learn the dag in presence of LVs with whitelist for LVs
  smallsample=FALSE   #if true, structure learning is performed on a reduced data sample
  nmax=dim(data_train)[1]
  if(smallsample){
    n=20000
  }else{
    n=nmax
  }
  
  #learn the graph with LVs by means of the Structural EM algorithm
  bn_hc=add_to_dag_wl(bn_hc,wl_Lifestyle) 
  start0=bn.fit(bn_hc, data_train,method="bayes") #use the dag without LV as a starting point for structural EM
  struct0=structural.em(data_train[sample(1:nmax,n),],maximize="hc",maximize.args=list(whitelist=wl_em,blacklist=bl_em), start=start0, max.iter=3, fit="bayes", return.all=TRUE)#, impute="exact")
  bn_hc_em_wl=struct0$dag
  fit_em0=struct0$fitted
  
  #plot the reduced graph with and without LV
  graph_final_file_name<- file.path("./Bootstrap_runs/plots", paste0('graph_LV', "_", idx, ".pdf"))
  pdf(file=graph_final_file_name, height=10, width=20)
  par(mfrow=c(1, 2))
  graphviz.plot(bn_hc,  shape="ellipse",layout="dot", main="graph without LV") 
  graphviz.plot(bn_hc_em_wl, shape="ellipse",layout="dot", main="graph with LV") 
  dev.off()
  #save arcs
  archs_noLV_file_name<- file.path("./Bootstrap_runs/tables", paste0('graph without LV', "_", idx, ".csv"))
  write.table(bn_hc$arcs,file=archs_noLV_file_name,sep=",",row.names=FALSE)
  archs_LV_WL_file_name<- file.path("./Bootstrap_runs/tables", paste0('graph with LV', "_", idx, ".csv"))
  write.table(bn_hc_em_wl$arcs,file=archs_LV_WL_file_name,sep=",",row.names=FALSE)
}