### averaged reduced causal model: averaged network obtained by considering only the most frequent arcs (n=100 repetitions of the learning process)
library("readxl")
library("bnlearn")
library("Rgraphviz")

#set working directory
setwd("...") 
# import auxiliary functions
source("./utils/add_arcs.R")
source("./utils/completeCPT.R")

#load data
data=read.table("./discretizedDataset.csv", sep=",", header=TRUE) 
data=data.frame(data)
summary(data)

#select subset of features to be considered in the reduced model
reduce_var=TRUE
if(reduce_var){
  reduced_var=c("age","sex", "BMI", "Pressure","FPG", "LDL", "HDL", "TG", "T2D")
  data=data[,reduced_var]
}

#add latent variables
num_LV=1
data[ "Lifestyle"] = factor(rep(NA, nrow(data)), levels = c("1", "2","3","4")) #the number of levels for each latent variable needs to be specified
#factorise data
for(i in 1:(dim(data)[2])){
  data[,i]=factor(data[,i])
}
data$Lifestyle=factor(data$Lifestyle,levels=c("1","2","3","4"))
summary(data)

#read bootstrapped networks file names as a list
filenames <- list.files("./Bootstrap_runs/graphs", pattern="*.RData", full.names=TRUE)
ldf <- lapply(filenames, load)
res <- lapply(ldf, summary)
names(res) <- substr(filenames, 6, 30)
list_of_networks <- list() #create empty list

for (i in filenames) { #loop through the files
  print(i)
  list_of_networks[[i]] <- get(load( i)) 
}
#compute strength of each arc (i.e., frequency of appearance in the bootstrapped networks)
strength_after_bootstrapping=custom.strength(list_of_networks, names(data), weights = NULL, cpdag = TRUE, debug = FALSE)
# Compute the averaged network structure by considering only archs with frequency above a certain threshold
threshold=0.9 # strength threshold for inclusion in the averaged network structure.
final_BN=averaged.network(strength_after_bootstrapping, threshold)
inclusion.threshold(strength_after_bootstrapping)

#plot the averaged reduced causal model
graph_final_file_name<- file.path("./Bootstrap_runs", "averaged_reduced_graph_thr0.9.pdf")
pdf(file=graph_final_file_name, height=10, width=10)
par(mfrow=c(1, 1))
graphviz.plot(final_BN,  shape="ellipse",layout="dot", main="averaged reduced graph with LV") 
dev.off()

#save arcs and corresponding strengths
archs_strength_file_name<- file.path("./Bootstrap_runs", 'arcs_strength_bootstrap_100.csv')
write.table(strength_after_bootstrapping,file=archs_strength_file_name,sep=",",row.names=FALSE)
