completeCPT=function(data, dag, node){
  dim_vec=length(levels(data[,node]))
  for(i in 1:length(dag$nodes[[node]]$parents)){
    dim_vec=c(dim_vec,length(levels(data[,dag$nodes[[node]]$parents[i]])))
  }
  v=rep(1/dim_vec[1],prod(dim_vec))
  CPT=array(v,dim=dim_vec)
  dim(CPT)=dim(CPT)
  dimnames(CPT)[[1]]=levels(data[,node])
  for(i in 1:length(dag$nodes[[node]]$parents)){
    dimnames(CPT)[[i+1]]=levels(data[,dag$nodes[[node]]$parents[i]])
  }
  names(dimnames(CPT))=c(node,dag$nodes[[node]]$parents)
  class(CPT)="table"
  class = ifelse(is(data[, node], "ordered"), "bn.fit.onode", "bn.fit.dnode")
  return(structure(list(node = node, parents = dag$nodes[[node]]$parents, children = dag$nodes[[node]]$children,
                                    prob = CPT), class = class))
}


completeCPT_LV=function(data, dag, node){
  dim_vec=length(levels(data[,node]))
  
  tab=matrix(runif(dim_vec),1,dim_vec)
  tab=tab/sum(tab)
  CPT=array(tab,dimnames=list(colnames(tab)))
  dim(CPT)=dim(CPT)
  dimnames(CPT)[[1]]=levels(data[,node])
  names(dimnames(CPT))=node
  class(CPT)="table"
  class = ifelse(is(data[, node], "ordered"), "bn.fit.onode", "bn.fit.dnode")
  return(structure(list(node = node, parents = NULL, children = dag$nodes[[node]]$children,
                        prob = CPT), class = class))
}

set_dag_LV=function(dag, LV, nodes){
  dag$nodes[[LV]]$mb=nodes
  dag$nodes[[LV]]$nbr=nodes
  dag$nodes[[LV]]$parents=character(0)
  dag$nodes[[LV]]$children=nodes
  for(node in nodes){
    dag$nodes[[node]]$mb=c(dag$nodes[[node]]$mb, LV)
    dag$nodes[[node]]$nbr=c(dag$nodes[[node]]$nbr, LV)
    dag$nodes[[node]]$parents=c(dag$nodes[[node]]$parents, LV)
  }
  return(dag)
}
