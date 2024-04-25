
add_to_dag_wl=function(dag, wl){
  LV=wl$from[1]
  for(node in wl$to){
    dag=set_dag_LV(dag, LV, node)
  }
  return(dag)
}


set_dag_LV=function(dag, LV, nodes){
  dag$nodes[[LV]]$mb=c(dag$nodes[[LV]]$mb,nodes)
  dag$nodes[[LV]]$nbr=c(dag$nodes[[LV]]$nbr,nodes)
  dag$nodes[[LV]]$parents=c(dag$nodes[[LV]]$parents,character(0))
  dag$nodes[[LV]]$children=c(dag$nodes[[LV]]$children,nodes)
  for(node in nodes){
    dag$nodes[[node]]$mb=c(dag$nodes[[node]]$mb, LV)
    dag$nodes[[node]]$nbr=c(dag$nodes[[node]]$nbr, LV)
    dag$nodes[[node]]$parents=c(dag$nodes[[node]]$parents, LV)
  }
  return(dag)
}
