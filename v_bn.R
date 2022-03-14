
v_bn<-function(batch = batch){
  
  
  bn_per<-c()
  
  for(i in 1:nrow(batch)){
    
    group_id_aux <- rep(0,nrow(batch))
    group_id_aux[i]<-1
    
    bn_per[i]<-bn(group_id_aux, data = batch)
  }
  
  return(var(bn_per))
  
  
}