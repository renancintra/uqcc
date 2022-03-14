
mpca_arma_fun=function(x,xnew,I,Inew,n,
                       confidence_alfa = 0.95){
  # n=10
  # I=100
  # Inew=500
  # x=matrix(ncol=n,nrow=I)
  # for(i in 1:I){
  #   x[i,] <- arima.sim(model=list(ar=0.2,ma=0.5),n=n)
  # }
  
  ##################PCA pra salvar os autovalores e autovetores
  
  av=eigen(cov(x))
  
  pos_av=av$values[which(av$values>0)]
  
  
  sav=c()
  sav[1]=0
  i=1
  while(sav[i]<0.90){
    sav[i+1]=sum(pos_av[1:i])/sum(pos_av)
    i=i+1
  }
  i=(i-1)
  
  
  va_pc=av$values[1:i]
  ve_pc=av$vectors[,1:i]
  
  
  ############Salvar xpc_mean e  xpc_cov pra estatostca de Hotelling
  
  xpc=x%*%ve_pc
  
  xpc_mean=colMeans(xpc)
  
  xpc_cov=cov(xpc)
  
  ##############################################################
  
  #T2 chart
  
  #########################################Limite usando a dist F
  confidence.level=confidence_alfa
  
  #T2_lim=((i*(I-1))/(I-i))*qf(confidence.level,i,(I-i))
  T2_lim = qchisq( df = i, p = confidence.level) 
  
  ################################Limite Empírico
  # T2=c()
  # for (k in 1:I){
  #   T2[k]= mahalanobis(
  #     xpc[k,],
  #     center = xpc_mean,
  #     cov = xpc_cov,
  #     inverted = FALSE)
  # }
  
  #T2_lim_emp=quantile(T2, confidence.level)
  ##############################################################
  
  #Q chart
  
  ########################################Limite usando a dist Normal 
  TT=length(pos_av)
  
  t1=sum(pos_av[(i+1):TT])
  t2=sum(pos_av[(i+1):TT]^2)
  t3=sum(pos_av[(i+1):TT]^3)
  ho=1-((2*t1*t3)/(3*t2^2))
  z=qnorm(confidence.level+(1-confidence.level)/2)
  
  a=z*sqrt(2*t2*ho^2)/t1
  
  b=t2*ho*(ho-1)/t1^2
  
  
  Q_lim=t1*(a+1+b)^(1/ho)
  
  
  #########################Limite usando a dist empírica
  # 
  # x_est=xpc%*%t(ve_pc)
  # 
  # QQo=x_est-x
  # 
  # QQ=c()
  # for (k in 1:I){
  #   QQ[k]= sum((QQo[k,])^2)
  # }
  # 
  # Q_lim_emp=quantile(QQ, confidence.level)
  
  ################################################
  #############################################################
  
  ### Fase 2
  
  # xnew=matrix(ncol=n,nrow=Inew)
  # for(j in 1:Inew){
  #   xnew[j,] <-  arima.sim(model=list(ar=0.2,ma=0.5),n=n)
  # }
  
  #########T2 new
  xpc_new=xnew%*%ve_pc
  
  T2_new=c()
  for (k in 1:Inew){
    T2_new[k]= mahalanobis(
      xpc_new[k,],
      center = xpc_mean,
      cov = xpc_cov,
      inverted = FALSE)
  }
  #
  
  #########Q new
  
  ############Calculando residuo usando somatorio de z^2 de i+1 até TT
  
  #TT=length(pos_av)
  
  ve_pc_q=av$vectors[,(i+1):TT]
  
  
  xpc_new_q=xnew%*%ve_pc_q
  
  
  Q_new=c()
  for (k in 1:Inew){
    Q_new[k]= sum((xpc_new_q[k,])^2)
  }
  
  #########################% fora e ARL
  ##T2 ARL
  
  #sum(T2_new>T2_lim)/length(T2_new); 
  
  
  arl_T2=ifelse(sum(T2_new>T2_lim)==0,NA,1/(sum(T2_new>T2_lim)/length(T2_new))) 
  
  ##T2
  
  #sum(T2_new>T2_lim_emp)/length(T2_new); 1/(sum(T2_new>T2_lim_emp)/length(T2_new))
  #arl_T2_emp=ifelse(sum(T2_new>T2_lim_emp)==0,NA,1/(sum(T2_new>T2_lim_emp)/length(T2_new))) 
  
  ## Q ARL
  #sum(Q_new>Q_lim)/length(Q_new); 1/(sum(Q_new>Q_lim)/length(Q_new))###Q
  
  arl_Q=ifelse(sum(Q_new>Q_lim)==0,NA, 1/(sum(Q_new>Q_lim)/length(Q_new))) 
  
  #sum(Q_new>Q_lim_emp)/length(Q_new); 1/(sum(Q_new>Q_lim_emp)/length(Q_new))###Q
  #arl_Q_emp=ifelse(sum(Q_new>Q_lim_emp)==0,NA, 1/(sum(Q_new>Q_lim_emp)/length(Q_new))) 
  
  return(list(arl_T2=arl_T2,
              #arl_T2_emp=arl_T2_emp,
              arl_Q=arl_Q#,
              #arl_Q_emp=arl_Q_emp
              ))
  
}

