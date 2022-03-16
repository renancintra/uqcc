library('uclust')
library("future")

source("utils.R")
source("mpca_arma_fun.R")

future::plan(multisession)

# numero de instantes
n_col <-1000 # t = 20 e 2000
# numero de bateladas
n_orig <- c(500)
# numero bateladas fase 2
n_fase2 <-200
# acf
lag_max_acf<- n_col-1
#alfa_h1<-c(-0.5,	-0.2,	-0.1,	0,	0.1,	0.2,	0.5)
#alfa_h1<-c(-2,-1,-0.5	,-0.2	,-0.1,0,0.1,	0.2,	0.5,1,2  )
alfa_h1<-c(0.5)

#alfa_h1<-c(-0.1)
# descontrole -0.1
lista_das_listas_acf<-list()
lista_das_listas_drift<-list()
lista_das_listas_mpa_T<-list()
lista_das_listas_mpa_Q<-list()

inicio_0<-Sys.time()
for( n_ite_norig in 1:length(n_orig)){
  
  lista_alfa_acf <- list()
  lista_alfa_drift <- list()
  mpa_list_T2<-list()
  mpa_list_Q<-list()
  print("norig")
  print(n_orig[n_ite_norig])
  
  for( alfas_utili in 1:length(alfa_h1)){
    print("alfa")
    print(alfa_h1[alfas_utili])
    
    resultados_replicacao_acf <- c()
    resultados_replicacao_drift <- c()
    arl_Q<-c()
    arl_T2<-c()
    inicio<-Sys.time()
    #for(repl in 1:100){
    for(repl in 1:100){
      inicio_rep <- Sys.time()
      # base crua
      # fase 1
      
      mdserh0 <-
        compute_mdser(nrows=n_orig[n_ite_norig],
                      ncols=n_col,
                      ar_param_list=list(ar = 0.5))

      mdserh1 <-
        compute_mdser(nrows=n_fase2,
                      ncols=n_col,
                      ar_param_list = list(ma = alfa_h1[alfas_utili]))
      
      #------------------
      # ACF -------------
      #------------------
      
      mdh0_acf <-
        compute_acf(nrows=n_orig[n_ite_norig],
                    lag_max_acf = lag_max_acf,
                    mdser = mdserh0)

      mdh1_acf <-
        compute_acf(nrows=n_fase2,
                    lag_max_acf = lag_max_acf,
                    mdser = mdserh1)
      
      #------------------
      # Drift------------
      #------------------
      
      mean_serieh0 <-colMeans(mdserh0)
      
      
      mdserh0_drift<-mdserh0[]-mean_serieh0[col(mdserh0)]
      mdserh1_drift<-mdserh1[]-mean_serieh0[col(mdserh1)]
      
      # calculos bn --------------------
      
      inicio_acf<-Sys.time()
      
      # v_bn_calc_acf<-v_bn(mdh0_acf)
      
      v_bn_calc_acf <- compute_v_bn(mdh0_acf)
      
      bn_new_pad_acf <- compute_bn_pad(md_h0 = mdh0_acf, md_h1=mdh1_acf, pad=sqrt(v_bn_calc_acf))
      
      fim_acf<-Sys.time()
      if(repl < 5){
        print("Acf:")
        print(fim_acf-inicio_acf)
      }
      
      inicio_drift<-Sys.time()

      
      v_bn_calc_drift<-compute_v_bn(mdserh0_drift)
      
      bn_new_pad_drift <- compute_bn_pad(md_h0=mdserh0_drift, md_h1 = mdserh1_drift, pad=sqrt(v_bn_calc_drift))
      
      fim_drift<-Sys.time()
      if(repl < 5){
        print("Drift:")
        print(fim_drift-inicio_drift)
      }
      resultados_replicacao_acf[repl]<-ifelse(sum(bn_new_pad_acf > -qnorm( p = 0.05) ) == 0,NA, 1/(sum(bn_new_pad_acf > -qnorm( p = 0.05)  )/nrow(mdh1_acf))) 
      resultados_replicacao_drift[repl]<-ifelse(sum(bn_new_pad_drift > -qnorm( p = 0.05) ) ==0,NA, 1/(sum(bn_new_pad_drift > -qnorm( p = 0.05)  )/nrow(mdserh1_drift))) 
      
      inicio_arl<-Sys.time()
      arl_T2[repl]<-mpca_arma_fun(x = mdserh0,
                                  xnew = mdserh1,
                                  I = n_orig[n_ite_norig],
                                  Inew = nrow(mdserh1),
                                  n=ncol(mdserh1),
                                  confidence_alfa = 0.95)$arl_T2
      
      
      arl_Q[repl]<-mpca_arma_fun(x = mdserh0,
                                 xnew = mdserh1,
                                 I = n_orig[n_ite_norig],
                                 Inew = nrow(mdserh1),
                                 n=ncol(mdserh1),
                                 confidence_alfa = 0.95)$arl_Q
      fim_arl<-Sys.time()
      if(repl <5){
        print("Arl")
        print(fim_arl - inicio_arl)
      }
      
      fim_rep <- Sys.time()
      print("tempo da replicacao:")
      print(fim_rep-inicio_rep)
    }
    
    # print("armazenando os resultados do alfa:")
    # print(alfa_h1[alfas_utili])
    # print(n_orig[n_ite_norig])
    print(round(Sys.time()-inicio))
    
    lista_alfa_acf[[alfas_utili]]<-resultados_replicacao_acf
    lista_alfa_drift[[alfas_utili]]<-resultados_replicacao_drift
    
    
    
    mpa_list_T2[[alfas_utili]]<-arl_T2
    mpa_list_Q[[alfas_utili]]<-arl_Q
    
    
  }
  lista_das_listas_mpa_T[[n_ite_norig]]<-mpa_list_T2
  lista_das_listas_mpa_Q[[n_ite_norig]]<-mpa_list_Q
  lista_das_listas_acf[[n_ite_norig]]<-lista_alfa_acf
  lista_das_listas_drift[[n_ite_norig]]<-lista_alfa_drift
}
fim<-Sys.time()
fim-inicio_0

future::plan(future::sequential)
# acf ------------------------------------------

#df_resultados_acf
df_resultados_acf <- compute_results(lista_das_listas_acf, lista_alfa_acf,
                                     n_orig, alfa_h1)
readr::write_rds(df_resultados_acf, 'df_resultados_acf.rds')



#  drift ------------------------------------------

#df_resultados_drift
df_resultados_drift <- compute_results(lista_das_listas_drift, lista_alfa_drift,
                                       n_orig, alfa_h1)

readr::write_rds(df_resultados_drift, 'df_resultados_drift.rds')

#  alt T ------------------------------------------

#df_resultados_alT
df_resultados_alT <- compute_results(lista_das_listas_mpa_T, lista_alfa_drift,
                                     n_orig, alfa_h1)

readr::write_rds(df_resultados_alT, 'df_resultados_alT.rds')

#  alt Q ------------------------------------------

#df_resultados_alQ
df_resultados_alQ <- compute_results(lista_das_listas_mpa_Q, lista_alfa_drift,
                                     n_orig, alfa_h1)

readr::write_rds(df_resultados_alQ, 'df_resultados_alQ.rds')
