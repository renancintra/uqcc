library('furrr')
library('future')

ma <- function(arr, n=10){
  res = arr
  for(i in (floor(n/2)+1):(length(arr)-floor(n/2))){
    res[i] = mean(arr[(i-floor(n/2)):(i+floor(n/2))])
  }
  res
}

# Dado o numero de linhas e colunas de uma matriz, gerar uma simulação
# do arima com os parametros passados como lista em `ar_param_list`
compute_mdser <-
    function(nrows, ncols, ar_param_list) {
        matrix_na <-
            matrix(rep(NA), 
                      nrow = nrows,
                      ncol = ncols)
        for(i in 1:nrows) {
            matrix_na[i, ] <- arima.sim(n=ncols, model=ar_param_list)
        }
        return(matrix_na)
    }

compute_acf <-
    function(nrows, lag_max_acf, mdser) {

        matrix_na <-
            matrix(rep(NA), 
                   nrow = nrows,
                   ncol = lag_max_acf+1)
        for(i in 1:nrows) {
            serie <- mdser[i, ]
            serie_ma <- serie-ma(serie)
            aux <- NULL
            aux <- acf(serie_ma, plot=FALSE, lag.max=lag_max_acf)$acf
            matrix_na[i, ] <- aux
        }
        return(matrix_na)
    }

compute_v_bn <-
    function(matrix_batch) {
  
        bn_per <- furrr::future_map(
            1:nrow(matrix_batch),
            function(idx) {
                group_id_aux <- rep(0, nrow(matrix_batch))
                group_id_aux[idx] <- 1

                bn_per_idx <- bn(group_id_aux, data=matrix_batch)

                return(bn_per_idx)
            },
            .options = furrr::furrr_options(seed=NULL)
        )
        return(var(unlist(bn_per), na.rm = T))
}

compute_bn_pad <-
    function(md_h0, md_h1, pad) {

        bn_new_pad <-
            furrr::future_map(
                1:nrow(md_h1),
                function(idx) {
                    bn_new <-
                        uclust::bn(group_id = c(rep(0,nrow(md_h0)),1),
                                   data=rbind(md_h0, md_h1[idx, ]))
                    return(bn_new/pad)
                }
            )
        return(unlist(bn_new_pad))
    }

compute_results <-
  function(lista_das_listas, lista_alfa, n_orig, alfa_h1) {
    
    df_results <- data.frame()
    
    for(j in 1:length(lista_das_listas)) {
      for(i in 1:length(lista_alfa)) {
        
        aux <- NULL
        media <- mean(lista_das_listas[[j]][[i]],na.rm =TRUE)
        desvpad <- sd(lista_das_listas[[j]][[i]],na.rm =TRUE)
        median <- median(lista_das_listas[[j]][[i]],na.rm =TRUE)
        n_origem <- n_orig[[j]]
        coeficiente <- alfa_h1[[i]]
        
        aux <- data.frame(media,
                          desvpad,
                          median,
                          n_origem,
                          coeficiente)
        df_results <- rbind(df_results, aux)
        
      }
    }
    
    return(df_results)
    
  }
