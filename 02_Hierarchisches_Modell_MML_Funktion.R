### Hierarchical Model of Marginal Maximal Likelihood
rm(list = ls())

#Funktionen
load("functions/functions.RData")


##### mml.hm-Funktion #####
mml.hm <- function(data_IR,
                   data_RT,
                   Anz_QPunkte = 20,
                   start_c = 0,
                   start_v = 1,
                   prune = NULL,
                   EM_conv_criterion = 0.001,
                   EM_iter_criterion = 500,
                   newtonraphson_conv_criterion = 0.001,
                   newtonraphson_iter_criterion = 1000,
                   D=1){

  
  ######## Vorbereitungen ########
  ##### Funktionen
  #PQ-Formel
  PQ_fun <- function(y, beta, alpha, theta, D) {
    z <- D*(alpha * theta + beta)
    P_theta <- 1 / (1 + exp(-z))
    PQ_fun <- P_theta ^ y * (1 - P_theta) ^ (1 - y)
  }
  
  #Log-Normalverteilung
  rtime <- function(t, ny, tau, omega) {
    if(is.na(t)){
      rtime <- NA
    } else if(t>0){
      z <- (ny - tau)
      rtime <- omega / (t* sqrt(2 * pi)) * 
        exp(-((omega) ^ 2 / 2) * (log(t) -z) ^ 2)
    } else if(any(t<=0)){
      stop('response time not greater than 0')
    }
  }
  
  ##### Gaussian quadrature
  SIGMA_G <- matrix(c(1, start_c, start_c, start_v), 2, 2)
  X <-
    mgauss.hermite(Anz_QPunkte,
                   mu = c(0, 0),
                   sigma = SIGMA_G,
                   prune = prune)
  
  ##### Anzahl Parameter
  I <- ncol(data_IR)
  P <- nrow(data_IR)
  Q <- nrow(X$points)
  
  
  ##### Data Spalten umbenennen
  colnames(data_IR) <- 1:I
  colnames(data_RT) <- 1:I
  
  ##### Startwerte für Parameter
  #IR-Modell
  alpha <- as.numeric(rep(1, I))
  ibeta <- apply(data_IR, 2, function(x) -log((length(x)-sum(x))/sum(x)))
  
  #RT-Modell
  ny <- apply(log(data_RT), 2, function(x) mean(x))
  omega <-  apply(log(data_RT), 2, function(x) 1/var(x))
  
  ##### Vorbereitung  Arrays und Dfs ####
  array_joint <- array(NA, dim = c(P, Q, I))
  
  posterior <- as.data.frame(matrix(NA, P, Q))
  
  r_iqr <- as.data.frame(matrix(NA, Q, I))
  
  se_itempara <- matrix(NA, I, 4)
  
  
  ######## EM-Algorithmus ########
  
  ## Vorbereitung EM-Werte
  alpha_old <- alpha
  ibeta_old <- ibeta
  ny_old <- ny
  omega_old <- omega
  SIGMA_G_old <- SIGMA_G
  EM <- T
  R <- 0
  
  while (EM) {
    R <- R + 1
    
    #### E-STEP ####
    print(paste("EM Iteration:", R))
    
    ## A-Posteriori-Verteilung
    #Gemeinsame Verteilung von y_pi und t_pi gegeben der Stützpunkte von lambda_p und der Itemparameter psi_i
    for (p in 1:P) {
      for (i in 1:I) {
        array_joint[p, , i] <-
          PQ_fun(data_IR[p, i], ibeta[i], alpha[i], X$points[, 1], D) * 
          rtime(data_RT[p, i], ny[i], X$points[, 2], omega[i])
      }
    }
    
    #Gemeinsame Wahrscheinlichkeit von y_p und t_p gegeben theta und itemparameter
    prod_PQr <-
      apply(array_joint, c(1, 2), function(i) prod(i))
    
    #A-Posteriori-Verteilung pro Person
    for (p in 1:P) {
      posterior[p, ] <-
        prod_PQr[p, ] * X$weights / sum(prod_PQr[p, ] * X$weights)
    }
    
    
    ##Brechnung von n_iqr
    n_iqr <- apply(posterior, 2, function(x) sum(x))
    
    
    ##Brechnung von r_iqr
    for (i in 1:I) {
      r_iqr[, i] <-
        apply(data_IR[, i] * posterior, 2, function(x) sum(x))
    }
    
    
    #### M-STEP ####
    
    #### Schätzung der Itemparameter ####
    ## Newton-Raphson-Verfahren
    for (i in 1:I) {
      newton_criterion <- T
      R2 <- 0
      
      #### Schätzung von alpha und ibeta ####
      while (newton_criterion) {
        R2 <- R2 + 1
        
        P_theta <- PQ_fun(1, ibeta[i], alpha[i], X$points[, 1], D)
        Q_theta <- PQ_fun(0, ibeta[i], alpha[i], X$points[, 1], D)
        
        ### Berechnunung von L1, L2, L3, L4 sowie L11, L12, L22, L33, L34, L44
        ## erste Ableitungen
        L1 <-
          D*sum((r_iqr[, i] - n_iqr * P_theta) * X$points[, 1])
        
        L2 <- D*sum(r_iqr[, i] - n_iqr * P_theta)
        
        
        ## zweite Ableitungen
        L11 <-
          -D^2*sum((n_iqr * P_theta * Q_theta) * (X$points[, 1]) ^ 2)
        
        L12 <-
          L21 <- -D^2*sum((n_iqr * P_theta * Q_theta) * X$points[, 1])
        
        L22 <- -D^2*sum(n_iqr * P_theta * Q_theta)
        
        #Hesse-Matrix erstellen
        Hes <- matrix(c(L11, L12,
                        L21, L22), 2, 2, byrow = T)
        #Newton-Raphson-Schritt anwenden
        res <- c(alpha[i], ibeta[i]) - solve(Hes) %*% c(L1, L2)
        #Erfüllung des Konvergenzkriteriums prüfen
        newton_criterion <-
          (abs(alpha[i] - res[1,]) >= newtonraphson_conv_criterion |
              abs(ibeta[i] - res[2,]) >= newtonraphson_conv_criterion) &
          newtonraphson_iter_criterion > R2
        #Update von Parametern
        alpha[i] <- res[1,]
        ibeta[i] <- res[2,]
        
      }
      #Warnung, Netwon-Raphson-Iterations-Kriterium erreicht wurde
      if(R2==newtonraphson_iter_criterion){
        warning("Newton-Raphson-Kriterium erreicht")
      }
      
      ## Standardfehler von alpha und beta
      se_itempara[i,1] <- sqrt(-solve(Hes)[1,1])
      se_itempara[i,2] <- sqrt(-solve(Hes)[2,2])
      
      
      #### Schätzung von ny und omega ####
      num_ny <- 0
      for (q in 1:Q) {
        num_ny <-
          num_ny + sum(posterior[, q] * (log(data_RT[, i]) + X$points[q, 2]))
      }
      ny[i] <- num_ny / P
      
      denom_omega <- 0
      for (q in 1:Q) {
        denom_omega <-
          denom_omega + sum(posterior[, q] * (log(data_RT[, i]) - (ny[i] - X$points[q, 2]))^2)
      }
      omega[i] <- sqrt(P / denom_omega)
      
      
      ## Standardfehler von omega und ny
      info_fisher <- matrix(c(omega[i]^-2*P,0,
                              0,omega[i]^2*P), 2, 2, byrow = T)
      
      se_itempara[i,3] <- sqrt(solve(info_fisher)[1,1])
      se_itempara[i,4] <- sqrt(solve(info_fisher)[2,2])
      
    }
    
    
    #### Schätzung von SIGMA ####
    
    num_sigma_12 <- num_sigma_22 <- 0
    for (q in 1:Q) {
      num_sigma_12 <-
        num_sigma_12 + sum(posterior[, q] * X$points[q, 1] * X$points[q, 2])
      num_sigma_22 <-
        num_sigma_22 + sum(posterior[, q] * X$points[q, 2] * X$points[q, 2])
    }
    SIGMA_G[2, 1] <-  SIGMA_G[1, 2] <- num_sigma_12 / P
    SIGMA_G[2, 2] <- num_sigma_22 / P
    
    
    ### Stützpunkte und Gewichtswerte anpassen
    X <-
      mgauss.hermite(Anz_QPunkte,
                     mu = c(0, 0),
                     sigma = SIGMA_G,
                     prune = prune)
    
    
    #### EM-Kriterium testen ####
    EM_alpha <- abs(alpha_old - alpha) >= EM_conv_criterion
    EM_ibeta <- abs(ibeta_old - ibeta) >= EM_conv_criterion
    EM_ny <- abs(ny_old - ny) >= EM_conv_criterion
    EM_omega <- abs(omega_old - omega) >= EM_conv_criterion
    EM_SIGMA_G <- abs(SIGMA_G_old - SIGMA_G) >= EM_conv_criterion
    
    alpha_old <- alpha
    ibeta_old <- ibeta
    ny_old <- ny
    omega_old <- omega
    SIGMA_G_old <- SIGMA_G
    
    EM <- (any(EM_alpha) | any(EM_ibeta) | any(EM_ny) |  
             any(EM_omega) |  any(EM_SIGMA_G)) & EM_iter_criterion > R
    #### M-Step Ende ####

  } #### EM-Algorithmus Ende ####
  
  # Itemparameter als Matrix
  itempara <- cbind(alpha, ibeta, omega, ny)
  colnames(itempara) <- c("Alpha", "Beta", "Omega", "Ny")
  
  mml.hm <- list("Item_Parameter" = itempara, "Standardfehler Itemparameter"=se_itempara, "Kovarianzmatrix" = SIGMA_G)
}


### mml.hm- Funktion speichern
save(mml.hm, file="functions/mml_hm.RData")

####################### End of Skript ###