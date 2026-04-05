#---------------------------Description---------------------------------------------------#
#  Input:                                                                                 #
#                X: High-frequency sample of fractional Gaussian Noise, for a > 0 that is #        
#                   X = a*(B_{1/n}^H,B_{2/n}^H - B_{1/n}^H,...,B_{n/n}^H - B_{(n-1)/n}^H) #
#                                                                                         #
#                K: Number of terms in Paxson's approximation,                            #
#                   Default: K=50 (not relevant for method "MLE")                         #            
#                                                                                         #
#           method: Estimation method to be used: "MLE", "Whittle", "OSMLE" or "OSSE"     #
#                                                                                         #
#          epsilon: truncation parameter for OSMLE and OSSE methods,                      #
#                   Default: epsilon = 10^{-3}                                            #
#                                                                                         #
#  Output:       Estimated Hurst parameter H, assuming a > 0 is unknown                   #
#                                                                                         #
#  Example:      Estimator(c(0.35,0.37,1.24,-0.75,-0.054), method="MLE")                  #
#                                                                                         #
#-----------------------------------------------------------------------------------------#
#             J.Reuß          04/2026                                                     #
#-----------------------------------------------------------------------------------------#

Estimator <- function(X, K=50, method, epsilon=1e-3) {
  
  # MLE
  if(method == "MLE") {
    n <- length(X) # Number of observations
    
    # Autocovariance function and its derivative
    gamma <- function(H) {
      k <- 0:(n-1)
      e <-2*H
      
      return((abs(k+1)^e - 2*abs(k)^e + abs(k-1)^e)/2)
    }
    
    dH_gamma <- function(H) {
      k <- 0:(n-1)
      e <-2*H
      
      t1 <- (k+1)^e * log(k+1)
      t2 <- ifelse(k>0, k^e * log(k), 0)
      t3 <- ifelse(abs(k-1)>0, abs(k-1)^e * log(abs(k-1)), 0)
      
      return(t1 - 2*t2 + t3)
    } 
    
    # Target function for root-finding
    l <- function(H) {
      Tn <- toeplitz(gamma(H))
      dTn <- toeplitz(dH_gamma(H))
      
      Tn_inv <- tryCatch({chol2inv(chol(Tn))}, error=function(e){solve(Tn)})
      
      w <- Tn_inv %*% X 
      
      return(sum(Tn_inv * dTn) - n*sum(w * (dTn %*% w))/sum(X * w))
    }
    
    # Find and return the MLE using uniroot
    return(uniroot(l, c(0.01,0.99), extendInt = "yes")$root)
  }
  
  # Whittle
  else if(method == "Whittle") {
    
    n <- length(X) # Number of observations
    
    # If n < 2 the estimator is not defined
    if(n < 2) stop("The Whittle estimator requires a minimum of two observations.")
    
    
    m <- floor(n/2)
    lambda <- 2*pi*(1:m)/n
    
    # Fast Fourier Transformation 
    I <- fft(X)
    I_n <- abs(I[2:(m+1)])^2/(2*pi*n)
    
    V <- 2*pi*(1:K)
    W <- 2*pi*(K+1)
    
    # Paxson's approximation of the spectral density f_H
    Pax_f_H <- function(H, x) {
      e1 <- 2*H+1
      e2 <- 2*H
      C_H <- gamma(e1)*sin(pi*H)/(2*pi)
      
      Sum <- colSums(outer(V, x, "+")^(-e1) + outer(V, x, "-")^(-e1))
      
      a_H_Sum <- ((V[K] + x)^(-e2) + (V[K] - x)^(-e2) +  
                    (W + x)^(-e2) + (W - x)^(-e2))/(8*pi*H)
      
      return(C_H*2*(1 - cos(x)) * (abs(x)^(-e1) + Sum + a_H_Sum))
    }
    
    # Target function for maximization
    r <- function(H) {
      f_lambda <- Pax_f_H(H, lambda)
      
      norm_fac <- exp(integrate(function(y) log(Pax_f_H(H, y)), lower=0, 
                                upper=pi, stop.on.error=FALSE)$value/pi)
      
      return(2*sum(I_n*norm_fac/f_lambda)/n)
    }
    
    # Locate the minimum by using optimize
    return(optimize(r, c(0.01,0.99))$minimum)
  }
  
  # OSMLE
  else if(method == "OSMLE") {
    
    n <- length(X) # Number of observations
    
    # If n < 4 the estimator is not defined
    if(n < 4) stop("The OSMLE requires a minimum of four observations.")
    
    Delta_n = 1/n
    
    # Calculate the initial guess estimators
    V_n1 <- sum(diff(X)^2)/(n-1)
    V_n2 <- sum((X[4:n] + X[3:(n-1)] - X[2:(n-2)] - X[1:(n-3)])^2)/(n-3)
    
    tilde_H <- max(epsilon, log2(V_n2/V_n1)/2) 
    tilde_sigma <- sqrt( V_n1/(Delta_n^(2*tilde_H) * (4 - 2^(2*tilde_H))) )
    
    e1 <- 2*tilde_H
    e2 <- e1 + 1
    
    # Calculate the covariance matrix, its inverse and derivative
    k <- 0:(n-1)
    Tn <- toeplitz((abs(k+1)^e1 - 2*abs(k)^e1 + abs(k-1)^e1)/2)
    
    t1 <- (k+1)^e1 * log(k+1)
    t2 <- ifelse(k>0, k^e1 * log(k), 0)
    t3 <- ifelse(abs(k-1)>0, abs(k-1)^e1 * log(abs(k-1)), 0)
    dTn <- toeplitz(t1 - 2*t2 + t3)
    
    Tn_inv <- tryCatch({chol2inv(chol(Tn))}, error=function(e){solve(Tn)})
    w <- Tn_inv %*% X 
    c <- tilde_sigma^2 * Delta_n^e1
    
    # A_n and B_n
    A_n <- sqrt(n) * (sum(X * w)/(n*c) - 1)
    B_n <- (sum(Tn_inv * dTn)/2 - sum(w * (dTn %*% w))/(2*c))/sqrt(n)
    
    # Paxson's approximation of f_H and its derivative
    C_H <- gamma(e2) * sin(pi*tilde_H)/(2*pi)
    C_H2 <- gamma(e2)/(2*pi) * (2*digamma(e2)*sin(pi*tilde_H) + pi*cos(pi*tilde_H))
    V <- (1:K)*2*pi
    W <- 2*pi*(K+1)
    
    h <- function(v, z) {
      -(log(v + z) + 1/e1)/(2*pi*tilde_H * (v + z)^(e1))
    }
    
    # Calculate the integrand
    Sums <- function(x) {
      m1 <- outer(V, x, "+")
      m2 <- outer(V, x, "-")
      
      Sum <- colSums(m1^(-e2) + m2^(-e2))
      Sum2 <- colSums(-2 * (log(m1) * m1^(-e2) + log(m2) * m2^(-e2)))
      
      aH_avg <- ((V[K] + x)^(-e1) + (V[K] - x)^(-e1) +  
                   (W + x)^(-e1) + (W - x)^(-e1))/(8*pi*tilde_H)
      
      aH2_avg <- (h(V[K], x) + h(V[K], -x) + h(W, x) + h(W, -x))/2
      
      list(Sum=Sum, Sum2=Sum2, aH_avg=aH_avg, aH2_avg=aH2_avg)
    }
    
    integrand <- function(x) {
      s <- Sums(x)
      fac <- 2*(1 - cos(x))
      G <- (abs(x)^(-e2) + s$Sum + s$aH_avg)
      
      f_H <- C_H * fac * G
      df_H <- C_H2 * fac * G + 
        C_H * fac * (-2*log(abs(x)) * abs(x)^(-e2) + s$Sum2 + s$aH2_avg)
      return(df_H/f_H)
    }
    
    # Calculate the integrals
    L_1 <- integrate(function(x) integrand(x)^2, lower=0, 
                     upper=pi, stop.on.error = FALSE)$value/pi
    L_2 <- integrate(function(x) integrand(x), lower=0, 
                     upper=pi, stop.on.error = FALSE)$value/pi
    
    # Return the OSMLE
    return(tilde_H - (2*B_n + A_n*L_2)/(sqrt(n) * L_1 - sqrt(n) * L_2^2))
  }
  
  # OSSE
  else if(method == "OSSE") {
    
    n <- length(X) # Number of observations
    
    # If n < 2 the estimator is not defined
    if(n < 2) stop("The OSSE requires a minimum of two observations.")
    
    Delta_n = 1/n
    
    # Calculate the initial guess estimators
    X_2 <- X[1:(n-1)] + X[2:n]
    
    tilde_H <- max(epsilon, 1/(2*log(2)) * log(sum(X_2^2)/sum(X^2)))
    tilde_sigma <- sqrt(n^(2*tilde_H - 1) * sum(X^2))
    
    e1 <- 2*tilde_H
    e2 <- e1 + 1
    
    # Calculate the covariance matrix, its inverse and derivative
    k <- 0:(n-1)
    Tn <- toeplitz((abs(k+1)^e1 - 2*abs(k)^e1 + abs(k-1)^e1)/2)
    
    t1 <- (k+1)^e1 * log(k+1)
    t2 <- ifelse(k>0, k^e1 * log(k), 0)
    t3 <- ifelse(abs(k-1)>0, abs(k-1)^e1 * log(abs(k-1)), 0)
    dTn <- toeplitz(t1 - 2*t2 + t3)
    
    Tn_inv <- tryCatch({chol2inv(chol(Tn))}, error=function(e){solve(Tn)})
    w <- Tn_inv %*% X 
    c <- tilde_sigma^2 * Delta_n^e1
    
    # A_n and B_n
    A_n <- sqrt(n) * (sum(X * w)/(n*c) - 1)
    B_n <- (sum(Tn_inv * dTn)/2 - sum(w * (dTn %*% w))/(2*c))/sqrt(n)
    
    # Paxson's approximation of f_H and its derivative
    C_H <- gamma(e2) * sin(pi*tilde_H)/(2*pi)
    C_H2 <- gamma(e2)/(2*pi) * (2*digamma(e2)*sin(pi*tilde_H) + pi*cos(pi*tilde_H))
    V <- (1:K)*2*pi
    W <- 2*pi*(K+1)
    
    h <- function(v, z) {
      -(log(v + z) + 1/e1)/(2*pi*tilde_H * (v + z)^(e1))
    }
    
    # Calculate the integrand
    Sums <- function(x) {
      m1 <- outer(V, x, "+")
      m2 <- outer(V, x, "-")
      
      Sum <- colSums(m1^(-e2) + m2^(-e2))
      Sum2 <- colSums(-2 * (log(m1) * m1^(-e2) + log(m2) * m2^(-e2)))
      
      aH_avg <- ((V[K] + x)^(-e1) + (V[K] - x)^(-e1) +  
                   (W + x)^(-e1) + (W - x)^(-e1))/(8*pi*tilde_H)
      
      aH2_avg <- (h(V[K], x) + h(V[K], -x) + h(W, x) + h(W, -x))/2
      
      list(Sum=Sum, Sum2=Sum2, aH_avg=aH_avg, aH2_avg=aH2_avg)
    }
    
    integrand <- function(x) {
      s <- Sums(x)
      fac <- 2*(1 - cos(x))
      G <- (abs(x)^(-e2) + s$Sum + s$aH_avg)
      
      f_H <- C_H * fac * G
      df_H <- C_H2 * fac * G + 
        C_H * fac * (-2*log(abs(x)) * abs(x)^(-e2) + s$Sum2 + s$aH2_avg)
      return(df_H/f_H)
    }
    
    # Calculate the integrals
    L_1 <- integrate(function(x) integrand(x)^2, lower=0, 
                     upper=pi, stop.on.error = FALSE)$value/pi
    L_2 <- integrate(function(x) integrand(x), lower=0, 
                     upper=pi, stop.on.error = FALSE)$value/pi
    
    # Return the OSSE
    return(tilde_H - (2*B_n + A_n*L_2)/(sqrt(n) * L_1 - sqrt(n) * L_2^2))
  }
  
  else print("There is no such method.")
}




#---------------------------Description---------------------------------------------------#
#  Input:                                                                                 #
#                H: Hurst parameter H for which the lower bound should be calculated      #
#                                                                                         #
#                K: Number of terms in Paxson's approximation                             #
#                                                                                         #
#  Output:       lower asymptotic bound of the second moment of sqrt(n)*(Estimator - H)   #
#                                                                                         #
#  Example:      asym_lowbound(H=0.3, K=1000)                                             #
#                                                                                         #
#-----------------------------------------------------------------------------------------#
#             J.Reuß          04/2026                                                     #
#-----------------------------------------------------------------------------------------#

asym_lowbound <- function(H, K) {
  # Hurwitz zeta function and its derivative with respect to s using the 
  # integral representation
  Hzeta <- function(s, a) {
    vapply(a, function(a_i) {
      
      g <- function(x) {
        return((x^(s-1) * exp(-a_i*x))/(1-exp(-x)))
      }
      result <- integrate(g, lower=0, upper=Inf, stop.on.error = FALSE)
      (1/gamma(s))*result$value
      
    }, numeric(1))
  }
  
  ds_hzeta <- function(s, a) {
    vapply(a, function(a_i) {
      
      g <- function(x) {
        return((exp(-a_i*x) * x^(s-1) * (log(x) - digamma(s)))/(1-exp(-x)))
      }
      result <- integrate(g, lower=0, upper=Inf, stop.on.error = FALSE)
      (1/gamma(s))*result$value
      
    }, numeric(1))
  }
  
  # Spectral density f_H and its derivative 
  f <- function(H,x) {
    gamma(2*H+1)*sin(pi*H)/(2*pi) * 2*(1-cos(x))*(2*pi)^(-1-2*H) * 
      (Hzeta(1+2*H,1-x/(2*pi)) + Hzeta(1+2*H,x/(2*pi)))
  }
  
  partialH_f_H <- function(H,x) {
    
    partial_tildeC_H <- function(H) {
      return((2*pi)^(-2*H-2) * gamma(2*H+1)*(2*digamma(2*H+1)*sin(pi*H) + 
                                               pi*cos(pi*H) - 2*log(2*pi)*sin(pi*H)))
    }
    
    tildeC_H <- function(H) {
      return(gamma(2*H+1) * sin(pi*H) * (2*pi)^(-2*H-2))
    }
    
    Z_H <- function(H,x) {
      Hzeta(2*H + 1, 1-x/(2*pi)) + Hzeta(2*H + 1, x/(2*pi))
    }
    
    partial_Z_H <- function(H,x) {
      2*ds_hzeta(2*H + 1, 1-x/(2*pi)) + 2*ds_hzeta(2*H + 1, x/(2*pi))
    }
    
    return(2*(1-cos(x)) * (partial_tildeC_H(H) * Z_H(H,x) + 
                             tildeC_H(H) * partial_Z_H(H,x) ))
  } 
  
  # Paxson's approximation of f_H and its derivative
  Pax_f_H <- function(H, x, K) {
    e1 <- 2*H+1
    e2 <- 2*H
    V <- 1:K
    
    C_H <- gamma(e1) * sin(pi*H)/(2*pi)
    
    sapply(x, function(y) {
      Sum <- sum( (2*V*pi + y)^(-e1) + (2*V*pi - y)^(-e1) )
      
      a_H <- function(u) ((2*u*pi + y)^(-e2) + (2*u*pi - y)^(-e2) )/(4*pi*H)
      
      return(C_H * 2*(1 - cos(y)) * ( abs(y)^(-e1) + Sum + (a_H(K) + a_H(K+1))/2 ))
    })
  }
  
  Pax_partial_f_H <- function(H, x, K) {
    e1 <- 2*H+1
    e2 <- 2*H
    V <- 1:K
    
    C_H <- gamma(e1) * sin(pi*H)/(2*pi)
    C_H2 <- gamma(e1)/(2*pi)*(2*digamma(e1)*sin(pi*H) + pi*cos(pi*H))
    
    sapply(x, function(y) {
      Sum <- sum( (2*V*pi + y)^(-e1) + (2*V*pi - y)^(-e1) )
      Sum2 <- sum( -2*(log(2*V*pi + y)*(2*V*pi + y)^(-e1) + 
                         log((2*V*pi - y))*(2*V*pi - y)^(-e1)) )
      
      a_H <- function(u) ((2*u*pi + y)^(-e2) + (2*u*pi - y)^(-e2) )/(4*pi*H)
      a_H2 <- function(u) {
        h <- function(u, z) -(log(2*u*pi + z) + (e2)^(-1))/(2*pi*H * (2*u*pi + z)^(e2))
        return(h(u, y) + h(u, -y))
      }
      
      S1 <- C_H2 * 2*(1-cos(y)) * ( abs(y)^(-e1) + Sum + (a_H(K) + a_H(K+1))/2 )
      S2 <- C_H * 2*(1-cos(y)) * ( -2*log(abs(y))*abs(y)^(-e1) + Sum2 + 
                                     (a_H2(K) + a_H2(K+1))/2 )
      
      return(S1 + S2)
    })
  }
  
  # Calculate the entries of the Fisher information by combining the analytic 
  # representation and Paxson's approximation
  JH22 <- function(H, K) {
    if(H < 0.006) {
      I <- integrate(function(x) (Pax_partial_f_H(H, x, K)/Pax_f_H(H, x, K))^2, 
                     lower=0, upper=pi, stop.on.error = FALSE)$value
      return(I/(2*pi))
    }
    else {
      I_1 <- integrate(function(x) (partialH_f_H(H, x)/f(H, x))^2, 
                       lower=1e-3, upper=pi, stop.on.error=FALSE)$value 
      I_2 <- integrate(function(x) (Pax_partial_f_H(H, x, K)/Pax_f_H(H, x, K))^2, 
                       lower=0, upper=1e-3, stop.on.error = FALSE)$value
      return((I_1 + I_2)/(2*pi))
    }
  }
  
  JH21 <- function(H, K) {
    if(H < 0.006) {
      I <- integrate(function(x) Pax_partial_f_H(H, x, K)/Pax_f_H(H, x, K), 
                     lower=0, upper=1e-3, stop.on.error = FALSE)$value
      return(-I/pi)
    }
    else {
      I_1 <- integrate(function(x) partialH_f_H(H, x)/f(H, x), 
                       lower=1e-3, upper=pi, stop.on.error=FALSE)$value  
      I_2 <- integrate(function(x) Pax_partial_f_H(H, x, K)/Pax_f_H(H, x, K), 
                       lower=0, upper=1e-3, stop.on.error = FALSE)$value
      return(-(I_1 + I_2)/pi)
    }
  }
  
  # Return the asymptotic lower bound
  return(2/(2*JH22(H, K) - JH21(H, K)^2))
}
