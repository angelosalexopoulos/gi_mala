rm(list=ls())



# Function to plot the ergodic mean
plot_ergodic_mean <- function(samples) {
  n <- length(samples)
  ergodic_mean <- cumsum(samples) / seq_len(n)
  
  plot(seq_len(n), ergodic_mean, type = "l", col = "blue",
       xlab = "Iteration", ylab = "Ergodic Mean",
       main = "")
}


log_pnorm <- function(z) {
  # Numerically stable computation of log(pnorm(z))
  result <- numeric(length(z))
  low <- z < -10
  result[low] <- pnorm(z[low], log.p = TRUE)
  result[!low] <- log(pnorm(z[!low]))
  return(result)
}

set.seed(121088)
mala_t <- function(n_samples, nu, gam_gi, x0 = 0, MALA = FALSE, precond,b) {
  # Target log-density and gradient
  log_p <- function(x) {
    return(-((nu + 1) / 2) * log(1 + (x^2 / nu)))
    #return(dnorm( x ,0, 1,log = T))
  }
  
  grad_log_p <- function(x) {
    return(-((nu + 1) * x) / (nu + x^2))
    #return(-x)
  }
  
  # Initialize
  samples <- numeric(n_samples)
  adjustments <- matrix(NA,n_samples,6)  # Store computed adjustments
  x <- x0
  accept_count <- 0  # Track accepted moves
  
  beta <- 1 - gam_gi  # Define beta based on gamma
  Kn <- matrix(0,2,2)
  H1 <- H2 <- rep(NA,n_samples)
  for (i in 1:n_samples) {
    # Compute the Langevin proposal
    if(!MALA){
      #precond <- 1+x^2
      precond <-  mynu/(mynu +1)
      #precond <- mynu/(mynu-2)
      #precond <-1
    }
    grad <- grad_log_p(x)
    mu <- x + gam_gi * grad * precond
    
    if (MALA) {
      prop_var <- 2 * gam_gi
    } else {
      prop_var <- 2 * gam_gi - gam_gi^2
    }
    
    x_prop <- rnorm(1, mean = mu, sd = sqrt(precond * prop_var))
    
    # Compute the acceptance probability
    log_p_x <- log_p(x)
    log_p_x_prop <- log_p(x_prop)
    
    grad_prop <- grad_log_p(x_prop)
    precondprop<-precond
    if(!MALA){
      #precondprop <- 1+x_prop^2
      precondprop <- mynu/(mynu +1)
      #precondprop <- mynu/(mynu-2)
      #precondprop <- 1
    }
    mu_prop <- x_prop + gam_gi * grad_prop * precondprop
    
    log_q_x_prop_given_x <- dnorm(x_prop, mean = mu, sd = sqrt(precond * prop_var), log = TRUE)
    log_q_x_given_x_prop <- dnorm(x, mean = mu_prop, sd = sqrt(precondprop * prop_var), log = TRUE)
    
    log_alpha <- (log_p_x_prop + log_q_x_given_x_prop) - (log_p_x + log_q_x_prop_given_x)

    acc_ratio <- if (log_alpha >= 0) 1 else if (log_alpha < -700) 0 else exp(log_alpha)
    acc_ratio <- min(1, max(0, acc_ratio))
    
    #cat(acc_ratio,'\n')
    if (log(runif(1)) < log_alpha) {
      x <- x_prop  # Accept
      grad <- grad_prop
      accept_count <- accept_count + 1
      precond <- precondprop
    }
    
    samples[i] <- x

    if(!MALA){
      beta_pow <- beta^(1:5)
      beta_pow2 <- beta^(2*(1:5))
      
      vPhi2 <- precond * (1 - beta_pow2)
      vPhi <- sqrt(pmax(vPhi2, .Machine$double.eps))
      
      #approx_prob = 1-pnorm(b/sqrt(precond))
      #approx_prob = 1 - pt(b, df = mynu)
      z_Gx <- (beta_pow * x-b) / vPhi
      z_Gy <- (beta_pow * x_prop-b) / vPhi
      Gx <- as.numeric(x > b)+ sum(exp(log_pnorm(z_Gx)))
      Gy <- as.numeric(x_prop > b) + sum(exp(log_pnorm(z_Gy)))
      #Gx <- as.numeric(x > b) + sum(exp(log_pnorm(z_Gx)))
      #Gy <- as.numeric(x_prop > b) + sum(exp(log_pnorm(z_Gx)))
      
      muq <- x + gam_gi * precond * grad
      s2q <- pmax(precond * (2 * gam_gi - gam_gi^2), .Machine$double.eps)
      sqrt_s2q <- sqrt(s2q)
      
      Eqa <- exp(log_pnorm((muq-b) / sqrt_s2q)) 
      
      denom2 <- vPhi2 + beta_pow2 * s2q
      denom <- sqrt(pmax(denom2, .Machine$double.eps))
      z_Eqb <- (beta_pow * muq-b) / denom
      Eqb <- sum(exp(log_pnorm(z_Eqb)))
      
      
      adjustments[i,] = ((beta^c(0,1:5))*muq-b)/sqrt(pmax(precond * (1 - beta^(2*c(0,1:5))) +   beta^(2*c(0,1:5))* s2q, .Machine$double.eps))
      H1[i] <- acc_ratio*(Gy-Gx)
      H2[i] <- Gy - Eqa - Eqb
      
      Kn[1,1] <- Kn[1,1] + H1[i]^2
      Kn[1,2] <- Kn[1,2] + H1[i]*H2[i]
      Kn[2,1] <- Kn[2,1] + H2[i]*H1[i]
      Kn[2,2] <- Kn[2,2] + H2[i]^2
      
      #adjustments[i] <- acc_ratio
    }
  }
  
  # Compute and return acceptance rate
  accept_rate <- accept_count / n_samples
  #cat("Acceptance Rate:", accept_rate, "\n")
  
  return(list(samples = samples,H1 = H1,H2 = H2,adjustments = adjustments, accept_rate = accept_rate,Kn=Kn/n_samples))
}

mynu_list <- c(1,2,5,30,100,300,1000)
c_list <- c(0, 1, 2, 3)

ratio_matrix <- matrix(NA, nrow = length(mynu_list), ncol = length(c_list),
                       dimnames = list(paste0("nu=", mynu_list), paste0("c=", c_list)))


for (i in seq_along(mynu_list)) {
  for (j in seq_along(c_list)) {
    mynu <- mynu_list[i]
    c_value <- c_list[j]
    myprecond <- mynu /(mynu +1)
    
    tail_prob_plain_mala <- rep(NA,100)
    tail_prob_plain_mala_precond <- rep(NA,100)
    tail_prob_gi <- rep(NA,100)
    tail_prob_gi_vr <- rep(NA,100)
    
    for(rep_iter in 1:100){
      #result_plain_mala <- mala_t(n_samples = 10000, nu = mynu, gam_gi = 2.2, MALA = TRUE, precond = 1, b = c_value)
      #result_plain_mala_precond <- mala_t(n_samples = 10000, nu = mynu, gam_gi = 1.7, MALA = TRUE, precond = myprecond, b = c_value)
      result_gi <- mala_t(n_samples = 10000, nu = mynu, gam_gi = 0.85, MALA = FALSE, precond = myprecond, b = c_value)
      
      #tail_prob_plain_mala[rep_iter] <- mean(result_plain_mala$samples > c_value)
      #tail_prob_plain_mala_precond[rep_iter] <- mean(result_plain_mala_precond$samples > c_value)
      tail_prob_gi[rep_iter] <- mean(result_gi$samples > c_value)
      
      vec1 <- c(mean((result_gi$samples > c_value) * result_gi$H1), mean((result_gi$samples > c_value) * result_gi$H2))
      vec2 <- c(mean(result_gi$samples > c_value) * mean(result_gi$H1), mean(result_gi$samples > c_value) * mean(result_gi$H2))
      betaVR <- solve(result_gi$Kn) %*% (vec1 - vec2)
      betaVR[2,1] =1
      betaVR[1,1] = -1 
      tail_prob_gi_vr[rep_iter] <- mean(result_gi$samples > c_value) +
        betaVR[2,1] * mean(result_gi$H1) +
        betaVR[1,1] * mean(result_gi$H2)
    }
    
    ratio_matrix[i, j] <- var(tail_prob_gi) / var(tail_prob_gi_vr)
    cat(sprintf("Done for nu = %d, c = %d\n", mynu, c_value))
  }
}

cat("\\begin{tabular}{l", paste(rep("c", length(c_list)), collapse = ""), "}\n", sep = "")
cat(" & ", paste0("c=", c_list, collapse = " & "), " \\\\\n\\hline\n")
for (i in seq_along(mynu_list)) {
  cat(paste0("nu=", mynu_list[i]), " & ",
      paste(sprintf("%.3f", ratio_matrix[i, ]), collapse = " & "), " \\\\\n")
}
cat("\\end{tabular}\n")

