## auxiliary functions ##
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol



#####-- Binary Response
#Wei's allocation
Rho_fun1 = function(z_RK) {
  # z_RK: a K-dim vector representing the parameters \theta_1, ...\theta_K
  temp = 1 / (1 - z_RK)
  return(temp / sum(temp))
}

#Neyman allocation
Rho_fun2 = function(z_RK) {
  # z_RK: a K-dim vector representing the parameters \theta_1, ...\theta_K
  temp = sqrt(z_RK * (1 - z_RK))
  return(temp / sum(temp))
}

#Rosenberger allocation
Rho_fun3 = function(z_RK) {
  # z_RK: a K-dim vector representing the parameters \theta_1, ...\theta_K
  temp = sqrt(z_RK)
  return(temp / sum(temp))
}

#####  -- Gaussian Response
#Zhang-Rosenberger's allocation
Rho_fun1_Gaussian = function(Mean_RK,SD_RK) {
  # Mean_RK: a K-dim vector representing the mean parameters
  # SD_RK: a K-dim vector representing the standaRK deviation parameters
  if(length(Mean_RK)!=2 | length(SD_RK)!=2) stop("Zhang-Rosenberger allocation is applicable to 2-arm design only.")
  if(any(Mean_RK<0)) return(c(0.5,0.5))
  m1s = sqrt(Mean_RK[1]); s1 = SD_RK[1]
  m2s = sqrt(Mean_RK[2]); s2 = SD_RK[2]
  r = s1*m2s/s2/m1s
  if((m1s<m2s & r>1)|(m1s>m2s & r<1)){
    return(c(m2s*s1,m1s*s2)/(m2s*s1+m1s*s2))
  }else{
    return(c(0.5,0.5))
  }
}


#Neyman allocation
Rho_fun2_Gaussian = function(Mean_RK,SD_RK) {
  # Mean_RK: a K-dim vector representing the mean parameters
  # SD_RK: a K-dim vector representing the standaRK deviation parameters
  temp = SD_RK
  return(temp / sum(temp))
}


#allocation function
g_fun = function(x_RK,
                 y_RK,
                 alpha = 2,
                 L = 1000) {
  temp = sapply(y_RK * (y_RK / x_RK) ^ alpha, function(x)
    min(x, L))
  return(temp / sum(temp))
}

## tests
test_oneside = function(S_RK, N_RK, sig_level) {
  # one-sided; test H0: p_1>p_k  for k=2,...,K
  K = length(S_RK)
  stat = rep(1, K - 1)

  # estimate = c(S_RK / N_RK)[2:K] - c(S_RK / N_RK)[1]

  for (k in 2:K) {
    pstar = (S_RK[1] * N_RK[1] + S_RK[k] * N_RK[k]) / (N_RK[1] + N_RK[k])
    stat[k - 1] = (S_RK[1] - S_RK[k]) /
      sqrt(pstar * (1 - pstar) * (1 /N_RK[1] + 1 / N_RK[k]))

  }
  pval = stats::pnorm(stat)
  # if(is.na(pval)) pval=1
  pval <- sapply(pval, function(x){if(is.na(x)) return(1) else return(x)})
  ### Will be needed in package development if deemed necessary
  # testmethod <- "Unadjusted One-sided Test"
  # names(estimate) <- "difference for treatment effect"
  # names(pval) <- "p.value"

  rval = list(p.value = pval)

  return(rval)
}

test_chisq = function(S_RK, N_RK, sig_level) {
  test = stats::chisq.test(x = rbind(round(S_RK * N_RK), N_RK - round(S_RK * N_RK)),
                    simulate.p.value = T)
  pval = test$p.value
  if(is.na(pval)) pval=1

  rval = list(p.value = pval)

  return(rval)
}

## tests - one-sided T
test_WelchT = function(X_new,sig_level,var.equal=F) {
  # one-sided; test H0: mu_1>mu_k  for k=2,...,K
  K = ncol(X_new)
  pval = rep(1, K - 1)

  for (k in 2:K) {
    pval[k-1] = stats::t.test(x=X_new[,k],y=X_new[,1],
                              alternative="greater",var.equal=F)$p.value

  }
  rval = list(p.value = pval)

  return(rval)
}


## tests - ANOVA
test_aov = function(X_new,sig_level) {
  K = ncol(X_new)
  data_long = na.omit(cbind(measure = as.vector(X_new),
                            arm = rep(1:K,each = nrow(X_new))))
  test = stats::anova(lm(measure~arm, data=as.data.frame(data_long)))
  pval = test$Pr[1]

  rval = list(p.value = pval)

  return(rval)
}
