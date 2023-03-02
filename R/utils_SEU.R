## SEU auxiliary functions ##
# fillNAgaps <- function(x, firstBack=FALSE) {
#   ## NA's in a vector or factor are replaced with last non-NA values
#   ## If firstBack is TRUE, it will fill in leading NA's with the first
#   ## non-NA value. If FALSE, it will not change leading NA's.
#
#   # If it's a factor, store the level labels and convert to integer
#   lvls <- NULL
#   if (is.factor(x)) {
#     lvls <- levels(x)
#     x    <- as.integer(x)
#   }
#
#   goodIdx <- !is.na(x)
#
#   # These are the non-NA values from x only
#   # Add a leading NA or take the first good value, depending on firstBack
#   if (firstBack)   goodVals <- c(x[goodIdx][1], x[goodIdx])
#   else             goodVals <- c(NA,            x[goodIdx])
#
#   # Fill the indices of the output vector with the indices pulled from
#   # these offsets of goodVals. Add 1 to avoid indexing to zero.
#   fillIdx <- cumsum(goodIdx)+1
#
#   x <- goodVals[fillIdx]
#
#   # If it was originally a factor, convert it back
#   if (!is.null(lvls)) {
#     x <- factor(x, levels=seq_along(lvls), labels=lvls)
#   }
#
#   x
# }


#####-- Binary Response
#randomized play-the-winner rule (RPW)
add_rule1 = function(X.df, arms) {
  #X.df: A data frame of two component: treatment arm and response value.
  #arms: a list of arm names

  K = length(arms)
  if(is.null(X.df)) return(rep(0,K))

  arm = factor(X.df[,1],levels=arms)
  response = X.df[,2]
  X_draw = model.matrix(~ arm - 1)
  colnames(X_draw) =  substr(colnames(X_draw),4, nchar(colnames(X_draw)))

  Dmat.func = function(resp) {
    if(resp==T) Dmat = diag(K) else Dmat = (1-diag(K))/(K-1)
    return(Dmat)
  }
  dot_func = function(a,b) a%*%b

  Dmat_list = lapply(response,Dmat.func)
  Xdraw_list = as.list(data.frame(t(X_draw)))
  Addition = t(mapply(dot_func, Xdraw_list, Dmat_list))


  return(colSums(Addition))
}


#randomized SEU rule: targeting Neyman
add_rule2 = function(X.df, arms) {
  #X.df: A data frame of two component: treatment arm and response value.
  #arms: a list of arm names

  K = length(arms)
  n = nrow(X.df)
  if(is.null(X.df)) return(rep(0,K))

  tmp = rbind(expand.grid(ARM=arms,RESPONSE=c(1/K)), X.df) #to adjust p_hat
  tmp = cbind(subject=seq(1:nrow(tmp)),tmp)
  X.df_wide = reshape2::dcast(tmp, subject ~ ARM, value.var="RESPONSE")[,-1]
  cumavg = function(x) vapply(seq_along(x), function(i) mean(x[1:i],na.rm = T), 1)

  p_hats = apply(X.df_wide,2,cumavg)[n+length(arms),]
  Addition_return = sqrt(p_hats*(1-p_hats))

  return(Addition_return)
}


#randomized SEU rule: targeting Neyman's allocation
add_rule3 = function(X.df, arms) {
  #X.df: A data frame of two component: treatment arm and response value.
  #arms: a list of arm names

  K = length(arms)
  n = nrow(X.df)
  if(is.null(X.df)) return(rep(0,K))

  tmp = rbind(expand.grid(ARM=arms,RESPONSE=c(1/K)), X.df) #to adjust p_hat
  tmp = cbind(subject=seq(1:nrow(tmp)),tmp)
  X.df_wide = reshape2::dcast(tmp, subject ~ ARM, value.var="RESPONSE")[,-1]
  cumavg = function(x) vapply(seq_along(x), function(i) mean(x[1:i],na.rm = T), 1)

  p_hats = apply(X.df_wide,2,cumavg)[n+length(arms),]
  Addition_return = sqrt(p_hats)

  return(Addition_return)
}


#randomized SEU rule: "randomized winner based on p"(???) allocation
add_rule4 = function(X.df, arms) {
  #X.df: A data frame of two component: treatment arm and response value.
  #arms: a list of arm names

  K = length(arms)
  n = nrow(X.df)
  if(is.null(X.df)) return(rep(0,K))

  tmp = rbind(expand.grid(ARM=arms,RESPONSE=c(1/K)), X.df) #to adjust p_hat
  tmp = cbind(subject=seq(1:nrow(tmp)),tmp)
  X.df_wide = reshape2::dcast(tmp, subject ~ ARM, value.var="RESPONSE")[,-1]
  cumavg = function(x) vapply(seq_along(x), function(i) mean(x[1:i],na.rm = T), 1)

  p_hats = apply(X.df_wide,2,cumavg)[n+length(arms),]
  Addition_return = rep(0.4/K,K) + 0.6*(p_hats==max(p_hats))

  return(Addition_return)
}


#randomized SEU rule: targeting Neyman
#It assumes observations from X.df are coming one by one, and the urn will
#update one by one. This is not the ideal case.
# add_rule2 = function(X.df, arms) {
#   #X.df: A data frame of two component: treatment arm and response value.
#   #arms: a list of arm names
#
#   K = length(arms)
#   n = nrow(X.df)
#   if(is.null(X.df)) return(rep(0,K))
#
#   tmp = rbind(X.df,expand.grid(ARM=arms,RESPONSE=c(0,1)))
#   tmp = cbind(subject=seq(1:nrow(tmp)),tmp)
#   X.df_wide = reshape2::dcast(tmp, subject ~ ARM, value.var="RESPONSE")[,-1]
#
#   X.df_cumcount = apply(X.df_wide,2,function(x) cumsum(ifelse(is.na(x), 0, 1)) + x*0)
#   X.df_cumcount = apply(rbind(0,X.df_cumcount), 2, zoo::na.locf)[-1,]
#   X.df_wide[is.na(X.df_wide)] = 0
#
#   p_hats = (apply(X.df_wide,2,cumsum) + 1/K) / (X.df_cumcount+1)
#   Addition = sqrt(p_hats*(1-p_hats))
#   if(n==1) Addition_return = Addition[1,] else Addition_return = colSums(Addition[1:n,])
#
#   return(Addition_return)
# }
#
# #randomized SEU rule: targeting Neyman's allocation
# add_rule3 = function(X.df, arms) {
#   #X.df: A data frame of two component: treatment arm and response value.
#   #arms: a list of arm names
#
#   K = length(arms)
#   n = nrow(X.df)
#   if(is.null(X.df)) return(rep(0,K))
#
#   tmp = rbind(X.df,expand.grid(ARM=arms,RESPONSE=c(0,1)))
#   tmp = cbind(subject=seq(1:nrow(tmp)),tmp)
#   X.df_wide = reshape2::dcast(tmp, subject ~ ARM, value.var="RESPONSE")[,-1]
#
#   X.df_cumcount = apply(X.df_wide,2,function(x) cumsum(ifelse(is.na(x), 0, 1)) + x*0)
#   X.df_cumcount = apply(rbind(0,X.df_cumcount), 2, zoo::na.locf)[-1,]
#   X.df_wide[is.na(X.df_wide)] = 0
#
#   p_hats = (apply(X.df_wide,2,cumsum) + 1/K) / (X.df_cumcount+1)
#   Addition = sqrt(p_hats)
#   if(n==1) Addition_return = Addition[1,] else Addition_return = colSums(Addition[1:n,])
#
#   return(Addition_return)
# }




#####  -- Gaussian Response

#Neyman allocation
add_rule1_Gaussian = function(X.df, arms) {
  #X.df: A data frame of two component: treatment arm and response value.
  #arms: a list of arm names

  K = length(arms)
  n = nrow(X.df)
  if(is.null(X.df)) return(rep(0,K))

  tmp = X.df
  tmp = cbind(subject=seq(1:nrow(tmp)),tmp)
  X.df_wide = reshape2::dcast(tmp, subject ~ ARM, value.var="RESPONSE")[,-1]
  cumsd = function(x) vapply(seq_along(x), function(i) sd(x[1:i],na.rm = T), 1)

  sd_hats = apply(X.df_wide,2,cumsd)[n,]
  Addition_return = sd_hats/sum(sd_hats)
  Addition_return[is.na(Addition_return)] = 0

  return(Addition_return)
}

#"randomized winner based on sd"(???) allocation
add_rule2_Gaussian = function(X.df, arms) {
  #X.df: A data frame of two component: treatment arm and response value.
  #arms: a list of arm names

  K = length(arms)
  n = nrow(X.df)
  if(is.null(X.df)) return(rep(0,K))

  tmp = X.df
  tmp = cbind(subject=seq(1:nrow(tmp)),tmp)
  X.df_wide = reshape2::dcast(tmp, subject ~ ARM, value.var="RESPONSE")[,-1]
  cumsd = function(x) vapply(seq_along(x), function(i) sd(x[1:i],na.rm = T), 1)

  sd_hats = apply(X.df_wide,2,cumsd)[n,]
  Addition_return = rep(0.4/K,K) + 0.6*(sd_hats==max(sd_hats))

  return(Addition_return)
}


## tests - one-sided T based on long data
test_WelchT_long = function(X_new,sig_level,var.equal=F) {
  # one-sided; test H0: mu_1>mu_k  for k=2,...,K
  K = length(levels(factor(X_new$ARM)))
  pval = rep(1, K - 1)

  for (k in 2:K) {
    pval[k-1] = stats::t.test(x=X_new[X_new$ARM==k,2],y=X_new[X_new$ARM==1,2],
                              alternative="greater",var.equal=F)$p.value

  }
  rval = list(p.value = pval)

  return(rval)
}

