library(tidyverse)
library(stringr)
library(reshape2)
library(splines)
library(nimble)
library(Matrix)
library(emulator)
library(MASS)
library(mvtnorm)
library(parallel)
rm(list = ls())
try(setwd("C:/Users/philaw/Box/Research/Joint_RT/code"))


dx_use_U = c(25,25,25,25)

model_output = mclapply(1:4,function(xx){
  
  dat = read.csv("data_submit.csv")[,-1]
  dat_use = dat[which(dat$family_GM %in% names(sort(table(dat$family_GM),decreasing = TRUE)[xx])),]
  # which(is.na(dat_use$Abundance))
  
  dat_use$family = dat_use$family_GM
  
  plot(dat_use$longitude,dat_use$latitude)
  
  dat_use = dat_use[order(dat_use$family_GM),]
  
  dim(dat_use)
  sort(table(dat_use$genus))
  sort(table(dat_use$family_GM))
  
  # dat_use = dat[dat$region %in% c("HTR","Cederberg") & dat$family_GM %in% c("ASTERACEAE"),]
  
  X = as.matrix(dat_use[,c("Elevation30m","Gmap","tminave01c","RFL_CONC","Abundance")])
  Y = as.matrix(log(dat_use[,c(which(colnames(dat_use) %in% c("lma","fwc","percent_N","succulence")),
                               grep("X",colnames(dat_use)))]))
  
  
  tmp1 = eigen(cor(scale(Y[,-(1:n_var)])))
  cumsum(tmp1$values / sum(tmp1$values))[1:10]
  
  tmp2 = eigen(cor((Y[,-(1:n_var)])))
  cumsum(tmp2$values / sum(tmp2$values))[1:10]
  
  
  r = Y[,-(1:n_var)]
  colnames(r) = 450:949
  rownames(r) = 1:nrow(r)
  
  r_temp = melt(r)
  colnames(r_temp) = c("Rep","Wavelength","Log Reflectance")
  r_temp$family = dat_use$family_GM[r_temp$Rep]
  r_temp$genus = dat_use$genus[r_temp$Rep]
  
  temp_use = r_temp %>% 
    group_by(family,Wavelength) %>% 
    summarize(mean = mean(`Log Reflectance`))
  
  # ggplot(data = temp_use) + 
  #   geom_line(aes(x = Wavelength,y = mean,group = family,color = family)) + 
  #   geom_vline(xintercept = c())
  # 
  # idx_use = which(r_temp$family == "ASTERACEAE" &
  #                   r_temp$genus %in% names(which(table(r_temp$genus)/500 >= 10)))
  # 
  # 
  # temp_use2 = r_temp[idx_use,] %>% 
  #   group_by(genus,Wavelength) %>% 
  #   summarize(mean = mean(exp(`Log Reflectance`)))
  # 
  # ggplot(data = temp_use2) + 
  #   geom_line(aes(x = Wavelength,y = mean,group = genus,color = genus))
  # X_use = X[idx_keep,] 
  # Y_use = Y[idx_keep,]  
  
  family_use = sort(unique(dat_use$family_GM))
  
  
  # idx_keep = which(dat_use$family_GM == "ASTERACEAE" &
  #                    dat_use$genus %in% names(which(table(dat_use$genus) >= 10)))
  # X_use = X[idx_keep,]
  # Y_use = Y[idx_keep,]
  X_use = scale(cbind(X))
  Y_use = Y
  
  
  # Refl = Y[,-(1:n_var)]
  # Tr = Y[,(1:n_var)]
  
  
  n_rep = nrow(Y_use)
  p = ncol(X_use)
  
  
  # n_family = length(family_use)
  # family_ind = as.numeric(as.factor(dat_use$family_GM))
  
  
  kern = function(x,v){
    dnorm(abs(x)/v)
  }
  
  dx_W = 10
  wave_knots_W = seq(450 ,950,by = dx_W)
  v_now_W = dx_W * 1.5
  dif_temp_W = outer(450:949,wave_knots_W,"-")
  K_W = cbind(1,kern(dif_temp_W,v_now_W))
  
  
  # K_W = cbind(1,bs(450:949,knots = seq(460,940,by = 10) ))
  n_basis_W = ncol(K_W)
  
  K_W_star = cbind(
    rbind(diag(n_var),matrix(0,n_wave,n_var)),
    rbind(matrix(0,n_var,n_basis_W),K_W)
  )
  n_basis_W = ncol(K_W)
  
  
  
  dx_beta = 100
  wave_knots_beta = seq(450,950,by = dx_beta)
  v_now_beta = dx_beta * 1.5
  dif_temp_beta = outer(450:949,wave_knots_beta,"-")
  K_beta = cbind(1,kern(dif_temp_beta,v_now_beta))
  
  # K_beta = cbind(1,bs(450:949,knots = c(seq(500,900,by = 100)) ))
  n_basis_beta = ncol(K_beta)
  
  K_beta_star = cbind(
    rbind(diag(n_var),matrix(0,n_wave,n_var)),
    rbind(matrix(0,n_var,n_basis_beta),K_beta)
  )
  
  
  # K_V = cbind(1,bs(450:949,knots = c(seq(475,675,by = 25),690,700,710,
  #                                    seq(725,750,by = 25),800,850,900) ))
  # n_basis_V = ncol(K_V)
  # 
  # K_V_star = cbind(
  #   rbind(diag(n_var),matrix(0,n_wave,n_var)),
  #   rbind(matrix(0,n_var,n_basis_V),K_V)
  # )
  # n_v = ncol(K_V_star)
  
  dx_U = dx_use_U[xx]
  wave_knots_U = seq(450,950,by = dx_U)
  v_now_U = dx_U * 1.5
  dif_temp_U = outer(450:949,wave_knots_U,"-")
  K_U = cbind(1,kern(dif_temp_U,v_now_U))
  
  
  # K_U = cbind(1,bs(450:949,knots = c(500,550,600,650,675,700,725,750,850) ))
  n_basis_U = ncol(K_U)
  
  K_U_star = rbind(matrix(0,n_var,n_basis_U),K_U)
  
  n_u = ncol(K_U_star)
  
  
  # K_psi = cbind(1,K_U)
  # n_basis_psi = ncol(K_psi)
  # K_alpha = K_psi
  # n_basis_alpha = ncol(K_alpha)
  
  
  
  K_eps = cbind(1,bs(450:949,knots = c(seq(475,925,by = 50)),degree = 1 ))
  n_eps = ncol(K_eps)
  
  temp_mod = lm(c(Y_use[,-(1:n_var)]) ~ 0 + K_W[rep(1:500,each = n_rep),] )
  summary(temp_mod)
  res_mat = matrix(resid(temp_mod),ncol = 500)
  plot(res_mat[1,])
  # 
  # temp_mod2 = lm(c(Y_use[,-(1:n_var)]) ~ 0 + K_W[rep(1:500,each = n_rep),]+ X_scale %x% K_beta)
  # summary(temp_mod2)
  # 
  # temp_mod3 = lm(c(Y_use[,-(1:n_var)]) ~ 0 + K_W[rep(1:500,each = n_rep),]+ cbind(X_scale,Y_use[,1:n_var]) %x% K_beta)
  # summary(temp_mod3)
  # 
  # 1- var(resid(temp_mod))/var(c(Y_use[,-(1:n_var)]))
  # 1- var(resid(temp_mod2))/var(c(Y_use[,-(1:n_var)]))
  # 1- var(resid(temp_mod3))/var(c(Y_use[,-(1:n_var)]))
  
  
  # anova(temp_mod,temp_mod2,temp_mod3)
  
  reps = 10e3
  burn = 10e4
  thin = 10
  tune = 100
  KX = K_beta_star %x% X_use
  
  XtX = t(X_use) %*% X_use
  
  tmp = coef(lm(apply(Y_use[,-(1:n_var)],2,mean) ~ 0 + K_beta ))
  # plot(K_beta %*% tmp)
  
  W_r = matrix(0,reps,n_basis_W); W_r_now = c(coef(temp_mod))
  W_t = matrix(0,reps,n_var); W_t_now = apply(Y_use[,(1:n_var)],2,mean)
  alp_now = c(K_W_star %*% c(W_t_now, W_r_now))
  
  beta_r = array(0,c(reps,n_basis_beta,p)); beta_r_now = matrix(0.0,n_basis_beta,p)
  
  # for(i in 1:n_family){
  #   beta_r_now[i,,1] = tmp
  # }
  
  # sig_pars_t = array(1,c(reps,n_var)); var_t_now = apply(Y_use[,1:n_var],2,var)
  log_sig_pars_R = matrix(0,reps,n_eps); log_sig_pars_R_now = rep(-1.5,n_eps)
  var_R_now = c(exp(K_eps %*% log_sig_pars_R_now))
  
  # gamma = array(0,c(reps,n_family,n_basis_gamma,n_var)); gamma_now = array(0,c(n_family,n_basis_gamma,n_var))
  # u_i = array(0,c(reps,n_family,n_u)); u_i_now = matrix(0,n_family,n_u)
  # psi_i_now =  t(K_psi_star %*% t(u_i_now))
  beta_t = array(0,c(reps,n_var,p)); beta_t_now = matrix(0,n_var,p)
  U_ij = array(0,c(reps,n_rep,n_basis_U)); U_ij_now = matrix(0,n_rep,n_basis_U)
  psi_ij_now = U_ij_now %*% t(K_U_star)
  
  # V_i = array(0,c(reps,n_family,n_v)); V_i_now = matrix(0,n_rep,n_v)
  # phi_i_now = V_i_now %*% t(K_V_star)
  
  sig2_W_save = numeric(reps); sig2_W_now = 1
  sig2_B_save = matrix(1,reps,p); sig2_B_now = rep(1,p)
  tau2_save = matrix(0,reps,n_wave); tau2_now = c(var_R_now)
  mu_save = array(0,c(reps,n_rep,n_var + n_wave))
  like_save = matrix(0,reps,2)
  Sigma_save = array(0,c(reps,n_u+n_var,n_var+n_u))
  # Sigma_inv_now = matrix(0,n_var+n_u,n_var+n_u)
  Sigma_now = diag(n_var+n_u)
  
  Sigma_now_proj = Sigma_now[-(1:n_var),(1:n_var)] %*% solve(Sigma_now[(1:n_var),(1:n_var)])
  
  
  V_now = Sigma_now[1:n_var,1:n_var]
  V_inv_now = solve(V_now)
  
  
  xb_r = X_use %*% t( K_beta %*% beta_r_now)
  xb_t = X_use %*% t(beta_t_now)
  
  xb = cbind(xb_t,xb_r)
  # xb = matrix(2.89607664,nrep,503)
  # ref_t = do.call("rbind",sapply(1:n_family,function(xx){
  #   Tr[family_ind_list[[xx]],] %*% t(K_gamma %*%  gamma_now[xx,,])
  # }))
  
  KtK_beta = bdiag(V_inv_now, quad.form(Matrix::diag(1/tau2_now),K_beta))
  KtK_U = quad.form(Matrix::diag(1/tau2_now),K_U)
  KtK_alpha = quad.form(Matrix::diag(1/tau2_now),K_W)
  
  
  
  mu_now = sweep(xb + psi_ij_now,2,alp_now,"+")
  
  like_func = function(mu,V,tau2){
    
    c(
      sum(dmvnorm(Y_use[,1:n_var]-mu[,1:n_var],sigma = V,log = TRUE)),
      sum(sapply(1:n_rep,function(jj){
        sum(dnorm(Y_use[jj,-(1:n_var)],mu[jj,-(1:n_var)],sqrt(tau2),log= TRUE))
      }))
    )
  }
  
  
  like_func_ind = function(mu,V,tau2){
    
    c(
      sum(dmvnorm(Y_use[,1:n_var]-mu[,1:n_var],sigma = V,log = TRUE)),
      sum(sapply(1:n_rep,function(jj){
        sum(dnorm(Y_use[jj,-(1:n_var)],mu[jj,-(1:n_var)],sqrt(tau2),log= TRUE))
      }))
    )
  }
  
  cand_v_tau = rep(0.2,n_eps)
  count_v_tau = rep(0,n_eps)
  
  
  like_now = sum(like_func(mu_now,V_now,tau2_now))
  
  st = proc.time()
  
  for(i in 2:(burn + reps * thin)){
    
    
    ### sample beta_r and beta_t
    mu_res = sweep(psi_ij_now,2,alp_now,"+")
    res_part = Y_use - mu_res
    res_adj = cbind(res_part[,1:n_var] %*% V_inv_now, 
                    sweep(res_part[,-(1:n_var)] ,2,tau2_now,"/"))
    
    
    v_beta = solve(KtK_beta %x% XtX + diag(c(rep(1e-3,n_var * (p) ),rep(1e-3,p), rep(1/sig2_B_now,times = (n_basis_beta - 1)))  ) ) 
    m_beta =    t(KX) %*% c(res_adj)
    beta_fam = mvrnorm(1,v_beta %*% m_beta,v_beta)
    
    beta_t_now = matrix(beta_fam[1:(n_var*p)],nrow = n_var,byrow = TRUE)
    beta_r_now = matrix(beta_fam[-(1:(n_var*p))],nrow = n_basis_beta,byrow = TRUE)
    
    xb_r = X_use %*% t( K_beta %*% beta_r_now)
    xb_t = X_use %*% t(beta_t_now)
    xb = cbind(xb_t,xb_r)
    
    ### sample alpha_r and alpha_t (This is W)
    
    mu_res = psi_ij_now + xb
    res_part = Y_use - mu_res
    res_adj = cbind(res_part[,1:n_var] %*% V_inv_now, 
                    sweep(res_part[,-(1:n_var)] ,2,tau2_now,"/"))
    
    
    v_alp = solve(n_rep * V_inv_now + 1e-3 * diag(n_var) )  
    m_alp = apply(res_adj[,1:n_var],2,sum)
    W_t_now = mvrnorm(1,v_alp %*% m_alp,v_alp)
    
    v_alp = solve(n_rep * KtK_alpha  +  diag(c(1e-3, rep(1/sig2_W_now,n_basis_W-1))))  
    m_alp = t(K_W) %*% apply(res_adj[,-(1:n_var)],2,sum)
    W_r_now = mvrnorm(1,v_alp %*% m_alp,v_alp)
    
    alp_now = c(K_W_star %*% c(W_t_now, W_r_now))
    
    mu_now = sweep(xb + psi_ij_now,2,alp_now,"+")
    
    res_now = Y_use - mu_now
    
    
    a_star = 1 + (n_basis_W - 1)/2
    b_star = 1 + 0.5 * sum(W_r_now[-1]^2)
    
    sig2_W_now = 1/rgamma(1,a_star,b_star)
    
    
    a_star_bet = 1 + (n_basis_beta - 1)/2
    
    sig2_B_now =  sapply(1:p, function(jj){
      b_star_bet  = 1 + 0.5 * sum(beta_r_now[-1,jj]^2)
      
      return(1/rgamma(1,a_star_bet,b_star_bet))
    })
    
    ### Sample U_{ij} for psi_{ij}
    
    mu_res = sweep(xb,2,alp_now,"+")
    res_part = Y_use - mu_res
    res_adj = cbind(res_part[,1:n_var] %*% V_inv_now, 
                    sweep(res_part[,-(1:n_var)] ,2,tau2_now,"/"))
    
    prior_inv = solve(Sigma_now[-(1:n_var),-(1:n_var)] - 
                        Sigma_now_proj %*% Sigma_now[(1:n_var),-(1:n_var)])
    prior_mean = Sigma_now_proj %*% t(res_now[,1:n_var])
    
    for(j in 1:n_rep){
      
      v_u = solve(KtK_U + prior_inv)
      m_u =  t(K_U) %*% c(res_adj[j,-(1:n_var)]) + prior_inv %*% prior_mean[,j]
      U_ij_now[j,] = mvrnorm(1,v_u %*% m_u,v_u)
      
    }
    
    psi_ij_now = U_ij_now %*% t(K_U_star)
    
    
    #### mean now 
    
    mu_now = sweep(xb + psi_ij_now,2,alp_now,"+")
    res_now = Y_use - mu_now
    
    ### Sigma_i
    U_star = cbind(res_now[,1:n_var],U_ij_now)
    
    nu_star = n_u + n_var+ 1 + n_rep
    V_sig = solve(diag(n_u + n_var)*1e-3  + crossprod(U_star))
    
    Sigma_now = solve(rWishart(1,nu_star,V_sig)[,,1])
    
    
    
    V_now = Sigma_now[1:n_var,1:n_var]
    V_inv_now = solve(V_now)
    Sigma_now_proj = Sigma_now[-(1:n_var),(1:n_var)] %*% V_inv_now
    
    ###  ###   ###  Sample Noise parameters 
    
    ### Deltas
    for(j in 1:n_eps){
      
      
      log_sig_pars_R_cand = log_sig_pars_R_now
      tau2_cand = tau2_now
      
      log_sig_pars_R_cand[j] = rnorm(1,log_sig_pars_R_now[j],cand_v_tau[j])
      
      tau2_cand = c(exp(K_eps %*% log_sig_pars_R_cand))
      
      if(j == 1){
        prior_dif = dnorm(log_sig_pars_R_cand[j],0,100,log = TRUE) - 
          dnorm(log_sig_pars_R_now[j],0,100,log = TRUE)
      } else{
        prior_dif = dnorm(log_sig_pars_R_cand[j],0,3,log = TRUE) - 
          dnorm(log_sig_pars_R_now[j],0,3,log = TRUE)
      }
      
      like_dif = sum(like_func(mu_now,V_now,tau2_cand) - 
                       like_func(mu_now,V_now,tau2_now))
      
      if(like_dif + prior_dif > log(runif(1))){
        
        tau2_now = tau2_cand
        log_sig_pars_R_now = log_sig_pars_R_cand
        
        count_v_tau[j] = count_v_tau[j] + 1
      }
      
    }
    
    # a_star = 1 + n_rep/2
    # b_star = 1 + 0.5 * apply(res_now,2,function(x){sum(x^2)})
    # 
    # tau2_now = sapply(1:n_wave,function(jj){
    #   1/rgamma(1,a_star,b_star[jj + n_var])
    # })
    
    KtK_beta = bdiag(V_inv_now, quad.form(Matrix::diag(1/tau2_now),K_beta))
    KtK_U = quad.form(Matrix::diag(1/tau2_now),K_U)
    KtK_alpha = quad.form(Matrix::diag(1/tau2_now),K_W)
    
    like_now = like_func(mu_now,V_now,tau2_now)
    
    if(i  %% tune == 0){
      
      if(i < burn){
        
        acc_v_tau = count_v_tau / tune; count_v_tau = rep(0,n_eps)
        
        cand_v_tau = ifelse(acc_v_tau > 0.6, cand_v_tau*2,
                            ifelse(acc_v_tau < 0.2, cand_v_tau/3, cand_v_tau) )
        
      } 
      
      time_its <- (proc.time() - st)[3] / (i - 1)
      time_used <- round((proc.time() - st)[3]/(60),digits=4)
      time_left <- round(time_its * (reps*thin + burn - i )/(60),digits=4)
      cat("\r", i, " of ", reps*thin + burn,"||| Time left: ",floor(time_left/60),
          " hours",time_left%%60," minutes")# |||| like = ", log_lik[i-1]) 
      flush.console()
      
      
    }
    
    
    
    if(i > burn & i %% thin == 0){
      beta_t[(i - burn)/thin,,] = beta_t_now
      beta_r[(i - burn)/thin,,] = beta_r_now
      Sigma_save[(i - burn)/thin,,] = Sigma_now
      log_sig_pars_R[(i - burn)/thin,] = log_sig_pars_R_now
      U_ij[(i - burn)/thin,,] = U_ij_now
      # mu_save[(i - burn)/thin,,] = mu_now
      like_save[(i - burn)/thin,] = like_now
      sig2_W_save[(i-burn)/thin] = sig2_W_now
      sig2_B_save[(i-burn)/thin,] = sig2_B_now
      
    }
    
  }
  
  plot(Y_use[1,-(1:4)])
  lines(mu_now[1,-(1:4)],col= "red")
  
  K_psi_star = cbind(
    rbind(diag(n_var),matrix(0,n_wave,n_var)),
    rbind(matrix(0,n_var,n_u),K_U)
  )
  
  
  rm(list = setdiff(ls(),c("beta_t","beta_r","Sigma_save","log_sig_pars_R",
                           "like_save","K_psi_star","sig2_W_save")))
  return(list(beta_t = beta_t,beta_r = beta_r,K_psi_star = K_psi_star,
              Sigma_save = Sigma_save,log_sig_pars_R = log_sig_pars_R,
              like_save = like_save,sig2_W_save = sig2_W_save))
  
},mc.preschedule = TRUE,mc.cores = 4)

save.image("four_fam_all.RData")
