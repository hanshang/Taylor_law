##############################
# long-run spatial covariance
##############################

# Version 1: S + T

long_run_spatial = function (dat, w_mat, C0 = 3, H = 3) 
{
  T = dim(dat)[2]
  N = dim(dat)[3]
  center_dat = array(NA, dim = dim(dat))
  for(i in 1:N)
  {
    center_dat[,,i] = t(scale(t(dat[,,i]), center = TRUE, scale = FALSE))
  }
  
  cov_l <- function(porder, band, nval, kern_type) {
    cov_sum = gamma_l(0, nval)
    if (kern_type == "FlatTop") {
      for (ik in 1:(nval - 1)) {
        cov_sum = cov_sum +  ftsa:::FlatTop(ik/band) * abs(ik)^(porder) * 
          (gamma_l(ik, nval) + t(gamma_l(ik, nval)))
      }
    }
    else {
      for (ik in 1:(nval - 1)) {
        cov_sum = cov_sum + ftsa:::kweights(ik/band, kernel = kern_type) * 
          abs(ik)^(porder) * (gamma_l(ik, nval) + t(gamma_l(ik, 
                                                            nval)))
      }
    }
    return(cov_sum)
  }
  gamma_l <- function(lag, T) {
    gamma_lag_sum = matrix(0, nrow = dim(dat)[1], ncol = dim(dat)[1])
    if (lag >= 0) {
      for (ij in 1:(T - lag)) 
      {
        gamma_lag_sum_cross = 0
        for(ik in 1:(N-1))
        {
          for(ih in (ik+1):N)
          {
            gamma_lag_sum_cross = gamma_lag_sum_cross + w_mat[ik,ih] * ifelse(is.na((dat[,ij,ik] - dat[,(ij+lag),ih])^2), 0, (dat[,ij,ik] - dat[,(ij+lag),ih])^2)
          }
        }
        
        gamma_lag_sum_auto = 0
        for(ik in 1:N)
        {
          gamma_lag_sum_auto = gamma_lag_sum_auto + ifelse(is.na(as.matrix(center_dat[,ij,ik]) %*% t(as.matrix(center_dat[, (ij + lag),ik]))), 0, as.matrix(center_dat[,ij,ik]) %*% t(as.matrix(center_dat[, (ij + lag),ik])))
        }
        
        gamma_lag_sum = gamma_lag_sum + diag(gamma_lag_sum_cross/(2*N*(N-1)*sum(w_mat))) + gamma_lag_sum_auto/(N-1)
      }
    }
    else {
      for (ij in 1:(T + lag)) {
        gamma_lag_sum_cross = 0
        for(ik in 1:(N-1))
        {
          for(ih in (ik+1):N)
          {
            gamma_lag_sum_cross = gamma_lag_sum_cross + w_mat[ik,ih] * ifelse(is.na((dat[,ij,ik] - dat[,(ij-lag),ih])^2), 0, (dat[,ij,ik] - dat[,(ij-lag),ih])^2)
          }
        }
        
        gamma_lag_sum_auto = 0
        for(ik in 1:N)
        {
          gamma_lag_sum_auto = gamma_lag_sum_auto + ifelse(is.na(as.matrix(center_dat[,ij,ik]) %*% t(as.matrix(center_dat[, (ij - lag),ik]))), 0, as.matrix(center_dat[,ij,ik]) %*% t(as.matrix(center_dat[, (ij - lag),ik])))
        }
        
        gamma_lag_sum = gamma_lag_sum + gamma_lag_sum_cross/(2*N*(N-1)*sum(w_mat)) + gamma_lag_sum_auto/(N-1)
      }
    }
    return(gamma_lag_sum/T)
  }
  gamma_hat = list()
  for (ik in 1:T) {
    gamma_hat[[ik]] = gamma_l(lag = ik - 1, T = T)
  }
  rho_hat <- function(gamma, ik) {
    a = sum(diag(gamma[[1]]))
    denom = a * a
    b = gamma[[ik + 1]] * gamma[[ik + 1]]
    num = sum(b)
    c = num/denom
    calc = sqrt(c)
    return(calc)
  }
  rho = vector(, T - 1)
  for (ik in 1:(T - 1)) {
    rho[ik] = rho_hat(gamma_hat, ik)
  }
  h_hat <- function(rho, N, H) {
    l <- -1
    c <- 1
    end <- N - H
    const <- C0 * sqrt(log10(N)/N)
    while (c <= end) {
      if (rho[c] <= const) {
        c2 <- 1
        while (c2 <= H - 1) {
          if (rho[c + c2] <= const) {
            c2 <- c2 + 1
            l <- c - 1
          }
          else {
            c2 <- H
            l <- -1
            c <- c + 1
          }
        }
        if (l != -1) {
          c <- end + 1
        }
      }
      else {
        c <- c + 1
      }
    }
    return(l)
  }
  h_adaptive_ini <- h_hat(rho = rho, N = T, H = H)
  h_opt_fun <- function(band_ini, T, kern_type, kern_type_ini) {
    kern_name_kweights_ini = switch(kern_type_ini, BT = "Bartlett", 
                                    PR = "Parzen", TH = "Tukey-Hanning", 
                                    QS = "Quadratic Spectral", FT = "FlatTop")
    kern_name_kweights = switch(kern_type, BT = "Bartlett", 
                                PR = "Parzen", TH = "Tukey-Hanning", 
                                QS = "Quadratic Spectral")
    w_weights = switch(kern_type, BT = 1, PR = 6, TH = pi^2/4, 
                       QS = 18 * pi * pi/125)
    q = switch(kern_type, BT = 1, PR = 2, TH = 2, QS = 2)
    C_0 = cov_l(porder = 0, band = band_ini, nval = T, kern_type = kern_name_kweights_ini)
    C_2 = w_weights * cov_l(porder = q, band = band_ini, 
                            nval = T, kern_type = kern_name_kweights_ini)
    first_part = (2 * q * sum(C_2^2))^(1/(1 + 2 * q))
    kernel_square_int = switch(kern_type, BT = 2/3, PR = 0.539285, 
                               TH = 3/4, QS = 1, FT = 4/3)
    second_part = ((sum(C_0^2) + sum(diag(C_0))^2) * kernel_square_int)^(-1/(1 + 
                                                                               2 * q))
    c_0 = first_part * second_part
    hat_h_opt = c_0 * (T^(1/(1 + 2 * q)))
    C_0_est = cov_l(porder = 0, band = hat_h_opt, nval = T, 
                    kern_type = kern_name_kweights)
    return(list(hat_h_opt = hat_h_opt, C_0_est = C_0_est))
  }
  h_opt_BT = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "BT", 
                       kern_type_ini = "BT")
  h_opt_PR = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "PR", 
                       kern_type_ini = "PR")
  h_opt_TH = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "TH", 
                       kern_type_ini = "TH")
  h_opt_QS = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "QS", 
                       kern_type_ini = "QS")
  h_opt_FT_BT = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "BT", 
                          kern_type_ini = "FT")
  h_opt_FT_PR = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "PR", 
                          kern_type_ini = "FT")
  h_opt_FT_TH = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "TH", 
                          kern_type_ini = "FT")
  h_opt_FT_QS = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "QS", 
                          kern_type_ini = "FT")
  h_adaptive_opt_FT_BT = h_opt_fun(band_ini = h_adaptive_ini, 
                                   T = T, kern_type = "BT", kern_type_ini = "FT")
  h_adaptive_opt_FT_PR = h_opt_fun(band_ini = h_adaptive_ini, 
                                   T = T, kern_type = "PR", kern_type_ini = "FT")
  h_adaptive_opt_FT_TH = h_opt_fun(band_ini = h_adaptive_ini, 
                                   T = T, kern_type = "TH", kern_type_ini = "FT")
  h_adaptive_opt_FT_QS = h_opt_fun(band_ini = h_adaptive_ini, 
                                   T = T, kern_type = "QS", kern_type_ini = "FT")
  return(list(BT_finite_fix_C0 = h_opt_BT$C_0_est, PR_finite_fix_C0 = h_opt_PR$C_0_est, 
              TH_finite_fix_C0 = h_opt_TH$C_0_est, QS_finite_fix_C0 = h_opt_QS$C_0_est, 
              BT_FT_fix_C0 = h_opt_FT_BT$C_0_est, PR_FT_fix_C0 = h_opt_FT_PR$C_0_est, 
              TH_FT_fix_C0 = h_opt_FT_TH$C_0_est, QS_FT_fix_C0 = h_opt_FT_QS$C_0_est, 
              BT_FT_adaptive_C0 = h_adaptive_opt_FT_BT$C_0_est, PR_FT_adaptive_C0 = h_adaptive_opt_FT_PR$C_0_est, 
              TH_FT_adaptive_C0 = h_adaptive_opt_FT_TH$C_0_est, QS_FT_adaptive_C0 = h_adaptive_opt_FT_QS$C_0_est))
}

female_long_run_spatial_var = long_run_spatial(dat = female_data[,,2:48], w_mat = weights_matrix_2)
male_long_run_spatial_var = long_run_spatial(dat = male_data[,,2:48], w_mat = weights_matrix_2)
total_long_run_spatial_var = long_run_spatial(dat = total_data[,,2:48], w_mat = weights_matrix_2)

# Version 2: c*S + (1-c)*T

long_run_spatial_combined = function (dat, w_mat, c_weight = 0, C0 = 3, H = 3) 
{
  T = dim(dat)[2]
  N = dim(dat)[3]
  center_dat = array(NA, dim = dim(dat))
  for(i in 1:N)
  {
    center_dat[,,i] = t(scale(t(dat[,,i]), center = TRUE, scale = FALSE))
  }
  
  cov_l <- function(porder, band, nval, kern_type) {
    cov_sum = gamma_l(0, nval)
    if (kern_type == "FlatTop") {
      for (ik in 1:(nval - 1)) {
        cov_sum = cov_sum +  ftsa:::FlatTop(ik/band) * abs(ik)^(porder) * 
          (gamma_l(ik, nval) + t(gamma_l(ik, nval)))
      }
    }
    else {
      for (ik in 1:(nval - 1)) {
        cov_sum = cov_sum + ftsa:::kweights(ik/band, kernel = kern_type) * 
          abs(ik)^(porder) * (gamma_l(ik, nval) + t(gamma_l(ik, 
                                                            nval)))
      }
    }
    return(cov_sum)
  }
  gamma_l <- function(lag, T) {
    gamma_lag_sum = matrix(0, nrow = dim(dat)[1], ncol = dim(dat)[1])
    if (lag >= 0) {
      for (ij in 1:(T - lag)) 
      {
        gamma_lag_sum_cross = 0
        for(ik in 1:(N-1))
        {
          for(ih in (ik+1):N)
          {
            gamma_lag_sum_cross = gamma_lag_sum_cross + w_mat[ik,ih] * ifelse(is.na((dat[,ij,ik] - dat[,(ij+lag),ih])^2), 0, (dat[,ij,ik] - dat[,(ij+lag),ih])^2)
          }
        }
        
        gamma_lag_sum_auto = 0
        for(ik in 1:N)
        {
          gamma_lag_sum_auto = gamma_lag_sum_auto + ifelse(is.na(as.matrix(center_dat[,ij,ik]) %*% t(as.matrix(center_dat[, (ij + lag),ik]))), 0, as.matrix(center_dat[,ij,ik]) %*% t(as.matrix(center_dat[, (ij + lag),ik])))
        }
        
        gamma_lag_sum = gamma_lag_sum + (c_weight*diag(gamma_lag_sum_cross/(2*N*(N-1)*sum(w_mat))) + (1-c_weight)*gamma_lag_sum_auto/(N-1)) 
      }
    }
    else {
      for (ij in 1:(T + lag)) {
        gamma_lag_sum_cross = 0
        for(ik in 1:(N-1))
        {
          for(ih in (ik+1):N)
          {
            gamma_lag_sum_cross = gamma_lag_sum_cross + w_mat[ik,ih] * ifelse(is.na((dat[,ij,ik] - dat[,(ij-lag),ih])^2), 0, (dat[,ij,ik] - dat[,(ij-lag),ih])^2)
          }
        }
        
        gamma_lag_sum_auto = 0
        for(ik in 1:N)
        {
          gamma_lag_sum_auto = gamma_lag_sum_auto + ifelse(is.na(as.matrix(center_dat[,ij,ik]) %*% t(as.matrix(center_dat[, (ij - lag),ik]))), 0, as.matrix(center_dat[,ij,ik]) %*% t(as.matrix(center_dat[, (ij - lag),ik])))
        }
        
        gamma_lag_sum = gamma_lag_sum + (c_weight*diag(gamma_lag_sum_cross/(2*N*(N-1)*sum(w_mat))) + (1-c_weight)*gamma_lag_sum_auto/(N-1)) 
      }
    }
    return(gamma_lag_sum/T)
  }
  gamma_hat = list()
  for (ik in 1:T) {
    gamma_hat[[ik]] = gamma_l(lag = ik - 1, T = T)
  }
  rho_hat <- function(gamma, ik) {
    a = sum(diag(gamma[[1]]))
    denom = a * a
    b = gamma[[ik + 1]] * gamma[[ik + 1]]
    num = sum(b)
    c = num/denom
    calc = sqrt(c)
    return(calc)
  }
  rho = vector(, T - 1)
  for (ik in 1:(T - 1)) {
    rho[ik] = rho_hat(gamma_hat, ik)
  }
  h_hat <- function(rho, N, H) {
    l <- -1
    c <- 1
    end <- N - H
    const <- C0 * sqrt(log10(N)/N)
    while (c <= end) {
      if (rho[c] <= const) {
        c2 <- 1
        while (c2 <= H - 1) {
          if (rho[c + c2] <= const) {
            c2 <- c2 + 1
            l <- c - 1
          }
          else {
            c2 <- H
            l <- -1
            c <- c + 1
          }
        }
        if (l != -1) {
          c <- end + 1
        }
      }
      else {
        c <- c + 1
      }
    }
    return(l)
  }
  h_adaptive_ini <- h_hat(rho = rho, N = T, H = H)
  h_opt_fun <- function(band_ini, T, kern_type, kern_type_ini) {
    kern_name_kweights_ini = switch(kern_type_ini, BT = "Bartlett", 
                                    PR = "Parzen", TH = "Tukey-Hanning", 
                                    QS = "Quadratic Spectral", FT = "FlatTop")
    kern_name_kweights = switch(kern_type, BT = "Bartlett", 
                                PR = "Parzen", TH = "Tukey-Hanning", 
                                QS = "Quadratic Spectral")
    w_weights = switch(kern_type, BT = 1, PR = 6, TH = pi^2/4, 
                       QS = 18 * pi * pi/125)
    q = switch(kern_type, BT = 1, PR = 2, TH = 2, QS = 2)
    C_0 = cov_l(porder = 0, band = band_ini, nval = T, kern_type = kern_name_kweights_ini)
    C_2 = w_weights * cov_l(porder = q, band = band_ini, 
                            nval = T, kern_type = kern_name_kweights_ini)
    first_part = (2 * q * sum(C_2^2))^(1/(1 + 2 * q))
    kernel_square_int = switch(kern_type, BT = 2/3, PR = 0.539285, 
                               TH = 3/4, QS = 1, FT = 4/3)
    second_part = ((sum(C_0^2) + sum(diag(C_0))^2) * kernel_square_int)^(-1/(1 + 
                                                                               2 * q))
    c_0 = first_part * second_part
    hat_h_opt = c_0 * (T^(1/(1 + 2 * q)))
    C_0_est = cov_l(porder = 0, band = hat_h_opt, nval = T, 
                    kern_type = kern_name_kweights)
    return(list(hat_h_opt = hat_h_opt, C_0_est = C_0_est))
  }
  h_opt_BT = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "BT", 
                       kern_type_ini = "BT")
  h_opt_PR = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "PR", 
                       kern_type_ini = "PR")
  h_opt_TH = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "TH", 
                       kern_type_ini = "TH")
  h_opt_QS = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "QS", 
                       kern_type_ini = "QS")
  h_opt_FT_BT = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "BT", 
                          kern_type_ini = "FT")
  h_opt_FT_PR = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "PR", 
                          kern_type_ini = "FT")
  h_opt_FT_TH = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "TH", 
                          kern_type_ini = "FT")
  h_opt_FT_QS = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "QS", 
                          kern_type_ini = "FT")
  h_adaptive_opt_FT_BT = h_opt_fun(band_ini = h_adaptive_ini, 
                                   T = T, kern_type = "BT", kern_type_ini = "FT")
  h_adaptive_opt_FT_PR = h_opt_fun(band_ini = h_adaptive_ini, 
                                   T = T, kern_type = "PR", kern_type_ini = "FT")
  h_adaptive_opt_FT_TH = h_opt_fun(band_ini = h_adaptive_ini, 
                                   T = T, kern_type = "TH", kern_type_ini = "FT")
  h_adaptive_opt_FT_QS = h_opt_fun(band_ini = h_adaptive_ini, 
                                   T = T, kern_type = "QS", kern_type_ini = "FT")
  return(list(BT_finite_fix_C0 = h_opt_BT$C_0_est, PR_finite_fix_C0 = h_opt_PR$C_0_est, 
              TH_finite_fix_C0 = h_opt_TH$C_0_est, QS_finite_fix_C0 = h_opt_QS$C_0_est, 
              BT_FT_fix_C0 = h_opt_FT_BT$C_0_est, PR_FT_fix_C0 = h_opt_FT_PR$C_0_est, 
              TH_FT_fix_C0 = h_opt_FT_TH$C_0_est, QS_FT_fix_C0 = h_opt_FT_QS$C_0_est, 
              BT_FT_adaptive_C0 = h_adaptive_opt_FT_BT$C_0_est, PR_FT_adaptive_C0 = h_adaptive_opt_FT_PR$C_0_est, 
              TH_FT_adaptive_C0 = h_adaptive_opt_FT_TH$C_0_est, QS_FT_adaptive_C0 = h_adaptive_opt_FT_QS$C_0_est))
}


# Version 3: c(S+T)+(1-c)(S*T)

long_run_spatial_combined_complex = function (dat, w_mat, c_weight = 0, C0 = 3, H = 3) 
{
  T = dim(dat)[2]
  N = dim(dat)[3]
  center_dat = array(NA, dim = dim(dat))
  for(i in 1:N)
  {
    center_dat[,,i] = t(scale(t(dat[,,i]), center = TRUE, scale = FALSE))
  }
  
  cov_l <- function(porder, band, nval, kern_type) {
    cov_sum = gamma_l(0, nval)
    if (kern_type == "FlatTop") {
      for (ik in 1:(nval - 1)) {
        cov_sum = cov_sum +  ftsa:::FlatTop(ik/band) * abs(ik)^(porder) * 
          (gamma_l(ik, nval) + t(gamma_l(ik, nval)))
      }
    }
    else {
      for (ik in 1:(nval - 1)) {
        cov_sum = cov_sum + ftsa:::kweights(ik/band, kernel = kern_type) * 
          abs(ik)^(porder) * (gamma_l(ik, nval) + t(gamma_l(ik, 
                                                            nval)))
      }
    }
    return(cov_sum)
  }
  gamma_l <- function(lag, T) {
    gamma_lag_sum = matrix(0, nrow = dim(dat)[1], ncol = dim(dat)[1])
    if (lag >= 0) {
      for (ij in 1:(T - lag)) 
      {
        gamma_lag_sum_cross = 0
        for(ik in 1:(N-1))
        {
          for(ih in (ik+1):N)
          {
            gamma_lag_sum_cross = gamma_lag_sum_cross + w_mat[ik,ih] * ifelse(is.na((dat[,ij,ik] - dat[,(ij+lag),ih])^2), 0, (dat[,ij,ik] - dat[,(ij+lag),ih])^2)
          }
        }
        
        gamma_lag_sum_auto = 0
        for(ik in 1:N)
        {
          gamma_lag_sum_auto = gamma_lag_sum_auto + ifelse(is.na(as.matrix(center_dat[,ij,ik]) %*% t(as.matrix(center_dat[, (ij + lag),ik]))), 0, as.matrix(center_dat[,ij,ik]) %*% t(as.matrix(center_dat[, (ij + lag),ik])))
        }
        
        gamma_lag_sum = gamma_lag_sum + c_weight*(diag(gamma_lag_sum_cross/(2*N*(N-1)*sum(w_mat))) + gamma_lag_sum_auto/(N-1)) + (1-c_weight)*(diag(gamma_lag_sum_cross/(2*N*(N-1)*sum(w_mat))) * gamma_lag_sum_auto/(N-1))
      }
    }
    else {
      for (ij in 1:(T + lag)) {
        gamma_lag_sum_cross = 0
        for(ik in 1:(N-1))
        {
          for(ih in (ik+1):N)
          {
            gamma_lag_sum_cross = gamma_lag_sum_cross + w_mat[ik,ih] * ifelse(is.na((dat[,ij,ik] - dat[,(ij-lag),ih])^2), 0, (dat[,ij,ik] - dat[,(ij-lag),ih])^2)
          }
        }
        
        gamma_lag_sum_auto = 0
        for(ik in 1:N)
        {
          gamma_lag_sum_auto = gamma_lag_sum_auto + ifelse(is.na(as.matrix(center_dat[,ij,ik]) %*% t(as.matrix(center_dat[, (ij - lag),ik]))), 0, as.matrix(center_dat[,ij,ik]) %*% t(as.matrix(center_dat[, (ij - lag),ik])))
        }
        
        gamma_lag_sum = gamma_lag_sum + c_weight*(diag(gamma_lag_sum_cross/(2*N*(N-1)*sum(w_mat))) + gamma_lag_sum_auto/(N-1)) + (1-c_weight)*(diag(gamma_lag_sum_cross/(2*N*(N-1)*sum(w_mat))) * gamma_lag_sum_auto/(N-1))
      }
    }
    return(gamma_lag_sum/T)
  }
  gamma_hat = list()
  for (ik in 1:T) {
    gamma_hat[[ik]] = gamma_l(lag = ik - 1, T = T)
  }
  rho_hat <- function(gamma, ik) {
    a = sum(diag(gamma[[1]]))
    denom = a * a
    b = gamma[[ik + 1]] * gamma[[ik + 1]]
    num = sum(b)
    c = num/denom
    calc = sqrt(c)
    return(calc)
  }
  rho = vector(, T - 1)
  for (ik in 1:(T - 1)) {
    rho[ik] = rho_hat(gamma_hat, ik)
  }
  h_hat <- function(rho, N, H) {
    l <- -1
    c <- 1
    end <- N - H
    const <- C0 * sqrt(log10(N)/N)
    while (c <= end) {
      if (rho[c] <= const) {
        c2 <- 1
        while (c2 <= H - 1) {
          if (rho[c + c2] <= const) {
            c2 <- c2 + 1
            l <- c - 1
          }
          else {
            c2 <- H
            l <- -1
            c <- c + 1
          }
        }
        if (l != -1) {
          c <- end + 1
        }
      }
      else {
        c <- c + 1
      }
    }
    return(l)
  }
  h_adaptive_ini <- h_hat(rho = rho, N = T, H = H)
  h_opt_fun <- function(band_ini, T, kern_type, kern_type_ini) {
    kern_name_kweights_ini = switch(kern_type_ini, BT = "Bartlett", 
                                    PR = "Parzen", TH = "Tukey-Hanning", 
                                    QS = "Quadratic Spectral", FT = "FlatTop")
    kern_name_kweights = switch(kern_type, BT = "Bartlett", 
                                PR = "Parzen", TH = "Tukey-Hanning", 
                                QS = "Quadratic Spectral")
    w_weights = switch(kern_type, BT = 1, PR = 6, TH = pi^2/4, 
                       QS = 18 * pi * pi/125)
    q = switch(kern_type, BT = 1, PR = 2, TH = 2, QS = 2)
    C_0 = cov_l(porder = 0, band = band_ini, nval = T, kern_type = kern_name_kweights_ini)
    C_2 = w_weights * cov_l(porder = q, band = band_ini, 
                            nval = T, kern_type = kern_name_kweights_ini)
    first_part = (2 * q * sum(C_2^2))^(1/(1 + 2 * q))
    kernel_square_int = switch(kern_type, BT = 2/3, PR = 0.539285, 
                               TH = 3/4, QS = 1, FT = 4/3)
    second_part = ((sum(C_0^2) + sum(diag(C_0))^2) * kernel_square_int)^(-1/(1 + 
                                                                               2 * q))
    c_0 = first_part * second_part
    hat_h_opt = c_0 * (T^(1/(1 + 2 * q)))
    C_0_est = cov_l(porder = 0, band = hat_h_opt, nval = T, 
                    kern_type = kern_name_kweights)
    return(list(hat_h_opt = hat_h_opt, C_0_est = C_0_est))
  }
  h_opt_BT = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "BT", 
                       kern_type_ini = "BT")
  h_opt_PR = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "PR", 
                       kern_type_ini = "PR")
  h_opt_TH = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "TH", 
                       kern_type_ini = "TH")
  h_opt_QS = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "QS", 
                       kern_type_ini = "QS")
  h_opt_FT_BT = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "BT", 
                          kern_type_ini = "FT")
  h_opt_FT_PR = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "PR", 
                          kern_type_ini = "FT")
  h_opt_FT_TH = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "TH", 
                          kern_type_ini = "FT")
  h_opt_FT_QS = h_opt_fun(band_ini = T^(0.2), T = T, kern_type = "QS", 
                          kern_type_ini = "FT")
  h_adaptive_opt_FT_BT = h_opt_fun(band_ini = h_adaptive_ini, 
                                   T = T, kern_type = "BT", kern_type_ini = "FT")
  h_adaptive_opt_FT_PR = h_opt_fun(band_ini = h_adaptive_ini, 
                                   T = T, kern_type = "PR", kern_type_ini = "FT")
  h_adaptive_opt_FT_TH = h_opt_fun(band_ini = h_adaptive_ini, 
                                   T = T, kern_type = "TH", kern_type_ini = "FT")
  h_adaptive_opt_FT_QS = h_opt_fun(band_ini = h_adaptive_ini, 
                                   T = T, kern_type = "QS", kern_type_ini = "FT")
  return(list(BT_finite_fix_C0 = h_opt_BT$C_0_est, PR_finite_fix_C0 = h_opt_PR$C_0_est, 
              TH_finite_fix_C0 = h_opt_TH$C_0_est, QS_finite_fix_C0 = h_opt_QS$C_0_est, 
              BT_FT_fix_C0 = h_opt_FT_BT$C_0_est, PR_FT_fix_C0 = h_opt_FT_PR$C_0_est, 
              TH_FT_fix_C0 = h_opt_FT_TH$C_0_est, QS_FT_fix_C0 = h_opt_FT_QS$C_0_est, 
              BT_FT_adaptive_C0 = h_adaptive_opt_FT_BT$C_0_est, PR_FT_adaptive_C0 = h_adaptive_opt_FT_PR$C_0_est, 
              TH_FT_adaptive_C0 = h_adaptive_opt_FT_TH$C_0_est, QS_FT_adaptive_C0 = h_adaptive_opt_FT_QS$C_0_est))
}

#############################################
# Compare long-run spatial variance settings
#############################################

long_run_spatial_parallel <- function(data, w_mat, c_all, ik, method = c("Simple", "Complex"))
{
  if(method == "Simple")
  {
    long_run_spatial_variance = long_run_spatial_combined(dat = data, w_mat = w_mat, c_weight = c_all[ik])
  }
  
  if(method == "Complex")
  {
    long_run_spatial_variance = long_run_spatial_combined_complex(dat = data, w_mat = w_mat, c_weight = c_all[ik])
  }
  
  return(list(long_run_spatial_variance = long_run_spatial_variance))
}

c_grid = seq(0, 1, by = 0.05)

# female

library(doParallel)
cl <- makeCluster(10, setup_strategy = "sequential")
registerDoParallel(cl)
female_long_run_spatial_combined = foreach(ij = 1:21) %dopar% long_run_spatial_parallel(data = female_data[,,2:48], w_mat = weights_matrix_2, method = "Simple", c_all = c_grid, ik = ij)
stopCluster(cl)
rm(cl)
gc()

cl <- makeCluster(10, setup_strategy = "sequential")
registerDoParallel(cl)
female_long_run_spatial_combined_complex = foreach(ij = 1:21) %dopar% long_run_spatial_parallel(data = female_data[,,2:48], w_mat = weights_matrix_2, method = "Complex", c_all = c_grid, ik = ij)
stopCluster(cl)
rm(cl)
gc()

save(female_long_run_spatial_combined, file = "female_long_run_spatial_combined.RData")
save(female_long_run_spatial_combined_complex, file = "female_long_run_spatial_combined_complex.RData")

# male
cl <- makeCluster(10, setup_strategy = "sequential")
registerDoParallel(cl)
male_long_run_spatial_combined = foreach(ij = 1:21) %dopar% long_run_spatial_parallel(data = male_data[,,2:48], w_mat = weights_matrix_2, method = "Simple", c_all = c_grid, ik = ij)
stopCluster(cl)
rm(cl)
gc()

cl <- makeCluster(10, setup_strategy = "sequential")
registerDoParallel(cl)
male_long_run_spatial_combined_complex = foreach(ij = 1:21) %dopar% long_run_spatial_parallel(data = male_data[,,2:48], w_mat = weights_matrix_2, method = "Complex", c_all = c_grid, ik = ij)
stopCluster(cl)
rm(cl)
gc()

save(male_long_run_spatial_combined, file = "male_long_run_spatial_combined.RData")
save(male_long_run_spatial_combined_complex, file = "male_long_run_spatial_combined_complex.RData")

# total
cl <- makeCluster(10, setup_strategy = "sequential")
registerDoParallel(cl)
total_long_run_spatial_combined = foreach(ij = 1:21) %dopar% long_run_spatial_parallel(data = total_data[,,2:48], w_mat = weights_matrix_2, method = "Simple", c_all = c_grid, ik = ij)
stopCluster(cl)
rm(cl)
gc()

cl <- makeCluster(10, setup_strategy = "sequential")
registerDoParallel(cl)
total_long_run_spatial_combined_complex = foreach(ij = 1:21) %dopar% long_run_spatial_parallel(data = total_data[,,2:48], w_mat = weights_matrix_2, method = "Complex", c_all = c_grid, ik = ij)
stopCluster(cl)
rm(cl)
gc()

save(total_long_run_spatial_combined, file = "total_long_run_spatial_combined.RData")
save(total_long_run_spatial_combined_complex, file = "total_long_run_spatial_combined_complex.RData")


# Long-run spatial variance version 1

# female
y = log(diag(female_long_run_spatial_var$BT_FT_fix_C0), base = 10)
x = log(rowMeans(Japan$rate$female), base = 10)

OLS_long_run_female = lm(y ~ x)
OLS_long_run_female_quadratic = lm(y ~ x + I(x^2))
OLS_long_run_female_cubic = lm(y ~ x + I(x^2) + I(x^3))

hampel_long_run_female = rlm(y ~ x, psi = psi.hampel, maxit = 200)
hampel_long_run_female_quadratic = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
hampel_long_run_female_cubic = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)

bisquare_long_run_female = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
bisquare_long_run_female_quadratic = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
bisquare_long_run_female_cubic = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)

# male
y = log(diag(male_long_run_spatial_var$BT_FT_fix_C0), base = 10)
x = log(rowMeans(Japan$rate$male), base = 10)

OLS_long_run_male = lm(y ~ x)
OLS_long_run_male_quadratic = lm(y ~ x + I(x^2))
OLS_long_run_male_cubic = lm(y ~ x + I(x^2) + I(x^3))

hampel_long_run_male = rlm(y ~ x, psi = psi.hampel, maxit = 200)
hampel_long_run_male_quadratic = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
hampel_long_run_male_cubic = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)

bisquare_long_run_male = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
bisquare_long_run_male_quadratic = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
bisquare_long_run_male_cubic = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)

# total
y = log(diag(total_long_run_spatial_var$BT_FT_fix_C0), base = 10)
x = log(rowMeans(Japan$rate$total), base = 10)

OLS_long_run_total = lm(y ~ x)
OLS_long_run_total_quadratic = lm(y ~ x + I(x^2))
OLS_long_run_total_cubic = lm(y ~ x + I(x^2) + I(x^3))

hampel_long_run_total = rlm(y ~ x, psi = psi.hampel, maxit = 200)
hampel_long_run_total_quadratic = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
hampel_long_run_total_cubic = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)

bisquare_long_run_total = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
bisquare_long_run_total_quadratic = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
bisquare_long_run_total_cubic = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)

# Long-run spatial variance version 2 & 3; search for the optimal c_weight

library(asbio)
library(MuMIn)
library(MASS)

# female

# Linear TL

OLS_list_female_combined = hampel_list_female_combined = bisquare_list_female_combined = list()
OLS_list_female_combined_complex = hampel_list_female_combined_complex = bisquare_list_female_combined_complex = list()

for(i in 1:21)
{
  y = log(diag(female_long_run_spatial_combined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_female_combined[[i]] = lm(y ~ x)
  
  
  hampel_list_female_combined[[i]] = rlm(y ~ x, psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_female_combined[[i]] = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
}

for(i in 1:21)
{
  y = log(diag(female_long_run_spatial_combined_complex[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_female_combined_complex[[i]] = lm(y ~ x)
  
  
  hampel_list_female_combined_complex[[i]] = rlm(y ~ x, psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_female_combined_complex[[i]] = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
}

OLS_female_fitting_combined_criteria = OLS_female_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_female_fitting_combined_criteria = hampel_female_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)
bisquare_female_fitting_combined_criteria = bisquare_female_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)


for(i in 1:21)
{
  # OLS
  OLS_female_fitting_combined_criteria[i, "R_squared"] = summary(OLS_list_female_combined[[i]])$adj.r.squared
  OLS_female_fitting_combined_criteria[i, "AIC"] = AIC(OLS_list_female_combined[[i]])
  OLS_female_fitting_combined_criteria[i, "AICc"] = AICc(OLS_list_female_combined[[i]])
  OLS_female_fitting_combined_criteria[i, "R_squared_pred"] = press(OLS_list_female_combined[[i]], as.R2 = TRUE)
  
  OLS_female_fitting_combined_complex_criteria[i, "R_squared"] = summary(OLS_list_female_combined_complex[[i]])$adj.r.squared
  OLS_female_fitting_combined_complex_criteria[i, "AIC"] = AIC(OLS_list_female_combined_complex[[i]])
  OLS_female_fitting_combined_complex_criteria[i, "AICc"] = AICc(OLS_list_female_combined_complex[[i]])
  OLS_female_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(OLS_list_female_combined_complex[[i]], as.R2 = TRUE)
  
  # hample
  hampel_female_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((hampel_list_female_combined[[i]]$residuals)^2)/sum((hampel_list_female_combined[[i]]$model[,1] - mean(hampel_list_female_combined[[i]]$model[,1]))^2))*(100)/(97)
  
  hampel_female_fitting_combined_criteria[i, "AIC"] = AIC(hampel_list_female_combined[[i]])
  hampel_female_fitting_combined_criteria[i, "AICc"] = AICc(hampel_list_female_combined[[i]])
  hampel_female_fitting_combined_criteria[i, "R_squared_pred"] = press(hampel_list_female_combined[[i]], as.R2 = TRUE)
  
  hampel_female_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((hampel_list_female_combined_complex[[i]]$residuals)^2)/sum((hampel_list_female_combined_complex[[i]]$model[,1] - mean(hampel_list_female_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  hampel_female_fitting_combined_complex_criteria[i, "AIC"] = AIC(hampel_list_female_combined_complex[[i]])
  hampel_female_fitting_combined_complex_criteria[i, "AICc"] = AICc(hampel_list_female_combined_complex[[i]])
  hampel_female_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(hampel_list_female_combined_complex[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_female_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_female_combined[[i]]$residuals)^2)/sum((bisquare_list_female_combined[[i]]$model[,1] - mean(bisquare_list_female_combined[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_female_fitting_combined_criteria[i, "AIC"] = AIC(bisquare_list_female_combined[[i]])
  bisquare_female_fitting_combined_criteria[i, "AICc"] = AICc(bisquare_list_female_combined[[i]])
  bisquare_female_fitting_combined_criteria[i, "R_squared_pred"] = press(bisquare_list_female_combined[[i]], as.R2 = TRUE)
  
  bisquare_female_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_female_combined_complex[[i]]$residuals)^2)/sum((bisquare_list_female_combined_complex[[i]]$model[,1] - mean(bisquare_list_female_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_female_fitting_combined_complex_criteria[i, "AIC"] = AIC(bisquare_list_female_combined_complex[[i]])
  bisquare_female_fitting_combined_complex_criteria[i, "AICc"] = AICc(bisquare_list_female_combined_complex[[i]])
  bisquare_female_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(bisquare_list_female_combined_complex[[i]], as.R2 = TRUE)
}

rownames(OLS_female_fitting_combined_criteria) = rownames(hampel_female_fitting_combined_criteria) = rownames(bisquare_female_fitting_combined_criteria) = rownames(OLS_female_fitting_combined_complex_criteria) = rownames(hampel_female_fitting_combined_complex_criteria) = rownames(bisquare_female_fitting_combined_complex_criteria) = c_grid

# Quadratic TL
OLS_list_female_quadratic_combined = hampel_list_female_quadratic_combined = bisquare_list_female_quadratic_combined = list()
OLS_list_female_quadratic_combined_complex = hampel_list_female_quadratic_combined_complex = bisquare_list_female_quadratic_combined_complex = list()

for(i in 1:21)
{
  y = log(diag(female_long_run_spatial_combined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_female_quadratic_combined[[i]] = lm(y ~ x + I(x^2))
  
  
  hampel_list_female_quadratic_combined[[i]] = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_female_quadratic_combined[[i]] = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
}

for(i in 1:21)
{
  y = log(diag(female_long_run_spatial_combined_complex[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_female_quadratic_combined_complex[[i]] = lm(y ~ x + I(x^2))
  
  
  hampel_list_female_quadratic_combined_complex[[i]] = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_female_quadratic_combined_complex[[i]] = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
}


OLS_female_quadratic_fitting_combined_criteria = OLS_female_quadratic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_female_quadratic_fitting_combined_criteria = hampel_female_quadratic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)
bisquare_female_quadratic_fitting_combined_criteria = bisquare_female_quadratic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)

for(i in 1:21)
{
  # OLS
  OLS_female_quadratic_fitting_combined_criteria[i, "R_squared"] = summary(OLS_list_female_quadratic_combined[[i]])$adj.r.squared
  OLS_female_quadratic_fitting_combined_criteria[i, "AIC"] = AIC(OLS_list_female_quadratic_combined[[i]])
  OLS_female_quadratic_fitting_combined_criteria[i, "AICc"] = AICc(OLS_list_female_quadratic_combined[[i]])
  OLS_female_quadratic_fitting_combined_criteria[i, "R_squared_pred"] = press(OLS_list_female_quadratic_combined[[i]], as.R2 = TRUE)
  
  OLS_female_quadratic_fitting_combined_complex_criteria[i, "R_squared"] = summary(OLS_list_female_quadratic_combined_complex[[i]])$adj.r.squared
  OLS_female_quadratic_fitting_combined_complex_criteria[i, "AIC"] = AIC(OLS_list_female_quadratic_combined_complex[[i]])
  OLS_female_quadratic_fitting_combined_complex_criteria[i, "AICc"] = AICc(OLS_list_female_quadratic_combined_complex[[i]])
  OLS_female_quadratic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(OLS_list_female_quadratic_combined_complex[[i]], as.R2 = TRUE)
  
  # hample
  hampel_female_quadratic_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((hampel_list_female_quadratic_combined[[i]]$residuals)^2)/sum((hampel_list_female_quadratic_combined[[i]]$model[,1] - mean(hampel_list_female_quadratic_combined[[i]]$model[,1]))^2))*(100)/(97)
  
  hampel_female_quadratic_fitting_combined_criteria[i, "AIC"] = AIC(hampel_list_female_quadratic_combined[[i]])
  hampel_female_quadratic_fitting_combined_criteria[i, "AICc"] = AICc(hampel_list_female_quadratic_combined[[i]])
  hampel_female_quadratic_fitting_combined_criteria[i, "R_squared_pred"] = press(hampel_list_female_quadratic_combined[[i]], as.R2 = TRUE)
  
  hampel_female_quadratic_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((hampel_list_female_quadratic_combined_complex[[i]]$residuals)^2)/sum((hampel_list_female_quadratic_combined_complex[[i]]$model[,1] - mean(hampel_list_female_quadratic_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  hampel_female_quadratic_fitting_combined_complex_criteria[i, "AIC"] = AIC(hampel_list_female_quadratic_combined_complex[[i]])
  hampel_female_quadratic_fitting_combined_complex_criteria[i, "AICc"] = AICc(hampel_list_female_quadratic_combined_complex[[i]])
  hampel_female_quadratic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(hampel_list_female_quadratic_combined_complex[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_female_quadratic_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_female_quadratic_combined[[i]]$residuals)^2)/sum((bisquare_list_female_quadratic_combined[[i]]$model[,1] - mean(bisquare_list_female_quadratic_combined[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_female_quadratic_fitting_combined_criteria[i, "AIC"] = AIC(bisquare_list_female_quadratic_combined[[i]])
  bisquare_female_quadratic_fitting_combined_criteria[i, "AICc"] = AICc(bisquare_list_female_quadratic_combined[[i]])
  bisquare_female_quadratic_fitting_combined_criteria[i, "R_squared_pred"] = press(bisquare_list_female_quadratic_combined[[i]], as.R2 = TRUE)
  
  bisquare_female_quadratic_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_female_quadratic_combined_complex[[i]]$residuals)^2)/sum((bisquare_list_female_quadratic_combined_complex[[i]]$model[,1] - mean(bisquare_list_female_quadratic_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_female_quadratic_fitting_combined_complex_criteria[i, "AIC"] = AIC(bisquare_list_female_quadratic_combined_complex[[i]])
  bisquare_female_quadratic_fitting_combined_complex_criteria[i, "AICc"] = AICc(bisquare_list_female_quadratic_combined_complex[[i]])
  bisquare_female_quadratic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(bisquare_list_female_quadratic_combined_complex[[i]], as.R2 = TRUE)
}

rownames(OLS_female_quadratic_fitting_combined_criteria) = rownames(hampel_female_quadratic_fitting_combined_criteria) = rownames(bisquare_female_quadratic_fitting_combined_criteria) = rownames(OLS_female_quadratic_fitting_combined_complex_criteria) = rownames(hampel_female_quadratic_fitting_combined_complex_criteria) = rownames(bisquare_female_quadratic_fitting_combined_complex_criteria) = c_grid


# Cubic TL
OLS_list_female_cubic_combined = hampel_list_female_cubic_combined = bisquare_list_female_cubic_combined = list()
OLS_list_female_cubic_combined_complex = hampel_list_female_cubic_combined_complex = bisquare_list_female_cubic_combined_complex = list()

for(i in 1:21)
{
  y = log(diag(female_long_run_spatial_combined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_female_cubic_combined[[i]] = lm(y ~ x + I(x^2) + I(x^3))
  
  
  hampel_list_female_cubic_combined[[i]] = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_female_cubic_combined[[i]] = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)
}

for(i in 1:21)
{
  y = log(diag(female_long_run_spatial_combined_complex[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_female_cubic_combined_complex[[i]] = lm(y ~ x + I(x^2) + I(x^3))
  
  
  hampel_list_female_cubic_combined_complex[[i]] = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_female_cubic_combined_complex[[i]] = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)
}


OLS_female_cubic_fitting_combined_criteria = OLS_female_cubic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_female_cubic_fitting_combined_criteria = hampel_female_cubic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)
bisquare_female_cubic_fitting_combined_criteria = bisquare_female_cubic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)

for(i in 1:21)
{
  # OLS
  OLS_female_cubic_fitting_combined_criteria[i, "R_squared"] = summary(OLS_list_female_cubic_combined[[i]])$adj.r.squared
  OLS_female_cubic_fitting_combined_criteria[i, "AIC"] = AIC(OLS_list_female_cubic_combined[[i]])
  OLS_female_cubic_fitting_combined_criteria[i, "AICc"] = AICc(OLS_list_female_cubic_combined[[i]])
  OLS_female_cubic_fitting_combined_criteria[i, "R_squared_pred"] = press(OLS_list_female_cubic_combined[[i]], as.R2 = TRUE)
  
  OLS_female_cubic_fitting_combined_complex_criteria[i, "R_squared"] = summary(OLS_list_female_cubic_combined_complex[[i]])$adj.r.squared
  OLS_female_cubic_fitting_combined_complex_criteria[i, "AIC"] = AIC(OLS_list_female_cubic_combined_complex[[i]])
  OLS_female_cubic_fitting_combined_complex_criteria[i, "AICc"] = AICc(OLS_list_female_cubic_combined_complex[[i]])
  OLS_female_cubic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(OLS_list_female_cubic_combined_complex[[i]], as.R2 = TRUE)
  
  # hample
  hampel_female_cubic_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((hampel_list_female_cubic_combined[[i]]$residuals)^2)/sum((hampel_list_female_cubic_combined[[i]]$model[,1] - mean(hampel_list_female_cubic_combined[[i]]$model[,1]))^2))*(100)/(97)
  
  hampel_female_cubic_fitting_combined_criteria[i, "AIC"] = AIC(hampel_list_female_cubic_combined[[i]])
  hampel_female_cubic_fitting_combined_criteria[i, "AICc"] = AICc(hampel_list_female_cubic_combined[[i]])
  hampel_female_cubic_fitting_combined_criteria[i, "R_squared_pred"] = press(hampel_list_female_cubic_combined[[i]], as.R2 = TRUE)
  
  hampel_female_cubic_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((hampel_list_female_cubic_combined_complex[[i]]$residuals)^2)/sum((hampel_list_female_cubic_combined_complex[[i]]$model[,1] - mean(hampel_list_female_cubic_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  hampel_female_cubic_fitting_combined_complex_criteria[i, "AIC"] = AIC(hampel_list_female_cubic_combined_complex[[i]])
  hampel_female_cubic_fitting_combined_complex_criteria[i, "AICc"] = AICc(hampel_list_female_cubic_combined_complex[[i]])
  hampel_female_cubic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(hampel_list_female_cubic_combined_complex[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_female_cubic_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_female_cubic_combined[[i]]$residuals)^2)/sum((bisquare_list_female_cubic_combined[[i]]$model[,1] - mean(bisquare_list_female_cubic_combined[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_female_cubic_fitting_combined_criteria[i, "AIC"] = AIC(bisquare_list_female_cubic_combined[[i]])
  bisquare_female_cubic_fitting_combined_criteria[i, "AICc"] = AICc(bisquare_list_female_cubic_combined[[i]])
  bisquare_female_cubic_fitting_combined_criteria[i, "R_squared_pred"] = press(bisquare_list_female_cubic_combined[[i]], as.R2 = TRUE)
  
  bisquare_female_cubic_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_female_cubic_combined_complex[[i]]$residuals)^2)/sum((bisquare_list_female_cubic_combined_complex[[i]]$model[,1] - mean(bisquare_list_female_cubic_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_female_cubic_fitting_combined_complex_criteria[i, "AIC"] = AIC(bisquare_list_female_cubic_combined_complex[[i]])
  bisquare_female_cubic_fitting_combined_complex_criteria[i, "AICc"] = AICc(bisquare_list_female_cubic_combined_complex[[i]])
  bisquare_female_cubic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(bisquare_list_female_cubic_combined_complex[[i]], as.R2 = TRUE)
}

rownames(OLS_female_cubic_fitting_combined_criteria) = rownames(hampel_female_cubic_fitting_combined_criteria) = rownames(bisquare_female_cubic_fitting_combined_criteria) = rownames(OLS_female_cubic_fitting_combined_complex_criteria) = rownames(hampel_female_cubic_fitting_combined_complex_criteria) = rownames(bisquare_female_cubic_fitting_combined_complex_criteria) = c_grid


# male

# Linear TL

OLS_list_male_combined = hampel_list_male_combined = bisquare_list_male_combined = list()
OLS_list_male_combined_complex = hampel_list_male_combined_complex = bisquare_list_male_combined_complex = list()

for(i in 1:21)
{
  y = log(diag(male_long_run_spatial_combined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_male_combined[[i]] = lm(y ~ x)
  
  
  hampel_list_male_combined[[i]] = rlm(y ~ x, psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_male_combined[[i]] = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
}

for(i in 1:21)
{
  y = log(diag(male_long_run_spatial_combined_complex[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_male_combined_complex[[i]] = lm(y ~ x)
  
  
  hampel_list_male_combined_complex[[i]] = rlm(y ~ x, psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_male_combined_complex[[i]] = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
}

OLS_male_fitting_combined_criteria = OLS_male_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_male_fitting_combined_criteria = hampel_male_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)
bisquare_male_fitting_combined_criteria = bisquare_male_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)


for(i in 1:21)
{
  # OLS
  OLS_male_fitting_combined_criteria[i, "R_squared"] = summary(OLS_list_male_combined[[i]])$adj.r.squared
  OLS_male_fitting_combined_criteria[i, "AIC"] = AIC(OLS_list_male_combined[[i]])
  OLS_male_fitting_combined_criteria[i, "AICc"] = AICc(OLS_list_male_combined[[i]])
  OLS_male_fitting_combined_criteria[i, "R_squared_pred"] = press(OLS_list_male_combined[[i]], as.R2 = TRUE)
  
  OLS_male_fitting_combined_complex_criteria[i, "R_squared"] = summary(OLS_list_male_combined_complex[[i]])$adj.r.squared
  OLS_male_fitting_combined_complex_criteria[i, "AIC"] = AIC(OLS_list_male_combined_complex[[i]])
  OLS_male_fitting_combined_complex_criteria[i, "AICc"] = AICc(OLS_list_male_combined_complex[[i]])
  OLS_male_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(OLS_list_male_combined_complex[[i]], as.R2 = TRUE)
  
  # hample
  hampel_male_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((hampel_list_male_combined[[i]]$residuals)^2)/sum((hampel_list_male_combined[[i]]$model[,1] - mean(hampel_list_male_combined[[i]]$model[,1]))^2))*(100)/(97)
  
  hampel_male_fitting_combined_criteria[i, "AIC"] = AIC(hampel_list_male_combined[[i]])
  hampel_male_fitting_combined_criteria[i, "AICc"] = AICc(hampel_list_male_combined[[i]])
  hampel_male_fitting_combined_criteria[i, "R_squared_pred"] = press(hampel_list_male_combined[[i]], as.R2 = TRUE)
  
  hampel_male_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((hampel_list_male_combined_complex[[i]]$residuals)^2)/sum((hampel_list_male_combined_complex[[i]]$model[,1] - mean(hampel_list_male_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  hampel_male_fitting_combined_complex_criteria[i, "AIC"] = AIC(hampel_list_male_combined_complex[[i]])
  hampel_male_fitting_combined_complex_criteria[i, "AICc"] = AICc(hampel_list_male_combined_complex[[i]])
  hampel_male_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(hampel_list_male_combined_complex[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_male_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_male_combined[[i]]$residuals)^2)/sum((bisquare_list_male_combined[[i]]$model[,1] - mean(bisquare_list_male_combined[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_male_fitting_combined_criteria[i, "AIC"] = AIC(bisquare_list_male_combined[[i]])
  bisquare_male_fitting_combined_criteria[i, "AICc"] = AICc(bisquare_list_male_combined[[i]])
  bisquare_male_fitting_combined_criteria[i, "R_squared_pred"] = press(bisquare_list_male_combined[[i]], as.R2 = TRUE)
  
  bisquare_male_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_male_combined_complex[[i]]$residuals)^2)/sum((bisquare_list_male_combined_complex[[i]]$model[,1] - mean(bisquare_list_male_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_male_fitting_combined_complex_criteria[i, "AIC"] = AIC(bisquare_list_male_combined_complex[[i]])
  bisquare_male_fitting_combined_complex_criteria[i, "AICc"] = AICc(bisquare_list_male_combined_complex[[i]])
  bisquare_male_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(bisquare_list_male_combined_complex[[i]], as.R2 = TRUE)
}

rownames(OLS_male_fitting_combined_criteria) = rownames(hampel_male_fitting_combined_criteria) = rownames(bisquare_male_fitting_combined_criteria) = rownames(OLS_male_fitting_combined_complex_criteria) = rownames(hampel_male_fitting_combined_complex_criteria) = rownames(bisquare_male_fitting_combined_complex_criteria) = c_grid

# Quadratic TL
OLS_list_male_quadratic_combined = hampel_list_male_quadratic_combined = bisquare_list_male_quadratic_combined = list()
OLS_list_male_quadratic_combined_complex = hampel_list_male_quadratic_combined_complex = bisquare_list_male_quadratic_combined_complex = list()

for(i in 1:21)
{
  y = log(diag(male_long_run_spatial_combined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_male_quadratic_combined[[i]] = lm(y ~ x + I(x^2))
  
  
  hampel_list_male_quadratic_combined[[i]] = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_male_quadratic_combined[[i]] = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
}

for(i in 1:21)
{
  y = log(diag(male_long_run_spatial_combined_complex[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_male_quadratic_combined_complex[[i]] = lm(y ~ x + I(x^2))
  
  
  hampel_list_male_quadratic_combined_complex[[i]] = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_male_quadratic_combined_complex[[i]] = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
}


OLS_male_quadratic_fitting_combined_criteria = OLS_male_quadratic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_male_quadratic_fitting_combined_criteria = hampel_male_quadratic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)
bisquare_male_quadratic_fitting_combined_criteria = bisquare_male_quadratic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)

for(i in 1:21)
{
  # OLS
  OLS_male_quadratic_fitting_combined_criteria[i, "R_squared"] = summary(OLS_list_male_quadratic_combined[[i]])$adj.r.squared
  OLS_male_quadratic_fitting_combined_criteria[i, "AIC"] = AIC(OLS_list_male_quadratic_combined[[i]])
  OLS_male_quadratic_fitting_combined_criteria[i, "AICc"] = AICc(OLS_list_male_quadratic_combined[[i]])
  OLS_male_quadratic_fitting_combined_criteria[i, "R_squared_pred"] = press(OLS_list_male_quadratic_combined[[i]], as.R2 = TRUE)
  
  OLS_male_quadratic_fitting_combined_complex_criteria[i, "R_squared"] = summary(OLS_list_male_quadratic_combined_complex[[i]])$adj.r.squared
  OLS_male_quadratic_fitting_combined_complex_criteria[i, "AIC"] = AIC(OLS_list_male_quadratic_combined_complex[[i]])
  OLS_male_quadratic_fitting_combined_complex_criteria[i, "AICc"] = AICc(OLS_list_male_quadratic_combined_complex[[i]])
  OLS_male_quadratic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(OLS_list_male_quadratic_combined_complex[[i]], as.R2 = TRUE)
  
  # hample
  hampel_male_quadratic_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((hampel_list_male_quadratic_combined[[i]]$residuals)^2)/sum((hampel_list_male_quadratic_combined[[i]]$model[,1] - mean(hampel_list_male_quadratic_combined[[i]]$model[,1]))^2))*(100)/(97)
  
  hampel_male_quadratic_fitting_combined_criteria[i, "AIC"] = AIC(hampel_list_male_quadratic_combined[[i]])
  hampel_male_quadratic_fitting_combined_criteria[i, "AICc"] = AICc(hampel_list_male_quadratic_combined[[i]])
  hampel_male_quadratic_fitting_combined_criteria[i, "R_squared_pred"] = press(hampel_list_male_quadratic_combined[[i]], as.R2 = TRUE)
  
  hampel_male_quadratic_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((hampel_list_male_quadratic_combined_complex[[i]]$residuals)^2)/sum((hampel_list_male_quadratic_combined_complex[[i]]$model[,1] - mean(hampel_list_male_quadratic_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  hampel_male_quadratic_fitting_combined_complex_criteria[i, "AIC"] = AIC(hampel_list_male_quadratic_combined_complex[[i]])
  hampel_male_quadratic_fitting_combined_complex_criteria[i, "AICc"] = AICc(hampel_list_male_quadratic_combined_complex[[i]])
  hampel_male_quadratic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(hampel_list_male_quadratic_combined_complex[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_male_quadratic_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_male_quadratic_combined[[i]]$residuals)^2)/sum((bisquare_list_male_quadratic_combined[[i]]$model[,1] - mean(bisquare_list_male_quadratic_combined[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_male_quadratic_fitting_combined_criteria[i, "AIC"] = AIC(bisquare_list_male_quadratic_combined[[i]])
  bisquare_male_quadratic_fitting_combined_criteria[i, "AICc"] = AICc(bisquare_list_male_quadratic_combined[[i]])
  bisquare_male_quadratic_fitting_combined_criteria[i, "R_squared_pred"] = press(bisquare_list_male_quadratic_combined[[i]], as.R2 = TRUE)
  
  bisquare_male_quadratic_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_male_quadratic_combined_complex[[i]]$residuals)^2)/sum((bisquare_list_male_quadratic_combined_complex[[i]]$model[,1] - mean(bisquare_list_male_quadratic_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_male_quadratic_fitting_combined_complex_criteria[i, "AIC"] = AIC(bisquare_list_male_quadratic_combined_complex[[i]])
  bisquare_male_quadratic_fitting_combined_complex_criteria[i, "AICc"] = AICc(bisquare_list_male_quadratic_combined_complex[[i]])
  bisquare_male_quadratic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(bisquare_list_male_quadratic_combined_complex[[i]], as.R2 = TRUE)
}

rownames(OLS_male_quadratic_fitting_combined_criteria) = rownames(hampel_male_quadratic_fitting_combined_criteria) = rownames(bisquare_male_quadratic_fitting_combined_criteria) = rownames(OLS_male_quadratic_fitting_combined_complex_criteria) = rownames(hampel_male_quadratic_fitting_combined_complex_criteria) = rownames(bisquare_male_quadratic_fitting_combined_complex_criteria) = c_grid


# Cubic TL
OLS_list_male_cubic_combined = hampel_list_male_cubic_combined = bisquare_list_male_cubic_combined = list()
OLS_list_male_cubic_combined_complex = hampel_list_male_cubic_combined_complex = bisquare_list_male_cubic_combined_complex = list()

for(i in 1:21)
{
  y = log(diag(male_long_run_spatial_combined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_male_cubic_combined[[i]] = lm(y ~ x + I(x^2) + I(x^3))
  
  
  hampel_list_male_cubic_combined[[i]] = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_male_cubic_combined[[i]] = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)
}

for(i in 1:21)
{
  y = log(diag(male_long_run_spatial_combined_complex[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_male_cubic_combined_complex[[i]] = lm(y ~ x + I(x^2) + I(x^3))
  
  
  hampel_list_male_cubic_combined_complex[[i]] = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_male_cubic_combined_complex[[i]] = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)
}


OLS_male_cubic_fitting_combined_criteria = OLS_male_cubic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_male_cubic_fitting_combined_criteria = hampel_male_cubic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)
bisquare_male_cubic_fitting_combined_criteria = bisquare_male_cubic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)

for(i in 1:21)
{
  # OLS
  OLS_male_cubic_fitting_combined_criteria[i, "R_squared"] = summary(OLS_list_male_cubic_combined[[i]])$adj.r.squared
  OLS_male_cubic_fitting_combined_criteria[i, "AIC"] = AIC(OLS_list_male_cubic_combined[[i]])
  OLS_male_cubic_fitting_combined_criteria[i, "AICc"] = AICc(OLS_list_male_cubic_combined[[i]])
  OLS_male_cubic_fitting_combined_criteria[i, "R_squared_pred"] = press(OLS_list_male_cubic_combined[[i]], as.R2 = TRUE)
  
  OLS_male_cubic_fitting_combined_complex_criteria[i, "R_squared"] = summary(OLS_list_male_cubic_combined_complex[[i]])$adj.r.squared
  OLS_male_cubic_fitting_combined_complex_criteria[i, "AIC"] = AIC(OLS_list_male_cubic_combined_complex[[i]])
  OLS_male_cubic_fitting_combined_complex_criteria[i, "AICc"] = AICc(OLS_list_male_cubic_combined_complex[[i]])
  OLS_male_cubic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(OLS_list_male_cubic_combined_complex[[i]], as.R2 = TRUE)
  
  # hample
  hampel_male_cubic_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((hampel_list_male_cubic_combined[[i]]$residuals)^2)/sum((hampel_list_male_cubic_combined[[i]]$model[,1] - mean(hampel_list_male_cubic_combined[[i]]$model[,1]))^2))*(100)/(97)
  
  hampel_male_cubic_fitting_combined_criteria[i, "AIC"] = AIC(hampel_list_male_cubic_combined[[i]])
  hampel_male_cubic_fitting_combined_criteria[i, "AICc"] = AICc(hampel_list_male_cubic_combined[[i]])
  hampel_male_cubic_fitting_combined_criteria[i, "R_squared_pred"] = press(hampel_list_male_cubic_combined[[i]], as.R2 = TRUE)
  
  hampel_male_cubic_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((hampel_list_male_cubic_combined_complex[[i]]$residuals)^2)/sum((hampel_list_male_cubic_combined_complex[[i]]$model[,1] - mean(hampel_list_male_cubic_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  hampel_male_cubic_fitting_combined_complex_criteria[i, "AIC"] = AIC(hampel_list_male_cubic_combined_complex[[i]])
  hampel_male_cubic_fitting_combined_complex_criteria[i, "AICc"] = AICc(hampel_list_male_cubic_combined_complex[[i]])
  hampel_male_cubic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(hampel_list_male_cubic_combined_complex[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_male_cubic_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_male_cubic_combined[[i]]$residuals)^2)/sum((bisquare_list_male_cubic_combined[[i]]$model[,1] - mean(bisquare_list_male_cubic_combined[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_male_cubic_fitting_combined_criteria[i, "AIC"] = AIC(bisquare_list_male_cubic_combined[[i]])
  bisquare_male_cubic_fitting_combined_criteria[i, "AICc"] = AICc(bisquare_list_male_cubic_combined[[i]])
  bisquare_male_cubic_fitting_combined_criteria[i, "R_squared_pred"] = press(bisquare_list_male_cubic_combined[[i]], as.R2 = TRUE)
  
  bisquare_male_cubic_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_male_cubic_combined_complex[[i]]$residuals)^2)/sum((bisquare_list_male_cubic_combined_complex[[i]]$model[,1] - mean(bisquare_list_male_cubic_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_male_cubic_fitting_combined_complex_criteria[i, "AIC"] = AIC(bisquare_list_male_cubic_combined_complex[[i]])
  bisquare_male_cubic_fitting_combined_complex_criteria[i, "AICc"] = AICc(bisquare_list_male_cubic_combined_complex[[i]])
  bisquare_male_cubic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(bisquare_list_male_cubic_combined_complex[[i]], as.R2 = TRUE)
}

rownames(OLS_male_cubic_fitting_combined_criteria) = rownames(hampel_male_cubic_fitting_combined_criteria) = rownames(bisquare_male_cubic_fitting_combined_criteria) = rownames(OLS_male_cubic_fitting_combined_complex_criteria) = rownames(hampel_male_cubic_fitting_combined_complex_criteria) = rownames(bisquare_male_cubic_fitting_combined_complex_criteria) = c_grid










###################################################


# female

OLS_model = lm(log(diag(female_long_run_spatial_var$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$female), base = 10))

female_TL_long_run_spatial = OLS_model$coefficients
female_TL_long_run_spatial_R_squared = summary(OLS_model)$r.squared

hampel_model = rlm(log(diag(female_long_run_spatial_var$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$female), base = 10), psi = psi.hampel, maxit = 200)
female_TL_hampel_long_run_spatial = hampel_model$coefficients
female_TL_hampel_long_run_spatial_R_squared = 1 - sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2)

bisquare_model = rlm(log(diag(female_long_run_spatial_var$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$female), base = 10), psi = psi.bisquare, maxit = 200)

female_TL_bisquare_long_run_spatial = bisquare_model$coefficients
female_TL_bisquare_long_run_spatial_R_squared = 1 - sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2)

# male

OLS_model = lm(log(diag(male_long_run_spatial_var$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$male), base = 10))

male_TL_long_run_spatial = OLS_model$coefficients
male_TL_long_run_spatial_R_squared = summary(OLS_model)$r.squared

hampel_model = rlm(log(diag(male_long_run_spatial_var$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$male), base = 10), psi = psi.hampel, maxit = 200)
male_TL_hampel_long_run_spatial = hampel_model$coefficients
male_TL_hampel_long_run_spatial_R_squared = 1 - sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2)

bisquare_model = rlm(log(diag(male_long_run_spatial_var$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$male), base = 10), psi = psi.bisquare, maxit = 200)

male_TL_bisquare_long_run_spatial = bisquare_model$coefficients
male_TL_bisquare_long_run_spatial_R_squared = 1 - sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2)

# total

OLS_model = lm(log(diag(total_long_run_spatial_var$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$total), base = 10))

total_TL_long_run_spatial = OLS_model$coefficients
total_TL_long_run_spatial_R_squared = summary(OLS_model)$r.squared

hampel_model = rlm(log(diag(total_long_run_spatial_var$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$total), base = 10), psi = psi.hampel, maxit = 200)
total_TL_hampel_long_run_spatial = hampel_model$coefficients
total_TL_hampel_long_run_spatial_R_squared = 1 - sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2)

bisquare_model = rlm(log(diag(total_long_run_spatial_var$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$total), base = 10), psi = psi.bisquare, maxit = 200)

total_TL_bisquare_long_run_spatial = bisquare_model$coefficients
total_TL_bisquare_long_run_spatial_R_squared = 1 - sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2)

# total

# Linear TL

OLS_list_total_combined = hampel_list_total_combined = bisquare_list_total_combined = list()
OLS_list_total_combined_complex = hampel_list_total_combined_complex = bisquare_list_total_combined_complex = list()

for(i in 1:21)
{
  y = log(diag(total_long_run_spatial_combined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_total_combined[[i]] = lm(y ~ x)
  
  
  hampel_list_total_combined[[i]] = rlm(y ~ x, psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_total_combined[[i]] = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
}

for(i in 1:21)
{
  y = log(diag(total_long_run_spatial_combined_complex[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_total_combined_complex[[i]] = lm(y ~ x)
  
  
  hampel_list_total_combined_complex[[i]] = rlm(y ~ x, psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_total_combined_complex[[i]] = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
}

OLS_total_fitting_combined_criteria = OLS_total_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_total_fitting_combined_criteria = hampel_total_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)
bisquare_total_fitting_combined_criteria = bisquare_total_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)


for(i in 1:21)
{
  # OLS
  OLS_total_fitting_combined_criteria[i, "R_squared"] = summary(OLS_list_total_combined[[i]])$adj.r.squared
  OLS_total_fitting_combined_criteria[i, "AIC"] = AIC(OLS_list_total_combined[[i]])
  OLS_total_fitting_combined_criteria[i, "AICc"] = AICc(OLS_list_total_combined[[i]])
  OLS_total_fitting_combined_criteria[i, "R_squared_pred"] = press(OLS_list_total_combined[[i]], as.R2 = TRUE)
  
  OLS_total_fitting_combined_complex_criteria[i, "R_squared"] = summary(OLS_list_total_combined_complex[[i]])$adj.r.squared
  OLS_total_fitting_combined_complex_criteria[i, "AIC"] = AIC(OLS_list_total_combined_complex[[i]])
  OLS_total_fitting_combined_complex_criteria[i, "AICc"] = AICc(OLS_list_total_combined_complex[[i]])
  OLS_total_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(OLS_list_total_combined_complex[[i]], as.R2 = TRUE)
  
  # hample
  hampel_total_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((hampel_list_total_combined[[i]]$residuals)^2)/sum((hampel_list_total_combined[[i]]$model[,1] - mean(hampel_list_total_combined[[i]]$model[,1]))^2))*(100)/(97)
  
  hampel_total_fitting_combined_criteria[i, "AIC"] = AIC(hampel_list_total_combined[[i]])
  hampel_total_fitting_combined_criteria[i, "AICc"] = AICc(hampel_list_total_combined[[i]])
  hampel_total_fitting_combined_criteria[i, "R_squared_pred"] = press(hampel_list_total_combined[[i]], as.R2 = TRUE)
  
  hampel_total_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((hampel_list_total_combined_complex[[i]]$residuals)^2)/sum((hampel_list_total_combined_complex[[i]]$model[,1] - mean(hampel_list_total_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  hampel_total_fitting_combined_complex_criteria[i, "AIC"] = AIC(hampel_list_total_combined_complex[[i]])
  hampel_total_fitting_combined_complex_criteria[i, "AICc"] = AICc(hampel_list_total_combined_complex[[i]])
  hampel_total_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(hampel_list_total_combined_complex[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_total_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_total_combined[[i]]$residuals)^2)/sum((bisquare_list_total_combined[[i]]$model[,1] - mean(bisquare_list_total_combined[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_total_fitting_combined_criteria[i, "AIC"] = AIC(bisquare_list_total_combined[[i]])
  bisquare_total_fitting_combined_criteria[i, "AICc"] = AICc(bisquare_list_total_combined[[i]])
  bisquare_total_fitting_combined_criteria[i, "R_squared_pred"] = press(bisquare_list_total_combined[[i]], as.R2 = TRUE)
  
  bisquare_total_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_total_combined_complex[[i]]$residuals)^2)/sum((bisquare_list_total_combined_complex[[i]]$model[,1] - mean(bisquare_list_total_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_total_fitting_combined_complex_criteria[i, "AIC"] = AIC(bisquare_list_total_combined_complex[[i]])
  bisquare_total_fitting_combined_complex_criteria[i, "AICc"] = AICc(bisquare_list_total_combined_complex[[i]])
  bisquare_total_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(bisquare_list_total_combined_complex[[i]], as.R2 = TRUE)
}

rownames(OLS_total_fitting_combined_criteria) = rownames(hampel_total_fitting_combined_criteria) = rownames(bisquare_total_fitting_combined_criteria) = rownames(OLS_total_fitting_combined_complex_criteria) = rownames(hampel_total_fitting_combined_complex_criteria) = rownames(bisquare_total_fitting_combined_complex_criteria) = c_grid

# Quadratic TL
OLS_list_total_quadratic_combined = hampel_list_total_quadratic_combined = bisquare_list_total_quadratic_combined = list()
OLS_list_total_quadratic_combined_complex = hampel_list_total_quadratic_combined_complex = bisquare_list_total_quadratic_combined_complex = list()

for(i in 1:21)
{
  y = log(diag(total_long_run_spatial_combined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_total_quadratic_combined[[i]] = lm(y ~ x + I(x^2))
  
  
  hampel_list_total_quadratic_combined[[i]] = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_total_quadratic_combined[[i]] = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
}

for(i in 1:21)
{
  y = log(diag(total_long_run_spatial_combined_complex[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_total_quadratic_combined_complex[[i]] = lm(y ~ x + I(x^2))
  
  
  hampel_list_total_quadratic_combined_complex[[i]] = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_total_quadratic_combined_complex[[i]] = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
}


OLS_total_quadratic_fitting_combined_criteria = OLS_total_quadratic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_total_quadratic_fitting_combined_criteria = hampel_total_quadratic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)
bisquare_total_quadratic_fitting_combined_criteria = bisquare_total_quadratic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)

for(i in 1:21)
{
  # OLS
  OLS_total_quadratic_fitting_combined_criteria[i, "R_squared"] = summary(OLS_list_total_quadratic_combined[[i]])$adj.r.squared
  OLS_total_quadratic_fitting_combined_criteria[i, "AIC"] = AIC(OLS_list_total_quadratic_combined[[i]])
  OLS_total_quadratic_fitting_combined_criteria[i, "AICc"] = AICc(OLS_list_total_quadratic_combined[[i]])
  OLS_total_quadratic_fitting_combined_criteria[i, "R_squared_pred"] = press(OLS_list_total_quadratic_combined[[i]], as.R2 = TRUE)
  
  OLS_total_quadratic_fitting_combined_complex_criteria[i, "R_squared"] = summary(OLS_list_total_quadratic_combined_complex[[i]])$adj.r.squared
  OLS_total_quadratic_fitting_combined_complex_criteria[i, "AIC"] = AIC(OLS_list_total_quadratic_combined_complex[[i]])
  OLS_total_quadratic_fitting_combined_complex_criteria[i, "AICc"] = AICc(OLS_list_total_quadratic_combined_complex[[i]])
  OLS_total_quadratic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(OLS_list_total_quadratic_combined_complex[[i]], as.R2 = TRUE)
  
  # hample
  hampel_total_quadratic_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((hampel_list_total_quadratic_combined[[i]]$residuals)^2)/sum((hampel_list_total_quadratic_combined[[i]]$model[,1] - mean(hampel_list_total_quadratic_combined[[i]]$model[,1]))^2))*(100)/(97)
  
  hampel_total_quadratic_fitting_combined_criteria[i, "AIC"] = AIC(hampel_list_total_quadratic_combined[[i]])
  hampel_total_quadratic_fitting_combined_criteria[i, "AICc"] = AICc(hampel_list_total_quadratic_combined[[i]])
  hampel_total_quadratic_fitting_combined_criteria[i, "R_squared_pred"] = press(hampel_list_total_quadratic_combined[[i]], as.R2 = TRUE)
  
  hampel_total_quadratic_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((hampel_list_total_quadratic_combined_complex[[i]]$residuals)^2)/sum((hampel_list_total_quadratic_combined_complex[[i]]$model[,1] - mean(hampel_list_total_quadratic_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  hampel_total_quadratic_fitting_combined_complex_criteria[i, "AIC"] = AIC(hampel_list_total_quadratic_combined_complex[[i]])
  hampel_total_quadratic_fitting_combined_complex_criteria[i, "AICc"] = AICc(hampel_list_total_quadratic_combined_complex[[i]])
  hampel_total_quadratic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(hampel_list_total_quadratic_combined_complex[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_total_quadratic_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_total_quadratic_combined[[i]]$residuals)^2)/sum((bisquare_list_total_quadratic_combined[[i]]$model[,1] - mean(bisquare_list_total_quadratic_combined[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_total_quadratic_fitting_combined_criteria[i, "AIC"] = AIC(bisquare_list_total_quadratic_combined[[i]])
  bisquare_total_quadratic_fitting_combined_criteria[i, "AICc"] = AICc(bisquare_list_total_quadratic_combined[[i]])
  bisquare_total_quadratic_fitting_combined_criteria[i, "R_squared_pred"] = press(bisquare_list_total_quadratic_combined[[i]], as.R2 = TRUE)
  
  bisquare_total_quadratic_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_total_quadratic_combined_complex[[i]]$residuals)^2)/sum((bisquare_list_total_quadratic_combined_complex[[i]]$model[,1] - mean(bisquare_list_total_quadratic_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_total_quadratic_fitting_combined_complex_criteria[i, "AIC"] = AIC(bisquare_list_total_quadratic_combined_complex[[i]])
  bisquare_total_quadratic_fitting_combined_complex_criteria[i, "AICc"] = AICc(bisquare_list_total_quadratic_combined_complex[[i]])
  bisquare_total_quadratic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(bisquare_list_total_quadratic_combined_complex[[i]], as.R2 = TRUE)
}

rownames(OLS_total_quadratic_fitting_combined_criteria) = rownames(hampel_total_quadratic_fitting_combined_criteria) = rownames(bisquare_total_quadratic_fitting_combined_criteria) = rownames(OLS_total_quadratic_fitting_combined_complex_criteria) = rownames(hampel_total_quadratic_fitting_combined_complex_criteria) = rownames(bisquare_total_quadratic_fitting_combined_complex_criteria) = c_grid


# Cubic TL
OLS_list_total_cubic_combined = hampel_list_total_cubic_combined = bisquare_list_total_cubic_combined = list()
OLS_list_total_cubic_combined_complex = hampel_list_total_cubic_combined_complex = bisquare_list_total_cubic_combined_complex = list()

for(i in 1:21)
{
  y = log(diag(total_long_run_spatial_combined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_total_cubic_combined[[i]] = lm(y ~ x + I(x^2) + I(x^3))
  
  
  hampel_list_total_cubic_combined[[i]] = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_total_cubic_combined[[i]] = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)
}

for(i in 1:21)
{
  y = log(diag(total_long_run_spatial_combined_complex[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10)
  x = log(rowMeans(Japan$rate$total), base = 10)
  
  OLS_list_total_cubic_combined_complex[[i]] = lm(y ~ x + I(x^2) + I(x^3))
  
  
  hampel_list_total_cubic_combined_complex[[i]] = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)
  
  
  bisquare_list_total_cubic_combined_complex[[i]] = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)
}


OLS_total_cubic_fitting_combined_criteria = OLS_total_cubic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_total_cubic_fitting_combined_criteria = hampel_total_cubic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)
bisquare_total_cubic_fitting_combined_criteria = bisquare_total_cubic_fitting_combined_complex_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA)

for(i in 1:21)
{
  # OLS
  OLS_total_cubic_fitting_combined_criteria[i, "R_squared"] = summary(OLS_list_total_cubic_combined[[i]])$adj.r.squared
  OLS_total_cubic_fitting_combined_criteria[i, "AIC"] = AIC(OLS_list_total_cubic_combined[[i]])
  OLS_total_cubic_fitting_combined_criteria[i, "AICc"] = AICc(OLS_list_total_cubic_combined[[i]])
  OLS_total_cubic_fitting_combined_criteria[i, "R_squared_pred"] = press(OLS_list_total_cubic_combined[[i]], as.R2 = TRUE)
  
  OLS_total_cubic_fitting_combined_complex_criteria[i, "R_squared"] = summary(OLS_list_total_cubic_combined_complex[[i]])$adj.r.squared
  OLS_total_cubic_fitting_combined_complex_criteria[i, "AIC"] = AIC(OLS_list_total_cubic_combined_complex[[i]])
  OLS_total_cubic_fitting_combined_complex_criteria[i, "AICc"] = AICc(OLS_list_total_cubic_combined_complex[[i]])
  OLS_total_cubic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(OLS_list_total_cubic_combined_complex[[i]], as.R2 = TRUE)
  
  # hample
  hampel_total_cubic_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((hampel_list_total_cubic_combined[[i]]$residuals)^2)/sum((hampel_list_total_cubic_combined[[i]]$model[,1] - mean(hampel_list_total_cubic_combined[[i]]$model[,1]))^2))*(100)/(97)
  
  hampel_total_cubic_fitting_combined_criteria[i, "AIC"] = AIC(hampel_list_total_cubic_combined[[i]])
  hampel_total_cubic_fitting_combined_criteria[i, "AICc"] = AICc(hampel_list_total_cubic_combined[[i]])
  hampel_total_cubic_fitting_combined_criteria[i, "R_squared_pred"] = press(hampel_list_total_cubic_combined[[i]], as.R2 = TRUE)
  
  hampel_total_cubic_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((hampel_list_total_cubic_combined_complex[[i]]$residuals)^2)/sum((hampel_list_total_cubic_combined_complex[[i]]$model[,1] - mean(hampel_list_total_cubic_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  hampel_total_cubic_fitting_combined_complex_criteria[i, "AIC"] = AIC(hampel_list_total_cubic_combined_complex[[i]])
  hampel_total_cubic_fitting_combined_complex_criteria[i, "AICc"] = AICc(hampel_list_total_cubic_combined_complex[[i]])
  hampel_total_cubic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(hampel_list_total_cubic_combined_complex[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_total_cubic_fitting_combined_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_total_cubic_combined[[i]]$residuals)^2)/sum((bisquare_list_total_cubic_combined[[i]]$model[,1] - mean(bisquare_list_total_cubic_combined[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_total_cubic_fitting_combined_criteria[i, "AIC"] = AIC(bisquare_list_total_cubic_combined[[i]])
  bisquare_total_cubic_fitting_combined_criteria[i, "AICc"] = AICc(bisquare_list_total_cubic_combined[[i]])
  bisquare_total_cubic_fitting_combined_criteria[i, "R_squared_pred"] = press(bisquare_list_total_cubic_combined[[i]], as.R2 = TRUE)
  
  bisquare_total_cubic_fitting_combined_complex_criteria[i, "R_squared"] = 1 - (sum((bisquare_list_total_cubic_combined_complex[[i]]$residuals)^2)/sum((bisquare_list_total_cubic_combined_complex[[i]]$model[,1] - mean(bisquare_list_total_cubic_combined_complex[[i]]$model[,1]))^2))*(100)/(97)
  bisquare_total_cubic_fitting_combined_complex_criteria[i, "AIC"] = AIC(bisquare_list_total_cubic_combined_complex[[i]])
  bisquare_total_cubic_fitting_combined_complex_criteria[i, "AICc"] = AICc(bisquare_list_total_cubic_combined_complex[[i]])
  bisquare_total_cubic_fitting_combined_complex_criteria[i, "R_squared_pred"] = press(bisquare_list_total_cubic_combined_complex[[i]], as.R2 = TRUE)
}

rownames(OLS_total_cubic_fitting_combined_criteria) = rownames(hampel_total_cubic_fitting_combined_criteria) = rownames(bisquare_total_cubic_fitting_combined_criteria) = rownames(OLS_total_cubic_fitting_combined_complex_criteria) = rownames(hampel_total_cubic_fitting_combined_complex_criteria) = rownames(bisquare_total_cubic_fitting_combined_complex_criteria) = c_grid


# female_long_run_spatial_combined_all = male_long_run_spatial_combined_all = total_long_run_spatial_combined_all = list()

library(doParallel)
cl <- makeCluster(30) 
registerDoParallel(cl)

# use grid search to determine the coefficient c

c_grid = seq(0, 01, by = 0.05)

female_long_run_spatial_combined_all = foreach(ij = 1:21) %dopar% long_run_spatial_parallel(data = female_data[,,2:48], w_mat = weights_matrix_2, method = "Simple", c_all = c_grid, ik = ij)
male_long_run_spatial_combined_all = foreach(ij = 1:21) %dopar% long_run_spatial_parallel(data = male_data[,,2:48], w_mat = weights_matrix_2, method = "Simple", c_all = c_grid, ik = ij)
total_long_run_spatial_combined_all = foreach(ij = 1:21) %dopar% long_run_spatial_parallel(data = total_data[,,2:48], w_mat = weights_matrix_2, method = "Simple", c_all = c_grid, ik = ij)

save(female_long_run_spatial_combined_all, file = "female_long_run_spatial_combined_all.RData")
save(male_long_run_spatial_combined_all, file = "male_long_run_spatial_combined_all.RData")
save(total_long_run_spatial_combined_all, file = "total_long_run_spatial_combined_all.RData")

# use grid search to determine the coefficient c

c_grid_refined = seq(0.001, 0.05, length.out = 20)

library(doParallel)
cl <- makeCluster(30) 
registerDoParallel(cl)

female_long_run_spatial_combined_all_refined = foreach(ij = 1:20) %dopar% long_run_spatial_parallel(data = female_data[,,2:48], w_mat = weights_matrix_2, method = "Simple", c_all = c_grid_refined, ik = ij)
male_long_run_spatial_combined_all_refined = foreach(ij = 1:20) %dopar% long_run_spatial_parallel(data = male_data[,,2:48], w_mat = weights_matrix_2, method = "Simple", c_all = c_grid_refined, ik = ij)
total_long_run_spatial_combined_all_refined = foreach(ij = 1:20) %dopar% long_run_spatial_parallel(data = total_data[,,2:48], w_mat = weights_matrix_2, method = "Simple", c_all = c_grid_refined, ik = ij)

save(female_long_run_spatial_combined_all_refined, file = "female_long_run_spatial_combined_all_refined.RData")
save(male_long_run_spatial_combined_all_refined, file = "male_long_run_spatial_combined_all_refined.RData")
save(total_long_run_spatial_combined_all_refined, file = "total_long_run_spatial_combined_all_refined.RData")


###########################
# Search for best c_weight
##########################


## female

# c_weight from 0 to 1 by 0.05

OLS_list_female = hampel_list_female = bisquare_list_female = list()

for(i in 1:21)
{
  OLS_list_female[[i]] = lm(log(diag(female_long_run_spatial_combined_all[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$female), base = 10))
  hampel_list_female[[i]] = rlm(log(diag(female_long_run_spatial_combined_all[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$female), base = 10), psi = psi.hampel, maxit = 200)
  bisquare_list_female[[i]] = rlm(log(diag(female_long_run_spatial_combined_all[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$female), base = 10), psi = psi.bisquare, maxit = 200)
}

OLS_female_fitting_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_female_fitting_criteria = data.frame(R_squared = NA, AIC = NA, R_squared_pred = NA)
bisquare_female_fitting_criteria = data.frame(R_squared = NA, AIC = NA, R_squared_pred = NA)

for(i in 1:21)
{
  # OLS
  OLS_female_fitting_criteria[i, "R_squared"] = summary(OLS_list_female[[i]])$adj.r.squared
  OLS_female_fitting_criteria[i, "AIC"] = AIC(OLS_list_female[[i]])
  OLS_female_fitting_criteria[i, "AICc"] = AICc(OLS_list_female[[i]])
  OLS_female_fitting_criteria[i, "R_squared_pred"] = press(OLS_list_female[[i]], as.R2 = TRUE)
  
  # hample
  hampel_female_fitting_criteria[i, "R_squared"] = 1 - sum((hampel_list_female[[i]]$residuals)^2)/sum((hampel_list_female[[i]]$model[,1] - mean(hampel_list_female[[i]]$model[,1]))^2)
  hampel_female_fitting_criteria[i, "AIC"] = AIC(hampel_list_female[[i]])
  hampel_female_fitting_criteria[i, "AICc"] = AICc(hampel_list_female[[i]])
  hampel_female_fitting_criteria[i, "R_squared_pred"] = press(hampel_list_female[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_female_fitting_criteria[i, "R_squared"] = 1 - sum((bisquare_list_female[[i]]$residuals)^2)/sum((bisquare_list_female[[i]]$model[,1] - mean(bisquare_list_female[[i]]$model[,1]))^2)
  bisquare_female_fitting_criteria[i, "AIC"] = AIC(bisquare_list_female[[i]])
  bisquare_female_fitting_criteria[i, "AICc"] = AICc(bisquare_list_female[[i]])
  bisquare_female_fitting_criteria[i, "R_squared_pred"] = press(bisquare_list_female[[i]], as.R2 = TRUE)
}

rownames(OLS_female_fitting_criteria) = rownames(hampel_female_fitting_criteria) = rownames(bisquare_female_fitting_criteria) = c_grid

# c_weight from 0.001 to 0.05

OLS_list_female_refined = hampel_list_female_refined = bisquare_list_female_refined = list()

for(i in 1:20)
{
  OLS_list_female_refined[[i]] = lm(log(diag(female_long_run_spatial_combined_all_refined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$female), base = 10))
  hampel_list_female_refined[[i]] = rlm(log(diag(female_long_run_spatial_combined_all_refined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$female), base = 10), psi = psi.hampel, maxit = 200)
  bisquare_list_female_refined[[i]] = rlm(log(diag(female_long_run_spatial_combined_all_refined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$female), base = 10), psi = psi.bisquare, maxit = 200)
}


OLS_female_fitting_criteria_refined = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_female_fitting_criteria_refined = data.frame(R_squared = NA, AIC = NA, R_squared_pred = NA)
bisquare_female_fitting_criteria_refined = data.frame(R_squared = NA, AIC = NA, R_squared_pred = NA)

for(i in 1:20)
{
  # OLS
  OLS_female_fitting_criteria_refined[i, "R_squared"] = summary(OLS_list_female_refined[[i]])$adj.r.squared
  OLS_female_fitting_criteria_refined[i, "AIC"] = AIC(OLS_list_female_refined[[i]])
  OLS_female_fitting_criteria_refined[i, "AICc"] = AICc(OLS_list_female_refined[[i]])
  OLS_female_fitting_criteria_refined[i, "R_squared_pred"] = press(OLS_list_female_refined[[i]], as.R2 = TRUE)
  
  # hample
  hampel_female_fitting_criteria_refined[i, "R_squared"] = 1 - sum((hampel_list_female_refined[[i]]$residuals)^2)/sum((hampel_list_female_refined[[i]]$model[,1] - mean(hampel_list_female_refined[[i]]$model[,1]))^2)
  hampel_female_fitting_criteria_refined[i, "AIC"] = AIC(hampel_list_female_refined[[i]])
  hampel_female_fitting_criteria_refined[i, "AICc"] = AICc(hampel_list_female_refined[[i]])
  hampel_female_fitting_criteria_refined[i, "R_squared_pred"] = press(hampel_list_female_refined[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_female_fitting_criteria_refined[i, "R_squared"] = 1 - sum((bisquare_list_female_refined[[i]]$residuals)^2)/sum((bisquare_list_female_refined[[i]]$model[,1] - mean(bisquare_list_female_refined[[i]]$model[,1]))^2)
  bisquare_female_fitting_criteria_refined[i, "AIC"] = AIC(bisquare_list_female_refined[[i]])
  bisquare_female_fitting_criteria_refined[i, "AICc"] = AICc(bisquare_list_female_refined[[i]])
  bisquare_female_fitting_criteria_refined[i, "R_squared_pred"] = press(bisquare_list_female_refined[[i]], as.R2 = TRUE)
}

rownames(OLS_female_fitting_criteria_refined) = rownames(hampel_female_fitting_criteria_refined) = rownames(bisquare_female_fitting_criteria_refined) = c_grid_refined

# OLS
OLS_female_fitting_criteria_all = rbind(OLS_female_fitting_criteria_refined, OLS_female_fitting_criteria[2:21,])
plot(x = as.numeric(rownames(OLS_female_fitting_criteria_all)), y = OLS_female_fitting_criteria_all$AICc, xlab = "C_weight", ylab = "AICc", main = "OLS estimation; female")

# hampel
hampel_female_fitting_criteria_all = rbind(hampel_female_fitting_criteria_refined, hampel_female_fitting_criteria[2:21,])
plot(x = as.numeric(rownames(hampel_female_fitting_criteria_all)), y = hampel_female_fitting_criteria_all$AICc, xlab = "C_weight", ylab = "AICc", main = "hampel estimation; female")

# bisquare
bisquare_female_fitting_criteria_all = rbind(bisquare_female_fitting_criteria_refined, bisquare_female_fitting_criteria[2:21,])
plot(x = as.numeric(rownames(bisquare_female_fitting_criteria_all)), y = bisquare_female_fitting_criteria_all$AICc, xlab = "C_weight", ylab = "AICc", main = "bisquare estimation; female")

## male

# c_weight from 0 to 1 by 0.05

OLS_list_male = hampel_list_male = bisquare_list_male = list()

for(i in 1:21)
{
  OLS_list_male[[i]] = lm(log(diag(male_long_run_spatial_combined_all[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$male), base = 10))
  hampel_list_male[[i]] = rlm(log(diag(male_long_run_spatial_combined_all[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$male), base = 10), psi = psi.hampel, maxit = 200)
  bisquare_list_male[[i]] = rlm(log(diag(male_long_run_spatial_combined_all[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$male), base = 10), psi = psi.bisquare, maxit = 200)
}

OLS_male_fitting_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_male_fitting_criteria = data.frame(R_squared = NA, AIC = NA, R_squared_pred = NA)
bisquare_male_fitting_criteria = data.frame(R_squared = NA, AIC = NA, R_squared_pred = NA)

for(i in 1:21)
{
  # OLS
  OLS_male_fitting_criteria[i, "R_squared"] = summary(OLS_list_male[[i]])$adj.r.squared
  OLS_male_fitting_criteria[i, "AIC"] = AIC(OLS_list_male[[i]])
  OLS_male_fitting_criteria[i, "AICc"] = AICc(OLS_list_male[[i]])
  OLS_male_fitting_criteria[i, "R_squared_pred"] = press(OLS_list_male[[i]], as.R2 = TRUE)
  
  # hample
  hampel_male_fitting_criteria[i, "R_squared"] = 1 - sum((hampel_list_male[[i]]$residuals)^2)/sum((hampel_list_male[[i]]$model[,1] - mean(hampel_list_male[[i]]$model[,1]))^2)
  hampel_male_fitting_criteria[i, "AIC"] = AIC(hampel_list_male[[i]])
  hampel_male_fitting_criteria[i, "AICc"] = AICc(hampel_list_male[[i]])
  hampel_male_fitting_criteria[i, "R_squared_pred"] = press(hampel_list_male[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_male_fitting_criteria[i, "R_squared"] = 1 - sum((bisquare_list_male[[i]]$residuals)^2)/sum((bisquare_list_male[[i]]$model[,1] - mean(bisquare_list_male[[i]]$model[,1]))^2)
  bisquare_male_fitting_criteria[i, "AIC"] = AIC(bisquare_list_male[[i]])
  bisquare_male_fitting_criteria[i, "AICc"] = AICc(bisquare_list_male[[i]])
  bisquare_male_fitting_criteria[i, "R_squared_pred"] = press(bisquare_list_male[[i]], as.R2 = TRUE)
}

rownames(OLS_male_fitting_criteria) = rownames(hampel_male_fitting_criteria) = rownames(bisquare_male_fitting_criteria) = c_grid

# c_weight from 0.001 to 0.05

OLS_list_male_refined = hampel_list_male_refined = bisquare_list_male_refined = list()

for(i in 1:20)
{
  OLS_list_male_refined[[i]] = lm(log(diag(male_long_run_spatial_combined_all_refined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$male), base = 10))
  hampel_list_male_refined[[i]] = rlm(log(diag(male_long_run_spatial_combined_all_refined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$male), base = 10), psi = psi.hampel, maxit = 200)
  bisquare_list_male_refined[[i]] = rlm(log(diag(male_long_run_spatial_combined_all_refined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$male), base = 10), psi = psi.bisquare, maxit = 200)
}


OLS_male_fitting_criteria_refined = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_male_fitting_criteria_refined = data.frame(R_squared = NA, AIC = NA, R_squared_pred = NA)
bisquare_male_fitting_criteria_refined = data.frame(R_squared = NA, AIC = NA, R_squared_pred = NA)

for(i in 1:20)
{
  # OLS
  OLS_male_fitting_criteria_refined[i, "R_squared"] = summary(OLS_list_male_refined[[i]])$adj.r.squared
  OLS_male_fitting_criteria_refined[i, "AIC"] = AIC(OLS_list_male_refined[[i]])
  OLS_male_fitting_criteria_refined[i, "AICc"] = AICc(OLS_list_male_refined[[i]])
  OLS_male_fitting_criteria_refined[i, "R_squared_pred"] = press(OLS_list_male_refined[[i]], as.R2 = TRUE)
  
  # hample
  hampel_male_fitting_criteria_refined[i, "R_squared"] = 1 - sum((hampel_list_male_refined[[i]]$residuals)^2)/sum((hampel_list_male_refined[[i]]$model[,1] - mean(hampel_list_male_refined[[i]]$model[,1]))^2)
  hampel_male_fitting_criteria_refined[i, "AIC"] = AIC(hampel_list_male_refined[[i]])
  hampel_male_fitting_criteria_refined[i, "AICc"] = AICc(hampel_list_male_refined[[i]])
  hampel_male_fitting_criteria_refined[i, "R_squared_pred"] = press(hampel_list_male_refined[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_male_fitting_criteria_refined[i, "R_squared"] = 1 - sum((bisquare_list_male_refined[[i]]$residuals)^2)/sum((bisquare_list_male_refined[[i]]$model[,1] - mean(bisquare_list_male_refined[[i]]$model[,1]))^2)
  bisquare_male_fitting_criteria_refined[i, "AIC"] = AIC(bisquare_list_male_refined[[i]])
  bisquare_male_fitting_criteria_refined[i, "AICc"] = AICc(bisquare_list_male_refined[[i]])
  bisquare_male_fitting_criteria_refined[i, "R_squared_pred"] = press(bisquare_list_male_refined[[i]], as.R2 = TRUE)
}

rownames(OLS_male_fitting_criteria_refined) = rownames(hampel_male_fitting_criteria_refined) = rownames(bisquare_male_fitting_criteria_refined) = c_grid_refined

# OLS
OLS_male_fitting_criteria_all = rbind(OLS_male_fitting_criteria_refined, OLS_male_fitting_criteria[2:21,])
plot(x = as.numeric(rownames(OLS_male_fitting_criteria_all)), y = OLS_male_fitting_criteria_all$AICc, xlab = "C_weight", ylab = "AICc", main = "OLS estimation; male")

# hampel
hampel_male_fitting_criteria_all = rbind(hampel_male_fitting_criteria_refined, hampel_male_fitting_criteria[2:21,])
plot(x = as.numeric(rownames(hampel_male_fitting_criteria_all)), y = hampel_male_fitting_criteria_all$AICc, xlab = "C_weight", ylab = "AICc", main = "hampel estimation; male")

# bisquare
bisquare_male_fitting_criteria_all = rbind(bisquare_male_fitting_criteria_refined, bisquare_male_fitting_criteria[2:21,])
plot(x = as.numeric(rownames(bisquare_male_fitting_criteria_all)), y = bisquare_male_fitting_criteria_all$AICc, xlab = "C_weight", ylab = "AICc", main = "bisquare estimation; male")

## total

# c_weight from 0 to 1 by 0.05

OLS_list_total = hampel_list_total = bisquare_list_total = list()

for(i in 1:21)
{
  OLS_list_total[[i]] = lm(log(diag(total_long_run_spatial_combined_all[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$total), base = 10))
  hampel_list_total[[i]] = rlm(log(diag(total_long_run_spatial_combined_all[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$total), base = 10), psi = psi.hampel, maxit = 200)
  bisquare_list_total[[i]] = rlm(log(diag(total_long_run_spatial_combined_all[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$total), base = 10), psi = psi.bisquare, maxit = 200)
}

OLS_total_fitting_criteria = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_total_fitting_criteria = data.frame(R_squared = NA, AIC = NA, R_squared_pred = NA)
bisquare_total_fitting_criteria = data.frame(R_squared = NA, AIC = NA, R_squared_pred = NA)

for(i in 1:21)
{
  # OLS
  OLS_total_fitting_criteria[i, "R_squared"] = summary(OLS_list_total[[i]])$adj.r.squared
  OLS_total_fitting_criteria[i, "AIC"] = AIC(OLS_list_total[[i]])
  OLS_total_fitting_criteria[i, "AICc"] = AICc(OLS_list_total[[i]])
  OLS_total_fitting_criteria[i, "R_squared_pred"] = press(OLS_list_total[[i]], as.R2 = TRUE)
  
  # hample
  hampel_total_fitting_criteria[i, "R_squared"] = 1 - sum((hampel_list_total[[i]]$residuals)^2)/sum((hampel_list_total[[i]]$model[,1] - mean(hampel_list_total[[i]]$model[,1]))^2)
  hampel_total_fitting_criteria[i, "AIC"] = AIC(hampel_list_total[[i]])
  hampel_total_fitting_criteria[i, "AICc"] = AICc(hampel_list_total[[i]])
  hampel_total_fitting_criteria[i, "R_squared_pred"] = press(hampel_list_total[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_total_fitting_criteria[i, "R_squared"] = 1 - sum((bisquare_list_total[[i]]$residuals)^2)/sum((bisquare_list_total[[i]]$model[,1] - mean(bisquare_list_total[[i]]$model[,1]))^2)
  bisquare_total_fitting_criteria[i, "AIC"] = AIC(bisquare_list_total[[i]])
  bisquare_total_fitting_criteria[i, "AICc"] = AICc(bisquare_list_total[[i]])
  bisquare_total_fitting_criteria[i, "R_squared_pred"] = press(bisquare_list_total[[i]], as.R2 = TRUE)
}

rownames(OLS_total_fitting_criteria) = rownames(hampel_total_fitting_criteria) = rownames(bisquare_total_fitting_criteria) = c_grid

# c_weight from 0.001 to 0.05

OLS_list_total_refined = hampel_list_total_refined = bisquare_list_total_refined = list()

for(i in 1:20)
{
  OLS_list_total_refined[[i]] = lm(log(diag(total_long_run_spatial_combined_all_refined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$total), base = 10))
  hampel_list_total_refined[[i]] = rlm(log(diag(total_long_run_spatial_combined_all_refined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$total), base = 10), psi = psi.hampel, maxit = 200)
  bisquare_list_total_refined[[i]] = rlm(log(diag(total_long_run_spatial_combined_all_refined[[i]]$long_run_spatial_variance$BT_FT_fix_C0), base = 10) ~ log(rowMeans(Japan$rate$total), base = 10), psi = psi.bisquare, maxit = 200)
}


OLS_total_fitting_criteria_refined = data.frame(R_squared = NA, AIC = NA, AICc = NA, R_squared_pred = NA) 
hampel_total_fitting_criteria_refined = data.frame(R_squared = NA, AIC = NA, R_squared_pred = NA)
bisquare_total_fitting_criteria_refined = data.frame(R_squared = NA, AIC = NA, R_squared_pred = NA)

for(i in 1:20)
{
  # OLS
  OLS_total_fitting_criteria_refined[i, "R_squared"] = summary(OLS_list_total_refined[[i]])$adj.r.squared
  OLS_total_fitting_criteria_refined[i, "AIC"] = AIC(OLS_list_total_refined[[i]])
  OLS_total_fitting_criteria_refined[i, "AICc"] = AICc(OLS_list_total_refined[[i]])
  OLS_total_fitting_criteria_refined[i, "R_squared_pred"] = press(OLS_list_total_refined[[i]], as.R2 = TRUE)
  
  # hample
  hampel_total_fitting_criteria_refined[i, "R_squared"] = 1 - sum((hampel_list_total_refined[[i]]$residuals)^2)/sum((hampel_list_total_refined[[i]]$model[,1] - mean(hampel_list_total_refined[[i]]$model[,1]))^2)
  hampel_total_fitting_criteria_refined[i, "AIC"] = AIC(hampel_list_total_refined[[i]])
  hampel_total_fitting_criteria_refined[i, "AICc"] = AICc(hampel_list_total_refined[[i]])
  hampel_total_fitting_criteria_refined[i, "R_squared_pred"] = press(hampel_list_total_refined[[i]], as.R2 = TRUE)
  
  # bisquare
  bisquare_total_fitting_criteria_refined[i, "R_squared"] = 1 - sum((bisquare_list_total_refined[[i]]$residuals)^2)/sum((bisquare_list_total_refined[[i]]$model[,1] - mean(bisquare_list_total_refined[[i]]$model[,1]))^2)
  bisquare_total_fitting_criteria_refined[i, "AIC"] = AIC(bisquare_list_total_refined[[i]])
  bisquare_total_fitting_criteria_refined[i, "AICc"] = AICc(bisquare_list_total_refined[[i]])
  bisquare_total_fitting_criteria_refined[i, "R_squared_pred"] = press(bisquare_list_total_refined[[i]], as.R2 = TRUE)
}

rownames(OLS_total_fitting_criteria_refined) = rownames(hampel_total_fitting_criteria_refined) = rownames(bisquare_total_fitting_criteria_refined) = c_grid_refined

# OLS
OLS_total_fitting_criteria_all = rbind(OLS_total_fitting_criteria_refined, OLS_total_fitting_criteria[2:21,])
plot(x = as.numeric(rownames(OLS_total_fitting_criteria_all)), y = OLS_total_fitting_criteria_all$AICc, xlab = "C_weight", ylab = "AICc", main = "OLS estimation; total")

# hampel
hampel_total_fitting_criteria_all = rbind(hampel_total_fitting_criteria_refined, hampel_total_fitting_criteria[2:21,])
plot(x = as.numeric(rownames(hampel_total_fitting_criteria_all)), y = hampel_total_fitting_criteria_all$AICc, xlab = "C_weight", ylab = "AICc", main = "hampel estimation; total")

# bisquare
bisquare_total_fitting_criteria_all = rbind(bisquare_total_fitting_criteria_refined, bisquare_total_fitting_criteria[2:21,])
plot(x = as.numeric(rownames(bisquare_total_fitting_criteria_all)), y = bisquare_total_fitting_criteria_all$AICc, xlab = "C_weight", ylab = "AICc", main = "bisquare estimation; total")



#############
# Make plots
############

## AICc

savepdf("total_AICC", width = 12, height = 10, toplines = 0.8)
par(mar = c(6,4,4,2))
plot(c_grid[2:20], OLS_total_cubic_fitting_combined_criteria$AICc[2:20], type = "l", ylim = c(-91.2, -80.5), ylab = "AICc", xlab = "", xaxt="n", cex.lab= 1.5)
axis(1, at = seq(0.05, 0.95, by = 0.05), las=2)
lines(c_grid[2:20], hampel_total_cubic_fitting_combined_criteria$AICc[2:20], col = 2)
lines(c_grid[2:20], bisquare_total_cubic_fitting_combined_criteria$AICc[2:20], col = 4)
legend("topleft", c("OLS", "Robust (Hampel)", "Robust (Bisquare)"), col = c(1, 2, 4), lty = c(1,1,1), cex = 1)
dev.off()

savepdf("female_AICC", width = 12, height = 10, toplines = 0.8)
par(mar = c(6,4,4,2))
plot(c_grid[2:20], OLS_female_cubic_fitting_combined_criteria$AICc[2:20], type = "l", ylim = c(-81.7, -77.5), ylab = "", xlab = "", xaxt="n", cex.lab= 1.5)
axis(1, at = seq(0.05, 0.95, by = 0.05), las=2)
lines(c_grid[2:20], hampel_female_cubic_fitting_combined_criteria$AICc[2:20], col = 2)
lines(c_grid[2:20], bisquare_female_cubic_fitting_combined_criteria$AICc[2:20], col = 4)
dev.off()

savepdf("male_AICC", width = 12, height = 10, toplines = 0.8)
par(mar = c(6,4,4,2))
plot(c_grid[2:20], OLS_male_cubic_fitting_combined_criteria$AICc[2:20], type = "l", ylim = c(-61.2, -42.1), ylab = "", xlab = "", xaxt="n", cex.lab= 1.5)
axis(1, at = seq(0.05, 0.95, by = 0.05), las=2)
lines(c_grid[2:20], hampel_male_cubic_fitting_combined_criteria$AICc[2:20], col = 2)
lines(c_grid[2:20], bisquare_male_cubic_fitting_combined_criteria$AICc[2:20], col = 4)
dev.off()

## predicted R-squared

savepdf("total_pred_r", width = 12, height = 10, toplines = 0.8)
par(mar = c(6,4,4,2))
plot(c_grid[2:20], OLS_total_cubic_fitting_combined_criteria$R_squared_pred[2:20], type = "l", ylim = c(0.9929, 0.994), ylab = "Predicted R-squared", xlab = "", xaxt="n", cex.lab= 1.5)
axis(1, at = seq(0.05, 0.95, by = 0.05), las=2)
title(xlab = expression(theta), line = 4, cex.lab=2)
lines(c_grid[2:20], hampel_total_cubic_fitting_combined_criteria$R_squared_pred[2:20], col = 2)
lines(c_grid[2:20], bisquare_total_cubic_fitting_combined_criteria$R_squared_pred[2:20], col = 4)
dev.off()

savepdf("female_pred_r", width = 12, height = 10, toplines = 0.8)
par(mar = c(6,4,4,2))
plot(c_grid[2:20], OLS_female_cubic_fitting_combined_criteria$R_squared_pred[2:20], type = "l", ylim = c(0.9931, 0.99385), ylab = "", xlab = "", xaxt="n", cex.lab= 1.5)
axis(1, at = seq(0.05, 0.95, by = 0.05), las=2)
title(xlab = expression(theta), line = 4, cex.lab=2)
lines(c_grid[2:20], hampel_female_cubic_fitting_combined_criteria$R_squared_pred[2:20], col = 2)
lines(c_grid[2:20], bisquare_female_cubic_fitting_combined_criteria$R_squared_pred[2:20], col = 4)
dev.off()

savepdf("male_pred_r", width = 12, height = 10, toplines = 0.8)
par(mar = c(6,4,4,2))
plot(c_grid[2:20], OLS_male_cubic_fitting_combined_criteria$R_squared_pred[2:20], type = "l", ylim = c(0.9852, 0.9906), ylab = "", xlab = "", xaxt="n", cex.lab=1.5)
axis(1, at = seq(0.05, 0.95, by = 0.05), las=2)
title(xlab = expression(theta), line = 4, cex.lab=2)
lines(c_grid[2:20], hampel_male_cubic_fitting_combined_criteria$R_squared_pred[2:20], col = 2)
lines(c_grid[2:20], bisquare_male_cubic_fitting_combined_criteria$R_squared_pred[2:20], col = 4)
dev.off()


###############
# Make a table
###############

long_run_spatial_result_female = long_run_spatial_result_male = long_run_spatial_result_total = matrix(NA, nrow = 6, ncol = 9)
rownames(long_run_spatial_result_female) = rownames(long_run_spatial_result_male) = rownames(long_run_spatial_result_total) = c("a", "b", "c", "d", "AICc", "Predicted R-squared")

## total series
# theta = 0.05
long_run_spatial_result_total[1:4,1] = OLS_list_total_cubic_combined[[2]]$coefficients
long_run_spatial_result_total[5,1] = OLS_total_cubic_fitting_combined_criteria$AICc[2]
long_run_spatial_result_total[6,1] = OLS_total_cubic_fitting_combined_criteria$R_squared_pred[2]

long_run_spatial_result_total[1:4,2] = hampel_list_total_cubic_combined[[2]]$coefficients
long_run_spatial_result_total[5,2] = hampel_total_cubic_fitting_combined_criteria$AICc[2]
long_run_spatial_result_total[6,2] = hampel_total_cubic_fitting_combined_criteria$R_squared_pred[2]

long_run_spatial_result_total[1:4,3] = bisquare_list_total_cubic_combined[[2]]$coefficients
long_run_spatial_result_total[5,3] = bisquare_total_cubic_fitting_combined_criteria$AICc[2]
long_run_spatial_result_total[6,3] = bisquare_total_cubic_fitting_combined_criteria$R_squared_pred[2]

# theta = 0.5
long_run_spatial_result_total[1:4,4] = OLS_list_total_cubic_combined[[11]]$coefficients
long_run_spatial_result_total[5,4] = OLS_total_cubic_fitting_combined_criteria$AICc[11]
long_run_spatial_result_total[6,4] = OLS_total_cubic_fitting_combined_criteria$R_squared_pred[11]

long_run_spatial_result_total[1:4,5] = hampel_list_total_cubic_combined[[11]]$coefficients
long_run_spatial_result_total[5,5] = hampel_total_cubic_fitting_combined_criteria$AICc[11]
long_run_spatial_result_total[6,5] = hampel_total_cubic_fitting_combined_criteria$R_squared_pred[11]

long_run_spatial_result_total[1:4,6] = bisquare_list_total_cubic_combined[[11]]$coefficients
long_run_spatial_result_total[5,6] = bisquare_total_cubic_fitting_combined_criteria$AICc[11]
long_run_spatial_result_total[6,6] = bisquare_total_cubic_fitting_combined_criteria$R_squared_pred[11]

# theta = 0.95
long_run_spatial_result_total[1:4,7] = OLS_list_total_cubic_combined[[20]]$coefficients
long_run_spatial_result_total[5,7] = OLS_total_cubic_fitting_combined_criteria$AICc[20]
long_run_spatial_result_total[6,7] = OLS_total_cubic_fitting_combined_criteria$R_squared_pred[20]

long_run_spatial_result_total[1:4,8] = hampel_list_total_cubic_combined[[20]]$coefficients
long_run_spatial_result_total[5,8] = hampel_total_cubic_fitting_combined_criteria$AICc[20]
long_run_spatial_result_total[6,8] = hampel_total_cubic_fitting_combined_criteria$R_squared_pred[20]

long_run_spatial_result_total[1:4,9] = bisquare_list_total_cubic_combined[[20]]$coefficients
long_run_spatial_result_total[5,9] = bisquare_total_cubic_fitting_combined_criteria$AICc[20]
long_run_spatial_result_total[6,9] = bisquare_total_cubic_fitting_combined_criteria$R_squared_pred[20]

library(xtable)
xtable(long_run_spatial_result_total, digits = 4)

## female series
# theta = 0.05
long_run_spatial_result_female[1:4,1] = OLS_list_female_cubic_combined[[2]]$coefficients
long_run_spatial_result_female[5,1] = OLS_female_cubic_fitting_combined_criteria$AICc[2]
long_run_spatial_result_female[6,1] = OLS_female_cubic_fitting_combined_criteria$R_squared_pred[2]

long_run_spatial_result_female[1:4,2] = hampel_list_female_cubic_combined[[2]]$coefficients
long_run_spatial_result_female[5,2] = hampel_female_cubic_fitting_combined_criteria$AICc[2]
long_run_spatial_result_female[6,2] = hampel_female_cubic_fitting_combined_criteria$R_squared_pred[2]

long_run_spatial_result_female[1:4,3] = bisquare_list_female_cubic_combined[[2]]$coefficients
long_run_spatial_result_female[5,3] = bisquare_female_cubic_fitting_combined_criteria$AICc[2]
long_run_spatial_result_female[6,3] = bisquare_female_cubic_fitting_combined_criteria$R_squared_pred[2]

# theta = 0.5
long_run_spatial_result_female[1:4,4] = OLS_list_female_cubic_combined[[11]]$coefficients
long_run_spatial_result_female[5,4] = OLS_female_cubic_fitting_combined_criteria$AICc[11]
long_run_spatial_result_female[6,4] = OLS_female_cubic_fitting_combined_criteria$R_squared_pred[11]

long_run_spatial_result_female[1:4,5] = hampel_list_female_cubic_combined[[11]]$coefficients
long_run_spatial_result_female[5,5] = hampel_female_cubic_fitting_combined_criteria$AICc[11]
long_run_spatial_result_female[6,5] = hampel_female_cubic_fitting_combined_criteria$R_squared_pred[11]

long_run_spatial_result_female[1:4,6] = bisquare_list_female_cubic_combined[[11]]$coefficients
long_run_spatial_result_female[5,6] = bisquare_female_cubic_fitting_combined_criteria$AICc[11]
long_run_spatial_result_female[6,6] = bisquare_female_cubic_fitting_combined_criteria$R_squared_pred[11]

# theta = 0.95
long_run_spatial_result_female[1:4,7] = OLS_list_female_cubic_combined[[20]]$coefficients
long_run_spatial_result_female[5,7] = OLS_female_cubic_fitting_combined_criteria$AICc[20]
long_run_spatial_result_female[6,7] = OLS_female_cubic_fitting_combined_criteria$R_squared_pred[20]

long_run_spatial_result_female[1:4,8] = hampel_list_female_cubic_combined[[20]]$coefficients
long_run_spatial_result_female[5,8] = hampel_female_cubic_fitting_combined_criteria$AICc[20]
long_run_spatial_result_female[6,8] = hampel_female_cubic_fitting_combined_criteria$R_squared_pred[20]

long_run_spatial_result_female[1:4,9] = bisquare_list_female_cubic_combined[[20]]$coefficients
long_run_spatial_result_female[5,9] = bisquare_female_cubic_fitting_combined_criteria$AICc[20]
long_run_spatial_result_female[6,9] = bisquare_female_cubic_fitting_combined_criteria$R_squared_pred[20]

library(xtable)
xtable(long_run_spatial_result_female, digits = 4)

## male series
# theta = 0.05
long_run_spatial_result_male[1:4,1] = OLS_list_male_cubic_combined[[2]]$coefficients
long_run_spatial_result_male[5,1] = OLS_male_cubic_fitting_combined_criteria$AICc[2]
long_run_spatial_result_male[6,1] = OLS_male_cubic_fitting_combined_criteria$R_squared_pred[2]

long_run_spatial_result_male[1:4,2] = hampel_list_male_cubic_combined[[2]]$coefficients
long_run_spatial_result_male[5,2] = hampel_male_cubic_fitting_combined_criteria$AICc[2]
long_run_spatial_result_male[6,2] = hampel_male_cubic_fitting_combined_criteria$R_squared_pred[2]

long_run_spatial_result_male[1:4,3] = bisquare_list_male_cubic_combined[[2]]$coefficients
long_run_spatial_result_male[5,3] = bisquare_male_cubic_fitting_combined_criteria$AICc[2]
long_run_spatial_result_male[6,3] = bisquare_male_cubic_fitting_combined_criteria$R_squared_pred[2]

# theta = 0.5
long_run_spatial_result_male[1:4,4] = OLS_list_male_cubic_combined[[11]]$coefficients
long_run_spatial_result_male[5,4] = OLS_male_cubic_fitting_combined_criteria$AICc[11]
long_run_spatial_result_male[6,4] = OLS_male_cubic_fitting_combined_criteria$R_squared_pred[11]

long_run_spatial_result_male[1:4,5] = hampel_list_male_cubic_combined[[11]]$coefficients
long_run_spatial_result_male[5,5] = hampel_male_cubic_fitting_combined_criteria$AICc[11]
long_run_spatial_result_male[6,5] = hampel_male_cubic_fitting_combined_criteria$R_squared_pred[11]

long_run_spatial_result_male[1:4,6] = bisquare_list_male_cubic_combined[[11]]$coefficients
long_run_spatial_result_male[5,6] = bisquare_male_cubic_fitting_combined_criteria$AICc[11]
long_run_spatial_result_male[6,6] = bisquare_male_cubic_fitting_combined_criteria$R_squared_pred[11]

# theta = 0.95
long_run_spatial_result_male[1:4,7] = OLS_list_male_cubic_combined[[20]]$coefficients
long_run_spatial_result_male[5,7] = OLS_male_cubic_fitting_combined_criteria$AICc[20]
long_run_spatial_result_male[6,7] = OLS_male_cubic_fitting_combined_criteria$R_squared_pred[20]

long_run_spatial_result_male[1:4,8] = hampel_list_male_cubic_combined[[20]]$coefficients
long_run_spatial_result_male[5,8] = hampel_male_cubic_fitting_combined_criteria$AICc[20]
long_run_spatial_result_male[6,8] = hampel_male_cubic_fitting_combined_criteria$R_squared_pred[20]

long_run_spatial_result_male[1:4,9] = bisquare_list_male_cubic_combined[[20]]$coefficients
long_run_spatial_result_male[5,9] = bisquare_male_cubic_fitting_combined_criteria$AICc[20]
long_run_spatial_result_male[6,9] = bisquare_male_cubic_fitting_combined_criteria$R_squared_pred[20]

library(xtable)
xtable(long_run_spatial_result_male, digits = 4)

### check significance

## OLS
summary(OLS_list_total_cubic_combined[[2]])
summary(OLS_list_female_cubic_combined[[2]])
summary(OLS_list_male_cubic_combined[[2]])

summary(OLS_list_total_cubic_combined[[11]])
summary(OLS_list_female_cubic_combined[[11]])
summary(OLS_list_male_cubic_combined[[11]])

summary(OLS_list_total_cubic_combined[[20]])
summary(OLS_list_female_cubic_combined[[20]])
summary(OLS_list_male_cubic_combined[[20]])

## Hampel
qt(0.975, df = 97, lower.tail = T)
summary(hampel_list_total_cubic_combined[[2]])
summary(hampel_list_female_cubic_combined[[2]])
summary(hampel_list_male_cubic_combined[[2]])

summary(hampel_list_total_cubic_combined[[11]])
summary(hampel_list_female_cubic_combined[[11]])
summary(hampel_list_male_cubic_combined[[11]])

summary(hampel_list_total_cubic_combined[[20]])
summary(hampel_list_female_cubic_combined[[20]])
summary(hampel_list_male_cubic_combined[[20]])

## Bisquare
qt(0.975, df = 97, lower.tail = T)
summary(bisquare_list_total_cubic_combined[[2]])
summary(bisquare_list_female_cubic_combined[[2]])
summary(bisquare_list_male_cubic_combined[[2]])

summary(bisquare_list_total_cubic_combined[[11]])
summary(bisquare_list_female_cubic_combined[[11]])
summary(bisquare_list_male_cubic_combined[[11]])

summary(bisquare_list_total_cubic_combined[[20]])
summary(bisquare_list_female_cubic_combined[[20]])
summary(bisquare_list_male_cubic_combined[[20]])

