
###############################
# Distance-weighted spatial TL
###############################

library(NipponMap)
library(spdep)
library(geosphere)

# binary weights matrix
shp = system.file("shapes/jpn.shp", package = "NipponMap")[1]
Japan_shp = read_sf(shp)
adjacency_weights = poly2nb(Japan_shp, row.names = Japan_shp$name)
weights_matrix = nb2mat(adjacency_weights, style= "B", zero.policy = TRUE)


# distance-related weight matrix

lonlat_mat = Japan_prefecture_city_population[2:48,c("Longitude", "Latitude")]
dist_mat = distm(lonlat_mat)
weights_matrix_2 = 1/dist_mat
diag(weights_matrix_2) = 0

spatial_var = function(y, w_mat)
{
  if(length(y) != nrow(w_mat))
  {
    warnings("uncomfortable matrix")
  }
  
  if(sum(is.na(y))>0)
  {
    ind = which(is.na(y))
    y = na.omit(y)
    w_mat = w_mat[-ind,-ind]
  }
  
  yi_yj = sapply(y, "-", y)
  weighted_cov = w_mat*yi_yj^2
  return(sum(weighted_cov/2)/sum(w_mat))
}


# Distance-weighted spatial covariance

spatial_dist_female_var = spatial_dist_male_var = spatial_dist_total_var = matrix(0, nrow = n_age, ncol = n_year)
for(ij in 1:n_age)
{
  for(it in 1:n_year)
  {
    spatial_dist_female_var[ij,it] = spatial_var(y = female_data[ij,it,2:48], w_mat = weights_matrix_2)
    spatial_dist_male_var[ij,it] = spatial_var(y = male_data[ij,it,2:48], w_mat = weights_matrix_2)
    spatial_dist_total_var[ij,it] = spatial_var(y = total_data[ij,it,2:48], w_mat = weights_matrix_2)
  }
}

# identify outliers

spatial_female_TL_mahalanobis = spatial_male_TL_mahalanobis = spatial_total_TL_mahalanobis = matrix(NA, nrow = 101, ncol = 44)
rownames(spatial_female_TL_mahalanobis) = rownames(spatial_male_TL_mahalanobis) = rownames(spatial_total_TL_mahalanobis) = 0:100
colnames(spatial_female_TL_mahalanobis) = colnames(spatial_male_TL_mahalanobis) = colnames(spatial_total_TL_mahalanobis) = 1975:2018

for(ij in 1:44)
{
  female_select = cbind(log(spatial_dist_female_var[,ij], base = 10), log(spatial_female_mean[,ij], base = 10))
  female_TL_mahalanobis[,ij] = mahalanobis(x = female_select, center = colMeans(female_select), cov = cov(female_select))
  
  male_select = cbind(log(spatial_dist_male_var[,ij], base = 10), log(spatial_male_mean[,ij], base = 10))
  male_TL_mahalanobis[!is.na(spatial_male_mean[,ij]),ij] = mahalanobis(x = na.omit(male_select), center = colMeans(na.omit(male_select)), cov = cov(na.omit(male_select)))
  
  total_select = cbind(log(spatial_dist_total_var[,ij], base = 10), log(spatial_total_mean[,ij], base = 10))
  total_TL_mahalanobis[,ij] = mahalanobis(x = total_select, center = colMeans(total_select), cov = cov(total_select))
}

spatial_female_TL_outliers = spatial_male_TL_outliers = spatial_total_TL_outliers = 
  spatial_female_TL_mahalanobis_outliers = spatial_male_TL_mahalanobis_outliers = spatial_total_TL_mahalanobis_outliers = 
  spatial_female_TL_rstudent_outliers = spatial_male_TL_rstudent_outliers = spatial_total_TL_rstudent_outliers = list()


for(ij in 1:44)
{
  female_select = data.frame(log_var = log10(spatial_female_mean[,ij]), log_mean = log10(spatial_dist_female_var[,ij]))
  spatial_female_TL_mahalanobis_outliers[[ij]] = female_select[spatial_female_TL_mahalanobis[,ij] > qchisq(p = 0.95, df = 2), , drop = FALSE]
  spatial_female_TL_rstudent_outliers[[ij]] = female_select[abs(rstudent(lm(log10(spatial_dist_female_var[,ij]) ~ log10(spatial_female_mean[,ij])))) > 2, , drop = FALSE]
  outlier_age = as.numeric(intersect(row.names(spatial_female_TL_mahalanobis_outliers[[ij]]), row.names(spatial_female_TL_rstudent_outliers[[ij]])))
  spatial_female_TL_outliers[[ij]] = female_select[outlier_age+1,]
  
  male_select = data.frame(log_var = log10(spatial_male_mean[,ij]), log_mean = log10(spatial_dist_male_var[,ij]))
  spatial_male_TL_mahalanobis_outliers[[ij]] = male_select[spatial_male_TL_mahalanobis[,ij] > qchisq(p = 0.95, df = 2), , drop = FALSE]
  spatial_male_TL_rstudent_outliers[[ij]] = male_select[abs(rstudent(lm(log10(spatial_dist_male_var[,ij]) ~ log10(spatial_male_mean[,ij])))) > 2, , drop = FALSE]
  outlier_age = as.numeric(intersect(row.names(spatial_male_TL_mahalanobis_outliers[[ij]]), row.names(spatial_male_TL_rstudent_outliers[[ij]])))
  spatial_male_TL_outliers[[ij]] = male_select[outlier_age+1,]
  
  total_select = data.frame(log_var = log10(spatial_total_mean[,ij]), log_mean = log10(spatial_dist_total_var[,ij]))
  spatial_total_TL_mahalanobis_outliers[[ij]] = total_select[spatial_total_TL_mahalanobis[,ij] > qchisq(p = 0.95, df = 2), , drop = FALSE]
  spatial_total_TL_rstudent_outliers[[ij]] = total_select[abs(rstudent(lm(log10(spatial_dist_total_var[,ij]) ~ log10(spatial_total_mean[,ij])))) > 2, , drop = FALSE]
  outlier_age = as.numeric(intersect(row.names(spatial_total_TL_mahalanobis_outliers[[ij]]), row.names(spatial_total_TL_rstudent_outliers[[ij]])))
  spatial_total_TL_outliers[[ij]] = total_select[outlier_age+1,]
}

names(spatial_female_TL_outliers) = names(spatial_male_TL_outliers) = names(spatial_total_TL_outliers) = 
  names(spatial_female_TL_mahalanobis_outliers) = names(spatial_male_TL_mahalanobis_outliers) = names(spatial_total_TL_mahalanobis_outliers) =
  names(spatial_female_TL_rstudent_outliers) = names(spatial_male_TL_rstudent_outliers) = names(spatial_total_TL_rstudent_outliers) = 1975:2018

sum(unlist(lapply(spatial_female_TL_mahalanobis_outliers, nrow)))
sum(unlist(lapply(spatial_female_TL_outliers, nrow)))

unlist(lapply(spatial_female_TL_outliers, nrow))
unlist(lapply(spatial_male_TL_outliers, nrow))
unlist(lapply(spatial_total_TL_outliers, nrow))

#########################################
# Taylor's law slope parameter estimates
#########################################
library(MASS)
library(tidyverse)

find_p_value = function(t_values, df)
{
  n = length(t_values)
  out = rep(0, n)
  for(i in 1:n)
  {
    if(t_values[i] < 0)
    {
      out[i] = pt(q = t_values[i], df = df, lower.tail = TRUE)
    } else {
      out[i] = pt(q = t_values[i], df = df, lower.tail = FALSE)
    }
  }
  return(out)
}

# female

# 1: benchmark of the classic OLS TL

spatial_female_classic = spatial_female_classic_hampel = spatial_female_classic_bisquare = matrix(NA, nrow = 2, ncol = 44)
spatial_female_classic_p = spatial_female_classic_hampel_p = spatial_female_classic_bisquare_p = matrix(NA, nrow = 2, ncol = 44)
spatial_female_classic_R_squared_prefecture = spatial_female_classic_hampel_R_squared_prefecture = spatial_female_classic_bisquare_R_squared_prefecture = vector("numeric", 44)

for(ij in 1:44)
{
  y = log(spatial_female_var[,ij], base = 10)
  x = log(spatial_female_mean[,ij], base = 10)
  
  
  OLS_model = lm(y ~ x)
  spatial_female_classic[,ij] = OLS_model$coefficients
  spatial_female_classic_p[,ij] = summary(OLS_model)$coef[,4]
  spatial_female_classic_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x, psi = psi.hampel, maxit = 200)
  spatial_female_classic_hampel[,ij] = hampel_model$coefficients
  spatial_female_classic_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  spatial_female_classic_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
  spatial_female_classic_bisquare[,ij] = bisquare_model$coefficients
  spatial_female_classic_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  spatial_female_classic_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# 2: benchmark of quadratic OLS TL

spatial_female_quadratic = spatial_female_quadratic_hampel = spatial_female_quadratic_bisquare = matrix(NA, nrow = 3, ncol = 44)
spatial_female_quadratic_p = spatial_female_quadratic_hampel_p = spatial_female_quadratic_bisquare_p = matrix(NA, nrow = 3, ncol = 44)
spatial_female_quadratic_R_squared_prefecture = spatial_female_quadratic_hampel_R_squared_prefecture = spatial_female_quadratic_bisquare_R_squared_prefecture = vector("numeric", 44)

for(ij in 1:44)
{
  y = log(spatial_female_var[,ij], base = 10)
  x = log(spatial_female_mean[,ij], base = 10)
  
  OLS_model = lm(y ~ x + I(x^2))
  spatial_female_quadratic[,ij] = OLS_model$coefficients
  spatial_female_quadratic_p[,ij] = summary(OLS_model)$coef[,4]
  spatial_female_quadratic_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
  spatial_female_quadratic_hampel[,ij] = hampel_model$coefficients
  spatial_female_quadratic_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  spatial_female_quadratic_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
  spatial_female_quadratic_bisquare[,ij] = bisquare_model$coefficients
  spatial_female_quadratic_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  spatial_female_quadratic_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# 3: cubic TL with OLS estimation

spatial_female_TL_coef = spatial_female_TL_hampel_coef = spatial_female_TL_bisquare_coef = matrix(NA, nrow = 4, ncol = 44)
spatial_female_TL_p = spatial_female_TL_hampel_p = spatial_female_TL_bisquare_p = matrix(NA, nrow = 4, ncol = 44)
spatial_female_TL_R_squared_prefecture = spatial_female_TL_hampel_R_squared_prefecture = spatial_female_TL_bisquare_R_squared_prefecture = vector("numeric", 44)

for(ij in 1:44)
{
  y = log(spatial_female_var[,ij], base = 10)
  x = log(spatial_female_mean[,ij], base = 10)
  
  OLS_model = lm(y ~ x + I(x^2) + I(x^3))
  spatial_female_TL_coef[,ij] = OLS_model$coefficients
  spatial_female_TL_p[,ij] = summary(OLS_model)$coef[,4]
  spatial_female_TL_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)
  spatial_female_TL_hampel_coef[,ij] = hampel_model$coefficients
  spatial_female_TL_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  spatial_female_TL_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)
  spatial_female_TL_bisquare_coef[,ij] = bisquare_model$coefficients
  spatial_female_TL_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  spatial_female_TL_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# male

# 1: benchmark of the classic OLS TL

spatial_male_classic = spatial_male_classic_hampel = spatial_male_classic_bisquare = matrix(NA, nrow = 2, ncol = 44)
spatial_male_classic_p = spatial_male_classic_hampel_p = spatial_male_classic_bisquare_p = matrix(NA, nrow = 2, ncol = 44)
spatial_male_classic_R_squared_prefecture = spatial_male_classic_hampel_R_squared_prefecture = spatial_male_classic_bisquare_R_squared_prefecture = vector("numeric", 44)

for(ij in 1:44)
{
  y = log(spatial_male_var[,ij], base = 10)
  x = log(spatial_male_mean[,ij], base = 10)
  
  
  OLS_model = lm(y ~ x)
  spatial_male_classic[,ij] = OLS_model$coefficients
  spatial_male_classic_p[,ij] = summary(OLS_model)$coef[,4]
  spatial_male_classic_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x, psi = psi.hampel, maxit = 200)
  spatial_male_classic_hampel[,ij] = hampel_model$coefficients
  spatial_male_classic_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  spatial_male_classic_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
  spatial_male_classic_bisquare[,ij] = bisquare_model$coefficients
  spatial_male_classic_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  spatial_male_classic_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# 2: benchmark of quadratic OLS TL

spatial_male_quadratic = spatial_male_quadratic_hampel = spatial_male_quadratic_bisquare = matrix(NA, nrow = 3, ncol = 44)
spatial_male_quadratic_p = spatial_male_quadratic_hampel_p = spatial_male_quadratic_bisquare_p = matrix(NA, nrow = 3, ncol = 44)
spatial_male_quadratic_R_squared_prefecture = spatial_male_quadratic_hampel_R_squared_prefecture = spatial_male_quadratic_bisquare_R_squared_prefecture = vector("numeric", 44)

for(ij in 1:44)
{
  y = log(spatial_male_var[,ij], base = 10)
  x = log(spatial_male_mean[,ij], base = 10)
  
  OLS_model = lm(y ~ x + I(x^2))
  spatial_male_quadratic[,ij] = OLS_model$coefficients
  spatial_male_quadratic_p[,ij] = summary(OLS_model)$coef[,4]
  spatial_male_quadratic_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
  spatial_male_quadratic_hampel[,ij] = hampel_model$coefficients
  spatial_male_quadratic_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  spatial_male_quadratic_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
  spatial_male_quadratic_bisquare[,ij] = bisquare_model$coefficients
  spatial_male_quadratic_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  spatial_male_quadratic_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# 3: cubic TL with OLS estimation

spatial_male_TL_coef = spatial_male_TL_hampel_coef = spatial_male_TL_bisquare_coef = matrix(NA, nrow = 4, ncol = 44)
spatial_male_TL_p = spatial_male_TL_hampel_p = spatial_male_TL_bisquare_p = matrix(NA, nrow = 4, ncol = 44)
spatial_male_TL_R_squared_prefecture = spatial_male_TL_hampel_R_squared_prefecture = spatial_male_TL_bisquare_R_squared_prefecture = vector("numeric", 44)

for(ij in 1:44)
{
  y = log(spatial_male_var[,ij], base = 10)
  x = log(spatial_male_mean[,ij], base = 10)
  
  OLS_model = lm(y ~ x + I(x^2) + I(x^3))
  spatial_male_TL_coef[,ij] = OLS_model$coefficients
  spatial_male_TL_p[,ij] = summary(OLS_model)$coef[,4]
  spatial_male_TL_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)
  spatial_male_TL_hampel_coef[,ij] = hampel_model$coefficients
  spatial_male_TL_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  spatial_male_TL_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)
  spatial_male_TL_bisquare_coef[,ij] = bisquare_model$coefficients
  spatial_male_TL_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  spatial_male_TL_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# total

# 1: benchmark of the classic OLS TL

spatial_total_classic = spatial_total_classic_hampel = spatial_total_classic_bisquare = matrix(NA, nrow = 2, ncol = 44)
spatial_total_classic_p = spatial_total_classic_hampel_p = spatial_total_classic_bisquare_p = matrix(NA, nrow = 2, ncol = 44)
spatial_total_classic_R_squared_prefecture = spatial_total_classic_hampel_R_squared_prefecture = spatial_total_classic_bisquare_R_squared_prefecture = vector("numeric", 44)

for(ij in 1:44)
{
  y = log(spatial_total_var[,ij], base = 10)
  x = log(spatial_total_mean[,ij], base = 10)
  
  
  OLS_model = lm(y ~ x)
  spatial_total_classic[,ij] = OLS_model$coefficients
  spatial_total_classic_p[,ij] = summary(OLS_model)$coef[,4]
  spatial_total_classic_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x, psi = psi.hampel, maxit = 200)
  spatial_total_classic_hampel[,ij] = hampel_model$coefficients
  spatial_total_classic_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  spatial_total_classic_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
  spatial_total_classic_bisquare[,ij] = bisquare_model$coefficients
  spatial_total_classic_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  spatial_total_classic_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# 2: benchmark of quadratic OLS TL

spatial_total_quadratic = spatial_total_quadratic_hampel = spatial_total_quadratic_bisquare = matrix(NA, nrow = 3, ncol = 44)
spatial_total_quadratic_p = spatial_total_quadratic_hampel_p = spatial_total_quadratic_bisquare_p = matrix(NA, nrow = 3, ncol = 44)
spatial_total_quadratic_R_squared_prefecture = spatial_total_quadratic_hampel_R_squared_prefecture = spatial_total_quadratic_bisquare_R_squared_prefecture = vector("numeric", 44)

for(ij in 1:44)
{
  y = log(spatial_total_var[,ij], base = 10)
  x = log(spatial_total_mean[,ij], base = 10)
  
  OLS_model = lm(y ~ x + I(x^2))
  spatial_total_quadratic[,ij] = OLS_model$coefficients
  spatial_total_quadratic_p[,ij] = summary(OLS_model)$coef[,4]
  spatial_total_quadratic_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
  spatial_total_quadratic_hampel[,ij] = hampel_model$coefficients
  spatial_total_quadratic_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  spatial_total_quadratic_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
  spatial_total_quadratic_bisquare[,ij] = bisquare_model$coefficients
  spatial_total_quadratic_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  spatial_total_quadratic_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# 3: cubic TL with OLS estimation

spatial_total_TL_coef = spatial_total_TL_hampel_coef = spatial_total_TL_bisquare_coef = matrix(NA, nrow = 4, ncol = 44)
spatial_total_TL_p = spatial_total_TL_hampel_p = spatial_total_TL_bisquare_p = matrix(NA, nrow = 4, ncol = 44)
spatial_total_TL_R_squared_prefecture = spatial_total_TL_hampel_R_squared_prefecture = spatial_total_TL_bisquare_R_squared_prefecture = vector("numeric", 44)

for(ij in 1:44)
{
  y = log(spatial_total_var[,ij], base = 10)
  x = log(spatial_total_mean[,ij], base = 10)
  
  OLS_model = lm(y ~ x + I(x^2) + I(x^3))
  spatial_total_TL_coef[,ij] = OLS_model$coefficients
  spatial_total_TL_p[,ij] = summary(OLS_model)$coef[,4]
  spatial_total_TL_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)
  spatial_total_TL_hampel_coef[,ij] = hampel_model$coefficients
  spatial_total_TL_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  spatial_total_TL_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)
  spatial_total_TL_bisquare_coef[,ij] = bisquare_model$coefficients
  spatial_total_TL_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  spatial_total_TL_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

