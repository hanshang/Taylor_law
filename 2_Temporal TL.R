##############
# Temporal TL
##############

n_age = 101
age_all = 0:100
year = 1975:2018
n_year = length(year)
n_prefecture = 47

total_data = male_data = female_data = total_data_smooth = male_data_smooth = female_data_smooth = array(NA, dim = c(n_age, n_year, n_prefecture+1))
for(ik in 1:(n_prefecture + 1))
{
    female_data[,,ik] = get(state[ik])$rate$female
    male_data[,,ik]   = get(state[ik])$rate$male
    total_data[,,ik]  = get(state[ik])$rate$total
}

##################################
# plot temporal mean and variance
##################################

library(MASS)
library(qpcR)
library(tidyverse)

# female

temporal_female_mean = apply(female_data[,,2:48], c(1, 3), mean)
temporal_female_var  = apply(female_data[,,2:48], c(1, 3), var)

rownames(temporal_female_mean) = rownames(temporal_female_var) = 0:100
colnames(temporal_female_mean) = colnames(temporal_female_var) = state[2:48]

plot(fts(0:100, log(temporal_female_mean, base = 10)), xlab = "Age", ylab = "Female mortality")
plot(fts(0:100, log(temporal_female_var,  base = 10)), xlab = "Age", ylab = "Female mortality")

# male

temporal_male_mean = apply(male_data[,,2:48], c(1, 3), mean)
temporal_male_var  = apply(male_data[,,2:48], c(1, 3), var)

rownames(temporal_male_mean) = rownames(temporal_male_var) = 0:100
colnames(temporal_male_mean) = colnames(temporal_male_var) = state[2:48]

plot(fts(0:100, log(temporal_male_mean, base = 10)), xlab = "Age", ylab = "Male mortality")
plot(fts(0:100, log(temporal_male_var,  base = 10)), xlab = "Age", ylab = "Male mortality")

# total

temporal_total_mean = apply(total_data[,,2:48], c(1, 3), mean)
temporal_total_var  = apply(total_data[,,2:48], c(1, 3), var)

rownames(temporal_total_mean) = rownames(temporal_total_var) = 0:100
colnames(temporal_total_mean) = colnames(temporal_total_var) = state[2:48]

plot(fts(0:100, log(temporal_total_mean, base = 10)), xlab = "Age", ylab = "Total mortality")
plot(fts(0:100, log(temporal_total_var,  base = 10)), xlab = "Age", ylab = "Total mortality")


# identify outliers

female_TL_mahalanobis = male_TL_mahalanobis = total_TL_mahalanobis = matrix(NA, nrow = 101, ncol = 47)
rownames(female_TL_mahalanobis) = rownames(male_TL_mahalanobis) = rownames(total_TL_mahalanobis) = 0:100
colnames(female_TL_mahalanobis) = colnames(male_TL_mahalanobis) = colnames(total_TL_mahalanobis) = state[-1]
  
for(ij in 1:47)
{
  female_select = cbind(log(temporal_female_var[,ij], base = 10), log(temporal_female_mean[,ij], base = 10))
  female_TL_mahalanobis[,ij] = mahalanobis(x = female_select, center = colMeans(female_select), cov = cov(female_select))
  
  male_select = cbind(log(temporal_male_var[,ij], base = 10), log(temporal_male_mean[,ij], base = 10))
  male_TL_mahalanobis[!is.na(temporal_male_mean[,ij]),ij] = mahalanobis(x = na.omit(male_select), center = colMeans(na.omit(male_select)), cov = cov(na.omit(male_select)))
  
  total_select = cbind(log(temporal_total_var[,ij], base = 10), log(temporal_total_mean[,ij], base = 10))
  total_TL_mahalanobis[,ij] = mahalanobis(x = total_select, center = colMeans(total_select), cov = cov(total_select))
}

female_TL_outliers = male_TL_outliers = total_TL_outliers = 
  female_TL_mahalanobis_outliers = male_TL_mahalanobis_outliers = total_TL_mahalanobis_outliers = 
  female_TL_rstudent_outliers = male_TL_rstudent_outliers = total_TL_rstudent_outliers = list()

for(ij in 1:47)
{
  female_select = data.frame(log_var = log10(temporal_female_mean[,ij]), log_mean = log10(temporal_female_var[,ij]))
  female_TL_mahalanobis_outliers[[ij]] = female_select[female_TL_mahalanobis[,ij] > qchisq(p = 0.95, df = 2), , drop = FALSE]
  female_TL_rstudent_outliers[[ij]] = female_select[abs(rstudent(lm(log10(temporal_female_var[,ij]) ~ log10(temporal_female_mean[,ij])))) > 2, , drop = FALSE]
  outlier_age = as.numeric(intersect(row.names(female_TL_mahalanobis_outliers[[ij]]), row.names(female_TL_rstudent_outliers[[ij]])))
  female_TL_outliers[[ij]] = female_select[outlier_age+1,]
  
  male_select = data.frame(log_var = log10(temporal_male_mean[,ij]), log_mean = log10(temporal_male_var[,ij]))
  male_TL_mahalanobis_outliers[[ij]] = male_select[male_TL_mahalanobis[,ij] > qchisq(p = 0.95, df = 2), , drop = FALSE]
  male_TL_rstudent_outliers[[ij]] = male_select[abs(rstudent(lm(log10(temporal_male_var[,ij]) ~ log10(temporal_male_mean[,ij])))) > 2, , drop = FALSE]
  outlier_age = as.numeric(intersect(row.names(male_TL_mahalanobis_outliers[[ij]]), row.names(male_TL_rstudent_outliers[[ij]])))
  male_TL_outliers[[ij]] = male_select[outlier_age+1,]
  
  total_select = data.frame(log_var = log10(temporal_total_mean[,ij]), log_mean = log10(temporal_total_var[,ij]))
  total_TL_mahalanobis_outliers[[ij]] = total_select[total_TL_mahalanobis[,ij] > qchisq(p = 0.95, df = 2), , drop = FALSE]
  total_TL_rstudent_outliers[[ij]] = total_select[abs(rstudent(lm(log10(temporal_total_var[,ij]) ~ log10(temporal_total_mean[,ij])))) > 2, , drop = FALSE]
  outlier_age = as.numeric(intersect(row.names(total_TL_mahalanobis_outliers[[ij]]), row.names(total_TL_rstudent_outliers[[ij]])))
  total_TL_outliers[[ij]] = total_select[outlier_age+1,]
}

names(female_TL_outliers) = names(male_TL_outliers) = names(total_TL_outliers) = 
  names(female_TL_mahalanobis_outliers) = names(male_TL_mahalanobis_outliers) = names(total_TL_mahalanobis_outliers) =
  names(female_TL_rstudent_outliers) = names(male_TL_rstudent_outliers) = names(total_TL_rstudent_outliers) = state[-1]

###########################
# Estimate TL coefficients
###########################


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

# benchmark of the classic OLS TL

female_classic = female_classic_hampel = female_classic_bisquare = matrix(NA, nrow = 2, ncol = 47)
female_classic_p = female_classic_hampel_p = female_classic_bisquare_p = matrix(NA, nrow = 2, ncol = 47)
female_classic_R_squared_prefecture = female_classic_hampel_R_squared_prefecture = female_classic_bisquare_R_squared_prefecture = vector("numeric", 47)

for(ij in 1:47)
{
  if(length(female_TL_outliers[[ij]]) > 0)
  {
    age_select = setdiff(age_all, as.numeric(rownames(female_TL_outliers[[ij]])))
    y = log(temporal_female_var[age_all %in% age_select,ij], base = 10)
    x = log(temporal_female_mean[age_all %in% age_select,ij], base = 10)
  } else {
    y = log(temporal_female_var[,ij], base = 10)
    x = log(temporal_female_mean[,ij], base = 10)
  }
  
  OLS_model = lm(y ~ x)
  female_classic[,ij] = OLS_model$coefficients
  female_classic_p[,ij] = summary(OLS_model)$coef[,4]
  female_classic_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x, psi = psi.hampel, maxit = 200)
  female_classic_hampel[,ij] = hampel_model$coefficients
  female_classic_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  female_classic_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
  female_classic_bisquare[,ij] = bisquare_model$coefficients
  female_classic_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  female_classic_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# benchmark of quadratic OLS TL

female_quadratic = female_quadratic_hampel = female_quadratic_bisquare = matrix(NA, nrow = 3, ncol = 47)
female_quadratic_p = female_quadratic_hampel_p = female_quadratic_bisquare_p = matrix(NA, nrow = 3, ncol = 47)
female_quadratic_R_squared_prefecture = female_quadratic_hampel_R_squared_prefecture = female_quadratic_bisquare_R_squared_prefecture = vector("numeric", 47)

for(ij in 1:47)
{
  if(length(female_TL_outliers[[ij]]) > 0)
  {
    age_select = setdiff(age_all, as.numeric(rownames(female_TL_outliers[[ij]])))
    y = log(temporal_female_var[age_all %in% age_select,ij], base = 10)
    x = log(temporal_female_mean[age_all %in% age_select,ij], base = 10)
  } else {
    y = log(temporal_female_var[,ij], base = 10)
    x = log(temporal_female_mean[,ij], base = 10)
  }
  
  OLS_model = lm(y ~ x + I(x^2))
  female_quadratic[,ij] = OLS_model$coefficients
  female_quadratic_p[,ij] = summary(OLS_model)$coef[,4]
  female_quadratic_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
  female_quadratic_hampel[,ij] = hampel_model$coefficients
  female_quadratic_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  female_quadratic_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
  female_quadratic_bisquare[,ij] = bisquare_model$coefficients
  female_quadratic_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  female_quadratic_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# cubic TL with OLS estimation

female_TL_coef = female_TL_hampel_coef = female_TL_bisquare_coef = matrix(NA, nrow = 4, ncol = 47)
female_TL_p = female_TL_hampel_p = female_TL_bisquare_p = matrix(NA, nrow = 4, ncol = 47)

female_TL_R_squared_prefecture = female_TL_hampel_R_squared_prefecture = female_TL_bisquare_R_squared_prefecture = vector("numeric", 47)

for(ij in 1:47)
{
  if(length(female_TL_outliers[[ij]]) > 0)
  {
    age_select = setdiff(age_all, as.numeric(rownames(female_TL_outliers[[ij]])))
    y = log(temporal_female_var[age_all %in% age_select,ij], base = 10)
    x = log(temporal_female_mean[age_all %in% age_select,ij], base = 10)
  } else {
    y = log(temporal_female_var[,ij], base = 10)
    x = log(temporal_female_mean[,ij], base = 10)
  }
  
  OLS_model = lm(y ~ x + I(x^2) + I(x^3))
  female_TL_coef[,ij] = OLS_model$coefficients
  female_TL_p[,ij] = summary(OLS_model)$coef[,4]
  female_TL_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared

  hampel_model = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)
  female_TL_hampel_coef[,ij] = hampel_model$coefficients
  female_TL_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  female_TL_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)
  female_TL_bisquare_coef[,ij] = bisquare_model$coefficients
  female_TL_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  female_TL_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# male

# benchmark of the classic OLS TL

male_classic = male_classic_hampel = male_classic_bisquare = matrix(NA, nrow = 2, ncol = 47)
male_classic_p = male_classic_hampel_p = male_classic_bisquare_p = matrix(NA, nrow = 2, ncol = 47)
male_classic_R_squared_prefecture = male_classic_hampel_R_squared_prefecture = male_classic_bisquare_R_squared_prefecture = vector("numeric", 47)

for(ij in 1:47)
{
  if(length(male_TL_outliers[[ij]]) > 0)
  {
    age_select = setdiff(age_all, as.numeric(rownames(male_TL_outliers[[ij]])))
    y = log(temporal_male_var[age_all %in% age_select,ij], base = 10)
    x = log(temporal_male_mean[age_all %in% age_select,ij], base = 10)
  } else {
    y = log(temporal_male_var[,ij], base = 10)
    x = log(temporal_male_mean[,ij], base = 10)
  }
  
  OLS_model = lm(y ~ x)
  male_classic[,ij] = OLS_model$coefficients
  male_classic_p[,ij] = summary(OLS_model)$coef[,4]
  male_classic_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x, psi = psi.hampel, maxit = 200)
  male_classic_hampel[,ij] = hampel_model$coefficients
  male_classic_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  male_classic_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
  male_classic_bisquare[,ij] = bisquare_model$coefficients
  male_classic_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  male_classic_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# benchmark of quadratic OLS TL

male_quadratic = male_quadratic_hampel = male_quadratic_bisquare = matrix(NA, nrow = 3, ncol = 47)
male_quadratic_p = male_quadratic_hampel_p = male_quadratic_bisquare_p = matrix(NA, nrow = 3, ncol = 47)
male_quadratic_R_squared_prefecture = male_quadratic_hampel_R_squared_prefecture = male_quadratic_bisquare_R_squared_prefecture = vector("numeric", 47)

for(ij in 1:47)
{
  if(length(male_TL_outliers[[ij]]) > 0)
  {
    age_select = setdiff(age_all, as.numeric(rownames(male_TL_outliers[[ij]])))
    y = log(temporal_male_var[age_all %in% age_select,ij], base = 10)
    x = log(temporal_male_mean[age_all %in% age_select,ij], base = 10)
  } else {
    y = log(temporal_male_var[,ij], base = 10)
    x = log(temporal_male_mean[,ij], base = 10)
  }
  
  OLS_model = lm(y ~ x + I(x^2))
  male_quadratic[,ij] = OLS_model$coefficients
  male_quadratic_p[,ij] = summary(OLS_model)$coef[,4]
  male_quadratic_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
  male_quadratic_hampel[,ij] = hampel_model$coefficients
  male_quadratic_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  male_quadratic_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
  male_quadratic_bisquare[,ij] = bisquare_model$coefficients
  male_quadratic_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  male_quadratic_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}


# cubic TL with OLS estimation

male_TL_coef = male_TL_hampel_coef = male_TL_bisquare_coef = matrix(NA, nrow = 4, ncol = 47)
male_TL_p = male_TL_hampel_p = male_TL_bisquare_p = matrix(NA, nrow = 4, ncol = 47)

male_TL_R_squared_prefecture = male_TL_hampel_R_squared_prefecture = male_TL_bisquare_R_squared_prefecture = vector("numeric", 47)

for(ij in 1:47)
{
  if(length(male_TL_outliers[[ij]]) > 0)
  {
    age_select = setdiff(age_all, as.numeric(rownames(male_TL_outliers[[ij]])))
    y = log(temporal_male_var[age_all %in% age_select,ij], base = 10)
    x = log(temporal_male_mean[age_all %in% age_select,ij], base = 10)
  } else {
    y = log(temporal_male_var[,ij], base = 10)
    x = log(temporal_male_mean[,ij], base = 10)
  }
  
  OLS_model = lm(y ~ x + I(x^2) + I(x^3))
  male_TL_coef[,ij] = OLS_model$coefficients
  male_TL_p[,ij] = summary(OLS_model)$coef[,4]
  male_TL_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)
  male_TL_hampel_coef[,ij] = hampel_model$coefficients
  male_TL_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  male_TL_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)
  male_TL_bisquare_coef[,ij] = bisquare_model$coefficients
  male_TL_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  male_TL_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# total

# benchmark of the classic OLS TL

total_classic = total_classic_hampel = total_classic_bisquare = matrix(NA, nrow = 2, ncol = 47)
total_classic_p = total_classic_hampel_p = total_classic_bisquare_p = matrix(NA, nrow = 2, ncol = 47)
total_classic_R_squared_prefecture = total_classic_hampel_R_squared_prefecture = total_classic_bisquare_R_squared_prefecture = vector("numeric", 47)

for(ij in 1:47)
{
  if(length(total_TL_outliers[[ij]]) > 0)
  {
    age_select = setdiff(age_all, as.numeric(rownames(total_TL_outliers[[ij]])))
    y = log(temporal_total_var[age_all %in% age_select,ij], base = 10)
    x = log(temporal_total_mean[age_all %in% age_select,ij], base = 10)
  } else {
    y = log(temporal_total_var[,ij], base = 10)
    x = log(temporal_total_mean[,ij], base = 10)
  }
  
  OLS_model = lm(y ~ x)
  total_classic[,ij] = OLS_model$coefficients
  total_classic_p[,ij] = summary(OLS_model)$coef[,4]
  total_classic_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x, psi = psi.hampel, maxit = 200)
  total_classic_hampel[,ij] = hampel_model$coefficients
  total_classic_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  total_classic_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x, psi = psi.bisquare, maxit = 200)
  total_classic_bisquare[,ij] = bisquare_model$coefficients
  total_classic_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  total_classic_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# benchmark of quadratic OLS TL

total_quadratic = total_quadratic_hampel = total_quadratic_bisquare = matrix(NA, nrow = 3, ncol = 47)
total_quadratic_p = total_quadratic_hampel_p = total_quadratic_bisquare_p = matrix(NA, nrow = 3, ncol = 47)
total_quadratic_R_squared_prefecture = total_quadratic_hampel_R_squared_prefecture = total_quadratic_bisquare_R_squared_prefecture = vector("numeric", 47)

for(ij in 1:47)
{
  if(length(total_TL_outliers[[ij]]) > 0)
  {
    age_select = setdiff(age_all, as.numeric(rownames(total_TL_outliers[[ij]])))
    y = log(temporal_total_var[age_all %in% age_select,ij], base = 10)
    x = log(temporal_total_mean[age_all %in% age_select,ij], base = 10)
  } else {
    y = log(temporal_total_var[,ij], base = 10)
    x = log(temporal_total_mean[,ij], base = 10)
  }
  
  OLS_model = lm(y ~ x + I(x^2))
  total_quadratic[,ij] = OLS_model$coefficients
  total_quadratic_p[,ij] = summary(OLS_model)$coef[,4]
  total_quadratic_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x + I(x^2), psi = psi.hampel, maxit = 200)
  total_quadratic_hampel[,ij] = hampel_model$coefficients
  total_quadratic_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  total_quadratic_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x + I(x^2), psi = psi.bisquare, maxit = 200)
  total_quadratic_bisquare[,ij] = bisquare_model$coefficients
  total_quadratic_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  total_quadratic_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

# cubic TL with OLS estimation

total_TL_coef = total_TL_hampel_coef = total_TL_bisquare_coef = matrix(NA, nrow = 4, ncol = 47)
total_TL_p = total_TL_hampel_p = total_TL_bisquare_p = matrix(NA, nrow = 4, ncol = 47)

total_TL_R_squared_prefecture = total_TL_hampel_R_squared_prefecture = total_TL_bisquare_R_squared_prefecture = vector("numeric", 47)

for(ij in 1:47)
{
  if(length(total_TL_outliers[[ij]]) > 0)
  {
    age_select = setdiff(age_all, as.numeric(rownames(total_TL_outliers[[ij]])))
    y = log(temporal_total_var[age_all %in% age_select,ij], base = 10)
    x = log(temporal_total_mean[age_all %in% age_select,ij], base = 10)
  } else {
    y = log(temporal_total_var[,ij], base = 10)
    x = log(temporal_total_mean[,ij], base = 10)
  }
  
  OLS_model = lm(y ~ x + I(x^2) + I(x^3))
  total_TL_coef[,ij] = OLS_model$coefficients
  total_TL_p[,ij] = summary(OLS_model)$coef[,4]
  total_TL_R_squared_prefecture[ij] = summary(OLS_model)$adj.r.squared
  
  hampel_model = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.hampel, maxit = 200)
  total_TL_hampel_coef[,ij] = hampel_model$coefficients
  total_TL_hampel_p[,ij] = find_p_value(t_values = summary(hampel_model)$coef[,3], df = summary(hampel_model)$df[2])
  total_TL_hampel_R_squared_prefecture[ij] = 1 - (sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
  
  bisquare_model = rlm(y ~ x + I(x^2) + I(x^3), psi = psi.bisquare, maxit = 200)
  total_TL_bisquare_coef[,ij] = bisquare_model$coefficients
  total_TL_bisquare_p[,ij] = find_p_value(t_values = summary(bisquare_model)$coef[,3], df = summary(bisquare_model)$df[2])
  total_TL_bisquare_R_squared_prefecture[ij] = 1 - (sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2))*(length(y)-1)/(length(y)-4)
}

