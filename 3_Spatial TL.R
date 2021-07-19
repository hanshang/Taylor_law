##############
# Spatial TL
##############

#################################
# plot spatial mean and variance
#################################

# female

spatial_female_mean = apply(female_data[,,2:48], c(1, 2), mean)
spatial_female_var  = apply(female_data[,,2:48], c(1, 2), var)

rownames(spatial_female_mean) = rownames(spatial_female_var) = 0:100
colnames(spatial_female_mean) = colnames(spatial_female_var) = year

plot(fts(0:100, log(spatial_female_mean, base = 10)), xlab = "Age", ylab = "Female mortality")
plot(fts(0:100, log(spatial_female_var,  base = 10)), xlab = "Age", ylab = "Female mortality")

# male

spatial_male_mean = apply(male_data[,,2:48], c(1, 2), mean)
spatial_male_var  = apply(male_data[,,2:48], c(1, 2), var)

rownames(spatial_male_mean) = rownames(spatial_male_var) = 0:100
colnames(spatial_male_mean) = colnames(spatial_male_var) = year

plot(fts(0:100, log(spatial_male_mean, base = 10)), xlab = "Age", ylab = "Male mortality")
plot(fts(0:100, log(spatial_male_var,  base = 10)), xlab = "Age", ylab = "Male mortality")

# total

spatial_total_mean = apply(total_data[,,2:48], c(1, 2), mean)
spatial_total_var  = apply(total_data[,,2:48], c(1, 2), var)

rownames(spatial_total_mean) = rownames(spatial_total_var) = 0:100
colnames(spatial_total_mean) = colnames(spatial_total_var) = year

plot(fts(0:100, log(spatial_total_mean, base = 10)), xlab = "Age", ylab = "Total mortality")
plot(fts(0:100, log(spatial_total_var,  base = 10)), xlab = "Age", ylab = "Total mortality")


###########################
# Estimate TL coefficients
###########################

# female ver 1

female_TL_slope_dynamic = female_TL_hampel_slope_dynamic = female_TL_bisquare_slope_dynamic = vector("numeric", 43)
female_TL_R_squared_dynamic = female_TL_hampel_R_squared_dynamic = female_TL_bisquare_R_squared_dynamic = vector("numeric", 43)

for(ij in 1:n_year)
{
  OLS_model = lm(log(spatial_female_var[,ij], base = 10) ~ log(spatial_female_mean[,ij], base = 10))
  female_TL_slope_dynamic[ij] = OLS_model$coefficients[2]
  female_TL_R_squared_dynamic[ij] = summary(OLS_model)$r.squared
  
  hampel_model = rlm(log(spatial_female_var[,ij], base = 10) ~ log(spatial_female_mean[,ij], base = 10), psi = psi.hampel, maxit = 200)
  female_TL_hampel_slope_dynamic[ij] = hampel_model$coefficients[2]
  female_TL_hampel_R_squared_dynamic[ij] = 1 - sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2)
  
  bisquare_model = rlm(log(spatial_female_var[,ij], base = 10) ~ log(spatial_female_mean[,ij], base = 10), psi = psi.bisquare, maxit = 200)
  female_TL_bisquare_slope_dynamic[ij] = bisquare_model$coefficients[2]
  female_TL_bisquare_R_squared_dynamic[ij] = 1 - sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2)
}

female_spatial_iid = lm(female_TL_slope_dynamic~year)
female_spatial_hampel_iid = lm(female_TL_hampel_slope_dynamic~year)
female_spatial_bisquare_iid = lm(female_TL_bisquare_slope_dynamic~year)

# female ver 2

female_TL_slope_dynamic_v2 = female_TL_hampel_slope_dynamic_v2 = female_TL_bisquare_slope_dynamic_v2 = vector("numeric", 43)
female_TL_R_squared_dynamic_v2 = female_TL_hampel_R_squared_dynamic_v2 = female_TL_bisquare_R_squared_dynamic_v2 = vector("numeric", 43)

for(ij in 1:n_year)
{
  OLS_model = lm(log(spatial_female_var[,ij], base = 10) ~ log(Japan$rate$female[,ij], base = 10))
  female_TL_slope_dynamic_v2[ij] = OLS_model$coefficients[2]
  female_TL_R_squared_dynamic_v2[ij] = summary(OLS_model)$r.squared
  
  hampel_model = rlm(log(spatial_female_var[,ij], base = 10) ~ log(Japan$rate$female[,ij], base = 10), psi = psi.hampel, maxit = 200)
  female_TL_hampel_slope_dynamic_v2[ij] = hampel_model$coefficients[2]
  female_TL_hampel_R_squared_dynamic_v2[ij] = 1 - sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2)
  
  bisquare_model = rlm(log(spatial_female_var[,ij], base = 10) ~ log(Japan$rate$female[,ij], base = 10), psi = psi.bisquare, maxit = 200)
  female_TL_bisquare_slope_dynamic_v2[ij] = bisquare_model$coefficients[2]
  female_TL_bisquare_R_squared_dynamic_v2[ij] = 1 - sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2)
}

# male ver 1

male_TL_slope_dynamic = male_TL_hampel_slope_dynamic = male_TL_bisquare_slope_dynamic = vector("numeric", 43)
male_TL_R_squared_dynamic = male_TL_hampel_R_squared_dynamic = male_TL_bisquare_R_squared_dynamic = vector("numeric", 43)

for(ij in 1:n_year)
{
  OLS_model = lm(log(spatial_male_var[,ij], base = 10) ~ log(spatial_male_mean[,ij], base = 10))
  male_TL_slope_dynamic[ij] = OLS_model$coefficients[2]
  male_TL_R_squared_dynamic[ij] = summary(OLS_model)$r.squared
  
  hampel_model = rlm(log(spatial_male_var[,ij], base = 10) ~ log(spatial_male_mean[,ij], base = 10), psi = psi.hampel, maxit = 200)
  male_TL_hampel_slope_dynamic[ij] = hampel_model$coefficients[2]
  male_TL_hampel_R_squared_dynamic[ij] = 1 - sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2)
  
  bisquare_model = rlm(log(spatial_male_var[,ij], base = 10) ~ log(spatial_male_mean[,ij], base = 10), psi = psi.bisquare, maxit = 200)
  male_TL_bisquare_slope_dynamic[ij] = bisquare_model$coefficients[2]
  male_TL_bisquare_R_squared_dynamic[ij] = 1 - sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2)
}

male_spatial_iid = lm(male_TL_slope_dynamic~year)
male_spatial_hampel_iid = lm(male_TL_hampel_slope_dynamic~year)
male_spatial_bisquare_iid = lm(male_TL_bisquare_slope_dynamic~year)

# male ver 2

male_TL_slope_dynamic_v2 = male_TL_hampel_slope_dynamic_v2 = male_TL_bisquare_slope_dynamic_v2 = vector("numeric", 43)
male_TL_R_squared_dynamic_v2 = male_TL_hampel_R_squared_dynamic_v2 = male_TL_bisquare_R_squared_dynamic_v2 = vector("numeric", 43)

for(ij in 1:n_year)
{
  OLS_model = lm(log(spatial_male_var[,ij], base = 10) ~ log(Japan$rate$male[,ij], base = 10))
  male_TL_slope_dynamic_v2[ij] = OLS_model$coefficients[2]
  male_TL_R_squared_dynamic_v2[ij] = summary(OLS_model)$r.squared
  
  hampel_model = rlm(log(spatial_male_var[,ij], base = 10) ~ log(Japan$rate$male[,ij], base = 10), psi = psi.hampel, maxit = 200)
  male_TL_hampel_slope_dynamic_v2[ij] = hampel_model$coefficients[2]
  male_TL_hampel_R_squared_dynamic_v2[ij] = 1 - sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2)
  
  bisquare_model = rlm(log(spatial_male_var[,ij], base = 10) ~ log(Japan$rate$male[,ij], base = 10), psi = psi.bisquare, maxit = 200)
  male_TL_bisquare_slope_dynamic_v2[ij] = bisquare_model$coefficients[2]
  male_TL_bisquare_R_squared_dynamic_v2[ij] = 1 - sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2)
}

# total ver 1

total_TL_slope_dynamic = total_TL_hampel_slope_dynamic = total_TL_bisquare_slope_dynamic = vector("numeric", 43)
total_TL_R_squared_dynamic = total_TL_hampel_R_squared_dynamic = total_TL_bisquare_R_squared_dynamic = vector("numeric", 43)


for(ij in 1:n_year)
{
  OLS_model = lm(log(spatial_total_var[,ij], base = 10) ~ log(spatial_total_mean[,ij], base = 10))
  total_TL_slope_dynamic[ij] = OLS_model$coefficients[2]
  total_TL_R_squared_dynamic[ij] = summary(OLS_model)$r.squared
  
  hampel_model = rlm(log(spatial_total_var[,ij], base = 10) ~ log(spatial_total_mean[,ij], base = 10), psi = psi.hampel, maxit = 200)
  total_TL_hampel_slope_dynamic[ij] = hampel_model$coefficients[2]
  total_TL_hampel_R_squared_dynamic[ij] = 1 - sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2)
  
  bisquare_model = rlm(log(spatial_total_var[,ij], base = 10) ~ log(spatial_total_mean[,ij], base = 10), psi = psi.bisquare, maxit = 200)
  total_TL_bisquare_slope_dynamic[ij] = bisquare_model$coefficients[2]
  total_TL_bisquare_R_squared_dynamic[ij] = 1 - sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2)
}

total_spatial_iid = lm(total_TL_slope_dynamic~year)
total_spatial_hampel_iid = lm(total_TL_hampel_slope_dynamic~year)
total_spatial_bisquare_iid = lm(total_TL_bisquare_slope_dynamic~year)

# total ver 2

total_TL_slope_dynamic_v2 = total_TL_hampel_slope_dynamic_v2 = total_TL_bisquare_slope_dynamic_v2 = vector("numeric", 43)
total_TL_R_squared_dynamic_v2 = total_TL_hampel_R_squared_dynamic_v2 = total_TL_bisquare_R_squared_dynamic_v2 = vector("numeric", 43)

for(ij in 1:n_year)
{
  OLS_model = lm(log(spatial_total_var[,ij], base = 10) ~ log(Japan$rate$total[,ij], base = 10))
  total_TL_slope_dynamic_v2[ij] = OLS_model$coefficients[2]
  total_TL_R_squared_dynamic_v2[ij] = summary(OLS_model)$r.squared
  
  hampel_model = rlm(log(spatial_total_var[,ij], base = 10) ~ log(Japan$rate$total[,ij], base = 10), psi = psi.hampel, maxit = 200)
  total_TL_hampel_slope_dynamic_v2[ij] = hampel_model$coefficients[2]
  total_TL_hampel_R_squared_dynamic_v2[ij] = 1 - sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2)
  
  bisquare_model = rlm(log(spatial_total_var[,ij], base = 10) ~ log(Japan$rate$total[,ij], base = 10), psi = psi.bisquare, maxit = 200)
  total_TL_bisquare_slope_dynamic_v2[ij] = bisquare_model$coefficients[2]
  total_TL_bisquare_R_squared_dynamic_v2[ij] = 1 - sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2)
}

