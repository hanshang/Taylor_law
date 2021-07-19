#############################
# Temporal long-run variance
#############################

library(sandwich)
library(MASS)

temporal_long_run_female_var = apply(female_data[,,2:48], c(1, 3), lrvar, type = "Andrews")
temporal_long_run_male_var = apply(male_data[,,2:48], c(1, 3), lrvar, type = "Andrews")
temporal_long_run_total_var = apply(total_data[,,2:48], c(1, 3), lrvar, type = "Andrews")


# female
female_long_run_TL_slope_prefecture = female_long_run_TL_hampel_slope_prefecture = female_long_run_TL_bisquare_slope_prefecture = vector("numeric", 47)
female_long_run_TL_R_squared_prefecture = female_long_run_TL_hampel_R_squared_prefecture = female_long_run_TL_bisquare_R_squared_prefecture = vector("numeric", 47)

for(ij in 1:47)
{
  OLS_model = lm(log(temporal_long_run_female_var[,ij], base = 10) ~ log(temporal_female_mean[,ij], base = 10))
  female_long_run_TL_slope_prefecture[ij] = OLS_model$coefficients[2]
  female_long_run_TL_R_squared_prefecture[ij] = summary(OLS_model)$r.squared
  
  hampel_model = rlm(log(temporal_long_run_female_var[,ij], base = 10) ~ log(temporal_female_mean[,ij], base = 10), psi = psi.hampel, maxit = 500)
  female_long_run_TL_hampel_slope_prefecture[ij] = hampel_model$coefficients[2]
  female_long_run_TL_hampel_R_squared_prefecture[ij] = 1 - sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2)
  
  bisquare_model = rlm(log(temporal_long_run_female_var[,ij], base = 10) ~ log(temporal_female_mean[,ij], base = 10), psi = psi.bisquare, maxit = 500)
  female_long_run_TL_bisquare_slope_prefecture[ij] = bisquare_model$coefficients[2]
  female_long_run_TL_bisquare_R_squared_prefecture[ij] = 1 - sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2)
}

# male

male_long_run_TL_slope_prefecture = male_long_run_TL_hampel_slope_prefecture = male_long_run_TL_bisquare_slope_prefecture = vector("numeric", 47)
male_long_run_TL_R_squared_prefecture = male_long_run_TL_hampel_R_squared_prefecture = male_long_run_TL_bisquare_R_squared_prefecture = vector("numeric", 47)

for(ij in 1:47)
{
  OLS_model = lm(log(temporal_long_run_male_var[,ij], base = 10) ~ log(temporal_male_mean[,ij], base = 10))
  male_long_run_TL_slope_prefecture[ij] = OLS_model$coefficients[2]
  male_long_run_TL_R_squared_prefecture[ij] = summary(OLS_model)$r.squared
  
  hampel_model = rlm(log(temporal_long_run_male_var[,ij], base = 10) ~ log(temporal_male_mean[,ij], base = 10), psi = psi.hampel, maxit = 200)
  male_long_run_TL_hampel_slope_prefecture[ij] = hampel_model$coefficients[2]
  male_long_run_TL_hampel_R_squared_prefecture[ij] = 1 - sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2)
  
  bisquare_model = rlm(log(temporal_long_run_male_var[,ij], base = 10) ~ log(temporal_male_mean[,ij], base = 10), psi = psi.bisquare, maxit = 200)
  male_long_run_TL_bisquare_slope_prefecture[ij] = bisquare_model$coefficients[2]
  male_long_run_TL_bisquare_R_squared_prefecture[ij] = 1 - sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2)
}

# total

total_long_run_TL_slope_prefecture = total_long_run_TL_hampel_slope_prefecture = total_long_run_TL_bisquare_slope_prefecture = vector("numeric", 47)
total_long_run_TL_R_squared_prefecture = total_long_run_TL_hampel_R_squared_prefecture = total_long_run_TL_bisquare_R_squared_prefecture = vector("numeric", 47)

for(ij in 1:47)
{
  OLS_model = lm(log(temporal_long_run_total_var[,ij], base = 10) ~ log(temporal_total_mean[,ij], base = 10))
  total_long_run_TL_slope_prefecture[ij] = OLS_model$coefficients[2]
  total_long_run_TL_R_squared_prefecture[ij] = summary(OLS_model)$r.squared
  
  hampel_model = rlm(log(temporal_long_run_total_var[,ij], base = 10) ~ log(temporal_total_mean[,ij], base = 10), psi = psi.hampel, maxit = 200)
  total_long_run_TL_hampel_slope_prefecture[ij] = hampel_model$coefficients[2]
  total_long_run_TL_hampel_R_squared_prefecture = 1 - sum((hampel_model$residuals)^2)/sum((hampel_model$model[,1] - mean(hampel_model$model[,1]))^2)
  
  bisquare_model = rlm(log(temporal_long_run_total_var[,ij], base = 10) ~ log(temporal_total_mean[,ij], base = 10), psi = psi.bisquare, maxit = 200)
  total_long_run_TL_bisquare_slope_prefecture[ij] = bisquare_model$coefficients[2]
  total_long_run_TL_bisquare_R_squared_prefecture[ij] = 1 - sum((bisquare_model$residuals)^2)/sum((bisquare_model$model[,1] - mean(bisquare_model$model[,1]))^2)
}


############################################
# K-means clustering of TL slope parameters
############################################

library(factoextra)

##############
# Temporal TL
##############

# female
state_prefecture_female_df = data.frame(a_linear = female_TL_coef[1,], b_linear = female_TL_coef[2,], c_linear = female_TL_coef[3,], d_linear = female_TL_coef[4,], row.names = state[-1])

fviz_nbclust_plot(state_prefecture_female_df, kmeans, method = "wss")
fviz_nbclust_plot(state_prefecture_female_df, kmeans, method = "silhouette")

kmeans_female = kmeans(state_prefecture_female_df, centers = 2, nstart = 25)
kmeans_female_plot = fviz_cluster(kmeans_female, data = state_prefecture_female_df, main = "")

# male
state_prefecture_male_df = data.frame(a_linear = male_TL_coef[1,], b_linear = male_TL_coef[2,], c_linear = male_TL_coef[3,], d_linear = male_TL_coef[4,], row.names = state[-1])

fviz_nbclust_plot(state_prefecture_male_df, kmeans, method = "wss")
fviz_nbclust(state_prefecture_male_df, kmeans, method = "silhouette")

kmeans_male = kmeans(state_prefecture_male_df, centers = 2, nstart = 25)
kmeans_male_plot = fviz_cluster(kmeans_male, data = state_prefecture_male_df, main = "")

# total
state_prefecture_total_df = data.frame(a_linear = total_TL_coef[1,], b_linear = total_TL_coef[2,], c_linear = total_TL_coef[3,], d_linear = total_TL_coef[4,], row.names = state[-1])

fviz_nbclust_plot(state_prefecture_total_df, kmeans, method = "wss")
fviz_nbclust(state_prefecture_total_df, kmeans, method = "silhouette")

kmeans_total = kmeans(state_prefecture_total_df, centers = 2, nstart = 25)
kmeans_total_plot = fviz_cluster(kmeans_total, data = state_prefecture_total_df, main = "")

## summary of clustering memberships
k_members_total = data.frame(Region = attributes(kmeans_total$cluster)$names, Cluster_total = as.vector(kmeans_total$cluster))
k_members_all = list(k_members_total, k_members_female, k_members_male) %>% reduce(left_join, by = "Region")

#######################
# Long-run temporal TL
#######################

# female
state_long_run_prefecture_female_df = data.frame( a_bisquare = female_TL_bisquare_coef[1,], b_bisquare = female_TL_bisquare_coef[2,], c_bisquare = female_TL_bisquare_coef[3,], d_bisquare = female_TL_bisquare_coef[4,], row.names = state[-1])

fviz_nbclust_plot(state_long_run_prefecture_female_df, kmeans, method = "wss")
fviz_nbclust(state_long_run_prefecture_female_df, kmeans, method = "silhouette")

kmeans_female_long_run = kmeans(state_long_run_prefecture_female_df, centers = 2, nstart = 25)
kmeans_female_plot_long_run = fviz_cluster(kmeans_female_long_run, data = state_long_run_prefecture_female_df, main = "", ylab = FALSE) + theme(legend.position="none") + theme(axis.title.x = element_text(size=20))

k_members_female_long_run = data.frame(Region = attributes(kmeans_female_long_run$cluster)$names, Cluster_female = as.vector(kmeans_female_long_run$cluster))

# male
state_long_run_prefecture_male_df = data.frame(a_bisquare = male_TL_bisquare_coef[1,], b_bisquare = male_TL_bisquare_coef[2,], c_bisquare = male_TL_bisquare_coef[3,], d_bisquare = male_TL_bisquare_coef[4,], row.names = state[-1])

fviz_nbclust_plot(state_long_run_prefecture_male_df, kmeans, method = "wss")
fviz_nbclust(state_long_run_prefecture_male_df, kmeans, method = "silhouette")

kmeans_male_long_run = kmeans(state_long_run_prefecture_male_df, centers = 2, nstart = 25)
kmeans_male_plot_long_run = fviz_cluster(kmeans_male_long_run, data = state_long_run_prefecture_male_df, main = "", ylab = FALSE)  + theme(legend.title = element_text(size=20), legend.text = element_text(size=15)) + theme(axis.title.x = element_text(size = 20))

k_members_male_long_run = data.frame(Region = attributes(kmeans_male_long_run$cluster)$names, Cluster_male = as.vector(kmeans_male_long_run$cluster))

# total

state_long_run_prefecture_total_df = data.frame(a_bisquare = total_TL_bisquare_coef[1,], b_bisquare = total_TL_bisquare_coef[2,], c_bisquare = total_TL_bisquare_coef[3,], d_bisquare = total_TL_bisquare_coef[4,], row.names = state[-1])

fviz_nbclust_plot(state_long_run_prefecture_total_df, kmeans, method = "wss")
fviz_nbclust(state_long_run_prefecture_total_df, kmeans, method = "silhouette")

kmeans_total_long_run = kmeans(state_long_run_prefecture_total_df, centers = 2, nstart = 25)
kmeans_total_plot_long_run = fviz_cluster(kmeans_total_long_run, data = state_long_run_prefecture_total_df, main = "") + theme(axis.title = element_text(size = 20)) + theme(legend.position="none")

k_members_total_long_run = data.frame(Region = attributes(kmeans_total_long_run$cluster)$names, Cluster_total = as.vector(kmeans_total_long_run$cluster))

k_members_all_long_run = list(k_members_total_long_run, k_members_female_long_run, k_members_male_long_run) %>% reduce(left_join, by = "Region")

################################
# Rangking of Prefecture by GDP
################################

GDP_prefecture <- read.csv("X:/Dropbox/Working projects/TL_JMD/TL project/Japan_prefecture.csv")
library(tidyverse)

GDP_prefecture[GDP_prefecture$Region == "Gumma", "Region"] = "Gunma"

# Set variables
all_measure = unique(GDP_prefecture$Measure)

# 2016 ranking
GDP_2016 = GDP_prefecture[,c("Year", "Region", "Value", "Measure")] %>% filter(Year == 2016, Region != "Japan", Measure == all_measure[1]) %>% arrange(Year, -Value, Region)
GDP_2016 = mutate(GDP_2016) %>% select(Region, Value)

GDP_2016_final = merge(GDP_2016, k_members_all_long_run, by = "Region") %>% arrange(-Value, Region, Cluster_total, Cluster_female, Cluster_male)

GDP_2016_ranking = cbind(GDP_2016, 1:47)
colnames(GDP_2016_ranking)[2:3] = c("GDP" , "GDP Ranking")

# GDP table output
GDP_table_13 = merge(prefecture_index, GDP_2016_ranking, by = "Region")
GDP_table_18 = merge(GDP_table_13, k_members_all_long_run, by = "Region")
GDP_table_1_12 = merge(GDP_table_18, k_members_all, by = "Region") %>% arrange(-GDP, Region, Index, Cluster_total, Cluster_female, Cluster_male, total_OLS, female_OLS, male_OLS)

