###############################
# functions used to save plots
###############################

savefig = function (filename, height=10, width = (1 + sqrt(5))/2*height, type=c("eps","pdf","jpg","png"), pointsize = 10, family = "Helvetica", sublines = 0, toplines = 0, leftlines = 0, res=300) 
{
    type <- match.arg(type)
    filename <- paste(filename, ".", type, sep = "")
    if(type=="eps")
    {
        postscript(file = filename, horizontal = FALSE, 
                width = width/2.54, height = height/2.54, pointsize = pointsize, 
                family = family, onefile = FALSE, print.it = FALSE)
    }
    else if(type=="pdf")
    {
        pdf(file = filename, width=width/2.54, height=height/2.54, pointsize=pointsize,
            family=family, onefile=TRUE)
    }
    else if(type=="jpg")
    {
        jpeg(filename=filename, width=width, height=height, res=res,quality=100, units="cm")#, pointsize=pointsize*50)
    }
    else if(type=="png")
    {
        png(filename=filename, width=width, height=height, res=res, units="cm")#, pointsize=pointsize*50)
    }
    else
        stop("Unknown file type")
    par(mgp = c(2.2, 0.45, 0), tcl = -0.4, mar = c(3.2 + sublines + 0.25 * (sublines > 0), 
         3.5 + leftlines, 1 + toplines, 1) + 0.1)
    par(pch = 1)
    invisible()
}

savepdf = function(...)
{
    savefig(...,type="pdf")
}

################
# Rainbow plots
################

Japan_log10 = Japan
Japan_log10$rate$female = log10(Japan$rate$female)
Japan_log10$rate$male = log10(Japan$rate$male)
Japan_log10$rate$total = log10(Japan$rate$total)

savepdf("Japan_female", width = 12, height = 10, toplines = 1)
plot(Japan_log10, series = "female", transform = FALSE, ylab = "Death rate (logarithm base 10)")
dev.off()

savepdf("Japan_male", width = 12, height = 10, toplines = 1)
plot(Japan_log10, series = "male", transform = FALSE, ylab = "Death rate (logarithm base 10)")
dev.off()

savepdf("Japan_total", width = 12, height = 10, toplines = 1)
plot(Japan_log10, series = "total", transform = FALSE, ylab = "Death rate (logarithm base 10)")
dev.off()

##############
# Image plots
##############

require(RColorBrewer)

# ratio between each prefecture and total 

prefecture_total = prefecture_female = prefecture_male = array(, dim = c(101, 44, 47))
for(iw in 2:48)
{
    gettotal <- get(state[iw])$rate$total
    gettotal[gettotal==0] <- NA
    getmale <- get(state[iw])$rate$male
    getmale[getmale==0] <- NA
    getfemale <- get(state[iw])$rate$female
    getfemale[getfemale==0] <- NA
    prefecture_total[,,iw-1]  = log(gettotal/Japan$rate$total)
    prefecture_female[,,iw-1] = log(getfemale/Japan$rate$female) 
    prefecture_male[,,iw-1]   = log(getmale/Japan$rate$male)
}  

# raw (averaged over age and prefecture)

total_rate_raw = apply(prefecture_total, c(1,3), mean, na.rm = TRUE)
female_rate_raw = apply(prefecture_female, c(1,3), mean, na.rm = TRUE)
male_rate_raw = apply(prefecture_male, c(1,3), mean, na.rm = TRUE)
colnames(male_rate_raw) = colnames(female_rate_raw) = colnames(total_rate_raw) = state[2:48]

total_rate_raw_arrange = total_rate_raw[,47:1]
female_rate_raw_arrange = female_rate_raw[,47:1]
male_rate_raw_arrange = male_rate_raw[,47:1]

# raw (averaged over year and prefecture)

total_rate_raw_year  = apply(prefecture_total, c(2,3), mean, na.rm = TRUE)
female_rate_raw_year = apply(prefecture_female, c(2,3), mean, na.rm = TRUE)
male_rate_raw_year   = apply(prefecture_male, c(2,3), mean, na.rm = TRUE)
colnames(male_rate_raw_year) = colnames(female_rate_raw_year) = colnames(total_rate_raw_year) = state[2:48]

total_rate_raw_arrange_year = total_rate_raw_year[,47:1]
female_rate_raw_arrange_year = female_rate_raw_year[,47:1]
male_rate_raw_arrange_year = male_rate_raw_year[,47:1]

# raw (averaged over year and age)
total_rate_raw_prefecture = apply(prefecture_total, c(1,2), mean, na.rm = TRUE)
female_rate_raw_prefecture = apply(prefecture_female, c(1,2), mean, na.rm = TRUE)
male_rate_raw_prefecture = apply(prefecture_male, c(1,2), mean, na.rm = TRUE)
colnames(total_rate_raw_prefecture) = colnames(female_rate_raw_prefecture) = colnames(male_rate_raw_prefecture) = 1975:2018

require(RColorBrewer)
par(mfrow = c(2, 3))
image(0:100, 1:47, total_rate_raw_arrange, col = brewer.pal(7, "RdBu"), xlab = "Age", ylab = "Prefecture", main = "Total", zlim = c(-1, 1)) 
box()
image(0:100, 1:47, female_rate_raw_arrange, col = brewer.pal(7, "RdBu"), xlab = "Age", ylab = "", main = "Female", zlim = c(-1, 1))
box()
image(0:100, 1:47, male_rate_raw_arrange, col = brewer.pal(7, "RdBu"), xlab = "Age", ylab = "", main = "Male", zlim = c(-1, 1))
box()

image(1975:2018, 1:47, total_rate_raw_arrange_year, col = brewer.pal(7, "RdBu"), xlab = "Year", ylab = "Prefecture", zlim = c(-1, 1))
box()
image(1975:2018, 1:47, female_rate_raw_arrange_year, col = brewer.pal(7, "RdBu"), xlab = "Year", ylab = "", zlim = c(-1, 1))
box()
image(1975:2018, 1:47, male_rate_raw_arrange_year, col = brewer.pal(7, "RdBu"), xlab = "Year", ylab = "", zlim = c(-1, 1))
box()

savepdf("image_plot_1", width = 22, height = 7, toplines = 0.8, pointsize = 12)
par(mfrow = c(1,3))
image(0:100, 1:47, total_rate_raw_arrange, col = brewer.pal(7, "RdBu"), xlab = "Age", ylab = "Prefecture", main = "Total", zlim = c(-1, 1))
image(0:100, 1:47, female_rate_raw_arrange, col = brewer.pal(7, "RdBu"), xlab = "Age", ylab = "", main = "Female", zlim = c(-1, 1))
image(0:100, 1:47, male_rate_raw_arrange, col = brewer.pal(7, "RdBu"), xlab = "Age", ylab = "", main = "Male", zlim = c(-1, 1))
dev.off()

savepdf("image_plot_2", width = 22, height = 7, toplines = 0.8, pointsize = 12)
par(mfrow = c(1,3))
image(1975:2018, 1:47, total_rate_raw_arrange_year, col = brewer.pal(7, "RdBu"), xlab = "Year", ylab = "Prefecture", zlim = c(-1, 1))
image(1975:2018, 1:47, female_rate_raw_arrange_year, col = brewer.pal(7, "RdBu"), xlab = "Year", ylab = "", zlim = c(-1, 1))
image(1975:2018, 1:47, male_rate_raw_arrange_year, col = brewer.pal(7, "RdBu"), xlab = "Year", ylab = "", zlim = c(-1, 1))
dev.off()

savepdf("image_plot_3", width = 22, height = 7, toplines = 0.8, pointsize = 12)
par(mfrow = c(1,3))
image(0:100, 1975:2018, total_rate_raw_prefecture, col = brewer.pal(7, "RdBu"), xlab = "Age", ylab = "Year", main = "Total", zlim = c(-1, 1))
image(0:100, 1975:2018, female_rate_raw_prefecture, col = brewer.pal(7, "RdBu"), xlab = "Age", ylab = "", main = "Female", zlim = c(-1, 1))
image(0:100, 1975:2018, male_rate_raw_prefecture, col = brewer.pal(7, "RdBu"), xlab = "Age", ylab = "", main = "Male", zlim = c(-1, 1))
dev.off()

###########################
# Temporal TL coefficients
###########################

# female
savepdf("female_TL_temporal_order1", width = 12, height = 10, toplines = 0.8)
plot(1:47, female_TL_coef[2,], xlab = "", 
     ylab = "", ylim = c(-0.9, 1.55))
     #,main = expression(paste(log[10], "(temporal variance) versus ", log[10], "(temporal mean)")))
points(1:47, female_TL_hampel_coef[2,], col = 2, pch = 2)
points(1:47, female_TL_bisquare_coef[2,], col = 4, pch = 4)
#legend("topright", c("LM", "RLM (Hampel weight)", "RLM (Bisquare weight)"), col = c(1, 2, 4), pch = c(1, 2, 4), cex = 1)
dev.off()

savepdf("female_TL_temporal_order2", width = 12, height = 10, toplines = 0.8)
plot(1:47, female_TL_coef[3,], xlab = "", 
     ylab = "", ylim = c(-1.6,-0.3))
#,main = expression(paste(log[10], "(temporal variance) versus ", log[10], "(temporal mean)")))
points(1:47, female_TL_hampel_coef[3,], col = 2, pch = 2)
points(1:47, female_TL_bisquare_coef[3,], col = 4, pch = 4)
dev.off()

savepdf("female_TL_temporal_order3", width = 12, height = 10, toplines = 0.8)
plot(1:47, female_TL_coef[4,], xlab = "Prefecture ordered geographically by North to South", 
     ylab = "", ylim = c(-0.26,-0.06))
#,main = expression(paste(log[10], "(temporal variance) versus ", log[10], "(temporal mean)")))
points(1:47, female_TL_hampel_coef[4,], col = 2, pch = 2)
points(1:47, female_TL_bisquare_coef[4,], col = 4, pch = 4)
dev.off()

# male
savepdf("male_TL_temporal_order1", width = 12, height = 10, toplines = 0.8)
plot(1:47, male_TL_coef[2,], xlab = "", 
     ylab = "", ylim = c(0.2,3.1))
#,main = expression(paste(log[10], "(temporal variance) versus ", log[10], "(temporal mean)")))
points(1:47, male_TL_hampel_coef[2,], col = 2, pch = 2)
points(1:47, male_TL_bisquare_coef[2,], col = 4, pch = 4)
#legend("topright", c("LM", "RLM (Hampel weight)", "RLM (Bisquare weight)"), col = c(1, 2, 4), pch = c(1, 2, 4), cex = 1)
dev.off()

savepdf("male_TL_temporal_order2", width = 12, height = 10, toplines = 0.8)
plot(1:47, male_TL_coef[3,], xlab = "", 
     ylab = "", ylim = c(-1,0.5))
#,main = expression(paste(log[10], "(temporal variance) versus ", log[10], "(temporal mean)")))
points(1:47, male_TL_hampel_coef[3,], col = 2, pch = 2)
points(1:47, male_TL_bisquare_coef[3,], col = 4, pch = 4)
dev.off()

savepdf("male_TL_temporal_order3", width = 12, height = 10, toplines = 0.8)
plot(1:47, male_TL_coef[4,], xlab = "Prefecture ordered geographically by North to South", 
     ylab = "", ylim = c(-0.17,0.05))
#,main = expression(paste(log[10], "(temporal variance) versus ", log[10], "(temporal mean)")))
points(1:47, male_TL_hampel_coef[4,], col = 2, pch = 2)
points(1:47, male_TL_bisquare_coef[4,], col = 4, pch = 4)
dev.off()

# total
savepdf("total_TL_temporal_order1", width = 12, height = 10, toplines = 0.8)
plot(1:47, total_TL_coef[2,], xlab = "", 
     ylab = expression(paste(log[10], "(Temporal Mean) estimates")), ylim = c(-0.5,3), main = "")
#,main = expression(paste(log[10], "(temporal variance) versus ", log[10], "(temporal mean)")))
points(1:47, total_TL_hampel_coef[2,], col = 2, pch = 2)
points(1:47, total_TL_bisquare_coef[2,], col = 4, pch = 4)
legend("topleft", c("LM", "RLM (Hampel weight)", "RLM (Bisquare weight)"), col = c(1, 2, 4), pch = c(1, 2, 4), cex = 1)
dev.off()

savepdf("total_TL_temporal_order2", width = 12, height = 10, toplines = 0.8)
plot(1:47, total_TL_coef[3,], xlab = "", 
     ylab = expression(paste(log[10]^2, "(Temporal Mean) estimates")), ylim = c(-1.3,0.1), main = "")
#,main = expression(paste(log[10], "(temporal variance) versus ", log[10], "(temporal mean)")))
points(1:47, total_TL_hampel_coef[3,], col = 2, pch = 2)
points(1:47, total_TL_bisquare_coef[3,], col = 4, pch = 4)
dev.off()

savepdf("total_TL_temporal_order3", width = 12, height = 10, toplines = 0.8)
plot(1:47, total_TL_coef[4,], xlab = "Prefecture ordered geographically by North to South", 
     ylab = expression(paste(log[10]^3, "(Temporal Mean) estimates")), ylim = c(-0.21, -0.03), main = "")
#,main = expression(paste(log[10], "(temporal variance) versus ", log[10], "(temporal mean)")))
points(1:47, total_TL_hampel_coef[4,], col = 2, pch = 2)
points(1:47, total_TL_bisquare_coef[4,], col = 4, pch = 4)
dev.off()


###########################
# Order 3 TL scatter plots
###########################

for(ij in 1:47)
{
  # female
  age_select = setdiff(age_all, as.numeric(rownames(female_TL_outliers[[ij]])))
  y_female = log(temporal_female_var[age_all %in% age_select,ij], base = 10)
  x_female = log(temporal_female_mean[age_all %in% age_select,ij], base = 10)
  
  OLS_model_female = lm(y_female ~ x_female + I(x_female^2) + I(x_female^3))
  newdata_female = data.frame(x_female = seq(from = min(x_female)-0.1, to = max(x_female)+0.2, by = 0.01))
  pred_female = predict(OLS_model_female, newdata_female, interval = "predict")
  
  savepdf(paste("prediction", ij, "female", sep = "_"), width = 20, height = 20, toplines = 0.8)
  par(mar = c(5,7,5,3))
  plot(x = log10(temporal_female_mean[,ij]), y = log10(temporal_female_var[,ij]), xlim = c(min(log10(temporal_female_mean[,ij]))-0.2, max(log10(temporal_female_mean[,ij]))+0.2), xlab = expression(paste(log[10], "(Mean)")), ylab = expression(paste(log[10], "(Variance)")), cex.lab = 3, col = rainbow(101, end = 0.9))
  if(nrow(female_TL_outliers[[ij]]) > 0)
  {
    points(female_TL_outliers[[ij]][,1], female_TL_outliers[[ij]][,2], pch = 21, bg = 3:(nrow(female_TL_outliers[[ij]])+2), col = "black", cex = 2)
    legend("topleft", col = 3:(nrow(female_TL_outliers[[ij]])+2), pch = 16, legend = rownames(female_TL_outliers[[ij]]), cex = 2)
  }
  title(main = state[ij+1], line = 1, cex.main = 3, font.main= 2)
  lines(x = newdata_female$x, y = pred_female[,1], col = 2, lwd = 2)
  polygon(c(rev(newdata_female$x), newdata_female$x), c(rev(pred_female[,3]), pred_female[,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)
  dev.off()
  
  # male
  age_select = setdiff(age_all, as.numeric(rownames(male_TL_outliers[[ij]])))
  y_male = log(temporal_male_var[age_all %in% age_select,ij], base = 10)
  x_male = log(temporal_male_mean[age_all %in% age_select,ij], base = 10)
  
  OLS_model_male = lm(y_male ~ x_male + I(x_male^2) + I(x_male^3))
  newdata_male = data.frame(x_male = seq(from = min(x_male, na.rm = T)-0.1, to = max(x_male, na.rm = T)+0.2, by = 0.01))
  pred_male = predict(OLS_model_male, newdata_male, interval = "predict")
  
  savepdf(paste("prediction", ij, "male", sep = "_"), width = 20, height = 20, toplines = 0.8)
  par(mar = c(5,7,5,3))
  plot(x = log10(temporal_male_mean[,ij]), y = log10(temporal_male_var[,ij]), xlim = c(min(log10(temporal_male_mean[,ij]), na.rm = T)-0.2, max(log10(temporal_male_mean[,ij]), na.rm = T)+0.2), xlab = expression(paste(log[10], "(Mean)")), ylab = expression(paste(log[10], "(Variance)")), cex.lab = 3, col = rainbow(101, end = 0.9))
  if(nrow(male_TL_outliers[[ij]]) > 0)
  {
    points(male_TL_outliers[[ij]][,1], male_TL_outliers[[ij]][,2], pch = 21, bg = 3:(nrow(male_TL_outliers[[ij]])+2), col = "black", cex = 2)
    legend("topleft", col = 3:(nrow(male_TL_outliers[[ij]])+2), pch = 16, legend = rownames(male_TL_outliers[[ij]]), cex = 2)
  }
  title(main = state[ij+1], line = 1, cex.main = 3, font.main= 2)
  lines(x = newdata_male$x, y = pred_male[,1], col = 2, lwd = 2)
  polygon(c(rev(newdata_male$x), newdata_male$x), c(rev(pred_male[,3]), pred_male[,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)
  dev.off()
  
  # total
  age_select = setdiff(age_all, as.numeric(rownames(total_TL_outliers[[ij]])))
  y_total = log(temporal_total_var[age_all %in% age_select,ij], base = 10)
  x_total = log(temporal_total_mean[age_all %in% age_select,ij], base = 10)
  
  OLS_model_total = lm(y_total ~ x_total + I(x_total^2) + I(x_total^3))
  newdata_total = data.frame(x_total = seq(from = min(x_total)-0.1, to = max(x_total)+0.2, by = 0.01))
  pred_total = predict(OLS_model_total, newdata_total, interval = "predict")
  
  savepdf(paste("prediction", ij, "total", sep = "_"), width = 20, height = 20, toplines = 0.8)
  par(mar = c(5,7,5,3))
  plot(x = log10(temporal_total_mean[,ij]), y = log10(temporal_total_var[,ij]), xlim = c(min(log10(temporal_total_mean[,ij]))-0.2, max(log10(temporal_total_mean[,ij]))+0.2), xlab = expression(paste(log[10], "(Mean)")), ylab = expression(paste(log[10], "(Variance)")), cex.lab = 3, col = rainbow(101, end = 0.9))
  if(nrow(total_TL_outliers[[ij]]) > 0)
  {
    points(total_TL_outliers[[ij]][,1], total_TL_outliers[[ij]][,2], pch = 21, bg = 3:(nrow(total_TL_outliers[[ij]])+2), col = "black", cex = 2)
    legend("topleft", col = 3:(nrow(total_TL_outliers[[ij]])+2), pch = 16, legend = rownames(total_TL_outliers[[ij]]), cex = 2)
  }
  title(main = state[ij+1], line = 1, cex.main = 3, font.main= 2)
  lines(x = newdata_total$x, y = pred_total[,1], col = 2, lwd = 2)
  polygon(c(rev(newdata_total$x), newdata_total$x), c(rev(pred_total[,3]), pred_total[,2]), col = rgb(0.7,0.7,0.7,0.4) , border = NA)
  dev.off()
  
}


##########################
# Weighted spatial TL coefficients
##########################

time_select = 1975:2018

# female
savepdf("female_TL_spatial_intercept", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(time_select, spatial_female_TL_coef[1,], xlab = "", ylab = "")
points(time_select, spatial_female_TL_hampel_coef[1,], col = 2, pch = 2)
points(time_select, spatial_female_TL_bisquare_coef[1,], col = 4, pch = 4)
dev.off()

savepdf("female_TL_spatial_order1", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(time_select, spatial_female_TL_coef[2,], xlab = "", ylab = "")
points(time_select, spatial_female_TL_hampel_coef[2,], col = 2, pch = 2)
points(time_select, spatial_female_TL_bisquare_coef[2,], col = 4, pch = 4)
dev.off()

savepdf("female_TL_spatial_order2", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(time_select, spatial_female_TL_coef[3,], xlab = "", ylab = "")
points(time_select, spatial_female_TL_hampel_coef[3,], col = 2, pch = 2)
points(time_select, spatial_female_TL_bisquare_coef[3,], col = 4, pch = 4)
dev.off()

savepdf("female_TL_spatial_order3", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(time_select, spatial_female_TL_coef[4,], xlab = "Year", ylab = "", cex.lab = 1.5)
points(time_select, spatial_female_TL_hampel_coef[4,], col = 2, pch = 2)
points(time_select, spatial_female_TL_bisquare_coef[4,], col = 4, pch = 4)
dev.off()

# male
savepdf("male_TL_spatial_intercept", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(time_select, (spatial_male_TL_coef[1,]), xlab = "", ylab = "", ylim = c(-1.7,0.8))
points(time_select, (spatial_male_TL_hampel_coef[1,]), col = 2, pch = 2)
points(time_select, (spatial_male_TL_bisquare_coef[1,]), col = 4, pch = 4)
dev.off()

savepdf("male_TL_spatial_order1", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(time_select, spatial_male_TL_coef[2,], xlab = "", ylab = "", ylim = c(3, 6.5))
points(time_select, spatial_male_TL_hampel_coef[2,], col = 2, pch = 2)
points(time_select, spatial_male_TL_bisquare_coef[2,], col = 4, pch = 4)
dev.off()

savepdf("male_TL_spatial_order2", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(time_select, spatial_male_TL_coef[3,], xlab = "", ylab = "", ylim = c(0.6, 2.3))
points(time_select, spatial_male_TL_hampel_coef[3,], col = 2, pch = 2)
points(time_select, spatial_male_TL_bisquare_coef[3,], col = 4, pch = 4)
dev.off()

savepdf("male_TL_spatial_order3", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(time_select, spatial_male_TL_coef[4,], xlab = "Year", ylab = "", ylim = c(0.06, 0.32), cex.lab = 1.5)
points(time_select, spatial_male_TL_hampel_coef[4,], col = 2, pch = 2)
points(time_select, spatial_male_TL_bisquare_coef[4,], col = 4, pch = 4)
dev.off()

# total
savepdf("total_TL_spatial_intercept", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(time_select, (spatial_total_TL_coef[1,]), xlab = "", ylab = expression(log[10](hat(a[3]))), main = "", cex.lab = 1.5)
points(time_select, (spatial_total_TL_hampel_coef[1,]), col = 2, pch = 2)
points(time_select, (spatial_total_TL_bisquare_coef[1,]), col = 4, pch = 4)
legend("topright", c("OLS", "Robust (Hampel)", "Robust (Bisquare)"), col = c(1, 2, 4), pch = c(1, 2, 4), cex = 1)
dev.off()

savepdf("total_TL_spatial_order1", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(time_select, spatial_total_TL_coef[2,], xlab = "", ylab = expression(hat(b)[3]), main = "", ylim = c(2.2, 5.5), cex.lab = 1.5)
points(time_select, spatial_total_TL_hampel_coef[2,], col = 2, pch = 2)
points(time_select, spatial_total_TL_bisquare_coef[2,], col = 4, pch = 4)
# legend("topright", c("OLS", "Robust (Hampel)", "Robust (Bisquare)"), col = c(1, 2, 4), pch = c(1, 2, 4), cex = 1)
dev.off()

savepdf("total_TL_spatial_order2", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(time_select, spatial_total_TL_coef[3,], xlab = "", ylab = expression(hat(c)[3]), main = "", cex.lab = 1.5)
points(time_select, spatial_total_TL_hampel_coef[3,], col = 2, pch = 2)
points(time_select, spatial_total_TL_bisquare_coef[3,], col = 4, pch = 4)
dev.off()

savepdf("total_TL_spatial_order3", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(time_select, spatial_total_TL_coef[4,], xlab = "Year", ylab = "",  main = "", cex.lab = 1.5)
title(ylab = expression(hat(d)[3]), cex.lab = 1.5)
points(time_select, spatial_total_TL_hampel_coef[4,], col = 2, pch = 2)
points(time_select, spatial_total_TL_bisquare_coef[4,], col = 4, pch = 4)
dev.off()

############################
# Spatial mean and variance
############################

savepdf("female_spatial_variance", width = 14, height = 10, toplines = 0.8)
par(mar = c(3,3,2,2))
matplot(log(spatial_female_var, base = 10), lty = 1, type = "l", col = "grey", xlab = "", ylab = "", ylim = c(-8.5,0))
lines(log(spatial_female_var[,1], base = 10), col = 2, lwd = 2)
lines(log(spatial_female_var[,37], base = 10), col = 3, lwd = 2)
lines(log(spatial_female_var[,44], base = 10), col = 4, lwd = 2)
dev.off()

savepdf("female_spatial_mean", width = 14, height = 10, toplines = 0.8)
par(mar = c(4,3,2,2))
matplot(log(spatial_female_mean, base = 10), lty = 1, type = "l", col = "grey", xlab = "Age", ylab = "", ylim = c(-4.5,0), cex.lab = 2)
lines(log(spatial_female_mean[,1], base = 10), col = 2, lwd = 2)
lines(log(spatial_female_mean[,37], base = 10), col = 3, lwd = 2)
lines(log(spatial_female_mean[,44], base = 10), col = 4, lwd = 2)
dev.off()

savepdf("male_spatial_variance", width = 14, height = 10, toplines = 0.8)
par(mar = c(3,3,2,2))
matplot(log(spatial_male_var, base = 10), lty = 1, type = "l", col = "grey", xlab = "", ylab = "", ylim = c(-8.5,0))
lines(log(spatial_male_var[,1], base = 10), col = 2, lwd = 2)
lines(log(spatial_male_var[,37], base = 10), col = 3, lwd = 2)
lines(log(spatial_male_var[,44], base = 10), col = 4, lwd = 2)
dev.off()

savepdf("male_spatial_mean", width = 14, height = 10, toplines = 0.8)
par(mar = c(4,3,2,2))
matplot(log(spatial_male_mean, base = 10), lty = 1, type = "l", col = "grey", xlab = "Age", ylab = "", ylim = c(-4.5,0), cex.lab = 2)
lines(log(spatial_male_mean[,1], base = 10), col = 2, lwd = 2)
lines(log(spatial_male_mean[,37], base = 10), col = 3, lwd = 2)
lines(log(spatial_male_mean[,44], base = 10), col = 4, lwd = 2)
dev.off()

savepdf("total_spatial_variance", width = 14, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
matplot(log(spatial_total_var, base = 10), lty = 1, type = "l", col = "grey", xlab = "", ylab = expression(paste(log[10], "(Spatial Variance)")), ylim = c(-8.5,0), cex.lab = 2)
lines(log(spatial_total_var[,1], base = 10), col = 2, lwd = 2)
lines(log(spatial_total_var[,37], base = 10), col = 3, lwd = 2)
lines(log(spatial_total_var[,44], base = 10), col = 4, lwd = 2)
legend("topleft", c("1985","2011","2018"), lty = c(1,1,1), col = c(2,3,4), lwd = c(2,2,2))
dev.off()

savepdf("total_spatial_mean", width = 14, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
matplot(log(spatial_total_mean, base = 10), lty = 1, type = "l", col = "grey", xlab = "Age", ylab = expression(paste(log[10], "(Spatial Mean)")), ylim = c(-4.5,0), cex.lab = 2)
lines(log(spatial_total_mean[,1], base = 10), col = 2, lwd = 2)
lines(log(spatial_total_mean[,37], base = 10), col = 3, lwd = 2)
lines(log(spatial_total_mean[,44], base = 10), col = 4, lwd = 2)
dev.off()

par(mfrow = c(1,3))
plot(log(spatial_female_mean[,1], base = 10), log(spatial_female_var[,1], base = 10))
points(log(spatial_male_mean[,1], base = 10), log(spatial_male_var[,1], base = 10), col = 2)
plot(log(spatial_female_mean[,37], base = 10), log(spatial_female_var[,37], base = 10))
points(log(spatial_male_mean[,37], base = 10), log(spatial_male_var[,37], base = 10), col = 2)
plot(log(spatial_female_mean[,44], base = 10), log(spatial_female_var[,44], base = 10))
points(log(spatial_male_mean[,44], base = 10), log(spatial_male_var[,44], base = 10), col = 2)

####################################
# Long-run temporal TL coefficients
####################################

# female
savepdf("female_long_run_TL_temporal", width = 14, height = 10, toplines = 0.8)
plot(1:47, female_long_run_TL_slope_prefecture, xlab = "Prefecture ordered geographically by North to South", 
     ylab = "Taylor's Law slope parameter estimates", ylim = c(1.6, 2.2))
     #,main = expression(paste(log[10], "(temporal long-run variance) versus ", log[10], "(temporal mean)")))
points(1:47, female_long_run_TL_hampel_slope_prefecture, col = 2, pch = 2)
points(1:47, female_long_run_TL_bisquare_slope_prefecture, col = 4, pch = 4)
legend("topright", c("TL", "RTL (Hampel weight)", "RTL (Bisquare weight)"), col = c(1, 2, 4), pch = c(1, 2, 4), cex = 1)
dev.off()

# male
savepdf("male_long_run_TL_temporal", width = 14, height = 10, toplines = 0.8)
plot(1:47, male_long_run_TL_slope_prefecture, xlab = "Prefecture ordered geographically by North to South", ylab = "Taylor's Law slope parameter estimates", ylim = c(1.6, 2.2))
    # , main = expression(paste(log[10], "(temporal long-run variance) versus ", log[10], "(temporal mean)")))
points(1:47, male_long_run_TL_hampel_slope_prefecture, col = 2, pch = 2)
points(1:47, male_long_run_TL_bisquare_slope_prefecture, col = 4, pch = 4)
#legend("topright", c("LM", "RLM (Hampel weight)", "RLM (Bisquare weight)"), col = c(1, 2, 4), pch = c(1, 2, 4), cex = 0.8)
dev.off()

# total
savepdf("total_long_run_TL_temporal", width = 14, height = 10, toplines = 0.8)
plot(1:47, total_long_run_TL_slope_prefecture, xlab = "Prefecture ordered geographically by North to South", ylab = "Taylor's Law slope parameter estimates", ylim = c(1.6, 2.2))
     #, main = expression(paste(log[10], "(temporal long-run variance) versus ", log[10], "(temporal mean)")))
points(1:47, total_long_run_TL_hampel_slope_prefecture, col = 2, pch = 2)
points(1:47, total_long_run_TL_bisquare_slope_prefecture, col = 4, pch = 4)
#legend("topright", c("LM", "RLM (Hampel weight)", "RLM (Bisquare weight)"), col = c(1, 2, 4), pch = c(1, 2, 4), cex = 0.8)
dev.off()


#############################################
# Plots of temporal coefficients against GDP
#############################################

library(tidyverse)

GDP_NtoS = GDP_2016 %>% arrange(factor(Region, levels = state[-1]))
GDPPP_NtoS = GDP_2016_per %>% arrange(factor(Region, levels = state[-1]))

## total
savepdf("total_temporal_GDP_intercept", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(log10(GDP_NtoS$Value), total_TL_coef[1,], xlab = "", ylab = expression(log[10](hat(a[3]))), main = "", ylim = c(-4.1, -2.3), cex.lab = 1.5)
points(log10(GDP_NtoS$Value), (total_TL_hampel_coef[1,]), col = 2, pch = 2)
points(log10(GDP_NtoS$Value), (total_TL_bisquare_coef[1,]), col = 4, pch = 4)
legend("topright", c("OLS", "Robust (Hampel)", "Robust (Bisquare)"), col = c(1, 2, 4), pch = c(1, 2, 4), cex = 1)
dev.off()

savepdf("total_temporal_GDP_order1", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(log10(GDP_NtoS$Value), total_TL_coef[2,], xlab = "", ylab = expression(hat(b)[3]),  main = "", cex.lab = 1.5, ylim = c(-2.3, 2.1))
points(log10(GDP_NtoS$Value), total_TL_hampel_coef[2,], col = 2, pch = 2)
points(log10(GDP_NtoS$Value), total_TL_bisquare_coef[2,], col = 4, pch = 4)
dev.off()

savepdf("total_temporal_GDP_order2", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(log10(GDP_NtoS$Value), total_TL_coef[3,], xlab = "", ylab = expression(hat(c)[3]), main = "", cex.lab = 1.5, ylim = c(-2.5, -0.1))
points(log10(GDP_NtoS$Value), total_TL_hampel_coef[3,], col = 2, pch = 2)
points(log10(GDP_NtoS$Value), total_TL_bisquare_coef[3,], col = 4, pch = 4)
dev.off()

savepdf("total_temporal_GDP_order3", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(log10(GDP_NtoS$Value), total_TL_coef[4,], xlab = expression(paste(log[10], "(Total GDP)")), 
     ylab = expression(hat(d)[3]), ylim = c(-0.4, -0.07), xlim = c(6.2, 8.2), main = "", cex.lab = 1.5)
points(log10(GDP_NtoS$Value), total_TL_hampel_coef[4,], col = 2, pch = 2)
points(log10(GDP_NtoS$Value), total_TL_bisquare_coef[4,], col = 4, pch = 4)
text_ind = c(which.max(total_TL_coef[4,]), which.min(total_TL_coef[4,]), which.max(log10(GDP_NtoS$Value)))
text(log10(GDP_NtoS$Value)[text_ind], total_TL_hampel_coef[4,text_ind], labels = state[-1][text_ind], cex = 0.8, pos = 4)
dev.off()

## female
savepdf("female_temporal_GDP_intercept", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(log10(GDP_NtoS$Value), female_TL_coef[1,], xlab = "", ylab = "", main = "", ylim = c(-4.3, -2.9), cex.lab = 1.5)
points(log10(GDP_NtoS$Value), (female_TL_hampel_coef[1,]), col = 2, pch = 2)
points(log10(GDP_NtoS$Value), (female_TL_bisquare_coef[1,]), col = 4, pch = 4)
dev.off()

savepdf("female_temporal_GDP_order1", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(log10(GDP_NtoS$Value), female_TL_coef[2,], xlab = "", ylab = "",  main = "", cex.lab = 1.5, ylim = c(-2.7, 0.1))
points(log10(GDP_NtoS$Value), female_TL_hampel_coef[2,], col = 2, pch = 2)
points(log10(GDP_NtoS$Value), female_TL_bisquare_coef[2,], col = 4, pch = 4)
dev.off()

savepdf("female_temporal_GDP_order2", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(log10(GDP_NtoS$Value), female_TL_coef[3,], xlab = "", ylab = "", main = "", ylim = c(-2.6, -1), cex.lab = 1.5)
points(log10(GDP_NtoS$Value), female_TL_hampel_coef[3,], col = 2, pch = 2)
points(log10(GDP_NtoS$Value), female_TL_bisquare_coef[3,], col = 4, pch = 4)
dev.off()

savepdf("female_temporal_GDP_order3", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(log10(GDP_NtoS$Value), female_TL_coef[4,], xlab = expression(paste(log[10], "(Total GDP)")),  ylab = "",  ylim = c(-0.41,-0.17), main = "", cex.lab = 1.5, xlim = c(6.2,8.3))
points(log10(GDP_NtoS$Value), female_TL_hampel_coef[4,], col = 2, pch = 2)
points(log10(GDP_NtoS$Value), female_TL_bisquare_coef[4,], col = 4, pch = 4)
text_ind = c(which.max(female_TL_coef[4,]), which.min(female_TL_coef[4,]), which.max(log10(GDP_NtoS$Value)))
text(log10(GDP_NtoS$Value)[text_ind], female_TL_coef[4,text_ind], labels = state[-1][text_ind], cex = 0.8, pos = 4)
dev.off()

## male
savepdf("male_temporal_GDP_intercept", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(log10(GDP_NtoS$Value), male_TL_coef[1,], xlab = "", ylab = "", main = "",  cex.lab = 1.5, ylim = c(-4.5, -2.1))
points(log10(GDP_NtoS$Value), (male_TL_hampel_coef[1,]), col = 2, pch = 2)
points(log10(GDP_NtoS$Value), (male_TL_bisquare_coef[1,]), col = 4, pch = 4)
dev.off()

savepdf("male_temporal_GDP_order1", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(log10(GDP_NtoS$Value), male_TL_coef[2,], xlab = "", ylab = "",  main = "", cex.lab = 1.5, ylim = c(-2.5,2.9))
points(log10(GDP_NtoS$Value), male_TL_hampel_coef[2,], col = 2, pch = 2)
points(log10(GDP_NtoS$Value), male_TL_bisquare_coef[2,], col = 4, pch = 4)
dev.off()

savepdf("male_temporal_GDP_order2", width = 12, height = 10, toplines = 0.8)
par(mar = c(3,4.5,2,2))
plot(log10(GDP_NtoS$Value), male_TL_coef[3,], xlab = "",  ylab = "", main = "", cex.lab = 1.5, ylim = c(-2.3,0.4))
points(log10(GDP_NtoS$Value), male_TL_hampel_coef[3,], col = 2, pch = 2)
points(log10(GDP_NtoS$Value), male_TL_bisquare_coef[3,], col = 4, pch = 4)
dev.off()

savepdf("male_temporal_GDP_order3", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(log10(GDP_NtoS$Value), male_TL_coef[4,], xlab = expression(paste(log[10], "(Total GDP)")), ylab = "", main = "", cex.lab = 1.5, ylim = c(-0.34,0.02), xlim = c(6.2, 8.3))
points(log10(GDP_NtoS$Value), male_TL_hampel_coef[4,], col = 2, pch = 2)
points(log10(GDP_NtoS$Value), male_TL_bisquare_coef[4,], col = 4, pch = 4)
text_ind = c(which.max(male_TL_coef[4,]), which.min(male_TL_coef[4,]), which.max(log10(GDP_NtoS$Value)), which.min(male_TL_bisquare_coef[4,]))
text(log10(GDP_NtoS$Value)[text_ind], male_TL_coef[4,text_ind], labels = state[-1][text_ind], cex = 0.8, pos = 4)
dev.off()


#################################
# Plots of Dimension1 against d3
#################################

## total
savepdf("total_dim1_temporal_intercept", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(x = kmeans_total_plot_long_run$plot_env$data$x, y = total_TL_coef[1,], xlab = "", ylab = expression(hat(a)[3]), cex.lab = 1.5, ylim = c(-4.2,-2.3))
points(x = kmeans_total_plot_long_run$plot_env$data$x, y = total_TL_hampel_coef[1,], col = 2, pch = 2)
points(x = kmeans_total_plot_long_run$plot_env$data$x, y = total_TL_bisquare_coef[1,], col = 4, pch = 4)
legend("topleft", c("OLS", "Robust (Hampel)", "Robust (Bisquare)"), col = c(1, 2, 4), pch = c(1, 2, 4), cex = 1)
dev.off()

savepdf("total_dim1_temporal_order1", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(x = kmeans_total_plot_long_run$plot_env$data$x, y = total_TL_coef[2,], xlab = "", ylab = expression(hat(b)[3]), cex.lab = 1.5,  ylim = c(-2.2, 2.2))
points(x = kmeans_total_plot_long_run$plot_env$data$x, y = total_TL_hampel_coef[2,], col = 2, pch = 2)
points(x = kmeans_total_plot_long_run$plot_env$data$x, y = total_TL_bisquare_coef[2,], col = 4, pch = 4)
dev.off()

savepdf("total_dim1_temporal_order2", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(x = kmeans_total_plot_long_run$plot_env$data$x, y = total_TL_coef[3,], xlab = "", ylab = expression(hat(c)[3]), cex.lab = 1.5, ylim = c(-2.55, -0.01))
points(x = kmeans_total_plot_long_run$plot_env$data$x, y = total_TL_hampel_coef[3,], col = 2, pch = 2)
points(x = kmeans_total_plot_long_run$plot_env$data$x, y = total_TL_bisquare_coef[3,], col = 4, pch = 4)
dev.off()

savepdf("total_dim1_temporal_order3", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(x = kmeans_total_plot_long_run$plot_env$data$x, y = total_TL_coef[4,], xlab = "Dimension 1", ylab = expression(hat(d)[3]), cex.lab = 1.5, ylim = c(-0.405, -0.05), xlim = c(-5.7,6.5))
points(x = kmeans_total_plot_long_run$plot_env$data$x, y = total_TL_hampel_coef[4,], col = 2, pch = 2)
points(x = kmeans_total_plot_long_run$plot_env$data$x, y = total_TL_bisquare_coef[4,], col = 4, pch = 4)
text_ind = c(which.max(total_TL_coef[4,]), which.min(total_TL_coef[4,]))
text(kmeans_total_plot_long_run$plot_env$data$x[text_ind], total_TL_coef[4,text_ind], labels = state[-1][text_ind], cex = 0.8, pos = 3)
dev.off()

## female
savepdf("female_dim1_temporal_intercept", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(x = kmeans_female_plot_long_run$plot_env$data$x, y = female_TL_coef[1,], xlab = "", ylab = "", cex.lab = 1.5, ylim = c(-4.3, -2.9), xlim = c(-5.6, 3.1))
points(x = kmeans_female_plot_long_run$plot_env$data$x, y = female_TL_hampel_coef[1,], col = 2, pch = 2)
points(x = kmeans_female_plot_long_run$plot_env$data$x, y = female_TL_bisquare_coef[1,], col = 4, pch = 4)
dev.off()

savepdf("female_dim1_temporal_order1", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(x = kmeans_female_plot_long_run$plot_env$data$x, y = female_TL_coef[2,], xlab = "", ylab = "", cex.lab = 1.5,  ylim = c(-2.7, 0.2))
points(x = kmeans_female_plot_long_run$plot_env$data$x, y = female_TL_hampel_coef[2,], col = 2, pch = 2)
points(x = kmeans_female_plot_long_run$plot_env$data$x, y = female_TL_bisquare_coef[2,], col = 4, pch = 4)
dev.off()

savepdf("female_dim1_temporal_order2", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(x = kmeans_female_plot_long_run$plot_env$data$x, y = female_TL_coef[3,], xlab = "", ylab = "", cex.lab = 1.5, ylim = c(-2.55, -1))
points(x = kmeans_female_plot_long_run$plot_env$data$x, y = female_TL_hampel_coef[3,], col = 2, pch = 2)
points(x = kmeans_female_plot_long_run$plot_env$data$x, y = female_TL_bisquare_coef[3,], col = 4, pch = 4)
dev.off()

savepdf("female_dim1_temporal_order3", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(x = kmeans_female_plot_long_run$plot_env$data$x, y = female_TL_coef[4,], xlab = "Dimension 1", ylab = "", cex.lab = 1.5, ylim = c(-0.405, -0.17), xlim = c(-5.8,3.5))
points(x = kmeans_female_plot_long_run$plot_env$data$x, y = female_TL_hampel_coef[4,], col = 2, pch = 2)
points(x = kmeans_female_plot_long_run$plot_env$data$x, y = female_TL_bisquare_coef[4,], col = 4, pch = 4)
text_ind = c(which.max(female_TL_bisquare_coef[4,]),which.min(female_TL_bisquare_coef[4,]))
text(kmeans_female_plot_long_run$plot_env$data$x[text_ind], female_TL_coef[4,text_ind], labels = state[-1][text_ind], cex = 0.8, pos = 3)
dev.off()

## male
savepdf("male_dim1_temporal_intercept", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(x = kmeans_male_plot_long_run$plot_env$data$x, y = male_TL_coef[1,], xlab = "", ylab = "", cex.lab = 1.5, ylim = c(-4.5,-2))
points(x = kmeans_male_plot_long_run$plot_env$data$x, y = male_TL_hampel_coef[1,], col = 2, pch = 2)
points(x = kmeans_male_plot_long_run$plot_env$data$x, y = male_TL_bisquare_coef[1,], col = 4, pch = 4)
dev.off()

savepdf("male_dim1_temporal_order1", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(x = kmeans_male_plot_long_run$plot_env$data$x, y = male_TL_coef[2,], xlab = "", ylab = "", cex.lab = 1.5,  ylim = c(-2.7, 3))
points(x = kmeans_male_plot_long_run$plot_env$data$x, y = male_TL_hampel_coef[2,], col = 2, pch = 2)
points(x = kmeans_male_plot_long_run$plot_env$data$x, y = male_TL_bisquare_coef[2,], col = 4, pch = 4)
dev.off()

savepdf("male_dim1_temporal_order2", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(x = kmeans_male_plot_long_run$plot_env$data$x, y = male_TL_coef[3,], xlab = "", ylab = "", cex.lab = 1.5, ylim = c(-2.3, 0.4))
points(x = kmeans_male_plot_long_run$plot_env$data$x, y = male_TL_hampel_coef[3,], col = 2, pch = 2)
points(x = kmeans_male_plot_long_run$plot_env$data$x, y = male_TL_bisquare_coef[3,], col = 4, pch = 4)
dev.off()

savepdf("male_dim1_temporal_order3", width = 12, height = 10, toplines = 0.8)
par(mar = c(4,4.5,2,2))
plot(x = kmeans_male_plot_long_run$plot_env$data$x, y = male_TL_coef[4,], xlab = "Dimension 1", ylab = "", cex.lab = 1.5, ylim = c(-0.35, 0.03), xlim = c(-5.5, 5.2))
points(x = kmeans_male_plot_long_run$plot_env$data$x, y = male_TL_hampel_coef[4,], col = 2, pch = 2)
points(x = kmeans_male_plot_long_run$plot_env$data$x, y = male_TL_bisquare_coef[4,], col = 4, pch = 4)
text_ind = c(which.max(male_TL_bisquare_coef[4,]), sort(male_TL_bisquare_coef[4,], decreasing = F, index.return = T)$ix[1:2])
text(kmeans_male_plot_long_run$plot_env$data$x[text_ind], male_TL_coef[4,text_ind], labels = state[-1][text_ind], cex = 0.8, pos = 3)
dev.off()

