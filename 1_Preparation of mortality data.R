###################################################
# read data from Japanese Human Mortality Database
###################################################

read.jpn_death <- function(region, label = region)
{
    path <- paste("http://www.ipss.go.jp/p-toukei/JMD/", region, "/STATS/",   "Deaths_1x1.txt", sep = "")
    txt <- RCurl::getURL(path)
    con <- textConnection(txt)
    mx <- try(read.table(con, skip = 2, header = TRUE, na.strings = "."), TRUE)
    close(con)
    if(class(mx) == "try-error") 
       stop("Connection error at www.mortality.org. Please check username, password and country label.")
    path <- paste("http://www.ipss.go.jp/p-toukei/JMD/", region, "/STATS/",   "Exposures_1x1.txt", sep = "")
    txt <- RCurl::getURL(path)
    con <- textConnection(txt)
    pop <- try(read.table(con, skip = 2, header = TRUE, na.strings = "."), TRUE)
    close(con)
    if(class(pop) == "try-error") 
        stop("Exposures file not found at www.mortality.org")
    return(list(Deaths = mx, pop = pop))
}

read.jpn <- function (region,  label = region) 
{
    path <- paste("http://www.ipss.go.jp/p-toukei/JMD/", region, "/STATS/",   "Mx_1x1.txt", sep = "")
    txt <- RCurl::getURL(path)
    con <- textConnection(txt)
    mx <- try(read.table(con, skip = 2, header = TRUE, na.strings = "."), 
              TRUE)
    close(con)
    if (class(mx) == "try-error") 
        stop("Connection error at www.mortality.org. Please check username, password and country label.")
    path <- paste("http://www.ipss.go.jp/p-toukei/JMD/", region, "/STATS/",   "Exposures_1x1.txt", sep = "")
    txt <- RCurl::getURL(path)
    con <- textConnection(txt)
    pop <- try(read.table(con, skip = 2, header = TRUE, na.strings = "."), 
               TRUE)
    close(con)
    if (class(pop) == "try-error") 
        stop("Exposures file not found at www.mortality.org")
    obj <- list(type = "mortality", label = label, lambda = 0)
    obj$year <- sort(unique(mx[, 1]))
    n <- length(obj$year)
    m <- length(unique(mx[, 2]))
    obj$age <- mx[1:m, 2]
    mnames <- names(mx)[-c(1, 2)]
    n.mort <- length(mnames)
    obj$rate <- obj$pop <- list()
    for (i in 1:n.mort) {
        obj$rate[[i]] <- matrix(mx[, i + 2], nrow = m, ncol = n)
        obj$rate[[i]][obj$rate[[i]] < 0] <- NA
        obj$pop[[i]] <- matrix(pop[, i + 2], nrow = m, ncol = n)
        obj$pop[[i]][obj$pop[[i]] < 0] <- NA
        dimnames(obj$rate[[i]]) <- dimnames(obj$pop[[i]]) <- list(obj$age, 
                                                                  obj$year)
    }
    names(obj$pop) = (names(obj$rate) <- tolower(mnames))
    obj$age <- as.numeric(as.character(obj$age))
    if (is.na(obj$age[m])) 
        obj$age[m] <- 2 * obj$age[m - 1] - obj$age[m - 2]
    return(structure(obj, class = "demogdata"))
}

state = c("Japan", "Hokkaido", "Aomori", "Iwate", "Miyagi", "Akita", "Yamagata", "Fukushima",
          "Ibaraki", "Tochigi", "Gunma", "Saitama", "Chiba", "Tokyo", "Kanagawa", "Niigata",
          "Toyama", "Ishikawa", "Fukui", "Yamanashi", "Nagano", "Gifu", "Shizuoka", "Aichi",
          "Mie", "Shiga", "Kyoto", "Osaka", "Hyogo", "Nara", "Wakayama", "Tottori", "Shimane",
          "Okayama", "Hiroshima", "Yamaguchi", "Tokushima", "Kagawa", "Ehime", "Kochi",
          "Fukuoka", "Saga", "Nagasaki", "Kumamoto", "Oita", "Miyazaki", "Kagoshima", "Okinawa")

state_smooth = c("Japan_smooth", "Hokkaido_smooth", "Aomori_smooth", "Iwate_smooth", 
                 "Miyagi_smooth", "Akita_smooth", "Yamagata_smooth", "Fukushima_smooth",
                 "Ibaraki_smooth", "Tochigi_smooth", "Gunma_smooth", "Saitama_smooth", 
                 "Chiba_smooth", "Tokyo_smooth", "Kanagawa_smooth", "Niigata_smooth",
                 "Toyama_smooth", "Ishikawa_smooth", "Fukui_smooth", "Yamanashi_smooth", 
                 "Nagano_smooth", "Gifu_smooth", "Shizuoka_smooth", "Aichi_smooth",
                 "Mie_smooth", "Shiga_smooth", "Kyoto_smooth", "Osaka_smooth", "Hyogo_smooth", 
                 "Nara_smooth", "Wakayama_smooth", "Tottori_smooth", "Shimane_smooth",
                 "Okayama_smooth", "Hiroshima_smooth", "Yamaguchi_smooth", "Tokushima_smooth", 
                 "Kagawa_smooth", "Ehime_smooth", "Kochi_smooth", "Fukuoka_smooth", "Saga_smooth", 
                 "Nagasaki_smooth", "Kumamoto_smooth", "Oita_smooth", "Miyazaki_smooth", 
                 "Kagoshima_smooth", "Okinawa_smooth")


#################################################
# full raw data (1975 to 2013) for ages 0 to 100
#################################################

#######################################
# precise death counts (no repetition)
######################################

dum = read.jpn_death("00", "Japan")$Deaths
Japan_count_F = matrix(dum[3109:7770,3], nrow=111)
Japan_count_M = matrix(dum[3109:7770,4], nrow=111)
Japan_count_T = matrix(dum[3109:7770,5], nrow=111)

Japan_count_female = rbind(Japan_count_F[1:100,], colSums(Japan_count_F[101:111,]))
Japan_count_male   = rbind(Japan_count_M[1:100,], colSums(Japan_count_M[101:111,]))
Japan_count_total  = rbind(Japan_count_T[1:100,], colSums(Japan_count_T[101:111,]))
rm(dum)

# Hokkaido

dum = read.jpn_death("01", "Hokkaido")$Deaths
Hokkaido_count_F = matrix(dum[,3], nrow=111)
Hokkaido_count_M = matrix(dum[,4], nrow=111)
Hokkaido_count_T = matrix(dum[,5], nrow=111)

Hokkaido_count_female = rbind(Hokkaido_count_F[1:100,], colSums(Hokkaido_count_F[101:111,]))
Hokkaido_count_male   = rbind(Hokkaido_count_M[1:100,], colSums(Hokkaido_count_M[101:111,]))
Hokkaido_count_total  = rbind(Hokkaido_count_T[1:100,], colSums(Hokkaido_count_T[101:111,])) 
rm(dum)

# Aomori

dum = read.jpn_death("02", "Aomori")$Deaths
Aomori_count_F = matrix(dum[,3], nrow=111)
Aomori_count_M = matrix(dum[,4], nrow=111)
Aomori_count_T = matrix(dum[,5], nrow=111)

Aomori_count_female = rbind(Aomori_count_F[1:100,], colSums(Aomori_count_F[101:111,]))
Aomori_count_male   = rbind(Aomori_count_M[1:100,], colSums(Aomori_count_M[101:111,]))
Aomori_count_total  = rbind(Aomori_count_T[1:100,], colSums(Aomori_count_T[101:111,])) 
rm(dum)

# Iwate

dum = read.jpn_death("03", "Iwate")$Deaths
Iwate_count_F = matrix(dum[,3], nrow=111)
Iwate_count_M = matrix(dum[,4], nrow=111)
Iwate_count_T = matrix(dum[,5], nrow=111)

Iwate_count_female = rbind(Iwate_count_F[1:100,], colSums(Iwate_count_F[101:111,]))
Iwate_count_male   = rbind(Iwate_count_M[1:100,], colSums(Iwate_count_M[101:111,]))
Iwate_count_total  = rbind(Iwate_count_T[1:100,], colSums(Iwate_count_T[101:111,])) 
rm(dum)

# Miyagi

dum = read.jpn_death("04", "Miyagi")$Deaths
Miyagi_count_F = matrix(dum[,3], nrow=111)
Miyagi_count_M = matrix(dum[,4], nrow=111)
Miyagi_count_T = matrix(dum[,5], nrow=111)

Miyagi_count_female = rbind(Miyagi_count_F[1:100,], colSums(Miyagi_count_F[101:111,]))
Miyagi_count_male   = rbind(Miyagi_count_M[1:100,], colSums(Miyagi_count_M[101:111,]))
Miyagi_count_total  = rbind(Miyagi_count_T[1:100,], colSums(Miyagi_count_T[101:111,])) 
rm(dum)

# Akita

dum = read.jpn_death("05", "Akita")$Deaths
Akita_count_F = matrix(dum[,3], nrow=111)
Akita_count_M = matrix(dum[,4], nrow=111)
Akita_count_T = matrix(dum[,5], nrow=111)

Akita_count_female = rbind(Akita_count_F[1:100,], colSums(Akita_count_F[101:111,]))
Akita_count_male   = rbind(Akita_count_M[1:100,], colSums(Akita_count_M[101:111,]))
Akita_count_total  = rbind(Akita_count_T[1:100,], colSums(Akita_count_T[101:111,])) 
rm(dum)

# Yamagata

dum = read.jpn_death("06", "Yamagata")$Deaths
Yamagata_count_F = matrix(dum[,3], nrow=111)
Yamagata_count_M = matrix(dum[,4], nrow=111)
Yamagata_count_T = matrix(dum[,5], nrow=111)

Yamagata_count_female = rbind(Yamagata_count_F[1:100,], colSums(Yamagata_count_F[101:111,]))
Yamagata_count_male   = rbind(Yamagata_count_M[1:100,], colSums(Yamagata_count_M[101:111,]))
Yamagata_count_total  = rbind(Yamagata_count_T[1:100,], colSums(Yamagata_count_T[101:111,])) 
rm(dum)

# Fukushima

dum = read.jpn_death("07", "Fukushima")$Deaths
Fukushima_count_F = matrix(dum[,3], nrow=111)
Fukushima_count_M = matrix(dum[,4], nrow=111)
Fukushima_count_T = matrix(dum[,5], nrow=111)

Fukushima_count_female = rbind(Fukushima_count_F[1:100,], colSums(Fukushima_count_F[101:111,]))
Fukushima_count_male   = rbind(Fukushima_count_M[1:100,], colSums(Fukushima_count_M[101:111,]))
Fukushima_count_total  = rbind(Fukushima_count_T[1:100,], colSums(Fukushima_count_T[101:111,])) 
rm(dum)

# Ibaraki

dum = read.jpn_death("08", "Ibaraki")$Deaths
Ibaraki_count_F = matrix(dum[,3], nrow=111)
Ibaraki_count_M = matrix(dum[,4], nrow=111)
Ibaraki_count_T = matrix(dum[,5], nrow=111)

Ibaraki_count_female = rbind(Ibaraki_count_F[1:100,], colSums(Ibaraki_count_F[101:111,]))
Ibaraki_count_male   = rbind(Ibaraki_count_M[1:100,], colSums(Ibaraki_count_M[101:111,]))
Ibaraki_count_total  = rbind(Ibaraki_count_T[1:100,], colSums(Ibaraki_count_T[101:111,])) 
rm(dum)

# Tochigi

dum = read.jpn_death("09", "Tochigi")$Deaths
Tochigi_count_F = matrix(dum[,3], nrow=111)
Tochigi_count_M = matrix(dum[,4], nrow=111)
Tochigi_count_T = matrix(dum[,5], nrow=111)

Tochigi_count_female = rbind(Tochigi_count_F[1:100,], colSums(Tochigi_count_F[101:111,]))
Tochigi_count_male   = rbind(Tochigi_count_M[1:100,], colSums(Tochigi_count_M[101:111,]))
Tochigi_count_total  = rbind(Tochigi_count_T[1:100,], colSums(Tochigi_count_T[101:111,])) 
rm(dum)

# Gunma

dum = read.jpn_death("10", "Gunma")$Deaths
Gunma_count_F = matrix(dum[,3], nrow=111)
Gunma_count_M = matrix(dum[,4], nrow=111)
Gunma_count_T = matrix(dum[,5], nrow=111)

Gunma_count_female = rbind(Gunma_count_F[1:100,], colSums(Gunma_count_F[101:111,]))
Gunma_count_male   = rbind(Gunma_count_M[1:100,], colSums(Gunma_count_M[101:111,]))
Gunma_count_total  = rbind(Gunma_count_T[1:100,], colSums(Gunma_count_T[101:111,])) 
rm(dum)

# Saitama

dum = read.jpn_death("11", "Saitama")$Deaths
Saitama_count_F = matrix(dum[,3], nrow=111)
Saitama_count_M = matrix(dum[,4], nrow=111)
Saitama_count_T = matrix(dum[,5], nrow=111)

Saitama_count_female = rbind(Saitama_count_F[1:100,], colSums(Saitama_count_F[101:111,]))
Saitama_count_male   = rbind(Saitama_count_M[1:100,], colSums(Saitama_count_M[101:111,]))
Saitama_count_total  = rbind(Saitama_count_T[1:100,], colSums(Saitama_count_T[101:111,])) 
rm(dum)

# Chiba

dum = read.jpn_death("12", "Chiba")$Deaths
Chiba_count_F = matrix(dum[,3], nrow=111)
Chiba_count_M = matrix(dum[,4], nrow=111)
Chiba_count_T = matrix(dum[,5], nrow=111)

Chiba_count_female = rbind(Chiba_count_F[1:100,], colSums(Chiba_count_F[101:111,]))
Chiba_count_male   = rbind(Chiba_count_M[1:100,], colSums(Chiba_count_M[101:111,]))
Chiba_count_total  = rbind(Chiba_count_T[1:100,], colSums(Chiba_count_T[101:111,])) 
rm(dum)

# Tokyo

dum = read.jpn_death("13", "Tokyo")$Deaths
Tokyo_count_F = matrix(dum[,3], nrow=111)
Tokyo_count_M = matrix(dum[,4], nrow=111)
Tokyo_count_T = matrix(dum[,5], nrow=111)

Tokyo_count_female = rbind(Tokyo_count_F[1:100,], colSums(Tokyo_count_F[101:111,]))
Tokyo_count_male   = rbind(Tokyo_count_M[1:100,], colSums(Tokyo_count_M[101:111,]))
Tokyo_count_total  = rbind(Tokyo_count_T[1:100,], colSums(Tokyo_count_T[101:111,])) 
rm(dum)

# Kanagawa

dum = read.jpn_death("14", "Kanagawa")$Deaths
Kanagawa_count_F = matrix(dum[,3], nrow=111)
Kanagawa_count_M = matrix(dum[,4], nrow=111)
Kanagawa_count_T = matrix(dum[,5], nrow=111)

Kanagawa_count_female = rbind(Kanagawa_count_F[1:100,], colSums(Kanagawa_count_F[101:111,]))
Kanagawa_count_male   = rbind(Kanagawa_count_M[1:100,], colSums(Kanagawa_count_M[101:111,]))
Kanagawa_count_total  = rbind(Kanagawa_count_T[1:100,], colSums(Kanagawa_count_T[101:111,])) 
rm(dum)

# Niigata

dum = read.jpn_death("15", "Niigata")$Deaths
Niigata_count_F = matrix(dum[,3], nrow=111)
Niigata_count_M = matrix(dum[,4], nrow=111)
Niigata_count_T = matrix(dum[,5], nrow=111)

Niigata_count_female = rbind(Niigata_count_F[1:100,], colSums(Niigata_count_F[101:111,]))
Niigata_count_male   = rbind(Niigata_count_M[1:100,], colSums(Niigata_count_M[101:111,]))
Niigata_count_total  = rbind(Niigata_count_T[1:100,], colSums(Niigata_count_T[101:111,])) 
rm(dum)

# Toyama

dum = read.jpn_death("16", "Toyama")$Deaths
Toyama_count_F = matrix(dum[,3], nrow=111)
Toyama_count_M = matrix(dum[,4], nrow=111)
Toyama_count_T = matrix(dum[,5], nrow=111)

Toyama_count_female = rbind(Toyama_count_F[1:100,], colSums(Toyama_count_F[101:111,]))
Toyama_count_male   = rbind(Toyama_count_M[1:100,], colSums(Toyama_count_M[101:111,]))
Toyama_count_total  = rbind(Toyama_count_T[1:100,], colSums(Toyama_count_T[101:111,])) 
rm(dum)

# Ishikawa

dum = read.jpn_death("17", "Ishikawa")$Deaths
Ishikawa_count_F = matrix(dum[,3], nrow=111)
Ishikawa_count_M = matrix(dum[,4], nrow=111)
Ishikawa_count_T = matrix(dum[,5], nrow=111)

Ishikawa_count_female = rbind(Ishikawa_count_F[1:100,], colSums(Ishikawa_count_F[101:111,]))
Ishikawa_count_male   = rbind(Ishikawa_count_M[1:100,], colSums(Ishikawa_count_M[101:111,]))
Ishikawa_count_total  = rbind(Ishikawa_count_T[1:100,], colSums(Ishikawa_count_T[101:111,])) 
rm(dum)

# Fukui

dum = read.jpn_death("18", "Fukui")$Deaths
Fukui_count_F = matrix(dum[,3], nrow=111)
Fukui_count_M = matrix(dum[,4], nrow=111)
Fukui_count_T = matrix(dum[,5], nrow=111)

Fukui_count_female = rbind(Fukui_count_F[1:100,], colSums(Fukui_count_F[101:111,]))
Fukui_count_male   = rbind(Fukui_count_M[1:100,], colSums(Fukui_count_M[101:111,]))
Fukui_count_total  = rbind(Fukui_count_T[1:100,], colSums(Fukui_count_T[101:111,])) 
rm(dum)

# Yamanashi

dum = read.jpn_death("19", "Yamanashi")$Deaths
Yamanashi_count_F = matrix(dum[,3], nrow=111)
Yamanashi_count_M = matrix(dum[,4], nrow=111)
Yamanashi_count_T = matrix(dum[,5], nrow=111)

Yamanashi_count_female = rbind(Yamanashi_count_F[1:100,], colSums(Yamanashi_count_F[101:111,]))
Yamanashi_count_male   = rbind(Yamanashi_count_M[1:100,], colSums(Yamanashi_count_M[101:111,]))
Yamanashi_count_total  = rbind(Yamanashi_count_T[1:100,], colSums(Yamanashi_count_T[101:111,])) 
rm(dum)

# Nagano

dum = read.jpn_death("20", "Nagano")$Deaths
Nagano_count_F = matrix(dum[,3], nrow=111)
Nagano_count_M = matrix(dum[,4], nrow=111)
Nagano_count_T = matrix(dum[,5], nrow=111)

Nagano_count_female = rbind(Nagano_count_F[1:100,], colSums(Nagano_count_F[101:111,]))
Nagano_count_male   = rbind(Nagano_count_M[1:100,], colSums(Nagano_count_M[101:111,]))
Nagano_count_total  = rbind(Nagano_count_T[1:100,], colSums(Nagano_count_T[101:111,])) 
rm(dum)

# Gifu

dum = read.jpn_death("21", "Gifu")$Deaths
Gifu_count_F = matrix(dum[,3], nrow=111)
Gifu_count_M = matrix(dum[,4], nrow=111)
Gifu_count_T = matrix(dum[,5], nrow=111)

Gifu_count_female = rbind(Gifu_count_F[1:100,], colSums(Gifu_count_F[101:111,]))
Gifu_count_male   = rbind(Gifu_count_M[1:100,], colSums(Gifu_count_M[101:111,]))
Gifu_count_total  = rbind(Gifu_count_T[1:100,], colSums(Gifu_count_T[101:111,])) 
rm(dum)

# Shizuoka

dum = read.jpn_death("22", "Shizuoka")$Deaths
Shizuoka_count_F = matrix(dum[,3], nrow=111)
Shizuoka_count_M = matrix(dum[,4], nrow=111)
Shizuoka_count_T = matrix(dum[,5], nrow=111)

Shizuoka_count_female = rbind(Shizuoka_count_F[1:100,], colSums(Shizuoka_count_F[101:111,]))
Shizuoka_count_male   = rbind(Shizuoka_count_M[1:100,], colSums(Shizuoka_count_M[101:111,]))
Shizuoka_count_total  = rbind(Shizuoka_count_T[1:100,], colSums(Shizuoka_count_T[101:111,])) 
rm(dum)

# Aichi

dum = read.jpn_death("23", "Aichi")$Deaths
Aichi_count_F = matrix(dum[,3], nrow=111)
Aichi_count_M = matrix(dum[,4], nrow=111)
Aichi_count_T = matrix(dum[,5], nrow=111)

Aichi_count_female = rbind(Aichi_count_F[1:100,], colSums(Aichi_count_F[101:111,]))
Aichi_count_male   = rbind(Aichi_count_M[1:100,], colSums(Aichi_count_M[101:111,]))
Aichi_count_total  = rbind(Aichi_count_T[1:100,], colSums(Aichi_count_T[101:111,])) 
rm(dum)

# Mie

dum = read.jpn_death("24", "Mie")$Deaths
Mie_count_F = matrix(dum[,3], nrow=111)
Mie_count_M = matrix(dum[,4], nrow=111)
Mie_count_T = matrix(dum[,5], nrow=111)

Mie_count_female = rbind(Mie_count_F[1:100,], colSums(Mie_count_F[101:111,]))
Mie_count_male   = rbind(Mie_count_M[1:100,], colSums(Mie_count_M[101:111,]))
Mie_count_total  = rbind(Mie_count_T[1:100,], colSums(Mie_count_T[101:111,])) 
rm(dum)

# Shiga

dum = read.jpn_death("25", "Shiga")$Deaths
Shiga_count_F = matrix(dum[,3], nrow=111)
Shiga_count_M = matrix(dum[,4], nrow=111)
Shiga_count_T = matrix(dum[,5], nrow=111)

Shiga_count_female = rbind(Shiga_count_F[1:100,], colSums(Shiga_count_F[101:111,]))
Shiga_count_male   = rbind(Shiga_count_M[1:100,], colSums(Shiga_count_M[101:111,]))
Shiga_count_total  = rbind(Shiga_count_T[1:100,], colSums(Shiga_count_T[101:111,])) 
rm(dum)

# Kyoto

dum = read.jpn_death("26", "Kyoto")$Deaths
Kyoto_count_F = matrix(dum[,3], nrow=111)
Kyoto_count_M = matrix(dum[,4], nrow=111)
Kyoto_count_T = matrix(dum[,5], nrow=111)

Kyoto_count_female = rbind(Kyoto_count_F[1:100,], colSums(Kyoto_count_F[101:111,]))
Kyoto_count_male   = rbind(Kyoto_count_M[1:100,], colSums(Kyoto_count_M[101:111,]))
Kyoto_count_total  = rbind(Kyoto_count_T[1:100,], colSums(Kyoto_count_T[101:111,])) 
rm(dum)

# Osaka

dum = read.jpn_death("27", "Osaka")$Deaths
Osaka_count_F = matrix(dum[,3], nrow=111)
Osaka_count_M = matrix(dum[,4], nrow=111)
Osaka_count_T = matrix(dum[,5], nrow=111)

Osaka_count_female = rbind(Osaka_count_F[1:100,], colSums(Osaka_count_F[101:111,]))
Osaka_count_male   = rbind(Osaka_count_M[1:100,], colSums(Osaka_count_M[101:111,]))
Osaka_count_total  = rbind(Osaka_count_T[1:100,], colSums(Osaka_count_T[101:111,])) 
rm(dum)

# Hyogo

dum = read.jpn_death("28", "Hyogo")$Deaths
Hyogo_count_F = matrix(dum[,3], nrow=111)
Hyogo_count_M = matrix(dum[,4], nrow=111)
Hyogo_count_T = matrix(dum[,5], nrow=111)

Hyogo_count_female = rbind(Hyogo_count_F[1:100,], colSums(Hyogo_count_F[101:111,]))
Hyogo_count_male   = rbind(Hyogo_count_M[1:100,], colSums(Hyogo_count_M[101:111,]))
Hyogo_count_total  = rbind(Hyogo_count_T[1:100,], colSums(Hyogo_count_T[101:111,])) 
rm(dum)

# Nara

dum = read.jpn_death("29", "Nara")$Deaths
Nara_count_F = matrix(dum[,3], nrow=111)
Nara_count_M = matrix(dum[,4], nrow=111)
Nara_count_T = matrix(dum[,5], nrow=111)

Nara_count_female = rbind(Nara_count_F[1:100,], colSums(Nara_count_F[101:111,]))
Nara_count_male   = rbind(Nara_count_M[1:100,], colSums(Nara_count_M[101:111,]))
Nara_count_total  = rbind(Nara_count_T[1:100,], colSums(Nara_count_T[101:111,])) 
rm(dum)

# Wakayama

dum = read.jpn_death("30", "Wakayama")$Deaths
Wakayama_count_F = matrix(dum[,3], nrow=111)
Wakayama_count_M = matrix(dum[,4], nrow=111)
Wakayama_count_T = matrix(dum[,5], nrow=111)

Wakayama_count_female = rbind(Wakayama_count_F[1:100,], colSums(Wakayama_count_F[101:111,]))
Wakayama_count_male   = rbind(Wakayama_count_M[1:100,], colSums(Wakayama_count_M[101:111,]))
Wakayama_count_total  = rbind(Wakayama_count_T[1:100,], colSums(Wakayama_count_T[101:111,])) 
rm(dum)

# Tottori

dum = read.jpn_death("31", "Tottori")$Deaths
Tottori_count_F = matrix(dum[,3], nrow=111)
Tottori_count_M = matrix(dum[,4], nrow=111)
Tottori_count_T = matrix(dum[,5], nrow=111)

Tottori_count_female = rbind(Tottori_count_F[1:100,], colSums(Tottori_count_F[101:111,]))
Tottori_count_male   = rbind(Tottori_count_M[1:100,], colSums(Tottori_count_M[101:111,]))
Tottori_count_total  = rbind(Tottori_count_T[1:100,], colSums(Tottori_count_T[101:111,])) 
rm(dum)

# Shimane

dum = read.jpn_death("32", "Shimane")$Deaths
Shimane_count_F = matrix(dum[,3], nrow=111)
Shimane_count_M = matrix(dum[,4], nrow=111)
Shimane_count_T = matrix(dum[,5], nrow=111)

Shimane_count_female = rbind(Shimane_count_F[1:100,], colSums(Shimane_count_F[101:111,]))
Shimane_count_male   = rbind(Shimane_count_M[1:100,], colSums(Shimane_count_M[101:111,]))
Shimane_count_total  = rbind(Shimane_count_T[1:100,], colSums(Shimane_count_T[101:111,])) 
rm(dum)

# Okayama

dum = read.jpn_death("33", "Okayama")$Deaths
Okayama_count_F = matrix(dum[,3], nrow=111)
Okayama_count_M = matrix(dum[,4], nrow=111)
Okayama_count_T = matrix(dum[,5], nrow=111)

Okayama_count_female = rbind(Okayama_count_F[1:100,], colSums(Okayama_count_F[101:111,]))
Okayama_count_male   = rbind(Okayama_count_M[1:100,], colSums(Okayama_count_M[101:111,]))
Okayama_count_total  = rbind(Okayama_count_T[1:100,], colSums(Okayama_count_T[101:111,])) 
rm(dum)

# Hiroshima

dum = read.jpn_death("34", "Hiroshima")$Deaths
Hiroshima_count_F = matrix(dum[,3], nrow=111)
Hiroshima_count_M = matrix(dum[,4], nrow=111)
Hiroshima_count_T = matrix(dum[,5], nrow=111)

Hiroshima_count_female = rbind(Hiroshima_count_F[1:100,], colSums(Hiroshima_count_F[101:111,]))
Hiroshima_count_male   = rbind(Hiroshima_count_M[1:100,], colSums(Hiroshima_count_M[101:111,]))
Hiroshima_count_total  = rbind(Hiroshima_count_T[1:100,], colSums(Hiroshima_count_T[101:111,])) 
rm(dum)

# Yamaguchi

dum = read.jpn_death("35", "Yamaguchi")$Deaths
Yamaguchi_count_F = matrix(dum[,3], nrow=111)
Yamaguchi_count_M = matrix(dum[,4], nrow=111)
Yamaguchi_count_T = matrix(dum[,5], nrow=111)

Yamaguchi_count_female = rbind(Yamaguchi_count_F[1:100,], colSums(Yamaguchi_count_F[101:111,]))
Yamaguchi_count_male   = rbind(Yamaguchi_count_M[1:100,], colSums(Yamaguchi_count_M[101:111,]))
Yamaguchi_count_total  = rbind(Yamaguchi_count_T[1:100,], colSums(Yamaguchi_count_T[101:111,])) 
rm(dum)

# Tokushima

dum = read.jpn_death("36", "Tokushima")$Deaths
Tokushima_count_F = matrix(dum[,3], nrow=111)
Tokushima_count_M = matrix(dum[,4], nrow=111)
Tokushima_count_T = matrix(dum[,5], nrow=111)

Tokushima_count_female = rbind(Tokushima_count_F[1:100,], colSums(Tokushima_count_F[101:111,]))
Tokushima_count_male   = rbind(Tokushima_count_M[1:100,], colSums(Tokushima_count_M[101:111,]))
Tokushima_count_total  = rbind(Tokushima_count_T[1:100,], colSums(Tokushima_count_T[101:111,])) 
rm(dum)

# Kagawa

dum = read.jpn_death("37", "Kagawa")$Deaths
Kagawa_count_F = matrix(dum[,3], nrow=111)
Kagawa_count_M = matrix(dum[,4], nrow=111)
Kagawa_count_T = matrix(dum[,5], nrow=111)

Kagawa_count_female = rbind(Kagawa_count_F[1:100,], colSums(Kagawa_count_F[101:111,]))
Kagawa_count_male   = rbind(Kagawa_count_M[1:100,], colSums(Kagawa_count_M[101:111,]))
Kagawa_count_total  = rbind(Kagawa_count_T[1:100,], colSums(Kagawa_count_T[101:111,])) 
rm(dum)

# Ehime

dum = read.jpn_death("38", "Ehime")$Deaths
Ehime_count_F = matrix(dum[,3], nrow=111)
Ehime_count_M = matrix(dum[,4], nrow=111)
Ehime_count_T = matrix(dum[,5], nrow=111)

Ehime_count_female = rbind(Ehime_count_F[1:100,], colSums(Ehime_count_F[101:111,]))
Ehime_count_male   = rbind(Ehime_count_M[1:100,], colSums(Ehime_count_M[101:111,]))
Ehime_count_total  = rbind(Ehime_count_T[1:100,], colSums(Ehime_count_T[101:111,])) 
rm(dum)

# Kochi

dum = read.jpn_death("39", "Kochi")$Deaths
Kochi_count_F = matrix(dum[,3], nrow=111)
Kochi_count_M = matrix(dum[,4], nrow=111)
Kochi_count_T = matrix(dum[,5], nrow=111)

Kochi_count_female = rbind(Kochi_count_F[1:100,], colSums(Kochi_count_F[101:111,]))
Kochi_count_male   = rbind(Kochi_count_M[1:100,], colSums(Kochi_count_M[101:111,]))
Kochi_count_total  = rbind(Kochi_count_T[1:100,], colSums(Kochi_count_T[101:111,])) 
rm(dum)

# Fukuoka

dum = read.jpn_death("40", "Fukuoka")$Deaths
Fukuoka_count_F = matrix(dum[,3], nrow=111)
Fukuoka_count_M = matrix(dum[,4], nrow=111)
Fukuoka_count_T = matrix(dum[,5], nrow=111)

Fukuoka_count_female = rbind(Fukuoka_count_F[1:100,], colSums(Fukuoka_count_F[101:111,]))
Fukuoka_count_male   = rbind(Fukuoka_count_M[1:100,], colSums(Fukuoka_count_M[101:111,]))
Fukuoka_count_total  = rbind(Fukuoka_count_T[1:100,], colSums(Fukuoka_count_T[101:111,])) 
rm(dum)

# Saga

dum = read.jpn_death("41", "Saga")$Deaths
Saga_count_F = matrix(dum[,3], nrow=111)
Saga_count_M = matrix(dum[,4], nrow=111)
Saga_count_T = matrix(dum[,5], nrow=111)

Saga_count_female = rbind(Saga_count_F[1:100,], colSums(Saga_count_F[101:111,]))
Saga_count_male   = rbind(Saga_count_M[1:100,], colSums(Saga_count_M[101:111,]))
Saga_count_total  = rbind(Saga_count_T[1:100,], colSums(Saga_count_T[101:111,])) 
rm(dum)

# Nagasaki

dum = read.jpn_death("42", "Nagasaki")$Deaths
Nagasaki_count_F = matrix(dum[,3], nrow=111)
Nagasaki_count_M = matrix(dum[,4], nrow=111)
Nagasaki_count_T = matrix(dum[,5], nrow=111)

Nagasaki_count_female = rbind(Nagasaki_count_F[1:100,], colSums(Nagasaki_count_F[101:111,]))
Nagasaki_count_male   = rbind(Nagasaki_count_M[1:100,], colSums(Nagasaki_count_M[101:111,]))
Nagasaki_count_total  = rbind(Nagasaki_count_T[1:100,], colSums(Nagasaki_count_T[101:111,])) 
rm(dum)

# Kumamoto

dum = read.jpn_death("43", "Kumamoto")$Deaths
Kumamoto_count_F = matrix(dum[,3], nrow=111)
Kumamoto_count_M = matrix(dum[,4], nrow=111)
Kumamoto_count_T = matrix(dum[,5], nrow=111)

Kumamoto_count_female = rbind(Kumamoto_count_F[1:100,], colSums(Kumamoto_count_F[101:111,]))
Kumamoto_count_male   = rbind(Kumamoto_count_M[1:100,], colSums(Kumamoto_count_M[101:111,]))
Kumamoto_count_total  = rbind(Kumamoto_count_T[1:100,], colSums(Kumamoto_count_T[101:111,])) 
rm(dum)

# Oita

dum = read.jpn_death("44", "Oita")$Deaths
Oita_count_F = matrix(dum[,3], nrow=111)
Oita_count_M = matrix(dum[,4], nrow=111)
Oita_count_T = matrix(dum[,5], nrow=111)

Oita_count_female = rbind(Oita_count_F[1:100,], colSums(Oita_count_F[101:111,]))
Oita_count_male   = rbind(Oita_count_M[1:100,], colSums(Oita_count_M[101:111,]))
Oita_count_total  = rbind(Oita_count_T[1:100,], colSums(Oita_count_T[101:111,])) 
rm(dum)

# Miyazaki

dum = read.jpn_death("45", "Miyazaki")$Deaths
Miyazaki_count_F = matrix(dum[,3], nrow=111)
Miyazaki_count_M = matrix(dum[,4], nrow=111)
Miyazaki_count_T = matrix(dum[,5], nrow=111)

Miyazaki_count_female = rbind(Miyazaki_count_F[1:100,], colSums(Miyazaki_count_F[101:111,]))
Miyazaki_count_male   = rbind(Miyazaki_count_M[1:100,], colSums(Miyazaki_count_M[101:111,]))
Miyazaki_count_total  = rbind(Miyazaki_count_T[1:100,], colSums(Miyazaki_count_T[101:111,])) 
rm(dum)

# Kagoshima

dum = read.jpn_death("46", "Kagoshima")$Deaths
Kagoshima_count_F = matrix(dum[,3], nrow=111)
Kagoshima_count_M = matrix(dum[,4], nrow=111)
Kagoshima_count_T = matrix(dum[,5], nrow=111)

Kagoshima_count_female = rbind(Kagoshima_count_F[1:100,], colSums(Kagoshima_count_F[101:111,]))
Kagoshima_count_male   = rbind(Kagoshima_count_M[1:100,], colSums(Kagoshima_count_M[101:111,]))
Kagoshima_count_total  = rbind(Kagoshima_count_T[1:100,], colSums(Kagoshima_count_T[101:111,])) 
rm(dum)

# Okinawa

dum = read.jpn_death("47", "Okinawa")$Deaths
Okinawa_count_F = matrix(dum[,3], nrow=111)
Okinawa_count_M = matrix(dum[,4], nrow=111)
Okinawa_count_T = matrix(dum[,5], nrow=111)

Okinawa_count_female = rbind(Okinawa_count_F[1:100,], colSums(Okinawa_count_F[101:111,]))
Okinawa_count_male   = rbind(Okinawa_count_M[1:100,], colSums(Okinawa_count_M[101:111,]))
Okinawa_count_total  = rbind(Okinawa_count_T[1:100,], colSums(Okinawa_count_T[101:111,])) 
rm(dum)


# rate

Japan = extract.years(extract.ages(read.jpn("00", "Japan"), 0:100), 1975:2018)
Hokkaido = extract.ages(read.jpn("01", "Hokkaido"), 0:100)
Aomori = extract.ages(read.jpn("02", "Aomori"), 0:100)
Iwate = extract.ages(read.jpn("03", "Iwate"), 0:100)
Miyagi = extract.ages(read.jpn("04", "Miyagi"), 0:100)
Akita = extract.ages(read.jpn("05", "Akita"), 0:100)
Yamagata = extract.ages(read.jpn("06", "Yamagata"), 0:100)
Fukushima = extract.ages(read.jpn("07", "Fukushima"), 0:100)
Ibaraki = extract.ages(read.jpn("08", "Ibaraki"), 0:100)
Tochigi = extract.ages(read.jpn("09", "Tochigi"), 0:100)
Gunma = extract.ages(read.jpn("10", "Gunma"), 0:100)
Saitama = extract.ages(read.jpn("11", "Saitama"), 0:100)
Chiba = extract.ages(read.jpn("12", "Chiba"), 0:100)
Tokyo = extract.ages(read.jpn("13", "Tokyo"), 0:100)
Kanagawa = extract.ages(read.jpn("14", "Kanagawa"), 0:100)
Niigata = extract.ages(read.jpn("15", "Niigata"), 0:100)
Toyama = extract.ages(read.jpn("16", "Toyama"), 0:100)
Ishikawa = extract.ages(read.jpn("17", "Ishikawa"), 0:100)
Fukui = extract.ages(read.jpn("18", "Fukui"), 0:100)
Yamanashi = extract.ages(read.jpn("19", "Yamanashi"), 0:100)
Nagano = extract.ages(read.jpn("20", "Nagano"), 0:100)
Gifu = extract.ages(read.jpn("21", "Gifu"), 0:100)
Shizuoka = extract.ages(read.jpn("22", "Shizuoka"), 0:100)
Aichi = extract.ages(read.jpn("23", "Aichi"), 0:100)
Mie = extract.ages(read.jpn("24", "Mie"), 0:100)
Shiga = extract.ages(read.jpn("25", "Shiga"), 0:100)
Kyoto = extract.ages(read.jpn("26", "Kyoto"), 0:100)
Osaka = extract.ages(read.jpn("27", "Osaka"), 0:100)
Hyogo = extract.ages(read.jpn("28", "Hyogo"), 0:100)
Nara = extract.ages(read.jpn("29", "Nara"), 0:100)
Wakayama = extract.ages(read.jpn("30", "Wakayama"), 0:100)
Tottori = extract.ages(read.jpn("31", "Tottori"), 0:100)
Shimane = extract.ages(read.jpn("32", "Shimane"), 0:100)
Okayama = extract.ages(read.jpn("33", "Okayama"), 0:100)
Hiroshima = extract.ages(read.jpn("34", "Hiroshima"), 0:100)
Yamaguchi = extract.ages(read.jpn("35", "Yamaguchi"), 0:100)
Tokushima = extract.ages(read.jpn("36", "Tokushima"), 0:100)
Kagawa    = extract.ages(read.jpn("37", "Kagawa"), 0:100)
Ehime     = extract.ages(read.jpn("38", "Ehime"), 0:100)
Kochi     = extract.ages(read.jpn("39", "Kochi"), 0:100)
Fukuoka   = extract.ages(read.jpn("40", "Fukuoka"), 0:100)
Saga      = extract.ages(read.jpn("41", "Saga"), 0:100)
Nagasaki  = extract.ages(read.jpn("42", "Nagasaki"), 0:100)
Kumamoto  = extract.ages(read.jpn("43", "Kumamoto"), 0:100)
Oita      = extract.ages(read.jpn("44", "Oita"), 0:100)
Miyazaki  = extract.ages(read.jpn("45", "Miyazaki"), 0:100)
Kagoshima = extract.ages(read.jpn("46", "Kagoshima"), 0:100)
Okinawa   = extract.ages(read.jpn("47", "Okinawa"), 0:100)

# check if the last year is 2017 for all states

year_store = vector("numeric",47)
for(ik in 1:47)
{
    year_store[ik] = tail(get(state[ik])$year, 1)
}
all(year_store == 2018)

# smoothed functional curves using penalized regression spline with monotonic constraint

for(ik in 1:48)
{
    assign(state_smooth[ik], smooth.demogdata(get(state[ik])))
    print(ik)
}

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

