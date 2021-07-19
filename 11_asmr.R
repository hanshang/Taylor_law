##########################################
# Japanese Mortality Database (1975-2018)
##########################################

Japan_extract = extract.years(Japan, 1975:2018)

age = 0:100
output_Japan = cbind(rep("Japan", 101*44), rep(1975:2018, each = 101), rep(age, 44),
                      as.numeric(Japan_extract$rate$female),
                    	as.numeric(Japan_extract$rate$male),
                    	as.numeric(Japan_extract$rate$total),
                    	as.numeric(Japan_extract$pop$female),
                    	as.numeric(Japan_extract$pop$male),
                    	as.numeric(Japan_extract$pop$total),
                      as.numeric(Japan_smooth$rate$female),
                      as.numeric(Japan_smooth$rate$male),
                      as.numeric(Japan_smooth$rate$total))

output_state = c("output_Hokkaido", "output_Aomori", "output_Iwate", "output_Miyagi", "output_Akita",
                 "output_Yamagata", "output_Fukushima", "output_Ibaraki", "output_Tochigi", "output_Gunma",
                 "output_Saitama",  "output_Chiba", "output_Tokyo", "output_Kanagawa",
                 "output_Niigata", "output_Toyama", "output_Ishikawa",
                 "output_Fukui", "output_Yamanashi", "output_Nagano", "output_Gifu",
                 "output_Shizuoka", "output_Aichi", "output_Mie", "output_Shiga", "output_Kyoto",
                 "output_Osaka",  "output_Hyogo", "output_Nara",
                 "output_Wakayama", "output_Tottori", "output_Shimane", "output_Okayama",
                 "output_Hiroshima", "output_Yamaguchi",
                 "output_Tokushima", "output_Kagawa", "output_Ehime", "output_Kochi",
                 "output_Fukuoka", "output_Saga",
                 "output_Nagasaki", "output_Kumamoto", "output_Oita", "output_Miyazaki",
                 "output_Kagoshima", "output_Okinawa")

for(iw in 1:47)
{
    assign(output_state[iw], cbind(rep(state[iw+1], 101*44), rep(1975:2018, each = 101), rep(age, 44),
                        as.numeric(get(state[iw+1])$rate$female),
                        as.numeric(get(state[iw+1])$rate$male),
                        as.numeric(get(state[iw+1])$rate$total),
                        as.numeric(get(state[iw+1])$pop$female),
                        as.numeric(get(state[iw+1])$pop$male),
                        as.numeric(get(state[iw+1])$pop$total),
                        as.numeric(get(state_smooth[iw+1])$rate$female),
                        as.numeric(get(state_smooth[iw+1])$rate$male),
                        as.numeric(get(state_smooth[iw+1])$rate$total)))
    print(iw); rm(iw)
}

output_comb = NULL
for(iw in 1:47)
{
    output_comb = rbind(output_comb, get(output_state[iw]))
    print(iw); rm(iw)
}
output_comb_final = rbind(output_Japan, output_comb)
colnames(output_comb_final) = c("Code", "Year", "Age", "MortFemale", "MortMale", "MortTotal", "ExpoFemale", "ExpoMale", "ExpoTotal",
                                "SmoothmortFemale", "SmoothmortMale", "SmoothmortTotal")

# export data in .csv format

write.csv(output_comb_final, "asmr.csv", quote=FALSE, row.names=FALSE)

# labels

cn0 = list()
cn0$Japan = noquote("Japan")
cn0$Hokkaido = noquote("Hokkaido")
cn0$Aomori = noquote("Aomori")
cn0$Iwate = noquote("Iwate")
cn0$Miyagi = noquote("Miyagi")
cn0$Akita = noquote("Akita")
cn0$Yamagata = noquote("Yamagata")
cn0$Fukushima = noquote("Fukushima")
cn0$Ibaraki = noquote("Ibaraki")
cn0$Tochigi = noquote("Tochigi")
cn0$Gunma = noquote("Gunma")
cn0$Saitama = noquote("Saitama")
cn0$Chiba = noquote("Chiba")
cn0$Tokyo = noquote("Tokyo")
cn0$Kanagawa = noquote("Kanagawa")
cn0$Niigata = noquote("Niigata")
cn0$Toyama = noquote("Toyama")
cn0$Ishikawa = noquote("Ishikawa")
cn0$Fukui = noquote("Fukui")
cn0$Yamanashi = noquote("Yamanashi")
cn0$Nagano = noquote("Nagano")
cn0$Gifu = noquote("Gifu")
cn0$Shizuoka = noquote("Shizuoka")
cn0$Aichi = noquote("Aichi")
cn0$Mie = noquote("Mie")
cn0$Shiga = noquote("Shiga")
cn0$Kyoto = noquote("Kyoto")
cn0$Osaka = noquote("Osaka")
cn0$Hyogo = noquote("Hyogo")
cn0$Nara = noquote("Nara")
cn0$Wakayama = noquote("Wakayama")
cn0$Tottori = noquote("Tottori")
cn0$Shimane = noquote("Shimane")
cn0$Okayama = noquote("Okayama")
cn0$Hiroshima = noquote("Hiroshima")
cn0$Yamaguchi = noquote("Yamaguchi")
cn0$Tokushima = noquote("Tokushima")
cn0$Kagawa = noquote("Kagawa")
cn0$Ehime = noquote("Ehime")
cn0$Kochi = noquote("Kochi")
cn0$Fukuoka = noquote("Fukuoka")
cn0$Saga = noquote("Saga")
cn0$Nagasaki = noquote("Nagasaki")
cn0$Kumamoto = noquote("Kumamoto")
cn0$Oita = noquote("Oita")
cn0$Miyazaki = noquote("Miyazaki")
cn0$Kagoshima = noquote("Kagoshima")
