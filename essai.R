library(class)
library(FNN)
file.path("eknn_uncertain_label.R")
file.path("my_lib_2.R")


#daily data

oid <- read.csv("owid-covid-data.csv")
min <- c(0.1,0.5,0.65)
max <- c(0.25,0.7,0.9)
clean_data <- oid[which(oid$location == "United Kingdom"),]
clean_data <- clean_data[,c(4,6,9)]
clean_data <- clean_data[complete.cases(clean_data),]

final_eval(clean_data,min,max,l,0.6,0.9,n,2,3,1,3,T)


#weekly data


min <- c(0.1,0.5,0.65)
max <- c(0.25,0.7,0.9)
clean_data <- oid[which(oid$location == "United Kingdom"),]
clean_data <- clean_data[,c(4,6,9)]
clean_data <- clean_data[complete.cases(clean_data),]
clean_data <- week_transform(clean_data)

final_eval(clean_data,min,max,l,0.6,0.9,n,2,3,1,3,T)




