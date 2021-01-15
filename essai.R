library(class)
library(FNN)
file.path("eknn_uncertain_label.R")
file.path("my_lib_3.R")

#daily data

oid <- read.csv("owid-covid-data.csv")
min <- c(0.1,0.5,0.65)
max <- c(0.25,0.7,0.9)
n <- c(3,4,5,6)
l <- c(1,2,3,4)
test_size <- 14
start <- 0.6
end <- 1
history_from <- 3
history_to <- 3
remove <- 1
target <- 3
optimize <- T
clean_data <- oid[which(oid$location == "France"),]
clean_data <- clean_data[,c(4,6,9)]
clean_data <- clean_data[complete.cases(clean_data),]

final_eval(clean_data,min,max,l,start,end,test_size,n,history_from,history_to,remove,target,optimize)





