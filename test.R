source("eknn_uncertain_regression.R")
source("evaluation_lib.R")
source("base_imputation.R")



oid <- read.csv("owid-covid-data.csv")
clean_data <- oid[which(oid$location == "France"),]
clean_data <- clean_data[,c(4,6,9)]
clean_data <- clean_data[complete.cases(clean_data),]


noised_data <- noise(clean_data,3,0.5)$noised_data

imputed_by_mean <- mean_imputation(noised_data,clean_data,3)$imputed_data
imputed_by_learning <- learning_imputation(noised_data,clean_data,1,3)$imputed_data 


n_neighbor <- 10
prediction_horizon <- 7
start_day <- 60
to_remove <- 1
target <- 3
reference <- clean_data


#KNN
knn_result_imputation_by_mean <- evaluation_unwknn(imputed_by_mean[,-to_remove],"new_deaths",reference[,-to_remove],start_day,prediction_horizon,n_neighbor)
mean(knn_result_imputation_by_mean)

knn_result_imputation_by_learning <- evaluation_unwknn(imputed_by_learning[,-to_remove],"new_deaths",reference[,-to_remove],start_day,prediction_horizon,n_neighbor)
mean(knn_result_imputation_by_learning)


#EKNN

eknn_result_imputation_by_mean <- evaluation_eknn(imputed_by_mean[,-to_remove],"new_deaths",reference[,-to_remove],start_day,prediction_horizon,n_neighbor)
mean(eknn_result_imputation_by_mean)

eknn_result_imputation_by_learning <- evaluation_eknn(imputed_by_learning[,-to_remove],"new_deaths",reference[,-to_remove],start_day,prediction_horizon,n_neighbor)
mean(eknn_result_imputation_by_learning)


