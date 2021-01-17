##########################################################
# MAIN FILE TO RUN FOR EVIDENTIAL IMPUTATION EXPERIMENTS #
##########################################################

# Variables and screen cleaning:
graphics.off();cat("\014");rm(list=ls());options(warn=-1); 

library(caret)

# INPUTS
n_run              <- 6
noise_level        <- 0.2
historic_length    <- 7
prediction_horizon <- 14 # in days
n_neighbours       <- 10
forest_size        <- 50
first_eval_date    <- as.Date("2020-02-23")
last_eval_date     <- as.Date("2020-11-02")
inputs <- list(n_run = n_run, noise_level = noise_level, historic_length = historic_length, 
               prediction_horizon = prediction_horizon, n_neighbours = n_neighbours, 
               first_eval_date = first_eval_date, last_eval_date = last_eval_date)

path <- "/Users/nsc/Google Drive/Recherche/Encadrements de these/Rayane Elimam/developpements/Imputation_credibiliste-main"
setwd(path)

source("my_lib_nico.R")

dataset <- read.csv("owid-covid-data.csv")
dataset <- dataset[which(dataset$location == "Algeria"),] # France, Algeria, 
dataset <- dataset[, c("date", "new_cases", "new_deaths")]

dataset$date <- as.Date(dataset$date)
dates        <- unique(dataset$date)
eval_dates   <- dates[dates >= first_eval_date & dates <= last_eval_date]

baseline                            <- c()
knn_rmses_no_imputation             <- c()
knn_rmses_mean_imputation           <- c()
knn_rmses_learning_imputation       <- c()
knn_rmses_temporal_imputation       <- c()
eknn_rmses_mean_imputation          <- c()
eknn_rmses_learning_imputation      <- c()
eknn_rmses_temporal_imputation      <- c()
eknn_rmses_mean_imputation_dist     <- c()
eknn_rmses_learning_imputation_dist <- c()
eknn_rmses_temporal_imputation_dist <- c()
for (i_run in 1 : n_run){
  cat(paste0("\nrun ", i_run, "/", n_run, ": "))
  # Noise injection
  noised_dataset <- noise_injection(dataset, noise_level)
  
  baseline_1run                            <- c()
  knn_rmses_no_imputation_1run             <- c()
  knn_rmses_mean_imputation_1run           <- c()
  knn_rmses_learning_imputation_1run       <- c()
  knn_rmses_temporal_imputation_1run       <- c()
  eknn_rmses_mean_imputation_1run          <- c()
  eknn_rmses_learning_imputation_1run      <- c()
  eknn_rmses_temporal_imputation_1run      <- c()
  eknn_rmses_mean_imputation_dist_1run     <- c()
  eknn_rmses_learning_imputation_dist_1run <- c()
  eknn_rmses_temporal_imputation_dist_1run <- c()
  cat(length(eval_dates), "dates -> ")
  i_date=100 # for debugging
  for (i_date in 1 : length(eval_dates)){
    cat(paste0(i_date, ", "))
    date <- eval_dates[i_date]
    
    # Split train/test
    train_dates <- dataset$date[dataset$date < date]
    test_dates  <- dataset$date[dataset$date >= date & dataset$date < date + prediction_horizon]
    train       <- noised_dataset[noised_dataset$date %in% train_dates, ]
    test        <- dataset[dataset$date %in% test_dates, ]
    
    #Separation attributes/labels
    train_X_raw <- subset(train, select = - new_deaths)
    test_X_raw  <- subset(test, select = - new_deaths)
    train_y <- train$new_deaths
    test_y  <- test$new_deaths
    
    # Past Attributes values computation
    train_X <- compute_past_attributes(train_X_raw, historic_length)
    test_X  <- compute_past_attributes(test_X_raw, historic_length)
    
    # Precise imputation
    train_y_complete_cases      <- train_y[complete.cases(train_y)]
    train_y_mean_imputated      <- impute(train_y, method = "mean")$train_y
    learning_imputation         <- impute(train_y, method = "learning", train_X = train_X)
    train_y_learning_imputated  <- learning_imputation$train_y
    # train_y_temporal_imputated  <- impute(train_y, method = "temporal", train_X = train_X)
    
    # Past labels computation
    train_y_complete_cases_past     <- compute_past_labels(train_y_complete_cases,     historic_length)
    train_y_mean_imputated_past     <- compute_past_labels(train_y_mean_imputated,     historic_length)
    train_y_learning_imputated_past <- compute_past_labels(train_y_learning_imputated, historic_length)
    # train_y_temporal_imputated_past <- compute_past_labels(train_y_temporal_imputated, historic_length)
    test_y_past                     <- compute_past_labels(test_y,                     historic_length)
    
    # Inclusion of past labels in training data
    train_X_complete_cases     <- cbind(train_X[!is.na(train_y), ], train_y_complete_cases_past)
    train_X_mean_imputated     <- cbind(train_X, train_y_mean_imputated_past)
    train_X_learning_imputated <- cbind(train_X, train_y_learning_imputated_past)
    # train_X_temporal_imputated <- cbind(train_X, train_y_temporal_imputated_past)
    test_X                     <- cbind(test_X, test_y_past)
    
    # Evidential modelling
    ### without taking into account distances between neighbours 
    # m_y_mean_imputation     <- mass_construction(train_y_mean_imputated,     method = "mean",     dist = FALSE)
    # m_y_learning_imputation <- mass_construction(train_y_learning_imputated, method = "learning", dist = FALSE, sds = learning_imputation$sds)
    # m_y_temporal_imputation <- mass_construction(train_y_temporal_imputated, method = "temporal", dist = FALSE, chronology = ...)
    ### with taking into account distances between neighbours 
    # m_y_mean_imputation_dist     <- mass_construction(train_y_mean_imputated,     method = "mean",     dist = TRUE)
    # m_y_learning_imputation_dist <- mass_construction(train_y_learning_imputated, method = "learning", dist = TRUE, sds = learning_imputation$sds)
    # m_y_temporal_imputation_dist <- mass_construction(train_y_temporal_imputated, method = "temporal", dist = TRUE, chronology = ...)
    
    # Date removing from attributes:
    train_X_complete_cases     <- subset(train_X_complete_cases, select = - date)
    train_X_mean_imputated     <- subset(train_X_mean_imputated, select = - date)
    train_X_learning_imputated <- subset(train_X_learning_imputated, select = - date)
    # train_X_temporal_imputated <- subset(train_X_temporal_imputated, select = - date)
    test_X                     <- subset(test_X, select = - date)
    
     # Attributes normalisation
    for (j in 1 : ncol(train_X_complete_cases)){
      if (max(train_X_complete_cases[, j]) > 0){
        train_X_complete_cases[, j]     <- train_X_complete_cases[, j]/max(train_X_complete_cases[, j])
        test_X[, j]                     <- test_X[, j]/max(train_X_complete_cases[, j])
      }
      if (max(train_X_mean_imputated[, j]) > 0){
        train_X_mean_imputated[, j]     <- train_X_mean_imputated[, j]/max(train_X_mean_imputated[, j])
      }
      if (max(train_X_learning_imputated[, j]) > 0){
        train_X_learning_imputated[, j] <- train_X_learning_imputated[, j]/max(train_X_learning_imputated[, j])
      }
      # if (max(train_X_temporal_imputated[, j]) > 0){
      #   train_X_temporal_imputated[, j] <- train_X_temporal_imputated[, j]/max(train_X_temporal_imputated[, j])
      # }
    }  
    
    # Prediction
    ### precise models
    preds_baseline          <- rep(mean(train_y_complete_cases), nrow(test_X))
    knn_complete_cases      <- knnreg(x = train_X_complete_cases, y = train_y_complete_cases, k = n_neighbours)
    knn_mean_imputation     <- knnreg(x = train_X_mean_imputated, y = train_y_mean_imputated, k = n_neighbours)
    knn_learning_imputation <- knnreg(x = train_X_learning_imputated, y = train_y_learning_imputated, k = n_neighbours)
    # knn_temporal_imputation <- knn(x = train_X_temporal_imputated,  y = train_y_temporal_imputated, k = n_neighbours)
    preds_complete_cases      <- predict(knn_complete_cases, test_X)
    preds_mean_imputation     <- predict(knn_mean_imputation, test_X)
    preds_learning_imputation <- predict(knn_learning_imputation, test_X)
    # preds_temporal_imputation <- predict(knn_temporal_imputation, test_X)
    ### uncertain models without taking into account distances between neighbours 
    # preds_mean_imputation_uncertain_labels     <- eknn(train_X_mean_imputated,     m_y_mean_imputation,     test_X)
    # preds_learning_imputation_uncertain_labels <- eknn(train_X_learning_imputated, m_y_learning_imputation, test_X)
    # preds_temporal_imputation_uncertain_labels <- eknn(train_X_temporal_imputated, m_y_temporal_imputation, test_X)
    ### uncertain models with taking into account distances between neighbours 
    # preds_mean_imputation_uncertain_labels_dist     <- eknn(train_X_mean_imputated,     m_y_mean_imputation_dist,     test_X)
    # preds_learning_imputation_uncertain_labels_dist <- eknn(train_X_learning_imputated, m_y_learning_imputation_dist, test_X)
    # preds_temporal_imputation_uncertain_labels_dist <- eknn(train_X_temporal_imputated, m_y_temporal_imputation_dist, test_X)
    
    # Evaluation
    ### precise models
    rmse_baseline                <- sum(sqrt(mean((preds_baseline - test_y)^2)))
    rmse_knn_complete_cases      <- sum(sqrt(mean((preds_complete_cases - test_y)^2)))
    rmse_knn_mean_imputation     <- sum(sqrt(mean((preds_mean_imputation - test_y)^2)))
    rmse_knn_learning_imputation <- sum(sqrt(mean((preds_learning_imputation - test_y)^2)))
    # rmse_knn_temporal_imputation <- sum(sqrt(mean((preds_temporal_imputation - test_y)^2)))
    ### uncertain models 
    # rmse_eknn_mean_imputation     <- sum(sqrt(mean((preds_mean_imputation_uncertain_labels - test_y)^2)))
    # rmse_eknn_learning_imputation <- sum(sqrt(mean((preds_learning_imputation_uncertain_labels - test_y)^2)))
    # rmse_eknn_temporal_imputation <- sum(sqrt(mean((preds_temporal_imputation_uncertain_labels - test_y)^2)))
    # rmse_eknn_mean_imputation_dist     <- sum(sqrt(mean((preds_mean_imputation_uncertain_labels_dist - test_y)^2)))
    # rmse_eknn_learning_imputation_dist <- sum(sqrt(mean((preds_learning_imputation_uncertain_labels_dist - test_y)^2)))
    # rmse_eknn_temporal_imputation_dist <- sum(sqrt(mean((preds_temporal_imputation_uncertain_labels_dist - test_y)^2)))
    
    # Incrementation
    baseline_1run                            <- c(baseline_1run,                            rmse_baseline)
    knn_rmses_no_imputation_1run             <- c(knn_rmses_no_imputation_1run,             rmse_knn_complete_cases)
    knn_rmses_mean_imputation_1run           <- c(knn_rmses_mean_imputation_1run,           rmse_knn_mean_imputation)
    knn_rmses_learning_imputation_1run       <- c(knn_rmses_learning_imputation_1run,       rmse_knn_learning_imputation)
    # knn_rmses_temporal_imputation_1run       <- c(knn_rmses_temporal_imputation_1run,       rmse_knn_temporal_imputation)
    # eknn_rmses_mean_imputation_1run          <- c(eknn_rmses_mean_imputation_1run,          rmse_eknn_mean_imputation)
    # eknn_rmses_learning_imputation_1run      <- c(eknn_rmses_learning_imputation_1run,      rmse_eknn_learning_imputation)
    # eknn_rmses_temporal_imputation_1run      <- c(eknn_rmses_temporal_imputation_1run,      rmse_eknn_temporal_imputation)
    # eknn_rmses_mean_imputation_dist_1run     <- c(eknn_rmses_mean_imputation_dist_1run,     rmse_eknn_mean_imputation_dist)
    # eknn_rmses_learning_imputation_dist_1run <- c(eknn_rmses_learning_imputation_dist_1run, rmse_eknn_learning_imputation_dist)
    # eknn_rmses_temporal_imputation_dist_1run <- c(eknn_rmses_temporal_imputation_dist_1run, rmse_eknn_temporal_imputation_dist)
  }
  
  # Incrementation
  baseline                            <- rbind(baseline,                            baseline_1run)
  knn_rmses_no_imputation             <- rbind(knn_rmses_no_imputation,             knn_rmses_no_imputation_1run)
  knn_rmses_mean_imputation           <- rbind(knn_rmses_mean_imputation,           knn_rmses_mean_imputation_1run)
  knn_rmses_learning_imputation       <- rbind(knn_rmses_learning_imputation,       knn_rmses_learning_imputation_1run)
  # knn_rmses_temporal_imputation       <- rbind(knn_rmses_temporal_imputation,       knn_rmses_temporal_imputation_1run)
  # eknn_rmses_mean_imputation          <- rbind(eknn_rmses_mean_imputation,          eknn_rmses_mean_imputation_1run)
  # eknn_rmses_learning_imputation      <- rbind(eknn_rmses_learning_imputation,      eknn_rmses_learning_imputation_1run)
  # eknn_rmses_temporal_imputation      <- rbind(eknn_rmses_temporal_imputation,      eknn_rmses_temporal_imputation_1run)
  # eknn_rmses_mean_imputation_dist     <- rbind(eknn_rmses_mean_imputation_dist,     eknn_rmses_mean_imputation_dist_1run)
  # eknn_rmses_learning_imputation_dist <- rbind(eknn_rmses_learning_imputation_dist, eknn_rmses_learning_imputation_dist_1run)
  # eknn_rmses_temporal_imputation_dist <- rbind(eknn_rmses_temporal_imputation_dist, eknn_rmses_temporal_imputation_dist_1run)
}

# Results aggregation
results <- data.frame(
  date = eval_dates,
  baseline = colMeans(baseline),
  knn_no_imp = colMeans(knn_rmses_no_imputation),
  knn_mean_imp = colMeans(knn_rmses_mean_imputation),
  knn_lear_imp = colMeans(knn_rmses_learning_imputation)
  # knn_temp_imp = colMeans(knn_rmses_temporal_imputation),
  # eknn_mean_imp = colMeans(eknn_rmses_mean_imputation),
  # eknn_lear_imp = colMeans(eknn_rmses_learning_imputation),
  # eknn_temp_imp = colMeans(eknn_rmses_temporal_imputation),
  # eknn_mean_imp_dist = colMeans(eknn_rmses_mean_imputation_dist),
  # eknn_lear_imp_dist = colMeans(eknn_rmses_learning_imputation_dist),
  # eknn_temp_imp_dist = colMeans(eknn_rmses_temporal_imputation_dist)
)

# Results plot
plots <- plot_results(results, dataset, eval_dates)
ggarrange(plots$results, plots$covid, ncol = 1, nrow = 2)

# Results saving
save_results(results, inputs, plots$results)




