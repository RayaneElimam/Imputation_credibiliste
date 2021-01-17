#####################
# FUNCTIONS LIBRARY #
#####################

library(randomForest)
library(ggplot2)
library(reshape2)

#########################################################################################################
noise_injection <- function(dataset, noise_level, max_length_noised_periods = 7){
  noised_dataset <- dataset
  noised_dataset <- noised_dataset[order(noised_dataset$date), ]
  labelled_dates <- noised_dataset$date
  
  nb_label_to_remove <- floor(noise_level * nrow(noised_dataset))
  nb_removed_labels  <- sum(is.na(noised_dataset$date))
  
  # Periods noising:
  while(nb_removed_labels < nb_label_to_remove - max_length_noised_periods){
    index_to_noise          <- sample(1 : length(labelled_dates), 1)
    first_day_noised_period <- labelled_dates[index_to_noise]
    v <- dataset[dataset$date >= first_day_noised_period & 
                   dataset$date <= first_day_noised_period + max_length_noised_periods - 1, "date"]
    max_days_to_noise <- sum(!is.na(v))
    nb_noised_days    <- sample(1 : min(max_length_noised_periods, max_days_to_noise), 1)
    date_to_noise     <- labelled_dates[index_to_noise : eval(index_to_noise + nb_noised_days - 1)]
    noised_dataset[noised_dataset$date %in% date_to_noise, "new_deaths"] <- NA # noise injection
    nb_removed_labels <- sum(is.na(noised_dataset$new_deaths))
  }
  
  # Single dates noising:
  labelled_dates <- noised_dataset[!is.na(noised_dataset$date), 'new_deaths']
  index_to_noise <- sample(1 : length(labelled_dates), nb_label_to_remove - nb_removed_labels)
  noised_dataset$date[index_to_noise] <- NA# noise injection
  
  return(noised_dataset)
}
#########################################################################################################
compute_past_attributes<- function(train_X_raw, historic_length){
  
  # New cases historic
  train_X <- train_X_raw
  for (past_length in 1 : historic_length){
    v <- c()
    for (i_date in 1 : nrow(train_X_raw)){
      if (i_date - past_length < 1){
        past_new_cases <- train_X_raw$new_cases[1] # for first dates, the first data is used
      } else {
        past_new_cases <- train_X_raw$new_cases[i_date - past_length]
      }
      v <- c(v, past_new_cases)
    }
    train_X <- cbind(train_X, v)
  }
  names(train_X)[3 : eval(2 + historic_length)] <- paste0("new_cases_past", 1 : historic_length)
  
  return(train_X)
}
#########################################################################################################
impute <- function(train_y, method, train_X = NULL){
  
  if (method == "mean"){
    train_y[is.na(train_y)] <- mean(train_y, na.rm = T)
    sds <- rep(sd(train_y, na.rm = T), sum(is.na(train_y)))
  }
  if (method == "learning"){
    X_train <- subset(train_X[complete.cases(train_y), ], select = - date)
    y_train <- train_y[complete.cases(train_y)]
    
    X_test <- subset(train_X[!complete.cases(train_y), ], select = - date)
    
    forest           <- randomForest(x = X_train, y = y_train, ntree = 50)
    preds_all_trees  <- predict(forest, X_test, predict.all = T)
    preds <- pmax(preds_all_trees$aggregate, 0)
    sds   <- rowMeans(preds_all_trees$individual^2)-rowMeans(preds_all_trees$individual)^2
    
    train_y[is.na(train_y)] <- preds
  }
  
  
  return(list(train_y = train_y, sds = sds))
}
#########################################################################################################
compute_past_labels <- function(train_y, historic_length){
  
  # New deaths historic
  past_labels <- c()
  for (past_length in 1 : historic_length){
    v <- c()
    for (i_date in 1 : length(train_y)){
      if (i_date - past_length < 1){
        past_label <- train_y[1] # for first dates, the first data is used
      } else {
        past_label <- train_y[i_date - past_length]
      }
      v <- c(v, past_label)
    }
    past_labels <- cbind(past_labels, v)
  }
  past_labels <- as.data.frame(past_labels)
  names(past_labels) <- paste0("new_deaths_past", 1 : historic_length)
  
  return(past_labels)
}
#########################################################################################################
mass_construction <- function(train_y, method, dist, sds = NULL){
  
  # Mass on Omega:
  m_omega <- c()
  for (i in 1 : length(train_y)){
    
  }
  
  return()
}
#########################################################################################################
plot_results <- function(results, dataset, eval_dates){
  
  df_results <- melt(results, id="date")
  names(df_results)[2 : 3] <- c('model', 'RMSE')
  df_covid   <- dataset[dataset$date %in% results$date, ]
  
  p_results <- ggplot(data=df_results, aes(x=date, y=RMSE, colour=model)) + geom_line()+ theme(legend.position="top")
  
  scaleFactor <- max(df_covid$new_cases) / max(df_covid$new_deaths)
  
  p_covid <- ggplot(df_covid, aes(x=date)) +
    geom_line(aes(y=new_cases), col="blue") +
    geom_line(aes(y=new_deaths * scaleFactor), col="red") +
    scale_y_continuous(name="new_cases", sec.axis=sec_axis(~./scaleFactor, name="new_deaths")) +
    theme(
      axis.title.y.left=element_text(color="blue"),
      axis.text.y.left=element_text(color="blue"),
      axis.title.y.right=element_text(color="red"),
      axis.text.y.right=element_text(color="red")
    )
  
  return(list(results = p_results, covid = p_covid))
}
#########################################################################################################
save_results <- function(results, inputs, my_plot){
  
  t <- Sys.time()
  write.csv2(results, paste0("./results/results ", t, ".csv"), row.names = F)
  write.csv2(inputs, paste0("./results/inputs ", t, ".csv"), row.names = F)
  
  ggsave(paste0("./results/plot ", t, '.pdf'), plot = my_plot)
  return()
}
#########################################################################################################















