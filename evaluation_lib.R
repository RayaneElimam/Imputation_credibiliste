library(kknn)
library(class)


evaluation_unwknn <- function(imputed_data,target_name,reference,start_day,prediction_horizon,n_neighbor){
  print("WEIGHTED KNN")
  set.seed(123)
  knn_acc <- c()
  for (i in start_day:(nrow(reference)-prediction_horizon)) {
    train_set <- imputed_data[1:i,]
    train_set <- train_set[,-which(names(train_set) == "uncertainty" )]
    
    test_set <- reference[(i+1):(i+prediction_horizon),]
    
    formula <- paste(target_name,"~.")
    formula <- eval(parse(text = formula))
    
    
    model <- kknn(formula,train_set, test_set,k = n_neighbor,kernel = "rectangular")
    knn_fit <- model$fitted.values
    acc <- sqrt(mean((knn_fit - test_set[,target_name])^2))
    
    knn_acc <- c(knn_acc,acc)
    nas <- sum(is.na(knn_fit))
    if(nas > 0){
      print("ISSUE KKNN")
      data <- NULL
    }
  }
  return(knn_acc)
}



evaluation_wknn <- function(imputed_data,target_name,reference,start_day,prediction_horizon,n_neighbor){
  print("WEIGHTED KNN")
  set.seed(123)
  knn_acc <- c()
  for (i in start_day:(nrow(reference)-prediction_horizon)) {
    train_set <- imputed_data[1:i,]
    train_set <- train_set[,-which(names(train_set) == "uncertainty" )]
    
    test_set <- reference[(i+1):(i+prediction_horizon),]
    
    formula <- paste(target_name,"~.")
    formula <- eval(parse(text = formula))
    
    
    model <- kknn(formula,train_set, test_set,k = n_neighbor,kernel = "triangular")
    knn_fit <- model$fitted.values
    acc <- sqrt(mean((knn_fit - test_set[,target_name])^2))
    
    knn_acc <- c(knn_acc,acc)
    nas <- sum(is.na(knn_fit))
    if(nas > 0){
      print("ISSUE KKNN")
      data <- NULL
    }
  }
  return(knn_acc)
}





evaluation_eknn <- function(imputed_data,target_name,reference,start_day,prediction_horizon,n_neighbor){
  print("EKNN DIS+LAB")
  set.seed(123)
  eknn_certain <- c()
  
  for (i in start_day:(nrow(reference)-prediction_horizon)) {
    train_set <- imputed_data[1:i,]
    test_set <- reference[(i+1):(i+prediction_horizon),]
    
    
    
    f_eknn <- EkNNval_uncertain_regression(train_set[,-which(names(train_set) %in% c(target_name,"uncertainty"))],
                                            train_set[,target_name],
                                            train_set[,"uncertainty"],
                                            test_set[,-which(names(test_set) %in% c(target_name,"uncertainty"))],
                                            n_neighbor,
                                            test_set[,target_name])
    
    eknn_certain <- c(eknn_certain,f_eknn$err)
    nas <- sum(is.na(f_eknn$ypred))
    if(nas > 0){
      print("ISSUE EKNN DIS+LAB")
      data <- NULL
    }
    
  }
  
  return(eknn_certain)
}
