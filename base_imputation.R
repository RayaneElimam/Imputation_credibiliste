library(rpart)
library(imputeTS)




noise <- function(clean_data,target,noise_level){
  set.seed(123)
  smp_size <- floor(noise_level * nrow(clean_data))
  index_to_noise <- sample(seq_len(nrow(clean_data)),size=smp_size)
  clean_data[index_to_noise,target] <- NA
  return(list(noise_id = index_to_noise,noised_data = clean_data ))
}


frame_noise <- function(noised_data,target,noise_id,max_noise_level) {
  set.seed(123)
  stop <- F
  j <- 1
  while(stop == F){
    frame_noise_size <- sample(0:(7-j),length(noise_id),replace = T)
    if((sum(frame_noise_size) + length(noise_id))/nrow(noised_data) <= max_noise_level) {
      stop <- T
    }
    j <- j+1
    
  }
  i <- 1
  while (i <= length(noise_id)) {
    frame <- noise_id[i]-frame_noise_size[i]
    if(frame > 0) {
      noised_data[frame:noise_id[i],target] <- NA
    }
    i <- i+1
  }
  return(noised_data)
}





mean_imputation <- function(noised_data,reference,target){
  
  noised_data$uncertainty <- 1
  only_noised <- noised_data[which(is.na(noised_data[,target])),]
  imputed <- na_mean(noised_data)
  mean <- mean(noised_data[complete.cases(noised_data),target])
  mean_vec <- rep(mean,nrow(noised_data[complete.cases(noised_data),]))
  err <- (imputed[which(rownames(imputed) %in% rownames(only_noised)),target] - reference[which(rownames(reference) %in% rownames(only_noised)),target])^2
  err_count <- (mean_vec - noised_data[complete.cases(noised_data),target])
  rmse_norm <- sqrt(mean(err_count))/sqrt(max(err_count))
  mean_un <- 1-rmse_norm
  noised_data[which(rownames(noised_data) %in% rownames(only_noised)),target] <- mean_un
  return(list(imputed_data = imputed,error = err))
}



learning_imputation <- function(noised_data,reference,to_remove,target){
  
  set.seed(123)
  
  base <- noised_data[complete.cases(noised_data),-to_remove]
  smp_size <- floor(0.8 * nrow(base))
  train_index <- sample(seq_len(nrow(base)),size=smp_size)
  train <- base[train_index,]
  test <- base[-train_index,]
  
  target2 <- grep(colnames(noised_data)[target],colnames(train))
  only_noised <- noised_data[which(is.na(noised_data[,target])),-to_remove]
  
  formula <- paste(colnames(train)[target2],"~")
  attr <- paste(colnames(train)[-target2],collapse = " + ")
  formula <- paste(formula,attr)
  formula <- eval(parse(text = formula))
  
  model <- rpart(formula, base)
  tstprd <- predict(model,base)
  rmse <- sqrt(mean( (tstprd - base[,target2])^2))
  max_rmse <- sqrt(max( (tstprd- base[,target2])^2)) 
  learning_un <- 1 - (rmse/max_rmse)
  prd <- predict(model,only_noised)
  prd <- as.data.frame(prd)
  noised_data$uncertainty <- 1
  for (i in 1:length(rownames(prd))) {
    noised_data[which(rownames(noised_data) == rownames(prd)[i]),target] <- prd[which(rownames(prd) == rownames(prd)[i]),1]
    
    noised_data[which(rownames(noised_data) %in% rownames(prd)[i]),"uncertainty"] <- learning_un
  }
  err <- (prd[,1] - reference[which(rownames(reference) %in% rownames(prd)),target] )^2
  
  return(list(imputed_data = noised_data, error = err))
  
}
