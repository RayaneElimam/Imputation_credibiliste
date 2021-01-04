library(FNN)
library(class)
library(caret) 
library(rpart)
library(stringr)
library(imputeTS)
library(evclass)
library(rpart)
library(gridExtra)
library(ggpubr)
library(kknn)
#Bruitage




Bruitage <- function(data,x,prc,end_eval){
  set.seed(123)
  smp_size <- floor(prc * end_eval)
  index <- sample(seq_len(end_eval),size=smp_size)
  data[index,x] <- NA
  return(list(id = index,data = data))
  }


Bruitage_fenetre <- function(data,x,id,max_prc,end_eval) {
  set.seed(123)
  ok <- F
  j <- 1
  while(ok == F){
    x2 <- sample(0:(6-j),length(id),replace = T)
    if((sum(x2) + length(id))/nrow(data) <= max_prc) {
      ok <- T
    }
    j <- j+1
    
  }
    i <- 1
  while (i <= length(id)) {
    b <- id[i]-x2[i]
    if(b > 0 && b < end_eval) {
      data[b:id[i],x] <- NA
      }
    i <- i+1
    }
  return(data)
  }


#Data manipulation

Discretisation <- function(data,ref,x){
  q <- quantile(ref[,x],probs = c(0.33,0.66,0.99))
  data$class <- ifelse(data[,x] <= q[1],"0",ifelse(data[,x] <= q[2] ,"1","2"))
  return(data)
  }

history_cr <- function(data,x,nlags) {
  col_num <- c()
  for (i in 1:nlags) {
    data[,ncol(data)+1] <- NA
    colnames(data)[ncol(data)] <- paste(paste(colnames(data)[x],"lags_"),-i,sep = "")
    col_num <- c(col_num,ncol(data))
    }
  return(list(data = data,cn = col_num))
  }




lags_fill <- function(data,x,col_num,nlags) {
  for (i in (nlags+1):nrow(data)) {
    for (j in 1:nlags) {
      data[i,col_num[j]] <- data[i-j,x]
      }
    }
  return(data[complete.cases(data),])
  }


history <- function(data,dep,arr,nl){
  
  for(i in dep:arr){
    data <- history_cr(data,i,nl)
    cn <- data$cn
    data <- lags_fill(data$data,i,cn,nl)
  }
  return(data)
}


#Imputation methods

imputation <- function(data,x,ref,k) {
  idx <- c(1:nrow(data))
  data <- cbind(data,idx)
  err <- data.frame(err = 0,dist = 0)
  err <- err[-1,]
  data$uncertainty <- 1
  truth <- data[complete.cases(data),]
  for (i in 1:nrow(data)) {
    if(is.na(data[i,x])){
      precedent <- truth[which(truth$idx < data[i,"idx"]),]
      futur <- truth[which(truth$idx > data[i,"idx"]),]
      precedent <- precedent[order(precedent$idx,decreasing = T),]
      futur <- futur[order(futur$idx),]
      prc <- data[i,"idx"] - precedent[1,"idx"]
      fut <- futur[1,"idx"] - data[i,"idx"]
      data[i,x] <- (precedent[1,x] + futur[1,x])/2
      data[i,"uncertainty"] <- k/sum(fut,prc)
      gap <- (ref[i,x] - data[i,x])^2
      dist <- futur[1,"idx"] - precedent[1,"idx"]
      err <- rbind(err,cbind(gap,dist))
      }
    
    }
  return(list(data_full = data[complete.cases(data),-which(names(data) %in% "idx")],error = err$gap))
  }




imputation_pond <- function(data,x,ref,k) {
  idx <- c(1:nrow(data))
  data <- cbind(data,idx)
  err <- data.frame(err = 0,dist = 0)
  err <- err[-1,]
  data$uncertainty <- 1
  truth <- data[complete.cases(data),]
  for (i in 1:nrow(data)) {
    if(is.na(data[i,x])){
      precedent <- truth[which(truth$idx < data[i,"idx"]),]
      futur <- truth[which(truth$idx > data[i,"idx"]),]
      precedent <- precedent[order(precedent$idx,decreasing = T),]
      futur <- futur[order(futur$idx),]
      tot <- futur[1,"idx"]-precedent[1,"idx"]
      prc <- data[i,"idx"] - precedent[1,"idx"]
      fut <- futur[1,"idx"] - data[i,"idx"]
      pp <- 1 - (prc/tot)
      pf <- 1 - (fut/tot)
      data[i,x] <- pp*(precedent[1,x]) +  pf*(futur[1,x])
      data[i,"uncertainty"] <- k/sum(fut,prc)
      gap <- (ref[i,x] - data[i,x] )^2
      dist <- futur[1,"idx"] - precedent[1,"idx"]
      err <- rbind(err,cbind(gap,dist))
      }
    }
  return(list(data_full = data[complete.cases(data),-which(names(data) %in% "idx")],error = err$gap))
}


tree_imputation <- function(data,ref,remove,x){
  
  set.seed(123)
  
  base <- data[complete.cases(data),-1]
  smp_size <- floor(0.8 * nrow(base))
  index <- sample(seq_len(nrow(base)),size=smp_size)
  app <- app[index,]
  test <- base[-index,]
  
  x2 <- grep(colnames(data)[x],colnames(app))
  miss <- data[which(is.na(data[,x])),-remove]
  formula <- paste(colnames(app)[x2],"~")
  attr <- paste(colnames(app)[-x2],collapse = " + ")
  formula <- paste(formula,attr)
  formula <- eval(parse(text = formula))
  model <- rpart(formula, app)
  tstprd <- predict(model,test)
  rmse <- sqrt(mean( (tstprd- test[,x2])^2))
  max_rmse <- sqrt(max( (tstprd- test[,x2])^2)) 
  tree_un <- 1 - (rmse/max_rmse)
  model <- rpart(formula, base)
  prd <- predict(model,miss)
  prd <- as.data.frame(prd)
  data$uncertainty <- 1
  for (i in 1:length(rownames(prd))) {
    data[which(rownames(data) == rownames(prd)[i]),x] <- prd[which(rownames(prd) == rownames(prd)[i]),1]
    
    data[which(rownames(data) %in% rownames(prd)[i]),"uncertainty"] <- tree_un
  }
  
  err <- (prd[,1] - ref[which(rownames(ref) %in% rownames(prd)),x] )^2
  
  return(list(data = data, error = err,f = formula))
  
}



mean_imputation <- function(data,ref,x){
  set.seed(123)
  data$uncertainty <- 1
  nas <- data[which(is.na(data[,x])),]
  imputed <- na_mean(data)
  imputed[which(rownames(imputed) %in% rownames(nas)),"uncertainty"] <- 1.01
  err <- (imputed[which(rownames(imputed) %in% rownames(nas)),x] - ref[which(rownames(ref) %in% rownames(nas)),x])^2  
  return(list(data = imputed,error = err))
}


locf_imputation <- function(data,ref,x){
  set.seed(123)
  data$uncertainty <- 1
  nas <- data[which(is.na(data[,x])),]
  imputed <- na_locf(data)
  err <- (imputed[which(rownames(imputed) %in% rownames(nas)),x] - ref[which(rownames(ref) %in% rownames(nas)),x])^2  
  imputed[which(rownames(imputed) %in% rownames(nas)),"uncertainty"] <- 1.01
  return(list(data = imputed,error = err))
  
  
}


inp_imputation <- function(data,ref,x){
  set.seed(123)
  data$uncertainty <- 1
  nas <- data[which(is.na(data[,x])),]
  imputed <- na_interpolation(data,option = "linear") 
  err <- (imputed[which(rownames(imputed) %in% rownames(nas)),x] - ref[which(rownames(ref) %in% rownames(nas)),x])^2  
  imputed[which(rownames(imputed) %in% rownames(nas)),"uncertainty"] <- 1.01
  return(list(data = imputed,error = err))
  
}


node_stat <- function(model,app ,x){
  
  w <- model$where
  w <- as.data.frame(w)
  w$row <- rownames(w)
  df <- data.frame(leaf = 0,
                   sd_leaf = 0,
                   mean_leaf = 0)
  leaf <- unique(w$w)
  df <- df[0,]
  print("boucle")
  for (i in 1:length(leaf)) {
    
    leaf_instance <- w[which(w$w == leaf[i]),]
    df[i,1] <- leaf[i]
    df[i,2] <- sd(app[which(rownames(app) %in% rownames(leaf_instance)),x])
    df[i,3] <- mean(app[which(rownames(app) %in% rownames(leaf_instance)),x])
    
  }
  return(df)  
}

#Evaluation

evaluation_knn <- function(data,start,end,neighbor){
  print("KNN")
  set.seed(123)
  knn_acc <- c()
  k <- c(start:end)
  for (i in k) {
    t_set <- data[1:i,]
    cl_tset <- t_set$class
    t_set <- t_set[,-which(names(data) %in% c("class","uncertainty") )]
    
    ts_set <- data[(i+1):nrow(data),]
    ts_set <- ts_set[which(ts_set$uncertainty == 1),]
    cl_tsset <- ts_set$class
    ts_set <- ts_set[,-which(names(data) %in% c("class","uncertainty"))]
    

    knn_fit <- knn(t_set, ts_set, cl_tset, k = neighbor)
    acc <- sum(knn_fit == cl_tsset)/nrow(ts_set)
    
    knn_acc <- c(knn_acc,acc)
    nas <- sum(is.na(knn_fit))
    if(nas > 0){
      print("ISSUE KNN")
      data <- NULL
    }
  }
  return(knn_acc)
}

evaluation_kknn <- function(data,start,end,neighbor){
  print("WEIGHTED KNN")
  set.seed(123)
  knn_acc <- c()
  k <- c(start:end)
  data$class <- as.factor(data$class)
  for (i in k) {
    t_set <- data[1:i,]
    t_set <- t_set[,-which(names(data) %in% c("uncertainty") )]
    
    ts_set <- data[(i+1):nrow(data),]
    ts_set <- ts_set[which(ts_set$uncertainty == 1),]
    cl_tsset <- ts_set$class
    ts_set <- ts_set[,-which(names(data) %in% c("uncertainty") )]
    

    model <- kknn(class~.,t_set, ts_set,k = neighbor)
    knn_fit <- model$fitted.values
    acc <- sum(knn_fit == cl_tsset)/nrow(ts_set)
    
    knn_acc <- c(knn_acc,acc)
    nas <- sum(is.na(knn_fit))
    if(nas > 0){
      print("ISSUE KKNN")
      data <- NULL
    }
  }
  return(knn_acc)
}


evaluation_eknn_certain <- function(data,start,end,neighbor,optimize){
  print("EKNN DIS+LAB")
  set.seed(123)
  eknn_certain <- c()
  k <- c(start:end)
  for (i in k) {
    t_set <- data[1:i,]
    ts_set <- data[(i+1):nrow(data),]
    ts_set <- ts_set[which(ts_set$uncertainty == 1),]
    ts_set <- ts_set[,-which(names(ts_set) %in% "uncertainty")]
    t_set[which(t_set$uncertainty == 1.01),"uncertainty"] <- 1
    param <- EkNNfit_uncertain_label(t_set[,-which(names(t_set) %in% c("class","uncertainty"))],
                                         t_set[,"class"],
                                         t_set[,"uncertainty"],
                                       neighbor,optimize = optimize)$param
    
    
    f_eknn <- EkNNval_uncertain_label(t_set[,-which(names(t_set) %in% c("class","uncertainty"))],
                                          t_set[,"class"],
                                          t_set[,"uncertainty"],
                                          ts_set[,-which(names(ts_set) %in% c("class","uncertainty"))],
                                          neighbor,
                                          ts_set[,"class"],param = param)
    f_eknn_perf <- 1 - f_eknn$err
    eknn_certain <- c(eknn_certain,f_eknn_perf)
    nas <- sum(is.na(f_eknn$ypred))
    if(nas > 0){
      print("ISSUE EKNN DIS+LAB")
      data <- NULL
    }
    
  }
  
  return(eknn_certain)
}


evaluation_eknn_dis <- function(data,start,end,neighbor,optimize){
  print("EKNN DIS")
  set.seed(123)
  eknn_certain <- c()
  k <- c(start:end)
  for (i in k) {
    t_set <- data[1:i,]
    ts_set <- data[(i+1):nrow(data),]
    ts_set <- ts_set[which(ts_set$uncertainty == 1),]
    ts_set <- ts_set[,-which(names(ts_set) %in% "uncertainty")]
    t_set[which(t_set$uncertainty == 1.01),"uncertainty"] <- 1
    param <- EkNNfit(t_set[,-which(names(t_set) %in% c("class","uncertainty"))],
                                     t_set[,"class"],
                                     neighbor,optimize = optimize)$param
    
    
    f_eknn <- EkNNval(t_set[,-which(names(t_set) %in% c("class","uncertainty"))],
                                      t_set[,"class"],
                                      ts_set[,-which(names(ts_set) %in% c("class","uncertainty"))],
                                      neighbor,
                                      ts_set[,"class"],param = param)
    f_eknn_perf <- 1 - f_eknn$err
    eknn_certain <- c(eknn_certain,f_eknn_perf)
    nas <- sum(is.na(f_eknn$ypred))
    if(nas > 0){
      print("ISSUE EKNN DIS")
      data <- NULL
    }
    
  }
  
  return(eknn_certain)
}



evaluation_eknn_uncertain <- function(data,start,end,neighbor,optimize){
  print("EKNN LAB")
  set.seed(123)
  eknn_uncertain <- c()
  k <- c(start:end)
  for (i in k) {
    t_set <- data[1:i,]
    ts_set <- data[(i+1):nrow(data),]
    ts_set <- ts_set[which(ts_set$uncertainty == 1),]
    ts_set <- ts_set[,-which(names(ts_set) %in% "uncertainty")]
    t_set[which(t_set$uncertainty == 1.01),"uncertainty"] <- 1
    #param <- EkNNfit_uncertain_label(t_set[,-which(names(t_set) %in% c("class","uncertainty"))],
     #                                     t_set[,"class"],
      #                                    t_set[,"uncertainty"],
       #                                   neighbor,optimize = optimize)$param
    
    
    f_eknn_un <- EkNNval_uncertain_label2(t_set[,-which(names(t_set) %in% c("class","uncertainty"))],
                                          t_set[,"class"],
                                          t_set[,"uncertainty"],
                                          ts_set[,-which(names(ts_set) %in% c("class","uncertainty"))],
                                          neighbor,
                                          ts_set[,"class"])
    f_eknn_un_perf <- 1 - f_eknn_un$err
    eknn_uncertain <- c(eknn_uncertain,f_eknn_un_perf)
    nas <- sum(is.na(f_eknn_un$ypred))
    if(nas > 0){
      print("ISSUE LAB")
      data <- NULL
    }
    
  }
  return(eknn_uncertain)
}


#Pipeline

pipeline_bruitage <- function(data,x,prc,max_prc) {
  set.seed(123)
  data <- Bruitage(data,x,prc)
  data <- Bruitage_fenetre(data$data,x,data$id,max_prc)
  return(list( data = data, na.rate = (sum(is.na(datawna))/nrow(data))))
}

imputation_evaluation <- function(clean_data,prc,max_prc,remove,x,end_eval){
  
  
  set.seed(123)
  ref <-  clean_data
  end_eval <- floor(end_eval * nrow(clean_data))
  data <- Bruitage(clean_data,x,prc,end_eval)
  data <- Bruitage_fenetre(data$data,x,data$id,max_prc,end_eval)
  na.rate <- sum(is.na(data))/nrow(data)
  
  treeimp <- tree_imputation(data,ref,remove,x) 
  evid_imp <- imputation(data,x,ref,1)
  evid_imp2 <- imputation_pond(data,x,ref,1)
  mean_imp <- mean_imputation(data,ref,x)
  locf_imp <- locf_imputation(data,ref,x)
  int_imp <- inp_imputation(data,ref,x)
  
  return(list(tree_imputed = treeimp$data,tree_error = treeimp$error,
              evid_data = evid_imp$data_full,evid_error = evid_imp$error,
              evid_data_pond = evid_imp2$data_full,evid_error2 = evid_imp2$error,
              mean_data = mean_imp$data,mean_error = mean_imp$error,
              locf_data = locf_imp$data,locf_error = locf_imp$error,
              inp_data = int_imp$data,inp_error = int_imp$error,
              na.rate = na.rate))
  
  
}


simple_extract <- function(imputed,i,){
  
   if(i == 1)  data <- imputed$tree_imputed
   if(i == 2)  data <- imputed$evid_data
   if(i == 3)  data <- imputed$evid_data_pond
   if(i == 4)  data <- imputed$mean_data
   if(i == 5)  data <- imputed$locf_data
   if(i == 6)  data <- imputed$inp_data 
   
   return(data)  
  
 }


error_extract <- function(imputed,i){
  
  if(i == 1)  err <- sqrt(mean(imputed$tree_error))
  if(i == 2)  err <- sqrt(mean(imputed$evid_error))
  if(i == 3)  err <- sqrt(mean(imputed$evid_error2))
  if(i == 4)  err <- sqrt(mean(imputed$mean_error))
  if(i == 5)  err <- sqrt(mean(imputed$locf_error))
  if(i == 6)  err <- sqrt(mean(imputed$inp_error))
  
  return(err)  
  
}



pipeline_data <- function(data,ref,x,nl,dep,arr) {
  if(nl > 0){
    data <- history(data,dep,arr,nl)
    data <- Discretisation(data,ref,x)
  }else{
    data <- Discretisation(data,ref,x)
  }
  return(data[complete.cases(data),-x])
}




pipeline_classifierold <- function(imputed,ref,x,sv,ev,nl,dep,arr,neighbor,remove,optimize){
  
  set.seed(123)
  df <- data.frame(matrix(ncol = 3, nrow = 6))
  colnames(df) <- c("knn","eknn certain","eknn incertain")
  rownames(df) <- c("tree_imp","evid_imp","evid_imp_pond","mean_imp",
                    'locf_imp',"inp_imp")
  
  for (i in 1:6) {
    
    data <- simple_extract(imputed,i)
    data <- pipeline_data(data,ref,x,nl,dep,arr)
    start <- floor(sv * nrow(data))
    end <- floor(ev * nrow(data))
    error <- error_extract(imputed,i)
    
    
    knn_test <- evaluation_knn(data[,-remove],start,end,neighbor)
    df[which(rownames(df) == rownames(df)[i]),"knn"] <- mean(knn_test)
    
    eknnc_test <- evaluation_eknn_certain(data[,-remove],start,end,neighbor,optimize)
    df[which(rownames(df) == rownames(df)[i]),"eknn certain"] <- mean(eknnc_test)
    
    eknnu_test <- evaluation_eknn_uncertain(data[,-remove],start,end,neighbor,optimize)
    df[which(rownames(df) == rownames(df)[i]),"eknn incertain"] <- mean(eknnu_test)
    
    df_res <- data.frame(       date = c(start:end),
                                 nn_res  = knn_test,
                                 eknnc_res = eknnc_test,
                                 eknnu_res = eknnu_test)
    
    ggplot(df_res, aes(date)) + 
      geom_line(aes(y = nn_res, colour = paste("KNN mean",round(mean(nn_res),digits = 2) ) )) + 
      geom_line(aes(y = eknnc_res, colour = paste("EKNN CERTAIN mean",round(mean(eknnc_res),digits = 2) ) )) +
      geom_line(aes(y = eknnu_res, colour = paste("EKNN INCERTAIN mean",round(mean(eknnu_res),digits = 2) ) )) +
      theme(legend.position="right") +
      ggtitle(paste(rownames(df)[i],paste("RMSE",round(error,digits = 2)))) +
      xlab("date") + ylab("Accuracy")
    ggsave(paste(rownames(df)[i],paste(imputed$na.rate,".png",sep = ""),sep = "_"),height = 7,width = 10)
    
    
    
  }
  return(df)
  
  
  
}

pipeline_classifier2old <- function(imputed,ref,x,sv,ev,nl,dep,arr,neighbor,remove,optimize){
  
  set.seed(123)
  df <- data.frame(matrix(ncol = 3, nrow = 6))
  colnames(df) <- c("knn","eknn certain","eknn incertain")
  rownames(df) <- c("tree_imp","evid_imp","evid_imp_pond","mean_imp",
                    'locf_imp',"inp_imp")
  
  for (i in 1:6) {
    
    data <- simple_extract(imputed,i)
    data <- pipeline_data(data,ref,x,nl,dep,arr)
    start <- floor(sv * nrow(data))
    end <- floor(ev * nrow(data))
    error <- error_extract(imputed,i)
    data[,-c(remove,(remove+2),ncol(data))] <- scale(data[,-c(remove,(remove+2),ncol(data))],center = T,scale = T)
    print(i)
    print(sum(is.na(data)))
    knn_test <- evaluation_knn(data[,-remove],start,end,neighbor)
    df[which(rownames(df) == rownames(df)[i]),"knn"] <- mean(knn_test)
    
    eknnc_test <- evaluation_eknn_certain(data[,-remove],start,end,neighbor,optimize)
    df[which(rownames(df) == rownames(df)[i]),"eknn certain"] <- mean(eknnc_test)
    
    eknnu_test <- evaluation_eknn_uncertain(data[,-remove],start,end,neighbor,optimize)
    df[which(rownames(df) == rownames(df)[i]),"eknn incertain"] <- mean(eknnu_test)
    
    df_res <- data.frame(       date = c(start:end),
                                nn_res  = knn_test,
                                eknnc_res = eknnc_test,
                                eknnu_res = eknnu_test)
    
    ggplot(df_res, aes(date)) + 
      geom_line(aes(y = nn_res, colour = paste("KNN mean",round(mean(nn_res),digits = 2) ) )) + 
      geom_line(aes(y = eknnc_res, colour = paste("EKNN CERTAIN mean",round(mean(eknnc_res),digits = 2) ) )) +
      geom_line(aes(y = eknnu_res, colour = paste("EKNN INCERTAIN mean",round(mean(eknnu_res),digits = 2) ) )) +
      theme(legend.position="right") +
      ggtitle(paste(rownames(df)[i],paste("RMSE",round(error,digits = 2)))) +
      xlab("date") + ylab("Accuracy")
    ggsave(paste(rownames(df)[i],paste(imputed$na.rate,".png",sep = ""),sep = "_"),height = 7,width = 10)
    
    
    
  }
  return(df)
  
  
  
}




pipeline_classifier2 <- function(imputed,ref,x,sv,ev,nl,dep,arr,neighbor,remove,optimize){
  
  set.seed(123)
  df <- data.frame(matrix(ncol = 5, nrow = 6))
  colnames(df) <- c("knn","kknn","eknn dis+lab","eknn lab","eknn dis")
  rownames(df) <- c("tree_imp","evid_imp","evid_imp_pond","mean_imp",
                    'locf_imp',"inp_imp")
  
  plot_i <- list()
  
  for (i in 1:6) {
    
    data <- simpl+e_extract(imputed,i)
    data <- pipeline_data(data,ref,x,nl,dep,arr)
    start <- floor(sv * nrow(data))
    end <- floor(ev * nrow(data))
    error <- error_extract(imputed,i)
    data[,-c(remove,(remove+2),ncol(data))] <- scale(data[,-c(remove,(remove+2),ncol(data))],center = T,scale = T)
    
    knn_test <- evaluation_knn(data[,-remove],start,end,neighbor)
    df[which(rownames(df) == rownames(df)[i]),"knn"] <- mean(knn_test)
    
    kknn_test <- evaluation_kknn(data[,-remove],start,end,neighbor)
    df[which(rownames(df) == rownames(df)[i]),"kknn"] <- mean(kknn_test)
    
    
    eknnc_test <- evaluation_eknn_certain(data[,-remove],start,end,neighbor,optimize)
    df[which(rownames(df) == rownames(df)[i]),"eknn dis+lab"] <- mean(eknnc_test)
    
    eknnu_test <- evaluation_eknn_uncertain(data[,-remove],start,end,neighbor,optimize)
    df[which(rownames(df) == rownames(df)[i]),"eknn lab"] <- mean(eknnu_test)
    
    eknnd_test <- evaluation_eknn_dis(data[,-remove],start,end,neighbor,optimize)
    df[which(rownames(df) == rownames(df)[i]),"eknn dis"] <- mean(eknnd_test)
    
    
    df_res <- data.frame(       date = c(start:end),
                                nn_res  = knn_test,
                                kknn_res = kknn_test,
                                eknnc_res = eknnc_test,
                                eknnu_res = eknnu_test,
                                eknnd_res = eknnd_test)
    
    
    mean_knn <- round(mean(df_res$nn_res),digits = 2)
    mean_kknn <- round(mean(df_res$kknn_res),digits = 2)
    mean_eknnc <- round(mean(df_res$eknnc_res),digits = 2)  
    mean_eknnu <- round(mean(df_res$eknnu_res),digits = 2)
    mean_eknnd <- round(mean(df_res$eknnd_res),digits = 2)  
    mk <- paste("KNN",mean_knn)
    mean_kknn <- paste("Weighted KNN",mean_kknn)
    mean_eknnc <- paste("EKNN DIS+LAB",mean_eknnc)
    mean_eknnu <- paste("EKNN LAB",mean_eknnu)
    mean_eknnd <- paste("EKNN DIS",mean_eknnd)
    
    k <- paste("K",neighbor)
    nlg <- paste("NLAGS",nl)
    paramk <- paste(k,nlg)
    
    plot_i[[i]] <- ggplot(df_res, aes(date)) + 
      geom_line(aes(y = nn_res, colour = mk )) +
      geom_line(aes(y = kknn_res, colour = mean_kknn )) + 
      geom_line(aes(y = eknnc_res, colour = mean_eknnc )) +
      geom_line(aes(y = eknnu_res, colour = mean_eknnu )) +
      geom_line(aes(y = eknnd_res, colour = mean_eknnd )) +
      ggtitle(paste(rownames(df)[i],paste("NA Rate",paste(round(imputed$na.rate,digits = 2),paste("RMSE",
                                                                                                  paste(round(error,digits = 2),paramk)))))) +
      theme(panel.background = element_rect(fill = 'gray20'),,
            plot.background=element_rect(fill = "gray20"),legend.position="bottom",plot.title = element_text(size = 9, face = "bold"),
            axis.text = element_text(size = 8,face = "bold",colour = "white")
            ,axis.title=element_text(size=8,face="bold",colour = "white"),
            text = element_text(size=4,face = "bold",colour = "white"),legend.direction = "horizontal",
            legend.background = element_rect(fill = "gray20", color = NA),
            legend.key = element_rect(color = "gray20", fill = "gray20"),
            legend.title = element_text(color = "white"),
            legend.text = element_text(color = "white",size = 4),
            legend.key.size = unit(0.3,"line"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "white", fill = NA )) + 
      xlab("date") + ylab("Accuracy") 

              

    ggplot(df_res, aes(date)) + 
      geom_line(aes(y = nn_res, colour = mk )) +
      geom_line(aes(y = kknn_res, colour = mean_kknn )) + 
      geom_line(aes(y = eknnc_res, colour = mean_eknnc )) +
      geom_line(aes(y = eknnu_res, colour = mean_eknnu )) +
      geom_line(aes(y = eknnd_res, colour = mean_eknnd )) +
      ggtitle(paste(rownames(df)[i],paste("NA Rate",paste(round(imputed$na.rate,digits = 2),paste("RMSE",
                                                                                                  paste(round(error,digits = 2),paramk)))))) +
      theme(panel.background = element_rect(fill = 'gray20'),
            plot.background=element_rect(fill = "gray20"),legend.position="bottom",plot.title = element_text(size = 13, face = "bold"),
            axis.text = element_text(size = 16,face = "bold",colour = "white")
            ,axis.title=element_text(size=16,face="bold",colour = "white"),
            text = element_text(size=16,face = "bold",colour = "white"),legend.direction = "horizontal",
            legend.background = element_rect(fill = "gray20", color = NA),
            legend.key = element_rect(color = "gray20", fill = "gray20"),
            legend.title = element_text(color = "white"),
            legend.text = element_text(color = "white",size = 12),
            legend.key.size = unit(0.3,"line"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "white", fill = NA )) + 
      xlab("date") + ylab("Accuracy") 
    
    
    ggsave(paste(rownames(df)[i],paste(imputed$na.rate,".png",sep = ""),sep = "_"),height = 7,width = 10)
    
    
  } 
  ggarrange(plotlist = plot_i,ncol = 2, nrow =3)
  ggsave(paste("full",paste(imputed$na.rate,".png",sep = ""),sep = "_"),height = 7,width = 10)
  
  return(plot_i)
  
}




full_evaluation <- function(clean_data,min_bruitage,max_bruitage,nlags,start_eval,end_eval,neighbors,dep,arr,remove,target,optimize) {
  path <- getwd()
  for (i in 1:length(min_bruitage)) {
    
    Test <- imputation_evaluation(clean_data,min_bruitage[i],max_bruitage[i],
                                  remove,target,end_eval)
    na.rate <- round(Test$na.rate,digits = 2)
    wd <- paste("NA",na.rate)
    dir.create(wd)
    setwd(wd)
    
    eval <- pipeline_classifier2(Test,clean_data,target,start_eval,end_eval,nlags
                                 ,dep,arr,target,remove,optimize)
    setwd(path)
  }
  
}


multiple_eval <- function(clean_data,min_bruitage,max_bruitage,nlags,start_eval,end_eval,neighbors,dep,arr,remove,target,optimize) {
  
  path <- getwd()
  for (neighbor in neighbors ) {
    wd <- paste("K",neighbor)
    dir.create(wd)
    setwd(wd)
    eval <- full_evaluation(clean_data,min_bruitage,max_bruitage,nlags,start_eval,end_eval,neighbor,dep,arr,remove,target,optimize)
    setwd(path)
  }
  
}



final_eval <- function(clean_data,min_bruitage,max_bruitage,nlags,start_eval,end_eval,neighbors,dep,arr,remove,target,optimize){
   path <- getwd()
   wd <- paste("all_eval")
   dir.create(wd)
   setwd(wd)
   transit <- getwd()
   for (nlag in nlags) {
     wd <- paste("nlags",nlag)
     dir.create(wd)
     setwd(wd)
     multiple_eval(clean_data,min_bruitage,max_bruitage,nlag,start_eval,end_eval,neighbors,dep,arr,remove,target,optimize)
     setwd(transit)
   }
  setwd(path)
}


week_transform <- function(clean_data){
  
  clean_data$week <- strftime(clean_data$date, format = "%V")
  reduced <- aggregate(. ~ week, data[,-1], sum)
  
}






