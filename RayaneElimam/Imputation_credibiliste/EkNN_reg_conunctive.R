library(FNN)
library(class)
library(kknn)

#EkNN regression Conjunctive fusion
EkNNval_uncertain_regression <- function(xtrain,ytrain,y_certainty,xtst,K,ytst=NULL,minset,maxset){
  
  omega <- c(minset:maxset)
  #Train/test samples to matrix - Labels to integer
  xtst<-as.matrix(xtst)
  
  xtrain<-as.matrix(xtrain)
  
  ytrain <- ytrain+1
  
  #If NULL initialize parameters Alpha gamma
  
  #Napp Train size
  Napp<-nrow(xtrain)
  
  #Class number
  M<-maxset+1
  
  #Test sample size
  N<-nrow(xtst)
  
  #nneighbors for test samples from train samples
  knn<-get.knnx(xtrain, xtst, k=K)
  
  #Index matrix
  is<-t(knn$nn.index)
  
  #Distance matrix
  ds<-t(knn$nn.dist)
  
  #Certainty matrix
  labels_mass <- t(knn$nn.dist)
  
  for (i in 1:nrow(labels_mass)) {
    for (j in 1:ncol(labels_mass)) {
      labels_mass[i,j] <- y_certainty[is[i,j]]
      
    }
    
  }
  #Mass matrix
  m = rbind(matrix(0,M,N),rep(1,N))
  
  #Mass for each neighbor, combination with dempster rule
  for(i in 1:N){
    for(j in 1:K){
      m1 <- rep(0,M+1)
      m1[ytrain[is[j,i]]] <- 0.999 * labels_mass[j,i]
      
      
      m1[M+1] <- 1 - m1[ytrain[is[j,i]]]
      m[1:M,i] <- m1[1:M]*m[1:M,i] + m1[1:M]*m[M+1,i] + m[1:M,i]*m1[M+1]
      m[M+1,i] <- m1[M+1] * m[M+1,i]
      m<-m/matrix(colSums(m),M+1,N,byrow=TRUE)
    }
  }
  
  #Mass matrix
  m<-t(m)
  #Pignistic Transformation
  for (i in 1:nrow(m)) {
    
    m[i,1:(ncol(m)-1)] <- m[i,1:(ncol(m)-1)] + ( m[i,ncol(m)] * 1/(ncol(m)-1))
    
  }
  
  ypred <- rep(0,N)
  
  #Pignistic expectation
  for (i in 1:N) {
    sum <- omega * m[i,1:M]
    ypred[i] <- sum(sum)
  }
  
  if(!is.null(ytst)) err<-sqrt(mean((ypred - ytst)^2)) else err<-NULL
  
  return(list(m=m[,-ncol(m)],ypred=ypred,err=err))
  
}



#EkNN regression Conjunctive fusion with missing values rule

EkNNval_uncertain_humble_regression <- function(xtrain,ytrain,y_certainty,xtst,K,ytst=NULL,minset,maxset){
  
  omega <- c(minset:maxset)
  #Train/test samples to matrix - Labels to integer
  xtst<-as.matrix(xtst)
  
  xtrain<-as.matrix(xtrain)
  
  ytrain <- ytrain+1
  
  
  #Napp Train size
  Napp<-nrow(xtrain)
  
  #Class number
  M<-maxset+1
  
  #Test sample size
  N<-nrow(xtst)
  
  #nneighbors for test samples from train samples
  knn<-get.knnx(xtrain, xtst, k=K)
  
  #Index matrix
  is<-t(knn$nn.index)
  
  #Distance matrix
  ds<-t(knn$nn.dist)
  
  #Certainty matrix
  labels_mass <- t(knn$nn.dist)
  for (i in 1:nrow(labels_mass)) {
    for (j in 1:ncol(labels_mass)) {
      labels_mass[i,j] <- y_certainty[is[i,j]]
      
    }
    
  }
  
  
  #Mass matrix
  m = rbind(matrix(0,M,N),rep(1,N))
  
  #Mass for each neighbor, combination with dempster rule
  for(i in 1:N){
    for(j in 1:K){
      if(!is.na(ytrain[is[j,i]]) ){
        m1 <- rep(0,M+1)
        m1[ytrain[is[j,i]]] <- 0.999 * labels_mass[j,i] 
        m1[M+1] <- 1 - m1[ytrain[is[j,i]]]
        
        m[1:M,i] <- m1[1:M]*m[1:M,i] + m1[1:M]*m[M+1,i] + m[1:M,i]*m1[M+1]
        m[M+1,i] <- m1[M+1] * m[M+1,i]
        m<-m/matrix(colSums(m),M+1,N,byrow=TRUE)
        
      } else{
        m1 <- rep(0,M+1)
        m1[M+1] <- 1  
        
        m[1:M,i] <- m1[1:M]*m[1:M,i] + m1[1:M]*m[M+1,i] + m[1:M,i]*m1[M+1]
        m[M+1,i] <- m1[M+1] * m[M+1,i]
        m<-m/matrix(colSums(m),M+1,N,byrow=TRUE)
        
      }
    }
  }
  #Mass matrix
  m<-t(m)
  
  #Pignistic transformation
  for (i in 1:nrow(m)) {
    
    m[i,1:(ncol(m)-1)] <- m[i,1:(ncol(m)-1)] + ( m[i,ncol(m)] * 1/(ncol(m)-1))
    
  }
  
  ypred <- rep(0,N)
  
  #Pignistic expectation
  for (i in 1:N) {
    sum <- omega * m[i,1:M]
    ypred[i] <- sum(sum)
  }
  
  
  if(!is.null(ytst)) err<-sqrt(mean((ypred - ytst)^2)) else err<-NULL
  
  return(list(m=m[,-ncol(m)],ypred=ypred,err=err))
  
}


#EkNN time regression Conjunctive fusion
EkNNval_uncertain_time_regression <- function(xtrain,ytrain,y_certainty,xtst,K,date,ytst=NULL,minset,maxset){
  
  omega <- c(minset:maxset)
  #Train/test samples to matrix - Labels to integer
  xtst<-as.matrix(xtst)
  
  xtrain<-as.matrix(xtrain)
  
  ytrain <- ytrain+1
  
  #If NULL initialize parameters Alpha gamma
  
  #Napp Train size
  Napp<-nrow(xtrain)
  
  #Class number
  M<-maxset+1
  
  #Test sample size
  N<-nrow(xtst)
  
  #nneighbors for test samples from train samples
  knn<-get.knnx(xtrain[,date], xtst[,date], k=K)
  
  #Index matrix
  is<-t(knn$nn.index)
  
  #Distance matrix
  ds<-t(knn$nn.dist)
  
  #Certainty matrix
  labels_mass <- t(knn$nn.dist)
  
  for (i in 1:nrow(labels_mass)) {
    for (j in 1:ncol(labels_mass)) {
      labels_mass[i,j] <- y_certainty[is[i,j]]
      
    }
    
  }
  #Mass matrix
  m = rbind(matrix(0,M,N),rep(1,N))
  
  #Mass for each neighbor, combination with dempster rule
  for(i in 1:N){
    for(j in 1:K){
      m1 <- rep(0,M+1)
      m1[ytrain[is[j,i]]] <- 0.999 * labels_mass[j,i]
      
      
      m1[M+1] <- 1 - m1[ytrain[is[j,i]]]
      m[1:M,i] <- m1[1:M]*m[1:M,i] + m1[1:M]*m[M+1,i] + m[1:M,i]*m1[M+1]
      m[M+1,i] <- m1[M+1] * m[M+1,i]
      m<-m/matrix(colSums(m),M+1,N,byrow=TRUE)
    }
  }
  
  #Mass matrix
  m<-t(m)
  #Pignistic Transformation
  for (i in 1:nrow(m)) {
    
    m[i,1:(ncol(m)-1)] <- m[i,1:(ncol(m)-1)] + ( m[i,ncol(m)] * 1/(ncol(m)-1))
    
  }
  
  ypred <- rep(0,N)
  
  #Pignistic expectation
  for (i in 1:N) {
    sum <- omega * m[i,1:M]
    ypred[i] <- sum(sum)
  }
  
  if(!is.null(ytst)) err<-sqrt(mean((ypred - ytst)^2)) else err<-NULL
  
  return(list(m=m[,-ncol(m)],ypred=ypred,err=err))
  
}



#EkNN time regression Conjunctive fusion with missing values rule

EkNNval_uncertain_humble_time_regression <- function(xtrain,ytrain,y_certainty,xtst,K,date,ytst=NULL,minset,maxset){
  
  omega <- c(minset:maxset)
  #Train/test samples to matrix - Labels to integer
  xtst<-as.matrix(xtst)
  
  xtrain<-as.matrix(xtrain)
  
  ytrain <- ytrain+1
  
  
  #Napp Train size
  Napp<-nrow(xtrain)
  
  #Class number
  M<-maxset+1
  
  #Test sample size
  N<-nrow(xtst)
  
  #nneighbors for test samples from train samples
  knn<-get.knnx(xtrain[,date], xtst[,date], k=K)
  
  #Index matrix
  is<-t(knn$nn.index)
  
  #Distance matrix
  ds<-t(knn$nn.dist)
  
  #Certainty matrix
  labels_mass <- t(knn$nn.dist)
  for (i in 1:nrow(labels_mass)) {
    for (j in 1:ncol(labels_mass)) {
      labels_mass[i,j] <- y_certainty[is[i,j]]
      
    }
    
  }
  
  
  #Mass matrix
  m = rbind(matrix(0,M,N),rep(1,N))
  
  #Mass for each neighbor, combination with dempster rule
  for(i in 1:N){
    for(j in 1:K){
      if(!is.na(ytrain[is[j,i]]) ){
        m1 <- rep(0,M+1)
        m1[ytrain[is[j,i]]] <- 0.999 * labels_mass[j,i] 
        m1[M+1] <- 1 - m1[ytrain[is[j,i]]]
        
        m[1:M,i] <- m1[1:M]*m[1:M,i] + m1[1:M]*m[M+1,i] + m[1:M,i]*m1[M+1]
        m[M+1,i] <- m1[M+1] * m[M+1,i]
        m<-m/matrix(colSums(m),M+1,N,byrow=TRUE)
        
      } else{
        m1 <- rep(0,M+1)
        m1[M+1] <- 1  
        
        m[1:M,i] <- m1[1:M]*m[1:M,i] + m1[1:M]*m[M+1,i] + m[1:M,i]*m1[M+1]
        m[M+1,i] <- m1[M+1] * m[M+1,i]
        m<-m/matrix(colSums(m),M+1,N,byrow=TRUE)
        
      }
    }
  }
  #Mass matrix
  m<-t(m)
  
  #Pignistic transformation
  for (i in 1:nrow(m)) {
    
    m[i,1:(ncol(m)-1)] <- m[i,1:(ncol(m)-1)] + ( m[i,ncol(m)] * 1/(ncol(m)-1))
    
  }
  
  ypred <- rep(0,N)
  
  #Pignistic expectation
  for (i in 1:N) {
    sum <- omega * m[i,1:M]
    ypred[i] <- sum(sum)
  }
  
  
  if(!is.null(ytst)) err<-sqrt(mean((ypred - ytst)^2)) else err<-NULL
  
  return(list(m=m[,-ncol(m)],ypred=ypred,err=err))
  
}




