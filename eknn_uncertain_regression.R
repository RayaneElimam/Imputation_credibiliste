library(FNN)


EkNNval_uncertain_regression <- function(xtrain,ytrain,y_certainty,xtst,K,ytst=NULL){
  
  all <- unique(ytrain)
  all <- all[order(all)]
  #Train/test samples to matrix - Labels to integer
  xtst<-as.matrix(xtst)
  
  xtrain<-as.matrix(xtrain)
  
  ytrain<-y<-as.integer(factor(ytrain,levels = all))
  
  
  #If NULL initialize parameters Alpha gamma
  
  #Napp Train size
  Napp<-nrow(xtrain)
  
  #Class number
  M<-max(ytrain)
  
  #Test sample size
  N<-nrow(xtst)
  
  #nneighbors for test samples from train samples
  knn<-get.knnx(xtrain, xtst, k=K)
  
  #Distance ^2 to give more weight to nearest neighbors
  knn$nn.dist<-knn$nn.dist^2
  
  #Index matrix
  is<-t(knn$nn.index)
  
  #Distance matrix
  ds<-t(knn$nn.dist)
  
  #Certainty matrix
  certainty <- t(knn$nn.dist)
  
  for (i in 1:nrow(certainty)) {
    for (j in 1:ncol(certainty)) {
      certainty[i,j] <- y_certainty[is[i,j]]
      
    }
    
  }
  
  #Mass matrix
  m = rbind(matrix(0,M,N),rep(1,N))
  
  #Mass for each neighbor, combination with dempster rule
  for(i in 1:N){
    for(j in 1:K){
      m1 <- rep(0,M+1)
      m1[ytrain[is[j,i]]] <- (0.95 * exp(-ds[j,i])) * certainty[j,i]
      
      
      m1[M+1] <- 1 - m1[ytrain[is[j,i]]]
      m[1:M,i] <- m1[1:M]*m[1:M,i] + m1[1:M]*m[M+1,i] + m[1:M,i]*m1[M+1]
      m[M+1,i] <- m1[M+1] * m[M+1,i]
      m<-m/matrix(colSums(m),M+1,N,byrow=TRUE)
    }
  }
  
  #Mass matrix
  m<-t(m)
  
  ypred <- rep(0,N)
  for (i in 1:N) {
    sum <- all * m[i,1:M]
    ypred[i] <- sum(sum)
  }
  
  
  if(!is.null(ytst)) err<-sqrt(mean((ypred - ytst)^2)) else err<-NULL
  
  return(list(m=m,ypred=ypred,err=err))
  
}
