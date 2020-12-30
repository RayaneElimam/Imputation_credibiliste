library(FNN)


EkNNfit_uncertain_label <-function(x,y,y_certainty,K,param=NULL,alpha=0.95,lambda=1/max(as.numeric(y)),optimize=TRUE,
                  options=list(maxiter=300,eta=0.1,gain_min=1e-6,disp=TRUE)){
  #Labels to integer
  y<-as.integer(as.factor(y))
  
  #train samples to matrix
  x<-as.matrix(x)
  
  #If NULL initialize parameters Alpha gamma
  if(is.null(param)) param<-EkNNinit(x,y,alpha)
  
  #Nneighbors for each instance
  knn<-get.knn(x,k=K)
  
  #Distance ^2 give more weights for the nearest neighbors from knn
  knn$nn.dist<-knn$nn.dist^2
  
  #Optimize parameters or not
  if(optimize){
    opt<-optimds_uncertain_label(x,y,y_certainty,param,knn,K,lambda,options)
    
    class <- classds_uncertain_label(opt$param,knn,y,y_certainty,K)
    
    return(list(param=opt$param,cost=opt$cost,err=class$err,ypred=class$ypred,m=class$m))
  }else{
    
    class <- classds_uncertain_label(param,knn,y,y_certainty,K)
    
    return(list(param=param,err=class$err,ypred=class$ypred,m=class$m))
  }
}


EkNNval_uncertain_label <- function(xtrain,ytrain,y_certainty,xtst,K,ytst=NULL,param=NULL){
  
  #Train/test samples to matrix - Labels to integer
  xtst<-as.matrix(xtst)
  
  xtrain<-as.matrix(xtrain)
  
  ytrain<-y<-as.integer(as.factor(ytrain))
  
  if(!is.null(ytst)) ytst<-y<-as.integer(as.factor(ytst))
  
  #If NULL initialize parameters Alpha gamma
  if(is.null(param)) param<-EkNNinit(xtrain,ytrain)
  
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
      m1[ytrain[is[j,i]]] <- (param$alpha * exp(-param$gamma[ytrain[is[j,i]]]^2 * ds[j,i])) * certainty[j,i]
      
      
      m1[M+1] <- 1 - m1[ytrain[is[j,i]]]
      m[1:M,i] <- m1[1:M]*m[1:M,i] + m1[1:M]*m[M+1,i] + m[1:M,i]*m1[M+1]
      m[M+1,i] <- m1[M+1] * m[M+1,i]
      m<-m/matrix(colSums(m),M+1,N,byrow=TRUE)
    }
  }
  
  #Mass matrix
  m<-t(m)
  
  m <- as.data.frame(m)
  #Retire mass on the entire space
  m <- m[,-ncol(m)]
  #Predict maximum mass for each sample
  ypred<-max.col(m[,1:M])
  
  print(max.col(m[,1:M]))
  
  if(!is.null(ytst)) err<-length(which(ypred != ytst))/N else err<-NULL
  
  return(list(m=m,ypred=ypred,err=err))
  
}


EkNNval_uncertain_label2 <- function(xtrain,ytrain,y_certainty,xtst,K,ytst=NULL,param=NULL){
  
  #Train/test samples to matrix - Labels to integer
  xtst<-as.matrix(xtst)
  
  xtrain<-as.matrix(xtrain)
  
  ytrain<-y<-as.integer(as.factor(ytrain))
  
  if(!is.null(ytst)) ytst<-y<-as.integer(as.factor(ytst))
  
  #If NULL initialize parameters Alpha gamma
  if(is.null(param)) param<-EkNNinit(xtrain,ytrain)
  
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
      m1[ytrain[is[j,i]]] <- certainty[j,i] * 0.99 
      
      
      m1[M+1] <- 1 - m1[ytrain[is[j,i]]]
      m[1:M,i] <- m1[1:M]*m[1:M,i] + m1[1:M]*m[M+1,i] + m[1:M,i]*m1[M+1]
      m[M+1,i] <- m1[M+1] * m[M+1,i]
      m<-m/matrix(colSums(m),M+1,N,byrow=TRUE)
    }
  }
  
  #Mass matrix
  m<-t(m)
  
  m <- as.data.frame(m)
  #Retire mass on the entire space
  m <- m[,-ncol(m)]
  #Predict maximum mass for each sample
  ypred<-max.col(m[,1:M])
  
  print(max.col(m[,1:M]))
  
  if(!is.null(ytst)) err<-length(which(ypred != ytst))/N else err<-NULL
  
  return(list(m=m,ypred=ypred,err=err))
  
}




classds_uncertain_label <-function(param,knn,y,y_certainty,K){
  
  #Sample size
  N<-length(y)
  #Class number
  M=max(y)
  
  L<-rep(0,N)
  
  #Mass matrix
  mk <- rbind(matrix(0,M,N),rep(1,N))
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
  
  #Mass attribution for the k neighbors of the N samples - combination with Dempster rule
  for(k in 1:K){
    Is <- is[k,]
    Is = y[Is]
    Tk <- matrix(0,M,N)
    for(j in 1:M){
      pos <- which(Is==j)
      if(length(pos) != 0) Tk[j,pos] <- rep(1,length(pos))
    }
    G <- matrix(param$gamm^2,M,N) * Tk
    
    gam <- apply(G,2,max)
    
    s <- (param$alpha*exp(-gam *ds[k,])) * certainty[k,]
    
    m <- rbind(Tk*matrix(s,M,N,byrow=TRUE),1-s)
    
    mk <- rbind( mk[1:M,]*(m[1:M,]+matrix(m[M+1,],M,N,byrow=TRUE))+
                   m[1:M,]*matrix(mk[M+1,],M,N,byrow=TRUE),mk[M+1,]*m[M+1,])
    Kn <- colSums(mk)
    
    mk <- mk/ matrix(Kn,M+1,N,byrow=TRUE)
  }
  mk<-t(mk)
  
  L<-max.col(mk[,1:M])
  
  err=length(which(L != y))/N
  
  return(list(m=mk,ypred=L,err=err))
}


#Optimization


gradientds_uncertain_label <-function(y,y_certainty,knn,P,alpha,K,lambda){
  
  N<-length(y)
  
  M <- max(y)
  
  Ident <- diag(M)
  
  T <- Ident[,y]
  
  gama <- P
  
  Dgama = rep(0,M)
  
  Ds <- rep(0,N)
  
  Dsgama <- matrix(0,M,N)
  
  mk<-rbind(matrix(0,M,N),rep(1,N))
  
  mm <- mk
  
  s <- matrix(0,K,N)
  
  is<-t(knn$nn.index)
  
  ds<-t(knn$nn.dist)
  
  certainty <- t(knn$nn.dist)
  for (i in 1:nrow(certainty)) {
    for (j in 1:ncol(certainty)) {
      certainty[i,j] <- y_certainty[is[i,j]]
      
    }
    
  }
  
  logK=0
  
  for(k in 1:K){
    Is <- is[k,]
    Is = y[Is]
    Tk <- matrix(0,M,N)
    for(j in 1:M){
      pos <- which(Is==j)
      if(length(pos) != 0) Tk[j,pos] <- rep(1,length(pos))
    }
    G <- matrix(gama^2,M,N) * Tk
    gam <- apply(G,2,max)
    s[k,] <- alpha*exp(-gam *ds[k,]) * certainty[k,]
    m <- rbind(Tk*matrix(s[k,],M,N,byrow=TRUE),1-s[k,])
    mk <- rbind( mk[1:M,]*(m[1:M,]+matrix(m[M+1,],M,N,byrow=TRUE))+
                   m[1:M,]*matrix(mk[M+1,],M,N,byrow=TRUE),mk[M+1,]*m[M+1,])
  }
  
  # Normalization
  Kn <- colSums(mk)
  
  mkn <- mk/ matrix(Kn,M+1,N,byrow=TRUE)
  
  Q <- mkn[1:M,]+lambda*matrix(mkn[M+1,],M,N,byrow=TRUE) - T
  
  ERR <- 0.5*mean(colSums(Q^2))
  
  #		gradient
  
  Dsm <- matrix(0,M+1,N)
  
  for(k in 1:K){
    Is <- is[k,]
    Is <- y[Is]
    Tk <- matrix(0,M,N)
    for(j in 1:M){
      pos <- which(Is==j)
      if(length(pos) != 0) Tk[j,pos] <- rep(1,length(pos))
    }
    m <- rbind(Tk*matrix(s[k,],M,N,byrow=TRUE),1-s[k,])
    if(length(which(m[M+1,] == 0)) > 1){
      Dsgama <- matrix(1e10,M,N)
      k <- K+1
    } else{
      mm[M+1,] <- mk[M+1,]/m[M+1,]
      H = matrix(mm[M+1,],M,N,byrow=TRUE)
      mm[1:M,] <- (mk[1:M,] - H*m[1:M,])/(m[1:M,]+ matrix(m[M+1,],M,N,byrow=TRUE))
      v<-  (mm[1:M,] + matrix(mm[M+1,],M,N,byrow=TRUE)) * Tk - mm[1:M,]
      DsK <- colSums(v) - mm[M+1,]
      Dsm[1:M,] <- (matrix(Kn,M,N,byrow=TRUE)*v -
                      mk[1:M,]*matrix(DsK,M,N,byrow=TRUE) ) / matrix(Kn^2,M,N,byrow=TRUE)
      Dsm[M+1,] <- (-Kn*mm[M+1,] - mk[M+1,]*DsK)/(Kn^2)
      Ds <- colSums(Q*(Dsm[1:M,]+lambda*matrix(Dsm[M+1,],M,N,byrow=TRUE)))
      Dsgama <- Dsgama -2*gama%*%t(Ds*ds[k,]*s[k,]) *Tk
    }
  }
  
  return(list(cost=ERR,grad=rowMeans(Dsgama)))
  
}


optimds_uncertain_label<- function(x,y,y_certainty,param,knn,K,lambda,options){
  
  M <- max(y)
  
  alpha <-param$alpha
  
  P<-param$gamma
  
  a <-1.2
  
  b<-0.8
  
  c <- 0.5
  
  mi <- 1e-4
  
  mx <- 1e6
  
  pas <- options$eta * rep(1,M)
  
  it <- 0
  
  gain <- 1
  
  costgrad<-gradientds_uncertain_label(y,y_certainty,knn,P,alpha,K,lambda)
  
  Errcou1<-costgrad$cost
  
  Dp <- costgrad$grad
  
  Pp <- P
  
  Errcop <- Errcou1 + 1
  
  while((gain >= options$gain_min) & (it <= options$maxiter)){
    it = it + 1;
    costgrad<-gradientds_uncertain_label(y,y_certainty,knn,P,alpha,K,lambda)
    Errco<-costgrad$cost
    D <- costgrad$grad
    if(is.nan(gain) | is.infinite(gain)) gain <- 1
    if(options$disp) print(c(it,Errco,gain))
    if(Errco > Errcop){
      P <- Pp
      D <- Dp
      pas <- pas * c
      P <- P - pas * D
    }else{
      gain <- .9*gain + .1*abs(Errcop - Errco)
      Pp <- P
      test <-  ((D * Dp) >= 0)
      pas <- ((test * a) + ((!test) * b)) * pas
      pas <- (pas <= mx) * pas + (pas > mx) * mx
      pas <- (pas >= mi) * pas + (pas < mi) * mi
      Dp <- D
      P <- P - pas * D
      Errcop <- Errco
    }
  }
  
  return(list(param=list(gamma=P,alpha=alpha),cost=Errco))
}




EkNNinit <- function(x,y,alpha=0.95){
  y<-as.integer(as.factor(y))
  
  x<-as.matrix(x)
  
  M<-max(y)
  
  gamm<-rep(0,M)
  
  for(k in 1:M){
    D<-dist(x[y==k,])
    gamm[k]<-1/sqrt(mean(D))
  }
  return(list(gamma=gamm,alpha=alpha))
}





                                                     train_evid_imputation[,"uncertainty"],
                                                     K = 5,
                                                     optimize = T)





  