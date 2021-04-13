library(ibelief)
library(FNN)
library(class)

EkNNval_uncertain_regression_disjunctive <- function(xtrain,ytrain,y_certainty,xtst,K,date,ytst=NULL,minset,maxset){
  #Omega space
  omega <- c(minset:maxset)
  
  #Train/test samples to matrix - Labels to integer
  xtst<-as.matrix(xtst)
  
  xtrain<-as.matrix(xtrain)
  
  #Test sample size
  N<-nrow(xtst)
  
  #nneighbors for test samples from train samples
  knn<-get.knnx(xtrain[,date], xtst[,date], k=K)
  
  #Index matrix
  is<-t(knn$nn.index)
  
  #Distance matrix
  ds <- t(knn$nn.dist)
  
  #Certainty matrix
  label_mass <- t(knn$nn.dist)
  
  for (i in 1:nrow(label_mass)) {
    for (j in 1:ncol(label_mass)) {
      label_mass[i,j] <- y_certainty[is[i,j]]
      
    }
    
  }
  #Mass matrix
  ps <- powerset(omega)
  m <- rep(0,length(ps))
  
  #Mass for each neighbor, combination with dempster rule
  for(i in 1:N){
    neighbormass <- matrix(nrow = length(ps),ncol = K)
    print(dim(neighbormass))
    neighbormass[,] <- 0 
    neighbormass[nrow(neighbormass),1:ncol(neighbormass)] <- 1
    for(j in 1:K){
      idx <- getVecIndex(ytrain[is[j,i]],ps)
      print(idx)
      neighbormass[idx,j] <- 0.999 * label_to_mass[j,i]
      neighbormass[nrow(neighbormass),j] <- neighbormass[nrow(neighbormass),j] - neighbormass[idx,j]
    }
    M_comb=DST(neighbormass,4)
    m <- cbind(m,M_comb)
  }
  
  #Mass matrix
  m <- m[,-1]
  m<-t(m)
  
  #Pignistic transformation
  
  for (i in 1:length(ps)) {
    ps[[i]] <- as.character(ps[[i]])
    
  }
  
  
  names <- as.vector(ps)
  colnames(m) <- names
  colnames(m) <- str_remove_all(colnames(m),"\"")
  colnames(m) <- str_remove_all(colnames(m)," ")
  
  pignistic <- matrix(nrow = nrow(m),ncol = length(omega),0)
  element <- as.character(omega)
  colnames(pignistic) <- element
  
  
  
  for (i in 1:nrow(pignistic)) {
    for (j in 1:ncol(pignistic)) {
      for (f in 1:ncol(m)) {
        if(grepl(colnames(pignistic)[j],colnames(m)[f])){
          cardinal <- str_count(colnames(m)[f],",") + 1
          pignistic[i,j] <- pignistic[i,j] + (m[i,f]/cardinal)
          print(i)
        }
      }
    }
  }
  
  #Pignistic expectation
  
  ypred <- rep(0,N)
  
  for (i in 1:N) {
    sum <- omega * pignistic[i,]
    ypred[i] <- sum(sum)
  }
  
  ypred <- as.integer(ypred)
  if(!is.null(ytst)) err<-sqrt(mean((ypred - ytst)^2)) else err<-NULL
  
  return(list(m=m[,-ncol(m)],ypred=ypred,err=err))
  
}

#####

disjunctive_fusion <- function(labels,labels_mass){
  #Mass_1 u Mass_2
  first_union <- union(labels[1],labels[2])
  mass_fusion <- labels_mass[1]*labels_mass[2]
  omega_mass <- 1 - mass_fusion 
  union <- first_union
  #Mass_1,2 U Mass_3 U....Mass_I
  for (i in 3:length(labels)) {
    union <- union(union,labels[i])
    mass_fusion <- mass_fusion * labels_mass[i]
    omega_mass <- 1 - mass_fusion
  }
  return(list(focal_elements = union,mass = mass_fusion,omega_mass = omega_mass ))
}

y_exemple <- c(1,2,5)
mass_exemple <- c(0.3,0.4,0.6)

result <- disjunctive_fusion(y_exemple,mass_exemple)
#Focal elements
result$focal_elements
#Mass for focal elements
result$mass
#Mass for omega
result$omega_mass



EkNNval_uncertain_regression_disjunctiveV2 <- function(xtrain,ytrain,y_certainty,xtst,K,date,ytst=NULL,minset,maxset){
  
  omega <- c(minset:maxset)
  #Train/test samples to matrix - Labels to integer
  xtst <- as.matrix(xtst)
  
  xtrain <- as.matrix(xtrain)
  
  #Test sample size
  N <- nrow(xtst)
  
  #nneighbors for test samples from train samples
  knn <- get.knnx(xtrain, xtst, k=K)
  
  #Index matrix
  is<-t(knn$nn.index)
  
  #Distance matrix
  ds<-t(knn$nn.dist)
  
  #Certainty matrix
  label_mass <- t(knn$nn.dist)
  ypred <- rep(0,N)
  
  for (i in 1:nrow(label_mass)) {
    for (j in 1:ncol(label_mass)) {
      label_mass[i,j] <- y_certainty[is[i,j]]
      
    }
    
  }
  #Mass matrix
  #Mass for each neighbor, combination with dempster rule
  
  
  for(i in 1:N){
    neighbor <- as.vector(ytrain[is[i,]])
    neighbors_mass <- label_mass[i,] * 0.99
    fusion <- disjunctive_fusion(neighbor,neighbors_mass)
    pignistic_vector <- rep(0,length(omega))
    pignistic_vector <- t(pignistic_vector)
    pignistic_vector <- as.matrix(pignistic_vector)
    for (j in 1:ncol(pignistic_vector)) {
      if(j %in% fusion$focal_elements){
        pignistic_vector[1,j] <- fusion$mass/length(fusion$focal_elements)
      }
      pignistic_vector[1,j] <- pignistic_vector[1,j] + (fusion$omega_mass/length(omega))
    }
    sum <- omega * pignistic_vector
    ypred[i] <- sum(sum)
  }
  
  #Mass matrix
  
  
  
  if(!is.null(ytst)) err<-sqrt(mean((ypred - ytst)^2)) else err<-NULL
  print(ytst)
  return(list(ypred=ypred,err=err))
  
}



