#flowdata: data frame of flowcytometry data
df <- flowdata
#dimentions of flowcytometry data
d <- 12
#transformed by pca
res.pca <- prcomp(asinh(df[,1:d]), scale = TRUE)
df2pcs <- data.frame(res.pca$x, df$species)
names(df2pcs)[length(names(df2pcs))] <- 'species'

#flowsom
library(FlowSOM)
#number of principle components
n_pcs <- 3
dftempmatrix <- as.matrix(df2pcs[, 1:n_pcs])
data_FlowSOM <- flowFrame(dftempmatrix)
res.flowsom <- ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
res.flowsom <- BuildSOM(res.flowsom, xdim = 10, ydim = 10, rlen = 1000)
res.flowsom <- BuildMST(res.flowsom, tSNE = T)

#res.flowsom is a flowsom object
#maximum number of clusters
k_max <- 6

#replace NA with node vector for nodes of 0 size
index_na <- unique(which(is.na(res.flowsom$map$medianValues),arr.ind=TRUE)[,1])
res.flowsom$map$medianValues[index_na,] <- res.flowsom$map$codes[index_na,]

#run k means in consensus clustering using the median
res.metacl <- ConsensusClusterPlus::ConsensusClusterPlus(t(res.flowsom$map$medianValues),clusterAlg = 'km', reps = 1000, maxK = k_max, innerLinkage="complete")

#function of predicting number of clusters based on CDF of consensus matrix
#input must be concensuscluster object and k_max
#implementation ref: https://www.biostars.org/p/198789/
k.hat <- function(result, k){
  Kvec <- 2:k
  x1 <- 0.3
  x2 <- 0.7 
  PAC <- rep(NA,length(Kvec)) 
  for(i in Kvec){
    M = result[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
  }
  optK = Kvec[which.min(PAC)]
  return(optK)
}


#performance metrics: accuracy, f1 score, v1 score
#function of accuracy
#defination of accuracy (purity): Manning, C., Raghavan, P., & Schütze, H. (2009)
purityfun <- function(clusters, classes) {
  tb <- table(as.character(clusters), classes)
  total <- sum(tb)
  sum <- 0
  
  for (i in 1: ncol(tb)){
    sum <- sum + max(tb[,i])
    
  }
  return(sum/total)
}

#function of f1 score
#defination of overall f1 score: Fung et al. 2003
f1.measure <- function(classes, clusters){
  tb <- table(as.character(classes), as.character(clusters))
  nk <- ncol(tb)
  nc <- nrow(tb)
  n <- sum(tb)
  sum <- 0
  for (i in 1: nk){
    a <- c()
    for (j in 1:nc){
      if (tb[j,i] != 0 ){
        p <- tb[j,i] / sum(tb[,i])
        r <-  tb[j,i] / sum(tb[j,])
        a <- c(a, 2*p*r/(p+r))
        
      }
      
    }
    sum <- sum + sum(tb[,i])/n * max(a)
    
    
    
  }
  return(sum)
  
}

#function of v1 score
#defination of v-measure: Rosenberg, A., & Hirschberg, J. (2007, June)
#implementation ref: https://gist.github.com/Nowosad/11c686e9b3b9d7b9d838f5de79a601e2
library(infotheo)
v.measure <- function(a, b) {
  mi <- mutinformation(a, b)
  entropy.a <- entropy(a)
  entropy.b <- entropy(b)
  if (entropy.a == 0.0) {
    homogeneity <- 1.0
  } else {
    homogeneity <- mi / entropy.a
  }
  if (entropy.b == 0.0) {
    completeness <- 1.0
  } else {
    completeness <- mi / entropy.b
  }
  if (homogeneity + completeness == 0.0) {
    v.measure.score <- 0.0
  } else {
    v.measure.score <- (2.0 * homogeneity * completeness
                        / (homogeneity + completeness))
  }
  # Can also return homogeneity and completeness if wanted
  return(list(v = v.measure.score, h = homogeneity, c = completeness))
}


#results
#classes:
#classes <- df2pcs$species
k_hat <- k.hat(res.metacl, k_max)
tb <- table(as.character(classes), res.metacl[[k_hat]]$consensusClass[res.flowsom$map$mapping[,1]])
acc <- purityfun(classes, res.metacl[[k_hat]]$consensusClass[res.flowsom$map$mapping[,1]])
f <- f1.measure(classes, res.metacl[[k_hat]]$consensusClass[res.flowsom$map$mapping[,1]])
v <- v.measure(classes, res.metacl[[k_hat]]$consensusClass[res.flowsom$map$mapping[,1]])