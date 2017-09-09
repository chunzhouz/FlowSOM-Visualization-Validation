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

#reverse pcs to fcs
#function that transforms flowsom model back from pcs to orignal scale (normalized)
sommodelpcstofcs <- function (res.flowsom, res.pca, dffcs, dim){
  data <- asinh(dffcs[,dim])
  dimpcs <- 1: ncol(res.flowsom$data)
  sizemap <- res.flowsom$map$nNodes
  a.flowsom <- res.flowsom
  a.flowsom$data <- data
  a.flowsom$map$codes <- (a.flowsom$map$codes %*% t(res.pca$rotation[,dimpcs]))[,dim]
  a.flowsom$map$colsUsed <- 1:length(dim)
  a.flowsom$prettyColnames <- colnames(data)
  mapping <- as.data.frame(res.flowsom$map$mapping)
  median <- matrix(0, nrow = sizemap, ncol = length(dim))
  sd <- matrix(0, nrow = sizemap, ncol = length(dim))
  colnames(median) <- colnames(data)
  colnames(sd) <- colnames(data)
  medianvalue <- matrix(NA, res.flowsom$map$nNodes, ncol(data))
  for (i in 1: res.flowsom$map$nNodes){
    temp <- matrix(data[a.flowsom$map$mapping[,1] == i,], ncol = ncol(data))
    median <- if (nrow(temp) == 0 | nrow(temp) == 1) a.flowsom$map$codes[i,] else apply(data[a.flowsom$map$mapping[,1] == i,], 2, median)
    medianvalue[i,] <- median
    rm(median)
    rm(temp)
  }

  a.flowsom$map$medianValues <- medianvalue
  colnames(a.flowsom$map$medianValues) <- a.flowsom$prettyColnames
  return(a.flowsom)
}

#for 12 dimensions, choise of variables: c(3,4,5,1,2)
asom1 <- sommodelpcstofcs(res.flowsom, res.pca, df, c(3,4,5,1,2))

#visualization of flowsom
PlotStars(asom1, view = 'grid')
#plot one node (e.g. node 68)
PlotNode(asom1, 68)