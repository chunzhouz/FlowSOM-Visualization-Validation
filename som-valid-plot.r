#validation plot of a flowsom object
#classes: external criteria
#clusters: optional clustering of som nodes
#
#plotflowsom.valid(flowsomobj, classes, clusters)
#
library(ggplot2)
library(scatterpie)
library(RColorBrewer)

plotflowsom.valid <- function (res.flowsom, classes, clusters = NA){
  classes <- factor(paste0(as.character(classes), '(species)'))
  dfplot <- data.frame(res.flowsom$map$grid)
  dfplot$radius <- NA
  for (sp in unique(classes)){
    temp <- rep(NA, nrow(dfplot))
    dfplot <- cbind(dfplot, temp)
    colnames(dfplot)[ncol(dfplot)] <- sp
    rm(temp)
  }
  for (i in 1: res.flowsom$map$nNodes){
    #radius
    dfplot$radius[i] <- sum(res.flowsom$map$mapping[,1] == i)
    dfplot[i ,(ncol(dfplot)-length(unique(classes))+1):ncol(dfplot)] <- table(classes[res.flowsom$map$mapping[,1] == i])
  }
  dfplot$clusters <- factor(paste0('clusters_',as.character(clusters)))
  if (is.na(clusters)){
    g <- ggplot() + 
      geom_scatterpie(data = dfplot, aes(x = Var1, y = Var2), cols = as.character(unique(classes))) + 
      coord_fixed() + 
      scale_fill_brewer(palette="Set2")
  } else {
    g <- ggplot() + 
      geom_rect(data = dfplot, aes(xmin = Var1-0.5, xmax = Var1 + 0.5, ymin = Var2-0.5, ymax = Var2+0.5, fill = as.factor(clusters)), alpha=0.8) + 
      geom_scatterpie(data = dfplot, aes(x = Var1, y = Var2), cols = as.character(unique(classes))) + 
      coord_fixed() + 
      scale_fill_brewer(palette="Set2")
  }
  
      
    
  
  return(g)
}