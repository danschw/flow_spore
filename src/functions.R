#~~~~~~~~~~~~#
###Packages###
#~~~~~~~~~~~~#

library(tidyverse)
library(viridis)
library(flowCore)
library(ggcyto)
library(cowplot)
library(ggridges)
library(mixtools)
library(clue)
library(mclust) #for GMM clustering


if(!require(dsHelper)){
  devtools::install_git("https://gitlab.com/drsudo/drsudo_helper.git")  
}
library(dsHelper)


#~~~~~~~~~~~~~#
###Functions###
#~~~~~~~~~~~~~#

#### GMM ####

  getEmMod<-function(inp.df,inp.type,inp.stain,inp.time,inp.channel,cutoff=5,lowPop=5,highPop=13){
    
    inp.df%>%
      dplyr::filter(type==inp.type,
                    stain==inp.stain,
                    time==inp.time)%>%
      dplyr::filter(channel==inp.channel)%>%
      dplyr::filter(value>cutoff)%>%
      with(value)%>%
      normalmixEM(lambda = .5, mu = c(lowPop, highPop),k = 2, sigma = 0.3)%>%
      return
  }
  
  
  getEmCutoff<-function(testmod){
    
    #helper function: sdnorm to adjust for mix
    sdnorm <- function(x, mean=0, sd=1, lambda=1){lambda*dnorm(x, mean=mean, sd=sd)}
    
    #function both normals
    f <- function(x) {
      sdnorm(x, m=testmod$mu[1], sd=testmod$sigma[1], lambda=testmod$lambda[1]) - #-
        sdnorm(x, m=testmod$mu[2], sd=testmod$sigma[2], lambda=testmod$lambda[2]) 
    }
    uniroot(f, interval=c(min(testmod$mu), max(testmod$mu)))$root
    
  }
  
  getEmPlot<-function(emcutoff,mix.mod,type="aa",stain=1,channel="FL1-A"){
    
    mix.mod$x%>%
      as.data.frame()%>%
      ggplot(aes(.))+
      geom_density()+
      geom_vline(aes(xintercept=emcutoff),col="red")+
      scale_y_continuous(expand=c(0,0))+
      scale_x_continuous(name=paste("stain:",type,
                                    " concentration:",stain,
                                    " channel:",channel))
  }

  getGroupParam <- function(inp.df=df2,inp.type = "SYBR1", inp.stain = 1,
                            inp.time = 30, inp.channel = "asinh.FL1.A",
                            popsep=8) {
    #subsetting
    dist1<-inp.df%>%
      dplyr::filter(type==inp.type,
                    stain==inp.stain,
                    time==inp.time)%>%
      dplyr::filter(channel==inp.channel)%>%
      with(value)
    
    res <- base::split(dist1, dist1 > popsep) %>%
      sapply(mean)
    
    res2 <- base::split(dist1, dist1 > popsep) %>%
      sapply(length)
    
    res3 <- base::split(dist1, dist1 > popsep) %>%
      sapply(sd)
    
    data.frame(
      sample = paste(inp.type, inp.stain, inp.channel, inp.time, sep = "_"),
      lowM = res[1],
      highM = res[2],
      lowSD = res3[1],
      highSD = res3[2],
      lowC = res2[1],
      highC = res2[2],
      popsep = popsep
    )
  }

#### K-means ####
#Getting reference matrix
ref.ds<-function(strn="Bs02003",time="24h",tripl=NA,stn=NA,df=NA){
  if(is.na(stn)){
    df%>%
      dplyr::filter(strain==strn)%>%
      dplyr::filter(tripl==tripl)%>%
      dplyr::filter(time==time)%>%
      dplyr::select(asinh.FSC.A,asinh.SSC.A,asinh.FL1.A)%>%
      as.matrix()
           }else{
     df%>%
             dplyr::filter(strain==strn)%>%
             dplyr::filter(time==time)%>%
             dplyr::filter(stain==stn)%>%
             dplyr::select(asinh.FSC.A,asinh.SSC.A,asinh.FL1.A)%>%
             as.matrix()
         }
}

#

  kmeans.predict<-function(model.kmeans,newdata){
    clusterB.pred<-lapply(newdata,function(x) {
      set.seed(1)
      
      cbind(x,cl_predict(newdata = x[,c(-4,-5,-6)],object = model.kmeans)%>%as.factor())})%>%
      do.call(rbind,.)
    
    colnames(clusterB.pred)[7]<-"cluster"
    
    clusterB.pred
  }
    

  kmeans.predict<-function(model.kmeans,newdata){
      set.seed(1)
      
    clusterB.pred<-cbind(newdata,
                         cluster=cl_predict(newdata = newdata[,c(-4,-5,-6)],object = model.kmeans)
                         %>%as.factor())
      #do.call(rbind,.)
    
    #colnames(clusterB.pred)[7]<-"cluster"
    
    as_tibble(clusterB.pred)
  }
    
#prediction function from package CLUE
predict.clusters2<-function (A, B, method = c("euclidean", "manhattan", "minkowski"),
                             ...) {
  method <- match.arg(method)
  FOO <- switch(method, euclidean = function(A, b) sqrt(rowSums(sweep(A,
                                                                      2, b)^2)), manhattan = function(A, b) rowSums(abs(sweep(A,
                                                                                                                              2, b))), minkowski = {
                                                                                                                                p <- list(...)[[1L]]
                                                                                                                                function(A, b) (rowSums(abs(sweep(A, 2, b))^p))^(1/p)
                                                                                                                              })
  out <- matrix(0, NROW(A), NROW(B))
  for (k in seq_len(NROW(B))) out[, k] <- FOO(A, B[k, ])

  max.col(-out)
}