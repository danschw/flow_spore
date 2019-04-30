#~~~~~~~~~~~~#
###Packages###
#~~~~~~~~~~~~#

library(flowCore)
library(tidyverse)
library(viridis) #color palette
library(ggcyto) #visualization
library(cowplot)
library(ggridges)

library(mclust) #for GMM clustering
library(beepr) #sounds

library(testthat) #testing functions
# library(mixtools)
# library(clue)

if(!require(dsHelper)){
  devtools::install_git("https://gitlab.com/drsudo/drsudo_helper.git")  
}
library(dsHelper)

installed.packages()

#~~~~~~~~~~~~~#
###Functions###
#~~~~~~~~~~~~~#

#' firstsecondElement
#'
#' @param x 
#'
#' @description If only one element exists in a vector, the respective element is returned. If it contains more than 2, the second element is returned.
#' @return 
#' @export
#'
#' @examples
firstsecondElement<-function(x){
  ifelse(is.na(x[2]),x[1],x[2])
}

expect_equal(firstsecondElement(c(1)),1)
expect_equal(firstsecondElement(c(1,2)),2)

#' firstsecondElement
#'
#' @param x 
#'
#' @description If only one element exists in a vector, a vector containing the respective element twice is returned. Otherwise, The vector itself is returned
#' @return 
#' @export
#'
#' @examples
firstsecondPair<-function(x){
  if(is.na(x[2])) {c(x[1],x[1])}else {x}
}

testthat::expect_equal(firstsecondPair(c(1,2)),c(1,2))
testthat::expect_equal(firstsecondPair(c(1)),c(1,1))


#' sdnorm
#'
#' @param x numeric
#' @param ... passed into dnorm
#' @description dnorm with additional scaling parameter lambda used as mixing component
#' @return 
#' @export
#'
#' @examples
#' 
sdnorm <- function(x, mean=0, sd=1, lambda=1){lambda*dnorm(x, mean=mean, sd=sd)}

#' ref.ds
#'
#' @description subsets dataframe extracted from fcs file
#' 
#' @param strn strain
#' @param time time
#' @param tripl triplicate number
#' @param stn stain type
#' @param df data frame
#'
#' @return
#' @export
#'
#' @examples

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



#### GMM ####
# 
#   getEmMod<-function(inp.df,inp.type,inp.stain,inp.time,inp.channel,cutoff=5,lowPop=5,highPop=13){
#     
#     inp.df%>%
#       dplyr::filter(type==inp.type,
#                     stain==inp.stain,
#                     time==inp.time)%>%
#       dplyr::filter(channel==inp.channel)%>%
#       dplyr::filter(value>cutoff)%>%
#       with(value)%>%
#       normalmixEM(lambda = .5, mu = c(lowPop, highPop),k = 2, sigma = 0.3)%>%
#       return
#   }
#   
#   
#   getEmCutoff<-function(testmod){
#     
#     #helper function: sdnorm to adjust for mix
#     sdnorm <- function(x, mean=0, sd=1, lambda=1){lambda*dnorm(x, mean=mean, sd=sd)}
#     
#     #function both normals
findCutoffs <-function(x_mu,x_sd,x_pro,y_mu,y_sd,y_pro){
      # findCutoff <- function(x_mu,x_sd,x_pro,y_mu,y_sd,y_pro) {
      #       cat(x_mu,x_sd,x_pro,y_mu,y_sd,y_pro,"\n")
      #       f<-function(x) sdnorm(x, m=x_mu, sd=x_sd, lambda=x_pro) - sdnorm(x, m=y_mu, sd=y_sd, lambda=y_pro)
      #       uniroot(f, interval=c(min(c(x_mu,y_mu)), max(c(x_mu,y_mu))))$root
      # }
      # findCutoff
      
      Vectorize(function(x_mu,x_sd,x_pro,y_mu,y_sd,y_pro) {
            # cat(x_mu,x_sd,x_pro,y_mu,y_sd,y_pro,"\n")
            f<-function(x) sdnorm(x, m=x_mu, sd=x_sd, lambda=x_pro) - sdnorm(x, m=y_mu, sd=y_sd, lambda=y_pro)
            uniroot(f, interval=c(min(c(x_mu,y_mu)), max(c(x_mu,y_mu))))$root
      })(x_mu,x_sd,x_pro,y_mu,y_sd,y_pro)
      # f<-Vectorize(findCutoff)
      # f()
}

findCutoffs(c(8.3,2),c(0.6,0.4),c(0.9,0.6),
           c(11.2,12.2),c(0.6,0.7),c(0.1,0.4))

#   
#   getEmPlot<-function(emcutoff,mix.mod,type="aa",stain=1,channel="FL1-A"){
#     
#     mix.mod$x%>%
#       as.data.frame()%>%
#       ggplot(aes(.))+
#       geom_density()+
#       geom_vline(aes(xintercept=emcutoff),col="red")+
#       scale_y_continuous(expand=c(0,0))+
#       scale_x_continuous(name=paste("stain:",type,
#                                     " concentration:",stain,
#                                     " channel:",channel))
#   }
# 
#   getGroupParam <- function(inp.df=df2,inp.type = "SYBR1", inp.stain = 1,
#                             inp.time = 30, inp.channel = "asinh.FL1.A",
#                             popsep=8) {
#     #subsetting
#     dist1<-inp.df%>%
#       dplyr::filter(type==inp.type,
#                     stain==inp.stain,
#                     time==inp.time)%>%
#       dplyr::filter(channel==inp.channel)%>%
#       with(value)
#     
#     res <- base::split(dist1, dist1 > popsep) %>%
#       sapply(mean)
#     
#     res2 <- base::split(dist1, dist1 > popsep) %>%
#       sapply(length)
#     
#     res3 <- base::split(dist1, dist1 > popsep) %>%
#       sapply(sd)
#     
#     data.frame(
#       sample = paste(inp.type, inp.stain, inp.channel, inp.time, sep = "_"),
#       lowM = res[1],
#       highM = res[2],
#       lowSD = res3[1],
#       highSD = res3[2],
#       lowC = res2[1],
#       highC = res2[2],
#       popsep = popsep
#     )
#   }