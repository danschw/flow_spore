#~~~~~~~~~~~~#
###Packages###
#~~~~~~~~~~~~#

####CRAN####
library(Biobase) #for annotated data frames
library(plyr)
library(tidyverse)
library(forcats)
library(xtable)
library(ggcyto)
library(xtable)
library(forcats)
library(mixtools)
library(plot3D)
library(ggridges)
library(kableExtra)
library(tidyverse)
library(cowplot)
library(ggridges)
library(kableExtra)
library(tidyverse)
library(flowCore)
library(viridis)
library(clue)
library(class)
library(ggforce)
library(scales)


####Bioconductor####
#source("https://bioconductor.org/biocLite.R")
#biocLite()



#~~~~~~~~~~~~~#
###Functions###
#~~~~~~~~~~~~~#



##### Loading FCS Files #####

createFlowset<-function(filepath,sample.variables=c("ident","type","stain","conc"),fact=1){
  #Extracting fcsset rownames from filenames
  fileset<-list.files(path =filepath,pattern = ".fcs") # 
  
  fcsnames<-fileset %>%
    as.data.frame()%>%
    separate(1,into = sample.variables,sep = "_")
  
  fcsnames$fact<-fact #adding additional information
  fcsnames<-fcsnames[,!names(fcsnames)=="trash"] #deleting trash column
  fcsnames$fullname<-apply(fcsnames, 1, paste, collapse="_") #add all names
  
  #Loading FCS Files
  fileloc<-paste0(filepath,fileset)
  fcsset<-read.flowSet(files=fileloc,transformation = "linearize",emptyValue = FALSE,ignore.text.offset = TRUE)
  
  sampleNames(fcsset)<-fcsnames$fullname #changing filename to something unique: fullname
  #creating annotated data frame to add information from filenames into phenoData slot of flowset
  adf<-AnnotatedDataFrame(fcsnames,
                          data.frame(labelDescription=colnames(fcsnames)))
  rownames(adf)<-sampleNames(fcsset) #changing rownames, so adf can be inserted into flowset
  phenoData(fcsset)<-adf
  
  #transformation
  fcsset<-transform(fcsset,
                    `asinh.FL1-A`=asinh(`FL1-A`),
                    `asinh.FL2-A`=asinh(`FL2-A`),
                    `asinh.FL3-A`=asinh(`FL3-A`),
                    `asinh.FL4-A`=asinh(`FL4-A`),
                    `asinh.FL5-A`=asinh(`FL5-A`),
                    `asinh.FSC-A`=asinh(`FSC-A`),
                    `asinh.FSC-H`=asinh(`FSC-H`),
                    `asinh.FSC-W`=asinh(`FSC-W`),
                    `asinh.SSC-A`=asinh(`SSC-A`),
                    `asinh.SSC-H`=asinh(`SSC-H`),
                    `asinh.SSC-W`=asinh(`SSC-W`)
  )
  
  fcsset
}



#### FCS Specific ####

fcstodf<-function(fcsset){
  i<-1
  fcsNames<-pData(phenoData(fcsset))#Sample names and descriptions (if available)
  numbc<-c(1:as.numeric(count(as.data.frame(fcsset@phenoData@data$name))))
  fcs.df<-NULL
  for (i in numbc){
    fcs.temp<-cbind(as.data.frame(exprs(fcsset[[i]])),File=as.factor(fcsNames$name[i]))
    fcs.df<-rbind(fcs.df,fcs.temp)
    print(paste("case number",i,"done"))
  }
  return(fcs.df)
}

samplefcstodf<-function(fcsset,nsamples=10000){
  i<-1
  set.seed(1)
  fcsNames<-pData(phenoData(fcsset))#Sample names and descriptions (if available)
  numbc<-c(1:as.numeric(count(as.data.frame(fcsset@phenoData@data$name)))) #number of samples
  fcs.df<-NULL
  for (i in numbc){
    fcs.temp<-cbind(as.data.frame(exprs(fcsset[[i]])),File=as.factor(fcsNames$name[i]))
    fcs.df<-rbind(fcs.df,fcs.temp[sample(nrow(fcs.temp),nsamples,replace = TRUE),])
    print(paste("case number",i,"done"))
  }
  row.names(fcs.df) <- 1:nrow(fcs.df)
  return(fcs.df)
}

#### Gaussian Mixture Modeling #### 
{
  getEmMod<-function(inp.df,inp.type,inp.stain,inp.time,inp.channel,cutoff=5,lowPop=5,highPop=13){
    
    inp.df%>%
      dplyr::filter(type==inp.type,
                    stain==inp.stain,
                    time==inp.time)%>%
      dplyr::filter(channel==inp.channel)%>%
      dplyr::filter(value>cutoff)%>%
      with(value)%>%
      normalmixEM(lambda = .5, mu = c(lowPop, highPop), sigma = 0.3,k = 2)%>%
      return
  }
  
  #helper function: sdnorm to adjust for the mix:
  #see: https://stackoverflow.com/questions/25313578/any-suggestions-for-how-i-can-plot-mixem-type-data-using-ggplot2
  sdnorm <- function(x, mean=0, sd=1, lambda=1){lambda*dnorm(x, mean=mean, sd=sd)}
  
  getEmCutoff<-function(testmod){
    
    #writing function for both normals
    f <- function(x) {
      sdnorm(x, m=testmod$mu[1], sd=testmod$sigma[1], lambda=testmod$lambda[1]) -
        sdnorm(x, m=testmod$mu[2], sd=testmod$sigma[2], lambda=testmod$lambda[2]) 
    }
    uniroot(f, interval=c(min(testmod$mu), max(testmod$mu)))$root
    
  }
  
  getEmPlot<-function(emcutoff,mix.mod,type="aa",stain=1,channel="FL1-AAA"){
    
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
}


### Group Parameters ####
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
  
  #
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
    #new:
    popsep = popsep
  )
}

#### K-means ####

#Getting reference matrix
ref.ds<-function(strn="Bs02018",tim="48h",stn="1",df=df3A){
  df%>%  
    #dplyr::filter(!stain=="unstained")%>%
    dplyr::filter(strain==strn)%>%
    dplyr::filter(time==tim)%>%
    dplyr::filter(stain==stn)%>%
    dplyr::select(asinh.FSC.A,asinh.SSC.A,asinh.FL1.A)%>%
    as.matrix()
}

#

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
