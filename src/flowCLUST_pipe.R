#clear environment
rm(list = ls())
graphics.off()
def.par <- par(no.readonly = TRUE)# save default, for resetting...
######################
#load libraries
library(tidyverse)
library(cowplot)
library(flowCore)
library(flowStats)
library(flowClust)
library(ggcyto)
######################

#################################
#         Functions             #
#################################

#function to read and check integer
###################################
read.int <- function (message, low, up){
  #check input is whole numbers
  chk <- sum(is.na(as.integer(c(low,up))))
  if (chk){
    message("bad input for 'low' or 'up'")
    return(NA)
  }
  
  #prompt user for number with input message
  num <- as.integer(readline(paste(message, 
                                   "between ", low, " and ", up, "\n")))
  
  ## check input from user
  #check for number
  if (!is.na(num))
    if (num >= low)
      if(num <= up)
        return (num)
  message ("not a whole number or out of bound")
  return (NA)
}
#####

# function to choose the number of clusters by BIC plot
#######################################################
# using the flowClust package
clust.num <- function (gate, flow.frame, var1, var2, k.low, k.up){
  #try several cluster numbers
  clusts <- flowClust(flow.frame, varNames=c(var1, var2), K=k.low:k.up, B=100)
  par (mfrow=c(2,3))
  for(i in 1:k.up)
    print(plot (clusts[[i]], data=flow.frame, subset=c(var1, var2), level=0.9))
  #plot BIC of bove atempts at clustering
  print(plot (k.low:k.up,criterion(clusts, "BIC"), type="b", main=paste("choose 'k' for",gate)))
  par(def.par)  # reset to default
  #prompt user for number of clusters
  note <- "look at BIC plot and choose the number of clusters to gate (the number where curve starts to saturate): "
  num <- read.int(note, k.low, k.up)
  
  return (num)
  
}
#####


# function to choose cluster(s) of interest
#############################################
choose.clust <- function (clust2choose, flow.frame, cluster, var1, var2, level){
  #store number of clusters
  n<- cluster@K
  
  #plot out clusters and histograms to inform the user when he makes choice
  layout(matrix(c(1,2,4,3), 2, 2, byrow = TRUE), heights = c(1,3))
  par(mar = c(1,3,1,1))
  print(hist(cluster, data = flow.frame, subset = var1, main=var1,col= c(2:(n+1))))
  print(hist(cluster, data = flow.frame, subset = var2, main=var2,col= c(2:(n+1))))
  par(mar = c(3,3,1,1))
  print(plot(cluster, data=flow.frame, subset=c( var1, var2), level=level))
  par (mar = def.par$mar)
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste0("choose the\n", clust2choose, "\ncluster(s)\n(how many and which)"))
  legend("topright", legend = c(1:n), col= c(2:(n+1)), lty=1, lwd=5)
  par(def.par)  # reset to default
  
  # ask user for the number of clusters of interest
  n2choose <- read.int("How many clusters will you choose?\n(choose 0 if no cluster is present)", low = 0, up = n)
  if (!n2choose)
    return(n2choose)
  
  # ask user for the number assignments of clusters of interest
  clust.oi <- rep(NA,n2choose)
  for (i in 1:n2choose){
    while (is.na(clust.oi[i])){
      clust.oi[i] <- read.int(paste(i, ") choose number assignment of cluster: "), low = 1, up = n)
      # check that the numbers enterd are different
      if (anyDuplicated(clust.oi[!is.na(clust.oi)])){
        message("same number chosen twice")
        clust.oi[i] <- NA
      }}}
  return(clust.oi)
}
#####

# save plots
#############################
save.png <- function (what, sample.name){
  f.name <- paste0(what, sample.name,".png")
  dev.copy(device = png,filename = f.name,res = 300,width = 120, height = 100, units = "mm")
  dev.off()
}
# transformation inspired by Karava et al.
Transform.Novocyte <- function(fcsset){
  flowCore::transform(fcsset,
                      `asinh.BL1.A`=asinh(`BL1.A`),
                      `asinh.BL1.H`=asinh(`BL1.H`),
                      `asinh.FSC.A`=asinh(`FSC.A`),
                      `asinh.FSC.H`=asinh(`FSC.H`),
                      `asinh.SSC.A`=asinh(`SSC.A`),
                      `asinh.SSC.H`=asinh(`SSC.H`)
  )
}
#####

#####################
#                   #
#   Main Code       #
#                   #
#####################

# select file and load it
######################
set <- "delta6_DSM_T4_x10"

setwd("C:/Users/danschw/GitHub/flow_spore/data/entrap_data/day1/")# data/entrap_data/day1")
fcs.files <- list.files(pattern = ".fcs", full = TRUE, recursive = T)
fcs.files <- fcs.files[grep(set,x = fcs.files)]
#for mean time I'll work on one file at a time
print(matrix(c(1:length(fcs.files),fcs.files), ncol=2 ))
fcs.files <- fcs.files[read.int (message = "choose file number", low = 1, up = length(fcs.files))]
sample.name <- gsub("./.*/","",fcs.files)
sample.name <- gsub(".fcs","",sample.name)

# save results if user wants to...
sv <- "a"
while (!(sv=="y" | sv=="n")){
  sv <- readline(prompt = "Do you want to save these results? [y/n] ")
  if (!(sv=="y" | sv=="n"))
    message("bad input, try again...")
}

x <- read.FCS(fcs.files, transformation=FALSE,alter.names = T)
# asinh transformation
x <- Transform.Novocyte(x)

n.initial <- as.numeric(x@description$`$TOT`)

#have a look at the file 
# summary(x)

setwd("../output")
#####


# filter out doublets
####################################
# this done by taking the main cluster of FSC height.area ratio (using single cluster: k=1)
singlets.clust <- flowClust(x, varNames=c("FSC.A", "FSC.H"), K=1, B=1000, level=0.995)
 n.singlets <- sum(!is.na(Map(singlets.clust)))

plot(singlets.clust, data=x,  level=0.995,
     main= paste("cluster of singlets ", 100*n.singlets/n.initial, "%"))

if (sv=="y")
  save.png(what = "singlets", sample.name)

#filter outliers
x.singlets <- x[x %in% singlets.clust]%>%
  Subset(rectangleGate("asinh.BL1.A" = c(0, 15), "asinh.SSC.A" = c(0, 15),"asinh.FSC.A" = c(0, 15))) # remove negatives

n.singlets <- nrow(x.singlets@exprs)


#####

# filter out noise
# #####################################
# # for now this will be simply the lower left corner cluster on log-log plot
# 
# 
# ## gate to filter noise
# # use the BIC plot to choose the number of clusters (the number where curve starts to saturate)
# k.noise<- clust.num( "noise",flow.frame = x.singlets, var1 = "asinh.SSC.A", var2 = "asinh.BL1.A", k.low = 1, k.up = 5)
# 
# if (sv=="y")
#   save.png(what = "noiseKchoice", sample.name)
# 
# # cluster using number chosen
# noise.clust <- flowClust(x.singlets, varNames=c("asinh.SSC.A", "asinh.BL1.A"), K=k.noise, B=1000)
# 
# # choose which cluster is noise
# noise.clust.num <- choose.clust("noise",flow.frame = x.singlets, cluster =noise.clust ,var1 ="asinh.SSC.A" ,var2 = "asinh.BL1.A",level = 0.8 )
# noise.clust.num <- sort(noise.clust.num)
# if (sv=="y")
#   save.png(what = "noiseGate", sample.name)
# 
# # filter out the noise cluster
# noise.filter <- tmixFilter("noise.filter",c("FSC.H", "SSC.H"), K=k.noise, B=1000, level=0.8)
# 
# # split data to cells and noise
# cells.clust.num <-c(1:k.noise)
# cells.clust.num <- cells.clust.num[!(cells.clust.num %in% noise.clust.num)]
# x.cells <- split(x.singlets, noise.filter, population=list(cells=cells.clust.num, noise=noise.clust.num))
# 
# n.cells <- nrow(x.cells$cells@exprs)
# 
# #####
# ggcyto(x.cells$cells, aes(asinh.SSC.A, asinh.BL1.A))+
#   geom_hex()

# seperate by SYBR fluoresence
##########################################

# use the BIC plot to choose the number of clusters (the number where curve starts to saturate)
k.gfp<- clust.num("GFP", flow.frame = x.singlets, var1 = "asinh.SSC.A", var2 = "asinh.BL1.A", k.low = 1, k.up = 5)
if (sv=="y")
  save.png(what = "GfpKchoice", sample.name)
# cluster using number chosen
gfp.clust <- flowClust(x.singlets, varNames=c("asinh.SSC.A", "asinh.BL1.A"), K=k.gfp, B=1000)

#choose GFP positive cluster
gfp.pos  <- choose.clust(clust2choose = "GFP positive\n(choose 0 if no GFP+ cluster is present)",flow.frame = x.singlets, cluster =gfp.clust ,var1 ="asinh.SSC.A" ,var2 = "asinh.BL1.A",level = 0.9 )
gfp.pos  <- sort(gfp.pos )
gfp.neg <- c(1:k.gfp)
gfp.neg <- gfp.neg[!(gfp.neg %in% gfp.pos)]
if (sv=="y")
  save.png(what = "GfpGate", sample.name)


#split the populations by GFP gate
gfp.filter <- tmixFilter("gfp.filter",c("asinh.SSC.A", "asinh.BL1.A"), K=k.gfp, B=100, level=0.9)
if (sum(gfp.pos))
  x.gfp <- split(x.singlets, gfp.filter, population=list(pos=gfp.pos, neg=gfp.neg))
if(!sum(gfp.pos)) # this deals with case of no GFP+ cluster
  x.gfp <- split(x.singlets, gfp.filter, population=list(neg=gfp.neg))


#####

# summarizing the results
###################################
#for my convience of plotting I'll extract data into dataframes
df.gfp.neg <- as.data.frame(exprs(x.gfp$neg))
if (sum(gfp.pos))  df.gfp.pos <- as.data.frame(exprs(x.gfp$pos))
if (!sum(gfp.pos)) df.gfp.pos <- df.gfp.neg[0,]
df.all <- as.data.frame(exprs(x))
df.sing <- as.data.frame(exprs(x.singlets))
df.noise <- as.data.frame(exprs(x.cells$noise))

#plot all points together colored by gates
par(def.par)
# par(mfrow=c(1,2))

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), heights = c(9,1))
# par(mar = c(1,3,1,1))
#project gated events on FSC X BL1
plot (log(df.all$asinh.SSC.A), log(df.all$asinh.BL1.A),col="black" , pch=20, cex=0.1,
      xlab="FSC-H (log)", ylab= "GFP-H (log)", main=sample.name)
points (log(df.sing$asinh.SSC.A), log(df.sing$asinh.BL1.A), pch=20 ,col="grey" , cex=0.1)
points (df.noise$asinh.SSC.A, df.noise$asinh.BL1.A, pch=20 ,col="blue" , cex=0.1)
points (df.gfp.pos$asinh.SSC.A, df.gfp.pos$asinh.BL1.A, pch=20 ,col="green" , cex=0.1)
points (df.gfp.neg$asinh.SSC.A, df.gfp.neg$asinh.BL1.A, pch=20 ,col="red", cex=0.1)
# legend("topleft", legend = c("doublets", "noise", "GFP+", "GFP-"),
#        col=c("black", "blue", "green", "red"), pch=20)
#project gated events on FSC X BL1
plot (log(df.all$asinh.SSC.A), log(df.all$SSC.H),col="black" , pch=20, cex=0.1,
      xlab="FSC-H (log)", ylab= "SSC-H (log)", main=sample.name)
points (log(df.sing$asinh.SSC.A), log(df.sing$SSC.H), pch=20 ,col="grey" , cex=0.1)
points (df.noise$asinh.SSC.A, df.noise$SSC.H, pch=20 ,col="blue" , cex=0.1)
points (df.gfp.pos$asinh.SSC.A, df.gfp.pos$SSC.H, pch=20 ,col="green" , cex=0.1)
points (df.gfp.neg$asinh.SSC.A, df.gfp.neg$SSC.H, pch=20 ,col="red", cex=0.1)
# legend("topleft", legend = c("doublets", "noise", "GFP+", "GFP-"),
#        col=c("black", "blue", "green", "red"), pch=20)
par(mar=c(0,0,0,0))
plot.new()
legend("center", legend = c("doublets", "noise", "GFP+", "GFP-"),
       col=c("black", "blue", "green", "red"), pch=20,horiz = T, text.width=c(0.2,0.2))
par(def.par)

if (sv=="y")
  save.png(what = "SumPlots", sample.name)

# write summary of gated events to file
n.events <- data.frame(sample=sample.name,
                       k.noise=k.noise,
                       noise.clust.num=noise.clust.num,
                       k.gfp=k.gfp,
                       gfp.pos.num=gfp.pos,
                       volume=x@description$`$VOL`,
                       initial=n.initial,
                       singlets= n.singlets,
                       noise=nrow(df.noise),
                       gfp.pos= nrow(df.gfp.pos),
                       gfp.neg= nrow(df.gfp.neg))
n.events$undefined <- n.events$singlets-(n.events$noise+n.events$gfp.pos+n.events$gfp.neg)

print(n.events)



if (sv=="y"){
  write.csv(colnames(n.events),file = "colnames.csv", row.names = F, col.names = F)
  write.table(n.events, "gatedByGFP.csv", sep=",",dec = ".",qmethod = "double"  ,append = T, col.names = F)
}  

setwd("../")

