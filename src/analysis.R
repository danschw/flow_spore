
#packages and functions
source("src/functions.R")

#### Figure 1 - Separation of Cells vs Spores possible, overview plot ####
#loading
sample.var<-c("strain","ident","stain","mode","run")
fcsset1<-createFlowset(filepath = "data/f1/",sample.variables=sample.var)
table(df1$ident,df1$mode)

#gating
fcsset1.g<-fcsset1%>%
  Subset(.,norm2Filter("asinh.FSC.H", "asinh.FSC.W",filterId="norm_ssc.fsc",scale=1))%>%
  Subset(.,norm2Filter("asinh.SSC.H", "asinh.SSC.W",filterId="norm_ssc.fsc",scale=1))

df1<-fcstodf(fcsset1.g)%>%
  separate(length(fcsset1@colnames)+1,into = sample.var,sep = "_")%>%
  mutate(stain=as.factor(strain),ident=as.factor(ident),stain=as.factor(stain),mode=as.factor(mode))

{
#SSC-Plot
SSC.plot<-df1%>%
  dplyr::filter(mode=="SSC")%>%
  select(asinh.SSC.A,ident,stain)%>% #,run
  #dplyr::filter(stain=="unstained",ident %in% c("spores","cells"))%>%
  ggplot(aes(asinh.SSC.A))+
  geom_density(aes(fill=ident,y=..scaled..),alpha=0.4)+
  ylab("norm. count")+
  guides(fill=FALSE)+theme(legend.position = c(.05,.9))+
  scale_fill_viridis(discrete = TRUE,alpha=0.7,begin = 0,end = 0.8,name="",direction = -1)

SSC.plot

#OR CHANGE TO FSC!

#FSC-Plot
FSC.plot<-df1%>%
  dplyr::filter(mode=="PI")%>%
  select(asinh.FSC.A,ident,stain,run)%>%
  ggplot(aes(asinh.FSC.A))+
  geom_density(aes(fill=ident,y=..scaled..),alpha=0.4)+
  ylab("")+
  guides(fill=FALSE)+
  scale_fill_viridis(discrete = TRUE,alpha=0.7,begin = 0,end = 0.8,direction = -1)
FSC.plot

#pi
PI.plot<-df1%>%
  dplyr::filter(mode=="PI")%>%
  select(asinh.FL3.A,ident,stain,run)%>%
  #dplyr::filter(stain=="PI",run=="2x")%>%
  ggplot(aes(asinh.FL3.A))+
  geom_density(aes(fill=ident,y=..scaled..),alpha=0.4)+
  guides(fill=FALSE)+
  scale_fill_viridis(discrete = TRUE,alpha=0.7,begin = 0,end = 0.8,direction = -1)+
  xlim(c(5,12))+
  ylab("")

PI.plot

#sybr1
S1.plot<-df1%>%
  dplyr::filter(mode=="SYBR1")%>%
  select(asinh.FL1.A,ident,stain,run)%>%
  ggplot(aes(asinh.FL1.A))+
  geom_density(aes(fill=ident,y=..scaled..),alpha=0.4)+
  guides(fill=FALSE)+
  ylab("norm. count")+xlim(c(5,14))+
  scale_fill_viridis(discrete = TRUE,alpha=0.7,begin = 0,end = 0.8,direction = -1)
S1.plot

#sybr2
S2.plot<-df1%>%
  dplyr::filter(mode=="SYBR2")%>%
  select(asinh.FL1.A,ident,stain,run)%>%
  ggplot(aes(asinh.FL1.A))+
  geom_density(aes(fill=ident,y=..scaled..),alpha=0.4)+
  guides(fill=FALSE)+
  ylab("norm. count")+xlim(c(5,14))+
  scale_fill_viridis(discrete = TRUE,alpha=0.7,begin = 0,end = 0.8,direction = -1)
S2.plot

# Legend Plot
leg.plot<-df1%>%
  select(asinh.FL1.A,ident,stain,run)%>%
  ggplot(aes(asinh.FL1.A))+geom_density(aes(fill=ident),alpha=0.4)+
  ylab("")+xlab("")+
  scale_fill_viridis(discrete = TRUE,alpha=0.7,begin = 0,end = 0.8,
                     labels=c("sample: cells","sample: spores"),
                     name="",direction = -1)+
  theme(legend.position = c(0.01,0.95))+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  ylim(0,100)

leg.plot
}

plot_grid(SSC.plot,FSC.plot,PI.plot,S1.plot,S2.plot,leg.plot,ncol = 3,labels = "AUTO")

ggsave("fig/f1_roc.pdf",width = 9, height = 5,units = "in",dpi=300)

#### Figure 2: - Measuring both together: ability to separate ####
{#loading
  rm(list=ls())
  source("src/functions.R")
  
  sample.var=c("ident","type","stain","time") 
  fcsset2.1<-createFlowset(filepath = "data/f2/set5_cells+spores_30/",
                         sample.variables=sample.var,fact=30)
  fcsset2.2<-createFlowset(filepath = "data/f2/set5_cells+spores_90/",
                           sample.variables=sample.var,fact=90)
  fcsset2.3<-createFlowset(filepath = "data/f2/set6_cells+spores_30/",
                           sample.variables=sample.var,fact=30)
  fcsset2.4<-createFlowset(filepath = "data/f2/set6_cells+spores_90/",
                           sample.variables=sample.var,fact=90)
  fcsset2.5<-rbind2(rbind2(fcsset2.1,fcsset2.2),rbind2(fcsset2.3,fcsset2.4))
}

"
  fcsset2.5<-rbind2(fcsset2.1,fcsset2.2)
  fcsset2.6<-rbind2(fcsset2.3,fcsset2.4)
  fcsset2.5
  fcsset2.6
  
  
  to do:
  * Do proper gating of singlets (CHECK!)
  * Check why stupid thresholding doesnt work
  * 
  "

fcsset2.5%>%
  autoplot("asinh.FSC.H","asinh.FSC.W",bins=200)+
  theme_minimal()

fcsset2.5%>%
  Subset(.,norm2Filter("FSC-H", "FSC-W",filterId="norm_ssc.fsc",scale.factor = 4))%>%
  autoplot("asinh.FSC.H","asinh.FSC.W",bins=200)+
  theme_minimal()

fcsset2.5g<-fcsset2.5%>%
  Subset(.,norm2Filter("FSC-H", "FSC-W",filterId="norm_ssc.fsc",scale=4))


#
df2<-fcstodf(fcsset2.5g)%>%
  separate(length(fcsset2.5g@colnames)+1,into = c(sample.var,"time2"),sep = "_")%>%
  dplyr::filter(ident=="Cells+spores")%>%
  select(type,stain,time2,asinh.FL1.A,asinh.FL3.A,asinh.FSC.A,asinh.SSC.A)%>%
  gather("channel","value",4:7)

colnames(df2)[3]<-"time"

str(df2)

#selected cases
{
#all cases
casesPre<-unique(with(df2,interaction(type,stain,time,channel,sep = "_")))%>%
  as.data.frame()%>%
  separate(col=1,into=c("type","stain","time"),sep="_")
casesPre

#selected cases
cases1<-data.frame(
  type=c("PI","PI","PI","SYBR1","SYBR1","SYBR1","SYBR2","SYBR2","SYBR2","unstained","unstained"),
  stain=c(1,2,4,1,2,4,1,2,4,0,0),
  time=c(rep(30,11)),
  channel=c(rep("asinh.FL3.A",3),rep("asinh.FL1.A",6),"asinh.FSC.A","asinh.SSC.A"),
  cutoff=c(5,5,5,5,5,5,5,5,5,5,5),
  lowPop=c(5,5,5,5,5,5,5,5,5,5,10),
  highPop=c(14,14,14,14,14,14,14,14,14,14,14)
  #lowPop=rep(6,11),
  #highPop=rep(14,11)
  )

cases<-rbind(cases1,cases1[1:9,])
cases$time[12:20]<-90
cases$lowPop[14]<-10;cases$highPop[14]<-12
}

##testing
casnumb=1
testmod<-getEmMod(inp.df = df2,inp.type = cases$type[casnumb],inp.stain = cases$stain[casnumb],
                  inp.time = cases$time[casnumb],inp.channel = cases$channel[casnumb],
                  cutoff = cases$cutoff[casnumb],lowPop = cases$lowPop[casnumb],
                  highPop = cases$highPop[casnumb])

casnumb<-9  

{
  par(mfrow=c(4,5),mar=c(1,1,1,1))
  for (casnumb in 1:nrow(cases)){
    
    getEmMod(inp.df=df2,cases$type[casnumb],cases$stain[casnumb],cases$time[casnumb],
             cases$channel[casnumb],cases$cutoff[casnumb],cases$lowPop[casnumb],
             cases$highPop[casnumb])%>%
      plot(whichplots=2)#;abline(v = cutoffs.list[casnumb])
    
  }
  par(mfrow=c(1,1))  
}

## for all samples:
models.list<-Map(function(type,stain,time,channel,cutoff,lowPop,highPop){
  mod.t<-getEmMod(df2,type,stain,time,channel,cutoff,lowPop,highPop)
},type=cases$type,stain=cases$stain,time=cases$time,channel=cases$channel,
cutoff=cases$cutoff,lowPop=cases$lowPop,highPop=cases$highPop)

#getting cutoffs
cutoffs.list<-sapply(models.list,function(x) getEmCutoff(x))

#plotting
getEmPlot(cutoffs.list[[1]],models.list[[1]])

## Supplemental 3
i=1
plot.list<-list()

for (i in 1:length(models.list)){
  plot.list[[i]]<-
    data.frame(x=models.list[[i]]$x,
               g1l=models.list[[i]]$posterior[,1],
               g2l=models.list[[i]]$posterior[,2])%>%
    gather("group","value",2:3)%>%
    #head%>%
    ggplot()+
    geom_histogram(aes(x=x,y=..ncount..),bins=30)+
    geom_line(aes(x,value,color=group))+
    geom_vline(aes_string(xintercept=cutoffs.list[i]),col="red")+ylab("")+
    #ggtitle(aes_string("aaa",))+
    xlab(paste("Stain:",cases$type[i],"Concentration:",cases$stain[i],
               "Time:",cases$time[i],"Channel:",cases$channel[i]))+
    scale_color_discrete(guide=FALSE)+
    #scale_y_continuous(labels = "",breaks="")+
    theme(axis.text.y=element_blank(),
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=10))+
    stat_function(fun = sdnorm, n = 999,
                  args = list(mean = models.list[[i]]$mu[1],
                              sd = models.list[[i]]$sigma[1],
                              lambda =models.list[[i]]$lambda[1]),
                  col="#F8766D",linetype="dashed")+
    stat_function(fun = sdnorm, n = 999,
                  args = list(mean = models.list[[i]]$mu[2],
                              sd = models.list[[i]]$sigma[2],
                              lambda =models.list[[i]]$lambda[2]),
                  col="#00BFC4",linetype="dashed")
  
}

plot.list[[4]]

do.call(plot_grid,c(plot.list,ncol = 4))

#plot_grid()
ggsave("fig/ps3_separation.pdf",width=16,height = 20)


## Get Group parameters
cases2<-cbind(cases,popsep=cutoffs.list)

getGroupParam(inp.df=df2,
              inp.type = cases$type[1],
              inp.stain = cases$stain[1],
              inp.time = cases$time[1],
              inp.channel = cases$channel[1],
              popsep = cases2$popsep[1])

allp.list.P<-Map(
  function(type,stain,time,channel,popsep) getGroupParam(
    inp.df=df2,inp.type=type, inp.stain=stain, inp.time=time, inp.channel=channel, popsep=popsep),
  type=cases$type,stain=cases$stain,time=cases$time,channel=cases$channel,popsep=cases2$popsep)

allp.list<-allp.list.P%>%do.call("rbind",.) 

#finding critical value
crit.val<-qnorm(0.001/2,lower.tail = TRUE)%>%
  abs()

tx<-allp.list%>%
  mutate(diffM=highM-lowM,
         #sdPx=((highC-1)*highSD^2)+((lowC-1)*lowSD^2)/(highC+lowC-2),
         sdP2=((highC-1)*highSD^2+(lowC-1)*lowSD^2)/(lowC+highC-2))%>%
  mutate(sem=sqrt((sdP2/highC)+(sdP2/lowC)))%>%
  mutate(CI95=crit.val*sem)


#Supplemental 2

s2<-tx%>%
  separate(1,into=c("type","stain","channel","time"),sep="_")%>%
  mutate(time=ifelse(time==30,30,240))%>%
  ggplot(aes(stain,diffM,fill=interaction(channel,type)))+
  geom_bar(stat="identity",position="dodge")+
  geom_errorbar(aes(ymin=diffM-CI95,ymax=diffM+CI95),width=0.5,position="dodge")+
  ylab("Difference spore/cell distributions")+
  facet_grid(time~.)+
  scale_fill_discrete(name="")


s2

viability <- read_csv("data/viability2.csv")%>%
  gather("tripl","value",5:7)

s233<-viability%>%
  #mutate(conc=fct_reorder(conc))
  dplyr::filter(!is.na(value))%>%
  #mutate()%>%
  ggplot(aes(fct_inorder(conc),as.numeric(value)/210,col=stain))+
  geom_point(aes(shape=type),position=position_dodge(width = 0.0))+
  geom_line(stat="summary",fun.y="mean",
            aes(as.numeric(fct_inorder(conc)),as.numeric(value)/210,
                linetype=type))+
  facet_grid(time~.)+
  ylab("% CFU")+xlab("")+
  scale_color_discrete(name="")+
  scale_shape_discrete(name="")+
  scale_linetype_discrete(name="")

plot_grid(s2,s233,nrow = 2,labels = c("A","B"))

ggsave("fig/ps2_opt.pdf",width = 10, height = 8,dpi=300,units = "in")











#Only: time 30 stain 2 and time 30 stain 0
#Figure 2

fig2.1<-tx%>%
  separate(1,into=c("type","stain","channel","time"),sep="_")%>%
  dplyr::filter(time==30&stain==2 |time==30&stain==0)%>%
  ggplot(aes(interaction(type,channel),diffM,fill=channel))+
  geom_bar(stat="identity",position="dodge",alpha=0.7)+
  geom_errorbar(aes(ymin=diffM-CI95,ymax=diffM+CI95),width=0.5)+xlab("")+
  ylab("Difference spore/cell distributions")+scale_color_discrete(name="")+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position = c(0.6,0.8),axis.text.x = element_text(angle=30,hjust=1))+
  scale_x_discrete(label=c("SYBR1","SYBR2","PI","FSC","SSC"))

fig2.1
##plot: ggridges with separator, translate x axis of distribution, so that the line is everywhere the same

{
  f2.2<-df1%>% #ex df1
    dplyr::filter(ident=="Cells+spores")%>%
    select(type,stain,time,asinh.FL1.A,asinh.FL3.A,asinh.FSC.A,asinh.SSC.A)%>%
    dplyr::filter(time==30&stain==2 |time==30&stain==0)
  #dplyr::select(asinh.FL1.A,asinh.FL3.A,asinh.FSC.A,asinh.SSC.A,type)
  #dplyr::filter(channel %in% c("asinh.FL1.A","asinh.FL3.A","asinh.FSC.A","asinh.SSC.A"))
  
  # NOT FURTHER FOLLOWED UP SINCE OTHER DISTRIBUTIONS ARE SHOWN!
  
  allp.list[5,]
  f2.2$asinh.FL1.A[f2.2$type=="SYBR1"]<-f2.2$asinh.FL1.A[f2.2$type=="SYBR1"]-
    allp.list$popsep[5]
  
  allp.list[8,]
  f2.2$asinh.FL1.A[f2.2$type=="SYBR2"]<-f2.2$asinh.FL1.A[f2.2$type=="SYBR2"]-
    allp.list$popsep[8]
  
  allp.list[10,]
  allp.list[11,]
  f2.2$asinh.FSC.A<-f2.2$asinh.FSC.A-allp.list$popsep[10]
  f2.2$asinh.SSC.A<-f2.2$asinh.SSC.A-allp.list$popsep[11]
  
  allp.list[2,]
  f2.2$asinh.FL3.A<-f2.2$asinh.FL3.A-allp.list$popsep[2]
  
}

fig2.2<-f2.2%>%
  #select(-"type")%>%
  #with(.,table("asinh.FSC.A"))%>%
  gather("channel","value",4:7)%>%
  dplyr::mutate(channel=paste(channel,type,sep="."))%>%
  dplyr::select(-type)%>%
  dplyr::filter(channel %in% c("asinh.FL1.A.SYBR1",
                               "asinh.FL1.A.SYBR2",
                               "asinh.FL3.A.PI",
                               "asinh.SSC.A.unstained",
                               "asinh.FSC.A.unstained"))%>%
  ggplot(aes(x=value,y=channel))+
  geom_density_ridges(alpha=0.5)+
  geom_vline(aes(xintercept=0),col="red")+
  xlim(c(-5,5))+
  xlab("translated scatter/fluorescence signal")+ylab("")


fig2.2

plot_grid(fig2.1,fig2.2,rel_widths = c(0.33,0.66))
ggsave("fig/f2_SepDist.pdf",width = 10, height = 4)


#### Figure 3: ####

{
  sample.var<-c("strain","time","tripl","trash")
  #simple plot
  fcsset<-createFlowset(filepath = "data/MKA_Spores_Separation_9/",
                        sample.variables=sample.var)
  
  fcssetG<-Subset(fcsset,norm2Filter("asinh.FSC.H","asinh.FSC.W",scale.factor = 2))
  
  autoplot(fcsset,"asinh.FSC.H","asinh.SSC.A")
  autoplot(fcssetG,"asinh.FSC.H","asinh.SSC.A")
  
  autoplot(fcsset,"asinh.FSC.H","asinh.FSC.W")
  autoplot(fcssetG,"asinh.FSC.H","asinh.FSC.W")
  
  autoplot(fcsset,"asinh.FL1.A")
  autoplot(fcssetG,"asinh.FL1.A")
  
  #to df
  df1<-fcstodf(fcssetG)%>%
    separate(length(fcsset@colnames)+1,into = sample.var,sep = "_")
  
  #preprocessing
  df2<-df1%>%
    dplyr::filter(asinh.SSC.A<14.4,
                  asinh.FL1.A>3,
                  asinh.SSC.A>8)
  
  df1%>%
    #dplyr::filter(!stain=="unstained")%>%
    #dplyr::filter(stain==1)%>%
    ggplot(aes(asinh.SSC.A,asinh.FL1.A))+
    #geom_density_ridges()+
    geom_hex(bins=300)+
    facet_grid(strain~time)+
    xlim(c(10,15))+ylim(c(2.5,14))+
    scale_fill_viridis()+
    theme_bw()
  
  df1%>%
    #dplyr::filter(!stain=="unstained")%>%
    #dplyr::filter(stain==1)%>%
    ggplot(aes(asinh.SSC.A,asinh.FL1.A))+
    #geom_density_ridges()+
    geom_hex(bins=300)+
    facet_grid(strain~time)+
    #xlim(c(10,15))+ylim(c(2.5,14))+
    scale_fill_viridis()+
    theme_bw()
  
  #ggsave("fig/S5-applic.pdf",width = 12, height = 7)
  
  #Create Reference
  dfap.2<-df1%>%  
    #dplyr::filter(!stain=="unstained")%>%
    dplyr::filter(strain=="Bs02003")%>%
    #dplyr::filter(time==tim)%>%
    #dplyr::filter(stain==stn)%>%
    dplyr::select(asinh.FSC.A,asinh.SSC.A,asinh.FL1.A)%>%
    as.matrix()
  
  matrix.start<-matrix(c(12,14,10,
                         10,13,11,
                         11,11.2,7.5),ncol=3)
  
  Appl.kmeans<-kmeans(dfap.2,matrix.start)
  
  data.frame(dfap.2,cluster=Appl.kmeans$cluster)%>%
    ggplot(aes(asinh.SSC.A,asinh.FL1.A))+
    #geom_point(aes(col=cluster))+
    geom_hex(aes(fill=as.factor(cluster)),bins=300)+
    xlim(c(10,15))+ylim(c(2.5,15))+
    scale_fill_viridis(discrete = TRUE,guide=FALSE,end=0.8)
  
  #ggsave("fig/F6S-Sep-Kmeans.pdf",width = 12, height = 7)
  #looks pretty good
  
  centers.list<-list(Appl.kmeans$centers)
  
  ###splitting and applying kmeans centers to predict
  
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
  
  df3<-df1%>%dplyr::select(asinh.FSC.A,asinh.SSC.A,asinh.FL1.A,strain,time,tripl)
  df3.list<-split(df3,df3$strain)
  
  clusterA.pred<-base::Map(function(x,y) {
    cbind(x,cluster=predict.clusters2(x[,c(-4,-5,-6)],y))%>%
      as.data.frame()
  },
  df3.list,centers.list)%>%
    do.call(rbind,.)
  
  clusterA.pred%>%
    #dplyr::filter(!stain=="unstained")%>%
    #dplyr::filter(stain==x)%>%
    ggplot(aes(asinh.SSC.A,asinh.FL1.A))+
    #geom_density_ridges()+
    geom_hex(aes(fill=as.factor(cluster)),bins=200)+
    facet_grid(strain~tripl)+
    xlim(c(10,15))+ylim(c(2.5,14))+
    #scale_fill_viridis()+
    theme_bw()
  
  #final figure
  c.count<-clusterA.pred%>%
    group_by(strain,time,cluster)%>%
    summarize(count=n())
  #dplyr::filter(stain!="unstained")
  #ungroup()%>%
  
  #group_by(strain,time,cluster)%>%
  #summarize(count.mean=mean(count),count.sd=sd(count))
}

#other:910

ccount2<-c.count%>%
  ggplot(aes(strain,count,fill=as.factor(cluster)))+
  #stat_summary(fun.data=mean_sdl,fun.args = list(mult=1),geom="errorbar", width=1)+
  geom_bar(stat="summary",fun.y="mean",size=1,position=position_fill())+
  #geom_line(stat="summary",fun.y="mean",size=0.5,position=position_stack())+
  ylab("Proportion / %")+xlab("time / h")+
  scale_y_continuous(expand = c(0,0),
                     trans = trans_new("percent",function(x) x/100,function(x) x*100))+
  scale_fill_viridis(labels=c("Cells","Endospores","Spores"),discrete=TRUE,end = c(0.8),direction = -1,
                     guide=FALSE,alpha = 0.8)+
  theme_minimal()+
  theme(panel.spacing = unit(1, "lines"),legend.title=element_blank())+
  xlab("")

ccount2
plot_grid(ccount1,ccount2,nrow = 1,rel_widths = c(0.7,0.3),labels = c("A","B"))

ggsave("fig/f4_applications.pdf",width = 8, height = 5)

#at t=24 hours

```

```{r}

clplot1<-data.frame(dfap.2,cluster=Appl.kmeans$cluster)%>%
  ggplot(aes(asinh.SSC.A,asinh.FL1.A))+
  #geom_point(aes(col=cluster))+
  geom_hex(aes(fill=as.factor(cluster)),bins=300)+
  xlim(c(10,15))+ylim(c(2.5,15))+
  scale_fill_viridis(discrete = TRUE,end=0.8,label=c("Spores","Endospores","Cells"),name="")+
  theme_bw()
clplot1
ggsave("fig/f4p1_class_fl.pdf",width = 5, height = 3)

clplot2<-data.frame(dfap.2,cluster=Appl.kmeans$cluster)%>%
  ggplot(aes(asinh.SSC.A,asinh.FSC.A))+
  #geom_point(aes(col=cluster))+
  geom_hex(aes(fill=as.factor(cluster)),bins=300)+
  xlim(c(10,15))+ylim(c(9,12))+
  scale_fill_viridis(discrete = TRUE,end=0.8,label=c("Spores","Endospores","Cells"),name="")+
  theme_bw()
clplot2
ggsave("fig/f4p2_class_fsc.pdf",width = 5, height = 3)
