source("src/functions.R")

#### Figure 1 - Separation of Cells vs Spores possible, overview plot ####
#loading
sample.var<-c("strain","ident","stain","mode","run")
fcsset1<-flowCreateFlowSet(filepath = "data/f1/",sample_variables =sample.var)

#gating
fcsset1.g<-fcsset1%>%
  Subset(.,norm2Filter("asinh.FSC.H", "asinh.FSC.W",filterId="norm_ssc.fsc",scale=1))%>%
  Subset(.,norm2Filter("asinh.SSC.H", "asinh.SSC.W",filterId="norm_ssc.fsc",scale=1))

df1<-flowFcsToDf(fcsset1.g)

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

#FSC-Plot
FSC.plot<-df1%>%
  dplyr::filter(mode=="FSC")%>%
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
  
  sample.var=c("ident","type","stain","t1","t2","t3","time") 
  fcsset2.1<-flowCreateFlowSet(filepath = "data/f2/set5_cells+spores_30/",
                               sample_variables=sample.var,additional_variable = 30,transformation = TRUE)
  fcsset2.2<-flowCreateFlowSet(filepath = "data/f2/set5_cells+spores_90/",
                               sample_variables =sample.var,additional_variable = 90,transformation = TRUE)
  fcsset2.5<-rbind2(fcsset2.1,fcsset2.2)
}

#gating check

fcsset2.5%>%
  autoplot("asinh.FSC.H","asinh.FSC.W",bins=200)+
  theme_minimal()

fcsset2.5%>%
  Subset(.,norm2Filter("FSC-H", "FSC-W",filterId="norm_ssc.fsc",scale.factor = 5))%>%
  autoplot("asinh.FSC.H","asinh.FSC.W",bins=200)+
  theme_minimal()

fcsset2.5g<-fcsset2.5%>%
  Subset(.,norm2Filter("FSC-H", "FSC-W",filterId="norm_ssc.fsc",scale=5))

#to df
df2<-flowFcsToDf(fcsset2.5g)%>%
  dplyr::filter(ident=="Cells+spores")%>%
  select(type,stain,time,asinh.FL1.A,asinh.FL3.A,asinh.FSC.A,asinh.SSC.A)%>%
  gather("channel","value",4:7)

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
  )

cases<-rbind(cases1,cases1[1:9,])
cases$time[12:20]<-90
cases$lowPop[14]<-10;cases$highPop[14]<-12
}

## gm for all models
models.list<-Map(function(type,stain,time,channel,cutoff,lowPop,highPop){
  mod.t<-getEmMod(df2,type,stain,time,channel,cutoff,lowPop,highPop)
},type=cases$type,stain=cases$stain,time=cases$time,channel=cases$channel,
cutoff=cases$cutoff,lowPop=cases$lowPop,highPop=cases$highPop)

#getting cutoffs
cutoffs.list<-sapply(models.list,function(x) getEmCutoff(x))

#plotting

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

do.call(plot_grid,c(plot.list,ncol = 4))

#plot_grid()
ggsave("suppl/ps3_separation.pdf",width=16,height = 24)
ggsave("suppl/ps3_separation.png",width=16,height = 24)

## Get Group parameters

{
cases2<-cbind(cases,popsep=cutoffs.list)

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
}


#Supplemental 1

ps2.1<-tx%>%
  separate(1,into=c("type","stain","channel","time"),sep="_")%>%
  mutate(time=ifelse(time==30,30,240))%>%
  ggplot(aes(stain,diffM,fill=interaction(channel,type)))+
  geom_bar(stat="identity",position=position_dodge(),alpha=0.8)+
  geom_errorbar(aes(ymin=diffM-CI95,ymax=diffM+CI95,group=interaction(channel,type)),
                position=position_dodge(),alpha=0.8)+
  ylab("Difference spore/cell distributions")+xlab("")+
  facet_grid(time~.)+
  scale_fill_discrete(name="")+
  scale_x_discrete(label=c("unstained","1x","2x","4x"))

ps2.1

viability <- read_csv("data/viability2.csv")%>%
  gather("tripl","value",5:7)

ps2.2<-viability%>%
  #mutate(conc=fct_reorder(conc))
  dplyr::filter(!is.na(value))%>%
  #mutate()%>%
  ggplot(aes(fct_inorder(conc),as.numeric(value)/210,col=stain))+
  geom_point(aes(shape=type),position=position_dodge(width = 0.0))+
  geom_line(stat="summary",fun.y="mean",
            aes(as.numeric(fct_inorder(conc)),as.numeric(value)/210,
                linetype=type),alpha=0.8)+
  facet_grid(time~.)+
  ylab("% CFU")+xlab("")+
  scale_color_discrete(name="")+
  scale_shape_discrete(name="")+
  scale_linetype_discrete(name="")

plot_grid(ps2.1,ps2.2,nrow = 2,labels = c("A","B"))

ggsave("suppl/ps1_opt.pdf",width = 8, height = 10,dpi=300,units = "in")
ggsave("suppl/ps1_opt.png",width = 8, height = 10,dpi=300,units = "in")

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

{
  f2.2<-df2%>%
    group_by(type,stain,time)%>%
    mutate(id = row_number())%>%
    tidyr::spread(channel,value,fill=NA)%>%
    ungroup()%>%
    dplyr::select(-id)
    
  #translate so that cutoff is at 0
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
  gather("channel","value",4:7)%>%
  dplyr::mutate(channel=paste(channel,type,sep="."))%>%
  dplyr::select(-type)%>%
  dplyr::filter(channel %in% c("asinh.FL1.A.SYBR1",
                               "asinh.FL1.A.SYBR2",
                               "asinh.FL3.A.PI",
                               "asinh.SSC.A.unstained",
                               "asinh.FSC.A.unstained"))%>%
  dplyr::filter(!is.na(value))%>%
  ggplot(aes(x=value,y=channel))+
  geom_density_ridges(alpha=0.5)+
  geom_vline(aes(xintercept=0),col="red")+
  xlim(c(-5,5))+
  xlab("translated scatter/fluorescence signal")+ylab("")


fig2.2

plot_grid(fig2.1,fig2.2,rel_widths = c(0.33,0.66),labels = c("A","B"))
ggsave("fig/f2_sepdist.pdf",width = 10, height = 5)

#### Figure 3+4: Clustering ####

## Set9: Figure 3 and 4B

{
  rm(list=ls())
  source("src/functions.R")
  sample.var<-c("strain","time","tripl")
  #simple plot
  fcsset3<-flowCreateFlowSet(filepath = "data/f3/set9/",sample_variables=sample.var,transformation = TRUE)
  
  fcsset3G<-fcsset3%>%
    Subset(norm2Filter("asinh.FSC.H","asinh.FSC.W",scale.factor = 2))%>%
    Subset(rectangleGate("asinh.FL1.A"=c(0,15),"asinh.FL3.A"=c(0,15)))
  
  #to df
  df3<-fcsset3G%>%
    flowFcsToDf()%>%
    separate(length(fcsset3G@colnames)+1,into = sample.var,sep = "_")
}
  #Get centers
{
  df3.ref<-df3%>%  
    dplyr::filter(strain=="Bs02003")%>%
    dplyr::select(asinh.FSC.A,asinh.SSC.A,asinh.FL1.A)%>%
    as.matrix()
  
  #finding centers: with random starting points
  set.seed(1)
  Appl.kmeans <- kmeans.rep(df3.ref,rep = 1000)
  
  #center positions
  centers.list<-list(Appl.kmeans$centers)
  centers.list.df<-as.data.frame(centers.list)
  
  write.csv(centers.list.df,"suppl/centers_f3.csv")

  center.locs<-factor(Appl.kmeans$cluster,levels=c(2,1,3))
  }

#Figure 3A
#at t=24 hours

#PCA
pca.res<-df3.ref%>%
  prcomp(center=TRUE,scale. = TRUE)

pca.res$x%>%
  as.tibble()%>%
  mutate(cluster=factor(Appl.kmeans$cluster,levels=c(2,1,3)))%>%
  ggplot(aes(PC1,PC2))+
  #geom_point(aes(col=cluster,alpha=0.6))+
  geom_hex(aes(fill=..count..),bins=200)+
  scale_color_viridis(discrete = TRUE,end=0.8,label=c("Cells","Endospores","Spores"),
                     name="",direction = -1)+
  scale_alpha_continuous(guide=FALSE)+
  scale_fill_viridis_c(guide=FALSE)+
  geom_density2d(aes(color=cluster),bins=3)

ggsave("fig/f3_pca.png",width = 6, height = 6,dpi = 900)

clplot1<-data.frame(df3.ref,cluster=factor(Appl.kmeans$cluster,levels=c(2,1,3)))%>%
  ggplot(aes(asinh.SSC.A,asinh.FL1.A))+
  #geom_point(aes(col=cluster))+
  geom_hex(aes(fill=as.factor(cluster)),bins=300)+ #,alpha=..ncount.. #order= ?
  xlim(c(10,15))+ylim(c(2.5,15))+
  scale_fill_viridis(discrete = TRUE,end=0.8,label=c("Cells","Endospores","Spores"),name="",direction = -1,
                     guide = FALSE)+
  scale_alpha_continuous(guide = FALSE)+
  theme_bw()+
  geom_point(aes(centers.list.df[1,2],centers.list.df[1,3]),col="red",size=1)+
  geom_point(aes(centers.list.df[2,2],centers.list.df[2,3]),col="red",size=1)+
  geom_point(aes(centers.list.df[3,2],centers.list.df[3,3]),col="red",size=1)

clplot1

clplot2<-data.frame(df3.ref,cluster=factor(Appl.kmeans$cluster,levels=c(2,1,3)))%>%
  ggplot(aes(asinh.SSC.A,asinh.FSC.A))+
  #geom_point(aes(col=cluster))+
  geom_hex(aes(fill=as.factor(cluster)),bins=300)+ #,alpha=..ncount..
  xlim(c(10,15))+ylim(c(9,12))+
  scale_fill_viridis(discrete = TRUE,end=0.8,label=c("Cells","Endospores","Spores"),name="",direction = -1,
                     guide=FALSE)+
  scale_alpha_continuous(guide = FALSE)+
  theme_bw()+
  geom_point(aes(centers.list.df[1,2],centers.list.df[1,1]),col="red",size=1)+
  geom_point(aes(centers.list.df[2,2],centers.list.df[2,1]),col="red",size=1)+
  geom_point(aes(centers.list.df[3,2],centers.list.df[3,1]),col="red",size=1)

clplot2

plot_grid(clplot1,clplot2,align = "h")
#ggsave("fig/f3_class.pdf",width = 8, height = 4)
#ggsave("fig/f3_class.svg",width = 8, height = 4)
ggsave("fig/f3_class.png",width = 8, height = 4,dpi = 900)


#Figure 3
#splitting according to kmeans centers
{
  df3.split<-df3%>%dplyr::select(asinh.FSC.A,asinh.SSC.A,asinh.FL1.A,strain,time,tripl)
  df3.list<-split(df3.split,df3$strain)
  
  clusterB.pred<-base::Map(function(x,y) {
    cbind(x,cluster=predict.clusters2(x[,c(-4,-5,-6)],y))%>%
      as.data.frame()
  },
  df3.list,centers.list)%>%
    do.call(rbind,.)
}

##Supplemental 4

clusterB.pred%>%
  #dplyr::filter(!stain=="unstained")%>%
  #dplyr::filter(stain==x)%>%
  ggplot(aes(asinh.SSC.A,asinh.FL1.A))+
  #geom_density_ridges()+
  geom_hex(aes(),bins=200)+
  facet_grid(tripl~strain)+
  xlim(c(10,15))+ylim(c(2.5,14))+
  scale_fill_viridis(guide=FALSE)+
  theme_bw()

ggsave("suppl/ps4_f4b.pdf",width = 4, height = 6)
#ggsave("suppl/ps4_f4b.png",width = 4, height = 6)

#FIGURE 4B

{
c.count.B<-clusterB.pred%>%
  group_by(strain,time,cluster,tripl)%>%
  summarize(count=n())%>%
  ungroup()%>%
  group_by(strain,time,cluster)%>%
  summarize(mean.count=mean(count),sd.count=sd(count))%>%
  ungroup()%>%
  group_by(strain,time)%>%
  mutate(perc.mean.count=100*mean.count/sum(mean.count),perc.sd.count=100*sd.count/sum(mean.count))

ccount2<-c.count.B%>%
  ggplot(aes(strain,perc.mean.count,fill=as.factor(cluster)))+
  geom_bar(size=1,position=position_dodge(),stat="identity")+
  geom_errorbar(aes(ymin=perc.mean.count-perc.sd.count,ymax=perc.mean.count+perc.sd.count,
                    group=as.factor(cluster)),position=position_dodge(),alpha=0.6)+
  ylab("Proportion / %")+xlab("time / h")+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_viridis(labels=c("Cells","Endospores","Spores"),discrete=TRUE,end = c(0.8),direction = -1,
                     alpha = 0.8)+
  theme_minimal()+
  theme(panel.spacing = unit(1, "lines"),legend.title=element_blank(),
        legend.position = c(0.44,0.82))+
  xlab("")

ccount2

}

#### Figure 4A

{
  sample.var <- c("strain","time","stain","tripl")
  fcsset3A<-flowCreateFlowSet(filepath = "data/f3/set8/",
                          sample_variables=sample.var,transformation = TRUE)
    
  fcsset3AG<-fcsset3A%>%
    Subset(norm2Filter("asinh.FSC.H","asinh.FSC.W",scale.factor = 2))%>%
    Subset(rectangleGate("asinh.FL1.A"=c(0,15),"asinh.FL3.A"=c(0,15)))
  
#to df
df3A<-fcsset3AG%>%
  flowFcsToDf()
}

{
  ## Use kmeans to identify subsets
  #selecting reference 
  Bs2003.reference<-ref.ds(strn="Bs02003",tim="48h",stn="1",df=df3A)
  Bs2018.reference<-ref.ds(strn="Bs02018",tim="56h",stn="2",df=df3A)
  Bs2020.reference<-ref.ds(strn="Bs02020",tim="56h",stn="1",df=df3A)
  
  set.seed(1)
  Bs2003.kmeans<-kmeans.rep(Bs2003.reference,rep = 1000)
  Bs2018.kmeans<-kmeans.rep(Bs2018.reference,rep = 1000)
  Bs2020.kmeans<-kmeans.rep(Bs2020.reference,rep = 1000)
  
}

#Supplemental ?

data.frame(Bs2003.reference,cluster=Bs2003.kmeans$cluster)%>%
    ggplot(aes(asinh.SSC.A,asinh.FL1.A))+
    geom_hex(aes(fill=log(..count..)),bins=300)+
    facet_grid(.~cluster)+xlim(c(10,15))+ylim(c(2.5,15))
  
data.frame(Bs2018.reference,cluster=Bs2018.kmeans$cluster)%>%
    ggplot(aes(asinh.SSC.A,asinh.FL1.A))+
    geom_hex(aes(fill=log(..count..)),bins=300)+
    facet_grid(.~cluster)+xlim(c(10,15))+ylim(c(2.5,15))
  
data.frame(Bs2020.reference,cluster=Bs2020.kmeans$cluster)%>%
    ggplot(aes(asinh.SSC.A,asinh.FL1.A))+
    geom_hex(aes(fill=log(..count..)),bins=300)+
    facet_grid(.~cluster)+xlim(c(10,15))+ylim(c(2.5,15))
  
  ###splitting and applying kmeans centers to predict
{
  centers.list<-list(Bs2003.kmeans$centers,
                     Bs2018.kmeans$centers,
                     Bs2020.kmeans$centers)
  
  df3A.list<-df3A%>%
    select(asinh.FSC.A,asinh.SSC.A,asinh.FL1.A,strain,stain,time)%>%
    split(.$strain)
  
  cluster.pred.A<-base::Map(function(x,y) {
    cbind(x,cluster=predict.clusters2(x[,c(-4,-5,-6)],y))%>%
      as.data.frame()
  },
  df3A.list,centers.list)%>%
    do.call(rbind,.)
  
}

#Supplemental ?

suppl6<-lapply(1:3,function(x) {
  cluster.pred.A%>%
    #dplyr::filter(!stain=="unstained")%>%
    dplyr::filter(stain==x)%>%
    ggplot(aes(asinh.SSC.A,asinh.FL1.A))+
    #geom_density_ridges()+
    geom_hex(aes(),bins=200)+
    scale_fill_viridis()+
    facet_grid(time~strain)+
    xlim(c(10,15))+ylim(c(2.5,14))+
    #scale_fill_viridis()+
    theme_bw()
})

#ggsave("suppl/ps4B_f4a.pdf",width = 4, height = 6)

lapply(1,function(x) {
  cluster.pred.A%>%
      #dplyr::filter(!stain=="unstained")%>%
      dplyr::filter(stain==x)%>%
      ggplot(aes(asinh.SSC.A,asinh.FL1.A))+
      #geom_density_ridges()+
      #geom_polygon()+
      geom_hex(aes(col=as.factor(cluster),fill=log(..count..)),bins=30)+
      facet_grid(strain~time)+
      xlim(c(10,15))+ylim(c(2.5,14))+
      #scale_fill_viridis()+
      theme_bw()+
      scale_fill_viridis()
  })

{
  c.count.A<-cluster.pred.A%>%
    group_by(strain,time,cluster,stain)%>%
    summarize(count=n())%>%
    ungroup()%>%
    group_by(strain,time,cluster)%>%
    summarize(mean.count=mean(count),sd.count=sd(count))%>%
    ungroup()%>%
    group_by(strain,time)%>%
    mutate(perc.mean.count=100*mean.count/sum(mean.count),perc.sd.count=100*sd.count/sum(mean.count))%>%
    separate(time,c("time","h"),2)%>%
    mutate(time=as.numeric(time))%>%
    #mutate(x=interaction(strain,time))%>%
    ungroup()%>%
    #switching group names
    mutate(cluster=replace(cluster,cluster==2&strain=="Bs02018",4),
           cluster=replace(cluster,cluster==2&strain=="Bs02020",4))%>%
    mutate(cluster=replace(cluster,cluster==3&strain=="Bs02018",2),
           cluster=replace(cluster,cluster==1&strain=="Bs02020",2))%>%
    mutate(cluster=replace(cluster,cluster==4&strain=="Bs02018",3),
           cluster=replace(cluster,cluster==3&strain=="Bs02020",1))%>%
    mutate(cluster=replace(cluster,cluster==4&strain=="Bs02020",3))
  
  c.count.A%>%write.csv(file = "suppl/population_count.csv")
  
  ccount1<-c.count.A%>%
    ggplot(aes(as.factor(time),perc.mean.count))+
    #geom_line(size=1,position=position_dodge(),stat="identity",aes(color=as.factor(cluster)))+
    geom_bar(size=1,position=position_dodge(),stat="identity",aes(fill=as.factor(cluster)))+
    geom_errorbar(aes(ymin=perc.mean.count-perc.sd.count,ymax=perc.mean.count+perc.sd.count,
                      group=as.factor(cluster)),position=position_dodge(),alpha=0.6)+
    ylab("Proportion / %")+xlab("time / h")+
    scale_y_continuous(expand = c(0,0))+
    facet_grid(strain~.)+
    scale_fill_viridis(labels=c("Cells","Endospores","Spores"),discrete=TRUE,end = c(0.8),direction = -1,
                       guide=FALSE,alpha = 0.8)+
    theme_minimal()+
    theme(panel.spacing = unit(1, "lines"),legend.title=element_blank())+
    xlab("time / h")
  
  ccount1
}

plot_grid(ccount1,ccount2,nrow = 1,rel_widths = c(0.6,0.4),labels = c("A","B"))

ggsave("fig/f4_applications.pdf",width = 8, height = 5)

#### Supplemental: Combinations ####

#GET: mka spores sepration 4 > combination stains!
source("src/functions.R")

sample.var=c("ident","type","stain","conc","trash")
fcsset5<-flowCreateFlowSet(filepath = "data/reg/",sample_variables=sample.var)

#Gating
fcsset5g<-fcsset5%>%
  Subset(.,norm2Filter("FSC-H", "FSC-W",filterId="norm_ssc.fsc",scale=1))%>%
  Subset(rectangleGate("asinh.FL1.A"=c(0,15),"asinh.FL3.A"=c(0,15)))

fcs.df1<-flowFcsToDf(fcsset5g)

fcs.df1%>%
  ggplot(aes(asinh.FL1.A,interaction(type,stain,conc)))+
  geom_density_ridges()+ylab("")

fcs.df2<-fcs.df1%>%
  group_by(type,stain,conc)%>%
  summarize(mean.fl1=mean(asinh.FL1.A),mean.fl3=mean(asinh.FL3.A))%>%
  ungroup()%>%
  mutate(PI=   c(0,2,1,1,0,0,0,0,2,1,1,0,0,0),
         SYBR1=c(0,0,1,0,2,1,0,0,0,1,0,2,1,0),
         SYBR2=c(0,0,0,1,0,1,2,0,0,0,1,0,1,2))

fcs.df3%>%
  ggplot(aes(as.numeric(SYBR1),as.numeric(mean.fl1)))+
  geom_point(aes(col=PI,shape=SYBR2))+
  geom_line(aes(col=PI),stat="summary",fun.y=mean)

fcs.df3%>%
  ggplot(aes(as.numeric(mean.fl1),as.numeric(mean.fl3)))+
  geom_point(aes(col=interaction(PI),shape=SYBR2))+
  geom_line(aes(col=SYBR1),stat="summary",fun.y=mean)+
  facet_grid(PI~.)

fcs.df3%>%
  dplyr::filter(PI==0)%>%
  lm(as.numeric(mean.fl1)~-1+as.numeric(SYBR2)*as.numeric(SYBR1)+
       as.numeric(SYBR2),data=.)%>%
  summary()

fcs.df3%>%
  dplyr::filter(SYBR2==0)%>%
  lm(as.numeric(mean.fl3)~-1+as.numeric(SYBR1)*as.numeric(PI),data=.)%>%
  summary()

fcs.df3%>%
  dplyr::filter(SYBR1==0)%>%
  lm(as.numeric(mean.fl3)~-1+as.numeric(SYBR2)*as.numeric(PI),data=.)%>%
  summary()

write.csv2(fcs.df3,"suppl/mixdesign.csv")
