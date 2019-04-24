source("src/functions.R")

#### Figure 1 - Separation of Cells vs Spores possible, overview plot ####
{
  #sample variables
  sample.var<-c("strain","ident","stain","mode","run")
  
  df1<-flowCreateFlowSet(filepath = "data/f1/",sample_variables =sample.var)%>%
    Subset(.,norm2Filter("asinh.FSC.H", "asinh.FSC.W",filterId="norm_ssc.fsc",scale=1))%>%
    Subset(.,norm2Filter("asinh.SSC.H", "asinh.SSC.W",filterId="norm_ssc.fsc",scale=1))%>%
    flowFcsToDf(.)
}

{
SSC.plot<-df1%>%
  dplyr::filter(mode=="SSC")%>%
  select(asinh.SSC.A,ident,stain)%>% #,run
  #dplyr::filter(stain=="unstained",ident %in% c("spores","cells"))%>%
  ggplot(aes(asinh.SSC.A))+
  geom_density(aes(fill=ident,y=..scaled..),alpha=0.4)+
  ylab("norm. count")+
  guides(fill=FALSE)+theme(legend.position = c(.05,.9))+
  scale_fill_viridis(discrete = TRUE,alpha=0.7,begin = 0,end = 0.8,name="",direction = -1)

FSC.plot<-df1%>%
  dplyr::filter(mode=="FSC")%>%
  select(asinh.FSC.A,ident,stain,run)%>%
  ggplot(aes(asinh.FSC.A))+
  geom_density(aes(fill=ident,y=..scaled..),alpha=0.4)+
  ylab("")+
  guides(fill=FALSE)+
  scale_fill_viridis(discrete = TRUE,alpha=0.7,begin = 0,end = 0.8,direction = -1)

PI.plot<-df1%>%
  dplyr::filter(mode=="PI")%>%
  select(asinh.FL3.A,ident,stain,run)%>%
  #dplyr::filter(stain=="PI",run=="2x")%>%
  ggplot(aes(asinh.FL3.A))+
  geom_density(aes(fill=ident,y=..scaled..),alpha=0.4)+
  guides(fill=FALSE)+
  scale_fill_viridis(discrete = TRUE,alpha=0.7,begin = 0,end = 0.8,direction = -1)+
  xlim(c(5,12))+
  # annotate(geom = "segment", x = 9, xend = 8.6, y = 0.3, yend = 0.1,
  #          arrow=arrow(),colour = "red",alpha=0.7,size=1)+
  ylab("")
 
PI.plot

S1.plot<-df1%>%
  dplyr::filter(mode=="SYBR1")%>%
  select(asinh.FL1.A,ident,stain,run)%>%
  ggplot(aes(asinh.FL1.A))+
  geom_density(aes(fill=ident,y=..scaled..),alpha=0.4)+
  guides(fill=FALSE)+
  ylab("norm. count")+xlim(c(5,14))+
  scale_fill_viridis(discrete = TRUE,alpha=0.7,begin = 0,end = 0.8,direction = -1)

S2.plot<-df1%>%
  dplyr::filter(mode=="SYBR2")%>%
  select(asinh.FL1.A,ident,stain,run)%>%
  ggplot(aes(asinh.FL1.A))+
  geom_density(aes(fill=ident,y=..scaled..),alpha=0.4)+
  guides(fill=FALSE)+
  ylab("norm. count")+xlim(c(5,14))+
  scale_fill_viridis(discrete = TRUE,alpha=0.7,begin = 0,end = 0.8,direction = -1)

legend.plot<-df1%>%
  select(asinh.FL1.A,ident,stain,run)%>%
  ggplot(aes(asinh.FL1.A))+geom_density(aes(fill=ident),alpha=0.4)+
  ylab("")+xlab("")+
  scale_fill_viridis(discrete = TRUE,alpha=0.7,begin = 0,end = 0.8,
                     labels=c("non-sporulating cells","purified spores"),
                     name="",direction = -1)+
  theme(legend.position = c(0.01,0.95))+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  ylim(0,100)+ 
  annotate(geom="text", x=-10, y=65, parse=TRUE,label = "bold('A')~Side~Scatter~(no~stain)",hjust=0,size=4)+
  annotate(geom="text", x=-10, y=55, parse=TRUE,label = "bold('B')~Forward~Scatter~(no~stain)",hjust=0,size=4)+
  annotate(geom="text", x=-10, y=45, parse=TRUE,label = "bold('C')~PI~stain",hjust=0,size=4)+
  annotate(geom="text", x=-10, y=35, parse=TRUE,label = "bold('D')~SYBR1~stain",hjust=0,size=4)+
  annotate(geom="text", x=-10, y=25, parse=TRUE,label = "bold('E')~SYBR2~stain",hjust=0,size=4)+
  geom_rect(aes(xmin = -10, xmax = 100, ymin = 0, ymax = 1),color = "white", size = 2, fill = "white")

plot_grid(SSC.plot,FSC.plot,PI.plot,S1.plot,S2.plot,legend.plot,ncol = 3,labels = "AUTO")
}

ggsave("fig/Figure_1.pdf",width = 9, height = 5,units = "in",dpi=300)

#### Figure 2: - Measuring cells and spores together: ability to separate with GMM ####

{#loading
  rm(list=ls())
  source("src/functions.R")
  
  sample.var=c("ident","type","stain","t1","t2","t3","time") 
  fcsset2.1<-flowCreateFlowSet(filepath = "data/f2/set5_cells+spores_30/",
                               sample_variables=sample.var,additional_variable = 30,transformation = TRUE)
  fcsset2.2<-flowCreateFlowSet(filepath = "data/f2/set5_cells+spores_90/",
                               sample_variables =sample.var,additional_variable = 90,transformation = TRUE)

  df2<-fcsset2.5<-rbind2(fcsset2.1,fcsset2.2)%>%
  Subset(.,norm2Filter("FSC-H", "FSC-W",filterId="norm_ssc.fsc",scale=5))%>%
  flowFcsToDf(.)%>%
  dplyr::filter(ident=="Cells+spores")%>%
  select(type,stain,time,asinh.FL1.A,asinh.FL3.A,asinh.FSC.A,asinh.SSC.A)%>%
  gather("channel","value",4:7)
}

{
#selected cases and set starting points for EMM depending on channel
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

## gmm for all samples
models.list<-Map(function(type,stain,time,channel,cutoff,lowPop,highPop){
  mod.t<-getEmMod(df2,type,stain,time,channel,cutoff,lowPop,highPop)
},type=cases$type,stain=cases$stain,time=cases$time,channel=cases$channel,
cutoff=cases$cutoff,lowPop=cases$lowPop,highPop=cases$highPop)

#getting cutoffs
cutoffs.list<-sapply(models.list,function(x) getEmCutoff(x))

## Supplemental figure 3
i=1;plot.list<-list()
sdnorm <- function(x, mean=0, sd=1, lambda=1){lambda*dnorm(x, mean=mean, sd=sd)}

for (i in 1:length(models.list)){
  plot.list[[i]]<-
    data.frame(x=models.list[[i]]$x,
               g1l=models.list[[i]]$posterior[,1],
               g2l=models.list[[i]]$posterior[,2])%>%
    gather("group","value",2:3)%>%
    ggplot()+
    geom_histogram(aes(x=x,y=..ncount..),bins=30)+
    geom_line(aes(x,value,color=group))+
    geom_vline(aes_string(xintercept=cutoffs.list[i]),col="red")+ylab("")+
    xlab(paste("Stain:",cases$type[i],"Concentration:",cases$stain[i],
               "Time:",cases$time[i],"Channel:",cases$channel[i]))+
    scale_color_discrete(guide=FALSE)+
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

ggsave("suppl/Supplemental3.pdf",width=16,height = 24)
ggsave("suppl/Supplemental3.png",width=16,height = 24)


#get distribution parameters
distr.values<-lapply(models.list,function(x) {
  cbind(t(x$mu),t(x$sigma))%>%
    return()
  })%>%
  do.call("rbind",.)%>%
  as_tibble()%>%
  magrittr::set_colnames(c("mu_1","mu_2","sd_1","sd_2"))%>%
  mutate(
    diffM=mu_2-mu_1, #difference between means
    pooledSD=sqrt((sd_1^2+sd_2^2)/2) #pooled standard deviation
  )%>%
  cbind(.,cases)
  
distr.values

PS2.A<-distr.values%>%
  ggplot(aes(stain,diffM,fill=interaction(channel,type)))+
  geom_bar(stat="identity",position=position_dodge(),alpha=0.8)+
  #geom_boxplot(stat="identity",position=position_dodge(),alpha=0.8)+
  geom_errorbar(aes(ymin=diffM-pooledSD,ymax=diffM+pooledSD,group=interaction(channel,type)),
                position=position_dodge(),alpha=0.8)+
  ylab("Difference spore/cell distributions")+xlab("")+
  facet_grid(time~.)+
  scale_fill_discrete(name="")+
  scale_x_discrete(label=c("unstained","1x","2x","4x"))

PS2.A

PS2.B<-read_csv("data/viability2.csv")%>%
  gather("tripl","value",5:7)%>%
  #mutate(conc=fct_reorder(conc))
  dplyr::filter(!is.na(value))%>%
  #mutate()%>%
  ggplot(aes(fct_inorder(conc),as.numeric(value)/210,col=stain))+
  geom_point(aes(shape=type),alpha=0.9,
             position=position_dodge(width = 0.3),
             size=3)+
  geom_line(stat="summary",fun.y="mean",
            aes(as.numeric(fct_inorder(conc)),as.numeric(value)/210,linetype=type),
            alpha=0.8,size=1,position=position_dodge(width = 0.3))+
  facet_grid(time~.)+
  ylab("% CFU")+xlab("")+
  scale_color_discrete(name="")+
  scale_shape_discrete(name="",labels=c("non-sporulating cells","purified spores"))+
  scale_linetype_discrete(name="",labels=c("non-sporulating cells","purified spores"))

PS2.B

plot_grid(PS2.A,PS2.B,nrow = 2,labels = c("A","B"))

ggsave("suppl/Supplemental1.pdf",width = 8, height = 10,dpi=300,units = "in")
ggsave("suppl/Supplemental1.png",width = 8, height = 10,dpi=300,units = "in")

FIG2.1<-distr.values%>%
  dplyr::filter(time==30&stain==2 |time==30&stain==0)%>%
  ggplot(aes(interaction(type,channel),diffM,fill=channel))+
  geom_bar(stat="identity",position="dodge",alpha=0.7)+
  geom_errorbar(aes(ymin=diffM-pooledSD,ymax=diffM+pooledSD),width=0.5)+xlab("")+
  ylab("Difference spore/cell distributions")+scale_color_discrete(name="")+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position = c(0.6,0.8),axis.text.x = element_text(angle=30,hjust=1))+
  scale_x_discrete(label=c("SYBR1","SYBR2","PI","FSC","SSC"))

FIG2.1

{
  f2.2<-df2%>%
    group_by(type,stain,time)%>%
    mutate(id = row_number())%>%
    tidyr::spread(channel,value,fill=NA)%>%
    ungroup()%>%
    dplyr::select(-id)
    
  #translate so that cutoff is at 0
  f2.2$asinh.FL1.A[f2.2$type=="SYBR1"]<-f2.2$asinh.FL1.A[f2.2$type=="SYBR1"]-cutoffs.list[5]
  f2.2$asinh.FL1.A[f2.2$type=="SYBR2"]<-f2.2$asinh.FL1.A[f2.2$type=="SYBR2"]-cutoffs.list[8]
  f2.2$asinh.FSC.A<-f2.2$asinh.FSC.A-cutoffs.list[10]
  f2.2$asinh.SSC.A<-f2.2$asinh.SSC.A-cutoffs.list[11]
  f2.2$asinh.FL3.A<-f2.2$asinh.FL3.A-cutoffs.list[2]
  
}

FIG2.2<-f2.2%>%
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

FIG2.2

plot_grid(FIG2.2,FIG2.1,rel_widths = c(0.60,0.4),labels = c("A","B"))
ggsave("fig/Figure_2.pdf",width = 10, height = 5)

#### Figure 3+4: Clustering ####
{
  rm(list=ls())
  source("src/functions.R")
  sample.var<-c("strain","time","tripl")
  fcsset3<-flowCreateFlowSet(filepath = "data/f3/set9/",sample_variables=sample.var,transformation = TRUE)
  
  df3<-fcsset3%>%
    Subset(norm2Filter("asinh.FSC.H","asinh.FSC.W",scale.factor = 2))%>%
    Subset(rectangleGate("asinh.FL1.A"=c(0,15),"asinh.FL3.A"=c(0,15)))%>%
    flowFcsToDf(.)
  
  #Load reference containing all subpopulations
  df3.ref<-ref.ds(strn="Bs02003",time="24h",tripl=2,df=df3) 
  set.seed(1)
  df3.mix<-mclust::Mclust(data = df3.ref,G = 3)
  
  #getting centers for visualization and export
  centers.list.df<-t(df3.mix$parameters$mean)
  write.csv(centers.list.df,"suppl/centers_f3.csv")
  center.locs<-factor(df3.mix$classification,levels=c(1,2,3))

# Prediction for other cases

  df3.list<-df3%>%dplyr::select(asinh.FSC.A,asinh.SSC.A,asinh.FL1.A,strain,time,tripl)%>%
    split(.,df3$strain)
  
  clusterB.predict<-c(
    predict.Mclust(df3.mix,newdata = df3.list[["Bs02003"]][,1:3])$classification,
    predict.Mclust(df3.mix,newdata = df3.list[["Bs02025"]][,1:3])$classification
  )
  
  clusterB.annot<-do.call(rbind,list(
    df3.list[["Bs02003"]],
    df3.list[["Bs02025"]])
  )
  
  clusterB.pred<-cbind(clusterB.annot,cluster=clusterB.predict)


#supplemental 4
clusterB.pred%>%
  ggplot(aes(asinh.SSC.A,asinh.FL1.A))+
  geom_hex(bins=200)+
  #geom_hex(aes(col=as.factor(cluster)),bins=50)+
  facet_grid(tripl~strain)+
  xlim(c(10,15))+ylim(c(2.5,14))+
  scale_fill_viridis(guide=FALSE)+
  theme_bw()

  ggsave("suppl/Supplemental_4.pdf",width = 4, height = 6);beep()
#ggsave("suppl/ps4_f4b.png",width = 4, height = 6)
}
# FIGURE 4B
{
c.count.B<-clusterB.pred%>%
  group_by(strain,time,cluster,tripl)%>%
  summarize(count=n())%>%
  ungroup()%>%
  group_by(strain,time,tripl)%>%
  mutate(perc.mean.count=100*count/sum(count))

  ccount2<-c.count.B%>%
  ggplot(aes(strain,perc.mean.count,
             col=factor(cluster,levels=c(2,1,3)),
             shape=factor(cluster,levels=c(2,1,3))))+
  geom_point(size=3,position=position_dodge(1),stat="identity")+
  #geom_dotplot(col="white",binaxis="y",position=position_dodge(),stackdir = "center")+
  # geom_errorbar(aes(ymin=perc.mean.count-perc.sd.count,ymax=perc.mean.count+perc.sd.count,
  #                   group=as.factor(cluster)),position=position_dodge(),alpha=0.6)+
  ylab("Proportion / %")+xlab("time / h")+
  #scale_y_continuous(expand = c(0,0))+
  scale_color_viridis(labels=c("Cells","Forespores","Spores"),discrete=TRUE,end = c(0.8),direction = -1,
                     alpha = 0.8)+
  scale_shape_discrete(labels=c("Cells","Forespores","Spores"))+
  theme_minimal()+
  theme(panel.spacing = unit(1, "lines"),legend.title=element_blank(),
        legend.position = c(0.20,0.9))+
  xlab("")

ccount2
}

#### Figure 4A

{
  source("src/functions.R")
  sample.var <- c("strain","time","stain","tripl")
  fcsset3A<-flowCreateFlowSet(filepath = "data/f3/set8/",
                          sample_variables=sample.var,transformation = TRUE)
    
  df3A<-fcsset3A%>%
    Subset(norm2Filter("asinh.FSC.H","asinh.FSC.W",scale.factor = 2))%>%
    Subset(rectangleGate("asinh.FL1.A"=c(0,15),"asinh.FL3.A"=c(0,15)))%>%
    flowFcsToDf()%>%
    dplyr::filter(stain!="unstained")

  ###splitting and applying gmm
  centers.list<-list(df3.mix$classification,
                     df3.mix$classification,
                     df3.mix$classification)
  
  df3A.list<-df3A%>%
    select(asinh.FSC.A,asinh.SSC.A,asinh.FL1.A,strain,stain,time)%>%
    split(.$strain)
  
  cluster.predict.A<-c(
    predict(df3.mix,newdata = df3A.list[["Bs02003"]][,1:3])$classification,
    predict(df3.mix,newdata = df3A.list[["Bs02018"]][,1:3])$classification,
    predict(df3.mix,newdata = df3A.list[["Bs02020"]][,1:3])$classification
  )
  
  cluster.annot.A<-do.call(rbind,list(
    df3A.list[["Bs02003"]],
    df3A.list[["Bs02018"]],
    df3A.list[["Bs02020"]])
  )
  
  cluster.pred.A<-cbind(cluster.annot.A,cluster=cluster.predict.A)
}

#Supplemental

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
    # group_by(strain,time,cluster)%>%
    # summarize(mean.count=mean(count),sd.count=sd(count))%>%
    # ungroup()%>%
    # group_by(strain,time)%>%
    # mutate(perc.mean.count=100*mean.count/sum(mean.count),perc.sd.count=100*sd.count/sum(mean.count))%>%
    separate(time,c("time","h"),2)%>%
    mutate(time=as.numeric(time))%>%
    group_by(strain,time,stain)%>%
    mutate(perc.count=100*count/sum(count))
  
  table(c.count.A$strain,c.count.A$time,c.count.A$cluster,c.count.A$stain)

  c.count.A%>%write.csv(file = "suppl/population_count.csv")
  
ccount1<-c.count.A%>%
    ggplot(aes(time,perc.count))+
     geom_point(aes(col=as.factor(cluster),shape=as.factor(cluster),group=time),
                size=2)+
     stat_summary(aes(col=as.factor(cluster)),
                fun.y = "mean",geom = "line",size=2)+
    #geom_dotplot(aes(x = as.factor(time),fill=as.factor(cluster)),binaxis="y",
                 # position=position_dodge(),stackdir = "center")+
    # stat_summary(aes(x = as.factor(time),y = perc.count,col=as.factor(cluster)),
    #              position=position_dodge(1),fun.y="mean",geom="point")+
    ylab("Proportion / %")+xlab("time / h")+
    scale_y_continuous(expand = c(0,0))+
    facet_grid(strain~.)+
    scale_color_viridis(labels=c("Cells","Forespores","Spores"),discrete=TRUE,end = c(0.8),direction = -1,
                       guide=FALSE,alpha = 0.8)+
    scale_shape_discrete(guide=FALSE)+
    theme_minimal()+
    theme(panel.spacing = unit(1, "lines"),legend.title=element_blank())+
    xlab("time / h")
  
  ccount1
  plot_grid(ccount1,ccount2,nrow = 1,rel_widths = c(0.6,0.4),labels = c("A","B"))

  ggsave("fig/Figure_4.pdf",width = 8, height = 5);beep()
  
}

#### Supplemental: Combinations ####

{
  remove(list=ls())
  source("src/functions.R")
  sample.var=c("ident","type","stain","conc","trash")
  fcs.df1<-flowCreateFlowSet(filepath = "data/reg/",sample_variables=sample.var)%>%
    Subset(.,norm2Filter("FSC-H", "FSC-W",filterId="norm_ssc.fsc",scale=1))%>%
    Subset(rectangleGate("asinh.FL1.A"=c(0,15),"asinh.FL3.A"=c(0,15)))%>%
    flowFcsToDf(.)


fcs.df1%>%
  select(type,stain,conc,asinh.FL1.A,asinh.FL3.A)%>%
  gather("channel","value",4:5)%>%
  ggplot(aes(value,interaction(stain,conc),fill=type))+
  geom_density_ridges(alpha=0.6)+ylab("")+
  scale_fill_discrete(name="",labels=c("non-sporulating cells","purified spores"))+
  theme(legend.position = c(0.45, 0.97))+
  facet_grid(.~channel)

ggsave("suppl/Supplemental_2.pdf",width = 10, height = 10)

fcs.df2<-fcs.df1%>%
  group_by(type,stain,conc)%>%
  summarize(mean.fl1=mean(asinh.FL1.A),mean.fl3=mean(asinh.FL3.A))%>%
  ungroup()%>%
  mutate(PI=   c(0,2,1,1,0,0,0,0,2,1,1,0,0,0),
         SYBR1=c(0,0,1,0,2,1,0,0,0,1,0,2,1,0),
         SYBR2=c(0,0,0,1,0,1,2,0,0,0,1,0,1,2))

fcs.df2%>%
  dplyr::filter(PI==0)%>%
  lm(as.numeric(mean.fl1)~-1+as.numeric(SYBR2)*as.numeric(SYBR1)+
       as.numeric(SYBR2),data=.)%>%
  summary()

fcs.df2%>%
  dplyr::filter(SYBR2==0)%>%
  lm(as.numeric(mean.fl3)~-1+as.numeric(SYBR1)*as.numeric(PI),data=.)%>%
  summary()

fcs.df2%>%
  dplyr::filter(SYBR1==0)%>%
  lm(as.numeric(mean.fl3)~-1+as.numeric(SYBR2)*as.numeric(PI),data=.)%>%
  summary()

write.csv2(fcs.df2,"suppl/mixdesign.csv")
}
