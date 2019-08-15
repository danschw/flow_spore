#analysis of SP-10 entrapment experiment using Karava's analysis.R as reference.
rm(list=ls())
source("src/functions.R")
source("src/DS_functions.R")

#### Load data, sample set ####
set.seed(1)
sample.var <- c("host","media","time","dilution","well","phage","rep")
fcsset <- flowCreateFlowSet(filepath = "data/entrap_data/day1/delta6_DSM_T4_x10/", sample_variables = sample.var, transformation = FALSE)
#transform with arcsine, recpmendded by Karava et al.
fcsset <- Transform.Novocyte(fcsset)

#### Gating for singlets with flowStats ####

# The gate function needs to be applied to each sample seperatly
# get number of samples
n.sample <- nrow(fcsset@phenoData@data)

# initialise daaframe to store noise cutoff values
noise.fsc <- data.frame(sample=fcsset@phenoData@data$name,well=fcsset@phenoData@data$well,cutoff=rep(NA,n.sample))

#initialise list to store singlet plots
plot.list <- vector('list', n.sample)

for (i in 1:n.sample){
      singlet_gate <- gate_singlet(fcsset[[i]], area = "FSC-A", height = "FSC-H", filterId = "Singlets",wider_gate = TRUE )

      #plot gate
      id <- fcsset[[i]]@description$GUID
      plot.list[[i]] <-
         as.ggplot(
            ggcyto(fcsset[[i]], aes(x = `FSC-A`, y =  `FSC-H`))+
               geom_hex(bins = 500)+
               geom_gate(singlet_gate, size=0.01)+
               geom_stats()+
               theme_cowplot(font_size = 5)+
               scale_y_continuous(labels =  function(x) format(x, scientific = TRUE))+
               scale_x_continuous(labels =  function(x) format(x, scientific = TRUE))+
               facet_null()+
               ggtitle(fcsset[[i]]@description$`$WELLID`)
            )

   
      #apply gate
      fcsset[[i]] <- fcsset[[i]] %>%
            Subset(singlet_gate)%>%
            Subset(rectangleGate("asinh.BL1.A" = c(0, 15), "asinh.FSC.A" = c(0, 15),"asinh.SSC.A" = c(0, 15))) # remove negatives

      #remove noise by FSC
      noise<- rangeGate(x = fcsset[[i]], "asinh.FSC.A", absolute = F)
      noise.fsc$cutoff[i] <- noise@min[[1]]

}

#save plot
ggsave2(paste0("fig/DS_figures/Singlet_gates/",fcsset[[i]]@description$`$SRC`,".pdf" ), plot_grid(plotlist = plot.list))


#### Gating for singlets with norm2Filter ####
# df.set <- fcsset %>%
#    Subset(norm2Filter("asinh.FSC.H", "asinh.FSC.A", scale.factor = 10)) %>%
#    Subset(rectangleGate("asinh.BL1.A" = c(0, 15))) %>%
#    flowFcsToDf(.)

# plotting noise gate

noise.plot <-   
ggcyto(fcsset,aes(asinh.FSC.A))+
   geom_density()+
   geom_vline(data = noise.fsc, aes(xintercept=cutoff), color="red")+
   facet_wrap(~well)+theme_cowplot()

ggsave2("fig/DS_figures/gates/noise.PDF", noise.plot)


#### transform to dataframe ####
df.set <- fcsset%>%
      flowFcsToDf(.)

df.set %>%
   group_by(well) %>%
   summarise(events=n())

# # plot
df.set %>%
      ggplot(aes(asinh.SSC.A, asinh.BL1.A))+
      geom_hex(bins=500)+
      facet_wrap(~well)

clean.df <- df.set[0,]

# remove noise
wells <- levels(as.factor(df.set$well))
for (wl in wells){
   clean.df <- rbind(clean.df,
   df.set %>%
      filter(well==wl)%>%
      filter(asinh.FSC.A>noise.fsc$cutoff[noise.fsc$well==wl]))
}

clean.df%>%
   ggplot(aes(asinh.SSC.A, asinh.BL1.A))+
   geom_hex(bins=500)+
   facet_wrap(~well)

clean.df %>%
   group_by(well) %>%
   summarise(events=n(), perc=100*n()/5e4)

complot <- 
   ggplot(df.set, aes(asinh.SSC.A, asinh.BL1.A))+
   geom_point(size=0.1, color="grey")+
   geom_point(data=clean.df,size=0.1, color="blue")+
   facet_wrap(~well)

ggsave2("fig/DS_figures/gates/com_noise.PNG", complot)
#### predict centers of sub-populations ####

# Load reference containing all subpopulations
# I will use no phage triplicate

df.mix <- clean.df%>%
      dplyr::filter(phage=="noPHI")%>%
      dplyr::select(asinh.FSC.A,asinh.BL1.A)%>%
      as.matrix() %>%
      mclust::Mclust(data = ., G = 2) #I  have only 2 clusters




# getting centers for visualization and export
centers.list.df <- t(df.mix$parameters$mean)
center.locs <- factor(df.mix$classification, levels = c(1, 2))

### Prediction of clusters for all samples ####
cluster.predict <- clean.df %>%
      dplyr::select( asinh.FSC.A, asinh.BL1.A) %>%
      predict.Mclust(df.mix,.)



# p <-       
#    data.frame(df.set, cluster = cluster.predict$classification) %>%
#    ggplot(aes(asinh.SSC.A, asinh.BL1.A)) +
#    geom_hex(aes(fill = factor(cluster,levels=c(2,1))), bins = 300) + # ,alpha=..ncount.. #order= ?
#    geom_density2d(col = "red", bins = 20, size = 0.5, alpha = 0.7) +
#    xlim(c(5, 15)) + ylim(c(2.5, 15)) +
#    scale_fill_viridis(
#       discrete = TRUE, end = 0.8, label = c("Cells", "Spores"), name = "", direction = -1,
#       guide = FALSE
#    ) +
#    scale_alpha_continuous(guide = FALSE) +
#    theme_bw() +
#    geom_point(aes(centers.list.df[1, 2], centers.list.df[1, 3]), col = "blue", size = 1) +
#    geom_point(aes(centers.list.df[2, 2], centers.list.df[2, 3]), col = "blue", size = 1) +
#    facet_wrap(~well)
# ggsave("fig/DS_figures/first/plot_statsGate.png", plot = p)


p <-       
   data.frame(clean.df, cluster = cluster.predict$classification) %>%
   ggplot(aes(asinh.FSC.A, asinh.BL1.A)) +
   geom_hex(aes(fill = factor(cluster,levels=c(2,1))), bins = 300) + # ,alpha=..ncount.. #order= ?
   # geom_density2d(col = "red", bins = 20, size = 0.5, alpha = 0.7) +
   # xlim(c(5, 15)) + ylim(c(2.5, 15)) +
   # scale_fill_viridis(
   #    discrete = TRUE, end = 0.8, label = c("Cells", "Spores"), name = "", direction = -1,
   #    guide = FALSE
   # ) +
   # scale_alpha_continuous(guide = FALSE) +
   theme_bw()+ 
   geom_point(aes(centers.list.df[1, 1], centers.list.df[1, 2]), col = "blue", size = 1) +
   geom_point(aes(centers.list.df[2, 1], centers.list.df[2, 2]), col = "blue", size = 1) +
   facet_wrap(~well)
ggsave("fig/DS_figures/gates/clust.png", plot = p)


#### Get the quantities ####
clust.count <- data.frame(clean.df, cluster = cluster.predict$classification) %>%
   group_by(well, cluster) %>%
   summarize(count = n()) %>%
   mutate(perc.mean.count = 100 * count / sum(count), total=sum(count))
