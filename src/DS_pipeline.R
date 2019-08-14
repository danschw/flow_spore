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

for (i in 1:n.sample){
      singlet_gate <- gate_singlet(fcsset[[i]], area = "FSC-A", height = "FSC-H", filterId = "Singlets",wider_gate = FALSE )

      #plot gate
      p.singlet <-
            ggcyto(fcsset[[i]], aes(x = `FSC-A`, y =  `FSC-H`))+
            geom_hex(bins = 500)+
            geom_gate(singlet_gate)+
            geom_stats()+
            ggtitle("Singlets")
      #save plot
      id <- fcsset[[i]]@description$GUID
      ggsave(paste0("fig/DS_figures/Singlet_gates/",id,".pdf" ), p.singlet)

      #apply gate
      fcsset[[i]] <- fcsset[[i]] %>%
            Subset(singlet_gate)%>%
            Subset(rectangleGate("asinh.BL1.A" = c(0, 15), "asinh.FSC.A" = c(0, 15))) # remove negatives

}

#### Gating for singlets with norm2Filter ####
# df.set <- fcsset %>%
#    Subset(norm2Filter("asinh.FSC.H", "asinh.FSC.A", scale.factor = 10)) %>%
#    Subset(rectangleGate("asinh.BL1.A" = c(0, 15))) %>%
#    flowFcsToDf(.)











#### transform to dataframe ####
df.set <- fcsset%>%
      flowFcsToDf(.)

df.set %>%
   group_by(well) %>%
   summarise(events=n())

# # plot
df.set %>%
      ggplot(aes(asinh.FSC.A, asinh.BL1.A))+
      geom_hex(bins=50)#+
#       facet_wrap(~well)

#### predict centers of sub-populations ####

# Load reference containing all subpopulations
# I will use no phage triplicate



df.mix <- df.set%>%
      # dplyr::filter(phage=="noPHI")%>%
      dplyr::select(asinh.FSC.A,asinh.SSC.A,asinh.BL1.A)%>%
      as.matrix() %>%
      mclust::Mclust(data = ., G = 2) #I  have only 2 clusters

# getting centers for visualization and export
centers.list.df <- t(df.mix$parameters$mean)
center.locs <- factor(df.mix$classification, levels = c(1, 2))

### Prediction of clusters for all samples ####
cluster.predict <- df.set %>%
      dplyr::select(asinh.FSC.A, asinh.SSC.A, asinh.BL1.A) %>%
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
   data.frame(df.set, cluster = cluster.predict$classification) %>%
   ggplot(aes(asinh.SSC.A, asinh.BL1.A)) +
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
   geom_point(aes(centers.list.df[3, 1], centers.list.df[3, 2]), col = "blue", size = 1) 
# +
#    facet_wrap(~well)
ggsave("fig/DS_figures/first/R1_3clust.png", plot = p)


#### Get the quantities ####
clust.count <- data.frame(df.set, cluster = cluster.predict$classification) %>%
   group_by(well, cluster) %>%
   summarize(count = n()) %>%
   mutate(perc.mean.count = 100 * count / sum(count), total=sum(count))
