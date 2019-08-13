#analysis of SP-10 entrap,ent experiment using Karava's analysis.R as reference.
rm(list=ls())
source("src/functions.R")
source("src/DS_functions.R")

###########################################
# PLUGGING MY DATA INOT FIGURE 2 ANALYSIS #
###########################################
set.seed(1)
# my understanding is that the "flowCreateFlowSet" function parse the file name at underscores, dashes, periods etc.
# In fact found this definition separators="[-_\\.]" 
# e.g. file name delta6_LB_T0_x100_E6_SP10dF_C.fcs
# host_media_time_dilution_well_phage_rep
sample.var <- c("host","media","time","dilution","well","phage","rep")
fcsset.1 <- flowCreateFlowSet(
      filepath = "data/entrap_data/day1/delta6_DSM_T4_x10/",
      sample_variables = sample.var, transformation = FALSE
)

fcsset.1 <- Transform.Novocyte(fcsset.1)


full.df <- fcsset.1 %>%
      Subset(rectangleGate(
            "asinh.FSC.A" = c(0, 15),
            "asinh.SSC.A" = c(0, 15),
            "asinh.BL1.A" = c(0, 15)
      )) 

subset.df <- full.df %>%
      Subset(. , norm2Filter("FSC-H", "FSC-A", filterId = "norm_ssc.fsc", scale = 2)) %>%
      flowFcsToDf(.) %>%
      # dplyr::filter(ident == "Cells+spores") %>%
      select(host,media,time,dilution,well,phage,rep, asinh.BL1.A, asinh.FSC.A, asinh.SSC.A) %>%
      gather("channel", "value", 8:10)

models.list <- subset.df %>%
      split(x = ., interaction(.$host, .$media, .$time, .$dilution, .$well, .$phage, .$rep, .$channel)) %>%
      Filter(function(x) nrow(x) != 0, .) %>%
      lapply(., function(x) {
            x %>%
                  dplyr::filter(value > 5) %>%
                  with(value) %>%
                  Mclust(., 2)
      })

# get distribution parameters
distr.values <- lapply(models.list, function(x) {
      c(x$parameters$mean, sqrt(firstsecondPair(x$parameters$variance$sigmasq)), x$parameters$pro)
}) %>%
      do.call("rbind", .) %>%
      as_tibble() %>%
      magrittr::set_colnames(c("mu_1", "mu_2", "sd_1", "sd_2", "pro_1", "pro_2")) %>%
      mutate(
            diffM = mu_2 - mu_1, # difference between means
            pooledSD = sqrt((sd_1^2 + sd_2^2) / 2), # pooled standard deviation
            # cutoff = findCutoffs(mu_1, sd_1, pro_1, mu_2, sd_2, pro_2) # cutoff values
      )

# %>%
#       cbind(., cases) %>%
#       separate(cases, into = c("type", "stain", "time", "z", "channel"))
#############################################
# PLUGGING MY DATA INOT FIGURE 3,4 ANALYSIS #
#############################################
set.seed(1)
sample.var <- c("host","media","time","dilution","well","phage","rep")
fcsset3 <- flowCreateFlowSet(filepath = "data/entrap_data/day1/delta6_DSM_T4_x10/", sample_variables = sample.var, transformation = FALSE)
fcsset3 <- Transform.Novocyte(fcsset3)

sub.fcsset3 <- fcsset3 %>%
      Subset(norm2Filter("asinh.FSC.H", "asinh.FSC.A", scale.factor = 6)) %>%
      Subset(rectangleGate("asinh.BL1.A" = c(0, 15), "asinh.FSC.A" = c(0, 15)))

df3 <- sub.fcsset3%>%
      flowFcsToDf(.)

full.df3 <- fcsset3 %>%
  flowFcsToDf(.)

########
#looking at the subsetting for singlets
  # ggplot(df3, aes(asinh.SSC.A, asinh.BL1.A))+
  ggplot(df3, aes(FSC.H, FSC.A))+
  geom_density2d(data=full.df3, color="red")+
    geom_density2d(color="blue")+
    facet_wrap(~well)

pt <- 
  ggplot(df3, aes(FSC.H, FSC.A))+
  geom_point(data=full.df3, color=rgb(1,0,0,0.5))+
  geom_point(color=rgb(0,1,0,0.5))+
  facet_wrap(~well)
ggsave("fig/DS_figures/first/outlierTest_sf6.png",plot = pt)

df3 %>%
  ggplot(aes(asinh.SSC.A, asinh.BL1.A))+
  geom_density2d()+
  facet_wrap(~well)
######

# Load reference containing all subpopulations
df3.ref <- df3%>%
      dplyr::filter(well=="C7")%>%
      dplyr::select(asinh.FSC.A,asinh.SSC.A,asinh.BL1.A)%>%
      as.matrix()
df3.mix <- mclust::Mclust(data = df3.ref, G = 2) #I mostly have only 2 clusters


# getting centers for visualization and export
centers.list.df <- t(df3.mix$parameters$mean)
# write.csv(centers.list.df, "suppl/centers_f3.csv")
center.locs <- factor(df3.mix$classification, levels = c(1, 2))


# Prediction for other cases
cluster.predict <- df3 %>%
      dplyr::select(asinh.FSC.A, asinh.SSC.A, asinh.BL1.A) %>%
      predict.Mclust(df3.mix,.)


p.ref <- 
    data.frame(df3.ref, cluster = center.locs) %>%
      ggplot(aes(asinh.SSC.A, asinh.BL1.A)) +
      geom_hex(aes(fill = factor(cluster,levels=c(2,1,3))), bins = 300) + # ,alpha=..ncount.. #order= ?
      geom_density2d(col = "red", bins = 20, size = 0.5, alpha = 0.7) +
      xlim(c(5, 15)) + ylim(c(2.5, 15)) +
      scale_fill_viridis(
            discrete = TRUE, end = 0.8, label = c("Cells", "Forespores", "Spores"), name = "", direction = -1,
            guide = FALSE
      ) +
      scale_alpha_continuous(guide = FALSE) +
      theme_bw() +
      geom_point(aes(centers.list.df[1, 2], centers.list.df[1, 3]), col = "blue", size = 1) +
      geom_point(aes(centers.list.df[2, 2], centers.list.df[2, 3]), col = "blue", size = 1) 
p <-       
  data.frame(df3, cluster = cluster.predict$classification) %>%
      ggplot(aes(asinh.SSC.A, asinh.BL1.A)) +
      geom_hex(aes(fill = factor(cluster,levels=c(2,1,3))), bins = 300) + # ,alpha=..ncount.. #order= ?
      geom_density2d(col = "red", bins = 20, size = 0.5, alpha = 0.7) +
      xlim(c(5, 15)) + ylim(c(2.5, 15)) +
      scale_fill_viridis(
            discrete = TRUE, end = 0.8, label = c("Cells", "Forespores", "Spores"), name = "", direction = -1,
            guide = FALSE
      ) +
      scale_alpha_continuous(guide = FALSE) +
      theme_bw() +
      geom_point(aes(centers.list.df[1, 2], centers.list.df[1, 3]), col = "blue", size = 1) +
      geom_point(aes(centers.list.df[2, 2], centers.list.df[2, 3]), col = "blue", size = 1) +
      facet_wrap(~well)
ggsave("fig/DS_figures/first/plot_sf6.png", plot = p)

#get the quantities
c.count.B <- data.frame(df3, cluster = cluster.predict$classification) %>%
      group_by(well, cluster) %>%
      summarize(count = n()) %>%
      mutate(perc.mean.count = 100 * count / sum(count), total=sum(count))


