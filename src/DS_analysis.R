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


f1.df <- fcsset.1 %>%
      Subset(rectangleGate(
            "asinh.FSC.A" = c(0, 15),
            "asinh.SSC.A" = c(0, 15),
            "asinh.BL1.A" = c(0, 15)
      )) %>%
      Subset(., norm2Filter("FSC-H", "FSC-A", filterId = "norm_ssc.fsc", scale = 2)) %>%
      flowFcsToDf(.) %>%
      # dplyr::filter(ident == "Cells+spores") %>%
      select(host,media,time,dilution,well,phage,rep, asinh.BL1.A, asinh.FSC.A, asinh.SSC.A) %>%
      gather("channel", "value", 8:10)

models.list <- f1.df %>%
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
###########################################
# PLUGGING MY DATA INOT FIGURE 2 ANALYSIS #
###########################################