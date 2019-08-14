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
      singlet_gate <- gate_singlet(fcsset[[i]], area = "FSC-A", height = "FSC-H", filterId = "Singlets",wider_gate = TRUE )
      
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


