# Additional Packages
if(!require(flowStats)){
      BiocManager::install("flowStats")  
}
library(flowStats) 


# additional functions for adapting analysis for my data

Transform.Novocyte <- function(fcsset){
      flowCore::transform(fcsset,
                          `asinh.BL1-A`=asinh(`BL1-A`),
                          `asinh.BL1-H`=asinh(`BL1-H`),
                          `asinh.FSC-A`=asinh(`FSC-A`),
                          `asinh.FSC-H`=asinh(`FSC-H`),
                          `asinh.SSC-A`=asinh(`SSC-A`),
                          `asinh.SSC-H`=asinh(`SSC-H`)
      )
}