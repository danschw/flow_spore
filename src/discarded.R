#####f1####

#
fcsset1<-createFlowset(filepath = "data/Mka_Spores_Separation_3/Setting 1/",
                       sample.variables=c("ident","type","stain","mode","run"))
fcsset1

fcsset2<-createFlowset(filepath = "data/MKA_Spores_Separation_4/Combination stains - 1/",
                       sample.variables=c("ident","type","stain","run"))
fcsset2

fcsset3<-createFlowset(filepath = "data/MKA_Spores_Separation_4/1fold_2fold_stain_30/",
                       sample.variables=c("ident","type","stain","run"))
fcsset3

#Set1:S3 effect of afterstaining!
fcsset1%>%autoplot("asinh.FSC.A")
fcsset1%>%autoplot("asinh.SSC.A")
fcsset1%>%autoplot("asinh.FSC.A","asinh.SSC.A")

fcsset1%>%
  Subset(.,norm2Filter("asinh.FSC.H", "asinh.FSC.W",filterId="norm_ssc.fsc",scale=1))%>%
  Subset(.,norm2Filter("asinh.SSC.H", "asinh.SSC.W",filterId="norm_ssc.fsc",scale=1))%>%
  autoplot("asinh.FSC.A")

#Set2:S4 combinations
#highly fragmented by staining procedure, not usable for scattering!
fcsset2%>%autoplot("asinh.FL1.A")
fcsset2%>%autoplot("asinh.FL3.A")

#Set3:S4 different concentrations
fcsset3%>%autoplot("asinh.FSC.A")
fcsset3%>%autoplot("asinh.SSC.A")


#AAAA



####f2####

####f3####

####f4####