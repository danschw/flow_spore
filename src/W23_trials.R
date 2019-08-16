ggcyto(fcsset, aes(asinh.SSC.A, asinh.BL1.A))+
      geom_hex(bins=300)

flowCore::transform(fcsset,`log.SSC-A`=logTransform(`SSC-A`))

logTrans <- biexponentialTransform()
trans <- transformList('SSC-A', logTrans)
dataTransform <- transform(fcsset, trans)

ggcyto(fcsset, aes(`SSC-A`, asinh.BL1.A))+
      geom_hex(bins=300)+
      scale_x_flowJo_biexp()
library(flowMeans)

res <- flowMeans(dplyr::filter(clean.df, well=="B1"), varNames=c("asinh.SSC.A","asinh.FSC.A","asinh.BL1.A"), NumC = 2, MaxN = 2)
# plot(clean.df[,c(14,9)], col=res@Labels[[1]], pch=20)

tmp <- dplyr::filter(clean.df, well=="B1")
tmp$clust <- factor(res@Labels[[1]])

tmp%>%
ggplot(aes(asinh.FSC.A,asinh.BL1.A))+
      geom_point(aes(color=clust), size=0.01)

library(flowClust)
res <- flowClust(dplyr::filter(clean.df, well=="B1"), varNames=c("asinh.BL1.A", "asinh.FSC.A"), K = 2)

tmp <- dplyr::filter(clean.df, well=="B1")
tmp$clust <- factor(res@label)

tmp%>%
      ggplot(aes(asinh.SSC.A,asinh.BL1.A))+
      geom_point(aes(color=clust), size=0.01)




res <- clean.df%>%
      # dplyr::filter( well=="B10")%>%
      dplyr::select(asinh.FSC.A,asinh.BL1.A)%>%
      as.matrix() %>%
      mclust::Mclust(data = ., G = 2) #I  have only 2 clusters


# tmp <- dplyr::filter(clean.df, well=="B10")
tmp <- clean.df
tmp$clust <- factor(res$classification)

tmp%>%
      ggplot(aes(asinh.SSC.A,asinh.BL1.A))+
      geom_point(aes(color=clust), size=0.01)+
      geom_density_2d(color="white", size=0.1)+
      facet_wrap(~well)

res <- list()
tmp <- clean.df
tmp$clust <- NA
for(i in 1:10){
      # res[[i]] <- clean.df%>%
      #       dplyr::filter( well==df.stats$well[i])%>%
      #       dplyr::select(asinh.FSC.A,asinh.BL1.A)%>%
      #       as.matrix() %>%
      #       mclust::Mclust(data = ., G = 2) #I  have only 2 clusters 
      tmp$clust[tmp$well==df.stats$well[i]] <- factor(res[[i]]$classification)
}

tmp%>%
      ggplot(aes(asinh.SSC.A,asinh.BL1.A))+
      geom_point(aes(color=clust), size=0.01)+
      geom_density_2d(color="white", size=0.1)+
      facet_wrap(~well)

   