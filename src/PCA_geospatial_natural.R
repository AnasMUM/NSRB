library(tidyverse)
library(vegan)
library(FactoMineR)

all<-read.csv("Data/All_var_comb_ed_220921.csv")
#all=na.omit(all)

sel.var<-distinct(all, OBJECTID, .keep_all = TRUE)
natural.var<-sel.var%>%select(12:14, 43:78)

natural.var<-natural.var%>%select(-TransD, -Cluster, -Natural.Region)


natural.var[,1]<-factor(natural.var[,1], levels =c("1", "2", "3", "4", "5")) #Natural region
natural.var[,2]<-factor(natural.var[,2], levels =c("1", "2", "3")) #Cluster
natural.var[,3]<-factor(natural.var[,3], levels =c("Y", "N")) # Sampling status

res.pca = PCA(natural.var[,3:ncol(natural.var)], scale.unit=T, ncp=4, graph=T)
summary(res.pca)
dimdesc(res.pca, axes = c(1: 4))

graph.var(res.pca, axes = c(1, 2),xlim = NULL, ylim = NULL, col.sup = "blue",col.var = "black", lim.cos2.var = 0.01,cex = 1, title = NULL)


##RDA based PCA----
natural.pca<- rda (natural.var[,3:ncol(natural.var)] , scale=T)
plot(natural.pca, display=c(  "species", "cn"))
plot(natural.pca, display=c(  "sites", "cn"))


#Multiple Factor Analysis

natural.var<-natural.var%>%select(-Cluster, -Natural.Region)

res.mfa = MFA(natural.var, group=c(1, 1, 3, 2, 4, 5, 6, 10, 5), type=c("n", rep("s",8)), ncp=4, 
              name.group=c( "categories", "protected_areas", "climate", "hydrologic", "vegetation", "wetlands", "topography","surficial_geology","riparian"),
              num.group.sup=c(1))

#Source Evplot script

evplot(res.mfa$eig[,1])

dimdesc(res.mfa, axes = c(1: 4))
summary(res.mfa)

plot(res.mfa,choix="ind",partial="all")

res.hcpc <- HCPC(res.mfa, nb.clust=0, conso=0, min=3, max=10)

res.hcpc$call$t$inert.gain

res.hcpc$data.clust

res.hcpc$desc.var

res.hcpc$desc.axes
