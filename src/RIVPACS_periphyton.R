require(adespatial)
library(tidyverse)
require (plyr)
require(dplyr)
require(ade4)
library(vegan)
library(analogue)
library(BiodiversityR)
library(MASS)
library(FactoMineR)

#####Combine data######
peri<-read.csv("Data/Periphyton_192021_070422.csv")

sampling<-read.csv("data/Sampling_sites_with_subbasin_ed.csv")

sampling$temp_ID <- paste(sampling$Station_Code, sampling$Year_sampled, sep="_") 

peri$temp_ID <- paste(peri$Station_Code, peri$Year, sep="_") 

Peri_OID<-merge(peri, sampling,  all.x=T, by="temp_ID")


geo<-read.csv("Data/All_var_comb_070422.csv")

comb<-merge(Peri_OID, geo, all.x=T, by="OBJECTID")



###Selelct natural variables#######

natural.var<-comb%>%dplyr::select(155:189)

natural.var.sel<-natural.var%>%dplyr::select(-Elev_Max, -Elev_Min, -Slope_Max, -Slope_Min, -MoraineD, -StMoraineD, 
                                      -ICTMoraineD, -ORGdepositD, -GlaciolacD, -LacustD, -BedrockD, -FluvialD,-EolianD, -ColluvialD )
#natural.var.sel<-natural.var%>%select( -Slope_Max, -Slope_Min)


###PCA

res.pca = PCA(natural.var.sel, scale.unit=T,ncp=4, graph=T)
evplot(res.pca$eig[,1])
summary(res.pca)
dimdesc(res.pca, axes = c(1: 3))

graph.var(res.pca, axes = c(1, 2))
plot(res.mfa,choix="ind", c(1,2))

graph.var(res.pca, axes = c(1, 3))
plot(res.mfa,choix="ind", c(1,3))

natural.PC<-cbind(res.pca$ind$coord[,1:3], comb["Ref_site"])

natural.PC.sel<-subset(natural.PC, comb$Ref_site =='Y')

##########

#Select taxa
comb.sel<-subset(comb, comb$Ref_site =='Y')

sp<-comb.sel[,5:94]


minocc=round(nrow(sp)*(5/100))
#zoop.sel <- chooseTaxa(cope, n.occ = minocc+1,  na.rm = TRUE)

sp.t<-decostand(sp, "total")

#zoop.x<- chooseTaxa(zoop.t, n.occ = minocc+1, max.abun = 0.01, type = "AND",  na.rm = TRUE)

#zoop.x<- chooseTaxa(zoop.t,  max.abun = 0.01, n.occ = minocc+1, type = "OR",  na.rm = TRUE)

sp.x<- chooseTaxa(sp.t,  max.abun = 0.01, n.occ = minocc+1, type = "AND",  na.rm = TRUE)

#zoop.x<- chooseTaxa(zoop.t,  max.abun = 0.01,  na.rm = TRUE)

#sp.names<-colnames(sp.x)
#sp.sel<-sp%>%select( all_of(sp.names))

sp.sel.pa<-decostand(sp.x, "pa")



#Cluster analysis

#D.jac = dist.ldc(sp,"jaccard")
D.sor = dist.ldc(sp.sel.pa,"sorensen")

D.sor.ed<-dist.zeroes(sp,D.sor)

D.sor.single<-hclust(D.sor.ed, method="single")


D.sor.complete<-hclust(D.sor.ed, method="complete")


D.sor.average<-hclust(D.sor.ed, method="average")


D.sor.ward<-hclust(D.sor.ed, method="ward.D")

D.sor.wardD2<-hclust(D.sor.ed, method="ward.D2")

D.sor.centroid <-hclust(D.sor.ed, method="centroid")

D.sor.median <-hclust(D.sor.ed, method="median")





D.sor.single.coph<-cophenetic(D.sor.single)
cor(D.sor.ed,  D.sor.single.coph)

D.sor.complete.coph<-cophenetic(D.sor.complete)
cor(D.sor.ed,  D.sor.complete.coph)

D.sor.average.coph<-cophenetic(D.sor.average)
cor(D.sor.ed,  D.sor.average.coph)

D.sor.ward.coph<-cophenetic(D.sor.ward)
cor(D.sor.ed,  D.sor.ward.coph)

D.sor.wardD2.coph<-cophenetic(D.sor.wardD2)
cor(D.sor.ed,  D.sor.wardD2.coph)

D.sor.centroid.coph<-cophenetic(D.sor.centroid)
cor(D.sor.ed,  D.sor.centroid.coph)

D.sor.median.coph<-cophenetic(D.sor.median)
cor(D.sor.ed,  D.sor.median.coph)

##Best clustering 
Best.clus<-D.sor.average

##Optimal number of clusters

grpdist <- function(X) {
  require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr, "gower")
  distgr
}


kt <- data.frame(k = 1:nrow(sp), r = 0)
for (i in 2:(nrow(sp) - 1)) {
  gr <- cutree(Best.clus, i)
  distgr <- grpdist(gr)
  mt <- cor(D.sor.ed, distgr, method = "spearman")
  #mt <- cor(D.jac, distgr, method = "spearman")
  kt[i, 2] <- mt
}
k.best <- which.max(kt$r)

kt.n<-kt[1:100,]
plot(
  kt.n$k,
  kt.n$r,
  type = "h",
  main = "Matrix correlation-optimal number of clusters",
  xlab = "k (number of clusters)",
  ylab = "Pearson's correlation"
)
axis(
  1,
  k.best,
  paste("optimum", k.best, sep = "\n"),
  col = "red",
  font = 2,
  col.axis = "red"
)
points(k.best,
       max(kt$r),
       pch = 16,
       col = "red",
       cex = 1.5)

asw<-numeric(nrow(sp))
for (k in 2:(nrow(sp)-1)) {
  sil<-silhouette(cutree(Best.clus, k=k), D.sor.ed)
  asw[k]<-summary(sil)$avg.width
}
k.best<-which.max(asw)
k.best

resize.win <- function(Width=6, Height=6){
  # works for windows
  dev.off() # dev.new(width=6, height=6)
  windows(record=TRUE, width=Width, height=Height)
}
resize.win(30,17)
par(mfrow = c(1,2), lwd=1.3, mar=c(3.6,3.4,.5,.5), col=1,  font=1, family="serif", cex=1.5, mgp=c(2.3, 1.2, 0))

plot(1:nrow(sp), asw, type="h", xlab="number of clusters", ylab="average silhouette width")
axis(1, k.best, paste("optimum",k.best,sep="\n"), col="lightgrey", font=2, family="serif",
     col.axis="lightgrey", mgp=c(1.9, 1.2, 0))
points(k.best, max(asw), pch=16, col="lightgrey", cex=1.5)

#Cluster summaries


cut_avg<-cutree(Best.clus, k=k.best)
#cut_avg<-cutree(Jac_average_hclust, k=14)

sp_avg<-cbind(cut_avg, sp)

cut_avg<-cutree(Best.clus, k=3)

######LDA####
#library(MASS)
library (Morpho)

linear <- lda(cut_avg~Dim.1+Dim.2+Dim.3, natural.PC.sel)

pred<-predict(linear, natural.PC)

fac<-pred$class

#fac<-as.factor(replace(fac, fac==3, 2))

sp.all<-comb[,5:94]

sp.all.t<-decostand(sp.all, "total")

sp.all.x<- chooseTaxa(sp.all.t,  max.abun = 0.05, n.occ = minocc+1, type = "AND",  na.rm = TRUE)

#zoop.x<- chooseTaxa(zoop.t,  max.abun = 0.01,  na.rm = TRUE)

#sp.names<-colnames(sp.x)
#sp.sel<-sp%>%select( all_of(sp.names))



sp.all.pa<-decostand(sp.all.x, "pa")



typClass2 <- typprobClass(sp.all.pa,groups=fac,method="c",cv=F)

Tprob<-typprobClass(sp.all.pa, groups = as.numeric(fac), small = T, method = "c",
                    outlier = 0.05, sep = F, cv = F)$probs




calculateOccurrenceProb <- function(occur, group){
  if (!is.matrix(occur)) {
    occur <- as.matrix(occur)
  }
  group.sum <- rowsum(occur, group)
  row.sum   <- as.numeric(table(group))
  group.sum / row.sum
}

calculateExpected <- function(x, occur, group) {
  prob <- calculateOccurrenceProb(occur, group)
  expect <- x %*% prob
  expect
}

Exp.sp<-calculateExpected(Tprob, sp.all.pa, fac)

#' Calculate the expected number of species from a null model.
#' 
#' Calculates the expected number of species based on a null model (all calibration data
#' is in 1 group).  These are just the probability of occurrence of each species
#' in the calibration data.
#' @inheritParams calculateExpected
#' @examples
#' x <- matrix(rbinom(100, size = 1, prob = 0.5), ncol = 10)
#' calculateNullExpected(x)
#' 
calculateNullExpected <- function(occur){
  colMeans(occur)
}

null.sp<-calculateNullExpected(sp.all.pa)

calculateRIVPACSMetrics <- function(O, E, N, cutoff = 0){
  N <- matrix(N, ncol = ncol(O), nrow = nrow(O), byrow = T)
  null.nonrare <- N >= cutoff
  On  <- O * null.nonrare
  En  <- N * null.nonrare
  Ons <- rowSums(On)
  Ens <- rowSums(En)
  BCn <- rowSums(abs(On - En)) / (Ons + Ens)
  
  nonrare <- E >= cutoff
  O  <- O * nonrare
  E  <- E * nonrare
  Os <- rowSums(O)
  Es <- rowSums(E)
  BC <- rowSums(abs(O - E)) / (Os + Es)
  
  cbind(O = Os, E = Es, OE = Os/Es, BC = BC, 
        Onull = Ons, Enull = Ens, OEnull = Ons / Ens, BCnull = BCn)
}

calculateRIVPACSMetrics(sp.all.pa, Exp.sp, null.sp, cutoff = 0)



###LCBD#####
nonrare <- sp.E >= 0.5

Osp  <- sp.all.pa * nonrare

Esp  <- Exp.sp * nonrare

#O.LCBD =  beta.div(Osp, "sorensen", nperm=999, sqrt.D = FALSE)$LCBD 
#E.LCBD =  beta.div(Esp, "sorensen", nperm=999, sqrt.D = FALSE)$LCBD 

#Observed

Ob.D.sor = dist.ldc(Osp,"sorensen")

Ob.D.sor<-dist.zeroes(Osp, Ob.D.sor)

O.LCBD = LCBD.comp(Ob.D.sor, sqrt.D = T, save.D = FALSE)$LCBD

#Expected

Ex.D.sor = dist.ldc(Esp,"sorensen")

Ex.D.sor<-dist.zeroes(Esp, Ex.D.sor)

E.LCBD = LCBD.comp(Ex.D.sor, sqrt.D = T, save.D = FALSE)$LCBD


LCBD.ratio = O.LCBD/E.LCBD


##Natural predictors

chem[,1:7]<-log10(chem[,1:7])
chem[,9:20]<-log10(chem[,9:20])


#chem.scale<-as.data.frame(apply(chem, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X))))

rda.chem.full<- rda(zoop.h ~.,   data=chem)
anova.cca(rda.chem.full, step=1000) #check if model is significant, proceed if only significant
(R2a.all<-RsquareAdj(rda.chem.full)$adj.r.squared)
(R2a.all<-RsquareAdj(rda.chem.full)$r.squared)
vif.cca(rda.chem.full) #Check multicolinearity

forward.chem<-forward.sel(zoop.h,chem,adjR2thresh = 0.17, nperm=999, alpha = 0.1)

#sel.env<-forward.sel(res.D$vectors,env.n, nperm=999, alpha = 0.05)
forward.chem

library(corrplot)
res<-cor(chem)
res1 <- cor.mtest(chem, conf.level = .95)


library(corrplot)
res<-cor(chem)
res1 <- cor.mtest(chem, conf.level = .95)

corrplot(res, type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45)

#chem.sel<-subset(chem, select = -c(Temp_C, Cond, TDS, TIC, SO4,  TN, TOC, Alkalinity))
#chem.sel<-subset(chem, select = -c(Cond, TDS, TIC, Alkalinity, TN, TOC))
chem.sel<-subset(chem, select = -c(Cond, TDS, TIC ))
rda.chem.full<- rda(zoop.pa ~.,   data=chem)

chem.sel<-chem[1:20,]
x<-predict(rda.chem.full, type="response", newdata=chem.sel, scaling="none")
x


##Predict
rda.db<-capscale(zoop.pa ~ .,chem, add =TRUE)
x<-predict(rda.db, type="response", newdata=chem.sel, scaling="none")
x