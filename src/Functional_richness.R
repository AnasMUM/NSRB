


species_trait<-read.table("Data/Periphyton_traits_ed.csv", row.names = 1, header=TRUE, sep=",")



##Select observed

common <- intersect(rownames(species_trait), colnames(Osp))

Trait_sel <- subset(species_trait, rownames(species_trait) %in% common) 

sp.names<-rownames(species_trait)

Osp_sel<-Osp%>%dplyr::select( all_of(common))


###Select expected

common <- intersect(rownames(species_trait), colnames(Esp))

Trait_sel <- subset(species_trait, rownames(species_trait) %in% common) 

sp.names<-rownames(species_trait)

Esp<-as.data.frame(Esp)

Esp_sel<-Esp%>%dplyr::select( all_of(common))

#Ob_species_sel<-Osp%>%dplyr::select( all_of(sp.names))

### Summarize species traits in PCoA on Gower's distance to compute FR and FD
# removing scientific name and family from trait matrix
temp<-species_traits_his %>% select (-NEW.SCINAME,-FAMILY)
traitgow.d<-gowdis(Trait_sel)
trait_pcoa<-pcoa(traitgow.d, correction="cailliez")

# extracting scores from the top 7 significant (broken-stick) and large PCoA axes
trait_pcoa_scores<-as.data.frame(cbind(rownames(Trait_sel),trait_pcoa$vectors.cor[,1:7]))

# a quick diagnostic PCoA ordination plot using just the first 2 PCoA Axes
ggplot(trait_pcoa_scores, aes(Axis.1, Axis.2)) + geom_point(aes(colour=V1),size=5,shape=16,alpha=.4) + theme(legend.position="bottom") 

# selecting only those watersheds (in both excepted and observed) with at least 3 species (needed to calculate accurate convex hulls for FR)
# note that increasing the number of PCoA trait axes used, necessarily increases the minumum number of species required in a watershed
Osp_red<- Osp_sel %>% filter(rowSums(Osp_sel)>=3)
Esp_red<- Esp_sel %>% filter(rowSums(Esp_sel)>=3)

# selecting those watersheds (n=2828) in common so that ultimately hist-current differences can be calculated
common <- intersect(rownames(Osp_red), rownames(Esp_red))
Osp_red2 <- subset(Osp_red, rownames(Osp_red) %in% common) 
Esp_red2 <- subset(Esp_red, rownames(Esp_red) %in% common) 

# FUNCTIONAL RICHNESS
# calculating taxonomic and functional richness (using just PCoA 1 and 2)
source("src/TR_FR.R")
a<-TR_FR(trait_pcoa_scores[,2:3],Osp_red2) 
b<-TR_FR(trait_pcoa_scores[,2:3],Esp_red2) 
Peri_TRandFR<-cbind(a,b)

colnames(Peri_TRandFR)<-c("TR_Obs","FR_Obs","TR_Exp","FR_Exp")
Peri_TRandFR<-as.data.frame(Peri_TRandFR)
Peri_TRandFR$OE_f<-Peri_TRandFR$FR_Obs/Peri_TRandFR$FR_Exp
Peri_TRandFR
