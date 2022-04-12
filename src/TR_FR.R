
TR_FR<-function(trait,abun)
{
  # library required for indices computation
  require (geometry)
  require(ape)
  
  ##############################################################################
  # names of assemblages
  nm_asb<-row.names(abun)
 
  # matrix to store results
  indices<-c( "TR", "FR")
  TFR<-matrix(NA, length(nm_asb), length(indices), dimnames=list(nm_asb,indices))
  
  # convex hull volume of the species pool
  FRic_pool<-convhulln(trait,"FA")$vol
  ##############################################################################

  # loop on assemblages for computing taxonomic and functional richness (TR,FR)
  for (k in nm_asb)
  {
    ###########################################################
    # names, number, abundance and trait of of species present
    abun_k<-abun[k,]
    nm_sp_k<-row.names(trait)[which(abun_k>0)]
    nb_sp_k<-length(nm_sp_k)
    abun_sp_k<-abun[k,nm_sp_k]
    trait_sp_k<-trait[nm_sp_k,]
    if(nb_sp_k==1) { trait_sp_k<-matrix(trait_sp_k,nrow=1,dimnames=list(nm_sp_k,colnames(trait)) ) } # matrix object
    
    # Taxonomic richness
    TFR[k,"TR"]<-nb_sp_k
   
    # Functional richness = convex hull volume
    TFR[k,"FR"]<-round(convhulln(trait_sp_k,"FA")$vol/FRic_pool,6)
   }

  return(TFR)  
}