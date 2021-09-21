library(tidyverse)
library(vegan)
library(FactoMineR)

all<-read.csv("Data/Geospatial_var_NSRB_200921.csv")
all=na.omit(all)
