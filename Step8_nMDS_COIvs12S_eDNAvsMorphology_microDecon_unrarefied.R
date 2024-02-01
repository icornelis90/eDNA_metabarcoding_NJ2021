#!/usr/bin/env Rscript
if(.Platform$OS.type == "unix"){home<-"/home/"} else{
  home<-"//192.168.236.131/"}

library(vegan)
library(seqRFLP)
library(dplyr)
library(tibble)
library(stringr)
library(here)

# make paths
proj.path.12S <- here("/home/genomics/icornelis/02_ZEROimpact/01_12S/NJ2021/MiFish-UE_run2")
proj.path.COI <- here("/home/genomics/icornelis/02_ZEROimpact/02_COI/NJ2021")

#upload data
table_12S <- readxl::read_excel(paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean.xlsx"))
table_morph_Fish <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Standerdized.xlsx"),sheet = "Fish - Standerdized")
table_COI <- readxl::read_excel(paste0(proj.path.COI,"/OWFvsCoastal_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean.xlsx"))
table_morph_Inv <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Standerdized.xlsx"),sheet = "Epi - Standerdized")
env_12S <- readRDS(paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/R_Environment/env_ordered_noNeg.rds"))
env_COI <- readRDS(paste0(proj.path.COI,"/OWFvsCoastal_concatenated/results_microDecon/R_Environment/env_ordered_noNeg.rds"))
env_morph <- read.csv(paste0(proj.path.12S,"/Step5_Statistics/environmental_data_morph.csv"),  header=TRUE, sep=";")

#add color and pch to environmental data to create the plot 
env_12S$pch <- ifelse(env_12S$Environment=="inside_OWF","1","21")
env_12S$pch <- as.integer(env_12S$pch)
env_COI$pch <- ifelse(env_COI$Environment=="inside_OWF","1","21")
env_COI$pch <- as.integer(env_COI$pch)
env_morph$Zone_color <- ifelse(env_morph$Zones=="Coastal","limegreen", 
                                      ifelse(env_morph$Zones=="zone1", "slateblue", 
                                             ifelse(env_morph$Zones=="zone2","darkorange","red")))
env_morph$pch <- ifelse(env_morph$Environment=="inside_OWF","1","21")
env_morph$pch <- as.integer(env_morph$pch)

#create table for morphological data
table_morph_Fish_2 <- table_morph_Fish[,2:ncol(table_morph_Fish)]
tab_morph_Fish <- as.data.frame(t(table_morph_Fish_2))
colnames(tab_morph_Fish) <- table_morph_Fish$...1
tab_morph_Fish <- round(tab_morph_Fish, digits = 0)

table_morph_Inv_2 <- table_morph_Inv[,2:ncol(table_morph_Inv)]
tab_morph_Inv <- as.data.frame(t(table_morph_Inv_2))
colnames(tab_morph_Inv) <- table_morph_Inv$...1
tab_morph_Inv <- round(tab_morph_Inv, digits = 0)

### Fish data
# 12S select ASVs assigned to fishes
#select Fish species and merge by species
fish_classes <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_classes.rds"))
freshwater_fish <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_Freshwater.rds"))
table_FishASVs <- as.data.frame(table_12S[table_12S$Class %in% fish_classes,])
table_FishASVs <- as.data.frame(table_FishASVs[!table_FishASVs$Species %in% c(freshwater_fish, "NA"),])
table_FishASVs[is.na(table_FishASVs)] <- 0
taxo <- 'Species'
table_Fish <- aggregate(table_FishASVs[,1:(ncol(table_FishASVs)-11)], 
                                         by= list(as.factor(table_FishASVs[,taxo])),FUN=sum)
fish_species <- c(table_Fish$Group.1)
rownames(table_Fish) <- table_Fish$Group.1
table_Fish$Group.1 <- NULL
seqtab_Fish <- as.data.frame(t(table_Fish))
Fish_empty <- c(rownames(seqtab_Fish[rowSums(seqtab_Fish, na.rm=T) == '0',]))
Fish_zero <- c(colnames(seqtab_Fish[colSums(seqtab_Fish, na.rm=T) == '0']))
seqtab_Fish <- seqtab_Fish[,!colnames(seqtab_Fish) %in% Fish_zero]
seqtab_Fish <- seqtab_Fish[!rownames(seqtab_Fish) %in% Fish_empty,]
env_12S <- env_12S[!env_12S$Niskin.sample %in% Fish_empty,]

### Invertebrate data
# COI select ASVs assigned to invertebrates
table_AnimaliaASVs <- as.data.frame(table_COI[table_COI$Kingdom %in% c("Animalia"),])
table_AnimaliaASVs <- as.data.frame(table_AnimaliaASVs[!table_AnimaliaASVs$Species %in% "NA",])
table_AnimaliaASVs <- as.data.frame(table_AnimaliaASVs[!table_AnimaliaASVs$Phylum %in% "Chordata",])
table_AnimaliaASVs[is.na(table_AnimaliaASVs)] <- 0
taxo <- "Species"
remove_taxonomy <- ncol(table_AnimaliaASVs)-11
table_Animalia <- aggregate(table_AnimaliaASVs[,1:(ncol(table_AnimaliaASVs)-11)], by= list(as.factor(table_AnimaliaASVs[,taxo])),FUN=sum)
rownames(table_Animalia) <- table_Animalia$Group.1
table_Animalia$Group.1 <- NULL
seqtab_Animalia <- as.data.frame(t(table_Animalia))
Animalia_empty <- c(rownames(seqtab_Animalia[rowSums(seqtab_Animalia, na.rm=T) == 0,]))
Animalia_zero  <- c(colnames(seqtab_Animalia[,colSums(seqtab_Animalia, na.rm=T) == 0]))
seqtab_Animalia <- seqtab_Animalia[,!colnames(seqtab_Animalia) %in% Animalia_zero]
seqtab_Animalia <- seqtab_Animalia[!rownames(seqtab_Animalia) %in% Animalia_empty,]

env_COI <- env_COI[!env_COI$Niskin.sample %in% Animalia_empty,]

####### create nMDS-plots #######
##### Coast vs Offshore
par(mar=c(5,5,3,3), mfrow= c(2,2))
### 12S-data
#data.nmds_12S <- decostand(seqtab_Fish, method="log")
#data.nmds_12S <- wisconsin(seqtab_Fish)
#data.nmds_12S <- decostand(seqtab_Fish, method="hellinger")
data.nmds_12S <- decostand(decostand(seqtab_Fish, method="tot"), method="max")
ord.NMDS_12S=metaMDS(data.nmds_12S, k=2, distace ="bray", trymax=100, )
#ord.NMDS=metaMDS(seqtab_Fish, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS_12S, type='n', xlim=c(-2, 2), ylim=c(-2, 2))
title(main = "eDNA - Fish species", adj = 0)
points(ord.NMDS_12S,display="sites",pch=21,col="black",bg=as.vector(env_12S$Zone_color),cex=1)
legend("bottomleft", legend = c("Coastal","Transition","Offshore"), pch=c(16,16,16), border = "black", col=c("limegreen", "slateblue", "darkorange"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS_12S$stress,4),sep=""), bty ="n", cex=0.8)
#text(ord.NMDS_12S, display = "sites", cex = 0.5, color = as.vector(env_12S$Zone_color))

### COI-data
#data.nmds_COI <- decostand(seqtab_AnimaliaASVs, method="log")
#data.nmds_COI <- decostand(seqtab_AnimaliaASVs_samples, method="hellinger")
#data.nmds_COI <- wisconsin(seqtab_AnimaliaASVs_samples)
data.nmds_COI <- decostand(decostand(seqtab_Animalia, method="tot"), method="max")
ord.NMDS_COI=metaMDS(data.nmds_COI, k=2, distace ="bray", trymax=100)
#ord.NMDS_COI=metaMDS(seqtab_AnimaliaASVs, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS_COI, type='n', ylim=c(-2,2), xlim=c(-1,1))
title(main = "eDNA - Invertebrate species", adj = 0)
points(ord.NMDS_COI,display="sites",pch=21,col="black",bg=as.vector(env_COI$Zone_color),cex=1)
legend("bottomleft", legend = c("Coastal","Transition","Offshore"), pch=c(16,16,16), border = "black", col=c("limegreen", "slateblue", "darkorange"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS_COI$stress,4),sep=""), bty ="n", cex=0.8)

### Morphology data
#create NMDS with morphological data on Fish
data.nmds_Fish <- decostand(tab_morph_Fish, method="hellinger")
ord.NMDS_Fish=metaMDS(data.nmds_Fish, k=2, distace ="bray", trymax=100)
#ord.NMDS=metaMDS(tab_morph_Fish, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS_Fish, type='n', xlim=c(-1, 1), ylim=c(-1, 1))
title(main = "Morphology - Fish species", adj = 0)
points(ord.NMDS_Fish,display="sites",pch=21,col="black",bg=as.vector(env_morph$Zone_color),cex=1)
legend("bottomleft", legend = c("Coastal","Transition","Offshore"), pch=c(16,16,16), border = "black", col=c("limegreen", "slateblue", "darkorange"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS_Fish$stress,4),sep=""), bty ="n", cex=0.8)

#create NMDS with morphological data on invertebrates
data.nmds_Inv <- decostand(tab_morph_Inv, method="hellinger")
ord.NMDS_Inv=metaMDS(data.nmds_Inv, k=2, distace ="bray", trymax=100)
#ord.NMDS=metaMDS(tab_morph_Inv, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS_Inv, type='n', xlim=c(-1, 1), ylim=c(-1, 1))
title(main = "Morphology - Invertebrate species", adj = 0)
points(ord.NMDS_Inv,display="sites",pch=21,col="black",bg=as.vector(env_morph$Zone_color),cex=1)
legend("bottomleft", legend = c("Coastal","Transition","Offshore"), pch=c(16,16,16), border = "black", col=c("limegreen", "slateblue", "darkorange"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS_Inv$stress,4),sep=""), bty ="n", cex=0.8)
#text(ord.NMDS, display="sites",cex=0.5)

##### Samples taken at OWFs only
par(mar=c(3,3,3,3), mfrow= c(2,2))
### 12S-data
#create NMDS only inside and outside OWF with ASVs assigned to fish
env_12S_OWF <- env_12S[!env_12S$Zone %in% "Coast",]
remove_coastal_12S <- c(env_12S_OWF$Niskin.sample)
seqtab_Fish_OWF <- seqtab_Fish[rownames(seqtab_Fish) %in% remove_coastal_12S,]
seqtab_Fish_OWF <- as.data.frame(seqtab_Fish_OWF[,!colSums(seqtab_Fish_OWF) == 0,])
data.nmds_OWF_12S  <- data.nmds_12S[rownames(data.nmds_12S)%in% remove_coastal_12S,]
ord.NMDS_12S_OWF=metaMDS(data.nmds_OWF_12S, k=2, distace ="bray", trymax=100)
#ord.NMDS=metaMDS(seqtab_FishASVs_OWF, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS_12S_OWF, type='n', xlim=c(-1.5,1.5), ylim=c(-1.5,1.5))
title(main = "eDNA - Fish species", adj = 0)
points(ord.NMDS_12S_OWF,display="sites",pch=as.vector(env_12S_OWF$pch),col=as.vector(env_12S_OWF$Zone_color),bg=as.vector(env_12S_OWF$Zone_color),cex=1,)
legend("bottomleft", legend = c("Transition","Offshore", "outside_OWF", "inside_OWF"), pch = c(16,16,16,1), col=c("slateblue", "darkorange", "black", "black"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS_12S_OWF$stress,4),sep=""), bty ="n", cex=0.8)
#text(ord.NMDS, display = "sites")

### COI-data
#create NMDS only inside and outside OWF with ASVs assigned to invertebrates
env_COI_OWF <- env_COI[!env_COI$Zone %in% "Coast",]
remove_coastal_COI <- c(env_COI_OWF$Niskin.sample)
seqtab_Animalia_OWF <- seqtab_Animalia[rownames(seqtab_Animalia) %in% remove_coastal_COI,]
data.nmds_OWF_COI <- data.nmds_COI[rownames(data.nmds_COI) %in% remove_coastal_COI,]
ord.NMDS_OWF_COI=metaMDS(data.nmds_OWF_COI, k=2, distace ="bray", trymax=100)
#ord.NMDS=metaMDS(seqtab_AnimaliaASVs_OWF, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS_OWF_COI, type='n', xlim=c(-3,5), ylim=c(-2,2))
title(main = "eDNA - Invertebrate species", adj = 0)
points(ord.NMDS_OWF_COI,display="sites",pch=as.vector(env_COI_OWF$pch),col=as.vector(env_COI_OWF$Zone_color),bg=as.vector(env_COI_OWF$Zone_color),cex=1,)
legend("bottomleft", legend = c("Transition","Offshore", "outside_OWF", "inside_OWF"), pch = c(16,16,16,1), col=c("slateblue", "darkorange", "black", "black"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS_OWF_COI$stress,4),sep=""), bty ="n", cex=0.8)
#text(ord.NMDS, display="sites",cex=0.5)

### Morphology data
#create NMDS with morphological data on Fish
env_Fish_OWF <- env_morph[!env_morph$Zones %in% "Coastal",]
remove_coastal_Fish <- c(env_Fish_OWF$Niskin.sample)
tab_morph_Fish_OWF <- tab_morph_Fish[rownames(tab_morph_Fish) %in% remove_coastal_Fish,]
data.nmds_OWF_Fish <- decostand(tab_morph_Fish_OWF, method="hellinger")
ord.NMDS_OWF_Fish=metaMDS(data.nmds_OWF_Fish, k=2, distace ="bray", trymax=100)
#ord.NMDS=metaMDS(tab_morph_Fish_OWF, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS_OWF_Fish, type='n', xlim=c(-0.75, 0.75), ylim=c(-0.75, 0.75))
title(main = "Morphology - Fish species", adj = 0)
points(ord.NMDS_OWF_Fish,display="sites",pch=as.vector(env_Fish_OWF$pch),col=as.vector(env_Fish_OWF$Zone_color),bg=as.vector(env_Fish_OWF$Zone_color),cex=1,)
legend("bottomleft", legend = c("Transition","Offshore", "outside_OWF", "inside_OWF"), pch = c(16,16,16,1), col=c("slateblue", "darkorange", "black", "black"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS_OWF_Fish$stress,4),sep=""), bty ="n", cex=0.8)

#create NMDS with morphological data on invertebrates
env_Inv_OWF <- env_morph[!env_morph$Zones %in% "Coastal",]
remove_coastal_Inv <- c(env_Inv_OWF$Niskin.sample)
tab_morph_Inv_OWF <- tab_morph_Inv[rownames(tab_morph_Inv) %in% remove_coastal_Inv,]
data.nmds_OWF_Inv <- decostand(tab_morph_Inv_OWF, method="hellinger")
ord.NMDS_OWF_Inv=metaMDS(data.nmds_OWF_Inv, k=2, distace ="bray", trymax=100)
#ord.NMDS=metaMDS(tab_morph_Inv_OWF, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS_OWF_Inv, type='n', xlim=c(-0.75, 0.75), ylim=c(-0.75, 0.75))
title(main = "Morphology - Invertebrate species", adj = 0)
points(ord.NMDS_OWF_Inv,display="sites",pch=as.vector(env_Inv_OWF$pch),col=as.vector(env_Inv_OWF$Zone_color),bg=as.vector(env_Inv_OWF$Zone_color),cex=1,)
legend("bottomleft", legend = c("Transition","Offshore", "outside_OWF", "inside_OWF"), pch = c(16,16,16,1), col=c("slateblue", "darkorange", "black", "black"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS_OWF_Inv$stress,4),sep=""), bty ="n", cex=0.8)

###Indicator species analysis 3 zones
library(indicspecies)
groups_12S <- c(env_12S_samples$Zones)
#indval_12S <- multipatt(seqtab_FishASVs_samples, groups_12S, func = "IndVal", control = how(nperm=9999), duleg = TRUE)
indval_12S <- multipatt(data.nmds_12S, groups_12S, func = "IndVal", control = how(nperm=9999), duleg = TRUE)
simper_12S <- simper(data.nmds_12S, groups_12S, permutations = 9999)
summary(indval_12S) 
summary(simper_12S) 
groups_COI <- c(env_COI_samples$Zones)
#indval_COI <- multipatt(seqtab_AnimaliaASVs_samples, groups_COI, func = "IndVal.g", control = how(nperm=9999), duleg = TRUE) 
indval_COI <- multipatt(data.nmds_COI, groups_COI, func = "IndVal.g", control = how(nperm=9999), duleg = TRUE) 
simper_COI <- simper(data.nmds_COI, groups_COI, permutations = 9999)
summary(indval_COI) 
summary(simper_COI)
groups_morph <- c(env_morph$Zones)
#indval_Fish <- multipatt(tab_morph_Fish, groups_morph, func = "IndVal.g", control = how(nperm=9999), duleg = TRUE)
indval_Fish <- multipatt(data.nmds_Fish, groups_morph, func = "IndVal.g", control = how(nperm=9999), duleg = TRUE)
simper_Fish <- simper(data.nmds_Fish, groups_morph, permutations = 9999)
summary(indval_Fish)
summary(simper_Fish)
#indval_Inv <- multipatt(tab_morph_Inv, groups_morph, func = "IndVal.g", control = how(nperm=9999), duleg = TRUE)
indval_Inv <- multipatt(data.nmds_Inv, groups_morph, func = "IndVal.g", control = how(nperm=9999), duleg = TRUE)
simper_Inv <- simper(data.nmds_Inv, groups_morph, permutations = 9999)
summary(indval_Inv)
summary(simper_Inv)

###Indicator species analysis OWF_only
env_12S_OWF$Env_Zone <- paste(env_12S_OWF$Environment,env_12S_OWF$Zones, sep="_")
groups_OWF_12S <- c(env_12S_OWF$Env_Zone)
indval_OWF_12S <- multipatt(data.nmds_OWF_12S, groups_OWF_12S, func = "IndVal.g", control = how(nperm=9999))
simper_OWF_12S <- simper(data.nmds_OWF_12S, groups_OWF_12S, permutations = 9999)
summary(indval_OWF_12S) 
summary(simper_OWF_12S)
env_COI_OWF$Env_Zone <- paste(env_COI_OWF$Environment,env_COI_OWF$Zones, sep="_")
groups_OWF_COI <- c(env_COI_OWF$Env_Zone)
indval_OWF_COI <- multipatt(data.nmds_OWF_COI, groups_OWF_COI, func = "IndVal.g", control = how(nperm=9999)) 
simper_OWF_COI <- simper(data.nmds_OWF_COI, groups_OWF_COI, permutations = 9999)
summary(indval_OWF_COI)
summary(simper_OWF_COI)
env_Fish_OWF$Env_Zone <- paste(env_Fish_OWF$Environment,env_Fish_OWF$Zones, sep="_")
groups_OWF_morph <- c(env_Fish_OWF$Env_Zone)
indval_OWF_Fish <- multipatt(data.nmds_OWF_Fish, groups_OWF_morph, func = "IndVal.g", control = how(nperm=9999))
simper_OWF_Fish <- simper(data.nmds_OWF_Fish, groups_OWF_morph, permutations = 9999)
summary(indval_OWF_Fish)
summary(simper_OWF_Fish)
indval_OWF_Inv <- multipatt(data.nmds_OWF_Inv, groups_OWF_morph, func = "IndVal.g", control = how(nperm=9999))
simper_OWF_Inv <- simper(data.nmds_OWF_Inv, groups_OWF_morph, permutations = 9999)
summary(indval_OWF_Inv)
summary(simper_OWF_Inv)
