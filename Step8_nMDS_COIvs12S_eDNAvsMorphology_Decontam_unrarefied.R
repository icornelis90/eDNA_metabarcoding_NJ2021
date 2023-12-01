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
table_unrarefied_12S <- readxl::read_excel(paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/table_unrarefied_concatenated_CleanedASVs_WithField_FullTaxonomicAssignment.xlsx"))
#table_morph_Fish <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Standerdized.xlsx"),sheet = "Fish")
table_morph_Fish <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Standerdized.xlsx"),sheet = "Fish - Standerdized")
table_unrarefied_COI <- readxl::read_excel(paste0(proj.path.COI,"/OWFvsCoastal_concatenated/results/Decontam/table_unrarefied_concatenated_CleanedASVs_FullTaxonomicAssignment_WithField.xlsx"))
table_morph_Inv <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Standerdized.xlsx"),sheet = "Epi - Standerdized")
env <- read.csv(paste0(proj.path.12S,"/Step5_Statistics/environmental_data.csv"),  header=TRUE, sep=";")
env_morph <- read.csv(paste0(proj.path.12S,"/Step5_Statistics/environmental_data_morph.csv"),  header=TRUE, sep=";")

#add color and pch to environmental data to create the plot 
env$Environment_color <- ifelse(env$Zones=="Coastal","limegreen", 
                                ifelse(env$Zones=="zone1", "slateblue", 
                                       ifelse(env$Zones=="zone2","darkorange","red")))
#env$Environment_color <- ifelse(env$Niskin.sample %in% c("ft230_1","ft230_2","ft230_3"),
#                                   "red", env$Environment_color )
env$pch <- ifelse(env$Environment=="inside_OWF","1","21")
#env$pch <- ifelse(env$Environment=="inside_OWF","24","21")
env$pch <- as.integer(env$pch)
env_morph$Environment_color <- ifelse(env_morph$Zones=="Coastal","limegreen", 
                                      ifelse(env_morph$Zones=="zone1", "slateblue", 
                                             ifelse(env_morph$Zones=="zone2","darkorange","red")))
#env_morph$Environment_color <- ifelse(env_morph$Niskin.sample %in% c("ft230"),
#                                "red", env_morph$Environment_color )
env_morph$pch <- ifelse(env_morph$Environment=="inside_OWF","1","21")
#env$pch <- ifelse(env$Environment=="inside_OWF","24","21")
env_morph$pch <- as.integer(env_morph$pch)

#create table for morphological data
table_morph_Fish_2 <- table_morph_Fish[,2:ncol(table_morph_Fish)]
#rownames(table_morph_Fish_2) <- table_morph_Fish$...1
tab_morph_Fish <- as.data.frame(t(table_morph_Fish_2))
#tab_morph_Fish[is.na(tab_morph_Fish)] <- 0
colnames(tab_morph_Fish) <- table_morph_Fish$...1
tab_morph_Fish <- round(tab_morph_Fish, digits = 0)

table_morph_Inv_2 <- table_morph_Inv[,2:ncol(table_morph_Inv)]
#rownames(table_morph_Inv_2) <- table_morph_Inv$...1
tab_morph_Inv <- as.data.frame(t(table_morph_Inv_2))
#tab_morph_Inv[is.na(tab_morph_Inv)] <- 0
colnames(tab_morph_Inv) <- table_morph_Inv$...1
tab_morph_Inv <- round(tab_morph_Inv, digits = 0)

#remove samples from environmental data that are absent from the unrarefied data
keep_samples_12S <- c(colnames(table_unrarefied_12S))
env_12S <- env[env$Niskin.sample %in% keep_samples_12S,]
env_12S <- env_12S[!env_12S$Niskin.sample %in% ("ftWT1bis_4"),]
colorder_12S <- c(env_12S$Niskin.sample)
keep_samples_COI <- c(colnames(table_unrarefied_COI))
env_COI <- env[env$Niskin.sample %in% keep_samples_COI,]
colorder_COI <- c(env_COI$Niskin.sample)

### Fish data
# 12S select ASVs assigned to fishes
#select Fish species and merge by species
fish_classes <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_classes.rds"))
freshwater_fish <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_Freshwater.rds"))
table_unrarefied_FishASVs <- as.data.frame(table_unrarefied_12S[table_unrarefied_12S$Class %in% fish_classes,])
table_unrarefied_FishASVs <- as.data.frame(table_unrarefied_FishASVs[!table_unrarefied_FishASVs$Species %in% c(freshwater_fish, "NA"),])
#table_unrarefied_FishASVs <- as.data.frame(table_unrarefied_FishASVs[,!colnames(table_unrarefied_FishASVs) %in% c("ftWT1bis_4")])
taxo <- 'Species'
merged_data_unrarefied_fish <- aggregate(table_unrarefied_FishASVs[,1:(ncol(table_unrarefied_FishASVs)-11)], 
                                         by= list(as.factor(table_unrarefied_FishASVs[,taxo])),FUN=sum)
fish_species <- c(merged_data_unrarefied_fish$Group.1)
rownames(merged_data_unrarefied_fish) <- merged_data_unrarefied_fish$Group.1
merged_data_unrarefied_fish$Group.1 <- NULL
seqtab_unrarefied_FishASVs <- as.data.frame(t(merged_data_unrarefied_fish[,1:(ncol(table_unrarefied_FishASVs)-11)]))
colnames(seqtab_unrarefied_FishASVs) <- rownames(merged_data_unrarefied_fish)
seqtab_unrarefied_FishASVs  <- seqtab_unrarefied_FishASVs[colorder_12S,]
#seqtab_unrarefied_FishASVs  <- seqtab_unrarefied_FishASVs[c(env_12S_samples$Niskin.sample),]

### Invertebrate data
# COI select ASVs assigned to invertebrates
table_unrarefied_AnimaliaASVs <- as.data.frame(table_unrarefied_COI[table_unrarefied_COI$Kingdom %in% c("Animalia"),])
table_unrarefied_AnimaliaASVs <- as.data.frame(table_unrarefied_AnimaliaASVs[!table_unrarefied_AnimaliaASVs$Full %in% "NA",])
#table_unrarefied_AnimaliaASVs <- as.data.frame(table_unrarefied_AnimaliaASVs[!table_unrarefied_AnimaliaASVs$Full %in% c("NA","Clathria prolifera"),])
table_unrarefied_AnimaliaASVs <- as.data.frame(table_unrarefied_AnimaliaASVs[!table_unrarefied_AnimaliaASVs$Phylum %in% "Chordata",])
taxo <- "Full"
remove_taxonomy <- ncol(table_unrarefied_AnimaliaASVs)-11
merged_data_unrarefied_Animalia<-aggregate(table_unrarefied_AnimaliaASVs[,1:remove_taxonomy], by= list(as.factor(table_unrarefied_AnimaliaASVs[,taxo])),FUN=sum)
rownames(merged_data_unrarefied_Animalia) <- merged_data_unrarefied_Animalia$Group.1
merged_data_unrarefied_Animalia$Group.1 <- NULL
seqtab_unrarefied_AnimaliaASVs <- as.data.frame(t(table_unrarefied_AnimaliaASVs[,1:remove_taxonomy]))
seqtab_unrarefied_AnimaliaASVs <- as.data.frame(t(merged_data_unrarefied_Animalia))
seqtab_unrarefied_AnimaliaASVs  <- seqtab_unrarefied_AnimaliaASVs[colorder_COI,]
remove_empty <- c(rownames(seqtab_unrarefied_AnimaliaASVs[rowSums(seqtab_unrarefied_AnimaliaASVs) > '0',]))
remove_zero <- c(colnames(seqtab_unrarefied_AnimaliaASVs[,colSums(seqtab_unrarefied_AnimaliaASVs) > '0']))
seqtab_unrarefied_AnimaliaASVs <- seqtab_unrarefied_AnimaliaASVs[,colnames(seqtab_unrarefied_AnimaliaASVs) %in% remove_zero]
seqtab_unrarefied_AnimaliaASVs <- seqtab_unrarefied_AnimaliaASVs[rownames(seqtab_unrarefied_AnimaliaASVs) %in% remove_empty,]
env_COI <- env_COI[env_COI$Niskin.sample %in% remove_empty,]

####### create nMDS-plots #######
##### Coast vs Offshore
par(mar=c(5,5,3,3), mfrow= c(2,2))
### 12S-data
#create NMDS without neg_controls and only ASVs assigned to Fish
env_12S_samples <- env_12S %>% filter(!grepl("neg", env_12S$Niskin.sample))
seqtab_unrarefied_FishASVs_samples <- seqtab_unrarefied_FishASVs %>% filter(!grepl("neg", rownames(seqtab_unrarefied_FishASVs)))
#data.nmds_12S <- decostand(seqtab_unrarefied_FishASVs_samples, method="log")
#data.nmds_12S <- wisconsin(seqtab_unrarefied_FishASVs_samples)
#data.nmds_12S <- decostand(seqtab_unrarefied_FishASVs_samples, method="hellinger")
data.nmds_12S <- decostand(decostand(seqtab_unrarefied_FishASVs_samples, method="tot"), method="max")
ord.NMDS=metaMDS(data.nmds_12S, k=2, distace ="bray", trymax=100)
#ord.NMDS=metaMDS(seqtab_unrarefied_FishASVs, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS, type='n', xlim=c(-2, 2), ylim=c(-2, 2))
title(main = "eDNA - Fish species", adj = 0)
points(ord.NMDS,display="sites",pch=21,col="black",bg=as.vector(env_12S_samples$Environment_color),cex=1)
legend("bottomleft", legend = c("Coastal","Transition","Offshore"), pch=c(16,16,16), border = "black", col=c("limegreen", "slateblue", "darkorange"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS$stress,4),sep=""), bty ="n", cex=0.8)
#text(ord.NMDS, display = "sites", cex = 0.5, color = as.vector(env_12S_samples$Environment_color))

### COI-data
# create NMDS unrarefied samples with ASVs assigned to Animalia
env_COI_samples <- env_COI %>% filter(!grepl("neg", env_COI$Niskin.sample))
seqtab_unrarefied_AnimaliaASVs_samples <- seqtab_unrarefied_AnimaliaASVs %>% filter(!grepl("neg", rownames(seqtab_unrarefied_AnimaliaASVs)))
#data.nmds_COI <- decostand(seqtab_unrarefied_AnimaliaASVs, method="log")
#data.nmds_COI <- decostand(seqtab_unrarefied_AnimaliaASVs_samples, method="hellinger")
#data.nmds_COI <- wisconsin(seqtab_unrarefied_AnimaliaASVs_samples)
data.nmds_COI <- decostand(decostand(seqtab_unrarefied_AnimaliaASVs_samples, method="tot"), method="max")
ord.NMDS=metaMDS(data.nmds_COI, k=2, distace ="bray", trymax=100)
#ord.NMDS=metaMDS(seqtab_unrarefied_AnimaliaASVs, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS, type='n', ylim=c(-3,3), xlim=c(-3,3))
title(main = "eDNA - Invertebrate species", adj = 0)
points(ord.NMDS,display="sites",pch=21,col="black",bg=as.vector(env_COI_samples$Environment_color),cex=1)
legend("bottomleft", legend = c("Coastal","Transition","Offshore"), pch=c(16,16,16), border = "black", col=c("limegreen", "slateblue", "darkorange"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS$stress,4),sep=""), bty ="n", cex=0.8)

### Morphology data
#create NMDS with morphological data on Fish
data.nmds_Fish <- decostand(tab_morph_Fish, method="hellinger")
ord.NMDS=metaMDS(data.nmds_Fish, k=2, distace ="bray", trymax=100)
#ord.NMDS=metaMDS(tab_morph_Fish, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS, type='n', xlim=c(-1, 1), ylim=c(-1, 1))
title(main = "Morphology - Fish species", adj = 0)
points(ord.NMDS,display="sites",pch=21,col="black",bg=as.vector(env_morph$Environment_color),cex=1)
legend("bottomleft", legend = c("Coastal","Transition","Offshore"), pch=c(16,16,16), border = "black", col=c("limegreen", "slateblue", "darkorange"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS$stress,4),sep=""), bty ="n", cex=0.8)

#create NMDS with morphological data on invertebrates
data.nmds_Inv <- decostand(tab_morph_Inv, method="hellinger")
ord.NMDS=metaMDS(data.nmds_Inv, k=2, distace ="bray", trymax=100)
#ord.NMDS=metaMDS(tab_morph_Inv, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS, type='n', xlim=c(-1, 1), ylim=c(-1, 1))
title(main = "Morphology - Invertebrate species", adj = 0)
points(ord.NMDS,display="sites",pch=21,col="black",bg=as.vector(env_morph$Environment_color),cex=1)
legend("bottomleft", legend = c("Coastal","Transition","Offshore"), pch=c(16,16,16), border = "black", col=c("limegreen", "slateblue", "darkorange"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS$stress,4),sep=""), bty ="n", cex=0.8)
#text(ord.NMDS, display="sites",cex=0.5)

##### Samples taken at OWFs only
par(mar=c(3,3,3,3), mfrow= c(2,2))
### 12S-data
#create NMDS only inside and outside OWF with ASVs assigned to fish
env_12S_OWF <- env_12S_samples[!env_12S_samples$Zones %in% "Coastal",]
remove_coastal_12S <- c(env_12S_OWF$Niskin.sample)
seqtab_unrarefied_FishASVs_OWF <- seqtab_unrarefied_FishASVs_samples[rownames(seqtab_unrarefied_FishASVs_samples) %in% remove_coastal_12S,]
seqtab_unrarefied_FishASVs_OWF <- as.data.frame(seqtab_unrarefied_FishASVs_OWF[,!colSums(seqtab_unrarefied_FishASVs_samples) == 0,])
data.nmds_OWF_12S  <- data.nmds_12S[rownames(data.nmds_12S)%in% remove_coastal_12S,]
ord.NMDS=metaMDS(data.nmds_OWF_12S, k=2, distace ="bray", trymax=100)
#ord.NMDS=metaMDS(seqtab_unrarefied_FishASVs_OWF, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS, type='n', xlim=c(-1.5,1.5), ylim=c(-1.5,1.5))
title(main = "eDNA - Fish species", adj = 0)
points(ord.NMDS,display="sites",pch=as.vector(env_12S_OWF$pch),col=as.vector(env_12S_OWF$Environment_color),bg=as.vector(env_12S_OWF$Environment_color),cex=1,)
legend("bottomleft", legend = c("Transition","Offshore", "outside_OWF", "inside_OWF"), pch = c(16,16,16,1), col=c("slateblue", "darkorange", "black", "black"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS$stress,4),sep=""), bty ="n", cex=0.8)
#text(ord.NMDS, display = "sites")

### COI-data
#create NMDS only inside and outside OWF with ASVs assigned to invertebrates
env_COI_OWF <- env_COI[!env_COI$Zones %in% c("Coastal", "neg_control"),]
remove_coastal_COI <- c(env_COI_OWF$Niskin.sample)
seqtab_unrarefied_AnimaliaASVs_OWF <- seqtab_unrarefied_AnimaliaASVs[rownames(seqtab_unrarefied_AnimaliaASVs) %in% remove_coastal_COI,]
data.nmds_OWF_COI <- data.nmds_COI[rownames(data.nmds_COI) %in% remove_coastal_COI,]
ord.NMDS=metaMDS(data.nmds_OWF_COI, k=2, distace ="bray", trymax=100)
#ord.NMDS=metaMDS(seqtab_unrarefied_AnimaliaASVs_OWF, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
show_axes <- c(1,2)
plot(ord.NMDS, type='n', xlim=c(-3,5), ylim=c(-2,2))
title(main = "eDNA - Invertebrate species", adj = 0)
points(ord.NMDS,display="sites",pch=as.vector(env_COI_OWF$pch),col=as.vector(env_COI_OWF$Environment_color),bg=as.vector(env_COI_OWF$Environment_color),cex=1,)
legend("bottomleft", legend = c("Transition","Offshore", "outside_OWF", "inside_OWF"), pch = c(16,16,16,1), col=c("slateblue", "darkorange", "black", "black"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS$stress,4),sep=""), bty ="n", cex=0.8)
#text(ord.NMDS, display="sites",cex=0.5)

### Morphology data
#create NMDS with morphological data on Fish
env_Fish_OWF <- env_morph[!env_morph$Zones %in% "Coastal",]
remove_coastal_Fish <- c(env_Fish_OWF$Niskin.sample)
tab_morph_Fish_OWF <- tab_morph_Fish[rownames(tab_morph_Fish) %in% remove_coastal_Fish,]
data.nmds_OWF_Fish <- decostand(tab_morph_Fish_OWF, method="hellinger")
ord.NMDS=metaMDS(data.nmds_OWF_Fish, k=2, distace ="bray", trymax=100)
#ord.NMDS=metaMDS(tab_morph_Fish_OWF, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS, type='n', xlim=c(-0.75, 0.75), ylim=c(-0.75, 0.75))
title(main = "Morphology - Fish species", adj = 0)
points(ord.NMDS,display="sites",pch=as.vector(env_Fish_OWF$pch),col=as.vector(env_Fish_OWF$Environment_color),bg=as.vector(env_Fish_OWF$Environment_color),cex=1,)
legend("bottomleft", legend = c("Transition","Offshore", "outside_OWF", "inside_OWF"), pch = c(16,16,16,1), col=c("slateblue", "darkorange", "black", "black"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS$stress,4),sep=""), bty ="n", cex=0.8)

#create NMDS with morphological data on invertebrates
env_Inv_OWF <- env_morph[!env_morph$Zones %in% "Coastal",]
remove_coastal_Inv <- c(env_Inv_OWF$Niskin.sample)
tab_morph_Inv_OWF <- tab_morph_Inv[rownames(tab_morph_Inv) %in% remove_coastal_Inv,]
data.nmds_OWF_Inv <- decostand(tab_morph_Inv_OWF, method="hellinger")
ord.NMDS=metaMDS(data.nmds_OWF_Inv, k=2, distace ="bray", trymax=100)
#ord.NMDS=metaMDS(tab_morph_Inv_OWF, k=2, distace ="bray", trymax=100)
#par(mar=c(5,5,3,3))
plot(ord.NMDS, type='n', xlim=c(-0.75, 0.75), ylim=c(-0.75, 0.75))
title(main = "Morphology - Invertebrate species", adj = 0)
points(ord.NMDS,display="sites",pch=as.vector(env_Inv_OWF$pch),col=as.vector(env_Inv_OWF$Environment_color),bg=as.vector(env_Inv_OWF$Environment_color),cex=1,)
legend("bottomleft", legend = c("Transition","Offshore", "outside_OWF", "inside_OWF"), pch = c(16,16,16,1), col=c("slateblue", "darkorange", "black", "black"), bty = "n", cex = 0.8)
legend("bottomright", legend = paste("NMDS, stress=",round(ord.NMDS$stress,4),sep=""), bty ="n", cex=0.8)

###Indicator species analysis 3 zones
library(indicspecies)
groups_12S <- c(env_12S_samples$Zones)
#indval_12S <- multipatt(seqtab_unrarefied_FishASVs_samples, groups_12S, func = "IndVal", control = how(nperm=9999), duleg = TRUE)
indval_12S <- multipatt(data.nmds_12S, groups_12S, func = "IndVal", control = how(nperm=9999), duleg = TRUE)
simper_12S <- simper(data.nmds_12S, groups_12S, permutations = 9999)
summary(indval_12S) 
summary(simper_12S) 
groups_COI <- c(env_COI_samples$Zones)
#indval_COI <- multipatt(seqtab_unrarefied_AnimaliaASVs_samples, groups_COI, func = "IndVal.g", control = how(nperm=9999), duleg = TRUE) 
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
