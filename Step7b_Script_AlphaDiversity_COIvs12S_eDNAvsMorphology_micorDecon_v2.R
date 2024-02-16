#set OS type for paths
if(.Platform$OS.type == "unix"){home<-"/home/"} else{
  home<-"//192.168.236.131/"}

#activate libraries
libraries <- c("BiocManager"
               , "ggplot2"
               , "tidyr"
               , "taxonomizr"
               , "naturalsort"
               , "stringr"
               , "phylotools"
               , "scales"
               , "ggpattern"
               , "ggh4x"
               , "reshape2"
               , "vegan"
               , "here"
               , "dplyr"
               , "phyloseq"
               , "ggpubr"
               , "metagMisc"
               , "iNEXT"
)

for (opties in libraries){
  
  if (opties %in% installed.packages()){
    
    library(opties,character.only = TRUE)
    
  } else {install.packages(opties,repos = "http://cran.us.r-project.org")
    
    library(opties,character.only = TRUE)
  }
}

# make paths
proj.path.12S <- here("/home/genomics/icornelis/02_ZEROimpact/01_12S/NJ2021/MiFish-UE_run2")
proj.path.COI <- here("/home/genomics/icornelis/02_ZEROimpact/02_COI/NJ2021")

#upload data
table_12S <- readxl::read_excel(paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean.xlsx"))
table_morph_Fish <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy.xlsx"),sheet = "Fish")

table_COI <- readxl::read_excel(paste0(proj.path.COI,"/OWFvsCoastal_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean.xlsx"))
table_morph_Inv <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy.xlsx"),sheet = "Epi")

env_12S_all <- readRDS(paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/R_Environment/env_ordered_noNeg.rds"))
env_COI_all <- readRDS(paste0(proj.path.COI,"/OWFvsCoastal_concatenated/results_microDecon/R_Environment/env_ordered_noNeg.rds"))
env_12S <- read.csv(paste0(proj.path.12S,"/Step5_Statistics/environmental_data_OWF3HighestReadNumbers_LocationsTrawl.csv"),  header=TRUE, sep=";")
env_COI <- read.csv(paste0(proj.path.COI,"/Step5_Statistics/environmental_data_OWF3HighestReadNumbers_LocationsTrawl.csv"),  header=TRUE, sep=";")
env_morph <- read.csv(paste0(proj.path.12S,"/Step5_Statistics/environmental_data_morph.csv"),  header=TRUE, sep=";")

#add color and pch to environmental data to create the plot 
env_morph$Zone_color <- ifelse(env_morph$Zone=="Coast","limegreen", 
                               ifelse(env_morph$Zone=="Transition", "slateblue", 
                                       ifelse(env_morph$Zone=="Offshore","darkorange","red")))

env_12S$Environment_color <- env_12S$Zone_color
env_12S$Zone_color <- ifelse(env_12S$Zone=="Coast","limegreen", 
                               ifelse(env_12S$Zone=="Transition", "slateblue", 
                                      ifelse(env_12S$Zone=="Offshore","darkorange","red")))


#morphological data
table_morph_Fish_2 <- as.data.frame(table_morph_Fish[,2:ncol(table_morph_Fish)])
table_morph_Fish_2 [is.na(table_morph_Fish_2 )] <- 0
rownames(table_morph_Fish_2) <- as.character(table_morph_Fish$...1)
table_morph_Fish_2 <- round(table_morph_Fish_2, digits = 0)
table_morph_Inv_2 <- as.data.frame(table_morph_Inv[,2:ncol(table_morph_Inv)])
table_morph_Inv_2 [is.na(table_morph_Inv_2)] <- 0
rownames(table_morph_Inv_2) <- as.character(table_morph_Inv$...1)
table_morph_Inv_2 <- round(table_morph_Inv_2, digits = 0)

#select Fish species and merge by species
fish_classes <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_classes.rds"))
freshwater_fish <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_Freshwater.rds"))
table_unrarefied_FishASVs <- as.data.frame(table_12S[table_12S$Class %in% fish_classes,])
table_unrarefied_FishASVs <- as.data.frame(table_unrarefied_FishASVs[
  !table_unrarefied_FishASVs$Species %in% c(freshwater_fish, "NA"),])
rownames(table_unrarefied_FishASVs) <- table_unrarefied_FishASVs$ASV

table_unrarefied_FishASVs_noT <- as.data.frame(table_unrarefied_FishASVs[,1:(ncol(table_unrarefied_FishASVs)-11)])
table_unrarefied_FishASVs_noT[is.na(table_unrarefied_FishASVs_noT)] <- 0
table_unrarefied_FishASVs_noT <- table_unrarefied_FishASVs_noT[!rowSums(table_unrarefied_FishASVs_noT) == 0,]
table_unrarefied_FishASVs_noT <- table_unrarefied_FishASVs_noT[,!colSums(table_unrarefied_FishASVs_noT) == 0]

taxo <- "Species"
table_unrarefied_Fish <- aggregate(table_unrarefied_FishASVs[,1:(ncol(table_unrarefied_FishASVs)-11)], 
                                   by= list(as.factor(table_unrarefied_FishASVs[,taxo])),FUN=sum)
rownames(table_unrarefied_Fish) <- as.character(table_unrarefied_Fish$Group.1)
table_unrarefied_Fish$Group.1 <- NULL
table_unrarefied_Fish[is.na(table_unrarefied_Fish)] <- 0
table_unrarefied_Fish <- table_unrarefied_Fish[!rowSums(table_unrarefied_Fish) == 0,]
table_unrarefied_Fish <- table_unrarefied_Fish[,!colSums(table_unrarefied_Fish) == 0]

#select Invertebrate species and merge by species
table_unrarefied_InvASVs <- as.data.frame(table_COI[table_COI$Kingdom %in% c("Animalia"),])
table_unrarefied_InvASVs <- as.data.frame(table_unrarefied_InvASVs[!table_unrarefied_InvASVs$Species %in% "NA",])
table_unrarefied_InvASVs <- as.data.frame(table_unrarefied_InvASVs[!table_unrarefied_InvASVs$Phylum %in% "Chordata",])

table_unrarefied_InvASVs_noT <- as.data.frame(table_unrarefied_InvASVs[,1:(ncol(table_unrarefied_InvASVs)-11)])
table_unrarefied_InvASVs_noT[is.na(table_unrarefied_InvASVs_noT)] <- 0
table_unrarefied_InvASVs_noT <- table_unrarefied_InvASVs_noT[!rowSums(table_unrarefied_InvASVs_noT) == 0,]
table_unrarefied_InvASVs_noT <- table_unrarefied_InvASVs_noT[,!colSums(table_unrarefied_InvASVs_noT) == 0]

taxo <- "Species"
table_unrarefied_Inv <- aggregate(table_unrarefied_InvASVs[,1:(ncol(table_unrarefied_InvASVs)-11)],
                                  by= list(as.factor(table_unrarefied_InvASVs[,taxo])),FUN=sum)
rownames(table_unrarefied_Inv) <-as.character(table_unrarefied_Inv$Group.1)
table_unrarefied_Inv$Group.1 <- NULL
table_unrarefied_Inv[is.na(table_unrarefied_Inv)] <- 0
table_unrarefied_Inv <- table_unrarefied_Inv[!rowSums(table_unrarefied_Inv) == 0,]
table_unrarefied_Inv <- table_unrarefied_Inv[,!colSums(table_unrarefied_Inv) == 0]

#Calculate the Coverage of Fish species level for all samples
#smpl_Fish <- env_12S_all[env_12S_all$Niskin.sample  %in%  colnames(table_unrarefied_Fish),]
smpl_Fish <- env_12S[env_12S$Niskin.sample  %in%  colnames(table_unrarefied_Fish),]
ps_Fish <- table_unrarefied_Fish[,colnames(table_unrarefied_Fish) %in% smpl_Fish$Niskin.sample]
rownames(smpl_Fish) <- colnames(ps_Fish)
Taxonomy_Fish <- as.matrix(rownames(ps_Fish))
rownames(Taxonomy_Fish) <- rownames(ps_Fish)

ps_unrarefied_Fish <- phyloseq(otu_table(ps_Fish, taxa_are_rows = TRUE), sample_data(smpl_Fish),tax_table(Taxonomy_Fish))
Coverage_Fish_all <- phyloseq_coverage(physeq = ps_unrarefied_Fish)
Coverage_Fish <- phyloseq_coverage(physeq = ps_unrarefied_Fish)

Fish_Diversity <- estimateD(ps_Fish, datatype="abundance",
                           base="coverage", level=0.9992178, conf=0.95)

#Calculate the Coverage of Invertebrate species level for all samples
#smpl_Inv <- env_COI_all[env_COI_all$Niskin.sample %in% colnames(table_unrarefied_Inv),]
smpl_Inv <- env_COI[env_COI$Niskin.sample %in% colnames(table_unrarefied_Inv),]
ps_Inv <- table_unrarefied_Inv[,colnames(table_unrarefied_Inv) %in% smpl_Inv$Niskin.sample]
rownames(smpl_Inv) <- colnames(ps_Inv)
Taxonomy_Inv <- as.matrix(rownames(ps_Inv))
rownames(Taxonomy_Inv) <- rownames(ps_Inv)

ps_unrarefied_Inv <- phyloseq(otu_table(ps_Inv, taxa_are_rows = TRUE), 
                              sample_data(smpl_Inv),tax_table(Taxonomy_Inv))
Coverage_Inv_all <- phyloseq_coverage(physeq = ps_unrarefied_Inv)
Coverage_Inv <- phyloseq_coverage(physeq = ps_unrarefied_Inv)
ps_rarefied_Inv <- phyloseq_coverage_raref(physeq = ps_unrarefied_Inv, 
                                            iter = 1, coverage = 0.7961672, drop_lowcoverage = T)
seqtab_rarefied_Inv <- as.data.frame(ps_rarefied_Inv@otu_table)

Inv_Diversity <- estimateD(ps_Inv, datatype="abundance",
          base="coverage", level=0.985, conf=0.95)

t <- c(1, 21, 41, 61, 81, 101, 150, 200, 250, 
       500, 750, 1000, 1500, 2000, 2500,
       5000, 7500, 10000, 15000, 20000, 25000, 30000)

iNext_Inv <- iNEXT(ps_Inv, q=c(0, 1), datatype="abundance", size = t)
g_Inv <- ggiNEXT(iNext_Inv, type=1, se=TRUE, grey=FALSE)
g_Inv_2 <- g_Inv + 
  scale_colour_manual(values=c(smpl_Inv$Zone_color)) +
  scale_fill_manual(values=c(smpl_Inv$Zone_color)) +
  scale_shape_manual(values=c(rep(20, times = 122))) +
  theme(legend.text = element_text(size = 9.5))
g_Inv_2

ps_Inv_Top3 <- ps_Inv[,colnames(ps_Inv) %in% c(env_COI$Niskin.sample)]
ps_unrarefied_Inv_Top3 <- phyloseq(otu_table(ps_Inv_Top3, 
                                        taxa_are_rows = TRUE))
Coverage_Inv_Top3 <- phyloseq_coverage(physeq = ps_unrarefied_Inv_Top3)


##Observed diversity and Shannon diversity index for all datasets seperatly
#eDNA-12S

#plot_12S <- plot_richness(ps_12S, x="Area", measures=c("Observed", "Shannon"), color="Area") + 
  #geom_boxplot(outlier.shape = NA)
#plot_12S

#eDNA-COI
ps_COI <- merged_data_unrarefied_Animalia
smpl_COI <- env_unrarefied_COI
#smpl_COI$Area <- paste(smpl_COI$Environment, smpl_COI$Zones, sep="_")
#smpl_COI$Area <- ifelse(smpl_COI$Area=="Coastal_Coastal","Coastal", smpl_COI$Area)
smpl_COI$Area <- ifelse(smpl_COI$Zones=="Coastal","Coast", 
                        ifelse(smpl_COI$Zones=="zone1", "Transition", 
                               ifelse(smpl_COI$Zones=="zone2","Offshore","red")))
rownames(smpl_COI) <- colnames(ps_COI)
Taxonomy_COI <- rownames(merged_data_unrarefied_Animalia)
Taxonomy_COI <- as.matrix(Taxonomy_COI)
rownames(Taxonomy_COI) <- rownames(ps_COI)
#ps_COI <- phyloseq(otu_table(ps_COI, taxa_are_rows = TRUE), sample_data(smpl_COI),tax_table(Taxonomy_COI))
#plot_COI <- plot_richness(ps_COI, x="Area", measures=c("Observed", "Shannon"), color="Area") + 
  #geom_boxplot(outlier.shape = NA)
#plot_COI

#Morphology-Fish
#smpl_morph_Fish <- env_morph
smpl_morph_Fish <- env_morph_EqZone
ps_morph_Fish <- as.data.frame(table_morph_Fish_2[,colnames(table_morph_Fish_2) %in% c(smpl_morph_Fish$Niskin.sample, 24, 3, 6)])
rownames(smpl_morph_Fish) <- colnames(ps_morph_Fish)
Taxonomy_morph_Fish <- rownames(table_morph_Fish_2)
Taxonomy_morph_Fish <- as.matrix(Taxonomy_morph_Fish)
rownames(Taxonomy_morph_Fish) <- rownames(ps_morph_Fish)
ps_morph_Fish_1 <- phyloseq(otu_table(ps_morph_Fish, taxa_are_rows = TRUE), sample_data(smpl_morph_Fish),tax_table(Taxonomy_morph_Fish))
Coverage_morph_Fish <- phyloseq_coverage(physeq = ps_morph_Fish_1)
#plot_richness(ps_morph_Fish, x="Area", measures=c("Observed", "Shannon"), color="Area") + 
  #geom_boxplot(outlier.shape = NA)

iNext_Fish_morph <- iNEXT(ps_morph_Fish, q=c(0,1), datatype="abundance", knot = 40)
g_Fish_morph <- ggiNEXT(iNext_Fish_morph, type=3, se=TRUE, grey=FALSE, facet.var="Order.q")
g_Fish_morph_2 <- g_Fish_morph + 
  scale_colour_manual(values=c(smpl_morph_Fish$Zone_color)) +
  scale_fill_manual(values=c(smpl_morph_Fish$Zone_color)) +
  scale_shape_manual(values=c(rep(20, times = 22))) +
  theme(legend.text = element_text(size = 9.5)) 
g_Fish_morph_2

Fish_Diversity_morph_morph <- as.data.frame(iNext_Fish_morph$iNextEst$coverage_based[
  which(iNext_Fish_morph$iNextEst$coverage_based$Order.q == 0),])
Fish_Diversity_morph_morph <- as.data.frame(Fish_Diversity_morph_morph[Fish_Diversity_morph_morph$SC > 0.9598690,])
Fish_Diversity_morph_morph_2 <- aggregate(Fish_Diversity_morph_morph,
                                    by= list(as.factor(Fish_Diversity_morph_morph[,"Assemblage"])), FUN=min)
#Fish_Diversity_morph_morph_2$Method <- ifelse(Fish_Diversity_morph_morph_2$Assemblage == "ftWG2_5", "Observed", "Rarefaction")

env_morph$Niskin.sample <- str_replace(env_morph$Niskin.sample, "_DYF_NJ2021", "") 
Coast_morph <- c(Fish_Diversity_morph_morph_2$qD[
  Fish_Diversity_morph_morph_2$Assemblage %in% env_morph$Niskin.sample[env_morph$Zones == "Coastal"]])
Transition_morph <- c(Fish_Diversity_morph_morph_2$qD[
  Fish_Diversity_morph_morph_2$Assemblage %in% env_morph$Niskin.sample[env_morph$Zones == "zone1"]])
Offshore_morph <- c(Fish_Diversity_morph_morph_2$qD[
  Fish_Diversity_morph_morph_2$Assemblage %in% env_morph$Niskin.sample[env_morph$Zones == "zone2"]])

boxplot(Coast_morph, Transition_morph, Offshore_morph,
        main="Average number of species per Zone", ylab ="Number of species", 
        names = c("Coast", "Transition", "Offshore"),
        col=c("limegreen", "slateblue" , "darkorange"))

#Morphology-Invertebrates
ps_morph_Inv <- table_morph_Inv_2
smpl_morph_Inv <- env_morph
#smpl_morph_Inv$Area <- paste(smpl_morph_Inv$Environment, smpl_morph_Inv$Zones, sep="_")
#smpl_morph_Inv$Area <- ifelse(smpl_morph_Inv$Area=="Coastal_Coastal","Coastal", smpl_morph_Inv$Area)
smpl_morph_Inv$Area <- ifelse(smpl_morph_Inv$Zones=="Coastal","Coast", 
                        ifelse(smpl_morph_Inv$Zones=="zone1", "Transition", 
                               ifelse(smpl_morph_Inv$Zones=="zone2","Offshore","red")))
rownames(smpl_morph_Inv) <- colnames(ps_morph_Inv)
Taxonomy_morph_Inv <- rownames(table_morph_Inv_2)
Taxonomy_morph_Inv <- as.matrix(Taxonomy_morph_Inv)
rownames(Taxonomy_morph_Inv) <- rownames(ps_morph_Inv)
ps_morph_Inv_1 <- phyloseq(otu_table(ps_morph_Inv, taxa_are_rows = TRUE), sample_data(smpl_morph_Inv),tax_table(Taxonomy_morph_Inv))
#plot_morph_Inv <- plot_richness(ps_morph_Inv, x="Area", measures=c("Observed", "Shannon")) + 
  #geom_boxplot(outlier.shape = NA)
Coverage_morph_Inv <- phyloseq_coverage(physeq = ps_morph_Inv_1)

iNext_Inv_morph <- iNEXT(ps_morph_Inv, q=c(0,1), datatype="abundance", knot = 40)
g_Inv_morph <- ggiNEXT(iNext_Inv_morph, type=1, se=TRUE, grey=FALSE, facet.var="Order.q")
g_Inv_morph_2 <- g_Inv_morph + 
  scale_colour_manual(values=c(smpl_morph_Inv$Zone_color)) +
  scale_fill_manual(values=c(smpl_morph_Inv$Zone_color)) +
  scale_shape_manual(values=c(rep(20, times = 22))) +
  theme(legend.text = element_text(size = 9.5)) 
g_Inv_morph_2


#Morphology-Standerdized
ps_morph_Fish_Std <- as.data.frame(table_morph_Fish_Std_2)
smpl_morph_Fish_Std <- env_morph
smpl_morph_Fish_Std$Area <- ifelse(smpl_morph_Fish_Std$Zones=="Coastal","Coast", 
                              ifelse(smpl_morph_Fish_Std$Zones=="zone1", "Transition", 
                                     ifelse(smpl_morph_Fish_Std$Zones=="zone2","Offshore","red")))
rownames(smpl_morph_Fish_Std) <- colnames(ps_morph_Fish_Std)
Taxonomy_morph_Fish_Std <- rownames(table_morph_Fish_Std_2)
Taxonomy_morph_Fish_Std <- as.matrix(Taxonomy_morph_Fish_Std)
rownames(Taxonomy_morph_Fish_Std) <- rownames(ps_morph_Fish_Std)

ps_morph_Inv_Std <- table_morph_Inv_Std_2
smpl_morph_Inv_Std <- env_morph
smpl_morph_Inv_Std$Area <- ifelse(smpl_morph_Inv_Std$Zones=="Coastal","Coast", 
                              ifelse(smpl_morph_Inv_Std$Zones=="zone1", "Transition", 
                                     ifelse(smpl_morph_Inv_Std$Zones=="zone2","Offshore","red")))
rownames(smpl_morph_Inv_Std) <- colnames(ps_morph_Inv_Std)
Taxonomy_morph_Inv_Std <- rownames(table_morph_Inv_Std_2)
Taxonomy_morph_Inv_Std <- as.matrix(Taxonomy_morph_Inv_Std)
rownames(Taxonomy_morph_Inv_Std) <- rownames(ps_morph_Inv_Std)


##Observed diversity and Shannon diversity index for all Fish data 
rownames(smpl_12S) <- paste(colnames(table_unrarefied_Fish), "eDNA_Fish", sep="_")
rownames(smpl_morph_Fish) <- paste(colnames(table_morph_Fish_2), "Morphology_Fish", sep="_")
rownames(smpl_morph_Fish_Std) <- paste(colnames(table_morph_Fish_Std_2), "Morphology_Standerdized_Fish", sep="_")
#smpl_Fish <- bind_rows(smpl_12S, smpl_morph_Fish, smpl_morph_Fish_Std)
#smpl_Fish$Method <- c(rep("eDNA", nrow(smpl_12S)), 
#                      rep("Morphology", nrow(smpl_morph_Fish)), 
#                      rep("Morphology - Standerdized", nrow(smpl_morph_Fish_Std)))
smpl_Fish <- bind_rows(smpl_12S, smpl_morph_Fish)
smpl_Fish$Method <- c(rep("eDNA", nrow(smpl_12S)), 
                      rep("Morphology", nrow(smpl_morph_Fish)))
smpl_Fish$Organism <- sub(".*_","",rownames(smpl_Fish))

colnames(ps_12S) <- paste(colnames(table_unrarefied_Fish), "eDNA_Fish", sep="_")
colnames(ps_morph_Fish) <- paste(colnames(table_morph_Fish_2), "Morphology_Fish", sep="_")
colnames(ps_morph_Fish_Std) <- paste(colnames(table_morph_Fish_Std_2), "Morphology_Standerdized_Fish", sep="_")

ps_Fish <- merge(ps_12S, ps_morph_Fish, by.x=0, by.y=0, all=T)
#ps_Fish <- merge(ps_Fish, ps_morph_Fish_Std, by.x=1, by.y=0, all=T)
rownames(ps_Fish) <- ps_Fish$Row.names
ps_Fish$Row.names <- NULL
ps_Fish[is.na(ps_Fish)] <- 0
colnames(ps_Fish) <- rownames(smpl_Fish)
Taxonomy_Fish <- as.matrix(rownames(ps_Fish))
rownames(Taxonomy_Fish) <- rownames(ps_Fish)

#ps_Fish <- phyloseq(otu_table(ps_Fish, taxa_are_rows = TRUE), sample_data(smpl_Fish),tax_table(Taxonomy_Fish))
#p <- plot_richness(ps_Fish, x="Area", measures=c("Observed", "Shannon")) + 
  #geom_boxplot(outlier.shape = 16, outlier.size = 2, aes(fill=Method)) + 
  #scale_fill_manual(values=c("lightblue3", "darkolivegreen3")) +
  #theme(axis.text.x=element_text(color = c("limegreen","slateblue1","darkorange", "slateblue", "darkorange3")))
#p$layers <- p$layers[-1]
#p

##Observed diversity and Shannon diversity index for all Invertebrate data 
rownames(smpl_COI) <- paste(colnames(merged_data_unrarefied_Animalia), "eDNA_Invertebrates", sep="_")
rownames(smpl_morph_Inv) <- paste(colnames(table_morph_Inv_2), "Morphology_Invertebrates", sep="_")
rownames(smpl_morph_Inv_Std) <- paste(colnames(table_morph_Inv_Std_2), "Morphology_Standerdized_Invertebrates", sep="_")
#smpl_Inv <- bind_rows(smpl_COI, smpl_morph_Inv, smpl_morph_Inv_Std)
#smpl_Inv$Method <- c(rep("eDNA", nrow(smpl_COI)), 
#                     rep("Morphology", nrow(smpl_morph_Inv)),
#                     rep("Morphology - Standerdized", nrow(smpl_morph_Inv_Std)))
smpl_Inv <- bind_rows(smpl_COI, smpl_morph_Inv)
smpl_Inv$Organism <- sub(".*_","",rownames(smpl_Inv))
smpl_Inv$Method <- c(rep("eDNA", nrow(smpl_COI)), 
                     rep("Morphology", nrow(smpl_morph_Inv)))

colnames(ps_COI) <- paste(colnames(merged_data_unrarefied_Animalia), "eDNA_Invertebrates", sep="_")
colnames(ps_morph_Inv) <- paste(colnames(table_morph_Inv_2), "morph_Invertebrates", sep="_")
colnames(ps_morph_Inv_Std) <- paste(colnames(table_morph_Inv_Std_2), "morph_Standerdized_Invertebrates", sep="_")

#ps_Inv <- bind_rows(ps_COI, ps_morph_Inv)
ps_Inv <- merge(ps_COI, ps_morph_Inv, by.x=0, by.y=0, all=T)
#ps_Inv <- merge(ps_Inv, ps_morph_Inv_Std, by.x=1, by.y=0, all=T)
rownames(ps_Inv) <- ps_Inv$Row.names
ps_Inv$Row.names <- NULL
ps_Inv[is.na(ps_Inv)] <- 0
colnames(ps_Inv) <- rownames(smpl_Inv)
Taxonomy_Inv <- as.matrix(rownames(ps_Inv))
rownames(Taxonomy_Inv) <- rownames(ps_Inv)

#ps_Inv <- phyloseq(otu_table(ps_Inv , taxa_are_rows = TRUE), sample_data(smpl_Inv ),tax_table(Taxonomy_Inv))
#p_Inv <- plot_richness(ps_Inv, x="Area", measures=c("Observed", "Shannon")) + 
  #geom_boxplot(outlier.shape = 16, outlier.size = 2, aes(fill=Method)) + 
  #scale_fill_manual(values=c("lightblue3", "darkolivegreen3")) +
  #theme(axis.text.x=element_text(color = c("limegreen","slateblue1","darkorange", "slateblue", "darkorange3")))
#p_Inv$layers <- p_Inv$layers[-1]
#p_Inv

##Observed diversity and Shannon diversity index for all data
smpl_all <- bind_rows(smpl_Fish, smpl_Inv)

ps_all <- merge(ps_Fish, ps_Inv, by.x=0, by.y=0, all=T)
rownames(ps_all) <- ps_all$Row.names
ps_all$Row.names <- NULL
#ps_all <- bind_rows(ps_Fish, ps_Inv)
ps_all[is.na(ps_all)] <- 0
ps_all <- ps_all[,rownames(smpl_all)]
#colnames(ps_all) <- rownames(smpl_all)
Taxonomy_all <- as.matrix(rownames(ps_all))
rownames(Taxonomy_all) <- rownames(ps_all)

ps_all_phylo <- phyloseq(otu_table(ps_all , taxa_are_rows = TRUE), sample_data(smpl_all),tax_table(Taxonomy_all))
p_all <- plot_richness(ps_all_phylo, x='Area', measures=c("Observed", "Shannon")) + 
  geom_boxplot(outlier.shape = 16, outlier.size = 2, aes(fill=Method)) + 
  scale_fill_manual(values=c("darkolivegreen3", "lightblue3")) +
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  facet_grid(variable~Organism, scales = "free")
p_all$layers <- p_all$layers[-1]
p_all$data$Area <- factor(p_all$data$Area, levels=unique(smpl_12S$Area))
p_all

##Statistical analysis - Fish
#Alpha diversity
library(car)
library(DescTools)
library(lsr)
DNA_species <- smpl_Fish
NumOfSp_Fish <- data.frame(colSums(ps_Fish>0))
DNA_species$NumOfSp <- colSums(ps_Fish>0)
NumOfSp_Fish$Area <- DNA_species$Area
colnames(NumOfSp_Fish) <- c("NumOfSp", "Area")

model_DNA <- lm(NumOfSp ~ Area, data=DNA_species) 
model_DNA2 <- lm(NumOfSp ~ Method, data=DNA_species) 
model_DNA3 <- lm(NumOfSp ~ Area*Method, data=DNA_species)
anova_DNA <- Anova(model_DNA, type=2)# SIGN for station , not sign for biological replicates
summary(model_DNA)
anova_DNA
anova_DNA2 <- Anova(model_DNA2, type=2)# SIGN for station , not sign for biological replicates
summary(model_DNA2)
anova_DNA2
anova_DNA3 <- Anova(model_DNA3, type=2)# SIGN for station , not sign for biological replicates
summary(model_DNA3)
anova_DNA3
PostHocTest(aov(NumOfSp ~ Area, data=DNA_species), method = "hsd")
PostHocTest(aov(NumOfSp ~ Method, data=DNA_species), method = "hsd")
PostHocTest(aov(NumOfSp ~ Area*Method, data=DNA_species), method = "hsd")

es <- etaSquared(model_DNA3, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) # 1

library(lsmeans)
library(multcomp)
cld(lsmeans(model_DNA, ~ Area), Letters=letters)
cld(lsmeans(model_DNA2, ~ Method), Letters=letters)
cld(lsmeans(model_DNA3, ~ Area*Method), Letters=letters)
#elk station in ander groep!
#lstrends(model_DNA, "Station", var="NR") #lstrends: moet continue zijn + volgorde belangrijk => DUS HIER NIET GEBRUIKEN
#DNA_species_2$NR <- as.numeric(as.factor(DNA_species_2$Biologic_replicate))

#check assumptions
#homogeneity
plot(model_DNA, 1)
plot(model_DNA2, 1)
plot(model_DNA3, 1)
#leveneTest(NumOfSp ~ Station, data=DNA_species) # SIGNIFICANT (but only on edges)
leveneTest(NumOfSp ~ Area, data=DNA_species) # SIGNIFICANT (but only on edges)
leveneTest(NumOfSp ~ Method, data=DNA_species)
leveneTest(NumOfSp ~ Area*Method, data=DNA_species)

#normality
resid <- residuals(model_DNA3)  # pull the residuals
hist(resid)
qqnorm(resid) 
plot(model_DNA, 2)
plot(model_DNA2, 2)
plot(model_DNA3, 2)
ggqqplot(DNA_species$NumOfSp)
shapiro.test(DNA_species$NumOfSp)  # NOT SIGNIFICANT
shapiro.test(resid)
shapiro.test(1/(DNA_species$NumOfSp))

library(MASS)
boxcox(model_DNA)
bc<-boxcox(model_DNA)
bc$x[which(bc$y==max(bc$y))]

#Shannon diversity
DNA_Shannon<- smpl_Fish
Shannon_Fish <- data.frame(diversity(ps_Fish, "shannon"))
Shannon_Fish$Area <- DNA_Shannon$Area
colnames(Shannon_Fish) <- c("Shannon", "Area")
DNA_Shannon$Shannon <- Shannon_Fish$Shannon

model_Shannon <- lm(Shannon ~ Area, data=DNA_Shannon)
model_Shannon2 <- lm(Shannon ~ Method, data=DNA_Shannon)
model_Shannon3 <- lm(Shannon ~ Area*Method, data=DNA_Shannon)
anova_Shannon <- Anova(model_Shannon, type=2)# SIGN for station , not sign for biological replicates
summary(model_Shannon)
anova_Shannon
anova_Shannon2 <- Anova(model_Shannon2, type=2)# SIGN for station , not sign for biological replicates
summary(model_Shannon2)
anova_Shannon2
anova_Shannon3 <- Anova(model_Shannon3, type=2)# SIGN for station , not sign for biological replicates
summary(model_Shannon3)
anova_Shannon3
PostHocTest(aov(Shannon ~ Area, data=DNA_Shannon), method = "hsd")
PostHocTest(aov(Shannon ~ Method, data=DNA_Shannon), method = "hsd")
PostHocTest(aov(Shannon ~ Area*Method, data=DNA_Shannon), method = "hsd")

es_Shannon <- etaSquared(model_Shannon3, type=2, anova=TRUE)
es_Shannon
sum(es_Shannon[,"eta.sq"]) # 1

cld(lsmeans(model_Shannon, ~ Area), Letters=letters)
cld(lsmeans(model_Shannon2, ~ Method), Letters=letters)
cld(lsmeans(model_Shannon3, ~ Area*Method), Letters=letters) #elk station in ander groep!

#check assumptions
#homogeneity
plot(model_Shannon3, 1)
leveneTest(Shannon ~ Area, data=DNA_Shannon) # SIGNIFICANT (but only on edges)
leveneTest(Shannon ~ Method, data=DNA_Shannon)
leveneTest(Shannon ~ Area*Method, data=DNA_Shannon)

#normality
resid_Shannon <- residuals(model_Shannon3)  # pull the residuals
hist(resid_Shannon)
qqnorm(resid_Shannon) 
plot(model_Shannon, 2)
plot(model_Shannon2, 2)
plot(model_Shannon3, 2)
ggqqplot(DNA_species$NumOfSp)
shapiro.test(DNA_species$NumOfSp)  # NOT SIGNIFICANT
shapiro.test(resid_Shannon)
shapiro.test(1/(DNA_species$NumOfSp))

boxcox(model_Shannon3)
bc<-boxcox(model_Shannon3)
bc$x[which(bc$y==max(bc$y))]

##Statistical analysis - Invertebrates
#Alpha diversity
DNA_species <- smpl_Inv
NumOfSp_Inv <- data.frame(colSums(ps_Inv>0))
DNA_species$NumOfSp <- colSums(ps_Inv>0)
NumOfSp_Inv$Area <- DNA_species$Area
colnames(NumOfSp_Inv) <- c("NumOfSp", "Area")

model_DNA <- lm(NumOfSp ~ Area, data=DNA_species) 
model_DNA2 <- lm(NumOfSp ~ Method, data=DNA_species) 
model_DNA3 <- lm(NumOfSp ~ Area*Method, data=DNA_species)
anova_DNA <- Anova(model_DNA, type=2)# SIGN for station , not sign for biological replicates
summary(model_DNA)
anova_DNA
anova_DNA2 <- Anova(model_DNA2, type=2)# SIGN for station , not sign for biological replicates
summary(model_DNA2)
anova_DNA2
anova_DNA3 <- Anova(model_DNA3, type=2)# SIGN for station , not sign for biological replicates
summary(model_DNA3)
anova_DNA3
PostHocTest(aov(NumOfSp ~ Area, data=DNA_species), method = "hsd")
PostHocTest(aov(NumOfSp ~ Method, data=DNA_species), method = "hsd")
PostHocTest(aov(NumOfSp ~ Area*Method, data=DNA_species), method = "hsd")

es <- etaSquared(model_DNA3, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) # 1

cld(lsmeans(model_DNA, ~ Area), Letters=letters)
cld(lsmeans(model_DNA2, ~ Method), Letters=letters)
cld(lsmeans(model_DNA3, ~ Area*Method), Letters=letters)
#elk station in ander groep!
#lstrends(model_DNA, "Station", var="NR") #lstrends: moet continue zijn + volgorde belangrijk => DUS HIER NIET GEBRUIKEN
#DNA_species_2$NR <- as.numeric(as.factor(DNA_species_2$Biologic_replicate))

#check assumptions
#homogeneity
plot(model_DNA, 1)
plot(model_DNA2, 1)
plot(model_DNA3, 1)
#leveneTest(NumOfSp ~ Station, data=DNA_species) # SIGNIFICANT (but only on edges)
leveneTest(NumOfSp ~ Area, data=DNA_species) # SIGNIFICANT (but only on edges)
leveneTest(NumOfSp ~ Method, data=DNA_species)
leveneTest(NumOfSp ~ Area*Method, data=DNA_species)

#normality
resid <- residuals(model_DNA3)  # pull the residuals
hist(resid)
qqnorm(resid) 
plot(model_DNA, 2)
plot(model_DNA2, 2)
plot(model_DNA3, 2)
ggqqplot(DNA_species$NumOfSp)
shapiro.test(DNA_species$NumOfSp)  # NOT SIGNIFICANT
shapiro.test(resid)
shapiro.test(1/(DNA_species$NumOfSp))

library(MASS)
boxcox(model_DNA)
bc<-boxcox(model_DNA)
bc$x[which(bc$y==max(bc$y))]

#Shannon diversity
DNA_Shannon <- smpl_Inv
Shannon <- data.frame(diversity(ps_Inv, "shannon"))
Shannon$Area <- DNA_Shannon$Area
colnames(Shannon) <- c("Shannon", "Area")
DNA_Shannon$Shannon <- Shannon$Shannon

model_Shannon <- lm(Shannon ~ Area, data=DNA_Shannon)
model_Shannon2 <- lm(Shannon ~ Method, data=DNA_Shannon)
model_Shannon3 <- lm(Shannon ~ Area*Method, data=DNA_Shannon)
anova_Shannon <- Anova(model_Shannon, type=2)# SIGN for station , not sign for biological replicates
summary(model_Shannon)
anova_Shannon
anova_Shannon2 <- Anova(model_Shannon2, type=2)# SIGN for station , not sign for biological replicates
summary(model_Shannon2)
anova_Shannon2
anova_Shannon3 <- Anova(model_Shannon3, type=2)# SIGN for station , not sign for biological replicates
summary(model_Shannon3)
anova_Shannon3
PostHocTest(aov(Shannon ~ Area, data=DNA_Shannon), method = "hsd")
PostHocTest(aov(Shannon ~ Method, data=DNA_Shannon), method = "hsd")
PostHocTest(aov(Shannon ~ Area*Method, data=DNA_Shannon), method = "hsd")

es_Shannon <- etaSquared(model_Shannon3, type=2, anova=TRUE)
es_Shannon
sum(es_Shannon[,"eta.sq"]) # 0.9924824: for only a small fraction of the variance we don’t know which of the effects in the model is responsible.

cld(lsmeans(model_Shannon, ~ Area), Letters=letters)
cld(lsmeans(model_Shannon2, ~ Method), Letters=letters)
cld(lsmeans(model_Shannon3, ~ Area*Method), Letters=letters) #elk station in ander groep!

#check assumptions
#homogeneity
plot(model_Shannon3, 1)
leveneTest(Shannon ~ Area, data=DNA_Shannon) # SIGNIFICANT (but only on edges)
leveneTest(Shannon ~ Method, data=DNA_Shannon)
leveneTest(Shannon ~ Area*Method, data=DNA_Shannon)

#normality
resid_Shannon <- residuals(model_Shannon3)  # pull the residuals
hist(resid_Shannon)
qqnorm(resid_Shannon) 
plot(model_Shannon, 2)
plot(model_Shannon2, 2)
plot(model_Shannon3, 2)
ggqqplot(DNA_species$NumOfSp)
shapiro.test(DNA_species$NumOfSp)  # NOT SIGNIFICANT
shapiro.test(resid_Shannon)
shapiro.test(1/(DNA_species$NumOfSp))

boxcox(model_Shannon3)
bc<-boxcox(model_Shannon3)
bc$x[which(bc$y==max(bc$y))]

## Calculate relative read abundances assigned
Clytia_hemisphaerica <- sum(merged_data_unrarefied_Animalia[rownames(merged_data_unrarefied_Animalia) %in% "Clytia hemisphaerica",])
Total_reads <- sum(merged_data_unrarefied_Animalia)

Whiting <- sum(table_unrarefied_Fish[rownames(table_unrarefied_Fish) %in% "Merlangius merlangus",])
Total_reads <- sum(table_unrarefied_Fish)

##Additional tests for coverage calculation based on ASV level instead of species level
#Calculate the Coverage of FishASV level
ps_FishASVs <- table_unrarefied_FishASVs_noT
smpl_FishASVs <- env_12S[env_12S$Niskin.sample %in% colnames(ps_FishASVs),]
rownames(smpl_FishASVs) <- colnames(ps_FishASVs)
Taxonomy_FishASVs <- as.matrix(table_unrarefied_FishASVs[,(ncol(table_unrarefied_FishASVs)-10):(ncol(table_unrarefied_FishASVs)-4)])

ps_unrarefied_FishASVs <- phyloseq(otu_table(ps_FishASVs, taxa_are_rows = TRUE),
                                   sample_data(smpl_FishASVs),
                                   tax_table(Taxonomy_FishASVs))
Coverage_FishASVs <- phyloseq_coverage(physeq = ps_unrarefied_FishASVs)
ps_rarefied_FishASVs <- phyloseq_coverage_raref(physeq = ps_unrarefied_FishASVs, 
                                                iter = 1, coverage = 0.99, drop_lowcoverage = T)

seqtab_unrarefied_FishASVs <- as.data.frame(t(table_unrarefied_FishASVs_noT))
seqtab_rarefied_FishASVs <- as.data.frame(ps_rarefied_FishASVs@otu_table)

rarecurve(seqtab_unrarefied_FishASVs, ylab = "Fish ASVs", 
          main = "Rarecurve of unrarefied samples after taxonomic assignment", 
          col = as.vector(env_12S$Environment_color), label = FALSE, step =100)
rarecurve(seqtab_rarefied_FishASVs, ylab = "Fish ASVs", 
          main = "Rarecurve of unrarefied samples after taxonomic assignment", 
          col = as.vector(smpl_12S$Environment_color), label = FALSE, step =100)

plot_FishASVs <- plot_richness(ps_rarefied_FishASVs, x="Zone", measures=c("Observed", "Shannon"), color="Zone") + 
  geom_boxplot(outlier.shape = NA)
plot_FishASVs

iNext_FishASVs <- iNEXT(ps_FishASVs, q=c(0,1,2), datatype="abundance")
ggiNEXT(iNext_FishASVs, type=1, se=TRUE, grey=FALSE)


#Calculate the Coverage of Invertebrate ASV level
ps_InvASVs <- table_unrarefied_InvASVs_noT
smpl_InvASVs <- env_COI[env_COI$Niskin.sample %in% colnames(ps_InvASVs),]
rownames(smpl_InvASVs) <- colnames(ps_InvASVs)
Taxonomy_InvASVs <- as.matrix(table_unrarefied_InvASVs[,(ncol(table_unrarefied_InvASVs)-10):(ncol(table_unrarefied_InvASVs)-4)])

ps_unrarefied_InvASVs <- phyloseq(otu_table(ps_InvASVs, taxa_are_rows = TRUE), sample_data(smpl_InvASVs),tax_table(Taxonomy_InvASVs))
Coverage_InvASVs <- phyloseq_coverage(physeq = ps_unrarefied_InvASVs)
ps_rarefied_InvASVs <- phyloseq_coverage_raref(physeq = ps_unrarefied_InvASVs, 
                                               iter = 1, coverage = 0.85, drop_lowcoverage = T)

seqtab_rarefied_InvASVs <- as.data.frame(ps_rarefied_InvASVs@otu_table)
plot_InvASVs <- plot_richness(ps_rarefied_InvASVs, x="Zone", measures=c("Observed", "Shannon"), color="Zone") + 
  geom_boxplot(outlier.shape = NA)
plot_InvASVs

rarecurve(t(table_unrarefied_InvASVs_noT), ylab = "Inv ASVs", 
          main = "Rarecurve of unrarefied samples after taxonomic assignment", 
          col = as.vector(env_COI$Environment_color), label = FALSE, step =100)
rarecurve(seqtab_rarefied_InvASVs, ylab = "Inv ASVs", 
          main = "Rarecurve of unrarefied samples after taxonomic assignment", 
          col = as.vector(smpl_COI$Environment_color), label = FALSE, step =100)

##Additional tests for coverage based rarefaction per Zone
#sum up the reads per Zone
seqtab_unrarefied_Fish <- as.data.frame(t(ps_Fish))
seqtab_unrarefied_Fish$Zone <- smpl_Fish$Zone
seqtab_unrarefied_Fish <- aggregate(seqtab_unrarefied_Fish[,1:(ncol(seqtab_unrarefied_Fish)-1)],
                                    by= list(as.factor(seqtab_unrarefied_Fish[,"Zone"])),FUN=mean)
rownames(seqtab_unrarefied_Fish) <- seqtab_unrarefied_Fish$Group.1
seqtab_unrarefied_Fish$Group.1 <- NULL
for (r in 1:nrow(seqtab_unrarefied_Fish)) {
  for (c in 1:ncol(seqtab_unrarefied_Fish)){
    seqtab_unrarefied_Fish[r,c] <- as.numeric(format(round(c(seqtab_unrarefied_Fish[r,c]), 0), nsmall = 0))
  }
}

iNext_Fish_Zone <- iNEXT(t(seqtab_unrarefied_Fish), q=c(0,1), datatype="abundance", size = t) #knot = 40
ggiNEXT(iNext_Fish_Zone, type=3, se=TRUE, grey=FALSE, facet.var = "Order.q")

ps_Fish_Zone <- phyloseq(otu_table(t(seqtab_unrarefied_Fish), taxa_are_rows = TRUE))
Coverage_Fish_Zone <- phyloseq_coverage(physeq = ps_Fish_Zone)

seqtab_unrarefied_Inv <- as.data.frame(t(ps_Inv))
seqtab_unrarefied_Inv$Zone <- smpl_Inv$Zone
seqtab_unrarefied_Inv <- aggregate(seqtab_unrarefied_Inv[,1:(ncol(seqtab_unrarefied_Inv)-1)], 
                                   by= list(as.factor(seqtab_unrarefied_Inv[,"Zone"])),FUN=sum)
rownames(seqtab_unrarefied_Inv) <- seqtab_unrarefied_Inv$Group.1
seqtab_unrarefied_Inv$Group.1 <- NULL
iNext_Inv_Zone <- iNEXT(t(seqtab_unrarefied_Inv), q=c(0), datatype="abundance")
ggiNEXT(iNext_Inv_Zone, type=1, se=TRUE, grey=FALSE)

ps_Inv_Zone <- phyloseq(otu_table(t(seqtab_unrarefied_Inv), taxa_are_rows = TRUE))
Coverage_Inv_Zone <- phyloseq_coverage(physeq = ps_Inv_Zone)

seqtab_morph_Fish <- as.data.frame(t(ps_morph_Fish[,colnames(ps_morph_Fish) %in% env_morph$Niskin.sample]))
seqtab_morph_Fish$Zone <- smpl_morph_Fish$Zone
seqtab_morph_Fish <- aggregate(seqtab_morph_Fish[,1:(ncol(seqtab_morph_Fish)-1)], 
                                   by= list(as.factor(seqtab_morph_Fish[,"Zone"])),FUN=mean)
rownames(seqtab_morph_Fish) <- seqtab_morph_Fish$Group.1
seqtab_morph_Fish$Group.1 <- NULL
for (r in 1:nrow(seqtab_morph_Fish)) {
  for (c in 1:ncol(seqtab_morph_Fish)){
    seqtab_morph_Fish[r,c] <- as.numeric(format(round(c(seqtab_morph_Fish[r,c]), 0), nsmall = 0))
  }
}

iNext_Fish_morph_Zone <- iNEXT(t(seqtab_morph_Fish), q=c(0,1), datatype="abundance", knot = 40) 
ggiNEXT(iNext_Fish_morph_Zone, type=1, se=TRUE, grey=FALSE, facet.var="Order.q")


seqtab_morph_Inv <- as.data.frame(t(ps_morph_Inv))
seqtab_morph_Inv$Zone <- smpl_morph_Inv$Zones
seqtab_morph_Inv <- aggregate(seqtab_morph_Inv[,1:(ncol(seqtab_morph_Inv)-1)], 
                               by= list(as.factor(seqtab_morph_Inv[,"Zone"])),FUN=mean)
rownames(seqtab_morph_Inv) <- seqtab_morph_Inv$Group.1
seqtab_morph_Inv$Group.1 <- NULL
for (r in 1:nrow(seqtab_morph_Inv)) {
  for (c in 1:ncol(seqtab_morph_Inv)){
    seqtab_morph_Inv[r,c] <- as.numeric(format(round(c(seqtab_morph_Inv[r,c]), 0), nsmall = 0))
  }
}

iNext_Inv_morph_Zone <- iNEXT(t(seqtab_morph_Inv), q=c(0,1), datatype="abundance", knot = 40)
ggiNEXT(iNext_Inv_morph_Zone, type=3, se=TRUE, grey=FALSE, facet.var="Order.q")



#keep the samples in the three zones separate
install.packages("remotes")
remotes::install_github("taowenmicro/ggClusterNet")
ps_Fish_Coast <- table_unrarefied_Fish[,
  colnames(table_unrarefied_Fish) %in% env_12S$Niskin.sample[
  which(env_12S$Zone == "Coast")]]
ps_Fish_Coast <- phyloseq(otu_table(ps_Fish_Coast, taxa_are_rows = TRUE))
veganComm_Fish_Coast <- vegan_otu(ps_Fish_Coast)
ps_Fish_Transition <- table_unrarefied_Fish[,
  colnames(table_unrarefied_Fish) %in% env_12S$Niskin.sample[
  which(env_12S$Zone == "Transition")]]
ps_Fish_Offshore <- table_unrarefied_Fish[,
  colnames(table_unrarefied_Fish) %in% env_12S$Niskin.sample[
  which(env_12S$Zone == "Offshore")]]

together_Fish <- list(t(ps_Fish_Coast), t(ps_Fish_Transition), t(ps_Fish_Offshore)) 
names(together_Fish)<-c("Coast","Transition","Offshore") # add the names of each habitat
str(together_Fish)

## ready for import!
output_Fish <- iNEXT(together_Fish, q=0, datatype="abundance", se=TRUE, knot = 20)
output_Fish

plot(output_Fish, type=1, se=T, show.legend=T, 
     title("a. Fish diversity - 12S",adj=0,line=1.5), 
     col= c("limegreen","slateblue","darkorange"))
plot(output_Fish, type=2, se=T, show.legend=T, 
     show.main=F,
     col= c("limegreen","slateblue","darkorange"))
plot(output_Fish, type=3, se=T, show.legend=T, 
     show.main=F, 
     col= c("limegreen","slateblue","darkorange"))

