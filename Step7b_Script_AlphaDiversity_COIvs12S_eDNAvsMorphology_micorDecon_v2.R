#### Isolde Cornelis
#### Alpha Diversity for eDNA metabarcoding data using coverage based rarefaction
#### 02/02/2024

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
               , "ggClusterNet"
)

for (opties in libraries){
  
  if (opties %in% installed.packages()){
    
    library(opties,character.only = TRUE)
    
  } else {install.packages(opties,repos = "http://cran.us.r-project.org")
    
    library(opties,character.only = TRUE)
  }
}

###make paths
proj.path.12S <- here("/home/genomics/icornelis/02_ZEROimpact/01_12S/NJ2021/MiFish-UE_run2")
proj.path.COI <- here("/home/genomics/icornelis/02_ZEROimpact/02_COI/NJ2021")

###upload data
table_12S <- readxl::read_excel(paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean.xlsx"))
table_morph_Fish <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Standerdized.xlsx"),sheet = "Fish - Standerdized")

table_COI <- readxl::read_excel(paste0(proj.path.COI,"/OWFvsCoastal_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean.xlsx"))
table_morph_Inv <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Standerdized.xlsx"),sheet = "Epi - Standerdized")

env_12S_all <- readRDS(paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/R_Environment/env_ordered_noNeg.rds"))
env_COI_all <- readRDS(paste0(proj.path.COI,"/OWFvsCoastal_concatenated/results_microDecon/R_Environment/env_ordered_noNeg.rds"))
env_12S <- read.csv(paste0(proj.path.12S,"/Step5_Statistics/environmental_data_OWF3HighestReadNumbers_LocationsTrawl.csv"),  header=TRUE, sep=";")
env_COI <- read.csv(paste0(proj.path.COI,"/Step5_Statistics/environmental_data_OWF3HighestReadNumbers_LocationsTrawl.csv"),  header=TRUE, sep=";")
env_morph <- read.csv(paste0(proj.path.12S,"/Step5_Statistics/environmental_data_morph.csv"),  header=TRUE, sep=";")

###add color and pch to environmental data to create the plot 
env_morph$Zone_color <- ifelse(env_morph$Zone=="Coast","limegreen", 
                               ifelse(env_morph$Zone=="Transition", "slateblue", 
                                       ifelse(env_morph$Zone=="Offshore","darkorange","red")))

env_12S$Environment_color <- env_12S$Zone_color
env_12S$Zone_color <- ifelse(env_12S$Zone=="Coast","limegreen", 
                               ifelse(env_12S$Zone=="Transition", "slateblue", 
                                      ifelse(env_12S$Zone=="Offshore","darkorange","red")))
env_COI$Environment_color <- env_COI$Zone_color
env_COI$Zone_color <- ifelse(env_COI$Zone=="Coast","limegreen", 
                             ifelse(env_COI$Zone=="Transition", "slateblue", 
                                    ifelse(env_COI$Zone=="Offshore","darkorange","red")))


###Prepare the morphological data
table_morph_Fish_2 <- as.data.frame(table_morph_Fish[,2:ncol(table_morph_Fish)])
table_morph_Fish_2 [is.na(table_morph_Fish_2 )] <- 0
rownames(table_morph_Fish_2) <- as.character(table_morph_Fish$...1)
table_morph_Fish_2 <- round(table_morph_Fish_2, digits = 0)
table_morph_Inv_2 <- as.data.frame(table_morph_Inv[,2:ncol(table_morph_Inv)])
table_morph_Inv_2 [is.na(table_morph_Inv_2)] <- 0
rownames(table_morph_Inv_2) <- as.character(table_morph_Inv$...1)
table_morph_Inv_2 <- round(table_morph_Inv_2, digits = 0)

###Prepare the 12S eDNA metabarcoding data
##select marine Fish species
fish_classes <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_classes.rds"))
freshwater_fish <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_Freshwater.rds"))
table_unrarefied_FishASVs <- as.data.frame(table_12S[table_12S$Class %in% fish_classes,])
table_unrarefied_FishASVs <- as.data.frame(table_unrarefied_FishASVs[
  !table_unrarefied_FishASVs$Species %in% c(freshwater_fish, "NA"),])
rownames(table_unrarefied_FishASVs) <- table_unrarefied_FishASVs$ASV
table_unrarefied_FishASVs[is.na(table_unrarefied_FishASVs)] <- 0

##merge the reads according on fish species level
taxo <- "Species"
table_unrarefied_Fish <- aggregate(table_unrarefied_FishASVs[,1:(ncol(table_unrarefied_FishASVs)-11)], 
                                   by= list(as.factor(table_unrarefied_FishASVs[,taxo])),FUN=sum)
rownames(table_unrarefied_Fish) <- as.character(table_unrarefied_Fish$Group.1)
table_unrarefied_Fish$Group.1 <- NULL
table_unrarefied_Fish <- table_unrarefied_Fish[!rowSums(table_unrarefied_Fish) == 0,]
table_unrarefied_Fish <- table_unrarefied_Fish[,!colSums(table_unrarefied_Fish) == 0]

###Prepare the COI eDNA metabarcoding data
##select Invertebrate species
table_unrarefied_InvASVs <- as.data.frame(table_COI[table_COI$Kingdom %in% c("Animalia"),])
table_unrarefied_InvASVs <- as.data.frame(table_unrarefied_InvASVs[!table_unrarefied_InvASVs$Species %in% "NA",])
table_unrarefied_InvASVs <- as.data.frame(table_unrarefied_InvASVs[!table_unrarefied_InvASVs$Phylum %in% "Chordata",])
table_unrarefied_InvASVs[is.na(table_unrarefied_InvASVs)] <- 0

##merge the reads according on invertebrate species level
taxo <- "Species"
table_unrarefied_Inv <- aggregate(table_unrarefied_InvASVs[,1:(ncol(table_unrarefied_InvASVs)-11)],
                                  by= list(as.factor(table_unrarefied_InvASVs[,taxo])),FUN=sum)
rownames(table_unrarefied_Inv) <-as.character(table_unrarefied_Inv$Group.1)
table_unrarefied_Inv$Group.1 <- NULL
table_unrarefied_Inv <- table_unrarefied_Inv[!rowSums(table_unrarefied_Inv) == 0,]
table_unrarefied_Inv <- table_unrarefied_Inv[,!colSums(table_unrarefied_Inv) == 0]

###Rarefy the 12S eDNA metabarcoding data
## Create a phyloseq object
ps_eDNA_Fish <- table_unrarefied_Fish[,colnames(table_unrarefied_Fish) %in% env_12S$Niskin.sample]
smpl_eDNA_Fish <- env_12S[env_12S$Niskin.sample  %in%  colnames(ps_eDNA_Fish),]
rownames(smpl_eDNA_Fish) <- colnames(ps_eDNA_Fish)
Taxonomy_eDNA_Fish <- as.matrix(rownames(ps_eDNA_Fish))
rownames(Taxonomy_eDNA_Fish) <- rownames(ps_eDNA_Fish)
ps_unrarefied_Fish <- phyloseq(otu_table(ps_eDNA_Fish, taxa_are_rows = TRUE),
                               sample_data(smpl_eDNA_Fish),
                               tax_table(Taxonomy_eDNA_Fish))

##Calculate the Coverage
Coverage_eDNA_Fish <- phyloseq_coverage(physeq = ps_unrarefied_Fish)

##Estimate the diversity based on the minimum coverage with a value > 0.5
eDNA_Fish_Diversity <- estimateD(ps_eDNA_Fish, datatype="abundance", base="coverage",
                            level=min(Coverage_eDNA_Fish$SampleCoverage[
                              which(Coverage_eDNA_Fish$SampleCoverage > 0.5)] ),
                            conf=0.95) #select the minimum coverage
eDNA_Fish_Diversity_100 <- estimateD(ps_eDNA_Fish, datatype="abundance",
                                base="coverage", level=1, conf=0.95)

##Rarefy the data based on a coverage just below the minimum coverage using the function phyloseq_coverage_raref,
##(using the minimum coverage will remove the sample from the rarefied dataset)
ps_rarefied_Fish <- phyloseq_coverage_raref(physeq = ps_unrarefied_Fish, 
                                           iter = 1, coverage = 0.998, drop_lowcoverage = T)
table_rarefied_Fish <- as.data.frame(ps_rarefied_Fish@otu_table)

##Optional: plot the rarefied and unrarefied data for visualization 
plot_Fish_rarefied <- plot_richness(ps_rarefied_Fish, x="Zone",
                                    measures=c("Observed", "Shannon"),
                                    color="Zone") + 
geom_boxplot(outlier.shape = NA)
plot_Fish_rarefied 

plot_Fish_unrarefied <- plot_richness(ps_unrarefied_Fish, x="Zone",
                                      measures=c("Observed", "Shannon"),
                                      color="Zone") + 
  geom_boxplot(outlier.shape = NA)
plot_Fish_unrarefied 

###Save Table containing the number of reads per species per Locations (only for locations where beam trawl data is available)
seqtab_Fish <- as.data.frame(t(ps_eDNA_Fish))
rownames(seqtab_Fish) <- str_replace(rownames(seqtab_Fish), "_DYF_NJ2021", "") # removes biological replicate from sample name
seqtab_Fish$names <- str_sub(rownames(seqtab_Fish), end=-3) # removes biological replicate from sample name
seqtab_concatenated_Fish <- aggregate(seqtab_Fish[,1:ncol(seqtab_Fish)-1], by= list(seqtab_Fish$names),FUN=sum)
rownames(seqtab_concatenated_Fish) <- seqtab_concatenated_Fish$Group.1
seqtab_concatenated_Fish$Group.1 <- NULL  
write.xlsx(as.data.frame(t(seqtab_concatenated_Fish)), 
           paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/table_Fish_BioRepConcatenates.xlsx"), 
           sheetName = "Fish Per Location", colNames = TRUE, rowNames = TRUE, append = FALSE)

###Rarefy the COI eDNA metabarcoding data
## Create a phyloseq object
ps_eDNA_Inv <- table_unrarefied_Inv[,colnames(table_unrarefied_Inv) %in% env_COI$Niskin.sample]
smpl_eDNA_Inv <- env_COI[env_COI$Niskin.sample %in% colnames(ps_eDNA_Inv),]
ps_eDNA_Inv <- ps_eDNA_Inv[!rowSums(ps_eDNA_Inv) == 0,]
ps_eDNA_Inv <- ps_eDNA_Inv[,!colSums(ps_eDNA_Inv) == 0]
rownames(smpl_eDNA_Inv) <- colnames(ps_eDNA_Inv)
Taxonomy_eDNA_Inv <- as.matrix(rownames(ps_eDNA_Inv))
rownames(Taxonomy_eDNA_Inv) <- rownames(ps_eDNA_Inv)
ps_unrarefied_Inv <- phyloseq(otu_table(ps_eDNA_Inv, taxa_are_rows = TRUE), 
                              sample_data(smpl_eDNA_Inv),tax_table(Taxonomy_eDNA_Inv))

##Calculate the Coverage
Coverage_eDNA_Inv <- phyloseq_coverage(physeq = ps_unrarefied_Inv)

##Estimate the diversity based on the minimum coverage with a value > 0.5
eDNA_Inv_Diversity <- estimateD(ps_eDNA_Inv, datatype="abundance", base="coverage",
                           level=min(Coverage_eDNA_Inv$SampleCoverage[
                             which(Coverage_eDNA_Inv$SampleCoverage > 0.5)] ),
                           conf=0.95) #select the minimum coverage
eDNA_Inv_Diversity_rare <- estimateD(ps_eDNA_Inv, datatype="abundance", base="coverage",
                                level=0.97, conf=0.95) #coverage used for rarefaction, with a lower coverage we get an error
eDNA_Inv_Diversity_100 <- estimateD(ps_eDNA_Inv, datatype="abundance", base="coverage",
                                level=1, conf=0.95)

##Rarefy the data based on a minimum coverage using the function phyloseq_coverage_raref,
##(using the minimum coverage that will not generate an error)
ps_rarefied_Inv <- phyloseq_coverage_raref(physeq = ps_unrarefied_Inv, 
                                           iter = 1, coverage = 0.97, drop_lowcoverage = T)
table_rarefied_Inv <- as.data.frame(ps_rarefied_Inv@otu_table)

##Optional: plot the rarefied and unrarefied data for visualization 
plot_Inv_rarefied <- plot_richness(ps_rarefied_Inv, x="Zone",
                                   measures=c("Observed", "Shannon"),
                                   color="Zone") + 
  geom_boxplot(outlier.shape = NA)
plot_Inv_rarefied 

plot_Inv_unrarefied <- plot_richness(ps_unrarefied_Inv, x="Zone",
                                     measures=c("Observed", "Shannon"),
                                     color="Zone") + 
  geom_boxplot(outlier.shape = NA)
plot_Inv_unrarefied 

###Save Table containing the number of reads per species per Locations (only for locations where beam trawl data is available)
seqtab_Inv <- as.data.frame(t(ps_eDNA_Inv))
rownames(seqtab_Inv) <- str_replace(rownames(seqtab_Inv), "_DYF_NJ2021", "") # removes biological replicate from sample name
seqtab_Inv$names <- str_sub(rownames(seqtab_Inv), end=-3) # removes biological replicate from sample name
seqtab_concatenated_Inv <- aggregate(seqtab_Inv[,1:ncol(seqtab_Inv)-1], by= list(seqtab_Inv$names),FUN=sum)
rownames(seqtab_concatenated_Inv) <- seqtab_concatenated_Inv$Group.1
seqtab_concatenated_Inv$Group.1 <- NULL  
write.xlsx(as.data.frame(t(seqtab_concatenated_Inv)), 
           paste0(proj.path.COI,"/OWFvsCoastal_concatenated/results_microDecon/table_Inv_BioRepConcatenates.xlsx"), 
           sheetName = "Inv Per Location", colNames = TRUE, rowNames = TRUE, append = FALSE)

###Prepare the Morphological data for the fish catch 
## Create a phyloseq object
ps_morph_Fish <- as.data.frame(table_morph_Fish_2)
smpl_morph_Fish_raw <- env_morph
rownames(smpl_morph_Fish_raw) <- colnames(ps_morph_Fish)
Taxonomy_morph_Fish <- as.matrix(rownames(table_morph_Fish_2))
rownames(Taxonomy_morph_Fish) <- rownames(ps_morph_Fish)
ps_morph_Fish_phylo <- phyloseq(otu_table(ps_morph_Fish, taxa_are_rows = TRUE),
                                sample_data(smpl_morph_Fish_raw),
                                tax_table(Taxonomy_morph_Fish))

##Calculate the Coverage
Coverage_morph_Fish <- phyloseq_coverage(physeq = ps_morph_Fish_phylo)

##Rarefy the data based on a minimum coverage using the function phyloseq_coverage_raref,
##(using the minimum coverage that will not generate an error)
ps_morph_rarefied_Fish <- phyloseq_coverage_raref(physeq = ps_morph_Fish_phylo, 
                                                 iter = 1, coverage = 0.998, drop_lowcoverage = T)
table_morph_rarefied_Fish <- as.data.frame(ps_morph_rarefied_Fish@otu_table)


###Prepare the Morphological data for the epibenthos catch 
## Create a phyloseq object
ps_morph_Inv <- table_morph_Inv_2
smpl_morph_Inv_raw <- env_morph
rownames(smpl_morph_Inv_raw) <- colnames(ps_morph_Inv)
Taxonomy_morph_Inv <- as.matrix(rownames(table_morph_Inv_2))
rownames(Taxonomy_morph_Inv) <- rownames(ps_morph_Inv)
ps_morph_Inv_phylo <- phyloseq(otu_table(ps_morph_Inv, taxa_are_rows = TRUE),
                               sample_data(smpl_morph_Inv_raw),
                               tax_table(Taxonomy_morph_Inv))

##Calculate the Coverage
Coverage_morph_Inv <- phyloseq_coverage(physeq = ps_morph_Inv_phylo)

##Rarefy the data based on a minimum coverage using the function phyloseq_coverage_raref,
ps_morph_rarefied_Inv <- phyloseq_coverage_raref(physeq = ps_morph_Inv_phylo, 
                                                  iter = 1, coverage = 0.97, drop_lowcoverage = T)
table_morph_rarefied_Inv <- as.data.frame(ps_morph_rarefied_Inv@otu_table)

###Analyze the observed diversity and Shannon diversity index for all data 
##Prepare the smpl matrix for the combined fish data
smpl_eDNA_Fish_unrarefied <- smpl_eDNA_Fish
smpl_eDNA_Fish_rarefied <- smpl_eDNA_Fish[smpl_eDNA_Fish$Niskin.sample %in% 
                                            colnames(table_rarefied_Fish),]
smpl_morph_Fish <- smpl_morph_Fish_raw

rownames(smpl_eDNA_Fish_unrarefied) <- paste(rownames(smpl_eDNA_Fish_unrarefied), 
                                  "eDNA_unrarefied_Fish", sep="_")
rownames(smpl_eDNA_Fish_rarefied) <- paste(rownames(smpl_eDNA_Fish_rarefied), 
                                  "eDNA_rarefied_Fish", sep="_")
rownames(smpl_morph_Fish) <- paste(rownames(smpl_morph_Fish),
                                   "Morphology_Fish", sep="_")

smpl_Fish <- bind_rows(smpl_eDNA_Fish_unrarefied,
                       smpl_eDNA_Fish_rarefied,
                       smpl_morph_Fish)
smpl_Fish$Method <- c(rep("eDNA_unrarefied", nrow(smpl_eDNA_Fish_unrarefied)), 
                      rep("eDNA_rarefied", nrow(smpl_eDNA_Fish_rarefied)), 
                      rep("Morphology", nrow(smpl_morph_Fish)))
smpl_Fish$Organism <- sub(".*_","",rownames(smpl_Fish))

##Prepare the abundance table for the combined fish data
ps_eDNA_Fish_unrarefied <- table_unrarefied_Fish
ps_eDNA_Fish_rarefied <- table_rarefied_Fish
ps_morph_Fish_2 <- ps_morph_Fish

colnames(ps_eDNA_Fish_unrarefied) <- paste(colnames(table_unrarefied_Fish), 
                                     "eDNA_unrarefied_Fish", sep="_")
colnames(ps_eDNA_Fish_rarefied) <- paste(colnames(table_rarefied_Fish),
                                   "eDNA_rarefied_Fish", sep="_")
colnames(ps_morph_Fish_2) <- paste(colnames(table_morph_Fish_2), "Morphology_Fish", sep="_")

ps_Fish <- merge(ps_eDNA_Fish_unrarefied, ps_eDNA_Fish_rarefied, by.x=0, by.y=0, all=T)
ps_Fish <- merge(ps_Fish, ps_morph_Fish_2, by.x=1, by.y=0, all=T)
rownames(ps_Fish) <- ps_Fish$Row.names
ps_Fish$Row.names <- NULL
ps_Fish[is.na(ps_Fish)] <- 0
Taxonomy_Fish <- as.matrix(rownames(ps_Fish))
rownames(Taxonomy_Fish) <- rownames(ps_Fish)

#ps_Fish <- phyloseq(otu_table(ps_Fish, taxa_are_rows = TRUE), sample_data(smpl_Fish),tax_table(Taxonomy_Fish))
#p <- plot_richness(ps_Fish, x="Area", measures=c("Observed", "Shannon")) + 
  #geom_boxplot(outlier.shape = 16, outlier.size = 2, aes(fill=Method)) + 
  #scale_fill_manual(values=c("lightblue3", "darkolivegreen3")) +
  #theme(axis.text.x=element_text(color = c("limegreen","slateblue1","darkorange", "slateblue", "darkorange3")))
#p$layers <- p$layers[-1]
#p

##Prepare the smpl matrix for the invertebrate fish data
smpl_eDNA_Inv_unrarefied <- smpl_eDNA_Inv
smpl_eDNA_Inv_rarefied <- smpl_eDNA_Inv[smpl_eDNA_Inv$Niskin.sample %in% 
                                            colnames(table_rarefied_Inv),]
smpl_morph_Inv <- smpl_morph_Inv_raw

rownames(smpl_eDNA_Inv_unrarefied) <- paste(rownames(smpl_eDNA_Inv_unrarefied), 
                                             "eDNA_unrarefied_Inv", sep="_")
rownames(smpl_eDNA_Inv_rarefied) <- paste(rownames(smpl_eDNA_Inv_rarefied), 
                                           "eDNA_rarefied_Inv", sep="_")
rownames(smpl_morph_Inv) <- paste(rownames(smpl_morph_Inv),
                                   "Morphology_Inv", sep="_")

smpl_Inv <- bind_rows(smpl_eDNA_Inv_unrarefied,
                       smpl_eDNA_Inv_rarefied,
                       smpl_morph_Inv)
smpl_Inv$Method <- c(rep("eDNA_unrarefied", nrow(smpl_eDNA_Inv_unrarefied)), 
                      rep("eDNA_rarefied", nrow(smpl_eDNA_Inv_rarefied)), 
                      rep("Morphology", nrow(smpl_morph_Inv)))
smpl_Inv$Organism <- "Invertebrates"

##Prepare the abundance table for the combined invertebrate data
ps_eDNA_Inv_unrarefied <- table_unrarefied_Inv
ps_eDNA_Inv_rarefied <- table_rarefied_Inv
ps_morph_Inv_2 <- ps_morph_Inv

colnames(ps_eDNA_Inv_unrarefied) <- paste(colnames(table_unrarefied_Inv), 
                                           "eDNA_unrarefied_Inv", sep="_")
colnames(ps_eDNA_Inv_rarefied) <- paste(colnames(table_rarefied_Inv),
                                         "eDNA_rarefied_Inv", sep="_")
colnames(ps_morph_Inv_2) <- paste(colnames(table_morph_Inv_2), "Morphology_Inv", sep="_")

ps_Inv <- merge(ps_eDNA_Inv_unrarefied, ps_eDNA_Inv_rarefied, by.x=0, by.y=0, all=T)
ps_Inv <- merge(ps_Inv, ps_morph_Inv_2, by.x=1, by.y=0, all=T)
rownames(ps_Inv) <- ps_Inv$Row.names
ps_Inv$Row.names <- NULL
ps_Inv[is.na(ps_Inv)] <- 0
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
p_all <- plot_richness(ps_all_phylo, x='Zone', measures=c("Observed", "Shannon")) + 
  geom_boxplot(outlier.shape = 16, outlier.size = 2, aes(fill=Method)) + 
  scale_fill_manual(values=c("darkolivegreen4", "darkolivegreen3", "lightblue3")) +
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  facet_grid(variable~Organism, scales = "free")
p_all$layers <- p_all$layers[-1]
p_all$data$Zone <- factor(p_all$data$Zone, levels=unique(smpl_eDNA_Fish$Zone))
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

### Calculate relative read abundances assigned
Clytia_hemisphaerica <- sum(table_unrarefied_Inv[rownames(table_unrarefied_Inv) %in% "Clytia hemisphaerica",])
Ch_Total_reads <- sum(table_unrarefied_Inv)
(Clytia_hemisphaerica/Ch_Total_reads)*100

Whiting <- sum(table_unrarefied_Fish[rownames(table_unrarefied_Fish) %in% "Merlangius merlangus",])
Wi_Total_reads <- sum(table_unrarefied_Fish)
(Whiting/Wi_Total_reads)*100

###Additional tests to determine the alpha diversity per Zone based on frequency samples
##12S eDNA metabarcoding on Fish species level
Coast_eDNA_Fish <- subset_samples(ps_unrarefied_Fish, Zone=="Coast")
Transition_eDNA_Fish <- subset_samples(ps_unrarefied_Fish, Zone=="Transition")
Offshore_eDNA_Fish <- subset_samples(ps_unrarefied_Fish, Zone=="Offshore")
seqtab_unrarefied_Fish <- as.data.frame(vegan_otu(ps_unrarefied_Fish))
veganComm_Coast_eDNA_Fish <- vegan_otu(Coast_eDNA_Fish)
veganComm_Transition_eDNA_Fish <- vegan_otu(Transition_eDNA_Fish)
veganComm_Offshore_eDNA_Fish <- vegan_otu(Offshore_eDNA_Fish)

mor_inc_eDNA_Fish <-list ("Coast" =t(veganComm_Coast_eDNA_Fish),
                "Transition"=t(veganComm_Transition_eDNA_Fish),
                "Offshore"=t(veganComm_Offshore_eDNA_Fish))
                
species_incidence_eDNA_Fish <- lapply(mor_inc_eDNA_Fish, as.incfreq)

#Convert to iNEXT format
t_eDNA <- seq(1, 40, by=1)
out.inc_eDNA_Fish <- iNEXT(species_incidence_eDNA_Fish, q=0, 
                  datatype="incidence_freq", size=t_eDNA)

ggiNEXT(out.inc_eDNA_Fish, type=3)

##COI eDNA metabarcoding on invertebrate species level
Coast_eDNA_Inv <- subset_samples(ps_unrarefied_Inv, Zone=="Coast")
Transition_eDNA_Inv <- subset_samples(ps_unrarefied_Inv, Zone=="Transition")
Offshore_eDNA_Inv <- subset_samples(ps_unrarefied_Inv, Zone=="Offshore")
seqtab_unrarefied_Inv <- as.data.frame(vegan_otu(ps_unrarefied_Inv))
veganComm_Coast_eDNA_Inv <- vegan_otu(Coast_eDNA_Inv)
veganComm_Transition_eDNA_Inv <- vegan_otu(Transition_eDNA_Inv)
veganComm_Offshore_eDNA_Inv <- vegan_otu(Offshore_eDNA_Inv)

mor_inc_eDNA_Inv <-list ("Coast" =t(veganComm_Coast_eDNA_Inv),
                          "Transition"=t(veganComm_Transition_eDNA_Inv),
                          "Offshore"=t(veganComm_Offshore_eDNA_Inv))

species_incidence_eDNA_Inv <- lapply(mor_inc_eDNA_Inv, as.incfreq)

#Convert to iNEXT format
out.inc_eDNA_Inv <- iNEXT(species_incidence_eDNA_Inv, q=0, 
                           datatype="incidence_freq", size=t_eDNA)

ggiNEXT(out.inc_eDNA_Inv, type=3)

##Morphological data on Fish species level
Coast_morph_Fish <- subset_samples(ps_morph_Fish_phylo, Zone=="Coast")
Transition_morph_Fish <- subset_samples(ps_morph_Fish_phylo, Zone=="Transition")
Offshore_morph_Fish <- subset_samples(ps_morph_Fish_phylo, Zone=="Offshore")
seqtab_unrarefied_Fish <- as.data.frame(vegan_otu(ps_morph_Fish_phylo))
veganComm_Coast_morph_Fish <- vegan_otu(Coast_morph_Fish)
veganComm_Transition_morph_Fish <- vegan_otu(Transition_morph_Fish)
veganComm_Offshore_morph_Fish <- vegan_otu(Offshore_morph_Fish)

mor_inc_morph_Fish <-list ("Coast" =t(veganComm_Coast_morph_Fish),
                          "Transition"=t(veganComm_Transition_morph_Fish),
                          "Offshore"=t(veganComm_Offshore_morph_Fish))

species_incidence_morph_Fish <- lapply(mor_inc_morph_Fish, as.incfreq)

#Convert to iNEXT format
t_morph <- seq(1, 40, by=1)
out.inc_morph_Fish <- iNEXT(species_incidence_morph_Fish, q=0, 
                           datatype="incidence_freq", size=t_morph)

ggiNEXT(out.inc_morph_Fish, type=3)

##Morphological data  on invertebrate species level
Coast_morph_Inv <- subset_samples(ps_morph_Inv_phylo, Zone=="Coast")
Transition_morph_Inv <- subset_samples(ps_morph_Inv_phylo, Zone=="Transition")
Offshore_morph_Inv <- subset_samples(ps_morph_Inv_phylo, Zone=="Offshore")
seqtab_unrarefied_Inv <- as.data.frame(vegan_otu(ps_morph_Inv_phylo))
veganComm_Coast_morph_Inv <- vegan_otu(Coast_morph_Inv)
veganComm_Transition_morph_Inv <- vegan_otu(Transition_morph_Inv)
veganComm_Offshore_morph_Inv <- vegan_otu(Offshore_morph_Inv)

mor_inc_morph_Inv <-list ("Coast" =t(veganComm_Coast_morph_Inv),
                         "Transition"=t(veganComm_Transition_morph_Inv),
                         "Offshore"=t(veganComm_Offshore_morph_Inv))

species_incidence_morph_Inv <- lapply(mor_inc_morph_Inv, as.incfreq)

#Convert to iNEXT format
out.inc_morph_Inv <- iNEXT(species_incidence_morph_Inv, q=0, 
                          datatype="incidence_freq", size=t_morph)

ggiNEXT(out.inc_morph_Inv, type=3)

