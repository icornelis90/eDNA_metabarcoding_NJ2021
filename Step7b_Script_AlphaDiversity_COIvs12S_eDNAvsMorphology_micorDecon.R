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
               , "iNext"
)

for (opties in libraries){
  
  if (opties %in% installed.packages()){
    
    library(opties,character.only = TRUE)
    
  } else {install.packages(opties,repos = "http://cran.us.r-project.org")
    
    library(opties,character.only = TRUE)
  }
}

install.packages("remotes")
remotes::install_github("vmikk/metagMisc")

# make paths
proj.path.12S <- here("/home/genomics/icornelis/02_ZEROimpact/01_12S/NJ2021/MiFish-UE_run2")
proj.path.COI <- here("/home/genomics/icornelis/02_ZEROimpact/02_COI/NJ2021")

#upload data
table_12S <- readxl::read_excel(paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean.xlsx"))
table_morph_12S <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Raw.xlsx"),sheet = "Fish")
table_morph_12S_Standerdized <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Standerdized.xlsx"),sheet = "Fish - Standerdized")

table_COI <- readxl::read_excel(paste0(proj.path.COI,"/OWFvsCoastal_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean.xlsx"))
table_morph_COI <- readxl::read_excel(paste0(proj.path.COI,"/Step5_Statistics/Morphology_Abundancy_Raw.xlsx"),sheet = "Epi")
table_morph_COI_Standerdized <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Standerdized.xlsx"),sheet = "Epi - Standerdized")

env_12S <- readRDS(paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/R_Environment/env_ordered_noNeg.rds"))
env_COI <- readRDS(paste0(proj.path.COI,"/OWFvsCoastal_concatenated/results_microDecon/R_Environment/env_ordered_noNeg.rds"))
env_morph <- read.csv(paste0(proj.path.12S,"/Step5_Statistics/environmental_data_morph.csv"),  header=TRUE, sep=";")

#add color and pch to environmental data to create the plot 
env_morph$Zone_color <- ifelse(env_morph$Zones=="Coastal","limegreen", 
                               ifelse(env_morph$Zones=="zone1", "slateblue", 
                                       ifelse(env_morph$Zones=="zone2","darkorange","red")))

#morphological data
table_morph_12S_2 <- as.data.frame(table_morph_12S[,2:ncol(table_morph_12S)])
table_morph_12S_2 [is.na(table_morph_12S_2 )] <- 0
rownames(table_morph_12S_2) <- as.character(table_morph_12S$...1)
table_morph_COI_2 <- as.data.frame(table_morph_COI[,2:ncol(table_morph_COI)])
table_morph_COI_2 [is.na(table_morph_COI_2 )] <- 0
rownames(table_morph_COI_2) <- as.character(table_morph_COI$...1)

table_morph_12S_Std_2 <- as.data.frame(table_morph_12S_Standerdized[,2:ncol(table_morph_12S_Standerdized)])
table_morph_12S_Std_2 [is.na(table_morph_12S_Std_2 )] <- 0
rownames(table_morph_12S_Std_2) <- as.character(table_morph_12S_Standerdized$...1)
table_morph_12S_Std_2 <- round(table_morph_12S_Std_2, digits = 0)
table_morph_COI_Std_2 <- as.data.frame(table_morph_COI_Standerdized[,2:ncol(table_morph_COI_Standerdized)])
table_morph_COI_Std_2 [is.na(table_morph_COI_Std_2 )] <- 0
rownames(table_morph_COI_Std_2) <- as.character(table_morph_COI_Standerdized$...1)
table_morph_COI_Std_2 <- round(table_morph_COI_Std_2, digits = 0)

#remove samples from environmental data that are absent from the unrarefied data
keep_samples_12S <- c(colnames(table_unrarefied_12S))
env_unrarefied_12S <- env_12S[env_12S$Niskin.sample %in% keep_samples_12S,]
keep_samples_COI <- c(colnames(table_unrarefied_COI))
env_unrarefied_COI <- env_COI[env_COI$Niskin.sample %in% keep_samples_COI,]

#select Fish species and merge by species
fish_classes <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_classes.rds"))
freshwater_fish <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_Freshwater.rds"))
table_unrarefied_FishASVs <- as.data.frame(table_12S[table_12S$Class %in% fish_classes,])
table_unrarefied_FishASVs <- as.data.frame(table_unrarefied_FishASVs[!table_unrarefied_FishASVs$Species %in% c(freshwater_fish, "NA"),])
rownames(table_unrarefied_FishASVs) <- table_unrarefied_FishASVs$ASV
table_unrarefied_FishASVs_noT <- as.data.frame(table_unrarefied_FishASVs[,1:(ncol(table_unrarefied_FishASVs)-11)])
table_unrarefied_FishASVs_noT[is.na(table_unrarefied_FishASVs_noT)] <- 0
table_unrarefied_FishASVs_noT <- table_unrarefied_FishASVs_noT[!rowSums(table_unrarefied_FishASVs_noT) == 0,]
table_unrarefied_FishASVs_noT <- table_unrarefied_FishASVs_noT[,!colSums(table_unrarefied_FishASVs_noT) == 0]

seqtab_unrarefied_FishASVs <- as.data.frame(t(table_unrarefied_FishASVs_noT))
seqtab_unrarefied_FishASVs <- seqtab_unrarefied_FishASVs[!rowSums(seqtab_unrarefied_FishASVs) < 1000,]
ps_12S <- table_unrarefied_FishASVs_noT
ps_12S <- seqtab_unrarefied_FishASVs
smpl_12S <- env_12S[env_12S$Niskin.sample %in% colnames(ps_12S),]
rownames(smpl_12S) <- colnames(ps_12S)
Taxonomy_12S <- as.matrix(table_unrarefied_FishASVs[,(ncol(table_unrarefied_FishASVs)-10):(ncol(table_unrarefied_FishASVs)-4)])

ps_unrarefied_FishASVs <- phyloseq(otu_table(ps_12S, taxa_are_rows = TRUE), sample_data(smpl_12S),tax_table(Taxonomy_12S))

ps_rarefied_FishASVs <- phyloseq_coverage_raref(physeq = ps_unrarefied_FishASVs, 
                                                   iter = 1, coverage = 0.99, drop_lowcoverage = T)
Coverage_FishASVs <- phyloseq_coverage(physeq = ps_unrarefied_FishASVs)

seqtab_rarefied_FishASVs <- as.data.frame(ps_rarefied_FishASVs@otu_table)
plot_12S <- plot_richness(ps_rarefied_FishASVs, x="Zone", measures=c("Observed", "Shannon"), color="Zone") + 
 geom_boxplot(outlier.shape = NA)
plot_12S

rarecurve(seqtab_unrarefied_FishASVs, ylab = "Fish ASVs", 
          main = "Rarecurve of unrarefied samples after taxonomic assignment", 
          col = as.vector(env_12S$Environment_color), label = FALSE, step =100)
rarecurve(seqtab_rarefied_FishASVs, ylab = "Fish ASVs", 
          main = "Rarecurve of unrarefied samples after taxonomic assignment", 
          col = as.vector(smpl_12S$Environment_color), label = FALSE, step =100)
merged_data_rarefied_12S <- phyloseq_coverage_raref(physeq = t(table_unrarefied_Fish),
                                                   iter = 1, coverage = 0.8)
seqtab_rarefied_FishASVs <- rrarefy(seqtab_unrarefied_FishASVs, 10000)


taxo <- "Species"
table_unrarefied_Fish <- aggregate(table_unrarefied_FishASVs[,1:(ncol(table_unrarefied_FishASVs)-11)], by= list(as.factor(table_unrarefied_FishASVs[,taxo])),FUN=sum)
rownames(table_unrarefied_Fish) <- as.character(table_unrarefied_Fish$Group.1)
table_unrarefied_Fish$Group.1 <- NULL
table_unrarefied_Fish[is.na(table_unrarefied_Fish)] <- 0
table_unrarefied_Fish <- table_unrarefied_Fish[!rowSums(table_unrarefied_Fish) == 0,]
table_unrarefied_Fish <- table_unrarefied_Fish[,!colSums(table_unrarefied_Fish) == 0]

#select Invertebrate species and merge by species
table_unrarefied_AnimaliaASVs <- as.data.frame(table_COI[table_COI$Kingdom %in% c("Animalia"),])
table_unrarefied_AnimaliaASVs <- as.data.frame(table_unrarefied_AnimaliaASVs[!table_unrarefied_AnimaliaASVs$Species %in% "NA",])
table_unrarefied_AnimaliaASVs <- as.data.frame(table_unrarefied_AnimaliaASVs[!table_unrarefied_AnimaliaASVs$Phylum %in% "Chordata",])

table_unrarefied_AnimaliaASVs_noT <- as.data.frame(table_unrarefied_AnimaliaASVs[,1:(ncol(table_unrarefied_AnimaliaASVs)-11)])
table_unrarefied_AnimaliaASVs_noT[is.na(table_unrarefied_AnimaliaASVs_noT)] <- 0
table_unrarefied_AnimaliaASVs_noT <- table_unrarefied_AnimaliaASVs_noT[!rowSums(table_unrarefied_AnimaliaASVs_noT) == 0,]
table_unrarefied_AnimaliaASVs_noT <- table_unrarefied_AnimaliaASVs_noT[,!colSums(table_unrarefied_AnimaliaASVs_noT) == 0]


ps_COI <- table_unrarefied_AnimaliaASVs_noT
smpl_COI <- env_COI[env_COI$Niskin.sample %in% colnames(ps_COI),]
rownames(smpl_COI) <- colnames(ps_COI)
Taxonomy_COI <- as.matrix(table_unrarefied_AnimaliaASVs[,(ncol(table_unrarefied_AnimaliaASVs)-10):(ncol(table_unrarefied_AnimaliaASVs)-4)])

ps_unrarefied_AnimaliaASVs <- phyloseq(otu_table(ps_COI, taxa_are_rows = TRUE), sample_data(smpl_COI),tax_table(Taxonomy_COI))

ps_rarefied_AnimaliaASVs <- phyloseq_coverage_raref(physeq = ps_unrarefied_AnimaliaASVs, 
                                                iter = 1, coverage = 0.85, drop_lowcoverage = T)
Coverage_AnimaliaASVs <- phyloseq_coverage(physeq = ps_unrarefied_AnimaliaASVs)

seqtab_rarefied_AnimaliaASVs <- as.data.frame(ps_rarefied_AnimaliaASVs@otu_table)
plot_COI <- plot_richness(ps_rarefied_AnimaliaASVs, x="Zone", measures=c("Observed", "Shannon"), color="Zone") + 
  geom_boxplot(outlier.shape = NA)
plot_COI

rarecurve(t(table_unrarefied_AnimaliaASVs_noT), ylab = "Animalia ASVs", 
          main = "Rarecurve of unrarefied samples after taxonomic assignment", 
          col = as.vector(env_COI$Environment_color), label = FALSE, step =100)
rarecurve(seqtab_rarefied_AnimaliaASVs, ylab = "Animalia ASVs", 
          main = "Rarecurve of unrarefied samples after taxonomic assignment", 
          col = as.vector(smpl_COI$Environment_color), label = FALSE, step =100)
merged_data_rarefied_COI <- phyloseq_coverage_raref(physeq = t(table_unrarefied_Animalia),
                                                    iter = 1, coverage = 0.8)
seqtab_rarefied_AnimaliaASVs <- rrarefy(seqtab_unrarefied_AnimaliaASVs, 10000)




taxo <- "Species"
table_unrarefied_Animalia <- aggregate(table_unrarefied_AnimaliaASVs[,1:(ncol(table_unrarefied_AnimaliaASVs)-11)], by= list(as.factor(table_unrarefied_AnimaliaASVs[,taxo])),FUN=sum)
rownames(table_unrarefied_Animalia) <-as.character(table_unrarefied_Animalia$Group.1)
table_unrarefied_Animalia$Group.1 <- NULL
table_unrarefied_Animalia[is.na(table_unrarefied_Animalia)] <- 0
table_unrarefied_Animalia <- table_unrarefied_Animalia[!rowSums(table_unrarefied_Animalia) == 0,]
table_unrarefied_Animalia <- table_unrarefied_Animalia[,!colSums(table_unrarefied_Animalia) == 0]

#remove negative controls and order colomuns according to zone
env_unrarefied_12S <- env_unrarefied_12S %>% filter(!grepl("neg", env_unrarefied_12S$Niskin.sample))
colorder_12S <- c(env_unrarefied_12S$Niskin.sample)
table_unrarefied_Fish <- table_unrarefied_Fish[,colorder_12S]
env_unrarefied_COI <- env_unrarefied_COI %>% filter(!grepl("neg", env_unrarefied_COI$Niskin.sample))
colorder_COI <- c(env_unrarefied_COI$Niskin.sample)
merged_data_unrarefied_Animalia <- merged_data_unrarefied_Animalia[,colorder_COI]

##Observed diversity and Shannon diversity index for all datasets seperatly
#eDNA-12S
ps_12S <- table_unrarefied_Fish
smpl_12S <- env_12S[env_12S$Niskin.sample %in% colnames(ps_12S),]
rownames(smpl_12S) <- colnames(ps_12S)
Taxonomy_12S <- rownames(table_unrarefied_Fish)
Taxonomy_12S <- as.matrix(Taxonomy_12S)
rownames(Taxonomy_12S) <- rownames(ps_12S)
ps_unrarefied_Fish <- phyloseq(otu_table(ps_12S, taxa_are_rows = TRUE), sample_data(smpl_12S),tax_table(Taxonomy_12S))

ps_rarefied_Fish <- phyloseq_coverage_raref(physeq = ps_unrarefied_Fish, 
                                                iter = 1, coverage = 0.98, drop_lowcoverage = T)
Coverage_Fish <- phyloseq_coverage(physeq = ps_unrarefied_Fish)
seqtab_rarefied_Fish <- as.data.frame(ps_rarefied_Fish@otu_table)
plot_12S <- plot_richness(ps_rarefied_Fish, x="Zone", measures=c("Observed", "Shannon"), color="Zone") + 
  geom_boxplot(outlier.shape = NA)
plot_12S

rarecurve(t(table_unrarefied_Fish), ylab = "Fish ASVs", 
          main = "Rarecurve of unrarefied samples after taxonomic assignment", 
          col = as.vector(env_12S$Environment_color), label = FALSE, step =100)
rarecurve(seqtab_rarefied_Fish, ylab = "Fish ASVs", 
          main = "Rarecurve of unrarefied samples after taxonomic assignment", 
          col = as.vector(env_12S$Environment_color), label = FALSE, step =100)

ps_Animalia <- table_unrarefied_Animalia
smpl_Animalia <- env_COI[env_COI$Niskin.sample %in% colnames(ps_Animalia),]
rownames(smpl_Animalia) <- colnames(ps_Animalia)
Taxonomy_Animalia <- rownames(table_unrarefied_Animalia)
Taxonomy_Animalia <- as.matrix(Taxonomy_Animalia)
rownames(Taxonomy_Animalia) <- rownames(ps_Animalia)
ps_unrarefied_Animalia <- phyloseq(otu_table(ps_Animalia, taxa_are_rows = TRUE), 
                                   sample_data(smpl_Animalia),tax_table(Taxonomy_Animalia))
ps_rarefied_Animalia <- phyloseq_coverage_raref(physeq = ps_unrarefied_Animalia, 
                                            iter = 1, coverage = 0.98, drop_lowcoverage = T)
Coverage_Animalia <- phyloseq_coverage(physeq = ps_unrarefied_Animalia)
Coverage_Animalia <- phyloseq_coverage(physeq = ps_rarefied_Animalia)

seqtab_rarefied_Animalia <- as.data.frame(ps_rarefied_Animalia@otu_table)
plot_12S <- plot_richness(ps_rarefied_Animalia, x="Zone", measures=c("Observed", "Shannon"), color="Zone") + 
  geom_boxplot(outlier.shape = NA)
plot_12S

rarecurve(t(table_unrarefied_Animalia), ylab = "Animalia ASVs", 
          main = "Rarecurve of unrarefied samples after taxonomic assignment", 
          col = as.vector(env_12S$Environment_color), label = FALSE, step =100)
rarecurve(seqtab_rarefied_Animalia, ylab = "Animalia ASVs", 
          main = "Rarecurve of unrarefied samples after taxonomic assignment", 
          col = as.vector(env_12S$Environment_color), label = FALSE, step =100)


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
ps_morph_12S <- as.data.frame(table_morph_12S_2)
smpl_morph_12S <- env_morph
#smpl_morph_12S$Area <- paste(smpl_morph_12S$Environment, smpl_morph_12S$Zones, sep="_")
#smpl_morph_12S$Area <- ifelse(smpl_morph_12S$Area=="Coastal_Coastal","Coastal", smpl_morph_12S$Area)
smpl_morph_12S$Area <- ifelse(smpl_morph_12S$Zones=="Coastal","Coast", 
                        ifelse(smpl_morph_12S$Zones=="zone1", "Transition", 
                               ifelse(smpl_morph_12S$Zones=="zone2","Offshore","red")))
rownames(smpl_morph_12S) <- colnames(ps_morph_12S)
Taxonomy_morph_12S <- rownames(table_morph_12S_2)
Taxonomy_morph_12S <- as.matrix(Taxonomy_morph_12S)
rownames(Taxonomy_morph_12S) <- rownames(ps_morph_12S)
#ps_morph_12S <- phyloseq(otu_table(ps_morph_12S, taxa_are_rows = TRUE), sample_data(smpl_morph_12S),tax_table(Taxonomy_morph_12S))
#plot_richness(ps_morph_12S, x="Area", measures=c("Observed", "Shannon"), color="Area") + 
  #geom_boxplot(outlier.shape = NA)

#Morphology-Invertebrates
ps_morph_COI <- table_morph_COI_2
smpl_morph_COI <- env_morph
#smpl_morph_COI$Area <- paste(smpl_morph_COI$Environment, smpl_morph_COI$Zones, sep="_")
#smpl_morph_COI$Area <- ifelse(smpl_morph_COI$Area=="Coastal_Coastal","Coastal", smpl_morph_COI$Area)
smpl_morph_COI$Area <- ifelse(smpl_morph_COI$Zones=="Coastal","Coast", 
                        ifelse(smpl_morph_COI$Zones=="zone1", "Transition", 
                               ifelse(smpl_morph_COI$Zones=="zone2","Offshore","red")))
rownames(smpl_morph_COI) <- colnames(ps_morph_COI)
Taxonomy_morph_COI <- rownames(table_morph_COI_2)
Taxonomy_morph_COI <- as.matrix(Taxonomy_morph_COI)
rownames(Taxonomy_morph_COI) <- rownames(ps_morph_COI)
#ps_morph_COI <- phyloseq(otu_table(ps_morph_COI, taxa_are_rows = TRUE), sample_data(smpl_morph_COI),tax_table(Taxonomy_morph_COI))
#plot_morph_COI <- plot_richness(ps_morph_COI, x="Area", measures=c("Observed", "Shannon")) + 
  #geom_boxplot(outlier.shape = NA)

#Morphology-Standerdized
ps_morph_12S_Std <- as.data.frame(table_morph_12S_Std_2)
smpl_morph_12S_Std <- env_morph
smpl_morph_12S_Std$Area <- ifelse(smpl_morph_12S_Std$Zones=="Coastal","Coast", 
                              ifelse(smpl_morph_12S_Std$Zones=="zone1", "Transition", 
                                     ifelse(smpl_morph_12S_Std$Zones=="zone2","Offshore","red")))
rownames(smpl_morph_12S_Std) <- colnames(ps_morph_12S_Std)
Taxonomy_morph_12S_Std <- rownames(table_morph_12S_Std_2)
Taxonomy_morph_12S_Std <- as.matrix(Taxonomy_morph_12S_Std)
rownames(Taxonomy_morph_12S_Std) <- rownames(ps_morph_12S_Std)

ps_morph_COI_Std <- table_morph_COI_Std_2
smpl_morph_COI_Std <- env_morph
smpl_morph_COI_Std$Area <- ifelse(smpl_morph_COI_Std$Zones=="Coastal","Coast", 
                              ifelse(smpl_morph_COI_Std$Zones=="zone1", "Transition", 
                                     ifelse(smpl_morph_COI_Std$Zones=="zone2","Offshore","red")))
rownames(smpl_morph_COI_Std) <- colnames(ps_morph_COI_Std)
Taxonomy_morph_COI_Std <- rownames(table_morph_COI_Std_2)
Taxonomy_morph_COI_Std <- as.matrix(Taxonomy_morph_COI_Std)
rownames(Taxonomy_morph_COI_Std) <- rownames(ps_morph_COI_Std)


##Observed diversity and Shannon diversity index for all Fish data 
rownames(smpl_12S) <- paste(colnames(table_unrarefied_Fish), "eDNA_Fish", sep="_")
rownames(smpl_morph_12S) <- paste(colnames(table_morph_12S_2), "Morphology_Fish", sep="_")
rownames(smpl_morph_12S_Std) <- paste(colnames(table_morph_12S_Std_2), "Morphology_Standerdized_Fish", sep="_")
#smpl_Fish <- bind_rows(smpl_12S, smpl_morph_12S, smpl_morph_12S_Std)
#smpl_Fish$Method <- c(rep("eDNA", nrow(smpl_12S)), 
#                      rep("Morphology", nrow(smpl_morph_12S)), 
#                      rep("Morphology - Standerdized", nrow(smpl_morph_12S_Std)))
smpl_Fish <- bind_rows(smpl_12S, smpl_morph_12S)
smpl_Fish$Method <- c(rep("eDNA", nrow(smpl_12S)), 
                      rep("Morphology", nrow(smpl_morph_12S)))
smpl_Fish$Organism <- sub(".*_","",rownames(smpl_Fish))

colnames(ps_12S) <- paste(colnames(table_unrarefied_Fish), "eDNA_Fish", sep="_")
colnames(ps_morph_12S) <- paste(colnames(table_morph_12S_2), "Morphology_Fish", sep="_")
colnames(ps_morph_12S_Std) <- paste(colnames(table_morph_12S_Std_2), "Morphology_Standerdized_Fish", sep="_")

ps_Fish <- merge(ps_12S, ps_morph_12S, by.x=0, by.y=0, all=T)
#ps_Fish <- merge(ps_Fish, ps_morph_12S_Std, by.x=1, by.y=0, all=T)
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
rownames(smpl_morph_COI) <- paste(colnames(table_morph_COI_2), "Morphology_Invertebrates", sep="_")
rownames(smpl_morph_COI_Std) <- paste(colnames(table_morph_COI_Std_2), "Morphology_Standerdized_Invertebrates", sep="_")
#smpl_Inv <- bind_rows(smpl_COI, smpl_morph_COI, smpl_morph_COI_Std)
#smpl_Inv$Method <- c(rep("eDNA", nrow(smpl_COI)), 
#                     rep("Morphology", nrow(smpl_morph_COI)),
#                     rep("Morphology - Standerdized", nrow(smpl_morph_COI_Std)))
smpl_Inv <- bind_rows(smpl_COI, smpl_morph_COI)
smpl_Inv$Organism <- sub(".*_","",rownames(smpl_Inv))
smpl_Inv$Method <- c(rep("eDNA", nrow(smpl_COI)), 
                     rep("Morphology", nrow(smpl_morph_COI)))

colnames(ps_COI) <- paste(colnames(merged_data_unrarefied_Animalia), "eDNA_Invertebrates", sep="_")
colnames(ps_morph_COI) <- paste(colnames(table_morph_COI_2), "morph_Invertebrates", sep="_")
colnames(ps_morph_COI_Std) <- paste(colnames(table_morph_COI_Std_2), "morph_Standerdized_Invertebrates", sep="_")

#ps_Inv <- bind_rows(ps_COI, ps_morph_COI)
ps_Inv <- merge(ps_COI, ps_morph_COI, by.x=0, by.y=0, all=T)
#ps_Inv <- merge(ps_Inv, ps_morph_COI_Std, by.x=1, by.y=0, all=T)
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
