#title: "Alpha Diversity for eDNA metabarcoding data using coverage based rarefaction"
#author: "Isolde Cornelis"
#date: "02/02/2024"

################################################################################

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
table_unrarefied_12S <- readxl::read_excel(paste0(
  proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean.xlsx"))
table_unrarefied_COI <- readxl::read_excel(paste0(
  proj.path.COI,"/OWFvsCoastal_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean.xlsx"))

table_morph_Fish <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Standerdized.xlsx"),
                                       sheet = "Fish")
table_morph_Inv <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Standerdized.xlsx"),
                                      sheet = "Epi")

env_12S <- read.csv(paste0(proj.path.12S,"/Step5_Statistics/environmental_data_OWF3HighestReadNumbers_LocationsTrawl.csv"),
                    header=TRUE, sep=";")
env_COI <- read.csv(paste0(proj.path.COI,"/Step5_Statistics/environmental_data_OWF3HighestReadNumbers_LocationsTrawl.csv"),
                    header=TRUE, sep=";")
env_morph <- read.csv(paste0(proj.path.12S,"/Step5_Statistics/environmental_data_morph.csv"),  header=TRUE, sep=";")

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
demersal_fish <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/R_Environment/Demersal_Fish.rds"))
table_unrarefied_FishASVs <- as.data.frame(table_unrarefied_12S[table_unrarefied_12S$Class %in% fish_classes,])
table_unrarefied_FishASVs <- as.data.frame(table_unrarefied_FishASVs[
  !table_unrarefied_FishASVs$Species %in% c(freshwater_fish, "NA"),])
rownames(table_unrarefied_FishASVs) <- table_unrarefied_FishASVs$ASV
table_unrarefied_FishASVs[is.na(table_unrarefied_FishASVs)] <- 0

###Rarefy the 12S eDNA metabarcoding data
## Create a phyloseq object
ps_eDNA_FishASVs <- as.data.frame(table_unrarefied_FishASVs[,colnames(table_unrarefied_FishASVs) %in% env_12S$Niskin.sample])
rownames(ps_eDNA_12S) <- table_unrarefied_12S$ASV
ps_eDNA_FishASVs[is.na(ps_eDNA_FishASVs)] <- 0
smpl_eDNA_FishASVs <- env_12S[env_12S$Niskin.sample  %in%  colnames(ps_eDNA_FishASVs),]
rownames(smpl_eDNA_FishASVs) <- colnames(ps_eDNA_FishASVs)
Taxonomy_eDNA_FishASVs <- as.matrix(table_unrarefied_FishASVs
                                [,(ncol(table_unrarefied_FishASVs)-10):ncol(table_unrarefied_FishASVs)])
rownames(Taxonomy_eDNA_FishASVs) <- table_unrarefied_FishASVs$ASV
ps_unrarefied_FishASVs <- phyloseq(otu_table(ps_eDNA_FishASVs, taxa_are_rows = TRUE),
                                   sample_data(smpl_eDNA_FishASVs),
                                   tax_table(Taxonomy_eDNA_FishASVs))

##Calculate the Coverage
Coverage_eDNA_FishASVs <- phyloseq_coverage(physeq = ps_unrarefied_FishASVs)
Coverage_eDNA_FishASVs_correct_singletons <- phyloseq_coverage(physeq = ps_unrarefied_FishASVs, correct_singletons = T)

##Rarefy the data based on a coverage just below the minimum coverage using the function phyloseq_coverage_raref,
##(using the minimum coverage will remove the sample from the rarefied dataset)
ps_rarefied_FishASVs <- phyloseq_coverage_raref(physeq = ps_unrarefied_FishASVs, correct_singletons = F, 
                                            iter = 1, coverage = 0.94, drop_lowcoverage = T)
table_rarefied_FishASVs <- as.data.frame(ps_rarefied_FishASVs@otu_table)
table_rarefied_FishASVs <- cbind(table_rarefied_FishASVs, Taxonomy_eDNA_FishASVs)

###Prepare the COI eDNA metabarcoding data
##select Invertebrate species
table_unrarefied_InvASVs <- as.data.frame(table_unrarefied_COI[table_unrarefied_COI$Kingdom %in% c("Animalia"),])
table_unrarefied_InvASVs <- as.data.frame(table_unrarefied_InvASVs[!table_unrarefied_InvASVs$Species %in% "NA",])
table_unrarefied_InvASVs <- as.data.frame(table_unrarefied_InvASVs[!table_unrarefied_InvASVs$Phylum %in% "Chordata",])
table_unrarefied_InvASVs[is.na(table_unrarefied_InvASVs)] <- 0

###Rarefy the COI eDNA metabarcoding data
## Create a phyloseq object
ps_eDNA_InvASVs <- as.data.frame(table_unrarefied_InvASVs[,colnames(table_unrarefied_InvASVs) %in% env_COI$Niskin.sample])
ps_eDNA_InvASVs <- ps_eDNA_InvASVs[,!colSums(ps_eDNA_InvASVs) == 0]
rownames(ps_eDNA_InvASVs) <- table_unrarefied_InvASVs$ASV
ps_eDNA_InvASVs[is.na(ps_eDNA_InvASVs)] <- 0
smpl_eDNA_InvASVs <- env_COI[env_COI$Niskin.sample  %in%  colnames(ps_eDNA_InvASVs),]
rownames(smpl_eDNA_InvASVs) <- colnames(ps_eDNA_InvASVs)
Taxonomy_eDNA_InvASVs <- as.matrix(table_unrarefied_InvASVs
                               [,(ncol(table_unrarefied_InvASVs)-10):ncol(table_unrarefied_InvASVs)])
rownames(Taxonomy_eDNA_InvASVs) <- table_unrarefied_InvASVs$ASV
ps_unrarefied_InvASVs <- phyloseq(otu_table(ps_eDNA_InvASVs, taxa_are_rows = TRUE), 
                              sample_data(smpl_eDNA_InvASVs),tax_table(Taxonomy_eDNA_InvASVs))

##Calculate the Coverage
Coverage_eDNA_InvASVs <- phyloseq_coverage(physeq = ps_unrarefied_InvASVs)
Coverage_eDNA_InvASVs_correct_singletons <- phyloseq_coverage(physeq = ps_unrarefied_InvASVs, correct_singletons = T)

##Rarefy the data based on a minimum coverage using the function phyloseq_coverage_raref,
##(using the minimum coverage that will not generate an error)
ps_rarefied_InvASVs <- phyloseq_coverage_raref(physeq = ps_unrarefied_InvASVs, correct_singletons = F,
                                           iter = 1, coverage = 0.94, drop_lowcoverage = T)
table_rarefied_InvASVs <- as.data.frame(ps_rarefied_InvASVs@otu_table)
table_rarefied_InvASVs <- cbind(table_rarefied_InvASVs, Taxonomy_eDNA_InvASVs)

##merge the reads according on fish species level
taxo <- "Species"
table_rarefied_Fish <- aggregate(table_rarefied_FishASVs[,1:(ncol(table_rarefied_FishASVs)-11)], 
                                   by= list(as.factor(table_rarefied_FishASVs[,taxo])),FUN=sum)
rownames(table_rarefied_Fish) <- as.character(table_rarefied_Fish$Group.1)
table_rarefied_Fish$Group.1 <- NULL
table_rarefied_Fish <- table_rarefied_Fish[!rowSums(table_rarefied_Fish) == 0,]
table_rarefied_Fish <- table_rarefied_Fish[,!colSums(table_rarefied_Fish) == 0]
table_rarefied_Demersal <-  table_rarefied_Fish[rownames(table_rarefied_Fish) %in% demersal_fish,]
table_rarefied_Demersal <- table_rarefied_Demersal[,!colSums(table_rarefied_Demersal) == 0]

##merge the reads according on invertebrate species level
taxo <- "Species"
table_rarefied_Inv <- aggregate(table_rarefied_InvASVs[,1:(ncol(table_rarefied_InvASVs)-11)],
                                  by= list(as.factor(table_rarefied_InvASVs[,taxo])),FUN=sum)
rownames(table_rarefied_Inv) <-as.character(table_rarefied_Inv$Group.1)
table_rarefied_Inv$Group.1 <- NULL
table_rarefied_Inv <- table_rarefied_Inv[!rowSums(table_rarefied_Inv) == 0,]
table_rarefied_Inv <- table_rarefied_Inv[,!colSums(table_rarefied_Inv) == 0]

###Determine alpha diversity
## Create a phyloseq object
smpl_eDNA_Fish <- as.data.frame(env_12S[env_12S$Niskin.sample %in% colnames(table_rarefied_Fish),])
ps_eDNA_Fish <- as.data.frame(table_rarefied_Fish)
rownames(smpl_eDNA_Fish) <- colnames(ps_eDNA_Fish)
Taxonomy_eDNA_Fish <- as.matrix(rownames(table_rarefied_Fish))
rownames(Taxonomy_eDNA_Fish) <- rownames(table_rarefied_Fish)
ps_rarefied_Fish <- phyloseq(otu_table(ps_eDNA_Fish, taxa_are_rows = TRUE),
                             sample_data(smpl_eDNA_Fish),
                             tax_table(Taxonomy_eDNA_Fish))

plot_Fish_rarefied <- plot_richness(ps_rarefied_Fish, x="Zone",
                                    measures=c("Observed", "Shannon"),
                                    color="Zone") + 
  geom_boxplot(outlier.shape = NA)
plot_Fish_rarefied 

## Create a phyloseq object
smpl_eDNA_Inv <- as.data.frame(env_COI[env_COI$Niskin.sample %in% colnames(table_rarefied_Inv),])
ps_eDNA_Inv <- as.data.frame(table_rarefied_Inv)
rownames(smpl_eDNA_Inv) <- colnames(ps_eDNA_Inv)
Taxonomy_eDNA_Inv <- as.matrix(rownames(table_rarefied_Inv))
rownames(Taxonomy_eDNA_Inv) <- rownames(table_rarefied_Inv)
ps_rarefied_Inv <- phyloseq(otu_table(ps_eDNA_Inv, taxa_are_rows = TRUE),
                             sample_data(smpl_eDNA_Inv),
                             tax_table(Taxonomy_eDNA_Inv))

plot_Inv_rarefied <- plot_richness(ps_rarefied_Inv, x="Zone",
                                    measures=c("Observed", "Shannon"),
                                    color="Zone") + 
  geom_boxplot(outlier.shape = NA)
plot_Inv_rarefied 


###Prepare the Morphological data for the fish catch 
## Create a phyloseq object
ps_morph_Fish <- as.data.frame(table_morph_Fish_2)
smpl_morph_Fish_raw <- env_morph
rownames(smpl_morph_Fish_raw) <- colnames(ps_morph_Fish)
Taxonomy_morph_Fish <- as.matrix(rownames(table_morph_Fish_2))
rownames(Taxonomy_morph_Fish) <- rownames(ps_morph_Fish)
ps_morph_Fish_phylo <- phyloseq(otu_table(ps_morph_Fish, taxa_are_rows = TRUE))#,
                                #sample_data(smpl_morph_Fish_raw),
                                #tax_table(Taxonomy_morph_Fish))

##Calculate the Coverage
Coverage_morph_Fish <- phyloseq_coverage(physeq = ps_morph_Fish_phylo)
Coverage_morph_Fish_correct_singletons <- phyloseq_coverage(physeq = ps_morph_Fish_phylo, correct_singletons = T)

##Rarefy the data based on a minimum coverage using the function phyloseq_coverage_raref,
##(using the minimum coverage that will not generate an error)
ps_morph_rarefied_Fish <- phyloseq_coverage_raref(physeq = ps_morph_Fish_phylo, correct_singletons = F, 
                                                  iter = 1, coverage = 0.94, drop_lowcoverage = T)
table_morph_rarefied_Fish <- as.data.frame(ps_morph_rarefied_Fish@.Data)
#table_morph_rarefied_Fish <- as.data.frame(table_morph_Fish_2)

plot_Fish_morph_rarefied <- plot_richness(ps_morph_rarefied_Fish, x="Zone",
                                    measures=c("Observed", "Shannon"),
                                    color="Zone") + 
  geom_boxplot(outlier.shape = NA)
plot_Fish_morph_rarefied 

##Select demersal fish species from the morphological data
table_morph_Demersal <- table_morph_rarefied_Fish[rownames(table_morph_rarefied_Fish) %in% demersal_fish,]
table_morph_Demersal <- table_morph_Demersal[,!colSums(table_morph_Demersal) == 0]

###Prepare the Morphological data for the epibenthos catch 
## Create a phyloseq object
ps_morph_Inv <- table_morph_Inv_2
smpl_morph_Inv_raw <- env_morph
rownames(smpl_morph_Inv_raw) <- colnames(ps_morph_Inv)
Taxonomy_morph_Inv <- as.matrix(rownames(table_morph_Inv_2))
rownames(Taxonomy_morph_Inv) <- rownames(ps_morph_Inv)
ps_morph_Inv_phylo <- phyloseq(otu_table(ps_morph_Inv, taxa_are_rows = TRUE))#,
                               #sample_data(smpl_morph_Inv_raw),
                               #tax_table(Taxonomy_morph_Inv))

##Calculate the Coverage
Coverage_morph_Inv <- phyloseq_coverage(physeq = ps_morph_Inv_phylo)
Coverage_morph_Inv_correct_singletons <- phyloseq_coverage(physeq = ps_morph_Inv_phylo, correct_singletons = T)

##Rarefy the data based on a minimum coverage using the function phyloseq_coverage_raref,
ps_morph_rarefied_Inv <- phyloseq_coverage_raref(physeq = ps_morph_Inv_phylo, correct_singletons = F,
                                                 iter = 1, coverage = 0.94, drop_lowcoverage = T)
table_morph_rarefied_Inv <- as.data.frame(ps_morph_rarefied_Inv@.Data)
#table_morph_rarefied_Inv <- table_morph_Inv_2

###create a boxplot for the coverage among samples
box_rarefied <- rbind(Coverage_eDNA_FishASVs, Coverage_eDNA_InvASVs, 
                      Coverage_morph_Fish, Coverage_morph_Inv)
box_rarefied$Method <- c(rep("eDNA", nrow(smpl_eDNA_FishASVs)), 
                         rep("eDNA", nrow(smpl_eDNA_InvASVs)),
                         rep("Morphology", nrow(smpl_morph_Fish_raw)),
                         rep("Morphology", nrow(smpl_morph_Inv_raw)))
box_rarefied$Organism <- c(rep("Fish", nrow(smpl_eDNA_FishASVs)), 
                           rep("Invertebrates", nrow(smpl_eDNA_InvASVs)),
                           rep("Fish", nrow(smpl_morph_Fish_raw)),
                           rep("Invertebrates", nrow(smpl_morph_Inv_raw)))
box_rarefied$Zone <- c(smpl_eDNA_FishASVs$Zone, 
                       smpl_eDNA_InvASVs$Zone,
                       smpl_morph_Fish_raw$Zone,
                       smpl_morph_Inv_raw$Zone)

box_plot_rarefied <- ggplot(box_rarefied, 
                            aes(x = Zone, y = SampleCoverage)) + 
  geom_boxplot(outlier.shape = 16, outlier.size = 2, aes(fill=Method)) +
  scale_fill_manual(values=c("darkolivegreen3", "lightblue3")) +
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  facet_wrap(~Organism, scales = "free")
box_plot_rarefied

###Create plot with the coverage based rarefied data
##Prepare the smpl matrix for the combined fish data
smpl_eDNA_Fish_Demersal <- smpl_eDNA_Fish[rownames(smpl_eDNA_Fish) %in% 
                                            colnames(table_rarefied_Demersal),]
smpl_morph_Fish_rarefied <- smpl_morph_Fish_raw[rownames(smpl_morph_Fish_raw) %in% 
                                               colnames(table_morph_rarefied_Fish),]
smpl_morph_Demersal <- smpl_morph_Fish_raw[rownames(smpl_morph_Fish_raw) %in% 
                                             colnames(table_morph_Demersal),]
smpl_eDNA_Inv_rarefied <- smpl_eDNA_InvASVs[smpl_eDNA_InvASVs$Niskin.sample %in% 
                                          colnames(table_rarefied_Inv),]
smpl_morph_Inv_rarefied <- smpl_morph_Inv_raw[rownames(smpl_morph_Inv_raw) %in% 
                                              colnames(table_morph_rarefied_Inv),]


smpl_all <- bind_rows(smpl_eDNA_Fish,
                      smpl_eDNA_Fish_Demersal,
                      smpl_eDNA_Inv_rarefied,
                      smpl_morph_Fish_rarefied,
                      smpl_morph_Demersal,
                      smpl_morph_Inv_rarefied)

smpl_all$Method <- c(rep("eDNA_All", nrow(smpl_eDNA_Fish)),
                     rep("eDNA_Demersal", nrow(smpl_eDNA_Fish_Demersal)),
                     rep("eDNA_All", nrow(smpl_eDNA_Inv_rarefied)), 
                     rep("Morphology_All", nrow(smpl_morph_Fish_rarefied)),
                     rep("Morphology_Demersal", nrow(smpl_morph_Demersal)),
                     rep("Morphology_All", nrow(smpl_morph_Inv_rarefied)))
smpl_all$Organism <- c(rep("Fish", nrow(smpl_eDNA_Fish)), 
                       rep("Fish", nrow(smpl_eDNA_Fish_Demersal)), 
                       rep("Invertebrates", nrow(smpl_eDNA_Inv_rarefied)), 
                       rep("Fish", nrow(smpl_morph_Fish_rarefied)),
                       rep("Fish", nrow(smpl_morph_Demersal)),
                       rep("Invertebrates", nrow(smpl_morph_Inv_rarefied)))
rownames(smpl_all) <- paste(smpl_all$Niskin.sample, smpl_all$Method, smpl_all$Organism, sep="_")

##Prepare the abundance table for the combined fish data
ps_eDNA_Fish_rarefied <- table_rarefied_Fish
ps_eDNA_Fish_Demersal <- table_rarefied_Demersal
ps_morph_Fish_rarefied <- table_morph_rarefied_Fish
ps_morph_Fish_Demersal <- table_morph_Demersal
ps_eDNA_Inv_rarefied <- table_rarefied_Inv
ps_morph_Inv_rarefied <- table_morph_rarefied_Inv

ps_all <- merge(ps_eDNA_Fish_rarefied, ps_eDNA_Fish_Demersal, by.x=0, by.y=0, all=T)
ps_all <- merge(ps_all, ps_eDNA_Inv_rarefied, by.x=1, by.y=0, all=T)
ps_all <- merge(ps_all, ps_morph_Fish_rarefied, by.x=1, by.y=0, all=T)
ps_all <- merge(ps_all, ps_morph_Fish_Demersal, by.x=1, by.y=0, all=T)
ps_all <- merge(ps_all, ps_morph_Inv_rarefied, by.x=1, by.y=0, all=T)
rownames(ps_all) <- ps_all$Row.names
ps_all$Row.names <- NULL
ps_all[is.na(ps_all)] <- 0
colnames(ps_all) <- rownames(smpl_all)
Taxonomy_all <- as.matrix(rownames(ps_all))
rownames(Taxonomy_all) <- rownames(ps_all)

ps_all_phylo <- phyloseq(otu_table(ps_all , taxa_are_rows = TRUE), sample_data(smpl_all),tax_table(Taxonomy_all))
p_all <- plot_richness(ps_all_phylo, x='Zone', measures=c("Observed", "Shannon")) + 
  geom_boxplot(outlier.shape = 16, outlier.size = 2, aes(fill=Method)) + 
  scale_fill_manual(values=c("darkolivegreen3", "darkolivegreen4", "lightblue3", "lightblue4")) +
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  facet_grid(variable~Organism, scales = "free")
p_all$layers <- p_all$layers[-1]
p_all$data$Zone <- factor(p_all$data$Zone, levels=unique(smpl_eDNA_Fish$Zone))
p_all

##Save Data
SaveData_eDNA_rarefied <- list(ps_micorDecon_rarefied_Fish = ps_eDNA_Fish_rarefied,
                               ps_micorDecon_rarefied_Inv = ps_eDNA_Inv_rarefied,
                               smpl_micorDecon_rarefied_Fish = smpl_eDNA_12S,
                               smpl_micorDecon_rarefied_Inv = smpl_eDNA_Inv_rarefied)
saveRDS(SaveData_eDNA_rarefied, paste0(proj.path.12S,
                                       "/MiFish_UE-S_concatenated/results_microDecon/R_Environment/microDecon_rarefied.rds"))
SaveData_morph_unrarefied <- list(ps_morph_unrarefied_Fish = ps_morph_Fish,
                                  ps_morph_unrarefied_Inv = ps_morph_Inv,
                                  smpl_morph_unrarefied_Fish = smpl_morph_Fish_raw,
                                  smpl_morph_unrarefied_Inv = smpl_morph_Inv_raw)
saveRDS(SaveData_morph_unrarefied, paste0(proj.path.12S,
                                          "/MiFish_UE-S_concatenated/results_microDecon/R_Environment/Morph_unrarefied.rds"))
SaveData_morph_rarefied <- list(ps_morph_rarefied_Fish = ps_morph_Fish_rarefied,
                                ps_morph_rarefied_Inv = ps_morph_Inv_rarefied,
                                smpl_morph_rarefied_Fish = smpl_morph_Fish_raw,
                                smpl_morph_rarefied_Inv = smpl_morph_Inv_raw)
saveRDS(SaveData_morph_rarefied, paste0(proj.path.12S,
                                        "/MiFish_UE-S_concatenated/results_microDecon/R_Environment/Morph_rarefied.rds"))


##Statistical analysis - Fish
#Alpha diversity
library(car)
library(DescTools)
library(lsr)
DNA_species <- smpl_all
DNA_species$NumOfSp <- colSums(ps_all>0)
DNA_species_Fish <- DNA_species[which(DNA_species$Organism == "Fish"),]
DNA_species_Demersal <- DNA_species_Fish[which(DNA_species_Fish$Method %in% c("eDNA_Demersal", "Morphology_Demersal")),]
DNA_species_AllFish <- DNA_species_Fish[which(DNA_species_Fish$Method %in% c("eDNA_All", "Morphology_All")),]
NumOfSp_AllFish <- as.data.frame(cbind(DNA_species_AllFish$NumOfSp,
                                       DNA_species_AllFish$Zone,
                                       DNA_species_AllFish$Method))
NumOfSp_Demersal <- as.data.frame(cbind(DNA_species_Demersal$NumOfSp,
                                       DNA_species_Demersal$Zone,
                                       DNA_species_Demersal$Method))

#model_DNA <- lm(NumOfSp ~ Zone, data=DNA_species, family=poisson)
#model_DNA_AllFish <- glm(NumOfSp ~ Zone*Method, data=DNA_species_AllFish, family=poisson)
#model_DNA_Demersal <- glm(NumOfSp ~ Zone*Method, data=DNA_species_Demersal, family=poisson)

model_DNA_AllFish <- lm(NumOfSp ~ Zone*Method, data=DNA_species_AllFish)
model_DNA_Demersal <- lm(NumOfSp ~ Zone*Method, data=DNA_species_Demersal) 

anova_DNA_AllFish <- Anova(model_DNA_AllFish, type=2)# SIGN for station , not sign for biological replicates
summary(model_DNA_AllFish)
anova_DNA_AllFish
anova_DNA_Demersal <- Anova(model_DNA_Demersal, type=2)# SIGN for station , not sign for biological replicates
summary(model_DNA_Demersal)
anova_DNA_Demersal
PostHocTest(aov(NumOfSp ~ Zone*Method, data=DNA_species_AllFish), method = "hsd")
PostHocTest(aov(NumOfSp ~ Zone*Method, data=DNA_species_Demersal), method = "hsd")

es <- etaSquared(model_DNA_AllFish, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) # 0.95539

es <- etaSquared(model_DNA_Demersal, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) # 0.872

library(lsmeans)
library(multcomp)
cld(lsmeans(model_DNA_AllFish, ~ Zone*Method), Letters=letters)
cld(lsmeans(model_DNA_Demersal, ~ Zone*Method), Letters=letters)
#elk station in ander groep!
#lstrends(model_DNA, "Station", var="NR") #lstrends: moet continue zijn + volgorde belangrijk => DUS HIER NIET GEBRUIKEN
#DNA_species_2$NR <- as.numeric(as.factor(DNA_species_2$Biologic_replicate))

#check assumptions
#homogeneity
plot(model_DNA_AllFish, 1)
plot(model_DNA_Demersal, 1)
#leveneTest(NumOfSp ~ Station, data=DNA_species) # SIGNIFICANT (but only on edges)
leveneTest(NumOfSp ~ Zone*Method, data=DNA_species_AllFish)
leveneTest(NumOfSp ~ Zone*Method, data=DNA_species_Demersal)

#normality
resid_all <- residuals(model_DNA_AllFish)  # pull the residuals
hist(resid_all)
qqnorm(resid_all) 
resid_demersal <- residuals(model_DNA_Demersal)  # pull the residuals
hist(resid_demersal)
qqnorm(resid_demersal) 
plot(model_DNA_AllFish, 2)
plot(model_DNA_Demersal, 2)
ggqqplot(DNA_species$NumOfSp)
shapiro.test(resid_all)
shapiro.test(resid_demersal)

library(MASS)
boxcox(model_DNA_AllFish)
bc<-boxcox(model_DNA_AllFish)
bc$x[which(bc$y==max(bc$y))] #0.4646465

boxcox(model_DNA_Demersal)
bc<-boxcox(model_DNA_Demersal)
bc$x[which(bc$y==max(bc$y))] #0.989899

#Shannon diversity
DNA_Shannon <- smpl_all
DNA_Shannon$NumOfSp <- (diversity(t(ps_all), "shannon"))
DNA_Shannon_Fish <- DNA_Shannon[which(DNA_Shannon$Organism == "Fish"),]
DNA_Shannon_Demersal <- DNA_Shannon_Fish[which(DNA_Shannon_Fish$Method %in% c("eDNA_Demersal", "Morphology_Demersal")),]
DNA_Shannon_AllFish <- DNA_Shannon_Fish[which(DNA_Shannon_Fish$Method %in% c("eDNA_All", "Morphology_All")),]
NumOfSp_AllFish <- as.data.frame(cbind(DNA_Shannon_AllFish$NumOfSp,
                                       DNA_Shannon_AllFish$Zone,
                                       DNA_Shannon_AllFish$Method))
NumOfSp_Demersal <- as.data.frame(cbind(DNA_Shannon_Demersal$NumOfSp,
                                        DNA_Shannon_Demersal$Zone,
                                        DNA_Shannon_Demersal$Method))

#model_Shannon_AllFish <- glm(NumOfSp ~ Zone*Method, data=DNA_Shannon_AllFish, family=poisson)
#model_Shannon_Demersal <- glm(NumOfSp ~ Zone*Method, data=DNA_Shannon_Demersal, family=poisson) 

model_Shannon_AllFish <- glm(NumOfSp ~ Zone*Method, data=DNA_Shannon_AllFish)
model_Shannon_Demersal <- glm(NumOfSp ~ Zone*Method, data=DNA_Shannon_Demersal) 

anova_Shannon_AllFish <- Anova(model_Shannon_AllFish, type=2)# SIGN for station , not sign for biological replicates
summary(model_Shannon_AllFish)
anova_Shannon_AllFish
anova_Shannon_Demersal <- Anova(model_Shannon_Demersal, type=2)# SIGN for station , not sign for biological replicates
summary(model_Shannon_Demersal)
anova_Shannon_Demersal
PostHocTest(aov(NumOfSp ~ Zone*Method, data=DNA_Shannon_AllFish), method = "hsd")
PostHocTest(aov(NumOfSp ~ Zone*Method, data=DNA_Shannon_Demersal), method = "hsd")

es <- etaSquared(model_Shannon_AllFish, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) # 0.9606682

es <- etaSquared(model_Shannon_Demersal, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) # 1.041221

cld(lsmeans(model_Shannon_AllFish, ~ Zone*Method), Letters=letters)
cld(lsmeans(model_Shannon_Demersal, ~ Zone*Method), Letters=letters)

#check assumptions
#homogeneity
plot(model_Shannon_AllFish, 1)
plot(model_Shannon_Demersal, 1)
leveneTest(NumOfSp ~ Zone*Method, data=DNA_Shannon_AllFish)
leveneTest(NumOfSp ~ Zone*Method, data=DNA_Shannon_Demersal)

#normality
resid_all <- residuals(model_Shannon_AllFish)  # pull the residuals
hist(resid_all)
qqnorm(resid_all) 
resid_demersal <- residuals(model_Shannon_Demersal)  # pull the residuals
hist(resid_demersal)
qqnorm(resid_demersal) 
plot(model_Shannon_AllFish, 2)
plot(model_Shannon_Demersal, 2)
ggqqplot(DNA_species$NumOfSp)
shapiro.test(resid_all)
shapiro.test(resid_demersal)

##Statistical analysis - Invertebrates
#Alpha diversity
DNA_species_Inv <- DNA_species[which(DNA_species$Organism == "Invertebrates"),]

model_DNA_Inv <- lm(NumOfSp ~ Zone*Method, data=DNA_species_Inv) 
model_DNA2_Inv <- glm(NumOfSp ~ Zone*Method, data=DNA_species_Inv, family = poisson) 
anova_DNA_Inv <- Anova(model_DNA_Inv, type=2)# SIGN for station , not sign for biological replicates
summary(model_DNA_Inv)
anova_DNA_Inv
anova_DNA2_Inv <- Anova(model_DNA2_Inv, type=2)# SIGN for station , not sign for biological replicates
summary(model_DNA2_Inv)
anova_DNA2_Inv
PostHocTest(aov(NumOfSp ~ Zone*Method, data=DNA_species_Inv), method = "hsd")

es <- etaSquared(model_DNA_Inv, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) # 1.000002

es <- etaSquared(model_DNA2_Inv, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) # 0.8668192

cld(lsmeans(model_DNA_Inv, ~ Zone*Method), Letters=letters)
cld(lsmeans(model_DNA2_Inv, ~ Zone*Method), Letters=letters)

#elk station in ander groep!
#lstrends(model_DNA, "Station", var="NR") #lstrends: moet continue zijn + volgorde belangrijk => DUS HIER NIET GEBRUIKEN
#DNA_species_2$NR <- as.numeric(as.factor(DNA_species_2$Biologic_replicate))

#check assumptions
#homogeneity
plot(model_DNA_Inv, 1)
plot(model_DNA2_Inv, 1)
leveneTest(NumOfSp ~ Zone*Method, data=DNA_species_Inv)

#normality
resid <- residuals(model_DNA_Inv)  # pull the residuals
resid2 <- residuals(model_DNA2_Inv)  # pull the residuals
hist(resid)
hist(resid2)
qqnorm(resid)
qqnorm(resid2) 
plot(model_DNA_Inv, 2)
plot(model_DNA2_Inv, 2)
ggqqplot(DNA_species_Inv$NumOfSp)
shapiro.test(resid)
shapiro.test(resid2)

boxcox(model_DNA_Inv)
bc<-boxcox(model_DNA_Inv)
bc$x[which(bc$y==max(bc$y))] #0.6666667

boxcox(model_DNA2_Inv)
bc<-boxcox(model_DNA2_Inv)
bc$x[which(bc$y==max(bc$y))] #0.6666667

#Shannon diversity
DNA_Shannon_Inv <- DNA_Shannon[which(DNA_Shannon$Organism == "Invertebrates"),]

model_Shannon_Inv <- lm(NumOfSp ~ Zone*Method, data=DNA_Shannon_Inv)
model_Shannon2_Inv <- glm(NumOfSp ~ Zone*Method, data=DNA_Shannon_Inv, family = poisson)
anova_Shannon_Inv <- Anova(model_Shannon_Inv, type=2)# SIGN for station , not sign for biological replicates
summary(model_Shannon_Inv)
anova_Shannon_Inv
anova_Shannon2_Inv <- Anova(model_Shannon2_Inv, type=2)# SIGN for station , not sign for biological replicates
summary(model_Shannon2_Inv)
anova_Shannon2_Inv
PostHocTest(aov(NumOfSp ~ Zone*Method, data=DNA_Shannon_Inv), method = "hsd")

es <- etaSquared(model_Shannon_Inv, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) #1.00072 

es <- etaSquared(model_Shannon2_Inv, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) #0.7441658

cld(lsmeans(model_Shannon_Inv, ~ Zone*Method), Letters=letters)
cld(lsmeans(model_Shannon2_Inv, ~ Zone*Method), Letters=letters)

#check assumptions
#homogeneity
plot(model_Shannon_Inv, 1)
plot(model_Shannon2_Inv, 1)
leveneTest(NumOfSp ~ Zone*Method, data=DNA_Shannon_Inv)

#normality
resid_Shannon_Inv <- residuals(model_Shannon_Inv)  # pull the residuals
resid_Shannon2_Inv <- residuals(model_Shannon2_Inv)  # pull the residuals
hist(resid_Shannon_Inv)
hist(resid_Shannon2_Inv)
qqnorm(resid_Shannon_Inv) 
qqnorm(resid_Shannon2_Inv) 
plot(model_Shannon_Inv, 2)
plot(model_Shannon2_Inv, 2)
ggqqplot(DNA_Shannon$Shannon)
shapiro.test(resid_Shannon_Inv)
shapiro.test(resid_Shannon2_Inv)

boxcox(model_Shannon_Inv)
bc<-boxcox(model_Shannon_Inv)
bc$x[which(bc$y==max(bc$y))]

boxcox(model_Shannon2_Inv)
bc<-boxcox(model_Shannon2_Inv)
bc$x[which(bc$y==max(bc$y))]

## Calculate relative read abundances assigned
Select_Transition <- c(smpl_eDNA_Inv$Niskin.sample[which(smpl_eDNA_Inv$Zone == "Transition")])
rarefied_eDNA_Inv_Transition <- table_rarefied_Inv[, colnames(table_rarefied_Inv) %in% Select_Transition]
Clytia_hemisphaerica <- sum(rarefied_eDNA_Inv_Transition[rownames(rarefied_eDNA_Inv_Transition) %in% "Clytia hemisphaerica",])
Total_reads <- sum(rarefied_eDNA_Inv_Transition)
Clytia_hemisphaerica/Total_reads

Select_Transition <- c(smpl_eDNA_Fish$Niskin.sample[which(smpl_eDNA_Fish$Zone == "Transition")])
rarefied_eDNA_Fish_Transition <- table_rarefied_Fish[, colnames(table_rarefied_Fish) %in% Select_Transition]
Whiting <- sum(rarefied_eDNA_Fish_Transition[rownames(rarefied_eDNA_Fish_Transition) %in% "Merlangius merlangus",])
Total_reads <- sum(rarefied_eDNA_Fish_Transition)
Whiting/Total_reads
Dab <- sum(rarefied_eDNA_Fish_Transition[rownames(rarefied_eDNA_Fish_Transition) %in% "Limanda limanda",])
Dab/Total_reads

Select_Transition <- c(smpl_morph_Fish_raw$Niskin.sample[which(smpl_morph_Fish_raw$Zone == "Transition")])
rarefied_morph_Transition <- table_morph_rarefied_Fish[, colnames(table_morph_rarefied_Fish) %in% Select_Transition]
Whiting <- sum(rarefied_morph_Transition[rownames(rarefied_morph_Transition) %in% "Merlangius merlangus",])
Total_reads <- sum(rarefied_morph_Transition)
Whiting/Total_reads
Dab <- sum(rarefied_morph_Transition[rownames(rarefied_morph_Transition) %in% "Limanda limanda",])
Dab/Total_reads




