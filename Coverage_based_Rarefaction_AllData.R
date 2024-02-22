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

###Rarefy the 12S eDNA metabarcoding data
## Create a phyloseq object
ps_eDNA_12S <- as.data.frame(table_unrarefied_12S[,colnames(table_unrarefied_12S) %in% env_12S$Niskin.sample])
rownames(ps_eDNA_12S) <- table_unrarefied_12S$ASV
ps_eDNA_12S[is.na(ps_eDNA_12S)] <- 0
smpl_eDNA_12S <- env_12S[env_12S$Niskin.sample  %in%  colnames(ps_eDNA_12S),]
rownames(smpl_eDNA_12S) <- colnames(ps_eDNA_12S)
Taxonomy_eDNA_12S <- as.matrix(table_unrarefied_12S
                                [,(ncol(table_unrarefied_12S)-10):ncol(table_unrarefied_12S)])
rownames(Taxonomy_eDNA_12S) <- table_unrarefied_12S$ASV
ps_unrarefied_12S <- phyloseq(otu_table(ps_eDNA_12S, taxa_are_rows = TRUE),
                               sample_data(smpl_eDNA_12S),
                               tax_table(Taxonomy_eDNA_12S))

##Calculate the Coverage
Coverage_eDNA_12S <- phyloseq_coverage(physeq = ps_unrarefied_12S)

##Rarefy the data based on a coverage just below the minimum coverage using the function phyloseq_coverage_raref,
##(using the minimum coverage will remove the sample from the rarefied dataset)
ps_rarefied_12S <- phyloseq_coverage_raref(physeq = ps_unrarefied_12S, 
                                            iter = 1, coverage = 0.94, drop_lowcoverage = T)
table_rarefied_12S <- as.data.frame(ps_rarefied_12S@otu_table)
table_rarefied_12S <- cbind(table_rarefied_12S, Taxonomy_eDNA_12S)

###Rarefy the COI eDNA metabarcoding data
## Create a phyloseq object
ps_eDNA_COI <- as.data.frame(table_unrarefied_COI[,colnames(table_unrarefied_COI) %in% env_COI$Niskin.sample])
rownames(ps_eDNA_COI) <- table_unrarefied_COI$ASV
ps_eDNA_COI[is.na(ps_eDNA_COI)] <- 0
smpl_eDNA_COI <- env_COI[env_COI$Niskin.sample  %in%  colnames(ps_eDNA_COI),]
rownames(smpl_eDNA_COI) <- colnames(ps_eDNA_COI)
Taxonomy_eDNA_COI <- as.matrix(table_unrarefied_COI
                               [,(ncol(table_unrarefied_COI)-10):ncol(table_unrarefied_COI)])
rownames(Taxonomy_eDNA_COI) <- table_unrarefied_COI$ASV
ps_unrarefied_COI <- phyloseq(otu_table(ps_eDNA_COI, taxa_are_rows = TRUE), 
                              sample_data(smpl_eDNA_COI),tax_table(Taxonomy_eDNA_COI))

##Calculate the Coverage
Coverage_eDNA_COI <- phyloseq_coverage(physeq = ps_unrarefied_COI)

##Rarefy the data based on a minimum coverage using the function phyloseq_coverage_raref,
##(using the minimum coverage that will not generate an error)
ps_rarefied_COI <- phyloseq_coverage_raref(physeq = ps_unrarefied_COI, 
                                           iter = 1, coverage = 0.92, drop_lowcoverage = T)
table_rarefied_COI <- as.data.frame(ps_rarefied_COI@otu_table)
table_rarefied_COI <- cbind(table_rarefied_COI, Taxonomy_eDNA_COI)

###Prepare the 12S eDNA metabarcoding data
##select marine Fish species
fish_classes <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_classes.rds"))
freshwater_fish <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_Freshwater.rds"))
table_rarefied_FishASVs <- as.data.frame(table_rarefied_12S[table_rarefied_12S$Class %in% fish_classes,])
table_rarefied_FishASVs <- as.data.frame(table_rarefied_FishASVs[
  !table_rarefied_FishASVs$Species %in% c(freshwater_fish, "NA"),])
rownames(table_rarefied_FishASVs) <- table_rarefied_FishASVs$ASV
table_rarefied_FishASVs[is.na(table_rarefied_FishASVs)] <- 0

##merge the reads according on fish species level
taxo <- "Species"
table_rarefied_Fish <- aggregate(table_rarefied_FishASVs[,1:(ncol(table_rarefied_FishASVs)-11)], 
                                   by= list(as.factor(table_rarefied_FishASVs[,taxo])),FUN=sum)
rownames(table_rarefied_Fish) <- as.character(table_rarefied_Fish$Group.1)
table_rarefied_Fish$Group.1 <- NULL
table_rarefied_Fish <- table_rarefied_Fish[!rowSums(table_rarefied_Fish) == 0,]
table_rarefied_Fish <- table_rarefied_Fish[,!colSums(table_rarefied_Fish) == 0]

###Prepare the COI eDNA metabarcoding data
##select Invertebrate species
table_rarefied_InvASVs <- as.data.frame(table_rarefied_COI[table_rarefied_COI$Kingdom %in% c("Animalia"),])
table_rarefied_InvASVs <- as.data.frame(table_rarefied_InvASVs[!table_rarefied_InvASVs$Species %in% "NA",])
table_rarefied_InvASVs <- as.data.frame(table_rarefied_InvASVs[!table_rarefied_InvASVs$Phylum %in% "Chordata",])
table_rarefied_InvASVs[is.na(table_rarefied_InvASVs)] <- 0

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

##Rarefy the data based on a minimum coverage using the function phyloseq_coverage_raref,
##(using the minimum coverage that will not generate an error)
ps_morph_rarefied_Fish <- phyloseq_coverage_raref(physeq = ps_morph_Fish_phylo, 
                                                  iter = 1, coverage = 0.94, drop_lowcoverage = T)
table_morph_rarefied_Fish <- as.data.frame(ps_morph_rarefied_Fish@.Data)

plot_Fish_morph_rarefied <- plot_richness(ps_morph_rarefied_Fish, x="Zone",
                                    measures=c("Observed", "Shannon"),
                                    color="Zone") + 
  geom_boxplot(outlier.shape = NA)
plot_Fish_morph_rarefied 


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

##Rarefy the data based on a minimum coverage using the function phyloseq_coverage_raref,
ps_morph_rarefied_Inv <- phyloseq_coverage_raref(physeq = ps_morph_Inv_phylo, 
                                                 iter = 1, coverage = 0.92, drop_lowcoverage = T)
table_morph_rarefied_Inv <- as.data.frame(ps_morph_rarefied_Inv@.Data)

###create a boxplot for the coverage among samples
box_rarefied <- rbind(Coverage_eDNA_12S, Coverage_eDNA_COI, 
                      Coverage_morph_Fish, Coverage_morph_Inv)
box_rarefied$Method <- c(rep("eDNA", nrow(smpl_eDNA_12S)), 
                         rep("eDNA", nrow(smpl_eDNA_COI)),
                         rep("Morphology", nrow(smpl_morph_Fish_raw)),
                         rep("Morphology", nrow(smpl_morph_Inv_raw)))
box_rarefied$Organism <- c(rep("Fish", nrow(smpl_eDNA_12S)), 
                           rep("Invertebrates", nrow(smpl_eDNA_COI)),
                           rep("Fish", nrow(smpl_morph_Fish_raw)),
                           rep("Invertebrates", nrow(smpl_morph_Inv_raw)))
box_rarefied$Zone <- c(smpl_eDNA_12S$Zone, 
                       smpl_eDNA_COI$Zone,
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
smpl_eDNA_Inv_rarefied <- smpl_eDNA_COI[smpl_eDNA_COI$Niskin.sample %in% 
                                          colnames(table_rarefied_Inv),]

smpl_all <- bind_rows(smpl_eDNA_Fish,
                      smpl_eDNA_Inv_rarefied,
                      smpl_morph_Fish_raw,
                      smpl_morph_Inv_raw)

rownames(smpl_all)[1:65] <- paste(rownames(smpl_eDNA_Fish), 
                                  "eDNA_Fish", sep="_")
rownames(smpl_all)[66:127] <- paste(rownames(smpl_eDNA_Inv_rarefied), 
                                    "eDNA_Inv", sep="_")
rownames(smpl_all)[128:149] <- paste(rownames(smpl_morph_Fish_raw),
                                     "Morphology_Fish", sep="_")
rownames(smpl_all)[150:171] <- paste(rownames(smpl_morph_Inv_raw),
                                     "Morphology_Inv", sep="_")

smpl_all$Method <- c(rep("eDNA", nrow(smpl_eDNA_Fish)), 
                      rep("eDNA", nrow(smpl_eDNA_Inv_rarefied)), 
                      rep("Morphology", nrow(smpl_morph_Fish_raw)),
                      rep("Morphology", nrow(smpl_morph_Inv_raw)))
smpl_all$Organism <- c(rep("Fish", nrow(smpl_eDNA_Fish)), 
                        rep("Invertebrates", nrow(smpl_eDNA_Inv_rarefied)), 
                        rep("Fish", nrow(smpl_morph_Fish_raw)),
                        rep("Invertebrates", nrow(smpl_morph_Inv_raw)))

##Prepare the abundance table for the combined fish data
ps_eDNA_Fish_rarefied <- table_rarefied_Fish
ps_morph_Fish_rarefied <- table_morph_rarefied_Fish
ps_eDNA_Inv_rarefied <- table_rarefied_Inv
ps_morph_Inv_rarefied <- table_morph_rarefied_Inv

colnames(ps_eDNA_Fish_rarefied) <- paste(colnames(table_rarefied_Fish),
                                          "eDNA_Fish", sep="_")
colnames(ps_morph_Fish_rarefied) <- paste(colnames(table_morph_rarefied_Fish), 
                                          "Morphology_Fish", sep="_")
colnames(ps_eDNA_Inv_rarefied) <- paste(colnames(table_rarefied_Inv),
                                         "eDNA_Inv", sep="_")
colnames(ps_morph_Inv_rarefied) <- paste(colnames(table_morph_rarefied_Inv), 
                                          "Morphology_Inv", sep="_")

ps_all <- merge(ps_eDNA_Fish_rarefied, ps_morph_Fish_rarefied, by.x=0, by.y=0, all=T)
ps_all <- merge(ps_all, ps_eDNA_Inv_rarefied, by.x=1, by.y=0, all=T)
ps_all <- merge(ps_all, ps_morph_Inv_rarefied, by.x=1, by.y=0, all=T)
rownames(ps_all) <- ps_all$Row.names
ps_all$Row.names <- NULL
ps_all[is.na(ps_all)] <- 0
Taxonomy_all <- as.matrix(rownames(ps_all))
rownames(Taxonomy_all) <- rownames(ps_all)

ps_all_phylo <- phyloseq(otu_table(ps_all , taxa_are_rows = TRUE), sample_data(smpl_all),tax_table(Taxonomy_all))
p_all <- plot_richness(ps_all_phylo, x='Zone', measures=c("Observed", "Shannon")) + 
  geom_boxplot(outlier.shape = 16, outlier.size = 2, aes(fill=Method)) + 
  scale_fill_manual(values=c("darkolivegreen3", "lightblue3")) +
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  facet_grid(variable~Organism, scales = "free")
p_all$layers <- p_all$layers[-1]
p_all$data$Zone <- factor(p_all$data$Zone, levels=unique(smpl_eDNA_Fish$Zone))
p_all




