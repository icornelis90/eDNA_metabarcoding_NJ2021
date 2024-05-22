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
               , "plotrix"
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
Coverage_eDNA_12S_correct_singletons <- phyloseq_coverage(physeq = ps_unrarefied_12S, correct_singletons = TRUE)

##Rarefy the data based on a coverage just below the minimum coverage using the function phyloseq_coverage_raref,
##(using the minimum coverage will remove the sample from the rarefied dataset)
ps_rarefied_12S <- phyloseq_coverage_raref(physeq = ps_unrarefied_12S, 
                                            iter = 1, coverage = 0.94, drop_lowcoverage = T)
table_rarefied_12S <- as.data.frame(ps_rarefied_12S@otu_table)
table_rarefied_12S <- cbind(table_rarefied_12S, Taxonomy_eDNA_12S)
table_unrarefied_12S_2 <- cbind(ps_eDNA_12S, Taxonomy_eDNA_12S)

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
Coverage_eDNA_COI_correct_singletons <- phyloseq_coverage(physeq = ps_unrarefied_COI, correct_singletons = TRUE)

##Rarefy the data based on a minimum coverage using the function phyloseq_coverage_raref,
##(using the minimum coverage that will not generate an error)
ps_rarefied_COI <- phyloseq_coverage_raref(physeq = ps_unrarefied_COI, 
                                           iter = 1, coverage = 0.92, drop_lowcoverage = T)
table_rarefied_COI <- as.data.frame(ps_rarefied_COI@otu_table)
table_rarefied_COI <- cbind(table_rarefied_COI, Taxonomy_eDNA_COI)
table_unrarefied_COI_2 <- cbind(ps_eDNA_COI, Taxonomy_eDNA_COI)

###Prepare the 12S eDNA metabarcoding data
##select marine Fish species
fish_classes <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_classes.rds"))
freshwater_fish <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_Freshwater.rds"))
demersal_fish <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/R_Environment/Demersal_Fish.rds"))

table_unrarefied_FishASVs <- as.data.frame(table_unrarefied_12S_2[table_unrarefied_12S_2$Class %in% fish_classes,])
table_unrarefied_FishASVs <- as.data.frame(table_unrarefied_FishASVs[
  !table_unrarefied_FishASVs$Species %in% c(freshwater_fish, "NA"),])
rownames(table_unrarefied_FishASVs) <- table_unrarefied_FishASVs$ASV
table_unrarefied_FishASVs[is.na(table_unrarefied_FishASVs)] <- 0

table_rarefied_FishASVs <- as.data.frame(table_rarefied_12S[table_rarefied_12S$Class %in% fish_classes,])
table_rarefied_FishASVs <- as.data.frame(table_rarefied_FishASVs[
  !table_rarefied_FishASVs$Species %in% c(freshwater_fish, "NA"),])
rownames(table_rarefied_FishASVs) <- table_rarefied_FishASVs$ASV
table_rarefied_FishASVs[is.na(table_rarefied_FishASVs)] <- 0

##merge the reads according on fish species level
taxo <- "Species"

table_unrarefied_Fish <- aggregate(table_unrarefied_FishASVs[,1:(ncol(table_unrarefied_FishASVs)-11)], 
                                   by= list(as.factor(table_unrarefied_FishASVs[,taxo])),FUN=sum)
rownames(table_unrarefied_Fish) <- as.character(table_unrarefied_Fish$Group.1)
table_unrarefied_Fish$Group.1 <- NULL
table_unrarefied_Fish <- table_unrarefied_Fish[!rowSums(table_unrarefied_Fish) == 0,]
table_unrarefied_Fish <- table_unrarefied_Fish[,!colSums(table_unrarefied_Fish) == 0]
table_unrarefied_Demersal <-  table_unrarefied_Fish[rownames(table_unrarefied_Fish) %in% demersal_fish,]
table_unrarefied_Demersal <- table_unrarefied_Demersal[,!colSums(table_unrarefied_Demersal) == 0]

table_rarefied_Fish <- aggregate(table_rarefied_FishASVs[,1:(ncol(table_rarefied_FishASVs)-11)], 
                                 by= list(as.factor(table_rarefied_FishASVs[,taxo])),FUN=sum)
rownames(table_rarefied_Fish) <- as.character(table_rarefied_Fish$Group.1)
table_rarefied_Fish$Group.1 <- NULL
table_rarefied_Fish <- table_rarefied_Fish[!rowSums(table_rarefied_Fish) == 0,]
table_rarefied_Fish <- table_rarefied_Fish[,!colSums(table_rarefied_Fish) == 0]
table_rarefied_Demersal <-  table_rarefied_Fish[rownames(table_rarefied_Fish) %in% demersal_fish,]
table_rarefied_Demersal <- table_rarefied_Demersal[,!colSums(table_rarefied_Demersal) == 0]

###Prepare the COI eDNA metabarcoding data
##select Invertebrate species
table_unrarefied_InvASVs <- as.data.frame(table_unrarefied_COI_2[table_unrarefied_COI_2$Kingdom %in% c("Animalia"),])
table_unrarefied_InvASVs <- as.data.frame(table_unrarefied_InvASVs[!table_unrarefied_InvASVs$Species %in% "NA",])
table_unrarefied_InvASVs <- as.data.frame(table_unrarefied_InvASVs[!table_unrarefied_InvASVs$Phylum %in% "Chordata",])
table_unrarefied_InvASVs[is.na(table_unrarefied_InvASVs)] <- 0

table_rarefied_InvASVs <- as.data.frame(table_rarefied_COI[table_rarefied_COI$Kingdom %in% c("Animalia"),])
table_rarefied_InvASVs <- as.data.frame(table_rarefied_InvASVs[!table_rarefied_InvASVs$Species %in% "NA",])
table_rarefied_InvASVs <- as.data.frame(table_rarefied_InvASVs[!table_rarefied_InvASVs$Phylum %in% "Chordata",])
table_rarefied_InvASVs[is.na(table_rarefied_InvASVs)] <- 0

##merge the reads according on invertebrate species level
table_unrarefied_Inv <- aggregate(table_unrarefied_InvASVs[,1:(ncol(table_unrarefied_InvASVs)-11)],
                                  by= list(as.factor(table_unrarefied_InvASVs[,taxo])),FUN=sum)
rownames(table_unrarefied_Inv) <-as.character(table_unrarefied_Inv$Group.1)
table_unrarefied_Inv$Group.1 <- NULL
table_unrarefied_Inv <- table_unrarefied_Inv[!rowSums(table_unrarefied_Inv) == 0,]
table_unrarefied_Inv <- table_unrarefied_Inv[,!colSums(table_unrarefied_Inv) == 0]

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
Coverage_morph_Fish_correct_singletons <- phyloseq_coverage(physeq = ps_morph_Fish_phylo, correct_singletons = TRUE)

##Rarefy the data based on a minimum coverage using the function phyloseq_coverage_raref,
##(using the minimum coverage that will not generate an error)
ps_morph_rarefied_Fish <- phyloseq_coverage_raref(physeq = ps_morph_Fish_phylo, 
                                                  iter = 1, coverage = 0.94, drop_lowcoverage = T)
table_morph_rarefied_Fish <- as.data.frame(ps_morph_rarefied_Fish@.Data)

##Select demersal fish species from the morphological data
table_morph_unrarefied_Demersal <- ps_morph_Fish[rownames(ps_morph_Fish) %in% demersal_fish,]
table_morph_unrarefied_Demersal <- table_morph_unrarefied_Demersal[,!colSums(table_morph_unrarefied_Demersal) == 0]

table_morph_rarefied_Demersal <- table_morph_rarefied_Fish[rownames(table_morph_rarefied_Fish) %in% demersal_fish,]
table_morph_rarefied_Demersal <- table_morph_rarefied_Demersal[,!colSums(table_morph_rarefied_Demersal) == 0]

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
  scale_fill_manual(values=c("lightblue3", "palevioletred2")) +
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  facet_wrap(~Organism, scales = "free")
box_plot_rarefied

###Create plot with the coverage based rarefied data
##Prepare the smpl matrix for the combined fish data
smpl_eDNA_Inv_rarefied <- smpl_eDNA_COI[smpl_eDNA_COI$Niskin.sample %in% 
                                          colnames(table_rarefied_Inv),]
smpl_morph_Demersal <- smpl_morph_Fish_raw[rownames(smpl_morph_Fish_raw) %in% 
                                       colnames(table_morph_rarefied_Demersal),]

smpl_all_rarefied <- bind_rows(smpl_eDNA_Fish,
                               smpl_eDNA_Fish,
                               smpl_eDNA_Inv_rarefied,
                               smpl_morph_Fish_raw,
                               smpl_morph_Demersal,
                               smpl_morph_Inv_raw)

rownames(smpl_all_rarefied)[1:66] <- paste(rownames(smpl_eDNA_Fish), 
                                           "eDNA_Fish_All", sep="_")
rownames(smpl_all_rarefied)[67:132] <- paste(rownames(smpl_eDNA_Fish), 
                                             "eDNA_Fish_Demersal", sep="_")
rownames(smpl_all_rarefied)[133:194] <- paste(rownames(smpl_eDNA_Inv_rarefied), 
                                              "eDNA_Inv", sep="_")
rownames(smpl_all_rarefied)[195:216] <- paste(rownames(smpl_morph_Fish_raw),
                                              "Morphology_Fish_All", sep="_")
rownames(smpl_all_rarefied)[217:237] <- paste(rownames(smpl_morph_Demersal),
                                              "Morphology_Fish_Demersal", sep="_")
rownames(smpl_all_rarefied)[238:259] <- paste(rownames(smpl_morph_Inv_raw),
                                              "Morphology_Inv", sep="_")

smpl_all_rarefied$Method <- c(rep("eDNA_All", nrow(smpl_eDNA_Fish)),
                              rep("eDNA_Demersal", nrow(smpl_eDNA_Fish)),
                              rep("eDNA_All", nrow(smpl_eDNA_Inv_rarefied)), 
                              rep("Morphology_All", nrow(smpl_morph_Fish_raw)),
                              rep("Morphology_Demersal", nrow(smpl_morph_Demersal)),
                              rep("Morphology_All", nrow(smpl_morph_Inv_raw)))
smpl_all_rarefied$Organism <- c(rep("Fish", nrow(smpl_eDNA_Fish)), 
                                rep("Fish", nrow(smpl_eDNA_Fish)), 
                                rep("Invertebrates", nrow(smpl_eDNA_Inv_rarefied)), 
                                rep("Fish", nrow(smpl_morph_Fish_raw)),
                                rep("Fish", nrow(smpl_morph_Demersal)),
                                rep("Invertebrates", nrow(smpl_morph_Inv_raw)))

##Prepare the abundance table for the combined fish data
ps_eDNA_Fish_rarefied <- table_rarefied_Fish
ps_eDNA_Fish_rarefied_nonpelagic <- table_rarefied_Demersal
ps_morph_Fish_rarefied <- table_morph_rarefied_Fish
ps_morph_Fish_rarefied_Demersal <- table_morph_rarefied_Demersal
ps_eDNA_Inv_rarefied <- table_rarefied_Inv
ps_morph_Inv_rarefied <- table_morph_rarefied_Inv

colnames(ps_eDNA_Fish_rarefied) <- paste(colnames(table_rarefied_Fish),
                                                  "eDNA_Fish_All", sep="_")
colnames(ps_eDNA_Fish_rarefied_nonpelagic) <- paste(colnames(ps_eDNA_Fish_rarefied_nonpelagic),
                                                    "eDNA_Fish_Demersal", sep="_")
colnames(ps_morph_Fish_rarefied) <- paste(colnames(table_morph_rarefied_Fish), 
                                          "Morphology_Fish_All", sep="_")
colnames(ps_morph_Fish_rarefied_Demersal) <- paste(colnames(ps_morph_Fish_rarefied_Demersal), 
                                                   "Morphology_Fish_Demersal", sep="_")
colnames(ps_eDNA_Inv_rarefied) <- paste(colnames(table_rarefied_Inv),
                                        "eDNA_Inv", sep="_")
colnames(ps_morph_Inv_rarefied) <- paste(colnames(table_morph_rarefied_Inv), 
                                         "Morphology_Inv", sep="_")

ps_all_rarefied <- merge(ps_eDNA_Fish_rarefied, ps_eDNA_Fish_rarefied_nonpelagic, by.x=0, by.y=0, all=T)
ps_all_rarefied <- merge(ps_all_rarefied, ps_morph_Fish_rarefied, by.x=1, by.y=0, all=T)
ps_all_rarefied <- merge(ps_all_rarefied, ps_morph_Fish_rarefied_Demersal, by.x=1, by.y=0, all=T)
ps_all_rarefied <- merge(ps_all_rarefied, ps_eDNA_Inv_rarefied, by.x=1, by.y=0, all=T)
ps_all_rarefied <- merge(ps_all_rarefied, ps_morph_Inv_rarefied, by.x=1, by.y=0, all=T)
rownames(ps_all_rarefied) <- ps_all_rarefied$Row.names
ps_all_rarefied$Row.names <- NULL
ps_all_rarefied[is.na(ps_all_rarefied)] <- 0
Taxonomy_all_rarefied <- as.matrix(rownames(ps_all_rarefied))
rownames(Taxonomy_all_rarefied) <- rownames(ps_all_rarefied)

ps_all_rarefied_phylo <- phyloseq(otu_table(ps_all_rarefied, taxa_are_rows = TRUE),
                                  sample_data(smpl_all_rarefied),tax_table(Taxonomy_all_rarefied))
#p_all_rarefied_box <- plot_richness(ps_all_rarefied_phylo, x='Zone', measures=c("Observed", "Shannon")) + 
p_all_rarefied_box <- plot_richness(ps_all_rarefied_phylo, x='Zone', measures=c("Observed")) + 
  geom_boxplot(outlier.shape = 16, outlier.size = 2, aes(fill=Method)) + 
  scale_fill_manual(values=c("lightblue3", "lightblue4", "palevioletred2", "palevioletred4")) +
  #scale_fill_manual(values=c("darkolivegreen3", "darkolivegreen4", "lightblue3", "lightblue4")) +
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  labs(y = "Richness")+
  facet_grid(variable~Organism, scales = "free")
p_all_rarefied_box$layers <- p_all_rarefied_box$layers[-1]
p_all_rarefied_box$data$Zone <- factor(p_all_rarefied_box$data$Zone, levels=unique(smpl_eDNA_Fish$Zone))
p_all_rarefied_box

p_all_rarefied_smooth <- plot_richness(ps_all_rarefied_phylo, x='Zone', measures=c("Observed")) + 
  #geom_point(aes(color=Method), size = 3) + 
  geom_smooth(method=NULL, aes(color=Method, fill=Method,as.numeric(Zone), value)) +
  scale_fill_manual(values=c("lightblue3", "lightblue4", "palevioletred2", "palevioletred4")) +
  scale_color_manual(values=c("lightblue3", "lightblue4", "palevioletred2", "palevioletred4")) +
  #scale_fill_manual(values=c("darkolivegreen3", "darkolivegreen4", "lightblue3", "lightblue4")) +
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  facet_grid(variable~Organism, scales = "free")
p_all_rarefied_smooth$layers <- p_all_rarefied_smooth$layers[-1]
p_all_rarefied_smooth$data$Zone <- factor(p_all_rarefied_smooth$data$Zone, levels=unique(smpl_eDNA_12S$Zone))
p_all_rarefied_smooth

p_all_rarefied_line <- estimate_richness(ps_all_rarefied_phylo, measures="Observed", split = T)
p_all_rarefied_line <- cbind(p_all_rarefied_line, ps_all_rarefied_phylo@sam_data)
p_all_rarefied_line_plot <- p_all_rarefied_line %>%
  group_by(Organism, Method, Zone) %>%
  summarise(
    n=n(),
    mean=mean(Observed),
    sd=sd(Observed)) %>%
  mutate(se=sd/sqrt(n))  %>%
  mutate(ci=se * qt((1-0.05)/2 + .5, n-1))

p_all_rarefied_line2 <- ggplot(p_all_rarefied_line_plot, aes(x=Zone, y=mean, group = Method)) + 
  geom_point(aes(color=Method), size = 3) + 
  #geom_errorbar(aes(color=Method, ymin = mean + sd, ymax = mean - sd), width = 0.2, size = 0.5) +
  #geom_errorbar(aes(color=Method, ymin = mean + se, ymax = mean - se), width = 0.2, size = 0.5) +
  geom_errorbar(aes(color=Method, ymin = mean + ci, ymax = mean - ci), width = 0.2, size = 0.5) +
  geom_line(linetype = "longdash", aes(color=Method), size = 1)+
  scale_color_manual(values=c("lightblue3", "lightblue4", "palevioletred2", "palevioletred4")) +
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme_bw()+
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  labs(y = "Richness")+
  facet_grid(~Organism, scales = "free")
p_all_rarefied_line2$data$Zone <- factor(p_all_rarefied_line2$data$Zone, levels=unique(smpl_eDNA_Fish$Zone))
p_all_rarefied_line2

###Create plot with the un-rarefied data
##Prepare the smpl matrix for the combined fish data
smpl_eDNA_Inv_unrarefied <- smpl_eDNA_COI[smpl_eDNA_COI$Niskin.sample %in% 
                                          colnames(table_unrarefied_Inv),]
smpl_morph_Demersal <- smpl_morph_Fish_raw[rownames(smpl_morph_Fish_raw) %in% 
                                             colnames(table_morph_unrarefied_Demersal),]

smpl_all_unrarefied <- bind_rows(smpl_eDNA_12S,
                                 smpl_eDNA_12S,
                                 smpl_eDNA_Inv_unrarefied,
                                 smpl_morph_Fish_raw,
                                 smpl_morph_Demersal,
                                 smpl_morph_Inv_raw)

rownames(smpl_all_unrarefied)[1:66] <- paste(rownames(smpl_eDNA_12S), 
                                             "eDNA_Fish_All", sep="_")
rownames(smpl_all_unrarefied)[67:132] <- paste(rownames(smpl_eDNA_12S), 
                                               "eDNA_Fish_Demersal", sep="_")
rownames(smpl_all_unrarefied)[133:197] <- paste(rownames(smpl_eDNA_Inv_unrarefied), 
                                                "eDNA_Inv", sep="_")
rownames(smpl_all_unrarefied)[198:219] <- paste(rownames(smpl_morph_Fish_raw),
                                                "Morphology_Fish_All", sep="_")
rownames(smpl_all_unrarefied)[220:241] <- paste(rownames(smpl_morph_Demersal),
                                                "Morphology_Fish_Demersal", sep="_")
rownames(smpl_all_unrarefied)[242:263] <- paste(rownames(smpl_morph_Inv_raw),
                                                "Morphology_Inv", sep="_")

smpl_all_unrarefied$Method <- c(rep("eDNA_All", nrow(smpl_eDNA_12S)),
                                rep("eDNA_Demersal", nrow(smpl_eDNA_12S)),
                                rep("eDNA_All", nrow(smpl_eDNA_Inv_unrarefied)), 
                                rep("Morphology_All", nrow(smpl_morph_Fish_raw)),
                                rep("Morphology_Demersal", nrow(smpl_morph_Demersal)),
                                rep("Morphology_All", nrow(smpl_morph_Inv_raw)))
smpl_all_unrarefied$Organism <- c(rep("Fish", nrow(smpl_eDNA_12S)), 
                                  rep("Fish", nrow(smpl_eDNA_12S)), 
                                  rep("Invertebrates", nrow(smpl_eDNA_Inv_unrarefied)), 
                                  rep("Fish", nrow(smpl_morph_Fish_raw)),
                                  rep("Fish", nrow(smpl_morph_Demersal)),
                                  rep("Invertebrates", nrow(smpl_morph_Inv_raw)))

##Prepare the abundance table for the combined fish data
ps_eDNA_Fish_unrarefied <- table_unrarefied_Fish
ps_eDNA_Fish_unrarefied_nonpelagic <- table_unrarefied_Demersal
ps_morph_Fish_unrarefied <- table_morph_Fish_2
ps_morph_Fish_unrarefied_Demersal <- table_morph_unrarefied_Demersal
ps_eDNA_Inv_unrarefied <- table_unrarefied_Inv
ps_morph_Inv_unrarefied <- table_morph_Inv_2

colnames(ps_eDNA_Fish_unrarefied) <- paste(colnames(table_unrarefied_Fish),
                                           "eDNA_Fish_All", sep="_")
colnames(ps_eDNA_Fish_unrarefied_nonpelagic) <- paste(colnames(ps_eDNA_Fish_unrarefied_nonpelagic),
                                                      "eDNA_Fish_Demersal", sep="_")
colnames(ps_morph_Fish_unrarefied) <- paste(colnames(ps_morph_Fish_unrarefied), 
                                            "Morphology_Fish_All", sep="_")
colnames(ps_morph_Fish_unrarefied_Demersal) <- paste(colnames(ps_morph_Fish_unrarefied_Demersal), 
                                                     "Morphology_Fish_Demersal", sep="_")
colnames(ps_eDNA_Inv_unrarefied) <- paste(colnames(table_unrarefied_Inv),
                                          "eDNA_Inv", sep="_")
colnames(ps_morph_Inv_unrarefied) <- paste(colnames(ps_morph_Inv_unrarefied), 
                                           "Morphology_Inv", sep="_")

ps_all_unrarefied <- merge(ps_eDNA_Fish_unrarefied, ps_eDNA_Fish_unrarefied_nonpelagic, by.x=0, by.y=0, all=T)
ps_all_unrarefied <- merge(ps_all_unrarefied, ps_eDNA_Inv_unrarefied, by.x=1, by.y=0, all=T)
ps_all_unrarefied <- merge(ps_all_unrarefied, ps_morph_Fish_unrarefied, by.x=1, by.y=0, all=T)
ps_all_unrarefied <- merge(ps_all_unrarefied, ps_morph_Fish_unrarefied_Demersal, by.x=1, by.y=0, all=T)
ps_all_unrarefied <- merge(ps_all_unrarefied, ps_morph_Inv_unrarefied, by.x=1, by.y=0, all=T)
rownames(ps_all_unrarefied) <- ps_all_unrarefied$Row.names
ps_all_unrarefied$Row.names <- NULL
ps_all_unrarefied[is.na(ps_all_unrarefied)] <- 0
Taxonomy_all_unrarefied <- as.matrix(rownames(ps_all_unrarefied))
rownames(Taxonomy_all_unrarefied) <- rownames(ps_all_unrarefied)

ps_all_unrarefied_phylo <- phyloseq(otu_table(ps_all_unrarefied , taxa_are_rows = TRUE), 
                                    sample_data(smpl_all_unrarefied),tax_table(Taxonomy_all_unrarefied))
#p_all_unrarefied_box <- plot_richness(ps_all_unrarefied_phylo, x='Zone', measures=c("Observed", "Shannon")) + 
p_all_unrarefied_box <- plot_richness(ps_all_unrarefied_phylo, x='Zone', measures=c("Observed")) + 
  geom_boxplot(outlier.shape = 16, outlier.size = 2, aes(fill=Method)) + 
  scale_fill_manual(values=c("lightblue3", "lightblue4", "palevioletred2", "palevioletred4")) +
  #scale_fill_manual(values=c("darkolivegreen3", "darkolivegreen4", "lightblue3", "lightblue4")) +
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  labs(y = "Richness")+
  facet_grid(variable~Organism, scales = "free")
p_all_unrarefied_box$layers <- p_all_unrarefied_box$layers[-1]
p_all_unrarefied_box$data$Zone <- factor(p_all_unrarefied_box$data$Zone, levels=unique(smpl_eDNA_12S$Zone))
p_all_unrarefied_box

p_all_unrarefied_smooth <- plot_richness(ps_all_unrarefied_phylo, x='Zone', measures=c("Observed")) + 
  #geom_point(aes(color=Method), size = 3) + 
  geom_smooth(method=NULL, aes(color=Method, fill=Method, as.numeric(Zone), value)) +
  scale_fill_manual(values=c("lightblue3", "lightblue4", "palevioletred2", "palevioletred4")) +
  scale_color_manual(values=c("lightblue3", "lightblue4", "palevioletred2", "palevioletred4")) +
  #scale_fill_manual(values=c("darkolivegreen3", "darkolivegreen4", "lightblue3", "lightblue4")) +
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  facet_grid(variable~Organism, scales = "free")
p_all_unrarefied_smooth$layers <- p_all_unrarefied_smooth$layers[-1]
p_all_unrarefied_smooth$data$Zone <- factor(p_all_unrarefied_smooth$data$Zone, levels=unique(smpl_eDNA_12S$Zone))
p_all_unrarefied_smooth

p_all_unrarefied_line <- estimate_richness(ps_all_unrarefied_phylo, measures="Observed", split = T)
p_all_unrarefied_line <- cbind(p_all_unrarefied_line, ps_all_unrarefied_phylo@sam_data)
p_all_unrarefied_line_plot <- p_all_unrarefied_line %>%
  group_by(Organism, Method, Zone) %>%
  summarise(
    n=n(),
    mean=mean(Observed),
    sd=sd(Observed)) %>%
  mutate(se=sd/sqrt(n))  %>%
  mutate(ci=se * qt((1-0.05)/2 + .5, n-1))

p_all_unrarefied_line2 <- ggplot(p_all_unrarefied_line_plot, aes(x=Zone, y=mean, group = Method)) + 
  geom_point(aes(color=Method), size = 3) + 
  geom_errorbar(aes(color=Method, ymin = mean + sd, ymax = mean - sd), width = 0.2, size = 0.5) +
  #geom_errorbar(aes(color=Method, ymin = mean + se, ymax = mean - se), width = 0.2, size = 0.5) +
  #geom_errorbar(aes(color=Method, ymin = mean + ci, ymax = mean - ci), width = 0.2, size = 0.5) +
  geom_line(linetype = "longdash", aes(color=Method), size = 1)+
  scale_color_manual(values=c("lightblue3", "lightblue4", "palevioletred2", "palevioletred4")) +
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme_bw()+
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  labs(y = "Richness")+
  facet_grid(~Organism, scales = "free")
p_all_unrarefied_line2$data$Zone <- factor(p_all_unrarefied_line2$data$Zone, levels=unique(smpl_eDNA_12S$Zone))
p_all_unrarefied_line2

##Save Data
SaveData_eDNA_rarefied <- list(ps_micorDecon_rarefied_Fish = ps_eDNA_Fish_rarefied,
                               ps_micorDecon_rarefied_Inv = ps_eDNA_Inv_rarefied,
                               smpl_micorDecon_rarefied_Fish = smpl_eDNA_12S,
                               smpl_micorDecon_rarefied_Inv = smpl_eDNA_Inv_rarefied)
saveRDS(SaveData_eDNA_rarefied, paste0(proj.path.12S,
                                       "/MiFish_UE-S_concatenated/results_microDecon/R_Environment/microDecon_rarefied.rds"))
SaveData_eDNA_unrarefied <- list(ps_micorDecon_unrarefied_Fish = ps_eDNA_Fish_unrarefied,
                                 ps_micorDecon_unrarefied_Inv = ps_eDNA_Inv_unrarefied,
                                 smpl_micorDecon_unrarefied_Fish = smpl_eDNA_12S,
                                 smpl_micorDecon_unrarefied_Inv = smpl_eDNA_COI)
saveRDS(SaveData_eDNA_unrarefied, paste0(proj.path.12S,
                                       "/MiFish_UE-S_concatenated/results_microDecon/R_Environment/microDecon_unrarefied.rds"))
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
DNA_species_rarefied <- smpl_all_rarefied
DNA_species_rarefied$NumOfSp <- colSums(ps_all_rarefied>0)
DNA_species_rarefied_Fish <- DNA_species_rarefied[which(DNA_species_rarefied$Organism == "Fish"),]
DNA_species_rarefied_Demersal <- DNA_species_rarefied_Fish[which(DNA_species_rarefied_Fish$Method 
                                                                 %in% c("eDNA_Demersal", "Morphology_Demersal")),]
DNA_species_rarefied_AllFish <- DNA_species_rarefied_Fish[which(DNA_species_rarefied_Fish$Method 
                                                                %in% c("eDNA_All", "Morphology_All")),]
NumOfSp_rarefied_AllFish <- as.data.frame(cbind(DNA_species_rarefied_AllFish$NumOfSp,
                                       DNA_species_rarefied_AllFish$Zone,
                                       DNA_species_rarefied_AllFish$Method))
NumOfSp_rarefied_Demersal <- as.data.frame(cbind(DNA_species_rarefied_Demersal$NumOfSp,
                                       DNA_species_rarefied_Demersal$Zone,
                                       DNA_species_rarefied_Demersal$Method))

DNA_species_unrarefied <- smpl_all_unrarefied
DNA_species_unrarefied$NumOfSp <- colSums(ps_all_unrarefied>0)
DNA_species_unrarefied_Fish <- DNA_species_unrarefied[which(DNA_species_unrarefied$Organism == "Fish"),]
DNA_species_unrarefied_Demersal <- DNA_species_unrarefied_Fish[which(DNA_species_unrarefied_Fish$Method 
                                                                 %in% c("eDNA_Demersal", "Morphology_Demersal")),]
DNA_species_unrarefied_AllFish <- DNA_species_unrarefied_Fish[which(DNA_species_unrarefied_Fish$Method 
                                                                %in% c("eDNA_All", "Morphology_All")),]
NumOfSp_unrarefied_AllFish <- as.data.frame(cbind(DNA_species_unrarefied_AllFish$NumOfSp,
                                                  DNA_species_unrarefied_AllFish$Zone,
                                                  DNA_species_unrarefied_AllFish$Method))
NumOfSp_unrarefied_Demersal <- as.data.frame(cbind(DNA_species_unrarefied_Demersal$NumOfSp,
                                                   DNA_species_unrarefied_Demersal$Zone,
                                                   DNA_species_unrarefied_Demersal$Method))

#model_DNA <- lm(NumOfSp ~ Zone, data=DNA_species, family=poisson)
model_DNA_rarefied_AllFish <- glm(NumOfSp ~ Zone*Method,
                                  data=DNA_species_rarefied_AllFish,
                                  family=poisson)
model_DNA_rarefied_AllFish_Z <- glm(NumOfSp ~ Zone,
                                    data=DNA_species_rarefied_AllFish,
                                    family=poisson)
model_DNA_rarefied_AllFish_M <- glm(NumOfSp ~ Method,
                                    data=DNA_species_rarefied_AllFish,
                                     family=poisson)
model_DNA_rarefied_Demersal <- glm(NumOfSp ~ Zone*Method,
                                   data=DNA_species_rarefied_Demersal,
                                   family=poisson)
model_DNA_rarefied_Demersal_Z <- glm(NumOfSp ~ Zone,
                                     data=DNA_species_rarefied_Demersal,
                                     family=poisson)
model_DNA_rarefied_Demersal_M <- glm(NumOfSp ~ Method,
                                     data=DNA_species_rarefied_Demersal,
                                     family=poisson)

summary(model_DNA_rarefied_AllFish)
summary(model_DNA_rarefied_Demersal)
summary(model_DNA_rarefied_AllFish_Z)
summary(model_DNA_rarefied_AllFish_M)
summary(model_DNA_rarefied_Demersal_Z)
summary(model_DNA_rarefied_Demersal_M)

anova_DNA_rarefied_AllFish <- Anova(model_DNA_rarefied_AllFish,
                                    type=2,
                                    test = "F")
anova_DNA_rarefied_AllFish
anova_DNA_rarefied_Demersal <- Anova(model_DNA_rarefied_Demersal,
                                     type=2,
                                     test = "F")
anova_DNA_rarefied_Demersal
anova_DNA_rarefied_AllFish_Z <- Anova(model_DNA_rarefied_AllFish_Z,
                                      type=2,
                                      test = "F")
anova_DNA_rarefied_AllFish_Z
anova_DNA_rarefied_AllFish_M <- Anova(model_DNA_rarefied_AllFish_M,
                                      type=2,
                                      test = "F")
anova_DNA_rarefied_AllFish_M
anova_DNA_rarefied_Demersal_Z <- Anova(model_DNA_rarefied_Demersal_Z,
                                       type=2,
                                       test = "F")
anova_DNA_rarefied_Demersal_Z
anova_DNA_rarefied_Demersal_M <- Anova(model_DNA_rarefied_Demersal_M,
                                       type=2,
                                       test = "F")
anova_DNA_rarefied_Demersal_M

PostHocTest(aov(NumOfSp ~ Zone*Method,
                data=DNA_species_rarefied_AllFish),
            method = "hsd")
PostHocTest(aov(NumOfSp ~ Zone,
                data=DNA_species_rarefied_AllFish),
            method = "hsd")
PostHocTest(aov(NumOfSp ~ Method,
                data=DNA_species_rarefied_AllFish),
            method = "hsd")
PostHocTest(aov(NumOfSp ~ Zone*Method,
                data=DNA_species_rarefied_Demersal),
            method = "hsd")
PostHocTest(aov(NumOfSp ~ Zone,
                data=DNA_species_rarefied_AllFish),
            method = "hsd")
PostHocTest(aov(NumOfSp ~ Method,
                data=DNA_species_rarefied_AllFish),
            method = "hsd")

model_DNA_unrarefied_AllFish <- glm(NumOfSp ~ Zone*Method,
                                    data=DNA_species_unrarefied_AllFish,
                                    family=poisson)
model_DNA_unrarefied_AllFish_Z <- glm(NumOfSp ~ Zone,
                                      data=DNA_species_unrarefied_AllFish,
                                      family=poisson)
model_DNA_unrarefied_AllFish_M <- glm(NumOfSp ~ Method,
                                      data=DNA_species_unrarefied_AllFish,
                                      family=poisson)
model_DNA_unrarefied_Demersal <- glm(NumOfSp ~ Zone*Method,
                                     data=DNA_species_unrarefied_Demersal,
                                     family=poisson)
model_DNA_unrarefied_Demersal_Z <- glm(NumOfSp ~ Zone,
                                       data=DNA_species_unrarefied_Demersal,
                                       family=poisson)
model_DNA_unrarefied_Demersal_M <- glm(NumOfSp ~ Method,
                                       data=DNA_species_unrarefied_Demersal,
                                       family=poisson)

summary(model_DNA_unrarefied_AllFish)
summary(model_DNA_unrarefied_Demersal)
summary(model_DNA_unrarefied_AllFish_Z)
summary(model_DNA_unrarefied_AllFish_M)
summary(model_DNA_unrarefied_Demersal_Z)
summary(model_DNA_unrarefied_Demersal_M)

anova_DNA_unrarefied_AllFish <- Anova(model_DNA_unrarefied_AllFish,
                                      type=2,
                                      test = "F")
anova_DNA_unrarefied_AllFish
anova_DNA_unrarefied_Demersal <- Anova(model_DNA_unrarefied_Demersal,
                                       type=2,
                                       test = "F")
anova_DNA_unrarefied_Demersal
anova_DNA_unrarefied_AllFish_Z <- Anova(model_DNA_unrarefied_AllFish_Z,
                                        type=2,
                                        test = "F")
anova_DNA_unrarefied_AllFish_Z
anova_DNA_unrarefied_Demersal_Z <- Anova(model_DNA_unrarefied_Demersal_Z,
                                         type=2,
                                         test = "F")
anova_DNA_unrarefied_Demersal_Z
anova_DNA_unrarefied_AllFish_M <- Anova(model_DNA_unrarefied_AllFish_M,
                                        type=2,
                                        test = "F")
anova_DNA_unrarefied_AllFish_M
anova_DNA_unrarefied_Demersal_M <- Anova(model_DNA_unrarefied_Demersal_M,
                                         type=2,
                                         test = "F")
anova_DNA_unrarefied_Demersal_M

PostHocTest(aov(NumOfSp ~ Zone*Method,
                data=DNA_species_unrarefied_AllFish),
            method = "hsd")
PostHocTest(aov(NumOfSp ~ Zone,
                data=DNA_species_unrarefied_AllFish),
            method = "hsd")
PostHocTest(aov(NumOfSp ~ Method,
                data=DNA_species_unrarefied_AllFish),
            method = "hsd")
PostHocTest(aov(NumOfSp ~ Zone*Method,
                data=DNA_species_unrarefied_Demersal),
            method = "hsd")
PostHocTest(aov(NumOfSp ~ Zone,
                data=DNA_species_unrarefied_Demersal),
            method = "hsd")
PostHocTest(aov(NumOfSp ~ Method,
                data=DNA_species_unrarefied_Demersal),
            method = "hsd")

es <- etaSquared(model_DNA_rarefied_AllFish, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) # 0.9554
es <- etaSquared(model_DNA_rarefied_Demersal, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) # 0.8728

es <- etaSquared(model_DNA_unrarefied_AllFish, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) # 0.9523
es <- etaSquared(model_DNA_unrarefied_Demersal, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) # 0.8972

library(lsmeans)
library(multcomp)
cld(lsmeans(model_DNA_unrarefied_AllFish, ~ Zone*Method), Letters=letters)
cld(lsmeans(model_DNA_unrarefied_Demersal, ~ Zone*Method), Letters=letters)

#check assumptions
#homogeneity
plot(model_DNA_unrarefied_AllFish, 1)
plot(model_DNA_unrarefied_Demersal, 1)
#leveneTest(NumOfSp ~ Station, data=DNA_species) # SIGNIFICANT (but only on edges)
leveneTest(NumOfSp ~ Zone*Method, data=DNA_species_rarefied_AllFish)
leveneTest(NumOfSp ~ Zone*Method, data=DNA_species_rarefied_Demersal)
leveneTest(NumOfSp ~ Zone*Method, data=DNA_species_unrarefied_AllFish)
leveneTest(NumOfSp ~ Zone*Method, data=DNA_species_unrarefied_Demersal)

#normality
resid_all_rarefied <- residuals(model_DNA_rarefied_AllFish)  # pull the residuals
hist(resid_all_rarefied)
qqnorm(resid_all_rarefied) 
plot(model_DNA_rarefied_AllFish, 2)
shapiro.test(resid_all_rarefied)

resid_rarefied_demersal <- residuals(model_DNA_rarefied_Demersal)  # pull the residuals
hist(resid_rarefied_demersal)
qqnorm(resid_rarefied_demersal) 
plot(model_DNA_rarefied_Demersal, 2)
shapiro.test(resid_rarefied_demersal)

ggqqplot(DNA_species_rarefied$NumOfSp)

resid_all_unrarefied <- residuals(model_DNA_unrarefied_AllFish)  # pull the residuals
hist(resid_all_unrarefied)
qqnorm(resid_all_unrarefied) 
plot(model_DNA_unrarefied_AllFish, 2)
shapiro.test(resid_all_unrarefied)

resid_unrarefied_demersal <- residuals(model_DNA_unrarefied_Demersal)  # pull the residuals
hist(resid_unrarefied_demersal)
qqnorm(resid_unrarefied_demersal) 
plot(model_DNA_unrarefied_Demersal, 2)
shapiro.test(resid_unrarefied_demersal)

ggqqplot(DNA_species_unrarefied$NumOfSp)

##Statistical analysis - Invertebrates
#Alpha diversity
DNA_species_rarefied_Inv <- DNA_species_rarefied[which(DNA_species_rarefied$Organism == "Invertebrates"),]
DNA_species_unrarefied_Inv <- DNA_species_unrarefied[which(DNA_species_unrarefied$Organism == "Invertebrates"),]

model_DNA_rarefied_Inv <- glm(NumOfSp ~ Zone*Method,
                              data=DNA_species_rarefied_Inv,
                              family = poisson)
model_DNA_rarefied_Inv_Z <- glm(NumOfSp ~ Zone,
                                data=DNA_species_rarefied_Inv,
                                family = poisson)

summary(model_DNA_rarefied_Inv)
summary(model_DNA_rarefied_Inv_Z)

anova_DNA_rarefied_Inv <- Anova(model_DNA_rarefied_Inv,
                                type=2,
                                test = "F")
anova_DNA_rarefied_Inv
anova_DNA_rarefied_Inv_Z <- Anova(model_DNA_rarefied_Inv_Z,
                                  type=2,
                                  test = "F")
anova_DNA_rarefied_Inv_Z
PostHocTest(aov(NumOfSp ~ Zone*Method,
                data=DNA_species_rarefied_Inv),
            method = "hsd")
PostHocTest(aov(NumOfSp ~ Zone*Method,
                data=DNA_species_rarefied_Inv_Z),
            method = "hsd")

model_DNA_unrarefied_Inv <- glm(NumOfSp ~ Zone*Method,
                                data=DNA_species_unrarefied_Inv,
                                family = poisson)
model_DNA_unrarefied_Inv_Z <- glm(NumOfSp ~ Zone,
                                  data=DNA_species_unrarefied_Inv,
                                  family = poisson)

summary(model_DNA_unrarefied_Inv)
summary(model_DNA_unrarefied_Inv_Z)

anova_DNA_unrarefied_Inv <- Anova(model_DNA_unrarefied_Inv,
                                  type=2,
                                  test = "F")
anova_DNA_unrarefied_Inv
anova_DNA_unrarefied_Inv_Z <- Anova(model_DNA_unrarefied_Inv_Z,
                                    type=2,
                                    test = "F")
anova_DNA_unrarefied_Inv_Z
PostHocTest(aov(NumOfSp ~ Zone*Method,
                data=DNA_species_unrarefied_Inv),
            method = "hsd")
PostHocTest(aov(NumOfSp ~ Zone*Method,
                data=DNA_species_unrarefied_Inv_Z),
            method = "hsd")

es <- etaSquared(model_DNA_rarefied_Inv, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) # 0.8668

es <- etaSquared(model_DNA_unrarefied_Inv, type=2, anova=TRUE)
es
sum(es[,"eta.sq"]) # 0.9467

cld(lsmeans(model_DNA_rarefied_Inv, ~ Zone*Method), Letters=letters)
cld(lsmeans(model_DNA_unrarefied_Inv, ~ Zone*Method), Letters=letters)

#check assumptions
#homogeneity
plot(model_DNA_rarefied_Inv, 1)
plot(model_DNA_unrarefied_Inv, 1)
leveneTest(NumOfSp ~ Zone*Method, data=DNA_species_rarefied_Inv)
leveneTest(NumOfSp ~ Zone*Method, data=DNA_species_unrarefied_Inv)

#normality
resid_rarefied <- residuals(model_DNA_rarefied_Inv)  # pull the residuals
hist(resid_rarefied)
qqnorm(resid_rarefied)
plot(model_DNA_rarefied_Inv, 2)
shapiro.test(resid_rarefied)

resid_unrarefied <- residuals(model_DNA_unrarefied_Inv)  # pull the residuals
hist(resid_unrarefied)
qqnorm(resid_unrarefied)
plot(model_DNA_unrarefied_Inv, 2)
shapiro.test(resid_unrarefied)

ggqqplot(DNA_species_rarefied_Inv$NumOfSp)
ggqqplot(DNA_species_unrarefied_Inv$NumOfSp)

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



