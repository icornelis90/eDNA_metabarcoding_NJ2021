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
#table_unrarefied_12S <- readxl::read_excel(paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/table_unrarefied_concatenated_CleanedASVs_WithField_FullTaxonomicAssignment.xlsx"))
table_unrarefied_12S <- readxl::read_excel(paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean.xlsx"))
#table_unrarefied_12S <- readxl::read_excel(paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_v2/table_rarefied_CleanedASVs_WithField_FullTaxonomicAssignment.xlsx"))
table_morph_12S <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Standerdized.xlsx"),sheet = "Fish")
table_morph_12S_Standerdized <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Standerdized.xlsx"),sheet = "Fish - Standerdized")
#table_unrarefied_COI <- readxl::read_excel(paste0(proj.path.COI,"/OWFvsCoastal_concatenated/results/Decontam/table_unrarefied_concatenated_CleanedASVs_FullTaxonomicAssignment_WithField.xlsx"))
table_unrarefied_COI <- readxl::read_excel(paste0(proj.path.COI,"/OWFvsCoastal_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean.xlsx"))
#table_unrarefied_COI <- readxl::read_excel(paste0(proj.path.COI,"/OWFvsCoastal_concatenated/results/table_rarefied_CleanedASVs_FullTaxonomicAssignment_WithField.xlsx"))
table_morph_COI <- readxl::read_excel(paste0(proj.path.COI,"/Step5_Statistics/Morphology_Abundancy_Raw.xlsx"),sheet = "Epi")
table_morph_COI_Standerdized <- readxl::read_excel(paste0(proj.path.12S,"/Step5_Statistics/Morphology_Abundancy_Standerdized.xlsx"),sheet = "Epi - Standerdized")
#env <- read.csv(paste0(proj.path.12S,"/Step5_Statistics/environmental_data.csv"),  header=TRUE, sep=";")
#env <- read.csv(paste0(proj.path.12S,"/Step5_Statistics/environmental_data_all3Samples.csv"),  header=TRUE, sep=";")
env_12S <- read.csv(paste0(proj.path.12S,"/Step5_Statistics/environmental_data_OWF3HighestReadNumbers_LocationsTrawl.csv"),  header=TRUE, sep=";")
env_COI <- read.csv(paste0(proj.path.COI,"/Step5_Statistics/environmental_data_OWF3HighestReadNumbers_LocationsTrawl.csv"),  header=TRUE, sep=";")
env_morph <- read.csv(paste0(proj.path.12S,"/Step5_Statistics/environmental_data_morph.csv"),  header=TRUE, sep=";")

#add color and pch to environmental data to create the plot 
env_12S$Environment_color <- ifelse(env_12S$Zone=="Coast","limegreen", 
                                    ifelse(env_12S$Zone=="Transition", "slateblue", 
                                           ifelse(env_12S$Zone=="Offshore","darkorange","red")))
env_COI$Environment_color <- ifelse(env_COI$Zone=="Coast","limegreen", 
                                    ifelse(env_COI$Zone=="Transition", "slateblue", 
                                           ifelse(env_COI$Zone=="Offshore","darkorange","red")))
env_morph$Environment_color <- ifelse(env_morph$Zone=="Coast","limegreen", 
                                      ifelse(env_morph$Zone=="Transition", "slateblue", 
                                             ifelse(env_morph$Zone=="Offshore","darkorange","red")))

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
demersal_fish <- readRDS(file = paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/R_Environment/Demersal_Fish.rds"))
table_unrarefied_FishASVs <- as.data.frame(table_unrarefied_12S[table_unrarefied_12S$Class %in% fish_classes,])
table_unrarefied_FishASVs <- as.data.frame(table_unrarefied_FishASVs[!table_unrarefied_FishASVs$Species %in% c(freshwater_fish, "NA"),])
rownames(table_unrarefied_FishASVs) <- table_unrarefied_FishASVs$ASV
table_unrarefied_FishASVs[is.na(table_unrarefied_FishASVs)] <- 0
taxo <- "Species"
merged_data_unrarefied_12S <- aggregate(table_unrarefied_FishASVs[,1:(ncol(table_unrarefied_FishASVs)-11)], by= list(as.factor(table_unrarefied_FishASVs[,taxo])),FUN=sum)
rownames(merged_data_unrarefied_12S) <- as.character(merged_data_unrarefied_12S$Group.1)
merged_data_unrarefied_12S$Group.1 <- NULL
merged_data_unrarefied_Demersal <- merged_data_unrarefied_12S[rownames(merged_data_unrarefied_12S) %in% demersal_fish,]
merged_data_unrarefied_Demersal  <- merged_data_unrarefied_Demersal[,!colSums(merged_data_unrarefied_Demersal) == 0]

#select Invertebrate species and merge by species
table_unrarefied_AnimaliaASVs <- as.data.frame(table_unrarefied_COI[table_unrarefied_COI$Kingdom %in% c("Animalia"),])
table_unrarefied_AnimaliaASVs <- as.data.frame(table_unrarefied_AnimaliaASVs[!table_unrarefied_AnimaliaASVs$Species %in% "NA",])
table_unrarefied_AnimaliaASVs <- as.data.frame(table_unrarefied_AnimaliaASVs[!table_unrarefied_AnimaliaASVs$Phylum %in% "Chordata",])
taxo <- "Species"
#taxo <- "Full"
merged_data_unrarefied_Animalia <- aggregate(table_unrarefied_AnimaliaASVs[,1:(ncol(table_unrarefied_AnimaliaASVs)-11)], by= list(as.factor(table_unrarefied_AnimaliaASVs[,taxo])),FUN=sum)
rownames(merged_data_unrarefied_Animalia) <-as.character(merged_data_unrarefied_Animalia$Group.1)
merged_data_unrarefied_Animalia$Group.1 <- NULL

#remove negative controls and order colomuns according to zone
env_unrarefied_12S <- env_unrarefied_12S %>% filter(!grepl("neg", env_unrarefied_12S$Niskin.sample))
colorder_12S <- c(env_unrarefied_12S$Niskin.sample)
merged_data_unrarefied_12S <- merged_data_unrarefied_12S[,colorder_12S]
env_unrarefied_COI <- env_unrarefied_COI %>% filter(!grepl("neg", env_unrarefied_COI$Niskin.sample))
colorder_COI <- c(env_unrarefied_COI$Niskin.sample)
merged_data_unrarefied_Animalia <- merged_data_unrarefied_Animalia[,colorder_COI]

##Observed diversity and Shannon diversity index for all datasets seperatly
#eDNA-12S
ps_12S <- merged_data_unrarefied_12S
smpl_12S <- env_unrarefied_12S
rownames(smpl_12S) <- colnames(ps_12S)
Taxonomy_12S <- rownames(merged_data_unrarefied_12S)
Taxonomy_12S <- as.matrix(Taxonomy_12S)
rownames(Taxonomy_12S) <- rownames(ps_12S)
#ps_12S <- phyloseq(otu_table(ps_12S, taxa_are_rows = TRUE), sample_data(smpl_12S),tax_table(Taxonomy_12S))
#plot_12S <- plot_richness(ps_12S, x="Area", measures=c("Observed", "Shannon"), color="Area") + 
  #geom_boxplot(outlier.shape = NA)
#plot_12S

#eDNA-COI
ps_COI <- merged_data_unrarefied_Animalia
smpl_COI <- env_unrarefied_COI
rownames(smpl_COI) <- colnames(ps_COI)
Taxonomy_COI <- rownames(merged_data_unrarefied_Animalia)
Taxonomy_COI <- as.matrix(Taxonomy_COI)
rownames(Taxonomy_COI) <- rownames(ps_COI)
#ps_COI <- phyloseq(otu_table(ps_COI, taxa_are_rows = TRUE), sample_data(smpl_COI),tax_table(Taxonomy_COI))
#plot_COI <- plot_richness(ps_COI, x="Area", measures=c("Observed", "Shannon"), color="Area") + 
  #geom_boxplot(outlier.shape = NA)
#plot_COI

SaveData <- list(ps_Decontam_unrarefied_Fish = ps_12S,
                 ps_Decontam_unrarefied_Inv = ps_COI,
                 smpl_Decontam_unrarefied_Fish = smpl_12S,
                 smpl_Decontam_unrarefied_Inv = smpl_COI)
saveRDS(SaveData, paste0(proj.path.12S,"/MiFish_UE-S_concatenated/results_microDecon/R_Environment/Decontam_unrarefied.rds"))


#Morphology-Fish
ps_morph_12S <- as.data.frame(table_morph_12S_2)
ps_morph_Demersal <- ps_morph_12S[rownames(ps_morph_12S) %in% demersal_fish,]
ps_morph_Demersal  <- ps_morph_Demersal[,!colSums(ps_morph_Demersal) == 0]
smpl_morph_12S <- env_morph
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
rownames(smpl_morph_12S_Std) <- colnames(ps_morph_12S_Std)
Taxonomy_morph_12S_Std <- rownames(table_morph_12S_Std_2)
Taxonomy_morph_12S_Std <- as.matrix(Taxonomy_morph_12S_Std)
rownames(Taxonomy_morph_12S_Std) <- rownames(ps_morph_12S_Std)

ps_morph_COI_Std <- table_morph_COI_Std_2
smpl_morph_COI_Std <- env_morph
rownames(smpl_morph_COI_Std) <- colnames(ps_morph_COI_Std)
Taxonomy_morph_COI_Std <- rownames(table_morph_COI_Std_2)
Taxonomy_morph_COI_Std <- as.matrix(Taxonomy_morph_COI_Std)
rownames(Taxonomy_morph_COI_Std) <- rownames(ps_morph_COI_Std)


##Observed diversity and Shannon diversity index for all Fish data 
#smpl_Fish <- bind_rows(smpl_12S, smpl_morph_12S, smpl_morph_12S_Std)
#smpl_Fish$Method <- c(rep("eDNA", nrow(smpl_12S)), 
#                      rep("Morphology", nrow(smpl_morph_12S)), 
#                      rep("Morphology - Standerdized", nrow(smpl_morph_12S_Std)))
smpl_Fish <- bind_rows(smpl_12S, smpl_12S, smpl_morph_12S, smpl_morph_12S)
smpl_Fish$Method <- c(rep("eDNA_All", nrow(smpl_12S)),
                      rep("eDNA_Demersal", nrow(smpl_12S)), 
                      rep("Morphology_All", nrow(smpl_morph_12S)),
                      rep("Morphology_Demersal", nrow(smpl_morph_12S)))
smpl_Fish$Organism <- "Fish"
rownames(smpl_Fish) <- paste(smpl_Fish$Niskin.sample, smpl_Fish$Method, smpl_Fish$Organism, sep="_")

ps_Fish <- merge(ps_12S, merged_data_unrarefied_Demersal, by.x=0, by.y=0, all=T)
ps_Fish <- merge(ps_Fish, ps_morph_12S, by.x=1, by.y=0, all=T)
ps_Fish <- merge(ps_Fish, ps_morph_Demersal, by.x=1, by.y=0, all=T)
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
#smpl_Inv <- bind_rows(smpl_COI, smpl_morph_COI, smpl_morph_COI_Std)
#smpl_Inv$Method <- c(rep("eDNA", nrow(smpl_COI)), 
#                     rep("Morphology", nrow(smpl_morph_COI)),
#                     rep("Morphology - Standerdized", nrow(smpl_morph_COI_Std)))
smpl_Inv <- bind_rows(smpl_COI, smpl_morph_COI)
smpl_Inv$Method <- c(rep("eDNA_All", nrow(smpl_COI)), 
                     rep("Morphology_All", nrow(smpl_morph_COI)))
smpl_Inv$Organism <- "Invertebrates"
rownames(smpl_Inv) <- paste(smpl_Inv$Niskin.sample, smpl_Inv$Method, smpl_Inv$Organism, sep="_")

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
p_all <- plot_richness(ps_all_phylo, x='Zone', measures=c("Observed", "Shannon")) + 
  geom_boxplot(outlier.shape = 16, outlier.size = 2, aes(fill=Method)) + 
  scale_fill_manual(values=c("darkolivegreen3", "darkolivegreen4", "lightblue3", "lightblue4")) +
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  facet_grid(variable~Organism, scales = "free")
p_all$layers <- p_all$layers[-1]
p_all$data$Zone <- factor(p_all$data$Zone, levels=unique(smpl_12S$Zone))
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

Whiting <- sum(merged_data_unrarefied_12S[rownames(merged_data_unrarefied_12S) %in% "Merlangius merlangus",])
Total_reads <- sum(merged_data_unrarefied_12S)
