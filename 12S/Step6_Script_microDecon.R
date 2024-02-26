#set OS type for paths
if(.Platform$OS.type == "unix"){home<-"/home/"} else{
  home<-"//192.168.236.131/"}

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(dplyr)
library(microbiome)
library(vegan)
library(openxlsx)
library(tibble)
library(tidyverse)
library(ade4)
library(here)
library(devtools)
library(microDecon)
#devtools::install_github("donaldtmcknight/microDecon")

# make paths
proj.path <- here("/home/genomics/icornelis/02_ZEROimpact/01_12S/NJ2021/MiFish-UE_run2")

# upload the dataset
table_raw <- readxl::read_excel(paste0(proj.path, "/MiFish_UE-S_concatenated/results_v2/table_unrarefied_raw_Full_TaxAss_WoRMS.xlsx"))

# Remove the taxonomic assignment
colnames(table_raw) <- str_replace(colnames(table_raw), "_1234", "")
colnames(table_raw) <- str_replace(colnames(table_raw), "_control", "_coastal")
colnames(table_raw) <- str_replace(colnames(table_raw), "_extraction_filter", "_DNA")
table_raw_noT <- as.data.frame(table_raw[,1:(ncol(table_raw)-11)])

#create environmental data
env <- as.data.frame(unique(str_sub(colnames(table_raw_noT), end=-4)))
colnames(env)[1] <- 'Niskin.sample'
Coast <- readRDS(file = here("/home/genomics/icornelis/03_RawScripts/Locations_Coast.rds"))
Coast <- c(Coast, "1", "13", "20", "23", "24", "3", "37", "38", "6", "92", "93")
Transition <- readRDS(file = here("/home/genomics/icornelis/03_RawScripts/Locations_Transition.rds"))
Transition <- c(Transition, "ftWT2bis")
Offshore <- readRDS(file = here("/home/genomics/icornelis/03_RawScripts/Locations_Offshore.rds"))
env$Location <- ifelse(grepl("neg", env$Niskin.sample), 
                       str_sub(env$Niskin.sample, end=-3),
                       gsub("\\_.*", "", env$Niskin.sample))
env$Biol_replicate <- ifelse(grepl("DYF", env$Niskin.sample), 
                             sub("^([^_]+_){1}([^_]+).*", "\\2", env$Niskin.sample),
                             gsub("^.*\\_", "", env$Niskin.sample))
env$Zone <- ifelse(env$Location %in% Coast, "Coast",
                   ifelse(env$Location %in% Transition, "Transition",
                          ifelse(env$Location %in% Offshore, "Offshore", "neg_control")))
env$Environment <- ifelse(env$Location %in% Coast, "Coast",
                          ifelse(grepl("Track|WBB05|WBB06|WBB07", env$Location), "inside_OWF",
                                 ifelse(grepl("neg", env$Location), "neg_control", "outside_OWF")))
env$Zone_color <- ifelse(env$Zone=="Coast","limegreen",
                         ifelse(env$Zone=="Transition", "slateblue",
                                ifelse(env$Zone=="Offshore","darkorange","red")))
for (r in 1:nrow(env)){
  if (env$Zone[r] == "Transition" && env$Environment[r] == "outside_OWF"){
    env$Environment_color[r] <- "slateblue4"
  }
  else if (env$Zone[r] =="Offshore" && env$Environment[r] == "outside_OWF"){
      env$Environment_color[r] <- "darkorange3"
  }
  else{
    env$Environment_color[r] <- env$Zone_color[r]
  }
}
env <- env[!env$Niskin.sample %in% c("ftWT2bis_1234_3", "ftWT2bis_5", "ftWT2bis_6", 
                                     "ftWT2bis_8", "ftWT2bis_9", "ft230_1_bis", "ft230_2_bis"),]
saveRDS(env, file = paste0(proj.path,"/MiFish_UE-S_concatenated/results_microDecon/R_Environment/env_AllSamples.rds"))

#Reorder the dataframe: ASV_names (always first column), negative controls, true samples, Taxa (always last column)
table_raw_neg <- table_raw_noT[, grepl("neg", colnames(table_raw_noT))]
table_raw_sample <- table_raw_noT[, !colnames(table_raw_noT) %in% c(colnames(table_raw_neg))]
table_raw_noT_order <- cbind(table_raw$ASV, table_raw_neg, table_raw_sample, table_raw$Species)
colnames(table_raw_noT_order)[c(1,ncol(table_raw_noT_order))] <- c("ASV", "Taxa")

#Select samples from the DYFS-campaign in September
DYFS_samples <- c(env$Niskin.sample[which(env$Zone == "Coast" & !env$Location == "ft230")])
DYFS <- c(DYFS_samples, env$Niskin.sample[which(env$Location %in% c("neg_coastal", "neg_PCR"))], 
          "filter_5", "filter_6", "neg_DNA_3", "ASV", "Taxa")

table_DYFS <- table_raw_noT_order[, grepl(paste(DYFS, collapse = "|"), 
                                          colnames(table_raw_noT_order))]
table_DYFS_num <- table_DYFS[,2:(ncol(table_DYFS)-1)]
remove_DYFS <- c(colnames(table_DYFS_num[colSums(table_DYFS_num) == 0]))
table_DYFS <- table_DYFS[, !colnames(table_DYFS) %in% remove_DYFS]
neg_DYFS <- length(c(colnames(table_DYFS[, grepl("neg", colnames(table_DYFS))])))
Locations_DYFS <- unique(env$Location[which(env$Zone == "Coast" & !env$Location == "ft230")])
numb_DYFS <- NULL
for (i in 1:length(Locations_DYFS)){
  numb_DYFS[i] <- sum(str_count(gsub("\\_.*", "", colnames(table_DYFS_num)), paste0("\\b", Locations_DYFS[i], "\\b"))) 
}

decontaminated_DYFS <- decon(data = table_DYFS, numb.blanks=neg_DYFS, numb.ind=numb_DYFS, taxa=T,
                             runs=2,thresh = 0.7,prop.thresh = 0.00005,regression=0,low.threshold=40,up.threshold=400)

decon_DYFS <- as.data.frame(decontaminated_DYFS$decon.table)
contaminant_DYFS <- as.data.frame(decontaminated_DYFS$OTUs.removed)
contaminant_reads_DYFS <- as.data.frame(decontaminated_DYFS$reads.removed)

table_DYFS_Field <- table_DYFS[,!grepl("filter|DNA|PCR", colnames(table_DYFS))]
neg_DYFS_Field <- length(c(colnames(table_DYFS_Field[, grepl("neg", colnames(table_DYFS_Field))])))

decontaminated_DYFS_Field <- decon(data = table_DYFS_Field, numb.blanks=neg_DYFS_Field, numb.ind=numb_DYFS, taxa=T,
                                   runs=2,thresh = 0.7,prop.thresh = 0.00005,regression=0,low.threshold=40,up.threshold=400)

decon_DYFS_Field <- as.data.frame(decontaminated_DYFS_Field$decon.table)
contaminant_DYFS_Field <- as.data.frame(decontaminated_DYFS_Field$OTUs.removed)
contaminant_reads_DYFS_Field <- as.data.frame(decontaminated_DYFS_Field$reads.removed)

table_DYFS_NoField <- table_DYFS[,!grepl("neg_control", colnames(table_DYFS))]
neg_DYFS_NoField <- length(c(colnames(table_DYFS_NoField[, grepl("neg", colnames(table_DYFS_NoField))])))

decontaminated_DYFS_NoField <- decon(data = table_DYFS_NoField, numb.blanks=neg_DYFS_NoField, numb.ind=numb_DYFS, taxa=T,
                                   runs=2,thresh = 0.7,prop.thresh = 0.00005,regression=0,low.threshold=40,up.threshold=400)

decon_DYFS_NoField <- as.data.frame(decontaminated_DYFS_NoField$decon.table)
contaminant_DYFS_NoField <- as.data.frame(decontaminated_DYFS_NoField$OTUs.removed)
contaminant_reads_DYFS_NoField <- as.data.frame(decontaminated_DYFS_NoField$reads.removed)

#Select samples from the GeoV-campaign in November
Geo_samples <- c(env$Niskin.sample[which(env$Zone %in% c("Transition", "Offshore") | env$Location == "ft230")])
Geo <- c(Geo_samples, env$Niskin.sample[which(env$Location %in% c("neg_inside_OWF", "neg_outside_OWF", "neg_PCR"))], 
         "filter_1", "filter_2", "filter_3", "filter_4", 
         "neg_DNA_1", "neg_DNA_2", "ASV", "Taxa")

table_Geo <- table_raw_noT_order[, grepl(paste(Geo, collapse = "|"), 
                                          colnames(table_raw_noT_order))]
table_Geo <- table_Geo[, !grepl("_bis", colnames(table_Geo))]
table_Geo_num <- table_Geo[,2:(ncol(table_Geo)-1)]
remove_Geo <- c(colnames(table_Geo_num[colSums(table_Geo_num) == 0]))
table_Geo <- table_Geo[, !colnames(table_Geo) %in% remove_Geo]
table_Geo <- table_Geo %>% relocate("ftTrack5_4_S2", .before = 'ftTrack5_4_S3')
table_Geo <- table_Geo %>% relocate("ftWT1bis_4_S2", .before = 'ftWT1bis_4_S3')
neg_Geo <- length(c(colnames(table_Geo[, grepl("neg", colnames(table_Geo))])))
Locations_Geo <- unique(env$Location[which(env$Zone %in% c("Transition", "Offshore") | env$Location == "ft230")])
numb_Geo <- NULL
for (i in 1:length(Locations_Geo)){
  numb_Geo[i] <- sum(str_count(gsub("\\_.*", "", colnames(table_Geo_num)), paste0("\\b", Locations_Geo[i], "\\b"))) 
}

decontaminated_Geo <- decon(data = table_Geo, numb.blanks=neg_Geo, numb.ind=numb_Geo, taxa=T,
                             runs=2,thresh = 0.7,prop.thresh = 0.00005,regression=0,low.threshold=40,up.threshold=400)

decon_Geo <- as.data.frame(decontaminated_Geo$decon.table)
contaminant_Geo <- as.data.frame(decontaminated_Geo$OTUs.removed)
contaminant_reads_Geo <- as.data.frame(decontaminated_Geo$reads.removed)

table_Geo_Field <- table_Geo[,!grepl("filter|DNA|PCR", colnames(table_Geo))]
neg_Geo_Field <- length(c(colnames(table_Geo_Field[, grepl("neg", colnames(table_Geo_Field))])))

decontaminated_Geo_Field <- decon(data = table_Geo_Field, numb.blanks=neg_Geo_Field, numb.ind=numb_Geo, taxa=T,
                                   runs=2,thresh = 0.7,prop.thresh = 0.00005,regression=0,low.threshold=40,up.threshold=400)

decon_Geo_Field <- as.data.frame(decontaminated_Geo_Field$decon.table)
contaminant_Geo_Field <- as.data.frame(decontaminated_Geo_Field$OTUs.removed)
contaminant_reads_Geo_Field <- as.data.frame(decontaminated_Geo_Field$reads.removed)

table_Geo_NoField <- table_Geo[,!grepl("OWF", colnames(table_Geo))]
neg_Geo_NoField <- length(c(colnames(table_Geo_NoField[, grepl("neg", colnames(table_Geo_NoField))])))

decontaminated_Geo_NoField <- decon(data = table_Geo_NoField, numb.blanks=neg_Geo_NoField, numb.ind=numb_Geo, taxa=T,
                                  runs=2,thresh = 0.7,prop.thresh = 0.00005,regression=0,low.threshold=40,up.threshold=400)

decon_Geo_NoField <- as.data.frame(decontaminated_Geo_NoField$decon.table)
contaminant_Geo_NoField <- as.data.frame(decontaminated_Geo_NoField$OTUs.removed)
contaminant_reads_Geo_NoField <- as.data.frame(decontaminated_Geo_NoField$reads.removed)

table_clean <- merge(decon_DYFS, decon_Geo, by = 1, all = T)
table_clean <- table_clean[, !grepl("Taxa|blank", colnames(table_clean))]
table_raw_clean <- merge(table_clean, table_raw[(ncol(table_raw)-10):ncol(table_raw)],
                         by.x = 1, by.y = "ASV", all.x = T)
colnames(table_raw_clean)[1] <- "ASV"
table_raw_clean <- table_raw_clean %>% relocate("ASV", .before = 'DADA2')

table_clean_Field <- merge(decon_DYFS_Field , decon_Geo_Field, by = 1, all = T)
table_clean_Field <- table_clean_Field[, !grepl("Taxa|blank", colnames(table_clean_Field))]
table_raw_clean_Field <- merge(table_clean_Field, table_raw[(ncol(table_raw)-10):ncol(table_raw)],
                         by.x = 1, by.y = "ASV", all.x = T)
colnames(table_raw_clean_Field)[1] <- "ASV"
table_raw_clean_Field <- table_raw_clean_Field %>% relocate("ASV", .before = 'DADA2')

#reorder the environmental data
env_order<- env %>% 
  arrange(factor(Location, levels = c("neg_coastal", "neg_inside_OWF", "neg_outside_OWF", 
                                      "neg_filter", "neg_DNA", "neg_PCR")))
env_order <- env_order %>% arrange(factor(Zone, levels = c("Coast", "Transition", "Offshore", "neg_control")), 
                                    Environment)
saveRDS(env_order, file = paste0(proj.path,"/MiFish_UE-S_concatenated/results_microDecon/R_Environment/env_ordered_AllSamples.rds"))
saveRDS(env_order[which(!env_order$Zone == "neg_control"),], file = paste0(proj.path,"/MiFish_UE-S_concatenated/results_microDecon/R_Environment/env_ordered_noNeg.rds"))

#concatenate the samples
seqtab_raw <- as.data.frame(t(table_raw_clean[,1:(ncol(table_raw_clean)-11)]))
colnames(seqtab_raw) <- table_raw_clean$ASV
seqtab_raw$names <- str_sub(rownames(seqtab_raw), end=-4) # removes PCR number (_S1, _S2 or _S3) from samplenames
seqtab_concatenated <- aggregate(seqtab_raw[,1:ncol(seqtab_raw)-1], by= list(seqtab_raw$names),FUN=sum)
rownames(seqtab_concatenated) <- seqtab_concatenated$Group.1
seqtab_concatenated <- seqtab_concatenated[,2:ncol(seqtab_raw)]  
table_clean_concatenated_noT <- as.data.frame(t(seqtab_concatenated))
table_clean_concatenated_noT <- table_clean_concatenated_noT[, c(env_order$Niskin.sample[which(!env_order$Zone == "neg_control")])]
table_unrarefied_concatenated <- merge(table_clean_concatenated_noT, table_raw[(ncol(table_raw)-10):ncol(table_raw)],
                                       by.x = 0, by.y = "ASV", all.x = T)
colnames(table_unrarefied_concatenated)[1] <- "ASV"
table_unrarefied_concatenated <- table_unrarefied_concatenated %>% relocate("ASV", .before = 'DADA2')

seqtab_raw_Field <- as.data.frame(t(table_raw_clean_Field[,1:(ncol(table_raw_clean_Field)-11)]))
colnames(seqtab_raw_Field) <- table_raw_clean_Field$ASV
seqtab_raw_Field$names <- str_sub(rownames(seqtab_raw_Field), end=-4) # removes PCR number (_S1, _S2 or _S3) from samplenames
seqtab_concatenated_Field <- aggregate(seqtab_raw_Field[,1:ncol(seqtab_raw_Field)-1], by= list(seqtab_raw_Field$names),FUN=sum)
rownames(seqtab_concatenated_Field) <- seqtab_concatenated_Field$Group.1
seqtab_concatenated_Field <- seqtab_concatenated_Field[,2:ncol(seqtab_raw_Field)]  
table_clean_concatenated_noT_Field <- as.data.frame(t(seqtab_concatenated_Field))
table_clean_concatenated_noT_Field <- table_clean_concatenated_noT_Field[, c(env_order$Niskin.sample[which(!env_order$Zone == "neg_control")])]
table_unrarefied_concatenated_Field <- merge(table_clean_concatenated_noT_Field, table_raw[(ncol(table_raw)-10):ncol(table_raw)],
                                       by.x = 0, by.y = "ASV", all.x = T)
colnames(table_unrarefied_concatenated_Field)[1] <- "ASV"
table_unrarefied_concatenated_Field <- table_unrarefied_concatenated_Field %>% relocate("ASV", .before = 'DADA2')

zero_ASVs <- table_unrarefied_concatenated$ASV[rowSums(
  table_unrarefied_concatenated[,1:(ncol(table_unrarefied_concatenated)-11)], na.rm = T) == 0]

table_unrarefied_concatenated <- table_unrarefied_concatenated[!table_unrarefied_concatenated$ASV %in% zero_ASVs,]

#Save Datasets
write.xlsx(table_raw_clean, 
           paste0(proj.path,"/MiFish_UE-S_concatenated/results_microDecon/table_raw_FullTaxonomicAssignment_clean.xlsx"), 
           sheetName = "FullTaxAss_CleanedASVs", colNames = TRUE, rowNames = FALSE, append = FALSE)

write.xlsx(table_unrarefied_concatenated, 
           paste0(proj.path,"/MiFish_UE-S_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean.xlsx"), 
           sheetName = "FullTaxAss_CleanedASVs", colNames = TRUE, rowNames = FALSE, append = FALSE)

write.xlsx(table_unrarefied_concatenated_Field, 
           paste0(proj.path,"/MiFish_UE-S_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean_Field.xlsx"), 
           sheetName = "FullTaxAss_CleanedASVs", colNames = TRUE, rowNames = FALSE, append = FALSE)

list_of_datasets <- list("DYFS" = contaminant_DYFS, 
                         "DYFS_Field" = contaminant_DYFS_Field, 
                         "DYFS_noField" = contaminant_DYFS_NoField,
                         "DYFS_reads" = contaminant_reads_DYFS, 
                         "DYFS_reads_Field" = contaminant_reads_DYFS_Field,
                         "DYFS_reads_noField" = contaminant_reads_DYFS_NoField,
                         "GeoV" = contaminant_Geo,
                         "GeoV_Field" = contaminant_Geo_Field,
                         "GeoV_NoField" = contaminant_Geo_NoField,
                         "GeoV_reads" = contaminant_reads_Geo,
                         "GeoV_reads_Field" = contaminant_reads_Geo_Field,
                         "GeoV_reads_NoField" = contaminant_reads_Geo_NoField)
write.xlsx(list_of_datasets, paste0(proj.path,"/MiFish_UE-S_concatenated/results_microDecon/contaminantASVs.xlsx"), colNames = T)

