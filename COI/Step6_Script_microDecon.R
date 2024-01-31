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
devtools::install_github("donaldtmcknight/microDecon")

# make paths
proj.path <- here("/home/genomics/icornelis/02_ZEROimpact/02_COI/NJ2021")

# upload the dataset
env <- read.csv(paste0(proj.path,"/Step5_Statistics/environmental_data.csv"),  header=TRUE, sep=";")
table_raw <- readxl::read_excel(paste0(proj.path,"/OWFvsCoastal_concatenated/results/table_unrarefied_unconcatenated_FULLTaxonomicAssignment.xlsx"),
                                sheet = "table_unrarefied_unconcatenated")

# Remove the taxonomic assignment
table_raw_noT <- as.data.frame(table_raw[,1:(ncol(table_raw)-11)])

#Reorder the dataframe: ASV_names (always first column), negative controls, true samples, Taxa (always last column)
table_raw_neg <- table_raw_noT[, grepl("neg", colnames(table_raw_noT))]
table_raw_sample <- table_raw_noT[, !colnames(table_raw_noT) %in% c(colnames(table_raw_neg))]
table_raw_noT_order <- cbind(table_raw$ASV, table_raw_neg, table_raw_sample, table_raw$Full)
colnames(table_raw_noT_order)[c(1,ncol(table_raw_noT_order))] <- c("ASV", "Taxa")

#Select samples from the DYFS-campaign in September
DYFS_samples <- c(env$Niskin.sample[which(env$Zones == "Coastal" & !env$Location == "ft230")])
DYFS <- c(DYFS_samples, env$Niskin.sample[which(env$Location %in% c("neg_Coastal", "neg_PCR"))], 
          "filter_5", "filter_6", "neg_DNA_3", "ASV", "Taxa")

table_DYFS <- table_raw_noT_order[, grepl(paste(DYFS, collapse = "|"), 
                                          colnames(table_raw_noT_order))]
table_DYFS_num <- table_DYFS[,2:(ncol(table_DYFS)-1)]
remove_DYFS <- c(colnames(table_DYFS_num[colSums(table_DYFS_num) == 0]))
table_DYFS <- table_DYFS[, !colnames(table_DYFS) %in% remove_DYFS]
table_DYFS <- table_DYFS %>% relocate("37_2_DYF_NJ2021_S3", .before = "37_3_DYF_NJ2021_S1")
neg_DYFS <- length(c(colnames(table_DYFS[, grepl("neg", colnames(table_DYFS))])))
numb_DYFS <- c(rep(9, times = 6), 7, rep(9, times=4))

decontaminated_DYFS <- decon(data = table_DYFS, numb.blanks=neg_DYFS, numb.ind=numb_DYFS, taxa=T,
                             runs=2,thresh = 0.7,prop.thresh = 0.00005,regression=0,low.threshold=40,up.threshold=400)

decon_DYFS <- as.data.frame(decontaminated_DYFS$decon.table)
contaminant_DYFS <- as.data.frame(decontaminated_DYFS$OTUs.removed)

table_DYFS_Field <- table_DYFS[,!grepl("filter|DNA|PCR", colnames(table_DYFS))]
neg_DYFS_Field <- length(c(colnames(table_DYFS_Field[, grepl("neg", colnames(table_DYFS_Field))])))

decontaminated_DYFS_Field <- decon(data = table_DYFS_Field, numb.blanks=neg_DYFS_Field, numb.ind=numb_DYFS, taxa=T,
                                   runs=2,thresh = 0.7,prop.thresh = 0.00005,regression=0,low.threshold=40,up.threshold=400)

decon_DYFS_Field <- as.data.frame(decontaminated_DYFS$decon.table)
contaminant_DYFS_Field <- as.data.frame(decontaminated_DYFS$OTUs.removed)

Geo_samples <- c(env$Niskin.sample[which(env$Zones %in% c("zone1", "zone2") | env$Location == "ft230")])
Geo <- c(Geo_samples, env$Location[which(env$Location %in% c("neg_inside_OWF", "neg_outside_OWF", "neg_PCR"))], 
         "filter_1", "filter_2", "filter_3", "filter_4", 
         "neg_DNA_1", "neg_DNA_2", "ASV", "Taxa")

table_Geo <- table_raw_noT_order[, grepl(paste(Geo, collapse = "|"), 
                                         colnames(table_raw_noT_order))]
table_Geo <- table_Geo[, !grepl("_bis", colnames(table_Geo))]
table_Geo_num <- table_Geo[,2:(ncol(table_Geo)-1)]
remove_Geo <- c(colnames(table_Geo_num[colSums(table_Geo_num) == 0]))
table_Geo <- table_Geo[, !colnames(table_Geo) %in% remove_Geo]
neg_Geo <- length(c(colnames(table_Geo[, grepl("neg", colnames(table_Geo))])))
numb_Geo <- c(9, rep(15, times=7), 12, rep(15, times=5), 12, 15, 15, 15, 12)

for (i in 1:length(Geo_samples)) {
  numb_Geo[i] <- length(grep(Geo_samples[i], colnames(table_Geo)))
}

decontaminated_Geo <- decon(data = table_Geo, numb.blanks=neg_Geo, numb.ind=numb_Geo, taxa=T,
                            runs=2,thresh = 0.7,prop.thresh = 0.00005,regression=0,low.threshold=40,up.threshold=400)

decon_Geo <- as.data.frame(decontaminated_Geo$decon.table)
contaminant_Geo <- as.data.frame(decontaminated_Geo$OTUs.removed)

table_Geo_Field <- table_Geo[,!grepl("filter|DNA|PCR", colnames(table_Geo))]
neg_Geo_Field <- length(c(colnames(table_Geo_Field[, grepl("neg", colnames(table_Geo_Field))])))

decontaminated_Geo_Field <- decon(data = table_Geo_Field, numb.blanks=neg_Geo_Field, numb.ind=numb_Geo, taxa=T,
                                  runs=2,thresh = 0.7,prop.thresh = 0.00005,regression=0,low.threshold=40,up.threshold=400)

decon_Geo_Field <- as.data.frame(decontaminated_Geo$decon.table)
contaminant_Geo_Field <- as.data.frame(decontaminated_Geo$OTUs.removed)

table_clean <- merge(decon_DYFS, decon_Geo, by = 1, all = T)
table_clean <- table_clean[, !grepl("Taxa|blank", colnames(table_clean))]
table_raw_clean <- merge(table_clean, table_raw[(ncol(table_raw)-10):ncol(table_raw)],
                         by.x = 1, by.y = "ASV", all.x = T)
colnames(table_raw_clean)[1] <- "ASV"
table_raw_clean <- table_raw_clean %>% relocate("ASV", .before = 'blastn_GenBank')

#concatenate the samples
seqtab_raw <- as.data.frame(t(table_raw_clean[,1:(ncol(table_raw_clean)-11)]))
colnames(seqtab_raw) <- table_raw_clean$ASV
seqtab_raw$names <- str_sub(rownames(seqtab_raw), end=-4) # removes PCR number (_S1, _S2 or _S3) from samplenames
seqtab_concatenated <- aggregate(seqtab_raw[,1:ncol(seqtab_raw)-1], by= list(seqtab_raw$names),FUN=sum)
rownames(seqtab_concatenated) <- seqtab_concatenated$Group.1
seqtab_concatenated <- seqtab_concatenated[,2:ncol(seqtab_raw)]  
NumberOfSequencesConcat <- as.data.frame(sort.default(rowSums(seqtab_concatenated[1:nrow(seqtab_concatenated),]))) #number of sequences per sample after concatenation, 8 samples <300 reads
NumberOfASVsConcat <- as.data.frame(rowSums(seqtab_concatenated !=0)) #number of ASVs per sample after concatenation
#write_tsv(as_tibble(NumberOfASVsConcat5, rownames="asv"), file=paste0(proj.path,"/MiFish_UE-S_concatenated/results2/NumberOfASVsPerSample_Concat_unrarefied_CleanedByNoField_05Threshold.tsv"))
#write_tsv(as_tibble(NumberOfSequencesConcat5, rownames="samples"), file=paste0(proj.path, "/MiFish_UE-S_concatenated/results2/NumberOfSequencesPerSamples_Concat_unrarefied_CleanedByNoField_05Threshold.tsv"))
table_clean_concatenated_noT <- as.data.frame(t(seqtab_concatenated))
table_unrarefied_concatenated <- merge(table_clean_concatenated_noT, table_raw[(ncol(table_raw)-10):ncol(table_raw)],
                                       by.x = 0, by.y = "ASV", all.x = T)
colnames(table_unrarefied_concatenated)[1] <- "ASV"
table_unrarefied_concatenated <- table_unrarefied_concatenated %>% relocate("ASV", .before = 'blastn_GenBank')

write.xlsx(table_unrarefied_concatenated, 
           paste0(proj.path,"/OWFvsCoastal_concatenated/results_microDecon/table_unrarefied_concatenated_FullTaxonomicAssignment_clean.xlsx"), 
           sheetName = "FullTaxAss_CleanedASVs", colNames = TRUE, rowNames = FALSE, append = FALSE)

