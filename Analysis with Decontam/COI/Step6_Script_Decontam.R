#set OS type for paths
if(.Platform$OS.type == "unix"){home<-"/home/"} else{
  home<-"//192.168.236.131/"}

library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(decontam); packageVersion("decontam")
library(dplyr)
library(microbiome)
library(vegan)
library(openxlsx)
library(tibble)
library(tidyverse)
library(ade4)
library(here)

# make paths
proj.path <- here("/home/genomics/icornelis/02_ZEROimpact/02_COI/NJ2021")

# upload the dataset
table_raw <- readxl::read_excel(paste0(proj.path,"/OWFvsCoastal_concatenated/results/table_unrarefied_unconcatenated_FULLTaxonomicAssignment.xlsx"),sheet = "table_unrarefied_unconcatenated")

# We us the unconcatenated, unrarefied data for DECONTAM cleaning. Steps to take;
# 1) calculate relative abundance 
# 2) create ASV, taxonomy and metadata files
# 3) create phloseq object
# 4) calculate decontam with combined method
# 5) concatenate PCR replicates
# 6) remove contaminants
# 7) rarefy data

### 1) calculate relative abundance of ASVs per sample; ASVs on column, samples on row
#Create ASV table, rows ASVs; columns samples; decontam will calculate the realtive abundance, then it should be samples are rows, relative abundances are column
table_raw_noT <- as.data.frame(table_raw[,1:(ncol(table_raw)-11)])
rownames(table_raw_noT) <- table_raw$ASV
table_raw_noT_noZero <- as.data.frame(table_raw_noT[,!colSums(table_raw_noT[, 1:ncol(table_raw_noT)]) == 0, ])
seqtab_raw_noZero <- as.data.frame(t(table_raw_noT_noZero))
seqtab_raw_noField <- seqtab_raw_noZero %>% filter(!grepl("neg_control", rownames(seqtab_raw_noZero)))
seqtab_raw_noField <- seqtab_raw_noField %>% filter(!grepl("OWF", rownames(seqtab_raw_noField)))
table_raw_noField <- as.data.frame(t(seqtab_raw_noField))

### 2) create ASV, taxonomy and metadata files
##ASV and Taxonomy table should be matrix for phyloseq
#Taxonomy table
table_raw_taxa <- as.data.frame(table_raw[, (ncol(table_raw) - 10):ncol(table_raw)])
rownames(table_raw_taxa) <- table_raw_taxa$ASV #you need to define rownames for phyloseq
table_raw_taxa_matrix <- as.matrix(table_raw_taxa)

#Create metadata table for phyloseq object
meta_raw <- as.data.frame(rownames(seqtab_raw_noZero))
colnames(meta_raw) <- "Sample"
rownames(meta_raw) <- meta_raw$Sample
meta_raw$Batch <- ifelse(grepl("_S1", meta_raw$Sample), "Batch1", 
                         ifelse(grepl("_S2", meta_raw$Sample), "Batch2", "Batch3"))
meta_raw$Sample_or_Control <- ifelse(grepl("neg_", meta_raw$Sample), "Control", "Sample")
meta_raw$is.neg <- meta_raw$Sample_or_Control== "Control"
meta_raw_noField <- meta_raw %>% filter(!grepl("neg_control", meta_raw$Sample))
meta_raw_noField <- meta_raw_noField %>% filter(!grepl("OWF", meta_raw_noField$Sample))

### 3) create the phyloseq object
phylo_raw_noField <- phyloseq(otu_table(seqtab_raw_noField, taxa_are_rows=FALSE),
                              sample_data(meta_raw_noField),
                              tax_table(table_raw_taxa_matrix))
               
phylo_raw <- phyloseq(otu_table(seqtab_raw_noZero, taxa_are_rows=FALSE),
                              sample_data(meta_raw),
                              tax_table(table_raw_taxa_matrix))

#inspect data to check for contaminants
df <- as.data.frame(sample_data(phylo_raw_noField)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phylo_raw_noField)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

df <- as.data.frame(sample_data(phylo_raw)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phylo_raw)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()


### 4) calculate decontam with prevalence method and the field negatives excluded
contamination_table_raw_prev <- isContaminant(phylo_raw_noField, method = "prevalence", 
                                                   neg = "is.neg", 
                                                   threshold = 0.5, #change the treshold to find the best division between true ASVs and contaminant ASVs
                                                   batch="Batch",
                                                   normalize = TRUE, # default TRUE, if FALSE data should be counts from equal-depth samples
                                                   detailed = TRUE) # if TRUE the return value is a data.frame

# plot histogram and prevalence of ASVs in positive vs negative samples
hist(contamination_table_raw_prev$p.prev, main="Histogram Prevalence Treshold 0.5 No Field negatives", xlab="p.prev values")

ps.pa_noField <- transform_sample_counts(phylo_raw_noField, function(abund) 1*(abund>0))
ps.pa.neg_noField <- prune_samples(sample_data(ps.pa_noField)$Sample_or_Control == "Control", ps.pa_noField)
ps.pa.pos_noField <- prune_samples(sample_data(ps.pa_noField)$Sample_or_Control == "Sample", ps.pa_noField)
df.pa_noField <- data.frame(pa.pos=taxa_sums(ps.pa.pos_noField), pa.neg=taxa_sums(ps.pa.neg_noField),
                            contaminant=contamination_table_raw_prev$contaminant)
ggplot(data=df.pa_noField, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#THRESHOLD 0.4 prevalence and field negatives removed
# Taxonomic assignments of contaminants
assign_freq <- subset_taxa(phylo_raw_noField, contamination_table_raw_prev$contaminant==TRUE)
contam_taxa <- data.frame(assign_freq@tax_table)

#The prevalence method with a threshold of 0.5 was without the Fiel negatives included was chosen based on the plots 
write.xlsx(as.data.frame(contam_taxa),
           paste0(proj.path, "/OWFvsCoastal_concatenated/results/Decontam/FullContaminationTable_raw_DecontamPrevalence_NoFieldNeg_Threshold05.xlsx"), sheetName= "Contamination Table", colNames = TRUE, rowNames = TRUE, append = FALSE)

### 4) calculate decontam with prevalence method and the field negatives included
contamination_table_raw_prev_2 <- isContaminant(phylo_raw, method = "prevalence", 
                                              neg = "is.neg", 
                                              threshold = 0.5, #change the treshold to find the best division between true ASVs and contaminant ASVs
                                              batch="Batch",
                                              normalize = TRUE, # default TRUE, if FALSE data should be counts from equal-depth samples
                                              detailed = TRUE) # if TRUE the return value is a data.frame

# plot histogram and prevalence of ASVs in positive vs negative samples
hist(contamination_table_raw_prev_2$p.prev, main="Histogram Prevalence Treshold 0.5 With Field negatives", xlab="p.prev values")

ps.pa_Field <- transform_sample_counts(phylo_raw, function(abund) 1*(abund>0))
ps.pa.neg_Field <- prune_samples(sample_data(ps.pa_Field)$Sample_or_Control == "Control", ps.pa_Field)
ps.pa.pos_Field <- prune_samples(sample_data(ps.pa_Field)$Sample_or_Control == "Sample", ps.pa_Field)
df.pa_Field <- data.frame(pa.pos=taxa_sums(ps.pa.pos_Field), pa.neg=taxa_sums(ps.pa.neg_Field),
                            contaminant=contamination_table_raw_prev_2$contaminant)
ggplot(data=df.pa_Field, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#THRESHOLD 0.6 prevalence and field negatives removed
# Taxonomic assignments of contaminants
assign_freq_2 <- subset_taxa(phylo_raw, contamination_table_raw_prev_2$contaminant==TRUE)
contam_taxa_Field <- data.frame(assign_freq_2@tax_table)

#contaminants assigned to fish taxa
contam_taxa_Field[is.na(contam_taxa_Field)] <- "NA" #125 taxa
contam_taxa_fish_Field <- as.data.frame(contam_taxa_Field[contam_taxa_Field$Class %in% fish_classes,])

#The prevalence method with a threshold of 0.5 was without the Field negatives included was chosen based on the plots 
write.xlsx(as.data.frame(contam_taxa_Field),
           paste0(proj.path, "/OWFvsCoastal_concatenated/results/Decontam/FullContaminationTable_raw_DecontamPrevalence_WithFieldNeg_Threshold05.xlsx"), sheetName= "Contamination Table", colNames = TRUE, rowNames = TRUE, append = FALSE)


### 5) concatenate PCR replicates
#concatenate PCR replicates from the unrarefied samples 
seqtab_raw <- as.data.frame(t(table_raw_noT))
seqtab_raw$names <- str_sub(rownames(seqtab_raw), end=-4) # removes PCR number (_S1, _S2 or _S3) from samplenames
seqtab_concatenated <- aggregate(seqtab_raw[,1:ncol(seqtab_raw)-1], by= list(seqtab_raw$names),FUN=sum)
rownames(seqtab_concatenated) <- seqtab_concatenated$Group.1
seqtab_concatenated <- seqtab_concatenated[,2:ncol(seqtab_raw)]  
NumberOfSequencesConcat <- as.data.frame(sort.default(rowSums(seqtab_concatenated[1:nrow(seqtab_concatenated),]))) #number of sequences per sample after concatenation, 8 samples <300 reads
NumberOfASVsConcat <- as.data.frame(rowSums(seqtab_concatenated !=0)) #number of ASVs per sample after concatenation
#write_tsv(as_tibble(NumberOfASVsConcat5, rownames="asv"), file=paste0(proj.path,"/MiFish_UE-S_concatenated/results2/NumberOfASVsPerSample_Concat_unrarefied_CleanedByNoField_05Threshold.tsv"))
#write_tsv(as_tibble(NumberOfSequencesConcat5, rownames="samples"), file=paste0(proj.path, "/MiFish_UE-S_concatenated/results2/NumberOfSequencesPerSamples_Concat_unrarefied_CleanedByNoField_05Threshold.tsv"))
table_unrarefied_concatenated_noT <- as.data.frame(t(seqtab_concatenated))
table_unrarefied_concatenated <- merge(table_unrarefied_concatenated_noT, table_raw_taxa, by = "row.names")
rownames(table_unrarefied_concatenated) <- table_unrarefied_concatenated[, 1]
table_unrarefied_concatenated[, 1] <- NULL

### 7) Remove contaminants
ASV_contaminants <- rownames(contamination_table_raw_prev_2[which(contamination_table_raw_prev_2$contaminant == TRUE),])
table_raw_clean <- table_raw[!table_raw$ASV %in% ASV_contaminants,]
table_unrarefied_concatenated_clean <- table_unrarefied_concatenated[!table_unrarefied_concatenated$ASV %in% ASV_contaminants,]
seqtab_concatenated_clean <- as.data.frame(t(table_unrarefied_concatenated_clean[,1:(ncol(table_unrarefied_concatenated_clean)-11)]))

### 6) rarefy data
#create environmental data
env_unrarefied <- as.data.frame(rownames(seqtab_concatenated))
colnames(env_unrarefied) <- 'Niskin.sample'
Coast <- readRDS(file = here("/home/genomics/icornelis/03_RawScripts/Locations_Coast.rds"))
Coast <- c(Coast, "1", "13", "20", "23", "24", "3", "37", "38", "6", "92", "93")
Transition <- readRDS(file = here("/home/genomics/icornelis/03_RawScripts/Locations_Transition.rds"))
Transition <- c(Transition, "ftWT2bis")
Offshore <- readRDS(file = here("/home/genomics/icornelis/03_RawScripts/Locations_Offshore.rds"))
env_unrarefied$Location <- gsub("\\_.*", "", env_unrarefied$Niskin.sample)
env_unrarefied$Biol_replicate <- gsub("^.*\\_", "", env_unrarefied$Niskin.sample)
env_unrarefied$Zone <- ifelse(env_unrarefied$Location %in% Coast, "Coast",
                              ifelse(env_unrarefied$Location %in% Transition, "Transition",
                                     ifelse(env_unrarefied$Location %in% Offshore, "Offshore", "neg_control")))
env_unrarefied$Environment_color <- ifelse(env_unrarefied$Zone=="Coast","limegreen",
                                           ifelse(env_unrarefied$Zone=="Transition", "slateblue",
                                                  ifelse(env_unrarefied$Zone=="Offshore","darkorange","red")))
saveRDS(env_unrarefied, file = paste0(proj.path,"/OWFvsCoastal_concatenated/results/Decontam/REnvironment/env_unrarefied.rds"))

#CREATE RARECURVES
raremax <- 25000
rarecurve(seqtab_concatenated_clean, ylab = "All ASVs", main = "Rarecurve of unrarefied samples after taxonomic assignment", col = as.vector(env_unrarefied$Environment_color), label = FALSE, step =1000)
rarecurve(seqtab_concatenated, ylab = "All ASVs", main = "Rarecurve of unrarefied samples before taxonomic assignment", col = as.vector(env_unrarefied$Environment_color), label = T, step =1000)
abline(v = raremax)
#title(sub = "Negative Field Controls are excluded to define the contaminants, threshold=0.4", line = -26)

#ASVs assigned to invertebrate taxa
table_unrarefied_concatenated_Inv <- as.data.frame(table_unrarefied_concatenated_clean[table_unrarefied_concatenated_clean$Kingdom %in% c("Animalia"),])
table_unrarefied_concatenated_Inv <- as.data.frame(table_unrarefied_concatenated_Inv[!table_unrarefied_concatenated_Inv$Full %in% "NA",])
#table_unrarefied_concatenated_Inv <- as.data.frame(table_unrarefied_concatenated_Inv[!table_unrarefied_concatenated_Inv$Full %in% c("NA","Clytia hemisphaerica"),])
table_unrarefied_concatenated_Inv <- as.data.frame(table_unrarefied_concatenated_Inv[!table_unrarefied_concatenated_Inv$Phylum %in% "Chordata",])

#remove samples from environmental data that are now absent from the data3
seqtab_Inv <- as.data.frame(t(table_unrarefied_concatenated_Inv[,1:(ncol(table_unrarefied_concatenated_Inv)-11)]))
#Put samples in same order as environmental data
seqtab_Inv <- seqtab_Inv[env_unrarefied$Niskin.sample,] #to match the colors correctly the samples in seqtab.concatenated must be in the same order as your environmental data
rarecurve(seqtab_Inv, ylab = "Invertebrate ASVs", main = "Rarecurve of unrarefied samples after taxonomic assignment", col = as.vector(env_unrarefied$Environment_color), label = T, step =100)
abline(v = raremax)

#rarefy the data and remove (reads < raremax), from cleaned and raw data
NotRarefiedConcat <- c(rownames(seqtab_concatenated_clean[rowSums(seqtab_concatenated_clean)<raremax,]))
seqtab_rarefied_clean <- as.data.frame(rrarefy(seqtab_concatenated_clean[1:nrow(seqtab_concatenated_clean),],raremax))
seqtab_rarefied_clean <- seqtab_rarefied_clean[,!colSums(seqtab_rarefied_clean) == 0]

table_rarefied_clean_allSamples<- as.data.frame(t(seqtab_rarefied_clean))
table_rarefied_clean_allSamples  <- merge(table_rarefied_clean_allSamples, table_raw_taxa, by= "row.names")
rownames(table_rarefied_clean_allSamples) <- table_rarefied_clean_allSamples[,1]
table_rarefied_clean_allSamples [,1] <- NULL
table_rarefied_clean <- table_rarefied_clean_allSamples[,!colnames(table_rarefied_clean_allSamples) %in% NotRarefiedConcat]
table_rarefied_clean <- table_rarefied_clean[!rowSums(table_rarefied_clean[, 1:(ncol(table_rarefied_clean) - 11)]) == 0,]

#you have to save here; we will use the cleaned and rarefied file with all samples still present for the heat map.
write.xlsx(as.data.frame(table_unrarefied_concatenated), 
           paste0(proj.path,"/OWFvsCoastal_concatenated/results/Decontam/table_unrarefied_AllASVs_FullTaxonomicAssignment.xlsx"), 
           sheetName = "FullTaxAss_AllASVs", colNames = TRUE, rowNames = FALSE, append = FALSE)
write.xlsx(as.data.frame(table_raw_clean), 
           paste0(proj.path,"/OWFvsCoastal_concatenated/results/table_raw_CleanedASVs_FullTaxonomicAssignment_NoField.xlsx"), 
           sheetName = "FullTaxAss_CleanedASVs_raw", colNames = TRUE, rowNames = FALSE, append = FALSE)
write.xlsx(as.data.frame(table_unrarefied_concatenated_clean), 
           paste0(proj.path,"/OWFvsCoastal_concatenated/results/Decontam/table_unrarefied_concatenated_CleanedASVs_FullTaxonomicAssignment_NoField.xlsx"), 
           sheetName = "FullTaxAss_CleanedASVs", colNames = TRUE, rowNames = FALSE, append = FALSE)
write.xlsx(as.data.frame(table_rarefied_clean_allSamples), 
           paste0(proj.path,"/OWFvsCoastal_concatenated/results/Decontam/table_rarefied_CleanedASVs_FullTaxonomicAssignment_allSamples_NoField.xlsx"), 
           sheetName = "FullTaxAss_CleanedASVs_rarefied", colNames = TRUE, rowNames = FALSE, append = FALSE)
write.xlsx(as.data.frame(table_rarefied_clean), 
           paste0(proj.path,"/OWFvsCoastal_concatenated/results/Decontam/table_rarefied_CleanedASVs_FullTaxonomicAssignment_NoField.xlsx"), 
           sheetName = "FullTaxAss_CleanesASVs_rarefied", colNames = TRUE, rowNames = FALSE, append = FALSE)

#keep the samples in env data according to rarefied data for rarecurve
keep_sampels_rarefied <- c(colnames(table_rarefied_clean))
env_rarefied <- env_unrarefied[env_unrarefied$Niskin.sample %in% keep_sampels_rarefied,]
seqtab_rarefied_clean <- seqtab_rarefied_clean[env_rarefied$Niskin.sample,]
rarecurve(seqtab_rarefied_clean, ylab = "All ASVs", 
          main = "Rarecurve of Rarefied samples after taxonomic assignment", 
          col = as.vector(env_rarefied$Environment_color), label = FALSE, step =100)

#create the rarecurve for invertebrate species
#ASVs assigned to invertebrate taxa
#table_rarefied_cleaned[is.na(table_rarefied_cleaned)] <- "NA"
table_rarefied_clean_Inv <- as.data.frame(table_rarefied_clean[table_rarefied_clean$Kingdom %in% c("Animalia"),])
table_rarefied_clean_Inv <- as.data.frame(table_rarefied_clean_Inv[!table_rarefied_clean_Inv$Full %in% "NA",])
table_rarefied_clean_Inv <- as.data.frame(table_rarefied_clean_Inv[!table_rarefied_clean_Inv$Phylum %in% "Chordata",])
seqtab_rarefied_clean_Inv <- as.data.frame(t(table_rarefied_clean_Inv[,1:(ncol(table_rarefied_clean_Inv)-11)]))
rarecurve(seqtab_rarefied_clean_Inv, ylab = "Fish ASVs", 
          main = "Rarecurve of Rarefied samples after taxonomic assignment", 
          col = as.vector(env_rarefied$Environment_color), label = FALSE, step =100)
