# make paths
proj.path <- here("/home/genomics/icornelis/02_ZEROimpact/02_COI/NJ2021")

# upload the dataset (rarfied data with samples and taxonomic assignment), blastn results from GenBank and blastn results from own reference list
table_unrarefied <- readxl::read_excel(paste0(proj.path,"/OWFvsCoastal_concatenated/results/table_unrarefied_concatenated.xlsx"),sheet = "table_unrarefied_concatenated")
table_rarefied_raw <- readxl::read_excel(paste0(proj.path,"/OWFvsCoastal_concatenated/results/Raremax_25000/table_rarefied_concatenated.xlsx"),sheet = "table_rarefied_concatenated")
#table_FullAssi <- readxl::read_excel(paste0(proj.path,"/OWFvsCoastal_concatenated/results/Raremax_25000/table_rarefied_FullTaxonomicAssignment_25000.xlsx"), sheet = "table_rarefied_FullTaxonomicAss")
blastn_Genbank <- readxl::read_excel(paste0(proj.path,"/Step4_blastn/blastn_GenBank_20221019.xlsx"),sheet = "blastn_GenBank_20221019")
blastn_ref <- readxl::read_excel(paste0(proj.path,"/Step4_blastn/blastn_own_references_tophits_20221019.xlsx"), sheet="blastn_own_reference_tophits_20")
asv <- read.table(file=paste0(proj.path,"/OWFvsCoastal_concatenated/results/asvs.tsv"), sep = '\t', header = TRUE)

# in the blastn data select ASV and species column 
TA_DADA2 <- data.frame(table_unrarefied$ASV, table_unrarefied$Species) 
blastn_Genbank_select <- blastn_Genbank[!blastn_Genbank$`pident` < 97000,]
blastn_Genbank <- data.frame(blastn_Genbank_select$ASV, blastn_Genbank_select$Species, blastn_Genbank_select$pident, blastn_Genbank_select$qcovs) 
blastn_ref_select <- blastn_ref[!blastn_ref$`pident` < 97000,]
blastn_ref <- data.frame(blastn_ref_select$ASV, blastn_ref_select$Species, blastn_ref_select$pident, blastn_ref_select$qcovs)
asv_seq <- data.frame(asv$names, asv$dnas)
#Rarefied_raw <- data.frame(table_rarefied_raw[,1:(ncol(table_rarefied_raw)-8)])
#Rarefied_raw$ASV <- table_rarefied_raw$ASV
#FullAss <- data.frame(table_FullAssi[,(ncol(table_FullAssi)-10):ncol(table_FullAssi)])

#merge 3 datasets according to ASV name
table_Species <- merge(TA_DADA2, blastn_Genbank , by.x='table_unrarefied.ASV', by.y='blastn_Genbank_select.ASV', all.x=TRUE)
table_Species_2 <- merge(table_Species, blastn_ref, by.x='table_unrarefied.ASV', by.y='blastn_ref_select.ASV', all.x=TRUE)
colnames(table_Species_2) <- c("ASV", "DADA2", "blastn_GenBank", "pident", "qcovs", "blastn_ref", "pident_ref", "qcovs_ref")
table_sequences <- merge(asv_seq, blastn_Genbank , by.x='asv.names', by.y='blastn_Genbank_select.ASV', all.x=TRUE)
write_tsv(as_tibble(table_sequences), file=paste0(proj.path,"/OWFvsCoastal_concatenated/results/blastn_Genbank.tsv"))
#table_Species_3 <- table_Species_2[!duplicated(table_Species_3$table_Species_2.ASV),]

#write full taxonomic assignment in seperate column Full
table_Species_2$Full <- table_Species_2$DADA2 #remove GenBAnk identifier from species name
#table_Species_2$blastn_ref <- str_extract(table_Species_2$blastn_ref, "^.*(?=(_))") #remove GenBank identifier from species name
table_Species_2[is.na(table_Species_2)] <- "NA"
table_Species_2$Full  <- ifelse(table_Species_2$DADA2=="NA", table_Species_2$blastn_ref, table_Species_2$Full) #add blastn results against custom reference database
table_Species_2$Full  <- ifelse(table_Species_2$Full=="Skeletonema pseudocostatum", table_Species_2$blastn_GenBank, 
                                       ifelse(table_Species_2$Full=="NA", table_Species_2$blastn_GenBank, table_Species_2$Full)) #add blastn results against GenBank
table_Species_2$Full <- word(table_Species_2$Full, 1,2, sep=" ") #keep only the first two words of the assigned species
table_Species_2[is.na(table_Species_2)] <- "NA"

#remove species with identical COI sequences 
table_Species_2$Full  <- ifelse(table_Species_2$Full=="Azadinium poporum", "NA", 
                                ifelse(table_Species_2$Full=="Azadinium obesum", "NA",
                                       ifelse(table_Species_2$Full=="Azadinium dalianense", "NA", 
                                              ifelse(table_Species_2$Full=="Trieres chinensis", "NA",
                                                     ifelse(table_Species_2$Full=="Leptocylindrus danicus", "NA", table_Species_2$Full)))))
table_Species_2$Full  <- ifelse(table_Species_2$Full=="Thoracosphaera heimii", "NA", 
                                ifelse(table_Species_2$Full=="Alexandrium tamiyavanichii", "NA",
                                       ifelse(table_Species_2$Full=="Alexandrium pseudogonyaulax", "NA", 
                                              ifelse(table_Species_2$Full=="Alexandrium catenella", "NA", table_Species_2$Full))))
table_Species_2$Full  <- ifelse(table_Species_2$Full=="Saprolegnia mixta", "NA", 
                                ifelse(table_Species_2$Full=="Saprolegnia delica", "NA",
                                       ifelse(table_Species_2$Full=="Saprolegnia bulbosa", "NA", table_Species_2$Full)))

write_tsv(as_tibble(table_Species_2), file=paste0(proj.path,"/OWFvsCoastal_concatenated/results/table_Species_2.tsv"))

#merge table_rarefied with final taxonomic assignment
#table_rarefied_raw$blastn_GenBank <- table_Species_2$blastn_GenBank
#table_rarefied_raw$blastn_ref <- table_Species_2$blastn_ref
#table_rarefied_raw$Full <- table_Species_2$Full
table_Species_3 <- merge(table_rarefied_raw, table_Species_2, by.x='ASV', by.y='ASV', all.x=TRUE)

#samples to be removed (< 25000 reads)
table_rarefied_noTA <- as.data.frame(table_rarefied_raw[,1:(ncol(table_rarefied_raw)-11)])
raremax <- 25000
rarefied <- c(colnames(table_rarefied_noTA[,colSums(table_rarefied_noTA) == raremax]))
Notrarefied <- c(colnames(table_rarefied_noTA[,colSums(table_rarefied_noTA) < raremax]))
table_rarefied <- table_rarefied_raw[,!colnames(table_rarefied_raw) %in% Notrarefied]
rownames(table_rarefied) <- table_rarefied_raw$ASV

#samples to be removed (< 25000 reads)
#table_rarefied_noTA <- as.data.frame(table_Species[,2:(ncol(table_Species)-10)])
#raremax <- 25000
#rarefied <- c(colnames(table_rarefied_noTA[,colSums(table_rarefied_noTA) == raremax]))
#Notrarefied <- c(colnames(table_rarefied_noTA[,colSums(table_rarefied_noTA) < raremax]))
#table_rarefied <- table_Species[,!colnames(table_Species) %in% Notrarefied]
#rownames(table_rarefied) <- table_rarefied_raw$ASV

#clean the rarefied table --> remove ASVs without reads after rarefaction, and with more than 10% presence in the negative control samples
neg_ctrl <- table_rarefied_raw[, grepl("neg", colnames(table_rarefied_raw))]
relative_abundance <- data.frame(matrix(nrow=nrow(neg_ctrl), ncol=ncol(neg_ctrl)
                                        , dimnames = list(c(rownames(neg_ctrl)), c(colnames(neg_ctrl)))))
for(i in 1:ncol(neg_ctrl)){relative_abundance[,i] <- (neg_ctrl[,i]/rowSums(table_rarefied_raw[1:(ncol(table_rarefied_raw)-11)]))*100}
rownames(relative_abundance) <- table_rarefied_raw$ASV
keep_ASVs <- c(rownames(relative_abundance[!rowSums(relative_abundance > 10),]))
#table_rarefied_raw_cleaned <- table_rarefied_raw[!rowSums(table_rarefied_raw[,1:(ncol(table_rarefied_raw)-11)])==0,]
#table_rarefied_raw_cleaned <- table_rarefied_raw_cleaned[table_rarefied_raw_cleaned$ASV %in% keep_ASVs,]
table_rarefied_raw_cleaned <- table_Species_3[!rowSums(table_Species_3[,2:(ncol(table_Species_3)-10)])==0,]
table_rarefied_raw_cleaned <- table_rarefied_raw_cleaned[table_rarefied_raw_cleaned$ASV %in% keep_ASVs,]
#table_rarefied_cleaned <- table_rarefied[!rowSums(table_rarefied[,1:(ncol(table_rarefied)-11)])==0,]
#table_rarefied_cleaned <- table_rarefied_cleaned[table_rarefied_cleaned$ASV %in% keep_ASVs,]
table_rarefied_cleaned <- table_rarefied[!rowSums(table_rarefied[,2:(ncol(table_rarefied)-10)])==0,]
table_rarefied_cleaned <- table_rarefied_cleaned[table_rarefied_cleaned$ASV %in% keep_ASVs,]
table_rarefied_contaminants <- table_rarefied_cleaned[!table_rarefied_cleaned$ASV %in% keep_ASVs,]
write_tsv(as_tibble(table_Species_3), 
          file=paste0(proj.path,"/OWFvsCoastal_concatenated/results/Raremax_25000/table_rarefied_AllASVs_FullTaxonomicAssignment_raw_97percent.tsv"))
write_tsv(as_tibble(table_rarefied_raw_cleaned), 
          file=paste0(proj.path,"/OWFvsCoastal_concatenated/results/Raremax_25000/table_rarefied_CleanedASVs_FullTaxonomicAssignment_raw_97percent.tsv"))
write_tsv(as_tibble(table_rarefied), 
          file=paste0(proj.path,"/OWFvsCoastal_concatenated/results/Raremax_25000/table_rarefied_AllASVs_FullTaxonomicAssignment_97percent.tsv"))
write_tsv(as_tibble(table_rarefied_cleaned), 
          file=paste0(proj.path,"/OWFvsCoastal_concatenated/results/Raremax_25000/table_rarefied_CleanedASVs_FullTaxonomicAssignment_97percent.tsv"))
write.xlsx(as.data.frame(table_rarefied_contaminants), paste0(proj.path,"/OWFvsCoastal_concatenated/results/table_rarefied_FullTaxonomicAssignment_contaminants.xlsx"))

#merge table_unrarefied with final taxonomic assignment
table_unrarefied$blastn_GenBank <- table_Species_2$blastn_GenBank
table_unrarefied$blastn_ref <- table_Species_2$blastn_ref
table_unrarefied$Full <- table_Species_2$Full
table_unrarefied_clean <- table_unrarefied[table_unrarefied$ASV %in% keep_ASVs,]
write.xlsx(as.data.frame(table_unrarefied_clean), paste0(proj.path,"/OWFvsCoastal_concatenated/results/table_unrarefied_FullTaxonomicAssignment.xlsx"))
write.xlsx(as.data.frame(table_rarefied_clean), paste0(proj.path,"/OWFvsCoastal_concatenated/results/table_unrarefied_FullTaxonomicAssignment_cleaned.xlsx"))

#Calculations
NumOfReads <- data.frame(colSums(table_unrarefied[10:43]))
meanReadsPerSample <- mean(NumOfReads$colSums.table_unrarefied.10.43..)
sdReadsPerSample <- sd(NumOfReads$colSums.table_unrarefied.10.43..)
NumOfReadsNov <- data.frame(colSums(table_unrarefied[c(1:5,47:76,89:150)]))
meanReadsPerSampleNov <- mean(NumOfReadsNov$colSums.table_unrarefied.c.1.5..47.76..89.150...)
sdReadsPerSampleNov <- se(NumOfReadsNov$colSums.table_unrarefied.c.1.5..47.76..89.150...)
f.x <- dnorm(NumOfReadsNov$colSums.table_unrarefied.c.1.5..47.76..89.150..., mean = 0, sd = 0)
plot(NumOfReadsNov$colSums.table_unrarefied.c.1.5..47.76..89.150..., f.x, type='l')
