#title: "Run DADA2 on demultiplexed MiSeq sequencing data"
#author: "Isolde Cornelis"
#date: "01/09/2022"

################################################################################

#set OS type for paths
if(.Platform$OS.type == "unix"){home<-"/home/"} else{
  home<-"//192.168.236.131/"}

#!/usr/bin/env Rscript
library(dada2)
library(vegan)
library(seqRFLP)
library(dplyr)
library(tibble)
library(stringr)


# report
writeLines("\nRunning dada2 denoising ...\n")

############## SET PARAMS ##############
############## SET PARAMS ##############

# load libs and funs
source("/home/genomics/icornelis/03_RawScripts/funs-libs.R")

# get args
option_list <- list( 
    make_option(c("-p","--primer"), type="character"),
    make_option(c("-l","--lib"), type="character")
    )

# set args
opt <- parse_args(OptionParser(option_list=option_list,add_help_option=FALSE))

# make paths
proj.path <- here("/home/genomics/icornelis/02_ZEROimpact/02_COI/NJ2021",paste0("OWFvsCoastal_concatenated"))
dir.sense <- paste0(proj.path,"/sense")
dir.antisense <- paste0(proj.path,"/antisense")

mkdir <- paste0(proj.path,"/results")

# confirm proj path 
writeLines(paste0("\n...\nOutput directory set to:\n",proj.path,"\n"))

# trucLens
# for Miya MiFish - truncLen 105 gives at least 29 bp overlap for the longest amplicons (e.g. Raja clavata @ 181 bp), and 40 bp for the regular 170 bp
trucVal <- c(105,105)


############## QUALITY TRIM TRUNCATE ##############
############## QUALITY TRIM TRUNCATE ##############

# report
writeLines("\n...\nQuality trimming and truncating\n")
Sys.sleep(3)

# quality trim Ns and truncate
#out.sense <- filterAndTrim(fwd=cpath("sense","trimmed","R1"), filt=cpath("sense","filtered","R1"), rev=cpath("sense","trimmed", "R2"), filt.rev=cpath("sense","filtered", "R2"), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=trucVal, multithread=TRUE, verbose=TRUE, matchIDs=TRUE)
#out.antisense <- filterAndTrim(fwd=cpath("antisense","trimmed","R1"), filt=cpath("antisense","filtered","R1"), rev=cpath("antisense","trimmed","R2"), filt.rev=cpath("antisense","filtered","R2"), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen=trucVal, multithread=TRUE, verbose=TRUE, matchIDs=TRUE)
out.sense <- filterAndTrim(fwd=cpath("sense","trimmed","R1"), filt=cpath("sense","filtered","R1"), rev=cpath("sense","trimmed", "R2"), filt.rev=cpath("sense","filtered", "R2"), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, multithread=TRUE, verbose=TRUE, matchIDs=TRUE)
out.antisense <- filterAndTrim(fwd=cpath("antisense","trimmed","R1"), filt=cpath("antisense","filtered","R1"), rev=cpath("antisense","trimmed","R2"), filt.rev=cpath("antisense","filtered","R2"), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, multithread=TRUE, verbose=TRUE, matchIDs=TRUE)


############## LEARN ERRORS ##############
############## LEARN ERRORS ##############

# report
writeLines("\n...\nLearning errors\n")
Sys.sleep(3)

# learn errors
set.seed(42)
sense.filt.R1.errs <- learnErrors(cpath("sense","filtered","R1"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)
set.seed(42)
sense.filt.R2.errs <- learnErrors(cpath("sense","filtered","R2"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)
set.seed(42)
antisense.filt.R1.errs <- learnErrors(cpath("antisense","filtered","R1"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)
set.seed(42)
antisense.filt.R2.errs <- learnErrors(cpath("antisense","filtered","R2"), multithread=TRUE, randomize=TRUE, nbases=1e+08, verbose=TRUE)


############## DENOISE ##############
############## DENOISE ##############

# report
writeLines("\n...\ndada2 denoising\n")
Sys.sleep(3)

# run dada denoising - takes time with pool=TRUE
sense.filt.R1.dada <- dada(cpath("sense","filtered","R1"), err=sense.filt.R1.errs, multithread=TRUE, pool=TRUE)
sense.filt.R2.dada <- dada(cpath("sense","filtered","R2"), err=sense.filt.R2.errs, multithread=TRUE, pool=TRUE)
antisense.filt.R1.dada <- dada(cpath("antisense","filtered","R1"), err=antisense.filt.R1.errs, multithread=TRUE, pool=TRUE)
antisense.filt.R2.dada <- dada(cpath("antisense","filtered","R2"), err=antisense.filt.R2.errs, multithread=TRUE, pool=TRUE)


############## DEREPLICATE ##############
############## DEREPLICATE ##############

# report
writeLines("\n...\nDereplication\n")
Sys.sleep(3)

# derep
sense.filt.R1.derep <- derepFastq(cpath("sense","filtered","R1"))
sense.filt.R2.derep <- derepFastq(cpath("sense","filtered","R2"))
antisense.filt.R1.derep <- derepFastq(cpath("antisense","filtered","R1"))
antisense.filt.R2.derep <- derepFastq(cpath("antisense","filtered","R2"))


############## MERGE ##############
############## MERGE ##############

# report
writeLines("\n...\nRead merging\n")
Sys.sleep(3)

# merge the R1 and R2
sense.merged <- mergePairs(dadaF=sense.filt.R1.dada, derepF=sense.filt.R1.derep, dadaR=sense.filt.R2.dada, derepR=sense.filt.R2.derep, verbose=TRUE, maxMismatch=0)
antisense.merged <- mergePairs(dadaF=antisense.filt.R1.dada, derepF=antisense.filt.R1.derep, dadaR=antisense.filt.R2.dada,  derepR=antisense.filt.R2.derep, verbose=TRUE, maxMismatch=0)

# make an OTU table
sense.seqtab <- makeSequenceTable(sense.merged)
antisense.seqtab <- makeSequenceTable(antisense.merged)

# reverse comp the antisense
colnames(antisense.seqtab) <- dada2::rc(colnames(antisense.seqtab))

# fix the names before merging
rownames(sense.seqtab) <- str_split_fixed(rownames(sense.seqtab),"\\.",4)[,1]
rownames(antisense.seqtab) <- str_split_fixed(rownames(antisense.seqtab),"\\.",4)[,1]

# merge the tables
merged.seqtab <- mergeSequenceTables(table1=sense.seqtab, table2=antisense.seqtab, repeats="sum")


############## REMOVE CHIMAERAS ##############
############## REMOVE CHIMAERAS ##############

# report
writeLines("\n...\nDetecting chimaeras\n")
Sys.sleep(3)

# remove chimaeras
merged.seqtab.nochim <- removeBimeraDenovo(merged.seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim <- as.data.frame(merged.seqtab.nochim)
seqtab.nochim.t <- as.data.frame(t(seqtab.nochim))

############## SAVE FILES ##############
############## SAVE FILES ##############

# report
writeLines("\n...\nSaving raw ASV files\n")
Sys.sleep(3)

# make df and fasta for IDs
otus.df <- tibble(names=paste0("MiFish_UE-","asv",str_pad(seq_along(colnames(merged.seqtab.nochim)),width=4,side="left",pad="0")), dnas=colnames(merged.seqtab.nochim)) %>% mutate(len=str_length(dnas))

# write out
write.FASTA(tab2fas(df=otus.df, seqcol="dnas", namecol="names"), file=paste0(proj.path,"/results/asvs.fa"))
write_tsv(as_tibble(otus.df, rownames="asv"), file=paste0(proj.path,"/results/asvs.tsv"))

# save the OTU table with sample names as df
colnames(seqtab.nochim) <- paste0("MiFish_UE-","asv",str_pad(seq_along(colnames(seqtab.nochim)),width=4,side="left",pad="0"))
write_tsv(as_tibble(t(seqtab.nochim), rownames="asv"), file=paste0(proj.path,"/results/asv-table.tsv"))

#count number of sequences per length category
library(Biostrings)
s <- readDNAStringSet(paste0(proj.path,"/results/asvs.fa"))
sum(width(s))/length(s)
sd(width(s))
sum(width(s)>186)
sum(width(s)<185)
sum(width(s)==185)
sum(width(s)<163)
min(width(s))
sort(width(s))
max(width(s))

# save runinfo as df
sum(merged.seqtab.nochim)/sum(merged.seqtab)
getN <- function(x) sum(getUniques(x))
trackdenoise.sense <- cbind(sapply(sense.filt.R1.dada, getN), sapply(sense.filt.R2.dada, getN))
trackdenoise.antisense <- cbind(sapply(antisense.filt.R1.dada, getN),sapply(antisense.filt.R2.dada, getN))
trackmerge.sense <- cbind(sapply(sense.merged, getN))
trackmerge.antisense <- cbind(sapply(antisense.merged, getN))
trackmerge <- cbind(rowSums(merged.seqtab), rowSums(merged.seqtab.nochim))

dada2_track <- merge(out.sense, trackdenoise.sense, by=0, all.x=TRUE)
colnames(dada2_track) <- c("samples","input sense", "filter sense", "denoised F", "denoised R")
dada2_track <- merge(dada2_track, trackmerge.sense, by.x=1, by.y=0, all.x=TRUE)
colnames(dada2_track) <- c("samples","input sense", "filter sense","denoised F", "denoised R",
                           "merge F+R")
dada2_track <- merge(dada2_track, out.antisense, by.x=1, by=0, all.x=TRUE)
colnames(dada2_track) <- c("samples","input sense", "filter sense", "denoised F", "denoised R",
                           "merge F+R",
                           "input antisense","filter antisense")
dada2_track <- merge(dada2_track, trackdenoise.antisense, by.x=1, by.y=0, all.x=TRUE)
colnames(dada2_track) <- c("samples","input sense", "filter sense", "denoised F", "denoised R",
                           "merge F+R",
                           "input antisense","filter antisense","denoised F AS", "denoised R AS")
dada2_track <- merge(dada2_track, trackmerge.antisense, by.x=1, by.y=0, all.x=TRUE)
colnames(dada2_track) <- c("samples","input sense", "filter sense","denoised F", "denoised R",
                           "merge F+R",
                           "input antisense","filter antisense", "denoised F AS", "denoised R AS",
                           "merge F+R AS")
dada2_track$samples <- gsub(".R1.fastq.gz", "", dada2_track$samples)
dada2_track <- merge(dada2_track, trackmerge, by.x=1, by.y=0, all.x=TRUE)
colnames(dada2_track) <- c("samples","input sense", "filter sense", "denoised F", "denoised R",
                           "merge F+R",
                           "input antisense","filter antisense", "denoised F AS", "denoised R AS",
                           "merge F+R AS", 
                           "remove chimeras", "#ASV")
dada2_track[is.na(dada2_track)] <- 0

write_tsv(as_tibble(dada2_track, rownames="names"), file=paste0(proj.path,"/results/dada2_track.tsv"))

# report
writeLines(paste0("\n...\nASVs written to:\n",paste0(proj.path,"/results/asvs.fna"),"\n"))
writeLines(paste0("\n...\nASV table written to:\n",paste0(proj.path,"/results/asv-table.tsv"),"\n"))

#save environment after merging in DADA2
save.image(paste0(proj.path,"/results/R_Environment/REnvironment_DADA2_merging.RData"))

#remove data that will not be used for further analysis
rm(antisense.filt.R1.dada, antisense.filt.R1.derep, antisense.filt.R1.errs, 
   antisense.filt.R2.dada, antisense.filt.R2.derep, antisense.filt.R2.errs, antisense.merged,
   sense.filt.R1.dada, sense.filt.R1.derep, sense.filt.R1.errs,
   sense.filt.R2.dada, sense.filt.R2.derep, sense.filt.R2.errs,
   out.antisense, out.sense, antisense.seqtab, sense.seqtab,
   trackfilter, trackdenoise.sense, trackdenoise.antisense, trackmerge.sense, trackmerge.antisense, trackmerge)

#count number of sequences per sample an sort from lowest to highest
NumberOfSequences <-as.data.frame(sort.default(rowSums(seqtab.nochim[1:nrow(seqtab.nochim),]))) #the total number of reads per sample
NumberOfRawReads <- as.data.frame(colSums(seqtab.nochim)) # the total number of reads per ASV
NumberOfASVs <- as.data.frame(rowSums(seqtab.nochim !=0)) # the total number of ASVs per sample
write_tsv(as_tibble(NumberOfASVs, rownames="asv"), file=paste0(proj.path,"/results/NumberOfASVs.tsv"))

# remove samples with low readdepth
noreads <- c(rownames(seqtab.nochim[rowSums(seqtab.nochim)==0,]))
lowreads <- c(rownames(seqtab.nochim[rowSums(seqtab.nochim)<300,]))
seqtab.nochim2<-seqtab.nochim[!rowSums(seqtab.nochim)<300,] #needed to add sample to the rarecurve

#Create rarecurve before concatenation of PCR replicates
rarecurve(seqtab.nochim, ylab = "ASVs", main = "Rarecurve of unconcatenated samples before taxonomic assignment", label = FALSE, step =1000) 
rarecurve(seqtab.nochim2, ylab = "ASVs", main = "Rarecurve of unconcatenated samples before taxonomic assignment", sample = 10000, label = FALSE, step =1000) #is raremax situated just after the steep part of the graph?
NotRarefied <- c(rownames(seqtab.nochim[rowSums(seqtab.nochim)<10000,]))

#concatenate PCR replicates from the unrarefied samples (samples with low read depth are already removed)
seqtab.concat <- as.data.frame(seqtab.nochim)
seqtab.concat$names <- str_sub(rownames(seqtab.concat), end=-3) # removes PCR number (_S1, _S2 or _S3) from samplenames
seqtab.concatenated<-aggregate(seqtab.concat[,1:ncol(seqtab.concat)-1], by= list(seqtab.concat$names),FUN=sum)
rownames(seqtab.concatenated) <- seqtab.concatenated$Group.1
seqtab.concatenated<-seqtab.concatenated[,2:ncol(seqtab.concat)]  
NumberOfSequencesConcat <-as.data.frame(sort.default(rowSums(seqtab.concatenated[1:nrow(seqtab.concatenated),])))
NumberOfASVsConcat <- as.data.frame(rowSums(seqtab.concatenated !=0))
write_tsv(as_tibble(NumberOfASVsConcat, rownames="asv"), file=paste0(proj.path,"/results/NumberOfASVsConcat.tsv"))

# remove samples without reads
noreadsConcat <- c(rownames(seqtab.concatenated[rowSums(seqtab.concatenated)==0,]))
lowreadsConcat <- c(rownames(seqtab.concatenated[rowSums(seqtab.concatenated)<300,]))
seqtab.concatenated2<-seqtab.concatenated[!rowSums(seqtab.concatenated)<300,] #needed to add sample to the rarecurve

#Create rarecurve after concatenation of PCR replicates and rarefy data
raremax <- 25000
rarecurve(seqtab.concatenated2, ylab = "ASVs", main = "Rarecurve of concatenated samples before taxonomic assignment", sample = raremax, label = FALSE, step =1000) #is raremax situated just after the steep part of the graph?
NotRarefiedConcat <- c(rownames(seqtab.concatenated[rowSums(seqtab.concatenated)<raremax,]))
seqtab.concatenated.rarefied <- rrarefy(seqtab.concatenated[1:nrow(seqtab.concatenated),],raremax)
#seqtab.nochim.rarefied <- rrarefy(seqtab.nochim[1:nrow(seqtab.nochim),],raremax)
seqtab.concatenated.rarefied2 <- seqtab.concatenated.rarefied[!rowSums(seqtab.concatenated.rarefied)<raremax,]
RemoveASVs <- c(colnames(seqtab.concatenated.rarefied[,colSums(seqtab.concatenated.rarefied)==0]))
seqtab.concatenated.rarefied3 <- seqtab.concatenated.rarefied[,!colSums(seqtab.concatenated.rarefied2)==0]
seqtab.concatenated.rarefied4 <- seqtab.concatenated.rarefied2[,!colSums(seqtab.concatenated.rarefied2)==0]
rarecurve(seqtab.concatenated.rarefied, ylab = "ASVs", main = "Rarecurve of rarefied samples before taxonomic assignment", label = FALSE, step =1000)
saveRDS(seqtab.concatenated.rarefied, "seqtab.concatenated.rarefied.rds")
saveRDS(seqtab.concatenated, "seqtab.concatenated.unrarefied.rds")
write_tsv(as_tibble(seqtab.concatenated.rarefied, rownames="asv"), file=paste0(proj.path,"/results/seqtab.concatenated.rarefied.tsv"))
write_tsv(as_tibble(seqtab.concatenated.rarefied3, rownames="asv"), file=paste0(proj.path,"/results/seqtab.concatenated.rarefied.ASVsRemoved.tsv"))
colnames(seqtab.concatenated.rarefied) <- otus.df$dnas

# add taxa
set.seed(100)
#tax_rare <- assignTaxonomy(seqtab.concatenated.rarefied3, "/home/genomics/icornelis/02_ZEROimpact/02_COI/NJ2021/Step3_DADA2/GEANS_ref_v4.fas", minBoot = 80, multithread = TRUE)
tax_rare <- assignTaxonomy(seqtab.concatenated.rarefied3, "/home/genomics/icornelis/03_RawScripts/Step4_COI_ReferenceDB_TaxonomicAssignment/WIDE_MIDORI_dada2_no_insects.fasta", minBoot = 80, multithread = TRUE)
#tax_rare2 <- addSpecies(tax_rare, "/home/genomics/MiFish-UE_run2/Step3_DADA2/icornelis/12S references_species-assignment.fa")
taxa <- assignTaxonomy(merged.seqtab.nochim, "/home/genomics/icornelis/03_RawScripts/Step4_COI_ReferenceDB_TaxonomicAssignment/WIDE_MIDORI_dada2_no_insects.fasta", minBoot = 80, multithread = TRUE)
write_tsv(as_tibble(tax_rare, rownames="ASV"), file=paste0(proj.path,"/results/taxonomic-assignment_rarefied.tsv"))
#write_tsv(as_tibble(tax_rare2, rownames="ASV"), file=paste0(proj.path,"/results/taxonomic-assignment_rarefied_2.tsv"))
write_tsv(as_tibble(taxa, rownames="ASV"), file=paste0(proj.path,"/results/taxonomic-assignment_unrarefied.tsv"))
taxa_rarefied <- data.frame(tax_rare)
#taxa_rarefied2 <- data.frame(tax_rare2)
taxa_unrarefied <- data.frame(taxa)
#write_tsv(as_tibble(taxa_rarefied), file=paste0(proj.path,"/results/TaxonomicAssignment_Macrobenthos.tsv"))
write_tsv(as_tibble(taxa_rarefied), file=paste0(proj.path,"/results/TaxonomicAssignment_MIDORI.tsv"))

#for unrarefied data
table_unrarefied<-as.data.frame(t(merged.seqtab.nochim))
#table_rarefied <- as.data.frame(t(merged.seqtab.nochim.rarefied))
table_unrarefied.concatenated <- as.data.frame(t(seqtab.concatenated))
table_rarefied.concatenated <- as.data.frame(t(seqtab.concatenated.rarefied3))

#remove rownames
rownames(table_unrarefied)<-NULL
#rownames(table_rarefied)<-NULL
rownames(table_unrarefied.concatenated)<-NULL
rownames(table_rarefied.concatenated)<-NULL

#add extra columns to table to add taxonomy, for this first transform the tax object to data frame
table_unrarefied$Kingdom<-taxa_unrarefied$Kingdom
table_unrarefied$Phylum<-taxa_unrarefied$Phylum
table_unrarefied$Class<-taxa_unrarefied$Class
table_unrarefied$Order<-taxa_unrarefied$Order
table_unrarefied$Family<-taxa_unrarefied$Family
table_unrarefied$Genus<-taxa_unrarefied$Genus
table_unrarefied$Species<-taxa_unrarefied$Species
table_unrarefied$ASV<-otus.df$names
write_tsv(as_tibble(table_unrarefied), file=paste0(proj.path,"/results/table_unrarefied.tsv"))

#table_rarefied$Kingdom<-taxa_rarefied$Kingdom
#table_rarefied$Phylum<-taxa_rarefied$Phylum
#table_rarefied$Class<-taxa_rarefied$Class
#table_rarefied$Order<-taxa_rarefied$Order
#table_rarefied$Family<-taxa_rarefied$Family
#table_rarefied$Genus<-taxa_rarefied$Genus
#table_rarefied$Species<-taxa_rarefied$Species
#table_rarefied$ASV<-otus.df$names
#write_tsv(as_tibble(table_rarefied), file=paste0(proj.path,"/results/table_rarefied.tsv"))

table_unrarefied.concatenated$Kingdom<-taxa_unrarefied$Kingdom
table_unrarefied.concatenated$Phylum<-taxa_unrarefied$Phylum
table_unrarefied.concatenated$Class<-taxa_unrarefied$Class
table_unrarefied.concatenated$Order<-taxa_unrarefied$Order
table_unrarefied.concatenated$Family<-taxa_unrarefied$Family
table_unrarefied.concatenated$Genus<-taxa_unrarefied$Genus
table_unrarefied.concatenated$Species<-taxa_unrarefied$Species
table_unrarefied.concatenated$ASV<-otus.df$names
write_tsv(as_tibble(table_unrarefied.concatenated), file=paste0(proj.path,"/results/table_unrarefied_concatenated.tsv"))

SelectASVs <- c(rownames(taxa_rarefied))
otus.df.rarefied <- otus.df[otus.df$dnas %in% SelectASVs,]

table_rarefied.concatenated$Kingdom<-taxa_rarefied$Kingdom
table_rarefied.concatenated$Phylum<-taxa_rarefied$Phylum
table_rarefied.concatenated$Class<-taxa_rarefied$Class
table_rarefied.concatenated$Order<-taxa_rarefied$Order
table_rarefied.concatenated$Family<-taxa_rarefied$Family
table_rarefied.concatenated$Genus<-taxa_rarefied$Genus
table_rarefied.concatenated$Species<-taxa_rarefied$Species
table_rarefied.concatenated$ASV<-otus.df.rarefied$names
write_tsv(as_tibble(table_rarefied.concatenated), file=paste0(proj.path,"/results/table_rarefied_concatenated.tsv"))

#get a list of the different families and species detected
Families_ur<-as.data.frame(unique(table_unrarefied$Family))
write_tsv(as_tibble(Families_ur), file=paste0(proj.path,"/results/Families_unrarefied.tsv"))
Families_r_co<-as.data.frame(unique(table_rarefied.concatenated$Family))
write_tsv(as_tibble(Families_r_co), file=paste0(proj.path,"/results/Families_rarefied.tsv"))
Genera<-as.data.frame(unique(table_rarefied.concatenated$Genus))
write_tsv(as_tibble(Genera), file=paste0(proj.path,"/results/Genera_rarefied.tsv"))
Species_ur<-as.data.frame(unique(table_unrarefied$Species))
write_tsv(as_tibble(Species_ur), file=paste0(proj.path,"/results/Species_unrarefied.tsv"))
Species_r_co<-as.data.frame(unique(table_rarefied.concatenated$Species))
write_tsv(as_tibble(Species_r_co), file=paste0(proj.path,"/results/Species_rarefied.tsv"))

#calculate number of samples, store sample names in list
samplenumber<-ncol(table_unrarefied)-8
samples<-colnames(table_unrarefied)[1:samplenumber]
samplenumber_rare<-ncol(table_rarefied.concatenated)-8
samples_rare<-colnames(table_rarefied.concatenated)[1:samplenumber_rare]

#add new colums with percentages
freq_table <- data.frame(matrix(nrow=nrow(table_unrarefied), ncol=samplenumber
                                , dimnames = list(c(rownames(taxa)), c(colnames(table_unrarefied[1:samplenumber])))))
freq_table_rare <- data.frame(matrix(nrow=nrow(table_rarefied.concatenated), ncol=samplenumber_rare
                                    , dimnames = list(c(rownames(tax_rare)), c(colnames(table_rarefied.concatenated[1:samplenumber_rare])))))
#calculate frequencies
for(i in 1:samplenumber){freq_table[,i] <- table_unrarefied[,i]/colSums(table_unrarefied[1:samplenumber])[i]}
for(i in 1:samplenumber_rare){freq_table_rare[,i] <- table_rarefied.concatenated[,i]/colSums(table_rarefied.concatenated[1:samplenumber_rare])[i]}

#add species information to frequency table
freq_table$Kingdom<-taxa[,1]
freq_table$Phylum<-taxa[,2]
freq_table$Class<-taxa[,3]
freq_table$Order<-taxa[,4]
freq_table$Family<-taxa[,5]
freq_table$Genus<-taxa[,6]
freq_table$Species<-taxa[,7]
freq_table$ASV<-otus.df$names
write_tsv(as_tibble(freq_table), file=paste0(proj.path,"/results/frequency_table.tsv"))

#add species information to frequency table
freq_table_rare$Kingdom<-tax_rare[,1]
freq_table_rare$Phylum<-tax_rare[,2]
freq_table_rare$Class<-tax_rare[,3]
freq_table_rare$Order<-tax_rare[,4]
freq_table_rare$Family<-tax_rare[,5]
freq_table_rare$Genus<-tax_rare[,6]
freq_table_rare$Species<-tax_rare[,7]
freq_table_rare$ASV<-otus.df.rarefied$names
write_tsv(as_tibble(freq_table_rare), file=paste0(proj.path,"/results/frequency_table_rare.tsv"))

##Check the diversity after rarefication for all ASVs, assigned ASVs and ASVs assigned to Fish
#remove negative controls
env <- read.csv("/home/genomics/icornelis/02_ZEROimpact/02_COI/NJ2021/Step3_DADA2/environmental_data.csv",  header=TRUE, sep=";")
#env$Environment_color <- ifelse(env$Environment=="Coastal","limegreen", 
#                                ifelse(env$Environment=="outside_OWF", "slateblue", 
#                                       ifelse(env$Environment=="inside_OWF","darkorange","red")))
env$Environment_color <- ifelse(env$Zones=="Coastal","limegreen", 
                                ifelse(env$Zones=="zone1", "slateblue", 
                                       ifelse(env$Zones=="zone2","darkorange","red")))

#remove samples with <raremax reads in total
table_rarefied.concatenated2 <- table_rarefied.concatenated[,!colSums(table_rarefied.concatenated[,1:(ncol(table_rarefied.concatenated)-8)])<raremax]
table_rarefied.concatenated2$Species <- table_rarefied.concatenated$Species

#remove unassigned ASVs
table_rarefied_assigned <- table_rarefied.concatenated[!table_rarefied.concatenated$Species %in% c(NA),]
merged.seqtab.nochim.assigned <- as.data.frame(t(table_rarefied_assigned[,1:(ncol(table_rarefied_assigned)-8)]))
table_rarefied_assigned2 <- table_rarefied.concatenated2[!table_rarefied.concatenated2$Species %in% c(NA),]
merged.seqtab.nochim.assigned2 <- as.data.frame(t(table_rarefied_assigned2[,1:(ncol(table_rarefied_assigned2)-8)]))

#remove non-Fish
#table_rarefied_fish <- table_rarefied.concatenated2[!table_rarefied.concatenated2$Order %in% c(NA,"Thalassiosirales"),]
#write_tsv(as_tibble(table_rarefied_fish), file=paste0(proj.path,"/results/seq_fish_order_rarefied.tsv"))
#merged.seqtab.nochim.fish <- as.data.frame(t(table_rarefied_fish[,1:(ncol(table_rarefied_fish)-8)]))
#merged.seqtab.nochim.fish <- merged.seqtab.nochim.fish %>% filter(!grepl("neg", rownames(merged.seqtab.nochim.fish)))
#colnames(merged.seqtab.nochim.fish) <- table_rarefied_fish$ASV

#remove samples from environmental data that are now absent from the data
keep_samples <- c(rownames(seqtab.concatenated2))
env<-env[env$Demultiplexing %in% keep_samples,]
keep_samples_rarefired <- c(rownames(seqtab.concatenated.rarefied3))
env_rarefied <- env[env$Demultiplexing %in% keep_samples_rarefired,]
#keep_samples_fish <- c(rownames(merged.seqtab.nochim.fish))
#env_fish <- env_rarefied[env_rarefied$Niskin.sample %in% keep_samples_fish,]

#Put samples in same order as environmental data
seqtab.concatenated2 <- seqtab.concatenated2[env$Demultiplexing,]
seqtab.concatenated.rarefied <- seqtab.concatenated.rarefied[env_rarefied$Demultiplexing,]
seqtab.concatenated.rarefied3 <- seqtab.concatenated.rarefied3[env_rarefied$Demultiplexing,]
merged.seqtab.nochim.assigned <- merged.seqtab.nochim.assigned[env_rarefied$Demultiplexing,]
merged.seqtab.nochim.assigned2 <- merged.seqtab.nochim.assigned2[env_rarefied$Demultiplexing,]
merged.seqtab.nochim.assigned2 <- na.omit(merged.seqtab.nochim.assigned2)
keep_samples_assigned <- c(rownames(merged.seqtab.nochim.assigned2))
env_assigned <- env[env$Demultiplexing %in% keep_samples_assigned,]

#rarecurve with color to check the biodiversity after rarefication
rarecurve(seqtab.concatenated.rarefied, ylab = "ASVs", main = "Rarecurve of rarefied samples before taxonomic assignment", label=FALSE, col = as.vector(env_rarefied$Environment_color), step = 1000)
rarecurve(seqtab.concatenated.rarefied3, ylab = "ASVs", main = "Rarecurve of rarefied samples before taxonomic assignment", label=FALSE, col = as.vector(env_rarefied$Environment_color), step = 1000)
rarecurve(merged.seqtab.nochim.assigned, ylab = "Assigned ASVs", main = "Rarecurve of rarefied samples after taxonomic assignment", col = as.vector(env_rarefied$Environment_color), label = FALSE, step =1000)
rarecurve(merged.seqtab.nochim.assigned2, ylab = "Assigned ASVs", main = "Rarecurve of rarefied samples after taxonomic assignment", col = as.vector(env_assigned$Environment_color), label = FALSE, step =1000)
#rarecurve(seqtab.concatenated2, ylab = "ASVs", main = "Rarecurve of concatenated samples before taxonomic assignment", col = as.vector(env$Environment_color), sample = raremax, label = FALSE, step =1000)
#rarecurve(merged.seqtab.nochim.fish, ylab = "ASVs assigned to fish", main = "Rarecurve of rarefied samples after taxonomix assignment", label=FALSE, col = as.vector(env_fish$Environment_color), step = 1000)  
#barplot(rowSums(merged.seqtab.nochim), main="Reads per sample", ylab="Number of reads", xlab="", las=2, col = as.vector(env$Environment_color))
#barplot(rowSums(merged.seqtab.nochim2>0), main="ASVs per sample", ylab="Number of reads", xlab="", las=2, col = as.vector(env$Environment_color))
#barplot(rowSums(merged.seqtab.nochim.samples>0), main="ASVs per sample", ylab="Number of reads", xlab="", las=2, col = as.vector(env.samples$Environment_color))
#barplot(rowSums(merged.seqtab.nochim.fish>0), main="ASVs per sample", ylab="Number of reads", xlab="", las=2, col = as.vector(env.samples$Environment_color))

#Create rarecurve without lines at raremax
out <- rarecurve(seqtab.concatenated2, ylab = "ASVs", main = "Rarecurve of concatenated samples before taxonomic assignment", col = as.vector(env$Environment_color), sample = raremax, label = FALSE, step =1000)
pars <- expand.grid(col = as.vector(env$Environment_color), lty = "solid", lwd = 1, 
                    stringsAsFactors = FALSE)
Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
Smax <- sapply(out, max)
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Sample Size",
     ylab = "ASVs", type = "n", main = "Rarecurve of concatenated samples before taxonomic assignment")
abline(v = raremax)
for (i in seq_along(out)) {
  N <- attr(out[[i]], "Subsample")
  with(pars, lines(N, out[[i]], col = col[i], lty = lty[i], lwd = lwd[i]))
}
