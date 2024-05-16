## eDNA metabarcoding Bioinformatics Offshore wind farm monitoring with eDNA
### Description
These scripts were used for the data anlysis of the seawater samples collected in autumn 2021. 
Amplification was performed with a one-step PCR-protocol using uniquely tagged primersets for each sample. Sequencing was done on a MiSeq illumina platform.
Detailed information about the scripts and how to use them can be found in the document **eDNAmetabarcoding_DataAnalysis_Workflow_MiSeq.docx**.\
The first steps in the processing include the quality control using FastQC,and the preparation of the datasets for demultiplexing.\
All samples were used for 12S and COI eDNA metabarcoding, all scripts used specifically for either of the two barcodes are stored seperatly in the 12S and COI folders, respectively. 

### Step2: Demultiplexing and concatenation of the PCR replicates per sample
In Step2a the raw reads from the paired R1 and R2 files will be demultiplexed and trimmed based on the sample barcodes. This generates separate sequencing files for each sample. In Step2b the demultiplexed files from the three technical replicates will be concatenated into one folder. To keep the files from the technical replicates separated within the same folder, the file name is extened with _S1, _S2 and _S3 for PCR1, PCR2 and PCR3, respectively.
#### Step2a – Demultiplexing
The raw reads will be demultiplexed based on the unique primer-tags, using Cutadapt. After demultiplexing the raw reads will be trimmed by removing the primer- and tag-sequences from the 5’-end and 3’-end of the read. In addition, trimmed sequences that are shorter than a set length will be removed from the file.
After running the script the demultiplexed and trimmed reads will be sorted into 4 folder: 2 folders in the sense folder, trimmed-R1 and trimmed-R2, and 2 folders in the antisense folder, trimmed-R1 and trimmed-R2. Each folder, contains one file per sample, that only contains the trimmed reads from that sample.
The script consits of four parts:\
1. The raw sequences in the paired R1 and R2 file are **orientated** into sense and antisense sequences, and are sorted into the respective **sense or antisense** folder.\
To increase the number of sequences that are kept for further processing, the number of mismatches allowed was set to **-e 0.15**. For the MiFish_UE primers (12S), this corresponds to a maximum of three mismatches in the forward primer (21 bp) and a maximum of four mismatches in the reverse primer (27 bp). For the Leray primers (COI) this corresponds to a maximum of 3 mismatches in the forward and reverse primer (24 bp).
2. Is a **visualization step**, to check for the presence of the forward or reverse primer in the first 20 sequences in the R1 and R2 files from the sense folder. The primer sequence is colored for visualization.
3. Is the **demultiplexing step**, the sequences will only be sorted if the tag is found in both paired sequences through the command --pair adapter, which means that for one sequence the tag is found in both the R1 and R2 file. During demultiplexing **no mismatches** are allowed, **-e 0.00**.
4. Is the **trimming step** to remove the primer sequences and their tag. This is first done at the 3’-end of the sequences using the reverse complement of the two primer sequences and then at the 5’-end of the sequences. In addition, all sequences with less nucleotides than the set TRUCVAL-value will be removed.
#### Step2b - Concatenate
After demultiplexing the sequences of the three technical replicates will be concatenated into one folder. This folder contains a processed-reads folder with the same structure as the folders for the separate technical replicates.

### Step3 - DADA2.R
Through seven different steps, the DADA2 pipeline will use the demultiplexed fastq files to produce amplicon sequence variants (ASVs). The output is a matrix containing the samples (rows) and ASVs (columns), in which the value of each entry is the number of times that ASV was observed in that sample (number of reads). After receiving the ASV table, the script will use the DADA2 function **assignTaxonomy** for the taxonomic assignment of the ASVs bases on the custom reference databases. Both references can be found in the folder Reference Databases (12S_references_20230705.fas and COI_references_WIDE_MIDORI_dada2_no_insects.fasta). the DADA2 function **assignTaxonomy** uses the naive Bayesian classifier method and makes taxonomic assignments by comparing the ASV sequences to a custom reference database. The function assignTaxonomy will only return ASVs at species level when more than 80 of 100 bootstraps result in the same species level. This is also true for the higher taxonomic levels.
More information about the different steps can be found in the Introduction to dada2 (https://benjjneb.github.io/dada2/tutorial.html).

### Step4 - nucleotide BLAST (BLASTn)
The BLASTn is used for the taxonomic assignment of the ASVs that remained unassigned after running DADA2. The cause may be the absence of the species from the custom reference database or a low taxonomic resolution between closely related species. BLASTn is run twice, once against the custom reference database, also used with assignTaxonomy in DADA2, and once against GenBank. Sequences were assigned to a species if there was ≥ 97% sequence identity and a query coverage of at least 75%. The chosen query coverage threshold was set to be less stringent in order to explore ASVs belonging to other kingdoms than Animalia for both marker genes.\ 
BLAST uses a heuristic method, which first starts to search for short matches, seeds, between two sequences. During the second step, local alignments are made between the matches found during the first search. However, this means that optimal alignments are not guaranteed and that some hits may be missed. Although the results from the nucleotide BLAST will only be used for the taxonomic assignment of the unassigned ASVs, blastn will be run on the asvs.fa file that contains all the ASVs. This results in an extra check for the results from assignTaxonomy in DADA2.

### Step5 - Get the Full Taxonomic Assignment
The script will be used to combine the taxonomic information gathered from DADA2 and BLASTn against the custom reference database and GenBank (tophits only).
During the first step the results from the taxonomic assignment in DADA2 are then copied into a new column Full. During the second step, the results from BLASTn against the custom reference database are added to all the ASVs that remained unassigned by the function assignTaxonomy in DADA2. During the final step, the results from BLASTn against GenBank will be used to add a taxonomic assignment to the ASVs that still remain unassigned after the second step.
Before getting the full taxonomic assignment all scientific names are check with WoRMS, and unaccepted scientific names are replaced by accepted scientific names. After getting the full taxonomic assignment the Classification received by the taxonomic assignment in DADA2 is replaced by the Classification from WoRMS and NCBI (for species that are not found in WoRMS).
For some ASVs taxonomic assignment to species level is not possible (due to it being absent from the reference database or due to identical barcode sequences). For these ASVs, the script uses a higher classification level (until Familiy-level) and adds sp. to it. If the Family is unknown the taxonomic assignment remains NA.
After running the Step5_Script_SequenceTable_FullTaxonomicAssignment.R the final ASV-table with the taxonomic assignment for all ASVs is received. 

### Step6 - Decontamination
Decontamination of the contaminant ASVs was performed on the raw ASV-table with the full taxonomic assignment (last 11 columns) using microDecon. The microDecon package will remove contaminant reads, rather than contaminant ASVs.\
After the decontamination step, the technical replicates (PCR replicates) of each sample are concatenated by summation. 
After running the script three files are provided, the decontaminated ASV-table with the technical replicates still separated, the decontaminated ASV-table with the concatenated technical replicates and a RDS-file containing the metadata of each sample, which will be used for data analysis.
 
### Step7: Continue with data analysis
- Step7a -  Create a stacked barplot with the origin of the ASVs and reads\
The input files are the ASV-table with full taxonomic assignment before decontamination and a file containing the metadata of each sample.\
The script will provide a stack barplot with the ASV-count and read-count in each sample. The colors of the bars represent the species origin of the ASVs/reads.
- Step7b - Create a Heatmap\
The input files are the decontaminated ASV-table with full taxonomic assignment and a file containing the metadata of each sample.\
The script will provide a heatmap representing the relative read abundance of each species per sample. The fish species are clustered using the ward.D method in the hclust function in the stats package. Before clustering the data was double transformed by converting the read counts to the relative abundance in each sample (using total in the decostand function) followed by scaling the relative abundances of each ASV to the highest observed relative abundance across all samples (using max in the decostand function)
- Step7c - Alpha diversity analysis\
The input files are the decontaminated ASV-table with full taxonomic assignment, the data from the morphological identification and the files containing the metadata of each sample.\
The Alpha diversity analysis is performed on the un-rarefied and rarefied data. The coverage based rarefaction is performed within the same script.\
The script will provide the plots, a box-plot and line-plot, with the alpha diversity analysis. In addition the statistical analysis using a generalized linear model is performed. 
- Step7d - Permanova analysis\
The input files are the decontaminated ASV-table with full taxonomic assignment, or the data from the morphological identification, and the file containing the metadata of each sample.\
The permanova analysis is performed in the un-rarefied data after double transformation. 
- Step7e - Species community structures (nMDS)\
The input files are the decontaminated ASV-table with full taxonomic assignment, the data from the morphological identification, and the file containing the metadata of each sample.\
The beta-diversity is analysed based on the un-rarefied data after double transformation.



