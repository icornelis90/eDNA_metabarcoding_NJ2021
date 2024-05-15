#title: "Full taxonomic assignment"
#author: "Joran Vanhollebeke & Isolde Cornelis"
#date: "27/02/2024"

################################################################################

#activate libraries
libraries <- c("here"
               , "taxize"
               , "curl"
               , "readr"
               , "xlsx"
               , "stringr"
               , "tibble"
               , "vegan"
               , "here"
               , "dplyr"
               , "tidyr"
)
for (opties in libraries){
  
  if (opties %in% installed.packages()){
    
    library(opties,character.only = TRUE)
    
  } else {install.packages(opties,repos = "http://cran.us.r-project.org")
    
    library(opties,character.only = TRUE)
  }
}

# make paths
proj.path <- here("/home/genomics/icornelis/02_ZEROimpact/01_12S/NJ2021/MiFish-UE_run2")

#Define blast genbank
{
  blastn_gen_com<-paste(
    "/usr/local/bin/blastn" 
    ,"-num_threads 5" 
    ,"-db /home/genomics/bioinf_databases/genbank/nt2/nt" 
    ,"-query ", paste0(proj.path,"/MiFish_UE-S_concatenated/results/asvs.fa")
    ,"-out "  , paste0(proj.path,"/MiFish_UE-S_concatenated/blastn_GenBank_",format(Sys.time(), "%Y%m%d"),".tsv")
    ,"-max_hsps 1 -max_target_seqs 1 -qcov_hsp_perc 75 -perc_identity 90" 
    ,"-outfmt '6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen evalue bitscore qcovs'"
  )
  # excute blastn
  system(blastn_gen_com)
  
}

#Define blast reference db
{
  blastn_ref_com<-paste(
    "/usr/local/bin/blastn" 
    ,"-num_threads 5" 
    ,"-db /home/genomics/icornelis/03_RawScripts/Step4_12S_ReferenceDB_TaxonomicAssignment/12S_references.fa" 
    ,"-query ", paste0(proj.path,"/MiFish_UE-S_concatenated/results/asvs.fa")
    ,"-out "  , paste0(proj.path,"/MiFish_UE-S_concatenated/blastn_own_reference_tophits_",format(Sys.time(), "%Y%m%d"),".tsv")
    ,"-max_hsps 1 -max_target_seqs 1 -qcov_hsp_perc 75 -perc_identity 90" 
    ,"-outfmt '6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen evalue bitscore qcovs'"
  )
  # excute blastn
  system(blastn_ref_com)
}

{
  blastn_ref_com<-paste(
    "/usr/local/bin/blastn" 
    ,"-num_threads 5" 
    ,"-db /home/genomics/icornelis/03_RawScripts/Step4_12S_ReferenceDB_TaxonomicAssignment/12S_references.fa" 
    ,"-query ", paste0(proj.path,"/MiFish_UE-S_concatenated/results/asvs.fa")
    ,"-out "  , paste0(proj.path,"/MiFish_UE-S_concatenated/blastn_own_reference_top10hits_",format(Sys.time(), "%Y%m%d"),".tsv")
    ,"-max_hsps 1 -max_target_seqs 10 -qcov_hsp_perc 75 -perc_identity 90" 
    ,"-outfmt '6 qseqid sseqid staxids stitle pident qlen length mismatch gapopen evalue bitscore qcovs'"
  )
  # excute blastn
  system(blastn_ref_com)
}

# Read results
blastn_ref <- read_tsv(file.path(proj.path, "MiFish_UE-S_concatenated", "blastn_own_reference_tophits_20230705.tsv"), 
                       col_names = c("ASV", "sseqid", "staxids", "Taxonomy", "pident", "qlen", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs"))

# Separate tax in each level
blastn_ref <- separate(blastn_ref, Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")

# Read Genbank results
#blastn_Genbank <- read_tsv(file.path(proj.path, "MiFish_UE-S_concatenated", "blastn_GenBank_20230427.tsv"), 
#                           col_names = c("ASV", "sseqid", "staxids", "Species", "pident", "qlen", "length", "mismatch", "gapopen", "evalue", "bitscore", "qcovs"))




# upload the dataset (rarefied data with samples and taxonomic assignment), blastn results from GenBank and blastn results from own reference list
table_raw <- readxl::read_excel(paste0(proj.path,"/MiFish_UE-S_concatenated/results/table_unrarefied.xlsx"),
                                sheet = "table_unrarefied")
blastn_Genbank <- readxl::read_excel(paste0(proj.path,"/Step4_blastn/blastn_GenBank_20221024.xlsx"),
                                     sheet = "blastn_GenBank_20221024")
blastn_ref <- readxl::read_excel(paste0(proj.path,"/Step4_blastn/blastn_own_reference_tophits_20221024.xlsx"), 
                                 sheet="blastn_own_reference_tophits_20")

# in the blastn data select ASV and species column 
TA_DADA2 <- data.frame(table_raw$ASV, table_raw$Species) 
TA_DADA2_fish <- TA_DADA2[!table_raw$Species %in% c("Skeletonema pseudocostatum_(MK372941)", "NA"),]
blastn_Genbank_select <- blastn_Genbank[!blastn_Genbank$pident < 97000,]
blastn_Genbank <- data.frame(blastn_Genbank_select$ASV, blastn_Genbank_select$Species)
blastn_Genbank_select$sseqidlong <- blastn_Genbank_select$sseqid
blastn_Genbank_select$sseqid <- str_extract(blastn_Genbank_select$sseqidlong, "([A-Z])\\w+\\d")
blastn_Genbank_select$gi_id <- str_extract(blastn_Genbank_select$sseqidlong, "\\d{7,10}")

#blastn_ref <- data.frame(blastn_ref$ASV, blastn_ref$Species)
blastn_ref_select <- blastn_ref[!blastn_ref$pident < 97000,]
blastn_ref <- data.frame(blastn_ref_select$ASV, blastn_ref_select$Species)

#scientific name from blastn_Genbank
unique_sseqid <- unique(blastn_Genbank_select$sseqid)
ENTREZ_KEY <- "8402b2366108c13f3fa86e1a7d1519efa308" #replace with your personal API_KEY on ncbi, 
Classification_ncbi <- classification(genbank2uid(id = unique(blastn_Genbank_select$gi_id), 
                                                  key = ENTREZ_KEY), db = 'ncbi') # can be run without ENTREZ_KEY but it will run slower
Classification_NCBI_1 <-cbind(Classification_ncbi)
Classification_NCBI_1$gi_id <- unique(blastn_Genbank_select$gi_id)

# select just which lineage level to keep
Classification_NCBI <- select(Classification_NCBI_1, gi_id, query, superkingdom, 
                              phylum, class, order, family, genus, species)
colnames(Classification_NCBI) <- c("gi_id", "query", "Kingdom", "Phylum", "Class", 
                                   "Order", "Family", "Genus" , "Species")

idx <- which(blastn_Genbank_select$gi_id %in% Classification_NCBI$gi_id)
blastn_taxa_top_lineage <-
  blastn_Genbank_select |> left_join(Classification_NCBI,
                                            by ="gi_id")
blastn_Genbank$blastn_Genbank_select.Species <- blastn_taxa_top_lineage$Species.y

Classification_NCBI <- Classification_NCBI[!duplicated(Classification_NCBI$Species),]

#merge 3 datasets according to ASV name
table_Species <- merge(TA_DADA2, blastn_Genbank, by.x='table_raw.ASV', by.y='blastn_Genbank_select.ASV', all.x=TRUE)
table_Species_2 <- merge(table_Species, blastn_ref, by.x='table_raw.ASV', by.y='blastn_ref_select.ASV', all.x=TRUE)
colnames(table_Species_2) <- c("ASV", "DADA2", "blastn_GenBank", "blastn_ref")
#table_Species_3 <- table_Species_2[!duplicated(table_Species_3$table_Species_2.ASV),]

#write full taxonomic assignment in seperate column Full
table_Species_2$Full <- str_extract(table_Species_2$DADA2, "^.*(?=(_))") #remove GenBank identifier from species name
table_Species_2$blastn_ref <- str_extract(table_Species_2$blastn_ref, "^.*(?=(_))") #remove GenBank identifier from species name
table_Species_2[is.na(table_Species_2)] <- "NA"
table_Species_2$Full  <- ifelse(table_Species_2$DADA2=="NA", table_Species_2$blastn_ref, table_Species_2$Full) #add blastn results against custom reference database
table_Species_2$Full  <- ifelse(table_Species_2$Full=="Skeletonema pseudocostatum", table_Species_2$blastn_GenBank, 
                                       ifelse(table_Species_2$Full=="NA", table_Species_2$blastn_GenBank, table_Species_2$Full)) #add blastn results against GenBank
table_Species_2[is.na(table_Species_2)] <- "NA"

# fix names so it will be recognized by worms
#Species_FULL<-data.frame("Original_full"=unique(table_Species_2$Full))
openxlsx::write.xlsx(as.data.frame(unique(table_Species_2$Full)), paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/Species_FULL_WORMS_Input.xlsx"))

# use https://www.marinespecies.org/aphia.php?p=match to match full taxonomy
worms_names <- readxl::read_excel(paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/Species_FULL_WORMS_Output.xlsx"),sheet = 'WoRMS match')
worms_names_exact <- worms_names[worms_names$`Match type` %in% "exact",]
worms_names_exact <- worms_names_exact %>% filter(str_detect(str_trim(ScientificName_accepted), "\\s+"))
Classification_WoRMS <- select(worms_names_exact, AphiaID_accepted, Kingdom, Phylum, Class, Order, Family, Genus, ScientificName_accepted)
colnames(Classification_WoRMS) <- c("species_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus" , "Species")
Classification_WoRMS <- Classification_WoRMS[!duplicated(Classification_WoRMS$Species),]

##Alternative if you forgot to activate Classification in Worms
#worms_names_exact <- worms_names[worms_names$`Match type` %in% "exact",]
#worms_names_exact <- worms_names_exact %>% filter(str_detect(str_trim(ScientificName_accepted), "\\s+"))
#Species <- unique(table_Species_2$Full)
#Species <- Species[!Species == 'NA']
#write_tsv(as.data.frame(Species),paste0(proj.path,"/MiFish_UE-S_concatenated/results/Species_FULL_WORMS_Input.tsv"),col_names = F)
#Worms_id <- get_wormsid(Species, marine_only = F)
#Worms_id <- unique(worms_names_exact$AphiaID_accepted)
#Worms_id <- na.omit(Worms_id)
#Classification_worms <- classification(Worms_id, db = 'worms')
#Classification_worms <- cbind(Classification_worms)
#Classification_WoRMS <- select(Classification_worms, species_id, kingdom, phylum, class, order, family, genus, species)
#colnames(Classification_WoRMS) <- c("species_id", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus" , "Species")
#Classification_WoRMS <- Classification_WoRMS[!duplicated(Classification_WoRMS$Species),]

#set unaccepted names to acepted names
for (r in 1:nrow(table_Species_2)){
  if(table_Species_2$Full[r] %in% worms_names_exact$`unique(table_Species_2$Full)`){
    table_Species_2$Full[r] <- worms_names_exact$ScientificName_accepted[which(worms_names_exact$`unique(table_Species_2$Full)`==table_Species_2$Full[r])]
  }
}

#remove species with identical 12S sequences 
Identical_12S <- c("Eutrigla gurnardus", "Chelidonichthys lucerna","Chelidonichthys cuculus","Chelidonichthys spinosus",
                   "Hyperoplus lanceolatus","Hyperoplus immaculatus", "Ammodytes tobianus", "Ammodytes marinus",
                   "Alosa alosa", "Alosa fallax")
for (r in 1:nrow(table_Species_2)){
  if(table_Species_2$Full[r] %in% Identical_12S){
    table_Species_2$Full[r] <- "NA"
  }
}

select_columns<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus" , "Species")

#add blastn results to the unconcatenated table
table_raw[table_raw=="NA"] = NA 
table_raw[table_raw=="Alosa fallax_(MT410940)"] = NA  
table_raw$DADA2 <- table_raw$Species
table_raw$blastn_ref <- table_Species_2$blastn_ref
table_raw$blastn_GenBank <- table_Species_2$blastn_GenBank
table_raw$Full <- table_Species_2$Full
table_raw[which(table_raw$DADA2=="Skeletonema pseudocostatum_(MK372941)"),select_columns]<-NA

# Find rows where all columns match the specified values to remove inaccurate classification during the dada2 assignment
matching_rows <- which(table_raw$Kingdom == "Animalia" &
                         (is.na(table_raw$Phylum) | table_raw$Phylum == "Chordata") &
                         (is.na(table_raw$Class) | table_raw$Class == "Actinopteri") &
                         is.na(table_raw$Order) &
                         is.na(table_raw$Family) &
                         is.na(table_raw$Genus) &
                         is.na(table_raw$Species))

# Set all matching rows in all columns to "NA"
table_raw[matching_rows,select_columns] <- NA


# Set new classification based on NCBI and WoRMS
for (r in 1:nrow(table_raw)){
  if(table_raw$Full[r] %in% Classification_NCBI$Species){
    table_raw[r,select_columns] <- Classification_NCBI[which(Classification_NCBI$Species==table_raw$Full[r]),select_columns]
  }
}

for (r in 1:nrow(table_raw)){
  if(table_raw$Full[r] %in% Classification_WoRMS$Species){
    table_raw[r,select_columns] <- Classification_WoRMS[which(Classification_WoRMS$Species==table_raw$Full[r]),select_columns]
  }
}

#Set Eukaryota to Animalia, Plantae, Chromista, ... 
Animalia_Phyla <- c("Arthropoda", "Mollusca", "Chordata", "Chordata",
                    "Nematoda", "Annelida", "Cnidaria", "Porifera", 
                    "Echinodermata", "Bryozoa", "Rotifera", "Nemertea", "Tardigrada")
for (r in 1:nrow(table_raw)){
  if(table_raw$Phylum[r] %in% Animalia_Phyla){
    table_raw$Kingdom[r] <- "Animalia"
  }
}

Chromista_Phyla <- c("Bacillariophyta", "Ochrophyta", "Cryptophyta", "Haptophyta",
               "Heliozoa", "Chromeridophyta", "Ciliophora", "Myzozoa",
               "Bigyra", "Ochrophyta", "Oomycota", "Cercozoa", "Foraminifera", "Radiozoa")
for (r in 1:nrow(table_raw)){
    if(table_raw$Phylum[r] %in% Chromista_Phyla){
      table_raw$Kingdom[r] <- "Chromista"
    }
}

Plantae_Phyla <- c("Chlorophyta", "Glaucophyta", "Rhodelphidia", "Rhodophyta",
                   "Prasinodermophyta", "Anthocerotophyta", "Bryophyta",
                   "Charophyta", "Marchantiophyta", "Tracheophyta",
                   "Streptophyta")
for (r in 1:nrow(table_raw)){
    if(table_raw$Phylum[r] %in% Plantae_Phyla){
      table_raw$Kingdom[r] <- "Plantae"
    }
}

#For values where the full taxonomy is empty but there's still a higher classification add the higher classification with SP to the full column
find_first_non_na_column <- function(df) {
  non_na_columns <- c()
  for (i in 1:nrow(df)) {
    if (df$Full[i] %in% c("NA",NA)) {
      for (column in names(df)[(ncol(df)-10):ncol(df)]) {
        if(is.na(df$Family[i])){next}
        if (column != "Full" && df[i, column] %in% c("NA",NA)) {
          NOT_NA_column<-which(colnames(df)==column)-1
          df[i,"Full"]<-paste0(df[i,NOT_NA_column]," sp.")
          print(paste0("On row ",i," Full was changed to the ",colnames(df)[NOT_NA_column]," column + sp."))
          break
        }
      }
    }
  }
  return(df)
}

table_raw <- find_first_non_na_column(table_raw)
table_raw[is.na(table_raw)] <- "NA"
table_raw$Species <- table_raw$Full
table_raw$Full <- NULL

openxlsx::write.xlsx(as.data.frame(table_raw), 
           paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/table_unrarefied_raw_Full_TaxAss_WoRMS.xlsx"))


Non_Fish_Chordata <- c("Appendicularia", "Ascidiacea", "Ascidiacea", "Thaliacea",
                       "Amphibia", "Aves", "Aves", "Reptilia", "Reptilia", "Mammalia")

table_raw_fish <- table_raw[table_raw$Kingdom %in% "Animalia",]
table_raw_fish <- table_raw_fish[table_raw_fish$Phylum %in% "Chordata",]
table_raw_fish <- table_raw_fish[!table_raw_fish$Class %in% Non_Fish_Chordata,]

fish_classes <- c(unique(table_raw_fish$Class))

#source https://www.britannica.com/animal/fish/Annotated-classification
fish_classes_Joran<-c( "Agnatha" #superclass
                 ,"Myxini" #class
                 ,"Petromyzonti" #class
                 ,"Elasmobranchii"#class
                 ,"Neoselachii" #subclass
                 ,"Batoidea" #infraclass
                 ,"Selachii" #infraclass
                 ,"Holocephali" #class
                 ,"Osteichthyes" #superclass
                 ,"Actinopterygii" #gigaclass
                 ,"Actinopteri" #superclass
                 ,"Chondrostei" #class
                 ,"Holostei"    #class   
                 ,"Teleostei"  #class
                 ,"Sarcopterygii" #gigaclass
                 ,"Coelacanthimorpha" #class
                 ,"Dipneusti" #class          
)

Pelagic_Fish <- c("Alosa fallax",
                  "Aphia minuta", 
                  "Atherina presbyter",
                  "Belone belone",
                  "Chelon auratus",
                  "Engraulis encrasicolus",
                  "Osmerus eperlanus",
                  "Salmo trutta",
                  "Sardina pilchardus",
                  "Scomber scombrus",
                  "Sprattus sprattus",
                  "Trachurus trachurus")

saveRDS(fish_classes, file = paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_classes.rds"))
saveRDS(fish_classes_Joran, file = paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/REnvironment/fish_classes_Joran.rds"))
saveRDS(Pelagic_Fish, file = paste0(proj.path,"/MiFish_UE-S_concatenated/results_microDecon/R_Environment/Pelagic_Fish.rds"))
