#title: "Creating Stacked Barplots with the un-rarefied data"
#author: "Isolde Cornelis"
#date: "16/11/2023"

##Stacked Barplot, x = Negative Control Samples, y = number of reads (unrarefied data), fill = human, Fish or others
#load data
proj.path <- here("/home/genomics/icornelis/02_ZEROimpact/01_12S/NJ2021/MiFish-UE_run2")
table_unrarefied_all <- readxl::read_excel(paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/table_unrarefied_AllASVs_FullTaxonomicAssignment.xlsx"))
#table_unrarefied_all <- readxl::read_excel("/home/genomics/icornelis/02_ZEROimpact/01_12S/NJ2021/MiFish-UE_run1/Step3_DADA2/Results/Taxonomic_assignment/seq_fish_order_v1.xlsx",sheet = "Sheet2")
env <- read.csv(paste0(proj.path,"/Step5_Statistics/environmental_data.csv"),  header=TRUE, sep=";")
env$Environment_color <- ifelse(env$Zones=="Coastal","limegreen", 
                                ifelse(env$Zones=="zone1", "slateblue", 
                                       ifelse(env$Zones=="zone2","darkorange","red")))

#update the env data to keep only the samples present in the ASV table with taxonomic assignment
keep_samples <- c(colnames(table_unrarefied_all))
env_all <- env[env$Niskin.sample %in% keep_samples,]
colorder_all <- env_all$Niskin.sample

#Select the ASVs assigned to Fish, Human and others
#select ASVs assigned to fish
fish_classes <- readRDS(file = paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_classes.rds"))
freshwater_fish <- readRDS(file = paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_Freshwater.rds"))
table_FishASVs <- as.data.frame(table_unrarefied_all[table_unrarefied_all$Class %in% fish_classes,])
table_FishASVs <- as.data.frame(table_FishASVs[!table_FishASVs$Species %in% c(freshwater_fish, "NA"),])
FishASVs <- c(table_FishASVs$ASV)
table_HumanASVs <- as.data.frame(table_unrarefied_all[table_unrarefied_all$Species %in% c("Homo sapiens"),])
HumanASVs <- c(table_HumanASVs$ASV)
table_NoTAASVs <- as.data.frame(table_unrarefied_all[table_unrarefied_all$Species %in% c("NA"),])
NoTAASVs <- c(table_NoTAASVs$ASV)
table_OtherASVs <-  as.data.frame(table_unrarefied_all[!table_unrarefied_all$ASV %in% c(FishASVs, HumanASVs, NoTAASVs),])

#merge
taxo <- 'Species'
merged_data_rarefied_fish <- aggregate(table_FishASVs[,1:(ncol(table_FishASVs)-11)], by= list(as.factor(table_FishASVs[,taxo])),FUN=sum)
merged_data_rarefied_other <- aggregate(table_OtherASVs[,1:(ncol(table_OtherASVs)-11)], by= list(as.factor(table_OtherASVs[,taxo])),FUN=sum)
rownames(merged_data_rarefied_fish) <- merged_data_rarefied_fish$Group.1
rownames(merged_data_rarefied_other) <- merged_data_rarefied_other$Group.1
merged_data_rarefied_fish <- merged_data_rarefied_fish[,colorder_all]
merged_data_rarefied_other <- merged_data_rarefied_other[,colorder_all]

#rarecurve unrarefied ASVs assigned to fish
seqtab.FishASVs <- as.data.frame(t(table_FishASVs[,1:(ncol(table_FishASVs)-11)]))
seqtab.FishASVs <- seqtab.FishASVs[colorder_all,]
raremax = 10000
out <- rarecurve(seqtab.FishASVs, ylab = "ASVs Assigned to Fish", main = "Rarecurve of unrarefied samples after taxonomic assignment", col = as.vector(env_all$Environment_color), label = FALSE, step =10)
pars <- expand.grid(col = as.vector(env_all$Environment_color), lty = "solid", lwd = 1, 
                    stringsAsFactors = FALSE)
Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
Smax <- sapply(out, max)
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Sample Size",
     ylab = "ASVs Assigned to Fish", type = "n",  main = "Rarecurve of unrarefied samples after taxonomic assignment")
abline(v = raremax)
for (i in seq_along(out)) {
  N <- attr(out[[i]], "Subsample")
  with(pars, lines(N, out[[i]], col = col[i], lty = lty[i], lwd = lwd[i]))
}

#create sequence table with only Negative Control Samples
seqtab_FishASVs <- as.data.frame(t(table_FishASVs[,1:(ncol(table_FishASVs)-11)]))
seqtab_HumanASVs <- as.data.frame(t(table_HumanASVs[,1:(ncol(table_HumanASVs)-11)]))
seqtab_NoTAASVs <- as.data.frame(t(table_NoTAASVs[,1:(ncol(table_NoTAASVs)-11)]))
seqtab_OtherASVs <- as.data.frame(t(table_OtherASVs[,1:(ncol(table_OtherASVs)-11)]))

seqtab_FishASVs <- seqtab_FishASVs[colorder_all,]
rownames(seqtab_FishASVs) <- ifelse(rownames(seqtab_FishASVs)=="ftWT2bis_4","ftWT2bis_3", 
                                    ifelse(rownames(seqtab_FishASVs)=="ftWT2bis_7", "ftWT2bis_4", 
                                           ifelse(rownames(seqtab_FishASVs)=="ftWT2bis_10", "ftWT2bis_5",  rownames(seqtab_FishASVs))))
seqtab_HumanASVs <- seqtab_HumanASVs[colorder_all,]
seqtab_NoTAASVs <- seqtab_NoTAASVs[colorder_all,]
seqtab_OtherASVs <- seqtab_OtherASVs[colorder_all,]

seqtab_ctrl_FishASVs <- seqtab_FishASVs %>% filter(grepl("neg", rownames(seqtab_FishASVs)))
seqtab_ctrl_HumanASVs <- seqtab_HumanASVs %>% filter(grepl("neg", rownames(seqtab_HumanASVs)))
seqtab_ctrl_NoTAASVs <- seqtab_NoTAASVs %>% filter(grepl("neg", rownames(seqtab_NoTAASVs)))
seqtab_ctrl_OtherASVs <- seqtab_OtherASVs %>% filter(grepl("neg", rownames(seqtab_OtherASVs)))

colnames(seqtab_ctrl_FishASVs ) <- table_FishASVs$ASV
colnames(seqtab_ctrl_HumanASVs ) <- table_HumanASVs$ASV
colnames(seqtab_ctrl_NoTAASVs ) <- table_NoTAASVs$ASV
colnames(seqtab_ctrl_OtherASVs ) <- table_OtherASVs$ASV

#create Stacked Barplot with total reads
NegCtrlSample <- rep(c(rownames(seqtab_ctrl_FishASVs)),4)
TypeOfASV <- c(rep("Fish", 22), rep("Human",22), rep("Other", 22), rep("NA", 22))
Reads_FishASVs <- rowSums(seqtab_ctrl_FishASVs)
Reads_HumanASVs <- rowSums(seqtab_ctrl_HumanASVs)
Reads_NoTAASVs <- rowSums(seqtab_ctrl_NoTAASVs)
Reads_OtherASVs <- rowSums(seqtab_ctrl_OtherASVs)
Reads <- c(Reads_FishASVs, Reads_HumanASVs, Reads_OtherASVs, Reads_NoTAASVs)
data_Reads <- data.frame(NegCtrlSample,TypeOfASV,Reads)

ggplot(data_Reads, aes(fill=TypeOfASV, y=Reads, x=factor(NegCtrlSample, level = c(rownames(seqtab_ctrl_FishASVs))))) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(size=10, angle=90))+
  labs(x = "", fill = "ASV origin", title = "Number of Reads in Negative Control Samples, unrarefied data")

#create Stacked Barplot with total ASVs
ASVs_FishASVs <- rowSums(seqtab_ctrl_FishASVs !=0)
ASVs_HumanASVs <- rowSums(seqtab_ctrl_HumanASVs !=0)
ASVs_NoTAASVs <- rowSums(seqtab_ctrl_NoTAASVs !=0)
ASVs_OtherASVs <- rowSums(seqtab_ctrl_OtherASVs !=0)
ASVs <- c(ASVs_FishASVs, ASVs_HumanASVs, ASVs_OtherASVs, ASVs_NoTAASVs)
data_ASVs <- data.frame(NegCtrlSample,TypeOfASV,ASVs)

ggplot(data_Reads, aes(fill=TypeOfASV, y=ASVs, x=factor(NegCtrlSample, level = c(rownames(seqtab_ctrl_FishASVs))))) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(size=10, angle=90))+
  labs(x = "", fill = "ASV origin", title = "Number of ASVs in Negative Control Samples, unrarefied data")

#create Stacked Barplot with all samples - Reads
Sample <- rep(c(rownames(seqtab_FishASVs)), 4)
TypeOfASV <- c(rep("Fish", 148), rep("Human", 148), rep("Other", 148), rep("NA", 148))
Reads_FishASVs <- rowSums(seqtab_FishASVs)
Reads_HumanASVs <- rowSums(seqtab_HumanASVs)
Reads_NoTAASVs <- rowSums(seqtab_NoTAASVs)
Reads_OtherASVs <- rowSums(seqtab_OtherASVs)
Reads <- c(Reads_FishASVs, Reads_HumanASVs, Reads_OtherASVs, Reads_NoTAASVs)
data_Reads <- data.frame(Sample,TypeOfASV,Reads)

ggplot(data_Reads, aes(fill=TypeOfASV, y=Reads, x=factor(Sample, level = c(rownames(seqtab_FishASVs))))) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(size=10, angle=90, colour = c(env_all$Zone_color)))+
  labs(x = "", fill = "ASV origin", title = "Number of Reads in All Samples, unrarefied data")

#create Stacked Barplot with all samples - ASVs
ASVs_FishASVs <- rowSums(seqtab_FishASVs !=0)
ASVs_HumanASVs <- rowSums(seqtab_HumanASVs !=0)
ASVs_NoTAASVs <- rowSums(seqtab_NoTAASVs !=0)
ASVs_OtherASVs <- rowSums(seqtab_OtherASVs !=0)
ASVs <- c(ASVs_FishASVs, ASVs_HumanASVs, ASVs_OtherASVs, ASVs_NoTAASVs)
data_ASVs <- data.frame(Sample,TypeOfASV,ASVs)

ggplot(data_Reads, aes(fill=TypeOfASV, y=ASVs, x=factor(Sample, level = c(rownames(seqtab_FishASVs))))) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(size=10, angle=90, colour = c(env_all$Zone_color)))+
  labs(x = "", fill = "ASV origin", title = "Number of ASVs in All Samples, unrarefied data")

#unique fish species and others
Fish_FishASVs <- colSums(merged_data_rarefied_fish[,1:ncol(merged_data_rarefied_fish)] !=0)
Other_OtherASVs <- colSums(merged_data_rarefied_other[,1:ncol(merged_data_rarefied_other)] !=0)

Samples <- c(rownames(seqtab_FishASVs))
data_summary <- data.frame(Reads_FishASVs, ASVs_FishASVs, Fish_FishASVs,
                              Reads_HumanASVs, ASVs_HumanASVs, 
                              Reads_OtherASVs, ASVs_OtherASVs, Other_OtherASVs,
                              Reads_NoTAASVs, ASVs_NoTAASVs)
write.xlsx(data_summary, paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/raw_data_reads_ASVs.xlsx"), rowNames =T)
write.xlsx(merged_data_rarefied_fish, paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/Fish_merged_unrarefied.xlsx", rowNames = T))
