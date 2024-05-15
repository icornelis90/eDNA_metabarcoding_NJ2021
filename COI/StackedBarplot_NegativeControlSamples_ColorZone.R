#title: "Creating Stacked Barplots with the un-rarefied data"
#author: "Isolde Cornelis"
#date: "01/03/2024"

################################################################################

##Stacked Barplot, x = Negative Control Samples, y = number of reads (un-rarefied data), fill = human, Fish or others
#load data
proj.path <- here("/home/genomics/icornelis/02_ZEROimpact/02_COI/NJ2021")
table_unrarefied_all <- readxl::read_excel(paste0(proj.path,"/OWFvsCoastal_concatenated/results/Decontam/table_unrarefied_AllASVs_FullTaxonomicAssignment.xlsx"))
env <- read.csv(paste0(proj.path,"/Step5_Statistics/environmental_data.csv"),  header=TRUE, sep=";")
env$Zone_color <- ifelse(env$Zone=="Coast","limegreen", 
                                ifelse(env$Zone=="Transition", "slateblue", 
                                       ifelse(env$Zone=="Offshore","darkorange","red")))

#update the env data to keep only the samples present in the ASV table with taxonomic assignment
keep_samples <- c(colnames(table_unrarefied_all))
env_all <- env[env$Niskin.sample %in% keep_samples,]
colorder_all <- env_all$Niskin.sample

#Select the ASVs assigned to Fish, Human and others
#select ASVs assigned to fish
table_AnimaliaASVs <- as.data.frame(table_unrarefied_all[table_unrarefied_all$Kingdom %in% c("Animalia"),])
table_InvASVs <- as.data.frame(table_AnimaliaASVs[!table_AnimaliaASVs$Phylum %in% c("Chordata"),])
table_InvASVsNoSp <- as.data.frame(table_AnimaliaASVs[!grepl("sp", table_AnimaliaASVs$Full),])
table_NoTAASVs <- as.data.frame(table_unrarefied_all[table_unrarefied_all$Kingdom %in% c("NA"),])
NoTAASVs <- c(table_NoTAASVs$ASV)
table_ChromistaASVs <- as.data.frame(table_unrarefied_all[table_unrarefied_all$Kingdom %in% c("Chromista"),])
ChromistaASVs <- c(table_ChromistaASVs$ASV)
table_OtherASVs <-  as.data.frame(table_unrarefied_all[!table_unrarefied_all$ASV %in% c(AnimaliaASVs,NoTAASVs,ChromistaASVs),])

#merge
taxo <- 'Full'
merged_data_rarefied_animalia <- aggregate(table_AnimaliaASVs[,1:(ncol(table_AnimaliaASVs)-11)], by= list(as.factor(table_AnimaliaASVs[,taxo])),FUN=sum)
merged_data_rarefied_other <- aggregate(table_OtherASVs[,1:(ncol(table_OtherASVs)-11)], by= list(as.factor(table_OtherASVs[,taxo])),FUN=sum)
merged_data_rarefied_chromista <- aggregate(table_ChromistaASVs[,1:(ncol(table_ChromistaASVs)-11)], by= list(as.factor(table_ChromistaASVs[,taxo])),FUN=sum)
rownames(merged_data_rarefied_animalia) <- merged_data_rarefied_animalia$Group.1
rownames(merged_data_rarefied_other) <- merged_data_rarefied_other$Group.1
rownames(merged_data_rarefied_chromista) <- merged_data_rarefied_chromista$Group.1
merged_data_rarefied_animalia <- merged_data_rarefied_animalia[,colorder_all]
merged_data_rarefied_other <- merged_data_rarefied_other[,colorder_all]
merged_data_rarefied_chromista <- merged_data_rarefied_chromista[,colorder_all]

#rarecurve unrarefied ASVs assigned to fish
seqtab.AnimaliaASVs <- as.data.frame(t(table_AnimaliaASVs[,1:(ncol(table_AnimaliaASVs)-11)]))
seqtab.AnimaliaASVs <- seqtab.AnimaliaASVs[colorder_all,]
seqtab.AnimaliaASVs2 <- seqtab.AnimaliaASVs[!rowSums(seqtab.AnimaliaASVs) > 150000,]
env_all2 <- env_all[env_all$Niskin.sample %in% c(rownames(seqtab.AnimaliaASVs2)),]

raremax = 25000
out <- rarecurve(seqtab.AnimaliaASVs, ylab = "ASVs Assigned to Animalia", main = "Rarecurve of unrarefied samples after taxonomic assignment", col = as.vector(env_all$Zone_color), label = TRUE, step =10)
out <- rarecurve(seqtab.AnimaliaASVs2, ylab = "ASVs Assigned to Animalia", main = "Rarecurve of unrarefied samples after taxonomic assignment", col = as.vector(env_all2$Zone_color), label = FALSE, step =10)
pars <- expand.grid(col = as.vector(env_all2$Environment_color), lty = "solid", lwd = 1, 
                    stringsAsFactors = FALSE)
Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
Smax <- sapply(out, max)
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Sample Size",
     ylab = "ASVs Assigned to Animalia", type = "n",  main = "Rarecurve of unrarefied samples after taxonomic assignment")
abline(v = raremax)
for (i in seq_along(out)) {
  N <- attr(out[[i]], "Subsample")
  with(pars, lines(N, out[[i]], col = col[i], lty = lty[i], lwd = lwd[i]))
}

#create sequence table with only Negative Control Samples
seqtab_AnimaliaASVs <- as.data.frame(t(table_AnimaliaASVs[,1:(ncol(table_AnimaliaASVs)-11)]))
seqtab_ChromistaASVs <- as.data.frame(t(table_ChromistaASVs[,1:(ncol(table_ChromistaASVs)-11)]))
seqtab_NoTAASVs <- as.data.frame(t(table_NoTAASVs[,1:(ncol(table_NoTAASVs)-11)]))
seqtab_OtherASVs <- as.data.frame(t(table_OtherASVs[,1:(ncol(table_OtherASVs)-11)]))

seqtab_AnimaliaASVs <- seqtab_AnimaliaASVs[colorder_all,]
rownames(seqtab_AnimaliaASVs) <- ifelse(rownames(seqtab_AnimaliaASVs)=="ftWT2bis_4","ftWT2bis_3", 
                                    ifelse(rownames(seqtab_AnimaliaASVs)=="ftWT2bis_7", "ftWT2bis_4", 
                                           ifelse(rownames(seqtab_AnimaliaASVs)=="ftWT2bis_10", "ftWT2bis_5",  rownames(seqtab_AnimaliaASVs))))
seqtab_ChromistaASVs <- seqtab_ChromistaASVs[colorder_all,]
seqtab_NoTAASVs <- seqtab_NoTAASVs[colorder_all,]
seqtab_OtherASVs <- seqtab_OtherASVs[colorder_all,]

seqtab_ctrl_AnimaliaASVs <- seqtab_AnimaliaASVs %>% filter(grepl("neg", rownames(seqtab_AnimaliaASVs)))
seqtab_ctrl_ChromistaASVs <- seqtab_ChromistaASVs %>% filter(grepl("neg", rownames(seqtab_ChromistaASVs)))
seqtab_ctrl_NoTAASVs <- seqtab_NoTAASVs %>% filter(grepl("neg", rownames(seqtab_NoTAASVs)))
seqtab_ctrl_OtherASVs <- seqtab_OtherASVs %>% filter(grepl("neg", rownames(seqtab_OtherASVs)))

colnames(seqtab_ctrl_AnimaliaASVs ) <- table_AnimaliaASVs$ASV
colnames(seqtab_ctrl_ChromistaASVs ) <- table_ChromistaASVs$ASV
colnames(seqtab_ctrl_NoTAASVs ) <- table_NoTAASVs$ASV
colnames(seqtab_ctrl_OtherASVs ) <- table_OtherASVs$ASV

#create Stacked Barplot with total reads
NegCtrlSample <- rep(c(rownames(seqtab_ctrl_AnimaliaASVs)),4)
TypeOfASV <- c(rep("Animalia", 19), rep("Chromista", 19), rep("Other", 19), rep("NA", 19))
Reads_AnimaliaASVs <- rowSums(seqtab_ctrl_AnimaliaASVs)
Reads_ChromistaASVs <- rowSums(seqtab_ctrl_ChromistaASVs)
Reads_NoTAASVs <- rowSums(seqtab_ctrl_NoTAASVs)
Reads_OtherASVs <- rowSums(seqtab_ctrl_OtherASVs)
Reads <- c(Reads_AnimaliaASVs, Reads_ChromistaASVs, Reads_OtherASVs, Reads_NoTAASVs)
data_Reads <- data.frame(NegCtrlSample,TypeOfASV,Reads)

ggplot(data_Reads, aes(fill=TypeOfASV, y=Reads, x=factor(NegCtrlSample, level = c(rownames(seqtab_ctrl_AnimaliaASVs))))) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(size=10, angle=90))+
  labs(x = "", fill = "ASV origin", title = "Number of Reads in Negative Control Samples, unrarefied data")

#create Stacked Barplot with total ASVs
ASVs_AnimaliaASVs <- rowSums(seqtab_ctrl_AnimaliaASVs !=0)
ASVs_ChromistaASVs <- rowSums(seqtab_ctrl_ChromistaASVs !=0)
ASVs_NoTAASVs <- rowSums(seqtab_ctrl_NoTAASVs !=0)
ASVs_OtherASVs <- rowSums(seqtab_ctrl_OtherASVs !=0)
ASVs <- c(ASVs_AnimaliaASVs, ASVs_ChromistaASVs, ASVs_OtherASVs, ASVs_NoTAASVs)
data_ASVs <- data.frame(NegCtrlSample,TypeOfASV,ASVs)

ggplot(data_Reads, aes(fill=TypeOfASV, y=ASVs, x=factor(NegCtrlSample, level = c(rownames(seqtab_ctrl_AnimaliaASVs))))) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(size=10, angle=90))+
  labs(x = "", fill = "ASV origin", title = "Number of ASVs in Negative Control Samples, unrarefied data")

#create Stacked Barplot with all samples - Reads
NegCtrlSample <- rep(c(rownames(seqtab_AnimaliaASVs)), 4)
TypeOfASV <- c(rep("Animalia", 142), rep("Chromista", 142), rep("Other", 142), rep("NA", 142))
Reads_AnimaliaASVs <- rowSums(seqtab_AnimaliaASVs)
Reads_ChromistaASVs <- rowSums(seqtab_ChromistaASVs)
Reads_NoTAASVs <- rowSums(seqtab_NoTAASVs)
Reads_OtherASVs <- rowSums(seqtab_OtherASVs)
Reads <- c(Reads_AnimaliaASVs, Reads_ChromistaASVs, Reads_OtherASVs, Reads_NoTAASVs)
data_Reads <- data.frame(NegCtrlSample,TypeOfASV,Reads)

ggplot(data_Reads, aes(fill=TypeOfASV, y=Reads, x=factor(NegCtrlSample, level = c(rownames(seqtab_AnimaliaASVs))))) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(size=10, angle=90, colour = c(env_all$Zone_color)))+
  labs(x = "", fill = "ASV origin", title = "Number of Reads in All Samples, unrarefied data")

#create Stacked Barplot with all samples - ASVs
ASVs_AnimaliaASVs <- rowSums(seqtab_AnimaliaASVs !=0)
ASVs_ChromistaASVs <- rowSums(seqtab_ChromistaASVs !=0)
ASVs_NoTAASVs <- rowSums(seqtab_NoTAASVs !=0)
ASVs_OtherASVs <- rowSums(seqtab_OtherASVs !=0)
ASVs <- c(ASVs_AnimaliaASVs, ASVs_ChromistaASVs, ASVs_OtherASVs, ASVs_NoTAASVs)
data_ASVs <- data.frame(NegCtrlSample,TypeOfASV,ASVs)

ggplot(data_Reads, aes(fill=TypeOfASV, y=ASVs, x=factor(NegCtrlSample, level = c(rownames(seqtab_AnimaliaASVs))))) + 
  geom_bar(position="stack", stat="identity") +
  theme(axis.text.x = element_text(size=10, angle=90, colour = c(env_all$Zone_color)))+
  labs(x = "", fill = "ASV origin", title = "Number of ASVs in All Samples, unrarefied data")

#unique animalia species and others
Animalia_AnimliaASVs <- colSums(merged_data_rarefied_animalia[,1:ncol(merged_data_rarefied_animalia)] !=0)
Other_OtherASVs <- colSums(merged_data_rarefied_other[,1:ncol(merged_data_rarefied_other)] !=0)
Chromista_ChromistaASVs <- colSums(merged_data_rarefied_other[,1:ncol(merged_data_rarefied_other)] !=0)

Samples <- c(rownames(seqtab_AnimaliaASVs))
data_summary <- data.frame(Reads_AnimaliaASVs, ASVs_AnimaliaASVs, Animalia_AnimliaASVs,
                           Reads_ChromistaASVs, ASVs_ChromistaASVs, Chromista_ChromistaASVs,
                           Reads_OtherASVs, ASVs_OtherASVs, Other_OtherASVs,
                           Reads_NoTAASVs, ASVs_NoTAASVs)
openxlsx::write.xlsx(data_summary, paste0(proj.path,"/OWFvsCoastal_concatenated/results/Decontam/raw_data_reads_ASVs_v2.xlsx"), rowNames = T)
openxlsx::write.xlsx(merged_data_rarefied_animalia, paste0(proj.path,"/OWFvsCoastal_concatenated/results/Decontam/Species_merged_unrarefied_v2.xlsx"), rowNames = T)
