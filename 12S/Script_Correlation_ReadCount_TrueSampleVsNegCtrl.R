#title: "Check fish species abundance in true samples vs negative control samples"
#author: "Isolde Cornelis"
#date: "15/05/2024"

################################################################################

# make paths
proj.path <- here("/home/genomics/icornelis/02_ZEROimpact/01_12S/NJ2021/MiFish-UE_run2")

#upload data
table_unrarefied <- readxl::read_excel(paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/table_unrarefied_AllASVs_FullTaxonomicAssignment.xlsx"))

#create rarecurve with ASVs assigned to Fish at Species level to check the readdepth
fish_classes <- readRDS(file = paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_classes.rds"))
freshwater_fish <- readRDS(file = paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_Freshwater.rds"))
table_FishASVs <- as.data.frame(table_unrarefied[table_unrarefied$Class %in% c(fish_classes),])
table_FishASVs <- as.data.frame(table_FishASVs[!table_FishASVs$Species %in% c(freshwater_fish, "NA"),])
rownames(table_FishASVs) <- table_FishASVs$ASV
table_FishASVs[is.na(table_FishASVs)] <- 0

#sort by species
taxo <- "Species"
table_FishASVs <- aggregate(table_FishASVs[,1:(ncol(table_FishASVs)-11)],
                         by= list(as.factor(table_FishASVs[,taxo])),FUN=sum)
rownames(table_FishASVs) <- table_FishASVs$Group.1
table_FishASVs$Group.1 <- NULL
seqtab.FishASVs.merged <- as.data.frame(t(table_FishASVs))

options(scipen = 999)

#Dot plot negative control vs samples (Fish species)
Neg_ctrl <- seqtab.FishASVs.merged %>% filter(grepl("neg", rownames(seqtab.FishASVs.merged)))
Sample <- seqtab.FishASVs.merged %>% filter(!grepl("neg", rownames(seqtab.FishASVs.merged)))

Abundance <- as.data.frame(rownames(table_FishASVs))
Abundance$Samples <- colSums(Sample)
Abundance$Neg_ctrl <- colSums(Neg_ctrl)
colnames(Abundance) <- c("Species", "Sample", "Neg_ctrl")

p <- ggplot(Abundance, aes(Neg_ctrl, Sample))
p + geom_point()+
  geom_smooth(method = 'lm')+
  theme(axis.text.y = element_text(size=10, colour = "black") #species size
        , plot.title = element_text(size=15) #size of the title of the plot
        , legend.text = element_text(size = 10)
        , strip.text.x = element_text(size = 15)
        , axis.text.x = element_text(size=10, colour = "black")
  )+
  labs(title = "Read counts of each fish species in the samples vs negative controls",
       x = "log10(Read counts in negative controls)",
       y = "log10(Read counts in all samples)"
  )+
  theme_bw() +
  scale_y_log10(labels = scales::comma) +
  scale_x_log10(labels = scales::comma) +
  stat_cor(aes(Neg_ctrl, Sample), label.x = 3,method = "spearman", p.accuracy = 0.001) 

cor.test(Abundance$Neg_ctrl, Abundance$Sample, 
         method = "spearman")
