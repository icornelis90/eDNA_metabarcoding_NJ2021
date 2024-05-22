-#title: "Compare the Alpha Diversity when using Decontam (un-rarefied) and microDecon"
#subtitle= "For microDecon we also compare non-rarefied vs coverage_based_rarefaction"
#author: "Isolde Cornelis"
#date: "28/02/2024"

################################################################################


proj.path.12S <- here("/home/genomics/icornelis/02_ZEROimpact/01_12S/NJ2021/MiFish-UE_run2")


microDecon_unrarefied <- readRDS(paste0(proj.path.12S,
                                        "/MiFish_UE-S_concatenated/results_microDecon/R_Environment/microDecon_unrarefied.rds"))
Decontam_unrarefied <- readRDS(paste0(proj.path.12S,
                                      "/MiFish_UE-S_concatenated/results_microDecon/R_Environment/Decontam_unrarefied.rds"))
microDecon_rarefied <- readRDS(paste0(proj.path.12S,
                                      "/MiFish_UE-S_concatenated/results_microDecon/R_Environment/microDecon_rarefied.rds"))
morph_unrarefied <- readRDS(paste0(proj.path.12S,
                                   "/MiFish_UE-S_concatenated/results_microDecon/R_Environment/Morph_unrarefied.rds"))
morph_rarefied <- readRDS(paste0(proj.path.12S,
                                 "/MiFish_UE-S_concatenated/results_microDecon/R_Environment/Morph_rarefied.rds"))


##create the sample dataframe
smpl_all <- bind_rows(Decontam_unrarefied$smpl_Decontam_unrarefied_Fish,
                      microDecon_unrarefied$smpl_microDecon_unrarefied_Fish,
                      microDecon_rarefied$smpl_micorDecon_rarefied_Fish,
                      morph_unrarefied$smpl_morph_unrarefied_Fish,
                      morph_rarefied$smpl_morph_rarefied_Fish,
                      Decontam_unrarefied$smpl_Decontam_unrarefied_Inv,
                      microDecon_unrarefied$smpl_microDecon_unrarefied_Inv,
                      microDecon_rarefied$smpl_micorDecon_rarefied_Inv,
                      morph_unrarefied$smpl_morph_unrarefied_Inv,
                      morph_rarefied$smpl_morph_rarefied_Inv)

smpl_all$Method <- c(rep("Decontam_non-rarefied", nrow(Decontam_unrarefied$smpl_Decontam_unrarefied_Fish)),
                     rep("microDecon_non-rarefied", nrow(microDecon_unrarefied$smpl_microDecon_unrarefied_Fish)),
                     rep("microDecon_rarefied", nrow(microDecon_rarefied$smpl_micorDecon_rarefied_Fish)), 
                     rep("Morphology_non-rarefied", nrow(morph_unrarefied$smpl_morph_unrarefied_Fish)),
                     rep("Morphology_rarefied", nrow(morph_rarefied$smpl_morph_rarefied_Fish)),
                     rep("Decontam_non-rarefied", nrow(Decontam_unrarefied$smpl_Decontam_unrarefied_Inv)),
                     rep("microDecon_non-rarefied", nrow(microDecon_unrarefied$smpl_microDecon_unrarefied_Inv)),
                     rep("microDecon_rarefied", nrow(microDecon_rarefied$smpl_micorDecon_rarefied_Inv)), 
                     rep("Morphology_non-rarefied", nrow(morph_unrarefied$smpl_morph_unrarefied_Inv)),
                     rep("Morphology_rarefied", nrow(morph_rarefied$smpl_morph_rarefied_Inv)))

repFish <- sum(nrow(Decontam_unrarefied$smpl_Decontam_unrarefied_Fish),
               nrow(microDecon_unrarefied$smpl_microDecon_unrarefied_Fish),
               nrow(microDecon_rarefied$smpl_micorDecon_rarefied_Fish),
               nrow(morph_unrarefied$smpl_morph_unrarefied_Fish),
               nrow(morph_rarefied$smpl_morph_rarefied_Fish))
repInv <- sum(nrow(Decontam_unrarefied$smpl_Decontam_unrarefied_Inv),
              nrow(microDecon_unrarefied$smpl_microDecon_unrarefied_Inv),
              nrow(microDecon_rarefied$smpl_micorDecon_rarefied_Inv),
              nrow(morph_unrarefied$smpl_morph_unrarefied_Inv),
              nrow(morph_rarefied$smpl_morph_rarefied_Inv))

smpl_all$Organism <- c(rep("Fish", repFish), rep("Invertebrates", repInv))
rownames(smpl_all) <- paste(smpl_all$Niskin.sample, smpl_all$Method, smpl_all$Organism, sep="_")

##Prepare the abundance table for the combined fish data
ps_Fish <- merge(Decontam_unrarefied$ps_Decontam_unrarefied_Fish, 
                 microDecon_unrarefied$ps_microDecon_unrarefied_Fish, by.x=0, by.y=0, all=T)
ps_Fish <- merge(ps_Fish, 
                 microDecon_rarefied$ps_micorDecon_rarefied_Fish, by.x=1, by.y=0, all=T)
ps_Fish <- merge(ps_Fish, 
                 morph_unrarefied$ps_morph_unrarefied_Fish, by.x=1, by.y=0, all=T)
ps_Fish <- merge(ps_Fish, 
                 morph_rarefied$ps_morph_rarefied_Fish, by.x=1, by.y=0, all=T)

ps_Inv <- merge(Decontam_unrarefied$ps_Decontam_unrarefied_Inv, 
                 microDecon_unrarefied$ps_micorDecon_unrarefied_Inv, by.x=0, by.y=0, all=T)
ps_Inv <- merge(ps_Inv, 
                 microDecon_rarefied$ps_micorDecon_rarefied_Inv, by.x=1, by.y=0, all=T)
ps_Inv <- merge(ps_Inv, 
                 morph_unrarefied$ps_morph_unrarefied_Inv, by.x=1, by.y=0, all=T)
ps_Inv <- merge(ps_Inv, 
                 morph_rarefied$ps_morph_rarefied_Inv, by.x=1, by.y=0, all=T)

ps_all <- merge(ps_Fish, ps_Inv, by.x=1, by.y=1, all=T)
rownames(ps_all) <- ps_all$Row.names
ps_all$Row.names <- NULL
ps_all[is.na(ps_all)] <- 0
colnames(ps_all) <- rownames(smpl_all)
Taxonomy_all <- as.matrix(rownames(ps_all))
rownames(Taxonomy_all) <- rownames(ps_all)

ps_all_phylo <- phyloseq(otu_table(ps_all, taxa_are_rows = TRUE), sample_data(smpl_all),tax_table(Taxonomy_all))
p_all <- plot_richness(ps_all_phylo, x='Zone', measures=c("Observed", "Shannon")) + 
  geom_boxplot(outlier.shape = 16, outlier.size = 2, aes(fill=Method)) + 
  scale_fill_manual(values=c("darkolivegreen3", "lightblue3", "lightblue4", "palevioletred2", "palevioletred4")) +
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  facet_grid(variable~Organism, scales = "free")
p_all$layers <- p_all$layers[-1]
p_all$data$Zone <- factor(p_all$data$Zone, levels=unique(smpl_all$Zone))
p_all$data$Method <- factor(p_all$data$Method, levels=unique(smpl_all$Method))
p_all

ps_phylo_microDecon <- phyloseq(otu_table(ps_all[,!grepl("Decontam", colnames(ps_all))], taxa_are_rows = TRUE),
                                sample_data(smpl_all[!grepl("Decontam", smpl_all$Method),]),
                                tax_table(Taxonomy_all))
p_microDecon <- plot_richness(ps_phylo_microDecon, x='Zone', measures=c("Observed")) + 
  geom_boxplot(outlier.shape = 16, outlier.size = 2, aes(fill=Method)) + 
  scale_fill_manual(values=c("lightblue3", "lightblue4", "palevioletred2", "palevioletred4"))+
  scale_x_discrete(labels=c("Coast", "Transition","Offshore")) +
  theme(axis.text.x=element_text(angle = 0, vjust = 0.5, hjust=0.5, color = c("limegreen","slateblue","darkorange"))) +
  labs(y = "Richness")+
  facet_grid(variable~Organism, scales = "free")
p_microDecon$layers <- p_microDecon$layers[-1]
p_microDecon$data$Zone <- factor(p_microDecon$data$Zone, levels=unique(smpl_all$Zone))
p_microDecon$data$Method <- factor(p_microDecon$data$Method, levels=unique(smpl_all$Method))
p_microDecon
