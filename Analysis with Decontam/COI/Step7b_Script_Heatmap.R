#title: "Heatmap of relative abundance of fish species"
#author: "Isolde Cornelis"
#date: "01/02/2024"

################################################################################

#set OS type for paths
if(.Platform$OS.type == "unix"){home<-"/home/"} else{
  home<-"//192.168.236.131/"}

libraries <- c("BiocManager"
               , "ggplot2"
               , "tidyr"
               , "taxonomizr"
               , "naturalsort"
               , "stringr"
               , "phylotools"
               , "scales"
               , "ggpattern"
               , "ggh4x"
               , "reshape2"
               , "grid"
               , "ggdendro"
               , "dendextend"
               , "here"
               , "dplyr"
               , "vegan"
)

for (opties in libraries){
  
  if (opties %in% installed.packages()){
    
    library(opties,character.only = TRUE)
    
  } else {install.packages(opties,repos = "http://cran.us.r-project.org")
    
    library(opties,character.only = TRUE)
  }
}

# make paths
proj.path <- here("/home/genomics/icornelis/02_ZEROimpact/02_COI/NJ2021")

#upload data
table_rarefied_all <- readxl::read_excel(paste0(proj.path,"/OWFvsCoastal_concatenated/results/Decontam/table_unrarefied_concatenated_CleanedASVs_FullTaxonomicAssignment_WithField.xlsx"))
env <- read.csv(paste0(proj.path,"/Step5_Statistics/environmental_data.csv"),  header=TRUE, sep=";")

#add color to environmental data to create the plot 
env$Environment_color <- ifelse(env$Zones=="Coastal","limegreen", 
                                ifelse(env$Zones=="zone1", "slateblue", 
                                       ifelse(env$Zones=="zone2","darkorange","red")))

#remove samples from environmental data that are absent from the rarefied data
keep_samples <- c(colnames(table_rarefied_all))
env_all <- env[env$Niskin.sample %in% keep_samples,]
env_all_samples <- env_all[!env_all$Zones == "neg_control",]
colorder_all <- c(env_all_samples$Niskin.sample)

#remove samples from environmental data that are absent from the rarefied data
keep_samples <- c(colnames(table_rarefied_all))
env_rarefied <- env[env$Niskin.sample %in% keep_samples,]
env_rarefied <- env_rarefied[!env_rarefied$Zones == "neg_control",]
colorder_rarefied <- c(env_rarefied$Niskin.sample)

#select ASVs that are assigned ASVs and ASVs assigned to Animalia
table_rarefied_allAnimaliaASVs <- as.data.frame(table_rarefied_all[table_rarefied_all$Kingdom %in% c("Animalia"),])
table_rarefied_allAnimaliaASVs <- as.data.frame(table_rarefied_allAnimaliaASVs[!table_rarefied_allAnimaliaASVs$Full %in% "NA",])
table_rarefied_ChordataASVs <- as.data.frame(table_rarefied_allAnimaliaASVs[table_rarefied_allAnimaliaASVs$Phylum %in% "Chordata",])
table_rarefied_allAnimaliaASVs <- as.data.frame(table_rarefied_allAnimaliaASVs[!table_rarefied_allAnimaliaASVs$Phylum %in% "Chordata",])

#3: sort by species
taxo <- "Full"
merged_data_rarefied_allAnimalia <- aggregate(table_rarefied_allAnimaliaASVs[,1:(ncol(table_rarefied_allAnimaliaASVs)-11)], by= list(as.factor(table_rarefied_allAnimaliaASVs[,taxo])),FUN=sum)
merged_data_rarefied_Chordata <- aggregate(table_rarefied_ChordataASVs[,1:(ncol(table_rarefied_ChordataASVs)-11)], by= list(as.factor(table_rarefied_ChordataASVs[,taxo])),FUN=sum)

#Set the tax names as rownames
rownames(merged_data_rarefied_allAnimalia)<-merged_data_rarefied_allAnimalia$Group.1
merged_data_rarefied_allAnimalia$Group.1<-NULL
unfound_species <- c(rownames(merged_data_rarefied_allAnimalia[rowSums(merged_data_rarefied_allAnimalia)==0,]))
remove_samples <- c(colnames(merged_data_rarefied_allAnimalia[,colSums(merged_data_rarefied_allAnimalia)==0]))
merged_data_rarefied_allAnimalia <- merged_data_rarefied_allAnimalia[!rownames(merged_data_rarefied_allAnimalia) %in% unfound_species,]
#merged_data_rarefied_allfish <- merged_data_rarefied_allfish[,!colnames(merged_data_rarefied_allassigned) %in% remove_samples]
merged_data_rarefied_allAnimalia<- merged_data_rarefied_allAnimalia[,colorder_all]
#a <- ifelse(rev(unique(rownames(merged_data_rarefied_allfish))) %in% unfound_species,"grey","black")

rownames(merged_data_rarefied_Chordata)<-merged_data_rarefied_Chordata$Group.1
merged_data_rarefied_Chordata$Group.1<-NULL
unfound_species <- c(rownames(merged_data_rarefied_Chordata[rowSums(merged_data_rarefied_Chordata)==0,]))
remove_samples <- c(colnames(merged_data_rarefied_Chordata[,colSums(merged_data_rarefied_Chordata)==0]))
merged_data_rarefied_Chordata <- merged_data_rarefied_Chordata[!rownames(merged_data_rarefied_Chordata) %in% unfound_species,]
#merged_data_rarefied_allfish <- merged_data_rarefied_allfish[,!colnames(merged_data_rarefied_allassigned) %in% remove_samples]
merged_data_rarefied_Chordata<- merged_data_rarefied_Chordata[,colorder_all]
#a <- ifelse(rev(unique(rownames(merged_data_rarefied_allfish))) %in% unfound_species,"grey","black")

#Order Species (Animalia) according to abundance
#sort according to abundance
merged_data_rarefied_allAnimalia_order <- merged_data_rarefied_allAnimalia
merged_data_rarefied_allAnimalia_order$Order <- rowSums(merged_data_rarefied_allAnimalia)
OrderAnimalia_all <- merged_data_rarefied_allAnimalia[order(merged_data_rarefied_allAnimalia_order$Order, decreasing = T),]
OrderAnimalia_all_top50 <- OrderAnimalia_all[1:50,]

#Location dendro based on transformed data
transformed_location <- decostand(decostand(merged_data_rarefied_allAnimalia, method="tot"), method="max")
Location_dendro <- as.dendrogram(hclust(d=dist(x = (t(transformed_location))), method = "ward.D"))
Location_order <- order.dendrogram(Location_dendro)
lable_color <- env_rarefied[Location_order,]
labels_colors(Location_dendro) <- c(lable_color$Zone_color)
dendro_plot <- ggdendrogram(data = Location_dendro, rotate = F)
par(mar = c(8, 3, 2, 2))
plot(Location_dendro)

#Animalia dendro
#merged_data_rarefied_allAnimalia_cluster <- (OrderAnimalia_all_top50)^0.25
OrderAnimalia_all_top50_cluster <- decostand(decostand(OrderAnimalia_all_top50, method="tot"), method="max")
Animalia_dendro_top50 <- as.dendrogram(hclust(d=dist(x = OrderAnimalia_all_top50_cluster),method = "ward.D"))
dendro_plot_top50 <- ggdendrogram(data = Animalia_dendro_top50, rotate = T)
print(dendro_plot_top50)

#Order Animalia according to dendrogram
Animalia_order_top50 <- order.dendrogram(Animalia_dendro_top50)
#OrderAnimalia_all_top50 <- merged_data_rarefied_allAnimalia[Animalia_order,]
OrderAnimalia_all_top50 <- OrderAnimalia_all_top50[Animalia_order_top50,]

#Calculate the relative abundance of each species per sample
#relative_abundance_all <- data.frame(matrix(nrow=nrow(OrderAnimalia_all), ncol=ncol(OrderAnimalia_all)
#                                            , dimnames = list(c(rownames(OrderAnimalia_all)), c(colnames(OrderAnimalia_all)))))
#for(i in 1:ncol(OrderAnimalia_all)){relative_abundance_all[,i] <- (OrderAnimalia_all[,i]/colSums(OrderAnimalia_all[1:ncol(OrderAnimalia_all)])[i])*100}
relative_abundance_all <- data.frame(matrix(nrow=nrow(OrderAnimalia_all_top50), ncol=ncol(OrderAnimalia_all_top50)
                                            , dimnames = list(c(rownames(OrderAnimalia_all_top50)), c(colnames(OrderAnimalia_all_top50)))))
for(i in 1:ncol(OrderAnimalia_all_top50)){relative_abundance_all[,i] <- (OrderAnimalia_all_top50[,i]/colSums(OrderAnimalia_all_top50[1:ncol(OrderAnimalia_all_top50)])[i])*100}

relative_abundance_all <- tibble::rownames_to_column(relative_abundance_all, "Species")
colnames(relative_abundance_all) <- gsub("X","",colnames(relative_abundance_all))
relative_species_alllocation <- melt(relative_abundance_all, na.rm = FALSE, value.name = "value", id='Species')
colnames(relative_species_alllocation) <- c("Species", "Location", "Relative_Abundance")
relative_species_alllocation[relative_species_alllocation == 0] <- NA

#Create heatmap rarefied data (all samples)
ggplot(data=relative_species_alllocation
       , aes(y=factor(Species, levels = rownames(OrderAnimalia_all_top50)),x=Location, fill= Relative_Abundance)
)+ 
  geom_tile(
  )+
  scale_fill_gradientn(colours=c( # For all the bad values <5%
    "darkgreen","green","yellow") # For the passed values
    , limits=c(0,100) # sets abundance range from 0 to 100 to get consitency when comparing samples
    # # , values=rescale(c(0,(5-(10^-10)) # For all the bad values <5%
    # # ,5,50,100)) # For the passed values
    , na.value="lightgrey"
    , breaks=c(25,50,75,100)
    , xlab("% \nRead\nAbundance\n")
  )+
  theme(legend.position = "right"
        , legend.title=element_text(size=10) #title text size
        , axis.text.y = element_text( face="italic",size=8, colour = "black") #species size
        , plot.title = element_text(size=15) #size of the title of the plot
        , legend.text = element_text(size = 10)
        , strip.text.x = element_text(size = 15)
        , axis.text.x = element_text(face="italic",size=8, colour = c(env_all$Zone_color), angle=90)
  )+
  guides(pattern = guide_legend(override.aes = list(fill = "green"),
                                title="")
  )+ 
  labs(title = "Relative read abundance % of each species per sample ", y=""
  )

heatmap_Animalia <- ggplot(data=relative_species_alllocation
                       , aes(y=factor(Species, levels = rownames(OrderAnimalia_all_top50)),x=Location, fill= Relative_Abundance)
)+ 
  geom_tile(
  )+
  scale_fill_gradientn(colours=c( # For all the bad values <5%
    "darkgreen","green","yellow") # For the passed values
    , limits=c(0,100) # sets abundance range from 0 to 100 to get consitency when comparing samples
    # # , values=rescale(c(0,(5-(10^-10)) # For all the bad values <5%
    # # ,5,50,100)) # For the passed values
    , na.value="lightgrey"
    , breaks=c(25,50,75,100)
    , xlab("% \nRead\nAbundance\n")
  )+
  theme(legend.position = "left"
        , legend.title=element_text(size=10) #title text size
        , axis.text.y = element_blank() #species size
        , axis.ticks.y = element_blank()
        , plot.title = element_text(size=15) #size of the title of the plot
        , legend.text = element_text(size = 10)
        , strip.text.x = element_text(size = 15)
        , axis.text.x = element_text(face="italic",size=8, colour = c(env_all$Zone_color), angle=90)
  )+
  guides(pattern = guide_legend(override.aes = list(fill = "green"),
                                title="")
  )+ 
  labs(title = "Relative read abundance % of each species per sample ", y=""
  )

grid.newpage()
print(heatmap_Animalia, 
      vp = viewport(x = 0.35, y = 0.5, width = 0.7, height = 1.0))
print(dendro_plot_top50,
      vp = viewport(x = 0.85, y = 0.545, width = 0.3, height = 0.90))



##Chordata

#Location dendro based on transformed data
transformed_location <- decostand(decostand(merged_data_rarefied_Chordata, method="tot"), method="max")
Location_dendro <- as.dendrogram(hclust(d=dist(x = (t(transformed_location))), method = "ward.D"))
Location_order <- order.dendrogram(Location_dendro)
lable_color <- env_rarefied[Location_order,]
labels_colors(Location_dendro) <- c(lable_color$Zone_color)
dendro_plot <- ggdendrogram(data = Location_dendro, rotate = F)
par(mar = c(8, 3, 2, 2))
plot(Location_dendro)

#Animalia dendro
merged_data_rarefied_Chordata_transformed <- decostand(decostand(merged_data_rarefied_Chordata, method="tot"), method="max")
Chordata_dendro <- as.dendrogram(hclust(d=dist(x = merged_data_rarefied_Chordata_transformed),method = "ward.D"))
dendro_plot_Chordata <- ggdendrogram(data = Chordata_dendro, rotate = T)
print(dendro_plot_Chordata)

#Order Animalia according to dendrogram
Chordata_order <- order.dendrogram(Chordata_dendro)
OrderChordata <- merged_data_rarefied_Chordata[Chordata_order,]

#Calculate the relative abundance of each species per sample
relative_abundance_all <- data.frame(matrix(nrow=nrow(OrderChordata), ncol=ncol(OrderChordata)
                                            , dimnames = list(c(rownames(OrderChordata)), c(colnames(OrderChordata)))))
for(i in 1:ncol(OrderChordata)){relative_abundance_all[,i] <- (OrderChordata[,i]/colSums(OrderChordata[1:ncol(OrderChordata)])[i])*100}

relative_abundance_all <- tibble::rownames_to_column(relative_abundance_all, "Species")
colnames(relative_abundance_all) <- gsub("X","",colnames(relative_abundance_all))
relative_species_alllocation <- melt(relative_abundance_all, na.rm = FALSE, value.name = "value", id='Species')
colnames(relative_species_alllocation) <- c("Species", "Location", "Relative_Abundance")
relative_species_alllocation[relative_species_alllocation == 0] <- NA

#Create heatmap rarefied data (all samples)
ggplot(data=relative_species_alllocation
       , aes(y=factor(Species, levels = rownames(OrderChordata)),x=Location, fill= Relative_Abundance)
)+ 
  geom_tile(
  )+
  scale_fill_gradientn(colours=c( # For all the bad values <5%
    "darkgreen","green","yellow") # For the passed values
    , limits=c(0,100) # sets abundance range from 0 to 100 to get consitency when comparing samples
    # # , values=rescale(c(0,(5-(10^-10)) # For all the bad values <5%
    # # ,5,50,100)) # For the passed values
    , na.value="lightgrey"
    , breaks=c(25,50,75,100)
    , xlab("% \nRead\nAbundance\n")
  )+
  theme(legend.position = "right"
        , legend.title=element_text(size=10) #title text size
        , axis.text.y = element_text( face="italic",size=8, colour = "black") #species size
        , plot.title = element_text(size=15) #size of the title of the plot
        , legend.text = element_text(size = 10)
        , strip.text.x = element_text(size = 15)
        , axis.text.x = element_text(face="italic",size=8, colour = c(env_all$Zone_color), angle=90)
  )+
  guides(pattern = guide_legend(override.aes = list(fill = "green"),
                                title="")
  )+ 
  labs(title = "Relative read abundance % of each species per sample ", y=""
  )
