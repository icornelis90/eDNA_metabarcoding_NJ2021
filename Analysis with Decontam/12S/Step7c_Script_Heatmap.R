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
proj.path <- here("/home/genomics/icornelis/02_ZEROimpact/01_12S/NJ2021/MiFish-UE_run2")

#upload data
table_rarefied <- readxl::read_excel(paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/table_unrarefied_concatenated_CleanedASVs_WithField_FullTaxonomicAssignment.xlsx"))
env <- read.csv(paste0(proj.path,"/Step5_Statistics/environmental_data.csv"),  header=TRUE, sep=";")

#add color and pch to environmental data to create the plot 
env$Environment_color <- ifelse(env$Zones=="Coastal","limegreen", 
                                ifelse(env$Zones=="zone1", "slateblue", 
                                       ifelse(env$Zones=="zone2","darkorange","red")))

#remove samples from environmental data that are absent from the rarefied data
keep_samples <- c(colnames(table_rarefied[,1:(ncol(table_rarefied)-11)]))
env_rarefied <- env[env$Niskin.sample %in% keep_samples,]
colorder_rarefied <- c(env_rarefied$Niskin.sample)
env_rarefied_samples <- env_rarefied %>% filter(!grepl("neg", env_rarefied$Niskin.sample))
colorder_samples <- c(env_rarefied_samples$Niskin.sample)

#select ASVs that are assigned to Fish
fish_classes <- readRDS(file = paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_classes.rds"))
freshwater_fish <- readRDS(file = paste0(proj.path,"/MiFish_UE-S_concatenated/results_v2/REnvironment/Fish_Freshwater.rds"))
table_rarefied_FishASVs <- as.data.frame(table_rarefied[table_rarefied$Class %in% c(fish_classes),])
table_rarefied_FishASVs <- as.data.frame(table_rarefied_FishASVs[!table_rarefied_FishASVs$Species %in% c(freshwater_fish, "NA"),])
rownames(table_rarefied_FishASVs) <- table_rarefied_FishASVs$ASV

#3: sort by species
taxo <- "Species"
merged_data_rarefied_fish<-aggregate(table_rarefied_FishASVs[,1:(ncol(table_rarefied_FishASVs)-11)], by= list(as.factor(table_rarefied_FishASVs[,taxo])),FUN=sum)
rownames(merged_data_rarefied_fish)<-merged_data_rarefied_fish$Group.1
merged_data_rarefied_fish$Group.1<-NULL
unfound_species <- c(rownames(merged_data_rarefied_fish[rowSums(merged_data_rarefied_fish)==0,]))
remove_samples <- c(colnames(merged_data_rarefied_fish[,colSums(merged_data_rarefied_fish)==0]))
merged_data_rarefied_fish <- merged_data_rarefied_fish[!rownames(merged_data_rarefied_fish) %in% unfound_species,]
#merged_data_rarefied_fish <- merged_data_rarefied_fish[,!colnames(merged_data_rarefied_fish) %in% remove_samples]
merged_data_rarefied_fish <- merged_data_rarefied_fish[,colorder_samples]
#b <- ifelse(rev(unique(rownames(merged_data_rarefied_fish))) %in% unfound_species,"grey","black")

#Location dendro based on transformed data
merged_data_rarefied_fish <- merged_data_rarefied_fish[,colorder_samples]
transformed_location <- decostand(decostand(merged_data_rarefied_fish, method="tot"), method="max")
Location_dendro <- as.dendrogram(hclust(d=dist(x = (t(transformed_location))), method = "ward.D"))
Location_order <- order.dendrogram(Location_dendro)
label_color <- env_rarefied[Location_order,]
labels_colors(Location_dendro) <- c(label_color$Zone_color)
plot(0,type='n',axes=FALSE,ann=FALSE)
dendro_plot <- ggdendrogram(data = Location_dendro, rotate = F)
par(mar = c(8, 3, 2, 2))
plot(Location_dendro)

#Fish dendro based on transformed data
transformed_fish <- decostand(merged_data_rarefied_fish, method="log")
Fish_dendro <- as.dendrogram(hclust(d=dist(x = transformed_fish), method = "ward.D"))
dendro_plot <- ggdendrogram(data = Fish_dendro, rotate = T) + 
  theme(axis.text.y = element_text(size = 8))
print(dendro_plot)

#Order Fish Species according to dendrogram
Fish_order <- order.dendrogram(Fish_dendro)
OrderFish <- merged_data_rarefied_fish[Fish_order, ]
env_order <- env_rarefied[Location_order,]

#Calculate the relative abundance of each fish species per sample
relative_abundance <- data.frame(matrix(nrow=nrow(OrderFish), ncol=ncol(OrderFish)
                                       , dimnames = list(c(rownames(OrderFish)), c(colnames(OrderFish)))))
for(i in 1:ncol(OrderFish)){relative_abundance[,i] <- (OrderFish[,i]/colSums(OrderFish[1:ncol(OrderFish)])[i])*100}
relative_abundance <- tibble::rownames_to_column(relative_abundance, "Species")
colnames(relative_abundance) <- gsub("X","",colnames(relative_abundance))
relative_species_location <- melt(relative_abundance, na.rm = FALSE, value.name = "value", id='Species')
colnames(relative_species_location) <- c("Species", "Location", "Relative_Abundance")
relative_species_location[relative_species_location == 0] <- NA

#Create heatmap rarefied data (samples < 10000 reads removed)
ggplot(data=relative_species_location
       , aes(y=factor(Species, levels = rownames(OrderFish)),x=Location, fill= Relative_Abundance)
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
        , axis.text.y = element_text(face="italic",size=8, colour = "black") #species size
        , plot.title = element_text(size=15) #size of the title of the plot
        , legend.text = element_text(size = 10)
        , strip.text.x = element_text(size = 15)
        , axis.text.x = element_text(face="italic",size=8, colour = c(env_rarefied$Zone_color), angle=90)
  )+
  guides(pattern = guide_legend(override.aes = list(fill = "green"),
                                title="")
  )+ 
  labs(title = "Relative read abundance % of each species per sample ", y=""
  )
