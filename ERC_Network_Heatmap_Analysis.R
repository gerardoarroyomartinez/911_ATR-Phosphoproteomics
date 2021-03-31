#0.1 Set working directory
setwd ("~/working_directory/")
setwd("C:/Users/gaam1/OneDrive/Documents/AZ20 Inhibitor/MouseMeiosis/species analysis")
#0.2 libraries that need to be loaded
library(qgraph) #Plot network
library(pheatmap) #make heatmap
library(dplyr)#data manipulation 
library(viridis)#change color in heatmap

################################################################################################
#                                Generate ERC Network Plot                                     #
################################################################################################

nw_erc_1 <- read.table("COMB-GO_0007127_meiosis_I.txt", sep = "\t", header = TRUE) #open ERC data
nw_erc <- nw_erc_1[2:6] #Select needed columns 
nw_erc[nw_erc == "ND"] <- 0 #Change missing values to 0
nw_erc[nw_erc == "N/A"] <- 0 #Change missing values to 0
nw_erc[nw_erc <= 0.4] <- 0 #Change  values below 0.3 ERC to 0
nw_erc <- mapply(nw_erc, FUN=as.numeric) #make columns to numeric 
nw_matrix <- t(matrix(unlist(nw_erc), ncol = 5, byrow = FALSE)) #transpose and make table into matrix
colnames(nw_matrix) <- nw_erc_1$X #rename column names into gene names 
row.names(nw_matrix) <- nw_erc_1[1:5,1] #rename rownames into 911 subunits 
empty_matrix <- matrix(0, nrow = 83, ncol = 88) #make empty matrix to have even columns and rows 
combined_matixes <- rbind(nw_matrix, empty_matrix) #combine data and empty matrix

subunits <- combined_matixes[1:38,1:5] #select 911 subunit columns 
meiosis_i <- combined_matixes[1:38,6:88] #select meiosis I genes 
meiosis_i_cero <- (colSums(meiosis_i, na.rm=T) != 0) #select columns with at least one ERC value for any of the subunits
meiosis_i_filtered <- meiosis_i[,meiosis_i_cero] # remove all the non-zero columns
comb <- cbind(subunits, meiosis_i_filtered) #combines filtered gene table and subunit table

groups <- list("RAD1"=1, "RAD9A"=2, "RAD1B"=3, "HUS1"=4, "HUS1B"=5) #make 911 subunit groups
#Plot network
qgraph(comb, shape="circle", layout="spring", vsize=5, threshold=0.3, theme="classic",
       label.scale.equal=TRUE,label.scale=FALSE, maximum=1, groups=groups,
       color = c("#a6cfcd", "#bc84b5", "#f1ee95", "#bdd697", "#f9a45e"), 
       border.color="black", 
       esize=FALSE,
       edge.color="black",
       arrows=FALSE, 
       bg="#bcbec0",
       legend=FALSE)

################################################################################################
#                                Generate ERC Heatmap Plot                                     #
################################################################################################

nw_erc_1 <- read.table("COMB-GO_0007127_meiosis_I.txt", sep = "\t", header = TRUE) #open ERC data 
nw_erc <- nw_erc_1[2:6] #select columns 
nw_erc[nw_erc == "ND"] <- 0 #Change missing values to 0
nw_erc[nw_erc == "N/A"] <- 0 #Change missing values to 0
nw_erc[nw_erc <= 0.4] <- 0 #Change 0.4 values to 0
nw_erc <- mapply(nw_erc, FUN=as.numeric) #make columns numeric 
nw_matrix <- (matrix(unlist(nw_erc), ncol = 5, byrow = FALSE)) #change table to matrix
colnames(nw_matrix) <- nw_erc_1[1:5,1] #rename columns 
row.names(nw_matrix) <- nw_erc_1$X #rename rows 

erc_meiosis <- nw_matrix #rename table
erc_meiosis_filter <- (rowSums(erc_meiosis, na.rm=T) != 0) #select rows with at least one ERC value for any of the subunits
erc_meiosis_filter_no_0 <- as.data.frame(erc_meiosis[erc_meiosis_filter,]) #make data frame
erc_meiosis_filter_no_0 <- cbind(X=rownames(erc_meiosis_filter_no_0), erc_meiosis_filter_no_0 ) #add column with gene names
#use gene list to clean ERC data sets to only plot having at least one ERC value >= 0.4 for any of the 911 subunit
filtered_table <- semi_join(nw_erc_1,erc_meiosis_filter_no_0, by=c("X")) 
filtered_table[filtered_table == "ND"] <- NA #Change ND values to NA
filtered_table[filtered_table == "N/A"] <- NA #Change N/A values to NA
filtered_table$X -> filtered_table_rownames #select row names 
filtered_table <- mapply(filtered_table[,2:6], FUN=as.numeric) #change columns into numeric
rownames(filtered_table) <- filtered_table_rownames #change row names

#plot heatmap
pheatmap(filtered_table[-19,], 
         scale = "none", 
         breaks = c(0,0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1, 1.1),
         legend_breaks = c(0,0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1),
         legend_labels = c(0,0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1),
         color = c("#e5e5e5", rev(magma(n=9))),
         cluster_cols  = FALSE,
         cluster_rows = TRUE,
         angle_col = "45",
         display_numbers = TRUE,
         number_color = "black",
         border_color = "black",
         cutree_rows=2, 
         cellwidth = 23, 
         cellheight = 13, 
         clustering_distance_rows = "euclidean",
         na_col = "white",
         filename = "hm_meiosis I.pdf"
)
dev.off()
