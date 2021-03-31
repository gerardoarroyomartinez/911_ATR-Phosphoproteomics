#Set working directory
setwd ("~/working_directory/")
setwd("C:/Users/gaam1/OneDrive/Documents/AZ20 Inhibitor/MouseMeiosis")
#0.2 libraries that need to be loaded
library(GOplot) #Plot GO Chord plot 


################################################################################################
#                                 Generate GO Chord Plot                                       #
################################################################################################

chr1 <- read.delim ("enrichment.Process_STQ.tsv", sep="\t", header =TRUE, fill = TRUE) #open GO analysis results
chr2 <- read.delim ("Chord_STQ_Q4.txt", sep="\t", header =TRUE, fill = TRUE) #open gene list used to obtain GO terms
plot_chord <- circle_dat(chr1, chr2) # combine GO analysis and gene list
chord <- chord_dat(plot_chord, chr2, chr1[1:10,3]) #select top 10 GO terms
chord[is.na(chord)] <- 0 #change missing values to 0 
pdf(file="chord2_GO.pdf", height = 15 , width = 13.5) # save plot to PDF
GOChord(chord, gene.space = 0.25, gene.size = 8, border.size = .5, space = 0.05, lfc.col = c("#f7392d", "white")) #plot GO Chord 
dev.off()

