#0.1 Set working directory
#setwd ("~/working_directory/")
setwd("C:/Users/gaam1/OneDrive/Documents/AZ20 Inhibitor/MouseMeiosis/RAD1_CKO_Scripts")
#0.2 libraries that need to be loaded
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggrepel) #label peptides
library(pheatmap) #make heatmap
library(ggExtra) #make ggplot and seqlogo grids
library(RColorBrewer)
library(tidyverse)

#########################################################################################
#                            Filter Phosphoproteomics Data                              #
#                           Localization Probability > 0.85                             #
#                    At least in 2 experiments for RAD1 CKO & AZ20                      #
#########################################################################################


pp <- 0 
cp <- 0
regulation <- 0.5

df.count <- data.frame(pp=0, cp = 0, regulation = 0,
                       my_factor=0, NAZ20 = 0, 
                       NRAD1= 0, 
                       BOTH=0, CT=0, NR= 0, 
                       Q1= 0, Q3= 0, 
                       ALL=0,
                       STQ_ALL=0,
                       Q2=0,
                       Q4=0,
                       az20=0,
                       rad1=0,
                       STQ_Q2=0,
                       STQ_Q4=0,
                       STQ_az20=0,
                       STQ_rad1=0
)


Loc.Prob <- c(0.85) # ENTER HERE, SEPARATED BY COMMA, ALL THE Localization.Probs TO BE TESTED
N.exp <- c(2) # ENTER HERE, SEPARATED BY COMMA, THE CUTTOF FOR # OF OBSERVATION IN EXPERIMENTS
Min.Ratio <- c(1) # ENTER HERE, SEPARATED BY COMMA, THE CUTOFF FOR REGULATION

for (v in 1:length(Loc.Prob)) {
  for (z in 1:length(N.exp)) {
    for (w in 1:length(Min.Ratio)) {
      
      experiment1 <- "MouseMeiosis.txt"  
      
      pp= Loc.Prob[v]
      cp= N.exp[z]
      regulation <- Min.Ratio[w]
      my_factor <- 4 # Bow-Tie "angle"
      
      Data <- read.delim (experiment1, sep="\t", header =TRUE, fill = TRUE) %>%
        filter (Localization.prob > Loc.Prob[v])%>%
        filter (AZ20 >= N.exp[z])%>%
        filter (RAD1 >= N.exp[z])%>%
        mutate (stq_sites = ifelse(str_detect(substr(Seq_window,8,8), "Q"), "STYQ", "STYX" ))%>%
        ungroup()
      
      #ADDING GEOMETRIC MEAN
      Mean_AZ20 <- 0
      Mean_RAD1 <- 0
      for(i in 1:length(Data$UNIPROT)){
        Mean_AZ20[i] <- mean(c(Data$AZ20_EXP1[i],Data$AZ20_EXP2[i], Data$AZ20_EXP3[i], Data$AZ20_EXP4[i]), na.rm = TRUE)
        Mean_RAD1[i] <- mean(c(Data$RAD1_EXP1[i],Data$RAD1_EXP2[i], Data$RAD1_EXP3[i]), na.rm = TRUE)
      }
      
      Data <- cbind(Data,Mean_AZ20)
      Data <- cbind(Data,Mean_RAD1)
      
      # 3.1 filter for matching data in both datasets
      DataComb.Filter <- Data%>%
        filter(Mean_AZ20 != "NaN" & Mean_RAD1 != "NaN"  )%>%
        #mutate (clust = ifelse(str_detect(cluster_pep, "\\("), "black", "red" ))%>%
        ungroup()  
      
      #filter for non-regulated proteins
      DataComb.center <- DataComb.Filter %>%
        filter((Mean_AZ20)^2 + (Mean_RAD1)^2 <= regulation) %>%
        mutate (quadr = "CENTER") %>%
        ungroup ()
      
      #filter for regulated proteins (BOW-TIE)
      DataComb.Q1 <- DataComb.Filter%>%
        filter((Mean_AZ20)^2 + (Mean_RAD1)^2 > regulation)%>% #treshold for exclusion of non-regulated proteins = 1
        filter (Mean_RAD1 > 0)%>% # remove Q3/Q4
        filter (Mean_AZ20 <= -1*((Mean_RAD1/my_factor)) & Mean_RAD1 >= abs(Mean_AZ20/my_factor)) %>%
        mutate (quadr = "Q1") %>%
        ungroup ()
     
      DataComb.Q2 <- DataComb.Filter%>%
        filter((Mean_AZ20)^2 + (Mean_RAD1)^2 > regulation)%>% #treshold for exclusion of non-regulated proteins = 1
        filter (Mean_AZ20 >0)%>% # remove Q1/Q4
        filter(Mean_AZ20 >= (Mean_RAD1/my_factor) & Mean_RAD1 >= (Mean_AZ20/my_factor)) %>%
        mutate (quadr = "Q2") %>%
        ungroup ()
      
      DataComb.Q3 <- DataComb.Filter%>%
        filter((Mean_AZ20)^2 + (Mean_RAD1)^2 > regulation)%>% #treshold for exclusion of non-regulated proteins = 1
        filter (Mean_RAD1 <0)%>% # remove Q1/Q2
        filter((Mean_AZ20 <= -1*((Mean_RAD1/my_factor)) & Mean_RAD1 >= abs(Mean_AZ20/my_factor)) | 
                 (Mean_AZ20 >= abs(Mean_RAD1/my_factor) & Mean_RAD1 <= -1*(Mean_AZ20/my_factor))) %>%
        mutate (quadr = "Q3") %>%
        ungroup ()
      
      DataComb.Q4 <- DataComb.Filter%>%
        filter((Mean_AZ20)^2 + (Mean_RAD1)^2 > regulation)%>% #treshold for exclusion of non-regulated proteins = 1
        filter (Mean_AZ20 <0)%>% # remove Q2/Q3
        filter(Mean_AZ20 <= (Mean_RAD1/my_factor) & Mean_RAD1 <= (Mean_AZ20/my_factor)) %>%
        mutate (quadr = "Q4")%>%
        ungroup ()
      
      DataComb.NotReg <- anti_join(DataComb.Filter, DataComb.center, by = "Phosphosite");
      DataComb.NotReg <- anti_join(DataComb.NotReg, DataComb.Q1, by = "Phosphosite");
      DataComb.NotReg <- anti_join(DataComb.NotReg, DataComb.Q2, by = "Phosphosite");
      DataComb.NotReg <- anti_join(DataComb.NotReg, DataComb.Q3, by = "Phosphosite");
      DataComb.NotReg <- anti_join(DataComb.NotReg, DataComb.Q4, by = "Phosphosite");
      
      DataComb.NotReg <- DataComb.NotReg%>%
        mutate (quadr = "NR")%>%
        ungroup ()
      
      DataBT <-  rbind(DataComb.Q1, DataComb.Q2, DataComb.Q3, DataComb.Q4, DataComb.NotReg, DataComb.center)

      DataComb.STQ.ALL <- DataBT %>%
        filter (stq_sites == "STYQ")%>%
        ungroup()
      DataComb.STQ.Q4 <- DataComb.Q4 %>%
        filter (stq_sites == "STYQ")%>%
        ungroup()
      #####az20#####
      DataComb.az20 <- DataComb.Filter%>%
        filter((Mean_AZ20)^2 + (Mean_RAD1)^2 > regulation)%>% #treshold for exclusion of non-regulated proteins = 1
        filter (Mean_AZ20 < 0 )%>%
        ungroup ()
      DataComb.STQ.az20 <- DataComb.az20 %>%
        filter (stq_sites == "STYQ")%>%
        ungroup()
      #####rad1#####
      DataComb.rad1 <- DataComb.Filter %>%
        filter((Mean_AZ20)^2 + (Mean_RAD1)^2 > regulation)%>% #treshold for exclusion of non-regulated proteins = 1
        filter (Mean_RAD1 <0)%>% 
        ungroup ()
      DataComb.STQ.rad1 <- DataComb.rad1 %>%
        filter (stq_sites == "STYQ")%>%
        ungroup()
      
      color_center <- "black"#333333"
      color_not_reg <- "#dbdbdb"
      color_13 <- "#28ACE3"   #"#4a9dc9"#"goldenrod1"
      color_12 <- "#0B0383"#"black" 
      point_size <- 2
      p <- ""
      p <- ggplot(DataComb.Filter, aes(y=(-Mean_RAD1), x=(-Mean_AZ20))) 
      p <- p + geom_point(data = DataComb.center, size=(point_size-0.5),  alpha=0.2, color=color_center)
      p <- p + expand_limits(x = 0, y = 0) +
        scale_y_continuous(expand = c(0,0), breaks=seq(0,5.5,1), minor_breaks = NULL, limits = c(0, 5.5)) +
        scale_x_continuous(expand = c(0,0), breaks=seq(0,5.5,1), minor_breaks = NULL, limits = c(0, 5.5)) +
        ylab(expression(Log[2](WT/RAD1~CKO))) +
        xlab(expression(Log[2](Vehicle/ATRi))) +
        ggtitle(paste("Loc.Prob >", pp, "; ", "N.exp =>", cp, "; ", "Min.Ratio >", regulation )) 
      #p <- p  
      p <- p + geom_point(data= DataComb.NotReg,  size=point_size, alpha=0.3, color=color_not_reg)
      p <- p + geom_point(data= DataComb.Q2,  size=point_size, alpha=0.6, color=color_not_reg) 
      p <- p + geom_point(data= DataComb.Q1,  size=point_size, alpha=0.3, color=color_not_reg)
      p <- p + geom_point(shape=21, data= DataComb.Q4,  size=point_size, fill=color_13, colour = color_12)
      p <- p + geom_point(data= DataComb.Q3,  size=point_size, alpha=0.3, color=color_not_reg)
      p.sqtq <- p + geom_point(shape = 21, data= DataComb.STQ.Q4, fill = "#ED2024",
                               colour = "#540500", size=2) + 
        theme_classic()+ theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                               axis.title.x = element_text(size=14, face="bold"),
                               axis.title.y = element_text(size=14, face="bold"), 
                               axis.text.x = element_text(face="bold", size=12, color = "black"),
                               axis.text.y = element_text(face="bold", size=12, color = "black"), 
                               axis.line = element_line(size = 1), 
        ) +
      geom_text_repel(aes(label = Phosphosite),
                      data = subset(DataBT, Phosphosite == "ATR_927" |
                                      Phosphosite == "ATR_440" |
                                      Phosphosite == "SMC3_787" |
                                      Phosphosite == "SMC1B_370" |
                                      Phosphosite == "MDC1_157" | 
                                      Phosphosite == "SYCP3_218" | 
                                      Phosphosite == "SYCP2_646" | 
                                      Phosphosite == "RAD50_635" | 
                                      Phosphosite == "CTIP_588" | 
                                      Phosphosite == "UIMC1_171" | 
                                      Phosphosite == "NBN_398" | 
                                      Phosphosite == "RAD9B_326" |
                                      Phosphosite == "RAD9B_361" | 
                                      Phosphosite == "MCM3_738" | 
                                      Phosphosite == "NEK1_750" | 
                                      Phosphosite == "TOPB1_1380" | 
                                      Phosphosite == "NIPBL_889"),
                      nudge_x = 5,
                      segment.size  = 0.2,
                      direction     = "y", 
                      hjust = 0)
      ggsave(paste(experiment1,pp,cp,regulation,"STQ_Q4_b", ".pdf", sep = "_"), height = 5, width = 5, dpi = 500)
      df.count[nrow(df.count) + 1,] = c(pp,
                                        cp, 
                                        regulation, my_factor, 
                                        NAZ20 = length(which(Data$AZ20 > 0)), 
                                        NRAD1= length(which(Data$RAD1 > 0)),
                                        BOTH=count(DataComb.Filter), 
                                        CT= count(DataComb.center), 
                                        NR= count(DataComb.NotReg), 
                                        Q1= count(DataComb.Q1), 
                                        Q3= count(DataComb.Q3), 
                                        ALL= count(Data), 
                                        STQ_ALL= count(DataComb.STQ.ALL),
                                        Q2= count(DataComb.Q2),
                                        Q4= count(DataComb.Q4),
                                        az20= count(DataComb.az20),
                                        rad1= count(DataComb.rad1),
                                        STQ_Q4= count(DataComb.STQ.Q4),
                                        STQ_az20= count(DataComb.STQ.az20),
                                        STQ_rad1= count(DataComb.STQ.rad1)
      )
      #Attach_protein_description
      write.table (DataBT, file= paste(experiment1,pp,cp,regulation,"ALL.xls", sep = "_"), sep="\t", row.names = FALSE)
      write.table (DataComb.Q4, file= paste(experiment1,pp,cp,regulation,"Q4.xls", sep = "_"), sep="\t", row.names = FALSE)
      write.table (DataComb.STQ.ALL, file= paste(experiment1,pp,cp,regulation,"STQ_ALL.xls", sep = "_"), sep="\t", row.names = FALSE)
      write.table (DataComb.STQ.Q4, file= paste(experiment1,pp,cp,regulation,"STQ_Q4.xls", sep = "_"), sep="\t", row.names = FALSE)
    }
  }
}

#########################################################################################
#                          Generate Amino Acid Frequency Function                       #
#########################################################################################

amino_acid_freq <- function(random){ 
  freq_data <- data.frame(random$Seq_window)
  colnames( freq_data) <- "Seq_window"
  data_frequency <-  freq_data %>% separate(Seq_window, sep = c(1,2,3,4,5,6,7,8,9,10,11,12,13),
                                            c("-6", "-5","-4", "-3", "-2","-1","P", "1", "2", "3", "4", "5", "6" )) %>%
    ungroup()
  `P-6` <- data.frame(table(data_frequency$`-6`))
  `P-5` <- data.frame(table(data_frequency$`-5`))
  `P-4` <- data.frame(table(data_frequency$`-4`))
  `P-3` <- data.frame(table(data_frequency$`-3`))
  `P-2` <- data.frame(table(data_frequency$`-2`))
  `P-1` <- data.frame(table(data_frequency$`-1`))
  P <- data.frame(table(data_frequency$P))
  P1 <- data.frame(table(data_frequency$`1`))
  P2 <- data.frame(table(data_frequency$`2`))
  P3 <- data.frame(table(data_frequency$`3`))
  P4 <- data.frame(table(data_frequency$`4`))
  P5 <- data.frame(table(data_frequency$`5`))
  P6 <- data.frame(table(data_frequency$`6`))
  
  r11 <- left_join(`P-6`, `P-5`, by = "Var1")
  r11 <- left_join(r11, `P-4`, by = "Var1")
  r11 <- left_join(r11, `P-3`, by = "Var1")
  r11 <- left_join(r11, `P-2`, by = "Var1")
  r11 <- left_join(r11, `P-1`, by = "Var1")
  r11 <- left_join(r11, P, by = "Var1")
  r11 <- left_join(r11, P1, by = "Var1")
  r11 <- left_join(r11, P2, by = "Var1")
  r11 <- left_join(r11, P3, by = "Var1")
  r11 <- left_join(r11, P4, by = "Var1")
  r11 <- left_join(r11, P5, by = "Var1")
  r11 <- left_join(r11, P6, by = "Var1")
  colnames(r11) <- c("Amino Acid",
                     "P-6", 
                     "P-5",
                     "P-4",
                     "P-3",
                     "P-2",
                     "P-1",
                     "P",
                     "P1",
                     "P2",
                     "P3",
                     "P4",
                     "P5",
                     "P6")
  r11 <- r11[-1,]
  freq_table <- as.data.frame(matrix(c("A","C",'D',"E",
                                       "F","G","H","I",
                                       "K","L","M","N",
                                       "P","Q","R","S",
                                       "T","V","W","Y"), ncol = 1))
  colnames(freq_table) <- c("Amino Acid")
  freq_table_data <- left_join(freq_table, r11)
  rownames(freq_table_data) <- freq_table$`Amino Acid`
  freq_table_data <- freq_table_data[,-1]
}

#########################################################################################
#                          Compare Amino Acid Frequency Function                        #
#########################################################################################
quadrant_analysis <- function(table_1, table_2){
  table_1_p <- data.frame(1:20)
  table_2_p <- data.frame(1:20)
  for(i in 1:length(table_1)){
    data <- table_1
    table_1_pro <- data[[i]]/sum(data[[i]])
    table_1_p[i] <- table_1_pro
  }
  for(i in 1:length(table_2)){
    data_2 <- table_2
    table_2_pro <- data_2[[i]]/sum(data_2[[i]])
    table_2_p[i] <- table_2_pro
  }
  am_comparison <- data.frame(log2(table_1_p/table_2_p))
}

#########################################################################################
#                         Plot Acid Frequency Comparison Function                       #
#########################################################################################
heatmap_am <- function(am_quadr, number_display,comparison, file_name){
  #jpeg(paste(experiment1,pp,cp,regulation,"AminoAcid_fdr",".jpeg", sep = "_"),width = 6, height = 6,units = "in", res = 300)
  breaks1 <- c(-2.5,-2,-1.5,-1,-0.5, 0.5, 1, 1.5, 2, 2.5) #seq(from= -1.5, to = 2, by = 0.5)
  pheatmap(am_quadr, 
           cluster_rows = FALSE, 
           cluster_cols = FALSE, 
           display_numbers = number_display,
           breaks = breaks1,
           labels_col = c("P-6","P-5","P-4","P-3","P-2","P-1","P","P1","P2","P3","P4","P5","P6"),
           labels_row =c("A","C",'D',"E",
                         "F","G","H","I",
                         "K","L","M","N",
                         "P","Q","R","S",
                         "T","V","W","Y"),
           legend_breaks =  c(-2.75,-2.25,-1.75,-1.25,-0.75, 0.75, 1.25, 1.75, 2.25, 2.75),
           legend_labels = breaks1, 
           color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(9), 
           number_color = "black",
           border_color = "black",
           angle_col = "0",
           main = paste("Loc.Prob >", pp, "; ", "N.exp >=", cp, "; ", "Min.Ratio >", regulation, ";", comparison ),
           filename = file_name, 
           fontsize = 12
           
  )
  
}
#########################################################################################
#                                  Plot Q4 vs Center                                    #
#########################################################################################
am_q4 <- amino_acid_freq(DataComb.Q4)
am_q4[is.na(am_q4)] <- 0
am_center <- amino_acid_freq(DataComb.center)
am_center[is.na(am_center)] <- 0
q4_analysis <- quadrant_analysis(am_q4, am_center)
q4_analysis[q4_analysis==-Inf] <- NA

heatmap_am(q4_analysis,TRUE,"Q4 vs Center", "q4_center_norm.pdf")

