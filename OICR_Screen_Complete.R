# OICR Membrane Screen Code
# This script takes the fluorescence data from DiOC2(3) and TO-PRO-3 Iodide dyes and plots the values relative to the untreated controls
# Author: Scott McAuley
# Last Updated: 2016-10-21

# Version: 1.1
# Edits: Changes function from running any calculation to just importing the data

# ================== SOURCES ================== 
source("http://bioconductor.org/biocLite.R") # Sources the biocLite.R installation script. 

# ================== PACKAGES ================== 
if (!require(ggplot2, quietly=TRUE)) {
  install.packages("ggplot2")
  library(ggplot2)
}

if (!require(compare, quietly=TRUE)) {
  biocLite("ChemmineR")
  library("ChemmineR")
}

if (!require(ggplot2, quietly=TRUE)) {
  biocLite("ChemmineOB")
  library("ChemmineOB") # Loads the package
}

if (!require(ggplot2, quietly=TRUE)) {
  install.packages("xlsx")
  library(xlsx) # Loads the package
}

# ================== GLOBAL CONSTANTS ================== 
PROJECTDIR <- "/Users/smcauley/Documents/Dropbox/nodwell_lab/Data/OICR screen/Membrane Pertubation/Membrane_Screen_Analysis"

# ================== FUNCTIONS ================== 

extract_data1 <- function(folder){
  plate1 <- c("A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10", "A11", "B02", "B03",
              "B04", "B05", "B06", "B07", "B08", "B09", "B10", "B11", "C02", "C03", "C04", 
              "Nisin", "CCCP")
  
  plate2 <- c("C05", "C06", "C07", "C08", "C09", "C10", "C11", "D02", "D03", "D04", "D05", "D06",
              "D07", "D08", "D09", "D10", "D11", "E02", "E03", "E04", "E05", "E06", "E07",
              "Nisin", "CCCP")
  
  plate3 <- c("E08", "E09", "E10", "E11", "F02", "F03", "F04", "F05", "F06", "F07", "F08", "F09", 
              "F10", "F11", "G02", "G03", "G04", "G05", "G06", "G07", "G08", "G09", "G10", 
              "Nisin", "CCCP")
  
  plate4 <- c("G11", "H02", "H03", "H04", "H05", "H06", "H07", "H07", "H09", "H10", "H11", "Nisin", 
              "CCCP", "N/A", "N/A", "N/A", "N/A","N/A", "N/A", "N/A", "N/A", "N/A", "N/A",  
              "N/A", "N/A")
  
  plateNames <- as.matrix(cbind(plate1, plate2, plate3, plate4))
  
  plates <- c("Plate1.txt", "Plate2.txt", "Plate3.txt", "Plate4.txt")
  
  full_data <- list()
  
  Final450 <- data.frame()
  Final600 <- data.frame()
  
  for(i in 1:4){
    # Read in data
    data <- read.table(plates[i], header = TRUE) 
    
    # Untreated Controls w/ culture
    culture_control600 <- vector()
    culture_control450 <- vector()
    UTControl600 <- as.numeric(data[93,3])
    UTControl450 <- as.numeric(data[93,2])
    
    # Untreated Controls w/o culture
    dye_control600 <- vector()
    dye_control450 <- vector()
    DControl600 <- as.numeric(data[94,3])
    DControl450 <- as.numeric(data[94,2])
    
    if (i != 4) {
      # Treated Samples at high concentration (with culture)
      highsamples450 <- as.numeric(data[c(seq(1,89,4),95,96),2])
      highsamples600 <- as.numeric(data[c(seq(1,89,4),95,96),3])
      
      # Treated Samples at low concentration (with culture)
      lowsamples450 <- as.numeric(data[c(seq(3,91,4),95,96),2])
      lowsamples600 <- as.numeric(data[c(seq(3,91,4),95,96),3])
      
      # Treated Controls at high concentration (no culture)
      highTControl450 <- as.numeric(data[c(seq(2,90,4),94,94),2])
      highTControl600 <- as.numeric(data[c(seq(2,90,4),94,94),3])
      
      # Treated Controls at low concentration (no culture)
      lowTControl450 <- as.numeric(data[c(seq(4,92,4),94,94),2])
      lowTControl600 <- as.numeric(data[c(seq(4,92,4),94,94),3])
      
      # Make Untreated culture the same length as samples
      culture_control600[1:25] <- UTControl600
      culture_control450[1:25] <- UTControl450
      
      # Make Untreated dye the same length as samples
      dye_control600[1:25] <- DControl600
      dye_control450[1:25] <- DControl450
      
    } else {
      # Treated Samples at high concentration (with culture)
      highsamples450 <- as.numeric(data[c(seq(1,41,4),95,96),2])
      highsamples600 <- as.numeric(data[c(seq(1,41,4),95,96),3])
      
      # Treated Samples at low concentration (with culture)
      lowsamples450 <- as.numeric(data[c(seq(3,43,4),95,96),2])
      lowsamples600 <- as.numeric(data[c(seq(3,43,4),95,96),3])
      
      # Treated Controls at high concentration (no culture)
      highTControl450 <- as.numeric(data[c(seq(2,42,4),94,94),2])
      highTControl600 <- as.numeric(data[c(seq(2,42,4),94,94),3])
      
      # Treated Controls at low concentration (no culture)
      lowTControl450 <- as.numeric(data[c(seq(4,44,4),94,94),2])
      lowTControl600 <- as.numeric(data[c(seq(4,44,4),94,94),3])
      
      # Make Untreated culture the same length as samples
      culture_control600[1:13] <- UTControl600
      culture_control450[1:13] <- UTControl450
      
      # Make Untreated dye the same length as samples
      dye_control600[1:13] <- DControl600
      dye_control450[1:13] <- DControl450
      
    }
    
    # Correction relative to untreated control
    # high concentration
    highcorrsample450 <- log2(highsamples450/UTControl450)
    highcorrsample600 <- log2(highsamples600/UTControl600)
    
    highcorrTControl450 <- log2(highTControl450/DControl450)
    highcorrTControl600 <- log2(highTControl600/DControl600)
    
    # low concentration
    lowcorrsample450 <- log2(lowsamples450/UTControl450)
    lowcorrsample600 <- log2(lowsamples600/UTControl600)
    
    lowcorrTControl450 <- log2(lowTControl450/DControl450)
    lowcorrTControl600 <- log2(lowTControl600/DControl600)
    
    # Correct for effect of compound on dye
    # high concentration
    Finalhighcorrsample450 <- highcorrsample450 - highcorrTControl450
    Finalhighcorrsample600 <- highcorrsample600 - highcorrTControl600
    
    Finallowcorrsample450 <- lowcorrsample450 - lowcorrTControl450
    Finallowcorrsample600 <- lowcorrsample600 - lowcorrTControl600
    
    # Create compound ID
    plate <- sapply(strsplit(as.character(folder), split = " "), "[", 2)
    
    if (i != 4){
      ID <- paste("OICR-CP", plate, "-", plateNames[1:25,i], sep = "")
      
      Finalcorrsample450 <- as.data.frame(cbind(ID, folder, plateNames[1:25,i], "DiOC2-3", highsamples450, lowsamples450, highTControl450, lowTControl450, culture_control450, dye_control450, Finalhighcorrsample450, Finallowcorrsample450))
      Finalcorrsample600 <- as.data.frame(cbind(ID, folder, plateNames[1:25,i], "TO-PRO-3", highsamples600, lowsamples600, highTControl600, lowTControl600, culture_control600, dye_control600, Finalhighcorrsample600, Finallowcorrsample600))
    } else {
      ID <- paste("OICR-CP", plate, "-", plateNames[1:13,i], sep = "")
      
      Finalcorrsample450 <- as.data.frame(cbind(ID, folder, plateNames[1:13,i], "DiOC2-3", highsamples450, lowsamples450, highTControl450, lowTControl450, culture_control450, dye_control450, Finalhighcorrsample450, Finallowcorrsample450))
      Finalcorrsample600 <- as.data.frame(cbind(ID, folder, plateNames[1:13,i], "TO-PRO-3", highsamples600, lowsamples600, highTControl600, lowTControl600, culture_control600, dye_control600, Finalhighcorrsample600, Finallowcorrsample600))
    }
    
    colnames(Finalcorrsample450) <- c("ID", "Plate", "Well", "Dye", "25_Sample", "5_Sample", "25_Control", "5_Control", "Untreated_culture", "Untreated_dye", "25_Ratio", "5_Ratio")
    colnames(Finalcorrsample600) <- c("ID", "Plate", "Well", "Dye", "25_Sample", "5_Sample", "25_Control", "5_Control", "Untreated_culture", "Untreated_dye", "25_Ratio", "5_Ratio")
    
    Final450 <- rbind(Final450, Finalcorrsample450)
    Final600 <- rbind(Final600, Finalcorrsample600)
    
    # Final <- list(Final450, Final600)
    Final <- rbind(Final450, Final600)
    
    # Forces values as numeric
    Final$`25_Sample` <- as.numeric(as.character(Final$`25_Sample`))
    Final$`5_Sample` <- as.numeric(as.character(Final$`5_Sample`))
    Final$`25_Control` <- as.numeric(as.character(Final$`25_Control`))
    Final$`5_Control` <- as.numeric(as.character(Final$`5_Control`))
    Final$`Untreated_culture` <- as.numeric(as.character(Final$`Untreated_culture`))
    Final$`Untreated_dye` <- as.numeric(as.character(Final$`Untreated_dye`))
    Final$`25_Ratio` <- as.numeric(as.character(Final$`25_Ratio`))
    Final$`5_Ratio` <- as.numeric(as.character(Final$`5_Ratio`))
  }
  
  return(Final)
}

extract_data2 <- function(folder){
  plate1 <- c("A02", "A03", "A04", "B02", "B03", "B04", "C02", "C03", "C04", "D02", "D03", "D04", 
              "E02", "E03", "E04", "F02", "F03", "F04", "G02", "G03", "G04", "Nisin", "CCCP")
  
  plate2 <- c("A05", "A06", "A07", "B05", "B06", "B07", "C05", "C06", "C07", "D05", "D06", "D07",
              "E05", "E06", "E07", "F05", "F06", "F07", "G05", "G06", "G07", "Nisin", "CCCP")
  
  plate3 <- c("A08", "A09", "A10", "B08", "B09", "B10", "C08", "C09", "C10", "D08", "D09", "D10", 
              "E08", "E09", "E10", "F08", "F09", "F10", "G08", "G09", "G10", "Nisin", "CCCP")
  
  plate4 <- c("A11", "H02", "H09", "B11", "H03", "H10", "C11", "H04", "H11", "D11", "H05", "E11",
              "H06", "F11", "H07", "G11", "H08", "Nisin", "CCCP", "N/A", "N/A", "N/A", "N/A")
  
  plateNames <- as.matrix(cbind(plate1, plate2, plate3, plate4))
  
  plates <- c("Plate1.txt", "Plate2.txt", "Plate3.txt", "Plate4.txt")
  
  full_data <- list()
  
  Final450 <- data.frame()
  Final600 <- data.frame()
  
  for(i in 1:4){
    # Read in data
    data <- read.table(plates[i], header = TRUE) 
    
    # Untreated Controls w/ culture
    culture_control600 <- vector()
    culture_control450 <- vector()
    UTControl600 <- as.numeric(data[93,3])
    UTControl450 <- as.numeric(data[93,2])
    
    # Untreated Controls w/o culture
    dye_control600 <- vector()
    dye_control450 <- vector()
    DControl600 <- as.numeric(data[94,3])
    DControl450 <- as.numeric(data[94,2])
    
    if (i != 4) {
      # Treated Samples at high concentration (with culture)
      highsamples450 <- as.numeric(data[c(seq(1,84,4),95,96),2])
      highsamples600 <- as.numeric(data[c(seq(1,84,4),95,96),3])
      
      # Treated Samples at low concentration (with culture)
      lowsamples450 <- as.numeric(data[c(seq(3,84,4),95,96),2])
      lowsamples600 <- as.numeric(data[c(seq(3,84,4),95,96),3])
      
      # Treated Controls at high concentration (no culture)
      highTControl450 <- as.numeric(data[c(seq(2,84,4),94,94),2])
      highTControl600 <- as.numeric(data[c(seq(2,84,4),94,94),3])
      
      # Treated Controls at low concentration (no culture)
      lowTControl450 <- as.numeric(data[c(seq(4,84,4),94,94),2])
      lowTControl600 <- as.numeric(data[c(seq(4,84,4),94,94),3])
      
      # Make Untreated culture the same length as samples
      culture_control600[1:23] <- UTControl600
      culture_control450[1:23] <- UTControl450
      
      # Make Untreated dye the same length as samples
      dye_control600[1:23] <- DControl600
      dye_control450[1:23] <- DControl450
      
    } else {
      # Treated Samples at high concentration (with culture)
      highsamples450 <- as.numeric(data[c(seq(1,44,4),seq(49,56,4),seq(61,68,4),seq(73,80,4),95,96),2])
      highsamples600 <- as.numeric(data[c(seq(1,44,4),seq(49,56,4),seq(61,68,4),seq(73,80,4),95,96),3])
      
      # Treated Samples at low concentration (with culture)
      lowsamples450 <- as.numeric(data[c(seq(3,46,4),seq(51,58,4),seq(63,70,4),seq(75,82,4),95,96),2])
      lowsamples600 <- as.numeric(data[c(seq(3,46,4),seq(51,58,4),seq(63,70,4),seq(75,82,4),95,96),3])
      
      # Treated Controls at high concentration (no culture)
      highTControl450 <- as.numeric(data[c(seq(2,45,4),seq(50,57,4),seq(62,69,4),seq(74,81,4),94,94),2])
      highTControl600 <- as.numeric(data[c(seq(2,45,4),seq(50,57,4),seq(62,69,4),seq(74,81,4),94,94),3])
      
      # Treated Controls at low concentration (no culture)
      lowTControl450 <- as.numeric(data[c(seq(4,47,4),seq(52,59,4),seq(64,71,4),seq(76,83,4),94,94),2])
      lowTControl600 <- as.numeric(data[c(seq(4,47,4),seq(52,59,4),seq(64,71,4),seq(76,83,4),94,94),3])
      
      # Make Untreated culture the same length as samples
      culture_control600[1:19] <- UTControl600
      culture_control450[1:19] <- UTControl450
      
      # Make Untreated dye the same length as samples
      dye_control600[1:19] <- DControl600
      dye_control450[1:19] <- DControl450
      
    }
    
    # Correction relative to untreated control
    # high concentration
    highcorrsample450 <- log2(highsamples450/UTControl450)
    highcorrsample600 <- log2(highsamples600/UTControl600)
    
    highcorrTControl450 <- log2(highTControl450/DControl450)
    highcorrTControl600 <- log2(highTControl600/DControl600)
    
    # low concentration
    lowcorrsample450 <- log2(lowsamples450/UTControl450)
    lowcorrsample600 <- log2(lowsamples600/UTControl600)
    
    lowcorrTControl450 <- log2(lowTControl450/DControl450)
    lowcorrTControl600 <- log2(lowTControl600/DControl600)
    
    # Correct for effect of compound on dye
    # high concentration
    Finalhighcorrsample450 <- highcorrsample450 - highcorrTControl450
    Finalhighcorrsample600 <- highcorrsample600 - highcorrTControl600
    
    Finallowcorrsample450 <- lowcorrsample450 - lowcorrTControl450
    Finallowcorrsample600 <- lowcorrsample600 - lowcorrTControl600
    
    # Create compound ID
    plate <- sapply(strsplit(as.character(folder), split = " "), "[", 2)
    
    if (i != 4){
      ID <- paste("OICR-CP", plate, "-", plateNames[1:23,i], sep = "")
      
      Finalcorrsample450 <- as.data.frame(cbind(ID, folder, plateNames[1:23,i], "DiOC2-3", highsamples450, lowsamples450, highTControl450, lowTControl450, culture_control450, dye_control450, Finalhighcorrsample450, Finallowcorrsample450))
      Finalcorrsample600 <- as.data.frame(cbind(ID, folder, plateNames[1:23,i], "TO-PRO-3", highsamples600, lowsamples600, highTControl600, lowTControl600, culture_control600, dye_control600, Finalhighcorrsample600, Finallowcorrsample600))
    } else {
      ID <- paste("OICR-CP", plate, "-", plateNames[1:19,i], sep = "")
      
      Finalcorrsample450 <- as.data.frame(cbind(ID, folder, plateNames[1:19,i], "DiOC2-3", highsamples450, lowsamples450, highTControl450, lowTControl450, culture_control450, dye_control450, Finalhighcorrsample450, Finallowcorrsample450))
      Finalcorrsample600 <- as.data.frame(cbind(ID, folder, plateNames[1:19,i], "TO-PRO-3", highsamples600, lowsamples600, highTControl600, lowTControl600, culture_control600, dye_control600, Finalhighcorrsample600, Finallowcorrsample600))
    }
    
    colnames(Finalcorrsample450) <- c("ID", "Plate", "Well", "Dye", "25_Sample", "5_Sample", "25_Control", "5_Control", "Untreated_culture", "Untreated_dye", "25_Ratio", "5_Ratio")
    colnames(Finalcorrsample600) <- c("ID", "Plate", "Well", "Dye", "25_Sample", "5_Sample", "25_Control", "5_Control", "Untreated_culture", "Untreated_dye", "25_Ratio", "5_Ratio")
    
    Final450 <- rbind(Final450, Finalcorrsample450)
    Final600 <- rbind(Final600, Finalcorrsample600)
    
    # Final <- list(Final450, Final600)
    Final <- rbind(Final450, Final600)
    
    # Forces values as numeric
    Final$`25_Sample` <- as.numeric(as.character(Final$`25_Sample`))
    Final$`5_Sample` <- as.numeric(as.character(Final$`5_Sample`))
    Final$`25_Control` <- as.numeric(as.character(Final$`25_Control`))
    Final$`5_Control` <- as.numeric(as.character(Final$`5_Control`))
    Final$`Untreated_culture` <- as.numeric(as.character(Final$`Untreated_culture`))
    Final$`Untreated_dye` <- as.numeric(as.character(Final$`Untreated_dye`))
    Final$`25_Ratio` <- as.numeric(as.character(Final$`25_Ratio`))
    Final$`5_Ratio` <- as.numeric(as.character(Final$`5_Ratio`))
  }
  
  return(Final)
}

plot_data <- function(data){
  for (i in 1:length(folders)){
    #Extract data
    plot_data <- data[grepl(folders[i], data$Plate),]
    
    plot_potential <- plot_data[grepl("DiOC2-3", plot_data$Dye),]
    plot_permeability <- plot_data[grepl("TO-PRO-3", plot_data$Dye),]
    
    #Generate plots
    potential <- ggplot() +
      geom_point(data = plot_potential, aes(x = plot_potential$Well, y = plot_potential$`25µM`), colour = 'black', size = 3) +
      geom_point(data = plot_potential, aes(x = plot_potential$Well, y = plot_potential$`5µM`), colour = 'red', size = 3) +
      scale_y_continuous(name="Corrected Fluorescence Ratio") +
      scale_x_discrete(name="") +
      ggtitle(paste("Membrane Potential ", plot_potential$Plate[1], sep = ""))
    
    permeability <- ggplot() +
      geom_point(data = plot_permeability, aes(x = plot_permeability$Well, y = plot_permeability$`25µM`), colour = 'black', size = 3) +
      geom_point(data = plot_permeability, aes(x = plot_permeability$Well, y = plot_permeability$`5µM`), colour = 'red', size = 3) +
      scale_y_continuous(name="Corrected Fluorescence Ratio") +
      scale_x_discrete(name="") +
      ggtitle(paste("Membrane Permeability ", plot_permeability$Plate[1], sep = ""))
    
    # Print plots as .pdf
    pdf(paste(getwd(),"/Potential ", plot_potential$Plate[1], ".pdf", sep = ""), width = 22, height = 5)
    print(potential)
    dev.off() # Required for printing PDF
    
    pdf(paste(getwd(),"/Permeability ", plot_permeability$Plate[1], ".pdf", sep = ""), width = 22, height = 5)
    print(permeability)
    dev.off() # Required for printing PDF
  
   }
}

# ================== DATA INPUT ================== 
# Set the working directory based on the plate number 

# Plate Layout 1
database1 <- data.frame()
hits <- data.frame()

folders <- list.files(paste(PROJECTDIR, "/Data/Layout 1/", sep = ""))
folders <- folders[grep("Plate ", folders)]

for (i in 1:length(folders)){
  setwd(paste(PROJECTDIR, "/Data/Layout 1/", folders[i], sep = ""))
  data <- extract_data1(folders[i])
  database1 <- rbind(database1, data)
}

# Plate Layout 2
database2 <- data.frame()
hits <- data.frame()

folders <- list.files(paste(PROJECTDIR, "/Data/Layout 2/", sep = ""))
folders <- folders[grep("Plate ", folders)]

for (i in 1:length(folders)){
  setwd(paste(PROJECTDIR, "/Data/Layout 2/", folders[i], sep = ""))
  data <- extract_data1(folders[i])
  database2 <- rbind(database2, data)
}

rm(data) # removes the data variable from the for loop
rm (i) # removes the i variable from the for loop

# Plot the variation in the control data
nisin <- database2[database2$Well == "Nisin",]
cccp <- database2[database2$Well == "CCCP",]

plot(nisin$`25_Sample`[nisin$Dye == "DiOC2-3"])
points(nisin$Untreated_culture[nisin$Dye == "DiOC2-3"], col = "red")
points(nisin$Untreated_dye[nisin$Dye == "DiOC2-3"], col = "blue")
summary(nisin$`25_Ratio`[nisin$Dye == "DiOC2-3"])

plot(nisin$`25_Sample`[nisin$Dye == "TO-PRO-3"])
points(nisin$Untreated_culture[nisin$Dye == "TO-PRO-3"], col = "red")
points(nisin$Untreated_dye[nisin$Dye == "TO-PRO-3"], col = "blue")
summary(nisin$`25_Ratio`[nisin$Dye == "TO-PRO-3"])

plot(cccp$`25_Sample`[cccp$Dye == "DiOC2-3"], ylim = c(0,50000))
points(cccp$Untreated_culture[cccp$Dye == "DiOC2-3"], col = "red")
points(cccp$Untreated_dye[cccp$Dye == "DiOC2-3"], col = "blue")
summary(cccp$`25_Ratio`[nisin$Dye == "DiOC2-3"])

plot(cccp$`25_Sample`[cccp$Dye == "TO-PRO-3"])
points(cccp$Untreated_culture[cccp$Dye == "TO-PRO-3"], col = "red")
points(cccp$Untreated_dye[cccp$Dye == "TO-PRO-3"], col = "blue")
summary(cccp$`25_Ratio`[nisin$Dye == "TO-PRO-3"])

# Is there a correlation between the increased untreated fluorescence and the treated fluorescence
plot(y = nisin$`25_Sample`[nisin$Dye == "DiOC2-3"], x = nisin$Untreated_culture[nisin$Dye == "DiOC2-3"])
plot(y = nisin$`25_Sample`[nisin$Dye == "TO-PRO-3"], x = nisin$Untreated_culture[nisin$Dye == "TO-PRO-3"])

plot(y = cccp$`25_Sample`[cccp$Dye == "DiOC2-3"], x = cccp$Untreated_culture[cccp$Dye == "DiOC2-3"])
plot(y = cccp$`25_Sample`[cccp$Dye == "TO-PRO-3"], x = cccp$Untreated_culture[cccp$Dye == "TO-PRO-3"])
# Perfect for CCCP and TO-PRO-3 but not for anythign else

# ================= Strip internal plate controls from dataset and separate into potential & permeability =================
potential_database2 <- database2[(database2$Dye == "DiOC2-3" & database2$Well != "Nisin" & database2$Well != "CCCP"),]
permeability_database2 <- database2[(database2$Dye == "TO-PRO-3" & database2$Well != "Nisin" & database2$Well != "CCCP"),]

# Plot the datasets
plot(potential_database2$`25_Ratio`)
plot(permeability_database2$`25_Ratio`)

# Potential
# Calculate just the log of the treated to untreated sample ratio
sam.pot.log <- log2(potential_database2$`25_Sample`/potential_database2$Untreated_culture)
plot(sam.pot.log)
hist(sam.pot.log)

# Which compounds interact with the dye?
comp.pot.log <- log2(potential_database2$`25_Control`/potential_database2$Untreated_dye)
plot(comp.pot.log)
hist(comp.pot.log, breaks = seq(-5, 1, 0.25))

# Filter out those that have larget effects (more negative than -0.25) on the dye alone from the total dataset
iRow <- -which(comp.pot.log < -0.25)
hist(comp.pot.log[iRow], breaks = seq(-0.5, 0.5, 0.25))

# Remove the compounds that interact with the dye from the samples
hist(sam.pot.log[iRow])

# Permeability
# Calculate just the log of the treated to untreated sample ratio
sam.per.log <- log2(permeability_database2$`25_Sample`/permeability_database2$Untreated_culture)
plot(sam.per.log)
hist(sam.per.log, breaks = seq(-6, 2.5, 0.25))

# Which compounds interact with the dye?
comp.per.log <- log2(permeability_database2$`25_Control`/permeability_database2$Untreated_dye)
plot(comp.per.log)
hist(comp.per.log, breaks = seq(-2, 2.5, 0.25))

# Filter out those that have larget effects (more negative than +/-0.25) on the dye alone from the total dataset
iRow <- -which(comp.per.log < -0.25)
iRow <- c(iRow, -which(comp.per.log > 0.25))

hist(comp.per.log[iRow], breaks = seq(-2, 2.5, 0.25))

# Remove the compounds that interact with the dye from the samples
plot(sam.per.log[iRow])
hist(sam.per.log[iRow])



# ================= Import and add data from B. subtilis inhibition screen =================
setwd("/Users/smcauley/Git/test/R Datasets")
load("Czarney_Screen.Rda")

# Creates matrix of the datablock section of the file
blockmatrix <- datablock2ma(datablocklist=datablock(sdfset))
blockframe <- as.data.frame(blockmatrix)

# Screen blockframe data for the values that are in the screen data
blockframe <- blockframe[(blockframe$`OICR CP ID` %in% potential_database$ID),]

# Ensure blockframe data is in the same order as the screen data
blockframe <- blockframe[with(blockframe, order(`OICR CP ID`)),]

# Reorder data to match the other database
potential_database <- potential_database[order(as.character(potential_database$ID)),]
permeability_database <- permeability_database[order(as.character(permeability_database$ID)),]

# Pull death results from screen for hits
Q23_squared1 <- as.double(as.character(blockframe$`B. sub Q23 Squared`))
Q23_squared2 <- as.double(as.character(blockframe$`B.sub Q23 Squared`))

# Add the B. subtilis inhibitory values to the hits database
potential_database <- cbind(potential_database, Q23_squared1, Q23_squared2)
permeability_database <- cbind(permeability_database, Q23_squared1, Q23_squared2)

# Pull hit structures as SMILES and add to database
#hit_cID <- blockframe$CdId
#structure_final <- vector()

#for (i in 1:length(hit_cID)){
#  structure_temp <- sdfset[as.integer(as.character(hit_cID[i]))]
#  structure_temp <- sdf2smiles(structure_temp)
#  structure_final <- c(structure_final, toString(structure_temp))
#  print(i)
#}

#plot(sdfset[641])

#toString(sdf2smiles(sdfset[2464]))

#as.integer(names(grepSDFset("OICR-CP31-G05", sdfset, mode="index")))

#datablock(sdfset[2464])
#datablock(sdfset[641])

# ================= Plot Results v Growth Inhibition =================
#ggplot() +
#  geom_point(data = potential_database, aes(x = Q23_squared1, y = `25_Ratio`), colour = 'black', size = 1) + 
#  geom_point(data = permeability_database, aes(x = Q23_squared1, y = `25_Ratio`), colour = 'red', size = 1)

# ================= Plot Distribution/Density of Ratios =================
#par(font.axis = 2, font.lab = 2, cex = 1.2) # bolds both the axis titles and the tick marks
#plot(density(potential_database$`25_Ratio`), 
     #main = "Density Function of Membrane Potential SSMD-scores", col = "black",
#     main = NA,
#     ylab = "Kernel Density Estimation",
#     lwd = 4,
#     xlab = "Ratio",
#     xlim = c(-0.5, 3),
#     ylim = c(0, 6))
#lines(density(potential_database$`5_Ratio`), pch=22, lty=2, lwd = 4, col = "red")

#par(font.axis = 2, font.lab = 2, cex = 1.2) # bolds both the axis titles and the tick marks
#plot(density(permeability_database$`25_Ratio`), 
     #main = "Density Function of Membrane Potential SSMD-scores", col = "black",
#     main = NA,
#     ylab = "Kernel Density Estimation",
#     lwd = 4,
#     xlab = "Ratio",
#     xlim = c(-2, 2),
#     ylim = c(0, 6))
#lines(density(permeability_database$`5_Ratio`), pch=22, lty=2, lwd = 4, col = "red")

# ================= Picking Hits =================
# This finds hits in the dataset based on a cut-off value of 0.3 for the membrane potential and 0.4 for the membrane permeability
hits_pot = potential_database[potential_database$`25_Ratio` > 0.3,]
hits_per = permeability_database[permeability_database$`25_Ratio` > 0.4,]

# Shows the results for the hits that are shared between the datasets
shared_per <- hits_per[(hits_per$ID %in% hits_pot$ID),]
shared_pot <- potential_database[unique(grep(paste(as.character(shared_per$ID),collapse="|"), potential_database$ID, value=FALSE)),]

# Shows the results for hits for only membran potential
only_pot <- hits_pot[!(hits_pot$ID %in% hits_per$ID),]
only_per <- permeability_database[unique(grep(paste(as.character(only_pot$ID),collapse="|"), permeability_database$ID, value=FALSE)),]

# Shows the results for hits for only membran permeability
difference_per <- hits_per[!(hits_per$ID %in% hits_pot$ID),]
difference_pot <- potential_database[unique(grep(paste(as.character(difference_per$ID),collapse="|"), potential_database$ID, value=FALSE)),]

# ================= Comparing with S.ven development screen =================
setwd("/Users/smcauley/Git/test")

Sven_Hits_Reg <- scan(file = "Sven_Hits.txt", what = character())
Sven_Hits <- blockframe[blockframe$Reg_Number %in% Sven_Hits_Reg,]
Sven_Hits <- database[database$ID %in% Sven_Hits$`OICR CP ID`,]

Sven_Hits_pot <- potential_database[potential_database$ID %in% Sven_Hits$ID,]
Sven_Hits_per <- permeability_database[permeability_database$ID %in% Sven_Hits$ID,]

hits_pot[hits_pot$ID %in% Sven_Hits_pot$ID,]
hits_per[hits_per$ID %in% Sven_Hits_per$ID,]

#write.xlsx(hits_pot[hits_pot$ID %in% Sven_Hits_pot$ID,], "Potential_Hits.xlsx")
#write.xlsx(hits_per[hits_per$ID %in% Sven_Hits_per$ID,], "Permeability_Hits.xlsx")
#write.xlsx(database, "Database.xlsx")

# ================= Comparing with S.ven development screen =================





ggplot() +
  geom_point(data = potential_database, aes(x = ID, y = potential_z_scores), colour = 'black', size = 3) +
  geom_point(data = potential_database[potential_database$Well == "Nisin",], aes(x = ID, y = potential_z_scores), colour = 'red', size = 3) +
  geom_point(data = potential_database[potential_database$Well == "CCCP",], aes(x = ID, y = potential_z_scores), colour = 'blue', size = 3) +
  #geom_point(data = potential_database[potential_database$Well == "Nisin"], aes(x = potential_database$ID, y = potential_database$potential_z_scores, colour = 'black', size = 3) +             
  scale_y_continuous(name="Z-Score") +
  scale_x_discrete(name="") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Permeability
permeability_database <- database[database$Dye == "TO-PRO-3",]
permeability_mean <- mean(permeability_database$`25µM`)
permeability_sd <- sd(permeability_database$`25µM`)

permeability_z_scores <- (permeability_database$`25µM`- permeability_mean)/permeability_sd

# Graph Z scores
ggplot() +
  geom_point(data = database, aes(x = c(1:nrow(database)), y = database$z_scores), colour = 'black', size = 3) +
  geom_point(data = database[database$Well], aes(x = c(1:nrow(database)), y = database$z_scores), colour = 'black', size = 3) +
  scale_y_continuous(name="Z-Score") +
  scale_x_discrete(name="")

# ================== REMOVE NEW VARIABLES ================== 
rm(list=setdiff(ls(), "PROJECTDIR"))

print("Process completed. Have a great day!")
