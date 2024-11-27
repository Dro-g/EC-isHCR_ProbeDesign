#Write the table file name (.txt) with an absolute path
RegionFile <- "/Users/Act5C_ProbeResionSence_45_55.txt"

#file name to save
nameL <- strsplit(RegionFile, "/")
str(nameL)
nameC <- as.vector(nameL[[1]])
str(nameC)
length(nameC)
dir <- dirname(RegionFile)
NameHead <- paste(dir, "/Probe_A161_", sep = "")
saveName <- paste0(NameHead, gsub(".txt", ".csv", nameC[length(nameC)]), sep = "")
str(saveName)

#load table
ProbeRegionFile <- read.delim(RegionFile, stringsAsFactors = F)
str(ProbeRegionFile)

#create dataframe
ProbeRegionFile2 <- data.frame(ProbeRegionFile$Probe_Region_Sense, ProbeRegionFile$StartBp_num)
str(ProbeRegionFile2)
colnames(ProbeRegionFile2) <- c("Probe_Region", "Start_Site")
str(ProbeRegionFile2)

#Add rows for probe
ProbeRegionFile2$P1 <- rep("NNN", nrow(ProbeRegionFile2))
str(ProbeRegionFile2)
ProbeRegionFile2$P2 <- rep("NNN", nrow(ProbeRegionFile2))
str(ProbeRegionFile2)

P1AddSeq <- "GGTACGCGAaa"
P2AddSeq <- "aaAGGTAGGTGTAA"


P1Region <- substr(ProbeRegionFile2$Probe_Region, 1, 25)
P2Region <- substr(ProbeRegionFile2$Probe_Region, 28, 52)
str(P1Region)
str(P2Region)

#RevCom and add initiator seq
P1Com <- chartr("ATGC", "TACG", P1Region)
str(P1Com)
P1RevCom <- sapply(lapply(strsplit(P1Com, NULL), rev), paste, collapse="")
str(P1RevCom)
P1Seq  <- paste(P1AddSeq, P1RevCom, sep = "")
str(P1Seq)
ProbeRegionFile2$P1 <- P1Seq

P2Com <- chartr("ATGC", "TACG", P2Region)
str(P2Com)
P2RevCom <- sapply(lapply(strsplit(P2Com, NULL), rev), paste, collapse="")
str(P2RevCom)
P2Seq  <- paste(P2RevCom, P2AddSeq, sep = "")
str(P2Seq)
ProbeRegionFile2$P2 <- P2Seq

str(ProbeRegionFile2)

write.csv(ProbeRegionFile2, file = saveName)
