#Install packages
#install.packages("stringr")

#Write the target mRNA sequence file as an absolute path in ""
LoadFileName <- "/Users/isHCR/Act5C.txt"

#Assign the number of ProbeRegions you want
CandidateNum <- 40

#Set GC content condition (%)
gc_min_per <- 40
gc_max_per <- 60


#################
str(RoadFileName)
BaseFileName <- gsub(".txt", "", LoadFileName)
BaseFileName
#load DNA sequence text files
dna_sequence_pre <- scan(file = LoadFileName, what = character())
str(dna_sequence_pre)

dna_sequence <- toupper(dna_sequence_pre)
str(dna_sequence)
seqLen <-  nchar(dna_sequence)
str(seqLen)


#Load the library
library(stringr)



subsequences <- matrix(nrow = CandidateNum, ncol = 6)
colnames(subsequences) <- c("Probe_Region_Sense", "PRS GC Content", "P1_GC Content", "P2_GC Content", "StartBp_num", "EndBp_num")

str(subsequences)

gc_min <- gc_min_per * 0.01
gc_max <- gc_max_per * 0.01

# User defined function for serch
find_sequence <- function(start) {
  end <- start + 51
  PutativeProbeRegion <- substr(dna_sequence, start, end)
  P1 <- substr(PutativeProbeRegion, 1, 25)
  P2 <- substr(PutativeProbeRegion, 28, 52)
  gc_content <- sum(str_count(PutativeProbeRegion, "G|C")) / 52
  P1_gc_content <- sum(str_count(P1, "G|C")) / 25
  P2_gc_content <- sum(str_count(P2, "G|C")) / 25
  InterBase1 <- substr(PutativeProbeRegion, 26, 26)
  InterBase2 <- substr(PutativeProbeRegion, 27, 27)
  return(list(PutativeProbeRegion, gc_content, P1_gc_content, P2_gc_content, start, end, InterBase1, InterBase2))
}


count <- 1
start <- 1

#Interval between probe region
interval <- 1


while (count <= CandidateNum && start <= seqLen - 51) {
  result <- find_sequence(start)
  if (result[[3]] >= gc_min & result[[3]] <= gc_max & result[[4]] >= gc_min & result[[4]] <= gc_max & result[[7]] != "T" & result[[8]] != "T") {
    subsequences[count, "Probe_Region_Sense"] <- result[[1]]
    subsequences[count, "PRS GC Content"] <- result[[2]]
    subsequences[count, "P1_GC Content"] <- result[[3]]
    subsequences[count, "P2_GC Content"] <- result[[4]]
    subsequences[count, "StartBp_num"] <- result[[5]]
    subsequences[count, "EndBp_num"] <- result[[6]]
    count <- count + 1
    start <- start + 52 + interval
  }	else {
  	start <- start + 1
  }
}


#User defined function of creating probe region seq in fasta format 
convert_to_fasta <- function(sequence, i) {
  fasta_string <- paste(">Probe_Region_Sense", i, "\n", sequence, "\n", sep = "")
  return(fasta_string)
}

#creating probe region seq in fasta format 
ProbeRegionSence <- convert_to_fasta(subsequences[1, 1], 1)

for (i in 2:CandidateNum) {
	a <- convert_to_fasta(subsequences[i, 1], i)
	ProbeRegionSence <- paste(ProbeRegionSence, a, sep = "")
}


subsequences


#output of the result
write.table(subsequences, file= paste(BaseFileName, "_ProbeResionSence_", gc_min_per, "_", gc_max_per, ".txt", sep = ""), quote=FALSE, sep="\t", row.names=FALSE)
writeLines(ProbeRegionSence, paste(BaseFileName, "_ProbeResionSence_", gc_min_per, "_", gc_max_per, ".fasta", sep = ""))


