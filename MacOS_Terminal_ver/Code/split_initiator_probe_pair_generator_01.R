#!/usr/bin/env Rscript

# Function to safely install and load packages
install_and_load <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  library(package, character.only = TRUE)
}

# Install and load required packages
install_and_load("optparse")

# Predefined sequences
sequence_dict <- list(
  S23 = "GGGTGGTCGTCGAAGTCGTAT",
  S41 = "GCTCGACGTTCCTTTGCAACA",
  S45 = "CCTCCACGTTCCATCTAAGCT",
  S72 = "CGGTGGAGTGGCAAGTAGGAT",
  S73 = "CGGTCAGGTGGCTAGTATGGA",
  A161 = "GGTACGCGAAGGTAGGTGTAA"
)

# Command line option parsing
option_list <- list(
  make_option(c("-i", "--id"), type="character", 
              help="Hairpin DNA ID (S23, S41, S45, S72, S73, and A161 are predefined sequences)", metavar="CHARACTER"),
  make_option(c("-s", "--initiator1_seq"), type="character", 
              help="Initiator1 sequence (optional)", metavar="CHARACTER")
)

opt_parser <- OptionParser(option_list=option_list)
args <- parse_args(opt_parser, positional_arguments = 1)

if (length(args$args) != 1) {
  cat("Usage: Rscript script.R -i <ID> [options] <input_file.fasta>\n")
  stop("Input FASTA file must be supplied as the last argument", call.=FALSE)
}

input_file <- args$args[1]
opt <- args$options

if (is.null(opt$id) && is.null(opt$initiator1_seq)) {
  stop("Either --id or --initiator1_seq must be specified", call.=FALSE)
}

# Function to split initiator sequence
split_initiator <- function(seq) {
  p1 <- paste0(substr(seq, 1, 9), "aa")
  p2 <- paste0("aa", substr(seq, 10, nchar(seq)))
  return(list(P1AddSeq = p1, P2AddSeq = p2))
}

# Determine the initiator sequence
if (!is.null(opt$id) && opt$id %in% names(sequence_dict)) {
  initiator_seq <- sequence_dict[[opt$id]]
} else if (!is.null(opt$initiator1_seq)) {
  initiator_seq <- opt$initiator1_seq
} else {
  stop("Invalid ID or initiator sequence not provided", call.=FALSE)
}

# Split the initiator sequence
initiator_parts <- split_initiator(initiator_seq)
P1AddSeq <- initiator_parts$P1AddSeq
P2AddSeq <- initiator_parts$P2AddSeq

# File naming
dir <- dirname(input_file)
NameHead <- paste0(dir, "/Probe_", ifelse(!is.null(opt$id), opt$id, "custom"), "_")
saveName <- paste0(NameHead, gsub(".txt", ".csv", basename(input_file)))

# Extract gene name from input file name
gene_name <- gsub("_ProbeRegionSense_.*", "", basename(input_file))
gene_name <- gsub("\\.txt$", "", gene_name)

# Load and process data
tryCatch({
  ProbeRegionFile <- read.delim(input_file, stringsAsFactors = FALSE)
  if (!"Probe_Region_Sense" %in% colnames(ProbeRegionFile) || !"StartBp_num" %in% colnames(ProbeRegionFile)) {
    stop("Required columns 'Probe_Region_Sense' and/or 'StartBp_num' not found in the input file")
  }
  ProbeRegionFile2 <- data.frame(
    Probe_Region = ProbeRegionFile$Probe_Region_Sense,
    Start_Site = ProbeRegionFile$StartBp_num,
    P1 = rep("NNN", nrow(ProbeRegionFile)),
    P2 = rep("NNN", nrow(ProbeRegionFile))
  )
}, error = function(e) {
  cat("Error reading the input file:", conditionMessage(e), "\n")
  cat("Please check if the file is tab-delimited and contains the required columns.\n")
  stop("File reading error", call. = FALSE)
})

P1Region <- substr(ProbeRegionFile2$Probe_Region, 1, 25)
P2Region <- substr(ProbeRegionFile2$Probe_Region, 28, 52)

# Reverse complement and add initiator sequences
revcom <- function(x) {
  chartr("ATGC", "TACG", sapply(lapply(strsplit(x, NULL), rev), paste, collapse=""))
}

ProbeRegionFile2$P1 <- paste0(P1AddSeq, revcom(P1Region))
ProbeRegionFile2$P2 <- paste0(revcom(P2Region), P2AddSeq)

# Write results
write.csv(ProbeRegionFile2, file = saveName, row.names = FALSE)

cat("Results saved to:", saveName, "\n")