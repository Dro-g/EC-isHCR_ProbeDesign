#!/usr/bin/env Rscript

# Function to safely install and load packages
install_and_load <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    repos <- "https://cloud.r-project.org"
    install.packages(package, repos = repos, quiet = TRUE)
  }
  library(package, character.only = TRUE)
}

# Install and load required packages
install_and_load("optparse")
install_and_load("stringr")

# Define command line options
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input FASTA file", metavar="FILE"),
  make_option(c("-m", "--min_gc"), type="integer", default=45, 
              help="Minimum GC content percentage [default= %default]", metavar="NUMBER"),
  make_option(c("-M", "--max_gc"), type="integer", default=55, 
              help="Maximum GC content percentage [default= %default]", metavar="NUMBER")
)

# Create usage message
usage <- "Usage: Rscript probe_region_finder.R [options] <input_file>\n\nOptions:"

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list, usage=usage)
args <- commandArgs(trailingOnly = TRUE)

# Check if no arguments are provided
if (length(args) == 0) {
  cat("Probe Region Finder\n\n")
  cat("This script finds probe regions in a DNA sequence based on GC content.\n\n")
  cat("Usage: Rscript probe_region_finder.R [options] <input_file>\n\n")
  cat("Options:\n")
  cat("  -i, --input FILE    Input FASTA file (required if not specified as last argument)\n")
  cat("  -m, --min_gc NUMBER Minimum GC content percentage (default: 45)\n")
  cat("  -M, --max_gc NUMBER Maximum GC content percentage (default: 55)\n")
  cat("  -h, --help          Show this help message and exit\n\n")
  cat("Example:\n")
  cat("  Rscript probe_region_finder.R input.fasta\n")
  cat("  Rscript probe_region_finder.R --input input.fasta --min_gc 40 --max_gc 60\n")
  cat("  Rscript probe_region_finder.R -m 40 -M 60 input.fasta\n")
  quit(status = 0)
}

# Parse options
opt <- parse_args(opt_parser, args = args, positional_arguments = TRUE)

# Check if input file is provided as a positional argument
if (length(opt$args) > 0) {
  opt$options$input <- opt$args[length(opt$args)]
} else if (is.null(opt$options$input)) {
  # If no input file is provided, use the last argument as input file
  opt$options$input <- args[length(args)]
}

# Check if input file is provided
if (is.null(opt$options$input)) {
  stop("Input file not specified. Use --input option or provide the file as the last argument.")
}

# Access parsed options
input_file <- opt$options$input
gc_min_per <- opt$options$min_gc
gc_max_per <- opt$options$max_gc

# Function to read DNA sequence from file and check FASTA format
read_dna_sequence <- function(file_path) {
  lines <- readLines(file_path)
  
  # Check if the file is in FASTA format
  if (!grepl("^>", lines[1])) {
    warning("WARNING: The input file may not be in FASTA format. Please check your input file.")
  }
  
  # Remove the header line(s)
  seq_lines <- lines[!grepl("^>", lines)]
  
  dna_sequence <- paste(seq_lines, collapse = "")
  return(toupper(dna_sequence))
}

# Read DNA sequence
dna_sequence <- read_dna_sequence(input_file)
seqLen <- nchar(dna_sequence)

cat("Input file:", input_file, "\n")
cat("Minimum GC content:", gc_min_per, "%\n")
cat("Maximum GC content:", gc_max_per, "%\n")
cat("Sequence length:", seqLen, "bp\n")

# Initialize variables
gc_min <- gc_min_per * 0.01
gc_max <- gc_max_per * 0.01
subsequences <- data.frame(Probe_Region_Sense = character(),
                           PRS_GC_Content = numeric(),
                           P1_GC_Content = numeric(),
                           P2_GC_Content = numeric(),
                           StartBp_num = integer(),
                           EndBp_num = integer(),
                           stringsAsFactors = FALSE)

# User defined function for search
find_sequence <- function(start) {
  end <- start + 51
  if (end > seqLen) return(NULL)  # Return NULL if we've reached the end of the sequence
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

# Search for probe regions
start <- 1
interval <- 1

while (start <= seqLen - 51) {
  result <- find_sequence(start)
  if (!is.null(result) &&
      result[[3]] >= gc_min && result[[3]] <= gc_max && 
      result[[4]] >= gc_min && result[[4]] <= gc_max && 
      result[[7]] != "T" && result[[8]] != "T") {
    subsequences <- rbind(subsequences, data.frame(
      Probe_Region_Sense = result[[1]],
      PRS_GC_Content = result[[2]],
      P1_GC_Content = result[[3]],
      P2_GC_Content = result[[4]],
      StartBp_num = result[[5]],
      EndBp_num = result[[6]],
      stringsAsFactors = FALSE
    ))
    start <- start + 52 + interval
  } else {
    start <- start + 1
  }
}

# User defined function of creating probe region seq in fasta format 
convert_to_fasta <- function(sequence, i) {
  fasta_string <- paste0(">Probe_Region_Sense_", i, "\n", sequence, "\n")
  return(fasta_string)
}

# Creating probe region seq in fasta format 
ProbeRegionSense <- sapply(1:nrow(subsequences), function(i) {
  convert_to_fasta(subsequences$Probe_Region_Sense[i], i)
})
ProbeRegionSense <- paste(ProbeRegionSense, collapse = "")

# Get base filename without extension
BaseFileName <- tools::file_path_sans_ext(input_file)

# Output of the result
write.table(subsequences, file = paste0(BaseFileName, "_ProbeRegionSense_", gc_min_per, "_", gc_max_per, ".txt"), 
            quote = FALSE, sep = "\t", row.names = FALSE)
writeLines(ProbeRegionSense, paste0(BaseFileName, "_ProbeRegionSense_", gc_min_per, "_", gc_max_per, ".fasta"))

cat("Probe regions have been successfully identified and saved.\n")
cat("Number of probe regions found:", nrow(subsequences), "\n")