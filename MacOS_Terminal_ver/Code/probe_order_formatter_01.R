#!/usr/bin/env Rscript

# Error handling function
handle_error <- function(msg) {
  cat("Error: ", msg, "\n")
  cat("Usage: Rscript probe_order_formatter.R <input_file.csv>\n")
  quit(status = 1)
}

# Check command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  handle_error("Please provide exactly one input file.")
}

input_file <- args[1]

# Check if input file exists
if (!file.exists(input_file)) {
  handle_error("Input file does not exist.")
}

# Generate output file name
output_file <- sub("\\.csv$", "_order.csv", input_file)

# Read input file
data <- tryCatch(
  read.csv(input_file, header = TRUE, stringsAsFactors = FALSE),
  error = function(e) {
    handle_error(paste("Failed to read input file:", e$message))
  }
)

# Check if required columns exist
required_columns <- c("Probe_Region", "Start_Site", "P1", "P2")
if (!all(required_columns %in% colnames(data))) {
  handle_error("Input file is missing required columns.")
}

# Extract hairpin DNA ID and gene name from the input file name
file_parts <- strsplit(basename(input_file), "_")[[1]]
hairpin_id <- file_parts[2]
gene_name <- file_parts[3]

# Prepare output data
output_data <- data.frame(
  Oligoname = c(paste0(gene_name, "_", hairpin_id, "_s", data$Start_Site, "_P1"),
                paste0(gene_name, "_", hairpin_id, "_s", data$Start_Site, "_P2")),
  sequence = c(data$P1, data$P2),
  stringsAsFactors = FALSE
)

# Sort the output data
output_data <- output_data[order(output_data$Oligoname), ]

# Write output file
tryCatch({
  write.table(data.frame(Oligoname = "Oligoname", sequence = "sequence"), 
              file = output_file, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(output_data, file = output_file, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}, error = function(e) {
  handle_error(paste("Failed to write output file:", e$message))
})

cat("Probe order file has been created:", output_file, "\n")
