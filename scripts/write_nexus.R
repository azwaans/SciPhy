#!/usr/bin/env Rscript

source("function.R")

# Get command-line arguments
args = commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 2) {
  stop("Usage: Rscript write_nexus.R <input_data_file> <output_file>")
}

# Read arguments
input_data_file = args[1]
output_file = args[2]

# Load the data (assuming it's in CSV format, but adjust as needed)
cell_by_site_dat = read.csv(input_data_file, row.names = 1)

write_nexus(cell_by_site_dat, output_file)

