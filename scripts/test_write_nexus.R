#!/usr/bin/env Rscript

# Source the function from the original script
source("function.R")


# create test input data and write to file
input_file = "example_data.csv"
dat = matrix(data = c(8,11,0,0,0,
                      8,1,0,0,0,
                      13,4,5,5,3), ncol = 5, byrow = T)

rownames(dat) = c("cell_1", "cell_2", "cell_3")
colnames(dat) = c("site_1", "site_2", "site_3", "site_4", "site_5")
write.csv(dat, file = input_file)


# Define output file
output_file = "example_output_test.sciphy"

# Call the function to generate the NEXUS file
write_nexus(read.csv(input_file, row.names = 1), output_file)

# Run a bash diff command to compare the generated file with the expected output
diff_result <- system(paste("diff", "example_output.sciphy", output_file), intern = TRUE)

# Check if there are any differences
if (length(diff_result) == 0) {
  cat("Test passed: Output NEXUS file matches expected output.\n")
} else {
  cat("Test failed: Output NEXUS file does not match expected output.\n")
  cat("Differences:\n")
  cat(diff_result, sep = "\n")
}
