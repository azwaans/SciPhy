write_nexus = function(cell_by_site_dat, output_file){
  
  # Basic checks on input data
  if (!is.matrix(cell_by_site_dat) && !is.data.frame(cell_by_site_dat)) {
    stop("Error: 'cell_by_site_dat' must be a matrix or a data frame.")
  }
  
  if (!all(sapply(cell_by_site_dat, is.numeric))) {
    stop("Error: All elements of 'cell_by_site_dat' must be numeric.")
  }
  
  # Test if we can write to the output file or if it exists
  if (file.exists(output_file) && file.access(output_file, mode = 2) != 0) {
    stop("Error: The output file is not writable.")
  }
  
  n_taxa = nrow(cell_by_site_dat)
  n_barcodes = ncol(cell_by_site_dat)
  
  tax_labels = seq(0, (n_taxa-1))
  
  #Write nexus header
  write(x = "#NEXUS", output_file)
  write(x = " ", output_file, append = T)
  
  # write taxa block
  write(x = "begin taxa;", output_file, append = T)
  write(x = paste0("dimensions ntax=", n_taxa, ";"), output_file, append = T)
  write(x= paste(c("taxlabels", tax_labels, ";")), output_file, sep = " ", append = T, ncolumns = (n_taxa+2))
  write(x = "end;", output_file, append = T)
  write(x = " ", output_file, append = T)
  
  # write data block
  write(x="begin characters;", output_file, append = T)
  write(x=paste0("dimensions nchar=", n_barcodes, ";"), output_file, append=T)
  write(x="format datatype=integer;", output_file, append = T)
  write(x="matrix", output_file, append = T) 
  
  for (taxon_nr in tax_labels){
    sequence_concatenated = paste(cell_by_site_dat[(taxon_nr+1), ], collapse = ",")
    write(x = paste(taxon_nr, sequence_concatenated), output_file, append = T)
  }
  write(x=";", output_file, append = T)
  write(x="end;", output_file, append = T)
}
