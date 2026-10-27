# README

This repository contains scripts to generate NEXUS input data for use with SciPhy.

## Input data
The starting point is the example_data.csv file, where each row represents a cell and contains its associated barcode information. In this example, we handle sequential editing data across 5 consecutive sites per cell.

## Converting to NEXUS Format
To convert the example_data.csv into NEXUS format, run the conversion script as shown below. The resulting file, here example_output.sciphy, can be loaded into BEAUti for setting up the analysis.

` Rscript write_nexus.R example_data.csv example_output.sciphy`

## Example xml
An example XML file generated from the example_data.csv file is available under examples/example.xml.


