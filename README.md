# Data & Code Repo for A discovery platform to identify inducible synthetic circuitry from varied microbial sources.

Relevent data & code for "A discovery platform to identify inducible synthetic circuitry from varied microbial sources."

## Description & Folder Structure: 

### SensorSearch: 
A Python script to search MistDB for paired two component systems &amp; extract their DNA sequences. Includes example folder outputs of 'searched' genomes for E. tarda and C. rodentium, with extracted TCS sequences sorted by restriction enzyme compatibility. 

#### Demo: 
The jupyter notebook file 'Example_Running_Sensor_Search.ipynb' can be run to use SensorSearch functions to extract TCS canonical sequences from specific inputs.

### Barcodes and Sequencing Analysis
All code and associated data to analyse sequencing results (note long read sequencing alignment was done in Geneious using Minimap software - alignment output data from geneious is included in this repository). Briefly sequencing analysis is done as follows: 
1. All library screen samples are analysed for putative barcodes using Dada2 (Short_read_sequencing_Library_screening folder, LB009 corresponds to in vitro library screening, data processing, ME4 corresponds to library colonisation, ME6 corresponds to in vivo library screening)
2. All inferred barcodes & barcode counts from all library screening experiments are combined & QCed (script: combining_all_barcodes.R, output: combined_bcs.csv)
3. Minimap alignments of long read sequencing are exported from Geneious & trimmed use Minimap_TCS_alignments.ipynb to get the barcode region of the alignments, which are exported as trimmed_barcode_alignments.csv. These are processed using TCS_to_barcode_assignment.R, and confirmed assigned barcodes are saved as assigned_barcodes.csv. 
4. assigned_barcodes.csv (from long read sequencing) and combined_bcs.csv (short read sequencing results from each library screen) are merged & processed to calculate odds ratios in Odds_ratio_plotting.R, which further performs QC & normalisation using positive controls. 

#### Demo: 
Once required dependencies are installed, all R scripts can be run to re-generate the results from barcode sequencing analysis. 

## Installation Guide: 
Programing language information & Software Dependencies: 
- SensorSearch functions are written in Python v3.8 or above, requires packages 'requests', 'pandas', 'Biopython', 'urllib3', 'reportlab' and 'os' 
- DNA Sequencing analysis is written in R v4.2.2 or above, packages required are as listed in scripts. 
Operating System information: 
- All scripts were run on MacOS 








