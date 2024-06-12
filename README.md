# Code Repo for A discovery platform to identify inducible synthetic circuitry from varied microbial sources.

Relevent data & code for "A discovery platform to identify inducible synthetic circuitry from varied microbial sources."
Raw sequencing data is not available on this github repo, it is deposited in and should be downloaded from ENA under project number: PRJEB66358. 

## Description & Folder Structure: 

### SensorSearch: 
A Python script to search MistDB for paired two component systems &amp; extract their DNA sequences. Includes example folder outputs of 'searched' genomes for E. tarda and C. rodentium, with extracted TCS sequences sorted by restriction enzyme compatibility. 

#### Demo: 
The jupyter notebook file 'Example_Running_Sensor_Search.ipynb' can be run to use SensorSearch functions to extract TCS canonical sequences from specific inputs.

### Barcodes and Sequencing Analysis
All code and associated data to analyse sequencing results (note long read sequencing alignment was done in Geneious using Minimap software - alignment output data from geneious is included in this repository). Briefly sequencing analysis is done as follows: 
1. All library screen samples are analysed for putative barcodes using Dada2 (Short_read_sequencing_Library_screening folder, LB009 corresponds to in vitro library screening, ME4 corresponds to library colonisation experiment, ME6 corresponds to in vivo Library 1 screening and ME9 corresponds to in vivo Library 2 screening)
2. All inferred barcodes & barcode counts from all library screening experiments are combined (output: combined_bcs.csv)
3. Long reads are aligned to known sensor regions using Minimap, and minimap alignments of long read sequencing are exported from Geneious & trimmed use Minimap_TCS_alignments.ipynb to get the barcode region of the alignments, which are exported as trimmed_barcode_alignments.csv. These are processed in R (Analysis_v2.R), and confirmed assigned barcodes are saved as assigned_bcs.csv. 
4. assigned_barcodes.csv (from long read sequencing) and combined_bcs.csv (short read sequencing results from each library screen - note new library screens can be run using the code as long as barcode abundances are formatted with barcode sequences in one column and abundances for each sample in adjacent columns) are merged & processed (using sample_codes.csv to specify sample specifications i.e. With or without spectinomycin etc.) to calculate fractional odds ratios, and perform QC & normalisation using specified positive controls. Various output csv results will be saved including csv with QC details, calculated OR and FOR values. The script is also compatible with pooled barcode analysis. 

#### Demo: 
Once required dependencies are installed, all R scripts can be run to re-generate the results from barcode sequencing analysis. Example lines to run analysis on experiment LB0012 are included in Analysis_v2.R script. 

#### Installation Guide: 
Programing language information & Software Dependencies: 
- SensorSearch functions are written in Python v3.8 or above, requires packages 'requests', 'pandas', 'Biopython', 'urllib3', 'reportlab' and 'os' 
- DNA Sequencing analysis is written in R v4.2.2 or above, packages required are as listed in scripts. 
- Operating System information: All scripts were run on MacOS, and should be portable to other operating systems given required coding languages & packages are installed, if required please contact the corresponding author of the linked publication to this repository for more information. 








