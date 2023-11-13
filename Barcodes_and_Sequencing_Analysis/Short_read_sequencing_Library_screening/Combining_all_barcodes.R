library(dplyr)

#For R Studio only: set workind directory to location of script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#read all barcodes determined by dada2 from each run
LB012_bcs <- read.csv('LB012/LB012_barcode_abundances.csv')
ME1_bcs <- read.csv('ME1/Barcode_abundances_ME1.csv')
ME4_bcs <- read.csv('ME4/Barcode_abundances_ME4.csv')
ME6_bcs <- read.csv('ME6/Barcode_abundances_ME6.csv')
LB009_bcs <- read.csv('LB009/Barcode_abundances_LB009.csv')

#read individually sanger sequenced barcodes that already have their assignment
ind_bcs <- read.csv('individually_sequenced_bcs.csv')

#make a list of all barcodes for all runs
df_list <- list(subset(LB012_bcs, select=-X),
                subset(ME1_bcs, select=-X),
                subset(ME4_bcs, select=-X),
                subset(ME6_bcs, select=-X),
                subset(LB009_bcs, select=-X),
                ind_bcs)

#and combine all barcodes & barcodes with known assignment 
combined_bcs <- df_list %>% purrr::reduce(dplyr::full_join, by='barcode')
combined_bcs <- combined_bcs %>%
  dplyr::select(sensor, everything()) #move sensor column to left side of df for easier viewing 

#cluster barcodes based on >95% similarity, using seq_cluster from bioseq package. 
barcodes <- bioseq::dna(c(combined_bcs$barcode))
combined_bcs$cluster <- bioseq::seq_cluster(barcodes,
                                    threshold = 0.05)

#get total number of reads for each Dada2 computed barcode ASV. 
combined_bcs$total_bc_reads <- rowSums(combined_bcs[ , 3:(ncol(combined_bcs)-1)],na.rm=TRUE)  

#select member of each cluster with highest number of reads and keep that DNA barcode. 
#get total reads per sample per cluster
com_bcs_clust_reads <- combined_bcs %>%
  dplyr::select(c(cluster,total_bc_reads), everything()) %>% #move cluster & total_bc_reads columns to left side of df for easier viewing. 
  dplyr::group_by(cluster) %>% 
  dplyr::arrange(desc(total_bc_reads), .by_group = TRUE) %>%
  dplyr::mutate_at(4:(ncol(combined_bcs)-1), ~tidyr::replace_na(.,0)) %>% #replace NAs with 0 so that they will be included in sums of reads for clusters of barcodes
  dplyr::summarise(across(where(is.numeric), sum)) #add all the reads from individual barcode clusters together per sample

#get dominant barcode DNA sequence per cluster
com_bcs_seq <- combined_bcs %>%
  dplyr::select(c(cluster,total_bc_reads), everything()) %>% #move cluster & total_bc_reads columns to front of df for easier viewing. 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(1, total_bc_reads) %>% #get DNA barcode with highest number of reads per 
  dplyr::select(sensor, barcode, cluster)

#combine total reads per cluster per sample with dominant (assumed correct) DNA barcode. 
com_bcs <- merge(com_bcs_seq,com_bcs_clust_reads, by='cluster')

#save barcodes, these are considered final barcodes present in the library. 
write.csv(com_bcs, file='combined_bcs.csv')

### Additional Analysis: Seeing how many reads are given to other barcodes in a cluster (other than main cluster barcode)
#get subtotals of reads from non main barcode: 
com_bcs_clust_reads2 <- combined_bcs %>%
  dplyr::select(c(cluster,total_bc_reads), everything()) %>% #move cluster & total_bc_reads columns to left side of df for easier viewing. 
  dplyr::group_by(cluster) %>% 
  dplyr::arrange(desc(total_bc_reads), .by_group = TRUE) %>%
  dplyr::slice(2:n()) %>% #get the barcodes read totals for all but the most abundant barcode
  dplyr::mutate_at(4:(ncol(combined_bcs)-1), ~tidyr::replace_na(.,0)) %>% #replace NAs with 0 so that they will be included in sums of reads for clusters of barcodes
  dplyr::summarise(across(where(is.numeric), sum)) #add all the reads from individual barcode clusters together per sample

#combine that with number of barcodes for that cluster (freq)
com_bcs_clust_reads2 <- cbind(com_bcs_clust_reads2, data.frame(table(combined_bcs$cluster)))%>%  #get number of barcodes per cluster, and merge that with the previous dataframe
  dplyr::select(c(Freq), everything()) %>%
  dplyr::rename(non_dom_bc_reads = total_bc_reads) %>%
  select(Freq, cluster, non_dom_bc_reads)

#comparison_df is a dataframe of the fraction of reads from a cluster that aren't from the 'dominant' barcode, the number of reads range from 4*10e-3 to 10e-5. 
comparison_df <- merge(com_bcs_clust_reads[c('cluster','total_bc_reads')], com_bcs_clust_reads2, by='cluster') %>%
  dplyr::mutate(fraction_non_dom = non_dom_bc_reads/total_bc_reads) %>%
  dplyr::group_by(cluster) %>% 
  dplyr::slice(1)%>%
  dplyr::filter(Freq!=1)
