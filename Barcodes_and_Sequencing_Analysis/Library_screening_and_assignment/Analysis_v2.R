# This code is to replicate analyses and figures for Robinson et al. 2024. Code developed by Clare Robinson. 


######### Loading required packages ###########
library(dplyr)
library(ggplot2)
library(purrr)

#For R Studio only: set working directory to location of script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######### Combining short reads from all Dada2 runs ###########
# This is used as a reference file against which long reads will be aligned. 
#read all barcodes determined by dada2 from each run
ME4_bcs <- read.csv('ME4/Barcode_abundances_ME4.csv')
ME6_bcs <- read.csv('ME6/Barcode_abundances_ME6.csv')
LB009_bcs <- read.csv('LB009/Barcode_abundances_LB009.csv')
ME9_bcs <- read.csv('ME9/Barcode_abundances_ME9.csv')

#make a list of all barcodes for all runs
df_list <- list(subset(LB012_bcs, select=-X),
                subset(ME1_bcs, select=-X),
                subset(ME4_bcs, select=-X),
                subset(ME6_bcs, select=-X),
                subset(LB009_bcs, select=-X))

#and combine all barcodes & barcodes with known assignment 
combined_bcs <- df_list %>% purrr::reduce(dplyr::full_join, by='barcode')
combined_bcs <- combined_bcs %>%
  dplyr::select(sensor, everything()) #move sensor column to left side of df for easier viewing 
write.csv(combined_bcs, 'combined_bcs.csv')

######### Long read assignment processing and clustering ###########

#This only needs to be run once (not repeatedly for every short read library screening experiment)
#Takes as input alignment csv of trimmed long reads aligned to Dada2 generated barcodes (alignment was done in Geneious), clusters barcodes & adds long reads within each cluster
#Assigns barcodes based on total number of long reads or proportional number of long reads & saves assigned bcs as 'assigned_bcs.csv' in 0_Long_read_assignment folder. 

Long_read_clustering_and_assignment <- function(geneious_alignments_csv, Dada2_sequences_aligned_to_in_geneious_csv, ind_seq_bcs){
  # Read in alignments (this is exported directly from Geneious)
  alignments <- read.csv(geneious_alignments_csv)
  
  # Read in DADA2 barcode sequences & merge to long read alignments (these are all Dada2 barcodes that were imported into Geneious to align long reads to)
  bc_sequences <- read.csv(Dada2_sequences_aligned_to_in_geneious_csv)
  alignments <- merge(alignments, bc_sequences, by='Name')
  
  #clean up sensor name from description column 
  alignments$sensor <- sapply(strsplit(alignments$Description, split='mapped'),'[',1)
  alignments$sensor <- gsub(".*from ","",alignments$sensor)
  alignments$sensor <- gsub(".fasta ","",alignments$sensor)
  alignments$sensor <- gsub("NZ_CP022412.2 -","",alignments$sensor)
  alignments$sensor <- gsub(" \\(reversed\\)","",alignments$sensor)
  alignments$sensor <- gsub("Filtered","",alignments$sensor)
  alignments <- alignments %>%
    dplyr::rename( n_seq = X..Sequences , 
                   per_of_refseq = X..Of.Ref.Seq ,
                   mean_coverage = Mean.Coverage)%>%
    mutate(per_of_refseq2 = as.numeric(gsub('%','',per_of_refseq)))%>%
    mutate(HQ = as.numeric(gsub('%','',X..HQ)))%>%
    mutate(LQ = as.numeric(gsub('%','',X..LQ)))
  alignments$library <- gsub('Filtered','',sapply(strsplit(alignments$sensor, split='_'),'[',1))
  
  #rename control sensors from plasmid names to control sensors names (for ex. plasmid CR21 pCR21 -> yeaR23)
  control_name_pattern <- c('pCR24 PCR Product','pCR23  PCR Product','pCR22  PCR Product','pCR21 PCR Product','pCR18  PCR Product','pCR17 PCR Product','pCR16  PCR Product','pCR15  PCR Product')
  control_name_replacement <- c('torC23','torC17','torC10','yeaR23','torCWT','hycAWT','ynfE17','ynfE15')
  alignments$sensor <- stringi::stri_replace_all_regex(alignments$sensor ,
                                                       pattern=control_name_pattern,
                                                       replacement =control_name_replacement,
                                                       vectorize=FALSE)
  
  # rename 'Sequence' column (which contains barcode DNA sequence) to 'barcode' 
  alignments <- alignments %>%
    dplyr::rename(barcode = Sequence)
  
  # read individually sanger sequenced barcodes that already have their assignment
  ind_bcs <- read.csv(ind_seq_bcs)%>%
    rename(sanger_seq_sensor = sensor)
  
  # combine long read assignment with sanger sequence assignment: 
  assignments <- list(alignments, ind_bcs) %>% purrr::reduce(dplyr::full_join, by='barcode') #join all 3 dfs (dada2 read counts, long read assignment & individually sequenced sensors)
  
  # Cluster barcodes based on >90% similarity, using seq_cluster from bioseq package -> the Dada2 runs used to barcodes generate 'ASVs' -> in practice observed this is likely too stringent -> merge barcodes with >90% similarity. 
  barcodes <- bioseq::dna(unique(c(assignments$barcode))) #Get bioseq format list of barcode DNA sequences 
  clustered_bcs <- bioseq::seq_cluster(barcodes, threshold = 0.1) #cluster on 10% similarity threshold
  clustered_bcs_df <- data.frame(barcode = barcodes, cluster = clustered_bcs)
  #write.csv(clustered_bcs_df, '2_Processing_intermediates_csv/lr_clustered_bcs_df.csv')
  
  # Combine clusters with original long read info df
  assignments <- inner_join(assignments, clustered_bcs_df, by='barcode')
  
  # Add subclusters column in data - this column is simply a count of the barcode number within the cluster (i.e. for cluster 1 if there are 30 barcodes, rows will be numbered within the cluster from 1 to 30)
  assignments2 <- assignments %>%
    dplyr::select(c(cluster),everything())%>%
    dplyr::group_by(cluster) %>% 
    mutate(subcluster_count = cumsum(barcode != lag(barcode, default = first(barcode)))+1) %>% #add a subcluster column with numbering of barcode's subcluster. 
    ungroup()%>% 
    select(cluster, subcluster_count, everything()) 
  
  # Pivot long read data wider based on subcluster (so that number of long reads can be totaled per sensor per cluster)
  Wide_subcluster_assigned_bcs <- assignments2 %>%
    select(c(cluster, subcluster_count, n_seq,sensor, sanger_seq_sensor))%>%
    tidyr::pivot_wider(names_from = subcluster_count, values_from = c(n_seq))
  
  # Totaling the long reads from a specific sensor (which are separated in subclusters by column) into a single column 
  cols_to_remove <- which(as.numeric(colnames(Wide_subcluster_assigned_bcs)) %in% 1:1000) # Identify any numeric columns named 1 to 1000 to remove
  Wide_subcluster_assigned_bcs$Long_reads_for_assignment <- rowSums(Wide_subcluster_assigned_bcs[, cols_to_remove], na.rm = TRUE) # Add the row_sums as a new column to the modified dataframe
  Wide_subcluster_assigned_bcs <- Wide_subcluster_assigned_bcs[, -cols_to_remove] # Remove these columns from the dataframe
  
  # Sensor barcode assignment by proportional top number of reads (diff TCS & barcodes have diff abundances, so will have different numbers of misassigned reads. Normalising per TCS abundance can help with barcode mis-assignment)
  top_proportional <- Wide_subcluster_assigned_bcs %>%
    dplyr::group_by(sensor) %>%
    dplyr::mutate(tot_long_rds_sensor = sum(Long_reads_for_assignment))%>% #get total reads per sensor
    dplyr::mutate(prop_long_rds_barcode=Long_reads_for_assignment/tot_long_rds_sensor)%>% #normalise reads for a barcode by it's total reads per assigned sensor
    dplyr::group_by(cluster)%>% #group by bc
    dplyr::slice_max(prop_long_rds_barcode)%>% #select TCS assigned with the higher proportional number of reads
    dplyr::mutate(top_prop=sensor)%>%
    dplyr::select(c(cluster,tot_long_rds_sensor,prop_long_rds_barcode, top_prop))
  
  # Sensor barcode assignment by top total number of reads
  top_by_reads_df <- Wide_subcluster_assigned_bcs%>%
    dplyr::group_by(cluster)%>% #group by bc
    dplyr::slice_max(Long_reads_for_assignment)%>% #select TCS with highest read assigned to that bc
    dplyr::mutate(top_by_reads=sensor)%>%
    dplyr::select(c(cluster,top_by_reads))
  
  # Merge alignments & assignment strategies
  alignments <- merge(merge(Wide_subcluster_assigned_bcs, top_by_reads_df, by='cluster', all=T), top_proportional, by='cluster', all=T)
  
  # Create column for if alignments agree or disagree 
  alignments <- alignments %>%
    mutate (dis = if_else(top_by_reads == top_prop, 'agreement','disagreement'))%>%
    mutate (dis = if_else(is.na(dis), 'manually_assigned', dis))
  
  # Get actual DNA sequence for highest number of reads barcode in each cluster
  barcode_cluster_seqs <- assignments2 %>%
    group_by(cluster) %>%
    slice_max(n_seq, with_ties = FALSE) %>% #select row with highest number of aligned reads in cluster
    select(cluster, barcode)
  
  # Merge DNA sequence of barcode with long read counts by cluster number: 
  Wide_subcluster_assigned_bcs2 <- inner_join(barcode_cluster_seqs, alignments, by='cluster')
  
  # Single replica df - this dataframe keeps only the main barcode per cluster. 
  assigned_bcs <- Wide_subcluster_assigned_bcs2%>%
    dplyr::filter(dis=='agreement'|dis=='manually_assigned')%>%
    group_by(cluster) %>%
    filter(!(is.na(top_prop) & n() > 1) | n() == 1|dis=='manually_assigned') %>%
    slice_max(Long_reads_for_assignment) %>% # Take the first row of each group - this is the assigned sensor alignment row
    ungroup()%>%
    dplyr::filter(Long_reads_for_assignment >=6|dis=='manually_assigned') #for strict assignment we've filtered this to be >=6 reads for assignment (Geneious includes reference sequence in the sequence count -> so actually filtering for 5 long reads aligned)
  
  # Consolidating sanger sequenced sensors & long read assigned sensors into a single 'sensor' column 
  assigned_bcs[is.na(assigned_bcs$sensor) & !is.na(assigned_bcs$sanger_seq_sensor), "sensor"] <- assigned_bcs[is.na(assigned_bcs$sensor) & !is.na(assigned_bcs$sanger_seq_sensor), "sanger_seq_sensor"]
  assigned_bcs <- assigned_bcs %>%
    select(-c(sanger_seq_sensor,top_by_reads, top_prop))%>%
    mutate_at(vars(sensor), ~ stringr::str_replace(., "Lib1_", ""))%>% #clean up long read assignment columns 
    mutate_at(vars(sensor), ~ stringr::str_replace(., "Lib2_", ""))
  
  write.csv(assigned_bcs, '0_Long_read_assignment/assigned_bcs.csv')
  return(assigned_bcs)
  
}

bc_assignments <- Long_read_clustering_and_assignment('0_Long_read_assignment/trimmed_alignments.csv',
                                            '0_Long_read_assignment/DADA2_bc_sequences.csv', 
                                            '0_Long_read_assignment/individually_sequenced_bcs.csv')

######### Short read assignment processing and clustering ###########

#This function takes as inputs barcode sequences & Dada2 calculated abundances, clusters barcodes & combines short reads within each cluster
# then assigns the barcodes using '0_Long_read_assignment/assigned_bcs.csv'. 

Cluster_Dada2_and_assign_barcodes <- function(Dada2_abudances_csv, assigned_bcs_csv) {
  
  # Read Dada2 barcode abundances for your specific experiment 
  all_bcs_reads <- read.csv(Dada2_abudances_csv)
  
  # Cluster barcodes based on >90% similarity, using seq_cluster from bioseq package -> each Dada2 runs generate 'ASVs' -> in practice observed this is likely too stringent -> merge barcodes with >90% similarity. 
  barcodes <- bioseq::dna(unique(c(all_bcs_reads$barcode)))
  clustered_bcs <- bioseq::seq_cluster(barcodes, threshold = 0.10) 
  clustered_bcs_df <- data.frame(barcode = barcodes, cluster = clustered_bcs)
  single_replica_clustered_bcs_df <- clustered_bcs_df%>%
    group_by(cluster)%>%
    slice_head(n = 1) %>%
    ungroup()
  #write.csv(single_replica_clustered_bcs_df, '2_Processing_intermediates_csv/single_replica_clustered_bcs_df_new.csv')
  
  # Combine clusters with all barcode dataframe 
  all_bcs_reads <- inner_join(all_bcs_reads, clustered_bcs_df, by='barcode')
  
  # Add subclusters column in data - this column is simply a count of the barcode number within the cluster (i.e. for cluster 1 if there are 30 barcodes, rows will be numbered within the cluster from 1 to 30)
  all_bcs_reads <- all_bcs_reads %>%
    dplyr::select(c(cluster),everything())%>%
    dplyr::group_by(cluster) %>% 
    mutate(subcluster_count = cumsum(barcode != lag(barcode, default = first(barcode)))+1) %>% #add a subcluster column with numbering of barcode's subcluster. 
    ungroup()%>% 
    select(cluster, subcluster_count, everything()) 
  
  # Add total number of short reads to dominant barcode within cluster: 
  all_bcs_reads <- all_bcs_reads %>%
    dplyr::group_by(cluster) %>%
    summarize(across(where(is.numeric), sum, na.rm = TRUE), .groups = 'drop')
  
  # Combine processed clusters DNA barcode sequence 
  all_bc_reads_new <- inner_join(single_replica_clustered_bcs_df, all_bcs_reads, by='cluster')
  
  # Combine processed clusters with barcode assignments from long read sequencing 
  #import csv of long read assignment alignment 
  assigned_bcs <- read.csv(assigned_bcs_csv)%>%
    dplyr::select(-c(X, cluster))
  
  #df with all data values: 
  all_bc_reads_new2 <- list(all_bc_reads_new, assigned_bcs) %>% purrr::reduce(dplyr::full_join, by='barcode')%>% #join all 2 dfs (dada2 read counts & long read assignment)
    filter(!is.na(cluster))%>% #filter only for barcodes appearing in the data for the experiment that was input into the function. 
    select(-c(X,Long_reads_for_assignment, tot_long_rds_sensor, prop_long_rds_barcode, dis)) #removing unecessary long read information
  
  all_bc_reads_new2$bc <- 1:nrow(all_bc_reads_new2) 
  all_bc_reads_new2$sensor <- gsub( ' ', '',all_bc_reads_new2$sensor)
  
  return(all_bc_reads_new2)
}

#example running the function for LB012 experiment. (takes ~30 seconds for small number of barcodes, but longer for higher numbers)
clust_assigned_Dada2_abundances <- Cluster_Dada2_and_assign_barcodes('1_Dada_bc_abundances_csv/LB012_barcode_abundances.csv', '0_Long_read_assignment/assigned_bcs.csv')

#example running the function for other experiments
clust_assigned_Dada2_abundances <- Cluster_Dada2_and_assign_barcodes('1_Reference_csv/ME9_ME4_ME6_LB009_bc_abundances.csv', '0_Long_read_assignment/assigned_bcs.csv')


######### Odds Ratio Analysis (not pooled barcodes) #########

#This function calculates odds ratios given specified positive control sensors (can also be input as 'fabR'), a specified sample_codes.csv, and the assigned barcodes & abundances
# from the Cluster_Dada2_and_assign_barcodes function. 
#Output is an odds ratio dataframe (output is not automatically saved as a csv)

odds_ratio <- function(positive_sensors,sample_codes_csv,assigned_bcs) {
  
  positive_bcs <- clust_assigned_Dada2_abundances[which(clust_assigned_Dada2_abundances$sensor %in% positive_sensors),]$bc
  
  ###Get positive_normalisation reads & calculate mean positive_normalisation reads per sample
  fabR_reads <- assigned_bcs[assigned_bcs$bc %in% positive_bcs,]%>% 
    dplyr::select(where(is.numeric))%>%
    select(-contains('bc_reads'))%>%
    tidyr::pivot_longer(cols=everything(), names_to = 'Sample', values_to = 'fabR_normed_reads')%>%
    group_by(Sample)%>%
    summarise(fabR_normed_reads = sum(fabR_normed_reads, na.rm = T))
  
  #read sample codes csv 
  Sample_codes <- read.csv(sample_codes_csv)
  
  #merge with sample codes, remove normalising barcodes that have fewer than 100 reads in none or spec sample. take geometric mean of these for positive normalisation factor 
  fabR_reads2 <- inner_join(Sample_codes,fabR_reads, by='Sample') %>%
    select(-c(Sample))%>%
    tidyr::pivot_wider(names_from = Antibiotic, values_from = c(fabR_normed_reads))%>%
    filter(none>=100)%>%
    filter(spec>=100)%>%
    group_by(Long_sample_name)%>%
    summarise_at(vars(none, spec),funs(exp(mean(log(.)))))%>%
    tidyr::pivot_longer(!c(Long_sample_name), names_to = 'Antibiotic',values_to = 'fabR_normalisation_factor')
  
  #Get sample data and format
  sample_data <- assigned_bcs%>%
    select(-tidyselect::any_of(c("barcode", "cluster", "subcluster_count")))%>%
    select(-contains('bc_reads'))%>%
    tidyr::pivot_longer(!c(sensor,bc), names_to = 'Sample', values_to = 'reads')
  
  #merge sample data with sample codes and positive normalisation reads, normalise and calculate odds ratio
  sample_data2 <-inner_join(Sample_codes,sample_data, by='Sample')
  Sample_data <- inner_join(sample_data2, fabR_reads2, by=c('Long_sample_name','Antibiotic'))%>%
    dplyr::mutate(fabR_normed_reads =reads/fabR_normalisation_factor)%>%
    select(-Sample)%>%
    tidyr::pivot_wider (names_from = Antibiotic, values_from = c(reads, fabR_normalisation_factor, fabR_normed_reads))%>%
    filter(reads_none!=0)%>% #only keep sensors with reads in the none sample
    dplyr::mutate_if(is.numeric, ~tidyr::replace_na(.,0))%>% #set NA values (NA reads) to 0 reads. 
    dplyr::mutate(odds_ratio = fabR_normed_reads_spec/fabR_normed_reads_none)%>% #calculate odds ratio
    dplyr::rename(none = fabR_normed_reads_none,
                  spec = fabR_normed_reads_spec)
  
  return (Sample_data)
}

### Normalising with positive bcs: 
Sample_data_positive_bcs <- odds_ratio (c('E_faecalis_TCS3', 'B_producta_TCS12','B_producta_TCS1'), #process all data & calculate odds ratio using specified positive control sensors
                                        '1_Dada_bc_abundances_csv/All_experiments_sample_codes.csv',
                                        clust_assigned_Dada2_abundances)

######### Sample QC Sample Exclusion Criteria (Spread exclusion, Rsquare and FC based) ###########

#This function performs a number of QC checks & filters normalised & odds ratio calculated data based parameters that can be specified:  
# Rsqare_threshold : Rsquare of the fit of positive controls
# FC_threshold : Fold change difference between average OR of positive and negative control sensors. 
# Upper_slope_threshold : maximum slope of fit of positive controls (slope should be ~1)
# Lower_slope_threshold : minimum slope of fit of positive controls (slope should be ~1)
# Fractional threshold : minimum fraction of total 'no spec' reads for a barcode to be included
# min_none_reads : minimum number of reads for a barcode in 'no spec' to be included 
# positive sensors : sensors used as positive controls (this is for plotting QC plots)
# negative sensors : sensors used as negative controls (this is for plotting QC plots)
# sample_codes_csv : specified sample_codes.csv

QC_data <- function(experiment_name,
                    Normalised_Sample_data_input, 
                    Dada2_data_input,
                    Rsquare_threshold, 
                    FC_threshold, 
                    Upper_slope_threshold, 
                    Lower_slope_threshold, 
                    Fractional_threshold, 
                    min_none_reads,
                    positive_sensors,
                    negative_sensors,
                    sample_codes_csv) {
  
  #First filter based on fraction of reads per bc for each sample -> to set threshold based on dilution
  #Get total number of reads per sample
  Total_n_reads_per_sample <- Normalised_Sample_data_input %>%
    group_by(Long_sample_name)%>%
    summarise_at(vars(reads_none, reads_spec),funs(sum(., na.rm=TRUE)))
  
  Sample_data_positive_bcs_no_QC <- Normalised_Sample_data_input
  
  Sample_data <- inner_join(Total_n_reads_per_sample, Normalised_Sample_data_input, by='Long_sample_name') %>%
    mutate(frac_reads_none = reads_none.y/reads_none.x)%>%
    mutate(frac_reads_spec = reads_spec.y/reads_spec.x)%>%
    filter(frac_reads_none >=Fractional_threshold)%>% #filter so that the fraction of library of the barcode isn't less than 1/5000
    filter(reads_none.y>=min_none_reads)%>% #filter to have at least 25 reads in the none sample. 
    mutate(reads_none = reads_none.y)%>%
    mutate(reads_spec = reads_spec.y)%>%
    select(-c(reads_none.y,reads_none.x, reads_spec.y,reads_spec.x,frac_reads_none,frac_reads_spec))
  
  Samples <- unique(Sample_data[which(Sample_data$Experiment != 'ME4'),]$Long_sample_name) 
  R_square_list <- vector("list", length = length(Samples))
  N_barcodes_list <- vector("list", length = length(Samples))
  N_sensors_list <- vector("list", length = length(Samples))
  FC_list <- vector("list", length = length(Samples))
  Slopes_list <- vector("list", length = length(Samples))
  
  output_dir<-paste0('3_plots/QC_plots_',experiment_name,'/')
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  
  #these are used for plotting QC plots 
  negative_bcs <- Normalised_Sample_data_input[which(Normalised_Sample_data_input$sensor %in% negative_sensors),]$bc
  positive_bcs <- Normalised_Sample_data_input[which(Normalised_Sample_data_input$sensor %in% positive_sensors),]$bc
  
  for (i in Samples) {
    
    #calculating difference in fold change between positive and negative sensors for that sample
    positive_data <- Sample_data[which(Sample_data$Long_sample_name == i & Sample_data$bc %in% positive_bcs),]
    positive_data["odds_ratio"][positive_data["odds_ratio"] <= 0.001] <- 0.001
    negative_data <- Sample_data[which(Sample_data$Long_sample_name == i & Sample_data$bc %in% negative_bcs),]
    negative_data["odds_ratio"][negative_data["odds_ratio"] <= 0.001] <- 0.001
    med_neg_or <- negative_data%>%
      group_by(Long_sample_name)%>%
      summarise(odds_ratio_neg = exp(mean(log(odds_ratio),na.rm=T)))
    med_pos_or <- positive_data%>%
      group_by(Long_sample_name)%>%
      summarise(odds_ratio_pos = mean(odds_ratio, na.rm=T))
    fc_neg_to_pos <-  med_pos_or$odds_ratio_pos[1]/med_neg_or$odds_ratio_neg[1]
    FC_list[[i]] <- fc_neg_to_pos
    
    #calculating total number of recovered barcodes for the sample & total number of recovered sensors 
    recovered_bcs <- length(unique(Sample_data[which(Sample_data$Long_sample_name == i),]$bc))
    recovered_sensors <- length(unique(Sample_data[which(Sample_data$Long_sample_name == i),]$sensor))
    N_barcodes_list[[i]] <- recovered_bcs
    N_sensors_list[[i]] <- recovered_sensors
    
    spread_plot_data <- Sample_data[which(Sample_data$Long_sample_name == i & (Sample_data$bc %in% positive_bcs)),]
    
    linear_reg<-lm(log10(spec) ~ log10(none), data=spread_plot_data)
    R_square<-summary(linear_reg)$r.squared
    R_square_list[[i]] <- R_square
    slope <- linear_reg$coefficients[[2]]
    Slopes_list[[i]] <- slope
    
    Sample_data["none"][Sample_data["none"] <= 0.0001] <- 0.0001
    min_line <- min(Sample_data[which(Sample_data$Long_sample_name == i),]$none) 
    Sample_data["spec"][Sample_data["spec"] <= 0.0001] <- 0.0001
    
    ###diagonal plot
    spread_plot <- ggplot(data = spread_plot_data, mapping=aes(x=none,y=spec), size=4)+
      geom_point(Sample_data_positive_bcs_no_QC[which(Sample_data_positive_bcs_no_QC$Long_sample_name == i),], mapping=aes(x=none,y=spec), colour='grey',alpha=0.3,size=1, shape=16)+
      geom_point(Sample_data[which(Sample_data$Long_sample_name == i),], mapping=aes(x=none,y=spec), colour='grey21',alpha=0.3,size=1, shape=16)+
      geom_point(spread_plot_data, mapping=aes(x=none,y=spec), colour='blue', size=2)+
      geom_point(Sample_data[which(Sample_data$Long_sample_name == i & Sample_data$bc %in% negative_bcs),], mapping=aes(x=none,y=spec), colour='red',size=2)+
      ggpmisc::stat_poly_line()+
      ggpmisc::stat_poly_eq()+
      annotate(geom='text',x=0.0001, y=100,hjust=0,label=paste('Barcodes:',recovered_bcs), size=2.5)+
      annotate(geom='text',x=0.0001, y=30,hjust=0,label=paste('Sensors:',recovered_sensors), size=2.5)+
      annotate(geom='text',x=0.0001, y=5,hjust=0,label=paste('FC:',round(fc_neg_to_pos,digits=1)), size=2.5)+
      scale_y_log10( labels = function(none) ifelse(none == 0.0001, "0", format(none, scientific = FALSE)), limits=c(0.0001,1000), breaks=c(0.0001,0.001,0.01,0.1,1,10,100,1000))+
      scale_x_log10( labels = function(spec) ifelse(spec == 0.0001, "0", format(spec, scientific = FALSE)), limits=c(0.0001,1000), breaks=c(0.0001,0.001,0.01,0.1,1,10,100,1000))+
      #scale_y_log10( labels = function(none) ifelse(none == 0.001, "0", format(none, scientific = FALSE)), limits=c(0.000001,1000))+
      #scale_x_log10( labels = function(spec) ifelse(spec == 0.001, "0", format(spec, scientific = FALSE)), limits=c(0.000001,1000))+
      geom_vline(xintercept = min_line, linetype = 'dashed')+
      ggtitle(gsub('_',' ',i))+
      xlab('-SPEC Normed Reads')+
      ylab('+SPEC Normed Reads')+
      theme_bw()+
      geom_abline(slope=1, linetype = 'dashed')+
      labs(color = NULL)+
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill='white'),
            text = element_text(size = 20),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(family = "Helvetica", size = 8, face = "bold"),
            axis.title = element_text(family = "Helvetica", size = 8),
            axis.text = element_text(family = "Helvetica", size = 8),
            legend.text = element_text(family = "Helvetica", size = 8),
            legend.title = element_text(family = "Helvetica", size = 8),
            legend.position ='none')
    spread_plot
    ggsave(spread_plot, file=paste('3_plots/QC_plots_',experiment_name,'/',i,'.jpg', sep=''),
           width = (210-30)/2, height = (297-30)/3, units = 'mm', dpi=300)
  }
  
  R_square_list <- as.data.frame(do.call(rbind, R_square_list))
  N_barcodes_list <-as.data.frame(do.call(rbind, N_barcodes_list))
  N_sensors_list <- as.data.frame(do.call(rbind, N_sensors_list ))
  FC_sensors_list <-as.data.frame(do.call(rbind, FC_list ))
  Slopes_list <- as.data.frame(do.call(rbind, Slopes_list ))
  
  QC_df <- cbind(R_square_list,N_barcodes_list,N_sensors_list,FC_sensors_list,Slopes_list)
  names(QC_df)<-c('Rsquare','N_barcodes','N_sensors', 'FC', 'slope')
  total_reads_per_sample <- as.data.frame(colSums(Filter(is.numeric, Dada2_data_input[, -which(names(Dada2_data_input) %in% c("bc","total_bc_reads",'X..Sequences.y', 'cluster','subcluster_count'))])))
  names(total_reads_per_sample)<-c('reads')
  total_reads_per_sample$Sample <- rownames(total_reads_per_sample)
  total_reads_per_sample <- inner_join(total_reads_per_sample, read.csv(sample_codes_csv), by=c('Sample'))%>%
    select(c(reads, Long_sample_name, Antibiotic, Experiment))%>%
    tidyr::pivot_wider(names_from = Antibiotic, values_from = c(reads))
  QC_df$Long_sample_name <- rownames(QC_df)
  QC_df <- inner_join(QC_df, total_reads_per_sample, by=c('Long_sample_name'))
  
  write.csv(QC_df,paste0(experiment_name, '_QC_df_details.csv'))
  
  QC_passing_samples <- QC_df[which((QC_df$Rsquare>=Rsquare_threshold)&(QC_df$FC>=FC_threshold)&(QC_df$slope<=Upper_slope_threshold)&(QC_df$slope>=Lower_slope_threshold)),]$Long_sample_name
  
  QC_passing_data<- Normalised_Sample_data_input[which(Normalised_Sample_data_input$Long_sample_name %in% QC_passing_samples),]
  write.csv(QC_passing_data, paste0(experiment_name,'_QC_passing_data_R',Rsquare_threshold,'_FC', FC_threshold,'.csv'))
  
  return(QC_passing_data)
}

# example running the function 
QC_passing_data <- QC_data('LB0012',
                           Sample_data_positive_bcs, 
                           clust_assigned_Dada2_abundances,
                           0.6,30,1.5,0.5,1/5000,25,
                           c('E_faecalis_TCS3','B_producta_TCS1','B_producta_TCS12'), 
                           c('L_gasseri_TCS1'),
                           '1_Dada_bc_abundances_csv/All_experiments_sample_codes.csv')

write.csv(QC_passing_data, 'LB0012_QC_passing_data.csv')

######### For samples that require pooling all Lib2 barcodes before calculation (with modified QC function) ############

pooled_odds_ratio <- function(positive_sensors, input_data, sample_codes_csv) {
  #group by sensor and pool all reads
  pooled_bc_assigned_bcs <- input_data%>%
    group_by(sensor)%>%
    summarise(across(where(is.numeric), sum))%>%
    select(-c(bc, cluster, subcluster_count))
  
  #assign new 'bc' numbers -> these aren't real barcodes (as all reads from all barcodes have been pooled) but allows the code to use the same odds_ratio function
  pooled_bc_assigned_bcs$bc <- 1:nrow(pooled_bc_assigned_bcs) 
  positive_bcs <- pooled_bc_assigned_bcs[which(pooled_bc_assigned_bcs$sensor %in% positive_sensors),]$bc

  ###Get positive_normalisation reads & calculate mean positive_normalisation reads per sample
  fabR_reads <- pooled_bc_assigned_bcs[pooled_bc_assigned_bcs$bc %in% positive_bcs,]%>% 
    dplyr::select(where(is.numeric))%>%
    select(-contains('bc_reads'))%>%
    tidyr::pivot_longer(cols=everything(), names_to = 'Sample', values_to = 'fabR_normed_reads')%>%
    group_by(Sample)%>%
    summarise(fabR_normed_reads = sum(fabR_normed_reads, na.rm = T))
  
  #read sample codes csv 
  Sample_codes <- read.csv(sample_codes_csv)
  
  #merge with sample codes, remove normalising barcodes that have fewer than 100 reads in none or spec sample. take geometric mean of these for positive normalisation factor 
  fabR_reads2 <- inner_join(Sample_codes,fabR_reads, by='Sample') %>%
    select(-c(Sample))%>%
    tidyr::pivot_wider(names_from = Antibiotic, values_from = c(fabR_normed_reads))%>%
    filter(none>=100)%>%
    filter(spec>=100)%>%
    group_by(Long_sample_name)%>%
    summarise_at(vars(none, spec),funs(exp(mean(log(.)))))%>%
    tidyr::pivot_longer(!c(Long_sample_name), names_to = 'Antibiotic',values_to = 'fabR_normalisation_factor')
  
  #Get sample data and format
  sample_data <- pooled_bc_assigned_bcs %>%
    select(-tidyselect::any_of(c("barcode", "cluster", "subcluster_count")))%>%
    select(-contains('bc_reads'))%>%
    tidyr::pivot_longer(!c(sensor,bc), names_to = 'Sample', values_to = 'reads')
  
  #merge sample data with sample codes and positive normalisation reads, normalise and calculate odds ratio
  sample_data2 <-inner_join(Sample_codes,sample_data, by='Sample')
  Sample_data <- inner_join(sample_data2, fabR_reads2, by=c('Long_sample_name','Antibiotic'))%>%
    dplyr::mutate(fabR_normed_reads =reads/fabR_normalisation_factor)%>%
    select(-Sample)%>%
    tidyr::pivot_wider (names_from = Antibiotic, values_from = c(reads, fabR_normalisation_factor, fabR_normed_reads))%>%
    filter(reads_none!=0)%>% #only keep sensors with reads in the none sample
    dplyr::mutate_if(is.numeric, ~tidyr::replace_na(.,0))%>% #set NA values (NA reads) to 0 reads. 
    dplyr::mutate(odds_ratio = fabR_normed_reads_spec/fabR_normed_reads_none)%>% #calculate odds ratio
    dplyr::rename(none = fabR_normed_reads_none,
                  spec = fabR_normed_reads_spec)

}

Sample_data_pooled_bcs <- pooled_odds_ratio(c('E_faecalis_TCS3','B_producta_TCS1','B_producta_TCS12'), 
                                            clust_assigned_Dada2_abundances,
                                            '1_Dada_bc_abundances_csv/All_experiments_sample_codes.csv')

QC_data_pooled_bcs <- function(experiment_name,
                               Sample_data_input,
                               Dada2_data_input,
                               Rsquare_threshold, 
                               FC_threshold, 
                               Upper_slope_threshold, 
                               Lower_slope_threshold, 
                               Fractional_threshold, 
                               min_none_reads, 
                               positive_sensors, 
                               negative_sensors, 
                               sample_codes_csv) {
  
  negative_bcs_ME9 <- Sample_data_input[which(Sample_data_input$sensor %in% negative_sensors),]$bc
  positive_bcs_ME9 <- Sample_data_input[which(Sample_data_input$sensor %in% positive_sensors),]$bc
  
  #getting fraction of reads per bc for each sample -> to set threshold based on dilution
  Total_n_reads_per_sample_pooled_bcs <- Sample_data_input %>%
    group_by(Long_sample_name)%>%
    summarise_at(vars(reads_none, reads_spec),funs(sum(., na.rm=TRUE)))
  
  Sample_data_pooled_bcs_no_QC <- Sample_data_input
  
  Sample_data_input <- inner_join(Total_n_reads_per_sample_pooled_bcs, Sample_data_input, by='Long_sample_name') %>%
    mutate(frac_reads_none = reads_none.y/reads_none.x)%>%
    mutate(frac_reads_spec = reads_spec.y/reads_spec.x)%>%
    filter(frac_reads_none >=Fractional_threshold)%>%
    filter(reads_none.y>=min_none_reads)%>%
    mutate(reads_none = reads_none.y)%>%
    mutate(reads_spec = reads_spec.y)%>%
    select(-c(reads_none.y,reads_none.x, reads_spec.y,reads_spec.x,frac_reads_none,frac_reads_spec))
  
  Sample_data <- Sample_data_input
  Samples <- unique(Sample_data[which(Sample_data$Experiment != 'ME4'),]$Long_sample_name) 
  R_square_list <- vector("list", length = length(Samples))
  N_barcodes_list <- vector("list", length = length(Samples))
  N_sensors_list <- vector("list", length = length(Samples))
  FC_list <- vector("list", length = length(Samples))
  Slopes_list <- vector("list", length = length(Samples))
  
  output_dir<-paste0('3_plots/QC_plots_pooled_bcs_',experiment_name,'/')
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  
  for (i in Samples) {
    
    #calculating difference in fold change between positive and negative sensors for that sample
    positive_data <- Sample_data[which(Sample_data$Long_sample_name == i & Sample_data$bc %in% positive_bcs_ME9),]
    positive_data["odds_ratio"][positive_data["odds_ratio"] <= 0.001] <- 0.001
    negative_data <- Sample_data[which(Sample_data$Long_sample_name == i & Sample_data$bc %in% negative_bcs_ME9),]
    negative_data["odds_ratio"][negative_data["odds_ratio"] <= 0.001] <- 0.001
    med_neg_or <- negative_data%>%
      group_by(Long_sample_name)%>%
      summarise(odds_ratio_neg = exp(mean(log(odds_ratio),na.rm=T)))
    med_pos_or <- positive_data%>%
      group_by(Long_sample_name)%>%
      summarise(odds_ratio_pos = mean(odds_ratio, na.rm=T))
    fc_neg_to_pos <-  med_pos_or$odds_ratio_pos[1]/med_neg_or$odds_ratio_neg[1]
    FC_list[[i]] <- fc_neg_to_pos
    
    #calculating total number of recovered barcodes for the sample & total number of recovered sensors 
    recovered_bcs <- length(unique(Sample_data[which(Sample_data$Long_sample_name == i),]$bc))
    recovered_sensors <- length(unique(Sample_data[which(Sample_data$Long_sample_name == i),]$sensor))
    N_barcodes_list[[i]] <- recovered_bcs
    N_sensors_list[[i]] <- recovered_sensors
    
    spread_plot_data <- Sample_data[which(Sample_data$Long_sample_name == i & (Sample_data$bc %in% positive_bcs_ME9)),]
    
    linear_reg<-lm(log10(spec) ~ log10(none), data=spread_plot_data)
    R_square<-summary(linear_reg)$r.squared
    R_square_list[[i]] <- R_square
    slope <- linear_reg$coefficients[[2]]
    Slopes_list[[i]] <- slope
    
    Sample_data["none"][Sample_data["none"] <= 0.0001] <- 0.0001
    min_line <- min(Sample_data[which(Sample_data$Long_sample_name == i),]$none) 
    Sample_data["spec"][Sample_data["spec"] <= 0.0001] <- 0.0001
    
    ###diagonal plot
    spread_plot <- ggplot(data = spread_plot_data, mapping=aes(x=none,y=spec), size=4)+
      geom_point(Sample_data_pooled_bcs_no_QC[which(Sample_data_pooled_bcs_no_QC$Long_sample_name == i),], mapping=aes(x=none,y=spec), colour='grey',alpha=0.3,size=1, shape=16)+
      geom_point(Sample_data[which(Sample_data$Long_sample_name == i),], mapping=aes(x=none,y=spec), colour='grey11',alpha=0.4,size=1, shape=16)+
      geom_point(spread_plot_data, mapping=aes(x=none,y=spec), colour='blue', size=2)+
      geom_point(Sample_data[which(Sample_data$Long_sample_name == i & Sample_data$bc %in% negative_bcs_ME9),], mapping=aes(x=none,y=spec), colour='red',size=2)+
      ggpmisc::stat_poly_line()+
      ggpmisc::stat_poly_eq()+
      annotate(geom='text',x=0.0001, y=100,hjust=0,label=paste('Barcodes:',recovered_bcs), size=2.5)+
      annotate(geom='text',x=0.0001, y=30,hjust=0,label=paste('Sensors:',recovered_sensors), size=2.5)+
      annotate(geom='text',x=0.0001, y=5,hjust=0,label=paste('FC:',round(fc_neg_to_pos,digits=1)), size=2.5)+
      scale_y_log10( labels = function(none) ifelse(none == 0.0001, "0", format(none, scientific = FALSE)), limits=c(0.0001,1000), breaks=c(0.0001,0.001,0.01,0.1,1,10,100,1000))+
      scale_x_log10( labels = function(spec) ifelse(spec == 0.0001, "0", format(spec, scientific = FALSE)), limits=c(0.0001,1000), breaks=c(0.0001,0.001,0.01,0.1,1,10,100,1000))+
      #scale_y_log10( labels = function(none) ifelse(none == 0.001, "0", format(none, scientific = FALSE)), limits=c(0.000001,1000))+
      #scale_x_log10( labels = function(spec) ifelse(spec == 0.001, "0", format(spec, scientific = FALSE)), limits=c(0.000001,1000))+
      geom_vline(xintercept = min_line, linetype = 'dashed')+
      ggtitle(gsub('_',' ',i))+
      xlab('-SPEC Normed Reads')+
      ylab('+SPEC Normed Reads')+
      theme_bw()+
      geom_abline(slope=1, linetype = 'dashed')+
      labs(color = NULL)+
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill='white'),
            text = element_text(size = 20),
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(family = "Helvetica", size = 8, face = "bold"),
            axis.title = element_text(family = "Helvetica", size = 8),
            axis.text = element_text(family = "Helvetica", size = 8),
            legend.text = element_text(family = "Helvetica", size = 8),
            legend.title = element_text(family = "Helvetica", size = 8),
            legend.position ='none')
    spread_plot
    ggsave(spread_plot, file=paste('3_plots/QC_plots_pooled_bcs_',experiment_name,'/',i,'.jpg', sep=''),
           width = (210-30)/2, height = (297-30)/3, units = 'mm', dpi=300)
  }
  
  R_square_list <- as.data.frame(do.call(rbind, R_square_list))
  N_barcodes_list <-as.data.frame(do.call(rbind, N_barcodes_list))
  N_sensors_list <- as.data.frame(do.call(rbind, N_sensors_list ))
  FC_sensors_list <-as.data.frame(do.call(rbind, FC_list ))
  Slopes_list <- as.data.frame(do.call(rbind, Slopes_list ))
  
  QC_df <- cbind(R_square_list,N_barcodes_list,N_sensors_list,FC_sensors_list,Slopes_list)
  names(QC_df)<-c('Rsquare','N_barcodes','N_sensors', 'FC', 'slope')
  total_reads_per_sample <- as.data.frame(colSums(Filter(is.numeric, Dada2_data_input[, -which(names(Dada2_data_input) %in% c("bc","total_bc_reads",'X..Sequences.y'))])))
  names(total_reads_per_sample)<-c('reads')
  total_reads_per_sample$Sample <- rownames(total_reads_per_sample)
  total_reads_per_sample <- inner_join(total_reads_per_sample, read.csv(sample_codes_csv), by=c('Sample'))%>%
    select(c(reads, Long_sample_name, Antibiotic, Experiment))%>%
    tidyr::pivot_wider(names_from = Antibiotic, values_from = c(reads))
  QC_df$Long_sample_name <- rownames(QC_df)
  QC_df <- inner_join(QC_df, total_reads_per_sample, by=c('Long_sample_name'))
  
  write.csv(QC_df,paste0(experiment_name,'_QC_df_details_pooled_bcs.csv'))
  
  QC_passing_samples <- QC_df[which((QC_df$Rsquare>=Rsquare_threshold)&(QC_df$FC>=FC_threshold)&(QC_df$slope<=Upper_slope_threshold)&(QC_df$slope>=Lower_slope_threshold)),]$Long_sample_name
  
  QC_passing_data<- Sample_data_input[which(Sample_data_input$Long_sample_name %in% QC_passing_samples),]
  write.csv(QC_passing_data, paste0(experiment_name,'_QC_passing_data_pooled_bcs_R',Rsquare_threshold,'_FC', FC_threshold,'.csv'))
  
  return(QC_passing_data)
}

QC_passing_data_pooled_bcs <- QC_data_pooled_bcs('LB0012',
                                                 Sample_data_pooled_bcs,
                                                 clust_assigned_Dada2_abundances,
                                                 0.6, 2, 2, 0.5,1/5000,25, 
                                                 c('E_faecalis_TCS3','B_producta_TCS1','B_producta_TCS12'), 
                                                 c('L_gasseri_TCS1'), 
                                                 '1_Dada_bc_abundances_csv/All_experiments_sample_codes.csv')

write.csv(QC_passing_data_pooled_bcs, 'QC_passing_data_pooled_bcs.csv')

######### Calculate fractional change: ######### 

# these functions calculate fractional change in OR (FOR) for each barcode and mean of all barcodes (for a specific sensor) and save the results as csvs. 
get_fc <- function(sampled_df, positive_bcs, negative_bcs, sensor_name) {
  
  positive_data <- sampled_df[which(sampled_df$bc %in% positive_bcs),]
  negative_data <- sampled_df[which(sampled_df$bc %in% negative_bcs),]
  
  #calculating split & odds ratio as distance between negative & positive mean of sensor medians. 
  med_neg_or <- negative_data%>%
    group_by(Long_sample_name)%>%
    summarise(odds_ratio_neg = exp(mean(log(odds_ratio),na.rm=T)))
  
  med_pos_or <- positive_data%>%
    group_by(Long_sample_name)%>%
    summarise(odds_ratio_pos = exp(mean(log(odds_ratio), na.rm=T)))
  
  sampled_df_fc <- inner_join(inner_join(sampled_df, med_neg_or, by=c('Long_sample_name')),med_pos_or, by=c('Long_sample_name'))%>%
    mutate(neg_pos_diff = odds_ratio_pos - odds_ratio_neg)%>%
    mutate(frac_sensor_OR = (odds_ratio - odds_ratio_neg)/neg_pos_diff)%>%
    filter(sensor == sensor_name)
  
  return(sampled_df_fc)
} #function to calculate fractional change
calculate_FOR <- function(experiment_name, input_data, positive_sensors, negative_sensors){
  
  positive_bcs <- input_data[which(input_data$sensor %in% positive_sensors),]$bc
  negative_bcs <- input_data[which(input_data$sensor %in% negative_sensors),]$bc
  
  List_of_sample_names <- unique(input_data$Long_sample_name)
  list_of_sensors <- unique(input_data$sensor)
  
  fc_distributions <- vector("list", length = length(List_of_sample_names)) #create vector in which to store fractional change values 
  
  #fractional change FOR is calculated iteratively for each sample in a loop 
  for(Chosen_sample in List_of_sample_names) {
    
    #get data for specific sample
    Sample_df <- input_data%>%
      filter(Long_sample_name == Chosen_sample)
    
    sensor_fcs <- vector("list")
    for (sensor_name in list_of_sensors) {
      ##get the actual fold change for all data points (will be plotted as red line)
      #test <- get_fc(Sample_df[which(Sample_df$sensor==sensor_name|Sample_df$bc %in% c(positive_bcs, negative_bcs)),], positive_bcs, negative_bcs, sensor_name)
      sensor_fcs[[sensor_name]] <- get_fc(Sample_df[which(Sample_df$sensor==sensor_name|Sample_df$bc %in% c(positive_bcs, negative_bcs)),], positive_bcs, negative_bcs, sensor_name)
     
    }
    
    #add sample fc distributions (of all sensors) to the vector of fold change distributions (all samples, all sensors)
    fc_distributions[[Chosen_sample]] <- sensor_fcs
  } 
  
  #combine FOR values into a data frame 
  fc_df <- imap_dfr(fc_distributions, function(sublist, top_index) {
    imap_dfr(sublist, function(df, sub_index) {
      df 
    })
  }) 
  
  #write csv with FOR for every barcode 
  write.csv(fc_df, paste(experiment_name,'_Fold_Change_QC_passing_data_all_barcodes.csv', sep=''))
  
  # repeat but with calculating MEAN FOR across multiple barcodes
  fc_distributions <- vector("list", length = length(List_of_sample_names))
  for(Chosen_sample in List_of_sample_names) {
    
    #get data for specific sample
    Sample_df <- input_data%>%
      filter(Long_sample_name == Chosen_sample)
    
    sensor_fcs <- vector("list")
    for (sensor_name in list_of_sensors) {
      sensor_fc <- get_fc(Sample_df[which(Sample_df$sensor==sensor_name|Sample_df$bc %in% c(positive_bcs, negative_bcs)),], positive_bcs, negative_bcs, sensor_name)
      sensor_fcs[[sensor_name]] <- mean(sensor_fc$frac_sensor_OR)
    }
    fc_distributions[[Chosen_sample]] <- sensor_fcs
  }
  
  Mean_Fc_distributions_df <- as.data.frame(do.call(cbind, fc_distributions))
  Mean_Fc_distributions_df$sensor <-rownames(Mean_Fc_distributions_df)
  Mean_Fc_distributions_df <- apply(Mean_Fc_distributions_df,2,as.character)
  write.csv(Mean_Fc_distributions_df, paste(experiment_name,'_Fold_Change_QC_passing_data_all_mean.csv',sep=''))
  
}

# example running the function 
#reading in QC passing data (from previous code block)
QC_passing_data <- read.csv('LB0012_QC_passing_data.csv')%>%
  select(-c(X))

calculate_FOR('LB0012',
              QC_passing_data,
              c('E_faecalis_TCS3','B_producta_TCS1','B_producta_TCS12'), 
              c('L_gasseri_TCS1'))

calculate_FOR('LB0012_pooled',
              QC_passing_data_pooled_bcs,
              c('E_faecalis_TCS3','B_producta_TCS1','B_producta_TCS12'), 
              c('L_gasseri_TCS1'))

######### ME4 N barcodes #######

#this is extra analysis specifically for the colonisation experiment ME4
#extra functions and libraries just for this block of code
nonzero <- function(x) sum(x != 0)
`%!in%` <- function(a,b) ! a %in% b
library(RColorBrewer)

#read in data 
assigned_bcs <- read.csv('2_Processing_intermediates_csv/assigned_bcs_and_counts.csv')%>%
  select(-c(X))
sample_codes <- read.csv('1_Reference_csv/All_experiments_sample_codes.csv')

#remove any fabR reads. 
no_fabR_data <- assigned_bcs%>%
  filter(sensor %!in% c('fabRbc1', 'fabRbc2', 'fabRbc3','fabRbc4'))

#calculate number of barcodes per sample
Number_of_barcodes_per_sample <- plyr::numcolwise(nonzero)(no_fabR_data) %>%
  select(-bc)%>%
  tidyr::pivot_longer(cols = everything(),names_to = 'Sample', values_to='N_barcodes')

#calculate number of sensors per sample
Number_of_sensors_per_sample <- no_fabR_data %>%
  tidyr::pivot_longer(!c(bc,sensor), names_to = 'Sample', values_to = 'reads')%>%
  filter(reads!=0)%>%
  group_by(Sample)%>%
  summarise(N_sensors=n_distinct(sensor))

sensor_bcs_counts <- full_join(sample_codes, full_join(Number_of_barcodes_per_sample, Number_of_sensors_per_sample, by='Sample'), by='Sample') 

write.csv(sensor_bcs_counts, '2_Processing_intermediates_csv/sensor_bcs_counts.csv')

#filter for ME4 data & format data for plotting 
no_fabR_data2 <- no_fabR_data%>%
  tidyr::pivot_longer(cols=!c(sensor, bc), names_to = 'Sample', values_to = 'reads')
no_fabR_data3 <- inner_join(sample_codes,no_fabR_data2, by='Sample') %>%
  filter(Experiment =='ME4')%>%
  filter(reads != 0)
reads_per_sample <- no_fabR_data3 %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(total_reads = sum(reads)) %>%
  dplyr::slice(1)
no_fabR_data3 <- merge(no_fabR_data3,y=reads_per_sample[,c('Sample','total_reads')], by='Sample')
no_fabR_data3$normed_Reads <- no_fabR_data3$reads/no_fabR_data3$total_reads

# Generate a palette of random colors
num_colors <- length(unique(no_fabR_data3$bc))
random_palette <- sample(colors(), num_colors)

#get days to plot & put in order 
day_order <- c(2, 4, 7, 12)
no_fabR_data3$Day <- factor(no_fabR_data3$Day, levels = day_order)
#save csv for future plotting
write.csv(no_fabR_data3, 'ME4_fractional_data.csv')

#create barcode fraction plot 
barcode_fraction_plot <- ggplot((no_fabR_data3[which(no_fabR_data3$Group=='5 mg strep'),]), aes(y=normed_Reads, x=Day, fill = as.factor(bc)))+
  geom_bar(stat = "identity")+
  facet_grid(rows=vars(Group), cols=vars(Mouse))+
  theme_bw()+
  ylab('Fraction of reads per bc')+
  xlab('Day')+
  scale_fill_manual(values = random_palette) +
  #scale_fill_viridis_b()+
  #scale_fill_brewer(palette = 'Dark2')+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.4),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA),
        strip.text.y = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) 
barcode_fraction_plot
ggsave(barcode_fraction_plot, file='barcode_fraction_plot_group1.pdf', width = 7, height = 3)

ME4_data_wider <- no_fabR_data3%>%
  select(c(Sample, sensor, bc, normed_Reads))%>%
  tidyr::pivot_wider(names_from = Sample, values_from = normed_Reads)

write.csv(ME4_data_wider, 'ME4_fractional_data_wider.csv')

