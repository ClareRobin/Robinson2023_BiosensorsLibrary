
######### Loading required packages ###########
library(extrafont)

# Load fonts for plotting
font_import(paths = NULL, recursive = TRUE, prompt = TRUE,pattern = NULL)
loadfonts()
'%!in%' <- function(x,y)!('%in%'(x,y))

library(dplyr)
library(ggplot2)

#For R Studio only: set workind directory to location of script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

######### Merging dada2 barcode count & barcode sensor assignments ###########

#Merge dada2 barcode count and barcode sensor assignments: 
merge_Dada2_and_Sensor_assignments <- function(com_bcs_csv, assigned_bcs_csv) {
  
  #import csv of dada2 barcodes & barcode read counts
  com_bcs <- read.csv(com_bcs_csv)%>%
    dplyr::select(-X)
  names(com_bcs)[names(com_bcs) == 'cluster'] <- 'bc'
  
  #import csv of barcode sensor assignments
  assigned_bcs <- read.csv(assigned_bcs_csv)%>%
    dplyr::select(-c(X, X..Sequences.y))
  names(assigned_bcs)[names(assigned_bcs) == 'barcode'] <- 'bc' #rename barcode column to bc for later merging with com_bcs
  
  #rename control sensors from plasmid names to control sensors names (for ex. plasmid CR21 pCR21 -> yeaR23)
  control_name_pattern <- c('pCR24 PCR Product','pCR23  PCR Product','pCR22  PCR Product','pCR21 PCR Product','pCR18  PCR Product','pCR17 PCR Product','pCR16  PCR Product','pCR15  PCR Product')
  control_name_replacement <- c('torC23','torC17','torC10','yeaR23','torCWT','hycAWT','ynfE17','ynfE15')
  assigned_bcs$top_by_reads <- stringi::stri_replace_all_regex(assigned_bcs$top_by_reads,
                                                               pattern=control_name_pattern,
                                                               replacement =control_name_replacement,
                                                               vectorize=FALSE)
  assigned_bcs$top_prop <- stringi::stri_replace_all_regex(assigned_bcs$top_prop,
                                                           pattern=control_name_pattern,
                                                           replacement =control_name_replacement,
                                                           vectorize=FALSE)
  
  ######### merge barcode counts & barcode sensor assignments & filter sensors with assignment disagreement #########
  assigned_bcs <- merge(com_bcs, assigned_bcs, by='bc', all=TRUE) %>%
    select(bc,top_by_reads,top_prop, everything()) %>%
    dplyr::mutate(dis = if_else(top_by_reads==top_prop,'agreement','disagreement'))%>% #add column for if assignment is in agreement or disagreement 
    dplyr::select(dis, everything()) %>% #move dis column to left side of df for easier viewing
    dplyr::mutate(dis= ifelse(is.na(dis)&!is.na(sensor),'manually assigned',dis)) #add manually assigned sensors
  
  # Remove assignments that disagree between top_by_reads and top_prop
  unassigned_bcs<-assigned_bcs[which(is.na(assigned_bcs$dis)),]
  disassigned_bcs<-assigned_bcs[which(assigned_bcs$dis=='disagreement'),]
  unassigned_bcs_over_500<-unassigned_bcs%>%
    filter(total_bc_reads >=100)
  assigned_bcs<-assigned_bcs[which(assigned_bcs$dis!='disagreement'),]
  
  # Consolidate sensor assignment into single column
  assigned_bcs[is.na(assigned_bcs$top_by_reads) & !is.na(assigned_bcs$sensor), "top_by_reads"] <- assigned_bcs[is.na(assigned_bcs$top_by_reads) & !is.na(assigned_bcs$sensor), "sensor"]
  assigned_bcs$sensor <- assigned_bcs$top_by_reads
  assigned_bcs <- assigned_bcs %>%
    select(-c(top_by_reads,top_prop,dis))

  return(assigned_bcs)
}
assigned_bcs<-merge_Dada2_and_Sensor_assignments('combined_bcs.csv','assigned_barcodes.csv')
write.csv(assigned_bcs, 'assigned_bcs_with_sequence.csv')

#number of barcodes per sensor
n_barcodes_per_sensor<-data.frame(table(assigned_bcs$sensor))
write.csv(n_barcodes_per_sensor, 'n_barcodes_per_sensor.csv')

######### Odds ratio & new bcs positive controls ###########

#Odds ratio calculation function: 
odds_ratio <- function(positive_bcs,sample_codes_csv,assigned_bcs) {

  ###Get positive_normalisation reads & calculate mean positive_normalisation reads per sample
  fabR_reads <- assigned_bcs[assigned_bcs$bc %in% positive_bcs,] %>% 
    dplyr::select(where(is.numeric)) %>%
    select(-c(total_bc_reads))%>%
    tidyr::pivot_longer(cols=!c(bc), names_to = 'Sample', values_to = 'fabR_normed_reads') 
  
  #read sample codes csv 
  Sample_codes <- read.csv(sample_codes_csv)
  
  #merge with sample codes, remove normalising barcodes that have fewer than 50 reads in none or 25 reads in spec sample. take geometric mean of these for fabR normalisation factor 
  fabR_reads2 <- merge(Sample_codes,fabR_reads, by='Sample') %>%
    select(-c(Sample))%>%
    tidyr::pivot_wider(names_from = Antibiotic, values_from = c(fabR_normed_reads))%>%
    filter(none>=50)%>%
    filter(spec>=25)%>%
    group_by(Long_sample_name)%>%
    summarise_at(vars(none, spec),funs(exp(mean(log(.)))))%>%
    tidyr::pivot_longer(!c(Long_sample_name), names_to = 'Antibiotic',values_to = 'fabR_normalisation_factor')
  
  #Get sample data and format
  sample_data <- assigned_bcs%>%
    select(-c(total_bc_reads))%>%
    tidyr::pivot_longer(!c(bc,sensor,barcode), names_to = 'Sample', values_to = 'reads')
  
  #merge sample data with sample codes and positive normalisation reads, normalise and calculate odds ratio
  sample_data2 <-merge(Sample_codes,sample_data, by='Sample')
  Sample_data <- merge(sample_data2, fabR_reads2, by=c('Long_sample_name','Antibiotic'))%>%
    dplyr::mutate(fabR_normed_reads =reads/fabR_normalisation_factor)%>%
    select(-Sample)%>%
    tidyr::pivot_wider (names_from = Antibiotic, values_from = c(reads, fabR_normalisation_factor, fabR_normed_reads))%>%
    filter(reads_none!=0)%>% #only keep sensors with reads in the none sample
    dplyr::mutate_if(is.numeric, ~tidyr::replace_na(.,0))%>% #set NA values (NA reads) to 0 reads. 
    dplyr::mutate(odds_ratio = fabR_normed_reads_spec/fabR_normed_reads_none)%>% #calculate odds ratio
    rename(none = fabR_normed_reads_none,
           spec = fabR_normed_reads_spec)
  
  return (Sample_data)
}

#Plotting functions
plot_conditions <- function(Sample_data, name, experiment, limits) {
  ######### Diagonal plotting all conditions ############
  list_of_conditions <- unique(Sample_data[which(Sample_data$Experiment == experiment  ),] $Long_sample_name)
  output_dir<-paste('odds_ratio_plots/',experiment,'diagonal_conditions_',name,'/')
  if (!dir.exists(output_dir)) {dir.create(output_dir)}
  
  for (i in list_of_conditions) {
    sensor_df <- Sample_data[which(Sample_data$Long_sample_name == i),] 
    
    sensor_df["odds_ratio"][sensor_df["odds_ratio"] <= 0.001] <- 0.001
    sensor_df["none"][sensor_df["none"] <= 0.001] <- 0.001
    sensor_df["spec"][sensor_df["spec"] <= 0.001] <- 0.001
    
    #diagonal plot
    odds_ratio_plot <- ggplot()+
      geom_point(sensor_df, mapping=aes(x=none,y=spec, colour='All Sensors'), colour='grey',alpha=0.5,size=1)+
      geom_point(sensor_df[which(sensor_df$sensor == 'fabRbc1'|
                                   sensor_df$sensor == 'fabRbc2'|
                                   sensor_df$sensor == 'fabRbc3'|
                                   sensor_df$sensor == 'fabRbc4'),], mapping=aes(x=none,y=spec, colour='test'), colour='black',size=3)+
      geom_point(sensor_df[which(sensor_df$sensor == 'E_faecalis_TCS3'|
                                   sensor_df$sensor == 'B_producta_TCS1'|
                                   sensor_df$sensor == 'B_producta_TCS12'),], mapping=aes(x=none,y=spec, colour=sensor),alpha=1,size=2)+
      scale_y_log10( labels = function(none) ifelse(none == 0.001, "0", format(none, scientific = FALSE)), limits=limits, breaks=c(0.001,0.01,0.1,1,10,100,1000))+
      scale_x_log10( labels = function(spec) ifelse(spec == 0.001, "0", format(spec, scientific = FALSE)), limits=limits, breaks=c(0.001,0.01,0.1,1,10,100,1000))+
      ggtitle(paste(gsub('_',' ',i),'',name,' normalised'))+
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
            legend.position = 'none')
    odds_ratio_plot
    ggsave(odds_ratio_plot, file=paste(output_dir,i,'.pdf', sep=''),
           width = (210-30)/2, height = (297-30)/3, units = 'mm', dpi=300)
  }
  
}

#Calculate odds ratio using fabR barcodes and plot: 
fabR_bcs <- assigned_bcs[which(assigned_bcs$sensor %in% c('fabRbc1','fabRbc2','fabRbc3','fabRbc4')),]$bc #get fabR barcodes
Sample_data_fabR_normalised <- odds_ratio (fabR_bcs, #process all data & calculate odds ratio using fabR barcodes
                           'All_experiments_sample_codes.csv',
                           assigned_bcs)%>%
  filter(reads_none>=10) #filter any barcodes with reads in none sample <10
plot_conditions(Sample_data_fabR_normalised, 'fabR','LB009',c(0.001,10000))

#Highest mean odds ratio sensors in Gavage sample from ME6 (these were cloned as individual strains)
Mean_odds_ratio <- Sample_data_fabR_normalised %>%
  filter(Experiment!='ME4')%>%
  group_by(Long_sample_name,sensor)%>%
  summarise_at(vars(odds_ratio),funs(mean(., na.rm=TRUE)))
mean_sensor_odds_ratio <- Mean_odds_ratio%>%
  filter(Long_sample_name=='Gavage')%>%
  group_by(sensor)%>%
  summarise_at(vars(odds_ratio),funs(mean(., na.rm=TRUE)))%>%
  filter(sensor!= 'fabRbc1' & sensor!= 'fabRbc2'& sensor!= 'fabRbc3'& sensor!= 'fabRbc4')%>%
  slice_max(order_by = odds_ratio, n = 3)

#Normalise using barcodes from these sensors instead of fabR
positive_bcs<-assigned_bcs[which(assigned_bcs$sensor %in% mean_sensor_odds_ratio$sensor),]$bc

Sample_data_positive_bcs <- odds_ratio (positive_bcs, #process all data & calculate odds ratio using positive barcodes
                                           'All_experiments_sample_codes.csv',
                                           assigned_bcs)
plot_conditions(Sample_data_positive_bcs, '3_pos_sensors','LB009',c(0.001,5000))
plot_conditions(Sample_data_positive_bcs, '3_pos_sensors','ME6',c(0.001,10000))
plot_conditions(Sample_data_positive_bcs, '3_pos_sensors','ME1',c(0.001,10000))

#Geometric means odds ratio for heat map: 
get_median_odds_ratio <- function(Sample_data,destination_file) {
  order <- c('ynfE15',
             'S_enterica_TCS_ttrSRRIG',
             'hycAWT',
             'torC10',
             'torC17',
             'torC23',
             'torCWT',
             'ynfE17',
             'yeaR23',
             'B_caccae_TCS4',
             'B_caccae_TCS5',
             'B_frag_TCS6',
             'B_theta_TCS2',
             'E_tarda_TCS4',
             'K_pneumonia_TCS8',
             'R_gnavus_TCS3',
             'S_enterica_TCS1',
             'S_enterica_TCS3',
             'Y_enterocolitica_TCS2',
             'C_rodentium_TCS1',
             'L_gasseri_TCS1',
             'E_tarda_TCS7',
             'C_rodentium_TCS3',
             'B_frag_TCS2',
             'C_rodentium_TCS10',
             'C_rodentium_TCS2',
             'R_gnavus_TCS5',
             'C_rodentium_TCS6',
             'K_pneumonia_TCS5',
             'Y_enterocolitica_TCS1',
             'E_tarda_TCS5',
             'V_cholera_TCS3',
             'C_rodentium_TCS7',
             'K_pneumonia_TCS3',
             'K_pneumoniae_TCS10',
             'B_producta_TCS6',
             'Y_enterocolitica_TCS5',
             'K_pneumonia_TCS4',
             'K_pneumonia_TCS9',
             'C_rodentium_TCS8',
             'E_tarda_TCS6',
             'S_enterica_TCS4',
             'C_rodentium_TCS4',
             'S_enterica_TCS5',
             'B_producta_TCS13',
             'C_rodentium_TCS9',
             'Y_enterocolitica_TCS3',
             'B_producta_TCS16',
             'L_plantarum_TCS2',
             'K_oxytoca_TCS1',
             'C_rodentium_TCS11',
             'C_rodentium_TCS5',
             'V_cholera_TCS4',
             'K_pneumonia_TCS7',
             'E_faecalis_TCS2',
             'E_tarda_TCS2',
             'B_frag_TCS3',
             'C_perfringens_TCS2',
             'V_cholera_TCS2',
             'B_producta_TCS12',
             'B_producta_TCS1',
             'E_faecalis_TCS3',
             'C_diff_TCS5')
  Log_odds_ratio <- Sample_data
  Log_odds_ratio["odds_ratio"][Log_odds_ratio["odds_ratio"] <= 0.001] <- 0.001 #set odds ratio below dilution factor a set minimum for LB009
  Log_odds_ratio["odds_ratio"][Log_odds_ratio["odds_ratio"] >= 1] <- 1
  Log_odds_ratio[which(Log_odds_ratio$Experiment =='ME6'),]["odds_ratio"][Log_odds_ratio[which(Log_odds_ratio$Experiment =='ME6'),]["odds_ratio"] <= 0.01] <- 0.01 #set odds ratio below dilution factor a set minimum for ME6
  Median_odds_ratio <- Log_odds_ratio %>%
    filter(Experiment!='ME4'&Experiment!='ME1')%>%
    group_by(Long_sample_name,sensor)%>%
    dplyr::mutate(median_OR = odds_ratio)%>%
    summarise_at(vars(median_OR),funs(median(., na.rm=TRUE)))
  
  Median_odds_ratio <- Median_odds_ratio%>%
    tidyr::pivot_wider(names_from = Long_sample_name, values_from=median_OR)%>%
    arrange(Gavage)%>%
    arrange(factor(sensor, levels = order))
  write.csv(Median_odds_ratio,destination_file)
  return(Median_odds_ratio)
}
get_geom_mean_odds_ratio <- function(Sample_data, destination_file) {
  order <- c('ynfE15',
             'S_enterica_TCS_ttrSRRIG',
             'hycAWT',
             'torC10',
             'torC17',
             'torC23',
             'torCWT',
             'ynfE17',
             'yeaR23',
             'B_caccae_TCS4',
             'B_caccae_TCS5',
             'B_frag_TCS6',
             'B_theta_TCS2',
             'E_tarda_TCS4',
             'K_pneumonia_TCS8',
             'R_gnavus_TCS3',
             'S_enterica_TCS1',
             'S_enterica_TCS3',
             'Y_enterocolitica_TCS2',
             'C_rodentium_TCS1',
             'L_gasseri_TCS1',
             'E_tarda_TCS7',
             'C_rodentium_TCS3',
             'B_frag_TCS2',
             'C_rodentium_TCS10',
             'C_rodentium_TCS2',
             'R_gnavus_TCS5',
             'C_rodentium_TCS6',
             'K_pneumonia_TCS5',
             'Y_enterocolitica_TCS1',
             'E_tarda_TCS5',
             'V_cholera_TCS3',
             'C_rodentium_TCS7',
             'K_pneumonia_TCS3',
             'K_pneumoniae_TCS10',
             'B_producta_TCS6',
             'Y_enterocolitica_TCS5',
             'K_pneumonia_TCS4',
             'K_pneumonia_TCS9',
             'C_rodentium_TCS8',
             'E_tarda_TCS6',
             'S_enterica_TCS4',
             'C_rodentium_TCS4',
             'S_enterica_TCS5',
             'B_producta_TCS13',
             'C_rodentium_TCS9',
             'Y_enterocolitica_TCS3',
             'B_producta_TCS16',
             'L_plantarum_TCS2',
             'K_oxytoca_TCS1',
             'C_rodentium_TCS11',
             'C_rodentium_TCS5',
             'V_cholera_TCS4',
             'K_pneumonia_TCS7',
             'E_faecalis_TCS2',
             'E_tarda_TCS2',
             'B_frag_TCS3',
             'C_perfringens_TCS2',
             'V_cholera_TCS2',
             'B_producta_TCS12',
             'B_producta_TCS1',
             'E_faecalis_TCS3',
             'C_diff_TCS5')
  Log_odds_ratio <- Sample_data
  Log_odds_ratio["odds_ratio"][Log_odds_ratio["odds_ratio"] <= 0.001] <- 0.001 #set odds ratio below dilution factor a set minimum for LB009
  Log_odds_ratio["odds_ratio"][Log_odds_ratio["odds_ratio"] >= 1] <- 1
  Log_odds_ratio[which(Log_odds_ratio$Experiment =='ME6'),]["odds_ratio"][Log_odds_ratio[which(Log_odds_ratio$Experiment =='ME6'),]["odds_ratio"] <= 0.01] <- 0.01 #set odds ratio below dilution factor a set minimum for ME6
  Mean_odds_ratio <- Log_odds_ratio %>%
    filter(Experiment!='ME4'&Experiment!='ME1')%>%
    group_by(Long_sample_name,sensor)%>%
    dplyr::mutate(log_odds_ratio= log10(odds_ratio))%>%
    summarise_at(vars(log_odds_ratio),funs(mean(., na.rm=TRUE)))%>%
    dplyr::mutate(geom_mean_odds_ratio= (log_odds_ratio))
  
  Mean_odds_ratio <- Mean_odds_ratio%>%
    select(-c(log_odds_ratio))%>%
    tidyr::pivot_wider(names_from = Long_sample_name, values_from=geom_mean_odds_ratio)%>%
    arrange(Gavage)%>%
    arrange(factor(sensor, levels = order))
  
  write.csv(Mean_odds_ratio,destination_file)
  
  return(Mean_odds_ratio)
}
get_range_odds_ratio <- function(Sample_data, destination_file) {
  order <- c('ynfE15',
             'S_enterica_TCS_ttrSRRIG',
             'hycAWT',
             'torC10',
             'torC17',
             'torC23',
             'torCWT',
             'ynfE17',
             'yeaR23',
             'B_caccae_TCS4',
             'B_caccae_TCS5',
             'B_frag_TCS6',
             'B_theta_TCS2',
             'E_tarda_TCS4',
             'K_pneumonia_TCS8',
             'R_gnavus_TCS3',
             'S_enterica_TCS1',
             'S_enterica_TCS3',
             'Y_enterocolitica_TCS2',
             'C_rodentium_TCS1',
             'L_gasseri_TCS1',
             'E_tarda_TCS7',
             'C_rodentium_TCS3',
             'B_frag_TCS2',
             'C_rodentium_TCS10',
             'C_rodentium_TCS2',
             'R_gnavus_TCS5',
             'C_rodentium_TCS6',
             'K_pneumonia_TCS5',
             'Y_enterocolitica_TCS1',
             'E_tarda_TCS5',
             'V_cholera_TCS3',
             'C_rodentium_TCS7',
             'K_pneumonia_TCS3',
             'K_pneumoniae_TCS10',
             'B_producta_TCS6',
             'Y_enterocolitica_TCS5',
             'K_pneumonia_TCS4',
             'K_pneumonia_TCS9',
             'C_rodentium_TCS8',
             'E_tarda_TCS6',
             'S_enterica_TCS4',
             'C_rodentium_TCS4',
             'S_enterica_TCS5',
             'B_producta_TCS13',
             'C_rodentium_TCS9',
             'Y_enterocolitica_TCS3',
             'B_producta_TCS16',
             'L_plantarum_TCS2',
             'K_oxytoca_TCS1',
             'C_rodentium_TCS11',
             'C_rodentium_TCS5',
             'V_cholera_TCS4',
             'K_pneumonia_TCS7',
             'E_faecalis_TCS2',
             'E_tarda_TCS2',
             'B_frag_TCS3',
             'C_perfringens_TCS2',
             'V_cholera_TCS2',
             'B_producta_TCS12',
             'B_producta_TCS1',
             'E_faecalis_TCS3',
             'C_diff_TCS5')
  Log_odds_ratio <- Sample_data
  Log_odds_ratio["odds_ratio"][Log_odds_ratio["odds_ratio"] <= 0.001] <- 0.001 #set odds ratio below dilution factor a set minimum for LB009
  Log_odds_ratio["odds_ratio"][Log_odds_ratio["odds_ratio"] >= 1] <- 1
  Log_odds_ratio[which(Log_odds_ratio$Experiment =='ME6'),]["odds_ratio"][Log_odds_ratio[which(Log_odds_ratio$Experiment =='ME6'),]["odds_ratio"] <= 0.01] <- 0.01 #set odds ratio below dilution factor a set minimum for ME6
  var_odds_ratio <- Log_odds_ratio %>%
    filter(Experiment!='ME4'&Experiment!='ME1')%>%
    group_by(Long_sample_name,sensor)%>%
    dplyr::mutate(log_odds_ratio= log10(odds_ratio))%>%
    summarise_at(vars(log_odds_ratio),funs(var(., na.rm=TRUE)))%>%
    dplyr::mutate(geom_mean_odds_ratio= (log_odds_ratio))
  
  var_odds_ratio <- var_odds_ratio%>%
    select(-c(log_odds_ratio))%>%
    tidyr::pivot_wider(names_from = Long_sample_name, values_from=geom_mean_odds_ratio)%>%
    arrange(Gavage)%>%
    arrange(factor(sensor, levels = order))
  
  write.csv(var_odds_ratio,destination_file)
  
  return(var_odds_ratio)
}

#lists for plotting loops
Sample_data<-Sample_data_positive_bcs
list_of_sensors <- unique(Sample_data$sensor)
LB009_list_of_conditions <- unique(Sample_data[which(Sample_data$Experiment == 'LB009'  ),] $Long_sample_name)
ME6_list_of_conditions <- unique(Sample_data[which(Sample_data$Experiment == 'ME6'  ),] $Long_sample_name)

#Functions for plotting formatting: 
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "10^", l)
  l <- gsub('\'1\'','',l)
  # return this as an expression
  parse(text=l)
}

######### QC Sample Exclusion Criteria (Spread exclusion) ###########

QC_data <- function(Sample_data_input) {
  Sample_data <- Sample_data_input
  Samples <- unique(Sample_data[which(Sample_data$Experiment != 'ME4' & Sample_data$Experiment != 'ME1'),]$Long_sample_name) 
  R_square_list <- vector("list", length = length(Samples))
  N_barcodes_list <- vector("list", length = length(Samples))
  N_sensors_list <- vector("list", length = length(Samples))
  
  for (i in Samples) {
    
    recovered_bcs <- length(unique(Sample_data[which(Sample_data$Long_sample_name == i),]$bc))
    recovered_sensors <- length(unique(Sample_data[which(Sample_data$Long_sample_name == i),]$sensor))
    N_barcodes_list[[i]] <- recovered_bcs
    N_sensors_list[[i]] <- recovered_sensors
    
    spread_plot_data <- Sample_data[which(Sample_data$Long_sample_name == i & (Sample_data$bc %in% positive_bcs)),]
    
    linear_reg<-lm(log10(spec) ~ log10(none), data=spread_plot_data)
    R_square<-summary(linear_reg)$r.squared
    R_square_list[[i]] <- R_square
    
    Sample_data["none"][Sample_data["none"] <= 0.0001] <- 0.0001
    Sample_data["spec"][Sample_data["spec"] <= 0.0001] <- 0.0001
    
    ###diagonal plot
    spread_plot <- ggplot(data = spread_plot_data, mapping=aes(x=none,y=spec), size=4)+
      geom_point(Sample_data[which(Sample_data$Long_sample_name == i),], mapping=aes(x=none,y=spec), colour='grey',alpha=0.5,size=1)+
      geom_point(spread_plot_data, mapping=aes(x=none,y=spec), colour='dark blue', size=4)+
      ggpmisc::stat_poly_line()+
      ggpmisc::stat_poly_eq()+
      annotate(geom='text',x=0.0001, y=100,hjust=0,label=paste('Barcodes:',recovered_bcs))+
      annotate(geom='text',x=0.0001, y=30,hjust=0,label=paste('Sensors:',recovered_sensors))+
      scale_y_log10( labels = function(none) ifelse(none == 0.0001, "0", format(none, scientific = FALSE)), limits=c(0.0001,5000), breaks=c(0.0001,0.001,0.01,0.1,1,10,100,1000))+
      scale_x_log10( labels = function(spec) ifelse(spec == 0.0001, "0", format(spec, scientific = FALSE)), limits=c(0.0001,5000), breaks=c(0.0001,0.001,0.01,0.1,1,10,100,1000))+
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
    ggsave(spread_plot, file=paste('odds_ratio_plots/spread_plots/',i,'.pdf', sep=''),
           width = (210-30)/2, height = (297-30)/3, units = 'mm', dpi=300)
  }
  R_square_list <- as.data.frame(do.call(rbind, R_square_list))
  N_barcodes_list <-as.data.frame(do.call(rbind, N_barcodes_list))
  N_sensors_list <- as.data.frame(do.call(rbind, N_sensors_list ))
  
  QC_df <- cbind(R_square_list,N_barcodes_list,N_sensors_list)
  names(QC_df)<-c('Rsquare','N_barcodes','N_sensors')
  write.csv(QC_df,'QC_df.csv')
  
  QC_df <- QC_df%>%
    filter(Rsquare>=0.7)
  QC_passing_samples <- rownames(QC_df)
  
  total_reads_per_sample <- as.data.frame(colSums(Filter(is.numeric, assigned_bcs[, -which(names(assigned_bcs) %in% c("bc","total_bc_reads",'X..Sequences.y'))])))
  names(total_reads_per_sample)<-c('reads')
  total_reads_per_sample <- total_reads_per_sample%>%
    filter(reads>=50000)
  
  total_reads_per_sample$Sample <- rownames(total_reads_per_sample)
  
  Sample_codes <- read.csv('All_experiments_sample_codes.csv')
  
  
  # Replace row names in total_reads_per_sample using sample_codes
  matched_rows <- match(total_reads_per_sample$Sample, Sample_codes$Sample)
  new_row_names <- Sample_codes$Long_sample_name[matched_rows]
  total_reads_per_sample$Sample <- new_row_names
  total_reads_per_sample <- as.data.frame(table(total_reads_per_sample$Sample ))%>%
    filter(Freq>1)
  
  QC_passing_data<- Sample_data_input[which(Sample_data_input$Long_sample_name %in% QC_passing_samples & Sample_data_input$Long_sample_name %in% total_reads_per_sample$Var1),]%>%
    filter(reads_none>=10) #filter any barcodes with reads in none sample <10
  write.csv(QC_passing_data, 'QC_passing_data.csv')
  
  return(QC_passing_data)
}

QC_passing_data <- QC_data(Sample_data_positive_bcs)
get_median_odds_ratio(QC_passing_data, 'QC_passing_data_median_odds_ratios.csv')
get_range_odds_ratio(QC_passing_data, 'QC_passing_data_range_odds_ratios.csv')


######### Average Number of reads per barcode per sample ###########

Ave_n_reads <- QC_passing_data %>%
  filter(Experiment!='ME1')%>%
  group_by(bc)%>%
  summarise_at(vars(reads_none),funs(sum(., na.rm=TRUE)))

Ave_n_reads2 <- Ave_n_reads %>%
  summarise_at(vars(reads_none),funs(mean(., na.rm=TRUE)))

Ave_n_reads3 <- Ave_n_reads %>%
  summarise_at(vars(reads_none),funs(median(., na.rm=TRUE)))

Ave_n_reads_plot <- ggplot(data = Ave_n_reads, mapping=aes(x=bc,y=reads_none))+
  geom_point()+
  scale_y_log10()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white'))
Ave_n_reads_plot

######### Plotting reads per barcode ######### 

#import csv of dada2 barcodes & barcode read counts
com_bcs <- read.csv('combined_bcs.csv')%>%
  dplyr::select(-X)
names(com_bcs)[names(com_bcs) == 'cluster'] <- 'bc'

#import csv of barcode sensor assignments
assigned_bcs <- read.csv('assigned_barcodes.csv')%>%
  dplyr::select(-X)
names(assigned_bcs)[names(assigned_bcs) == 'barcode'] <- 'bc' #rename barcode column to bc for later merging with com_bcs

#rename control sensors from plasmid names to control sensors names (for ex. plasmid CR21 pCR21 -> yeaR23)
control_name_pattern <- c('pCR24 PCR Product','pCR23  PCR Product','pCR22  PCR Product','pCR21 PCR Product','pCR18  PCR Product','pCR17 PCR Product','pCR16  PCR Product','pCR15  PCR Product')
control_name_replacement <- c('torC23','torC17','torC10','yeaR23','torCWT','hycAWT','ynfE17','ynfE15')
assigned_bcs$top_by_reads <- stringi::stri_replace_all_regex(assigned_bcs$top_by_reads,
                                                             pattern=control_name_pattern,
                                                             replacement =control_name_replacement,
                                                             vectorize=FALSE)
assigned_bcs$top_prop <- stringi::stri_replace_all_regex(assigned_bcs$top_prop,
                                                         pattern=control_name_pattern,
                                                         replacement =control_name_replacement,
                                                         vectorize=FALSE)

# merge barcode counts & barcode sensor assignments & filter sensors with assignment disagreement 
assigned_bcs <- merge(com_bcs, assigned_bcs, by='bc', all=TRUE) %>%
  select(bc,top_by_reads,top_prop, everything()) %>%
  dplyr::mutate(dis = if_else(top_by_reads==top_prop,'agreement','disagreement'))%>% #add column for if assignment is in agreement or disagreement 
  dplyr::select(dis, everything()) %>% #move dis column to left side of df for easier viewing
  dplyr::mutate(dis= ifelse(is.na(dis)&!is.na(sensor),'manually assigned',dis)) #add manually assigned sensors

samples_list <-c('bc','dis','Gavage', 'AM9', 'ANSOC')
norm_com_bcs <- assigned_bcs[,samples_list]
norm_com_bcs[,c(-1,-2)] <- data.frame(lapply(norm_com_bcs[,c(-1,-2)], function(x) x/sum(x, na.rm=TRUE))) #normalise all barcodes by total number of reads for that sample
norm_com_bcs <- norm_com_bcs[order(-norm_com_bcs$"Gavage"),]

norm_com_bcs <- norm_com_bcs %>%
  tidyr::pivot_longer(!c(bc, dis), names_to = 'sample', values_to ='reads')%>%
  tidyr::replace_na(list(dis = 'none'))

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

reads_per_bc_plot <- ggplot()+
  geom_point(data = norm_com_bcs[which(norm_com_bcs$dis=='agreement'|norm_com_bcs$dis=='none'),], mapping = aes(x=forcats::fct_inorder(as.factor(bc)), y=reads), colour='grey', size=0.8, alpha=0.5, shape=16)+
  scale_y_log10('Fraction of Reads',
                breaks = scales::trans_breaks("log10", function(x) 10^x, n=5),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  xlab('Barcode')+
  labs(fill = "Sample")+
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white'),
        text = element_text(size = 8),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(family = "Helvetica", size = 10, face = "bold"),
        axis.title = element_text(family = "Helvetica", size = 8),
        axis.text = element_text(family = "Helvetica", size = 8),
        legend.text = element_text(family = "Helvetica", size = 10),
        legend.title = element_text(family = "Helvetica", size = 10),
        strip.text = element_blank())+
        facet_grid(cols = vars(dis))
reads_per_bc_plot
ggsave(reads_per_bc_plot, file=paste('reads_per_bc_plot_assignment.pdf'),
       width = (210-30)*(2/3), height = (297-30)/6, units = 'mm', dpi=300)


######### T tests #######

#create list to store kruskal test vals. 
Kruskal_test_lists <- vector("list", length = length(list_of_sensors))

Sample_data <- QC_passing_data
#do T test between gavage and 5 mg strep groups for each sensor
for (i in list_of_sensors) {
  sensor_df <- Sample_data[which(Sample_data$sensor == i & Sample_data$Experiment == 'ME6'&
                                   (Sample_data$Day == 0|Sample_data$Day == 2)&
                                   (Sample_data$Group == 'gavage'|Sample_data$Group == '5 mg strep')),] %>%
    dplyr::mutate(Group = forcats::fct_relevel(Group,'gavage','5 mg strep'))
  sensor_df["odds_ratio"][sensor_df["odds_ratio"] <= 0.01] <- 0.01
  write.csv(sensor_df,paste(i,'.csv',sep = ''))
  
  if (length(unique(sensor_df$Group))>1) {
    P_val<-kruskal.test(odds_ratio ~ Group, data = sensor_df)$p.value
    Kruskal_test_lists[[i]] <- P_val
    
    wilcox.test(odds_ratio ~ Group, data = sensor_df) #Mann Whitney U Test (gives similar results to Kruskal Wallis)
  }
}

#convert results to dataframe and format
T_test_lists2 <- as.data.frame(do.call(rbind, Kruskal_test_lists))
T_test_lists2$sensor <- rownames(T_test_lists2)
T_test_lists2$Correct_P <- T_test_lists2$V1*nrow(T_test_lists2)


write.csv(T_test_lists2,'T_test_sensors.csv')

#get median odds ratios in order to filter for OFF sensors 
Median_odds_ratio <- Sample_data_positive_bcs %>%
  filter(Experiment!='ME4')%>%
  group_by(Long_sample_name,sensor)%>%
  summarise_at(vars(odds_ratio),funs(median(.)))%>%
  filter(Long_sample_name=='Gavage')
Median_odds_ratio["odds_ratio"][Median_odds_ratio["odds_ratio"] <= 0.01] <- 0.01

off_sensors <- Median_odds_ratio %>%
  filter(odds_ratio<=0.4)
off_sensors <- unique(off_sensors$sensor)

#filter T tests for those with p<0.05 (with multiple test comparison correction)
T_test_lists3<- T_test_lists2[which(T_test_lists2$sensor %in% off_sensors),] %>%
  filter(as.numeric(V1)<=0.05/(63)) #multiple test correction



######### ME4 barcodes per sample #########
ME4 <- assigned_bcs%>%
  filter(sensor %!in% c('fabRbc1','fabRbc2','fabRbc3','fabRbc4'))%>%
  select(ME4_D12_G10,
         ME4_D12_G11,
         ME4_D2_G10,
         ME4_D2_G11,
         ME4_D2_G20,
         ME4_D2_G21,
         ME4_D2_G30,
         ME4_D2_G31,
         ME4_D4_G10,
         ME4_D4_G11,
         ME4_D7_G10,
         ME4_D7_G11,
         ME4_gavage)
nonzero <- function(x) sum(x != 0)
Number_of_barcodes_per_sample <- plyr::numcolwise(nonzero)(ME4)%>%
  tidyr::pivot_longer(cols = everything(),names_to = 'Sample', values_to='N_barcodes')

ME4 <- assigned_bcs %>%
  select(-total_bc_reads)%>%
  tidyr::pivot_longer(!c(bc,sensor,barcode), names_to = 'Sample', values_to = 'reads')

#read sample codes csv 
Sample_codes <- read.csv('All_experiments_sample_codes.csv')

#merge sample data with sample codes and positive normalisation reads, normalise and calculate odds ratio
ME42 <-merge(Sample_codes,ME4, by='Sample')%>%
  filter(Experiment=='ME4')%>%
  filter(reads!=0)%>%
  group_by(Sample)%>%
  summarise(count=n_distinct(sensor))
  

ME6 <-merge(Sample_codes,ME4, by='Sample')%>%
  filter(Experiment=='ME6')%>%
  filter(reads!=0)%>%
  group_by(Sample)%>%
  summarise(count=n_distinct(sensor))

ME6 <-merge(Sample_codes,ME4, by='Sample')%>%
  filter(Experiment=='ME6')%>%
  filter(reads!=0)%>%
  group_by(Sample)%>%
  summarise(count=n_distinct(bc))

