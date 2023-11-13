
library(ggplot2)
library(tidyr)
library(dplyr)

#For R Studio only: set workind directory to location of script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#read in alignments 
alignments <- read.csv('trimmed_barcode_alignments.csv')
#get barcode from Name column
alignments$barcode <- alignments$Name
#clean up sensor name from description column 
alignments$sensor <- sapply(strsplit(alignments$Description, split='mapped'),'[',1)
alignments$sensor <- gsub(".*from ","",alignments$sensor)
alignments$sensor <- gsub(".fasta ","",alignments$sensor)
alignments$sensor <- gsub(" \\(reversed\\)","",alignments$sensor)

#TCS assignment to bc by proportional top number of reads (diff TCS & barcodes have diff abundances, so will have different numbers of misassigned reads. Normalising per TCS abundance can help with barcode mis-assignment)
top_proportional <- alignments %>%
  dplyr:: group_by(sensor) %>%
  dplyr:: mutate(tot_rds_sensor = sum(X..Sequences))%>% #get total reads per sensor
  dplyr::mutate(prop_rds_barcode=X..Sequences/tot_rds_sensor)%>% #normalise reads for a barcode by it's total reads per assigned sensor
  dplyr::group_by(barcode)%>% #group by bc
  dplyr::slice_max(prop_rds_barcode)%>% #select TCS assigned with the higher proportional number of reads
  dplyr::mutate(top_prop=sensor)%>%
  dplyr::select(c(barcode, top_prop))

#TCS assignment to bc by top number of reads
top_by_reads <- alignments%>%
  dplyr::group_by(barcode)%>% #group by bc
  dplyr::slice_max(X..Sequences)%>% #select TCS with highest read assigned to that bc
  dplyr::mutate(top_by_reads=sensor)%>%
  dplyr::select(-c(sensor, Name, Description))

alignments <- merge(alignments, top_proportional, by='barcode')
alignments <- merge(alignments, top_by_reads, by='barcode')

#any disagreements
disagrements <- unique(alignments[which(alignments$top_by_reads != alignments$top_prop),]$barcode)

#all assigned barcodes
assigned_barcodes <- alignments %>%
  distinct(barcode, .keep_all= TRUE)%>%
  select(c('barcode','top_by_reads','top_prop','X..Sequences.y')) %>%
  filter(X..Sequences.y >=4) #filter barcodes that have at least 4 trimmed long reads for their assignment (read counts are from Geneious alignment that include reference -> increases count by 1)
write.csv(assigned_barcodes, 'assigned_barcodes.csv')

#barcodes per sensor: by net total reads assignment
barcode_count_net <- top_by_reads %>% 
  group_by(top_by_reads) %>% 
  mutate(barcode_count = n())%>%
  distinct(top_by_reads, .keep_all= TRUE)

#barcodes per sensor: by proportional assignment
barcode_count_prop <- top_proportional %>% 
  group_by(top_prop) %>% 
  mutate(barcode_count = n())%>%
  distinct(top_prop, .keep_all= TRUE)

barcode_count <-merge(barcode_count_net,barcode_count_prop, by='barcode')
names(barcode_count)[names(barcode_count) == 'barcode_count.x'] <- 'barcode_count_net'
names(barcode_count)[names(barcode_count) == 'barcode_count.y'] <- 'barcode_count_prop'
names(barcode_count)[names(barcode_count) == 'top_prop'] <- 'sensor'

barcode_count<-barcode_count%>%
  tidyr::pivot_longer(cols = c('barcode_count_net','barcode_count_prop'),
               names_to = 'assignment_type',
               values_to = 'n_assigned_bcs')%>%
  dplyr::select(c('sensor','assignment_type','n_assigned_bcs'))%>% 
  group_by(assignment_type)

#barcodes per sensor: proportional vs net assignment
barcode_count <- barcode_count[order(-barcode_count$n_assigned_bcs),]
bc_per_sensor <- ggplot(data=barcode_count)+
  geom_point(mapping=aes(x=forcats::fct_inorder(gsub('_',' ',sensor)), y=n_assigned_bcs,colour=assignment_type), alpha=0.5)+
  xlab("Sensor")+
  ylab("Number of Assigned Barcodes")+
  ggtitle('Barcode assignment')+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 75, hjust=1),
        plot.margin = margin(20, 20, 20, 20),
        legend.position = 'right')
  #scale_color_manual(labels = c('Net assignment','Proportional assignment'))
bc_per_sensor
ggsave(bc_per_sensor, file='barcodes_per_sensor.png', width = 9, height = 5)




