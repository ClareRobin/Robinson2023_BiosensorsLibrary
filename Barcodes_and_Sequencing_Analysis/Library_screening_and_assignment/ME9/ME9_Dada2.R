### DADA2 Processing for barcode identification and relative abundance calculations for ME6 ###

library(dada2)

experiment <- 'ME9'

#path of forward and reverse reads for experiment 
pathF <-'/Users/clare/Desktop/PhD/Lab/1_TCS_Library/2023_11_17_ME9_second_lib_screen/Sequencing_data_and_analysis/Sequencing_Data/Forward_reads'
pathR <-'/Users/clare/Desktop/PhD/Lab/1_TCS_Library/2023_11_17_ME9_second_lib_screen/Sequencing_data_and_analysis/Sequencing_Data/Reverse_reads'
  
#assign directory and file names
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...
fastqFs <- sort(list.files(pathF, pattern="fq.gz"))
fastqRs <- sort(list.files(pathR, pattern="fq.gz"))

#make sure all forward reverse read pairs are present
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

##Filter and trim (Takes ~2 hours for 66 samples ~2 million reads each)
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(106,106),#truncate at the expected barcode length (remove the read through into adaptors from this being PE150 sequencing)
              maxEE=1, #max expected errors 
              truncQ=2, #truncate reads at the first instance of a quality score less than or equal to truncQ
              rm.phix=TRUE,#remove reads that match against the phiX genome
              compress=TRUE, 
              verbose=TRUE,
              multithread=TRUE)

####### infer sequence variants ####### 
#get filtered file paths
filtFs <- list.files(filtpathF, pattern="fq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fq.gz", full.names = TRUE)

#get sample names
sample.names <-  gsub('_EKDL230052792-1A_22GWC7LT3_L1_1.fq.gz','',list.files(filtpathF, pattern="fq.gz", full.names = F))
sample.namesR <- gsub('_EKDL230052792-1A_22GWC7LT3_L1_2.fq.gz','',list.files(filtpathR, pattern="fq.gz", full.names = F))

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#unfiltered file paths & names
Fs <- list.files(pathF, pattern="fq.gz", full.names = TRUE)
Rs <- list.files(pathR, pattern="fq.gz", full.names = TRUE)
names(Fs) <- sample.names
names(Rs) <- sample.names


set.seed(100)

# Learn forward error rates
errF <- learnErrors(filtFs, 
                    nbases=1e8, #minimum number of bases for error rate learning
                    multithread=TRUE)


# Learn forward error rates
errF <- learnErrors('/Users/clare/Desktop/PhD/Lab/1_TCS_Library/2023_11_17_ME9_second_lib_screen/Sequencing_data_and_analysis/Sequencing_Data/Forward_reads/filtered_sample_subset_for_Dada2_learning/subset_F.fastq.gz', 
                    nbases=1e8, #minimum number of bases for error rate learning
                    multithread=TRUE)

# Learn reverse error rates
errR <- learnErrors('/Users/clare/Desktop/PhD/Lab/1_TCS_Library/2023_11_17_ME9_second_lib_screen/Sequencing_data_and_analysis/Sequencing_Data/Reverse_reads/filtered_subset_for_error/subset_R.fastq.gz', 
                    nbases=1e8,
                    multithread=TRUE)

# Learn forward error rates
errF <- learnErrors(filtRs, 
                    nbases=1e8, #minimum number of bases for error rate learning
                    multithread=TRUE)


#visualise error rates: 
ggplot2::ggsave(plotErrors(errF, nominalQ=TRUE),file=paste(experiment,'_error_rate.png',sep=''))

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

tracks <- vector('list',length(sample.names))
names(tracks) <- sample.names

#function for count number of reads in filtered and unfiltered fastq files (without unzipping) - this pastes a bash script and runs it with system()
count_reads <- function(file) {
  bash_input <- paste('echo $(zcat ',file,'| wc -l) / 4 | bc', sep='')
  read_count <- system(bash_input, intern = TRUE)
  return(read_count)
}

for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]]) #read in and dereplicate fastq forward
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]]) #read in and dereplicate fastq reverse
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
  
  # Construct sequence table and remove chimeras
  seqtab <- makeSequenceTable(merger)
  # Remove chimeras
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)
  
  #count number of reads in filtered and unfiltered files (without unzipping) - this pastes a bash script and runs it with system()
  filtF_N <- count_reads(filtFs[[sam]])
  Fs_N <- count_reads(Fs[[sam]])
  
  #track reads throughout the pipeline & save as csv 
  track <- cbind(Fs_N,filtF_N,sum(getUniques(ddF)), sum(getUniques(ddR)), sum(getUniques(merger)), rowSums(seqtab.nochim))
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sam
  tracks[[sam]] <- track
}

#remove the derep objects to not take up R memory space
rm(derepF)
rm(derepR)

#save the pipeline read number data as csv 
tracks <- do.call(rbind, tracks)
write.csv(tracks, file = paste(experiment,'_dada2_read_tracking_through_pipeline.csv',sep = ''))


# Construct final sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

#only allow barcodes that are the expected 106bp 
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 106:106]

#save ASVs (barcodes) and their abundances in each sample
seqtab_barcodes <- tibble::rownames_to_column(as.data.frame(t(seqtab)), "barcode")
write.csv(seqtab_barcodes, file=paste('Barcode_abundances_',experiment,'.csv',sep=''))

