library(dada2)
library(dplyr)
library(assertthat)

#https://benjjneb.github.io/dada2/tutorial.html

read.table("metadata/SRP141397.metadata",sep="\t",header = TRUE) %>%
  dplyr::filter(str_detect(library_name,'16S')) %>% dplyr::mutate(pair1=paste0(paste("raw",study_accession,experiment_accession,run_accession,sep="/"),"_1.fastq")) %>%
  dplyr::filter(str_detect(library_name,'16S')) %>% dplyr::mutate(pair2=paste0(paste("raw",study_accession,experiment_accession,run_accession,sep="/"),"_2.fastq")) ->
  sra_metadata

sapply(sra_metadata$pair1,function(x){assert_that(file.exists(x))})
sapply(sra_metadata$pair2,function(x){assert_that(file.exists(x))})

fnFs <- sra_metadata$pair1
fnRs <- sra_metadata$pair2

sample.names <- sra_metadata$library_name

plotQualityProfile(fnFs[1:2])

plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# dada2 `filterAndTrim` parameters
#truncQ	
#(Optional). Default 2. Truncate reads at the first instance of a quality score less than or equal to truncQ.

#truncLen	
#(Optional). Default 0 (no truncation). Truncate reads after truncLen bases. Reads shorter than this are discarded.

#trimLeft	
#(Optional). Default 0. The number of nucleotides to remove from the start of each read. If both truncLen and trimLeft are provided, filtered reads will have length truncLen-trimLeft.

#trimRight	
#(Optional). Default 0. The number of nucleotides to remove from the end of each read. If both truncLen and trimRight are provided, truncation will be performed after trimRight is enforced.

#maxLen	
#(Optional). Default Inf (no maximum). Remove reads with length greater than maxLen. maxLen is enforced before trimming and truncation.

#minLen	
#(Optional). Default 20. Remove reads with length less than minLen. minLen is enforced after trimming and truncation.

#maxN	
#(Optional). Default 0. After truncation, sequences with more than maxN Ns will be discarded. Note that dada does not allow Ns.

#minQ	
#(Optional). Default 0. After truncation, reads contain a quality score less than minQ will be discarded.

#maxEE	
#(Optional). Default Inf (no EE filtering). After truncation, reads with higher than maxEE "expected errors" will be discarded. Expected errors are calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10))

#rm.phix	
#(Optional). Default TRUE. If TRUE, discard reads that match against the phiX genome, as determined by isPhiX.

#rm.lowcomplex	
#(Optional). Default 0. If greater than 0, reads with an effective number of kmers less than this value will be removed. The effective number of kmers is determined by seqComplexity using a Shannon information approximation. The default kmer-size is 2, and therefore perfectly random sequences will approach an effective kmer number of 16 = 4 (nucleotides) ^ 2 (kmer size).

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

errF <- learnErrors(filtFs, multithread=TRUE)
#105662880 total bases in 440262 reads from 31 samples will be used for learning the error rates

errR <- learnErrors(filtRs, multithread=TRUE)


dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)