library(dada2)
library(dplyr)
library(assertthat)
library(stringr)
library(readxl)
library(ggplot2)
library(ggbeeswarm)
#https://benjjneb.github.io/dada2/tutorial.html

read.table("metadata/SRP141397.metadata",sep="\t",header = TRUE) %>%
  dplyr::filter(str_detect(experiment_title,'16S')) %>% dplyr::mutate(full1=paste0(paste("raw",study_accession,experiment_accession,run_accession,sep="/"),"_1.fastq")) %>%
  dplyr::filter(str_detect(experiment_title,'16S')) %>% dplyr::mutate(full2=paste0(paste("raw",study_accession,experiment_accession,run_accession,sep="/"),"_2.fastq")) %>%
  dplyr::filter(str_detect(experiment_title,'16S')) %>% dplyr::mutate(pair1=paste0(paste("sratofastq",run_accession,sep="/"),"_1.fastq.gz")) %>%
  dplyr::filter(str_detect(experiment_title,'16S')) %>% dplyr::mutate(pair2=paste0(paste("sratofastq",run_accession,sep="/"),"_2.fastq.gz"))->
  sra_metadata

sapply(sra_metadata$pair1,function(x){assert_that(file.exists(x))})
sapply(sra_metadata$pair2,function(x){assert_that(file.exists(x))})

fnFs <- sra_metadata$pair1
fnRs <- sra_metadata$pair2

sample.names <- sra_metadata$experiment_title

plotQualityProfile(fnFs[1:2])

plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names


#Sequence data was processed using DADA2 [42]. Reads were trimmed from 251 bases to 240. Dereplication, error modeling, denoising, pair merging, and chimera removal were performed using default parameters.


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

#merge overlapping pairs?
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

#construct sequence table
seqtab <- makeSequenceTable(mergers)
taxa <- assignTaxonomy(seqtab, "SILVA_DADA/silva_nr_v123_train_set.fa.gz", multithread=TRUE)


#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


#assign taxonomy to nochim
taxaNochim <- assignTaxonomy(seqtab.nochim, "SILVA_DADA/silva_nr_v123_train_set.fa.gz", multithread=TRUE)

#this has experimental 
tabl11<-read_excel("metadata/table1.xls",sheet = 11)
names(tabl11)<-c("SampleID","Group","Type","Case_Control","Delivery","numReads","nonhost_shotgun_reads")
tabl11$SampleID<-str_replace(tabl11$SampleID,'16S GeneBlock ','POS.control')

trseqtabdfnoseq<-as.data.frame(trseqtab)
abundance<-as.data.frame(colSums(trseqtabdfnoseq))
abundance$sample<-str_replace(row.names(abundance),'_16S','')
names(abundance)<-c("dadareads","SampleID")
abundance<-merge(abundance,tabl11,by="SampleID")
abundanceMelted<-reshape2::melt(abundance)

trseqtabdf$seq<-row.names(trseqtabdfnoseq)
taxadf<-as.data.frame(taxa)
taxadf$seq<-row.names(taxadf)
trseqtabdf %>% merge(taxadf,by="seq") %>% group_by(Phylum,Genus) %>% select(-c(seq,Kingdom,Class,Order,Family)) %>% summarize_all(sum) -> genusCnts

genusCntsNrml<-data.frame(Phylum=genusCnts$Phylum,Genus=genusCnts$Genus,sweep(genusCnts[,-c(1,2)], 2, colSums(genusCnts[,-c(1,2)]), FUN ="/" ))
colnames(genusCntsNrml)<-str_replace(colnames(genusCntsNrml),'_16S','')
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 0.01, .3, 1), c("white", "blue", "green", "red"))
col_fun(seq(0, 1))


keepRows<-rowSums(genusCntsNrml[,-c(1,2)],na.rm = TRUE)>0.001
subset(genusCntsNrml,keepRows)->genusCntsNrmlAboveThreshold
cols<-tabl11[sort(sample(1:nrow(tabl11),50)),"SampleID"]

Type<-tabl11[match(colnames(genusCntsNrmlAboveThreshold[cols]),tabl11$SampleID),"Type"]
Delivery<-tabl11[match(colnames(genusCntsNrmlAboveThreshold[cols]),tabl11$SampleID),"Delivery"]

ha = HeatmapAnnotation(Type = Type, Delivery=Delivery,annotation_name_side = "left")
Heatmap(as.matrix(genusCntsNrmlAboveThreshold[,cols]),
        col=col_fun,
        top_annotation = ha,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_dend = FALSE) +
  Heatmap(genusCntsNrmlAboveThreshold[,"Phylum"], name = "Phylum", 
          top_annotation = HeatmapAnnotation(summary = anno_summary(height = unit(2, "cm"))),
          width = unit(15, "mm"))

