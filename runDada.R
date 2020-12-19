library(dada2)
library(dplyr)
library(assertthat)
library(stringr)
library(readxl)
library(ggplot2)
library(ggbeeswarm)
library(ComplexHeatmap)
library(circlize)
library(phyloseq)
library(phangorn)
library(msa)
#https://benjjneb.github.io/dada2/tutorial.html

read.table("metadata/SRP141397.metadata",sep="\t",header = TRUE) %>%
  dplyr::filter(str_detect(experiment_title,'16S')) %>% dplyr::mutate(full1=paste0(paste("raw",study_accession,experiment_accession,run_accession,sep="/"),"_1.fastq")) %>%
  dplyr::filter(str_detect(experiment_title,'16S')) %>% dplyr::mutate(full2=paste0(paste("raw",study_accession,experiment_accession,run_accession,sep="/"),"_2.fastq")) %>%
  dplyr::filter(str_detect(experiment_title,'16S')) %>% dplyr::mutate(pair1=paste0(paste("intermediates","fastq",run_accession,sep="/"),"_1.fastq.gz")) %>%
  dplyr::filter(str_detect(experiment_title,'16S')) %>% dplyr::mutate(pair2=paste0(paste("intermediates","fastq",run_accession,sep="/"),"_2.fastq.gz"))->
  sra_metadata

sapply(sra_metadata$pair1,function(x){assert_that(file.exists(x))})
sapply(sra_metadata$pair2,function(x){assert_that(file.exists(x))})

fnFs <- sra_metadata$pair1
fnRs <- sra_metadata$pair2

sample.names <- sra_metadata$experiment_title

plotQualityProfile(fnFs[1:2])

plotQualityProfile(fnRs[1:2])

# Place intermediates/filtered files in intermediates/filtered/ subdirectory
filtFs <- file.path("intermediates/filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("intermediates/filtered", paste0(sample.names, "_R_filt.fastq.gz"))
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
seqtab.chim <- makeSequenceTable(mergers)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab.chim, method="consensus", multithread=TRUE, verbose=TRUE)

noChimeras<-TRUE
if(noChimeras){
  seqtab<-seqtab.nochim
}else{
  seqtab<-seqtab.chim
}

trseqtab<-t(seqtab)
#assign taxonomy to nochim
taxa <- assignTaxonomy(seqtab, "SILVA_DADA/silva_nr_v123_train_set.fa.gz", multithread=TRUE)



trseqtabdfnoseq<-as.data.frame(trseqtab)
trseqtabdf$seq<-row.names(trseqtabdfnoseq)
abundance<-as.data.frame(colSums(trseqtabdfnoseq))


#this has experimental metadata
tabl11<-read_excel("metadata/table1.xls",sheet = 11)
names(tabl11)<-c("SampleID","Group","Type","Case_Control","Delivery","numReads","nonhost_shotgun_reads")
tabl11$SampleID<-str_replace(tabl11$SampleID,'16S GeneBlock ','POS.control')
tabl11$Type<-factor(tabl11$Type,levels=c("Air Swab","Blank","H2O","Placenta Fetal Side Delivery Biopsy","Placenta Maternal Side Delivery Biopsy","Positive Control (Gene Block)","Positive Control (Shotgun))","Maternal Saliva Microbiome Enrollment","Vaginal Swab Microbiome Enrollment"))
tabl11$Case_Control<-factor(tabl11$Case_Control,levels=c("Preterm","Term","n/a"))
tabl11$Delivery<-factor(tabl11$Delivery,levels=c("SVD","C-Section","n/a"))
tabl11 %>% arrange(Type,Case_Control,Delivery)->meta.sorted

#useful for qiime
meta.sorted %>% 
  mutate(Type=str_replace_all(pattern = '\\(',replacement='',Type)) %>%
  mutate(Type=str_replace_all(pattern = '\\)',replacement='',Type)) %>%
  mutate(Description=paste(Type,Case_Control,Delivery,sep="__"))  %>% 
  mutate(Description=str_replace_all(pattern = ' ',replacement = '_',string = Description)) %>%
  mutate(BarcodeSequence="") %>%
  mutate(LinkerPrimerSequence = "") %>%
  select(SampleID,BarcodeSequence,LinkerPrimerSequence,Type,Case_Control,Delivery,Description) %>% filter(SampleID != 'Shotgun Positive Control') %>%
  dplyr::rename("#SampleID"="SampleID") -> qiimetable
write.table(qiimetable,"intermediates/mapfile.txt",quote=FALSE,sep="\t",row.names=FALSE)



abundance$sample<-str_replace(row.names(abundance),'_16S','')
names(abundance)<-c("dadareads","SampleID")
abundance<-merge(abundance,meta.sorted,by="SampleID")
abundanceMelted<-reshape2::melt(abundance)
ggplot(abundanceMelted,aes(Type,value))+geom_beeswarm()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+facet_wrap("variable")


taxadf<-as.data.frame(taxa)
taxadf$seq<-row.names(taxadf)
trseqtabdf %>% merge(taxadf,by="seq") %>% group_by(Phylum,Genus) %>% select(-c(seq,Kingdom,Class,Order,Family)) %>% summarize_all(sum) -> genusCnts

genusCntsNrml<-data.frame(Phylum=genusCnts$Phylum,Genus=genusCnts$Genus,sweep(genusCnts[,-c(1,2)], 2, colSums(genusCnts[,-c(1,2)]), FUN ="/" ))
colnames(genusCntsNrml)<-str_replace(colnames(genusCntsNrml),'_16S','')




keepRows<-rowSums(genusCntsNrml[,-c(1,2)],na.rm = TRUE)>0.001
subset(genusCntsNrml,keepRows)->genusCntsNrmlAboveThreshold
cols<-meta.sorted[sort(sample(1:nrow(meta.sorted),50)),"SampleID"] %>% dplyr::pull('SampleID')


Type<-meta.sorted[match(colnames(genusCntsNrmlAboveThreshold[cols]),meta.sorted$SampleID),"Type"]  %>% dplyr::pull('Type')
Delivery<-meta.sorted[match(colnames(genusCntsNrmlAboveThreshold[cols]),meta.sorted$SampleID),"Delivery"] %>% dplyr::pull('Delivery')
Term<-meta.sorted[match(colnames(genusCntsNrmlAboveThreshold[cols]),meta.sorted$SampleID),"Case_Control"] %>% dplyr::pull('Case_Control')

set.seed(3)
col_fun = colorRamp2(c(0, 0.01, .3, 1), c("white", "blue", "green", "red"))
col_fun(seq(0, 1))
ha = HeatmapAnnotation(Type = Type, Delivery=Delivery,Term=Term,annotation_name_side = "left")
Heatmap(as.matrix(genusCntsNrmlAboveThreshold[,cols]),
        col=col_fun,
        top_annotation = ha,
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        show_row_dend = FALSE) +
  Heatmap(genusCntsNrmlAboveThreshold[,"Phylum"], name = "Phylum", 
          top_annotation = HeatmapAnnotation(summary = anno_summary(height = unit(2, "cm"))),
          width = unit(15, "mm"))

meta.sorted.inseq<-data.frame(meta.sorted[meta.sorted$SampleID %in% row.names(seqtab),])
rownames(meta.sorted.inseq)<-meta.sorted.inseq$SampleID
meta.sorted.inseq %>% select(-c(SampleID,numReads,nonhost_shotgun_reads)) -> meta.sorted.inseq.samp
rownames(seqtab)<-str_replace(rownames(seqtab),'_16S','')

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(meta.sorted.inseq.samp), 
               tax_table(taxa))



#seqs <- dada2::getSequences(seqtab)
#names(seqs) <- seqs # This propagates to the tip labels of the tree
#mult <- msa::msa(seqs, method="ClustalW", type="dna", order="input")
#save(mult,file="16s.msa")

#https://github.com/benjjneb/dada2/issues/88
library(DECIPHER)
dseqs<-DNAStringSet(seqs)
#about an hour
#2762+45+273+568+31+125+599
alignment <- AlignSeqs(dseqs, anchor=NA)
save(alignment,file="decipher.msa")
#https://compbiocore.github.io/metagenomics-workshop/assets/DADA2_tutorial.html

#Construct Phylogenetic Tree
#Extract sequences from DADA2 output

sequences<-getSequences(seqtab.nochim)
names(sequences)<-sequences
# Run Sequence Alignment (MSA) using DECIPHER

alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
# Change sequence alignment output into a phyDat structure

phang.align <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
#Create distance matrix

dm <- phangorn::dist.ml(phang.align)

#Perform Neighbor joining, 20 minutes

treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order

#Internal maximum likelihood

fit = pml(treeNJ, data=phang.align)
#negative edges length changed to 0!
  
#  fitGTR <- update(fit, k=4, inv=0.2)
#fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                    rearrangement = "stochastic", control = pml.control(trace = 0))
#alignment <- read.dna("alignment.fasta",format="fasta",as.matrix=TRUE)
alignment.rax.gtr <- raxml(alignment,
                           m="GTRGAMMAIX", # model
                           f="a", # best tree and bootstrap
                           p=1234, # random number seed
                           x=2345, # random seed for rapid bootstrapping
                           N=100, # number of bootstrap replicates
                           file="alignment", # name of output files
                           exec="raxmlHPC-PTHREADS-SSE3", # name of executable
                           threads=20
)

#https://f1000research.com/articles/5-1492/v1
phang.align <- as.phyDat(mult, type="DNA", names=getSequence(seqtab))
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

#fitGTR <- update(fit, k=4, inv=0.2)
#fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
     #               rearrangement = "stochastic", control = pml.control(trace = 0))







#pslog <- transform_sample_counts(ps, function(x) log(1 + x))
#out.wuf.log <- ordinate(pslog, method = "MDS", distance = "unifrac")
