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
                                                     dplyr::mutate(full2=paste0(paste("raw",study_accession,experiment_accession,run_accession,sep="/"),"_2.fastq")) %>%
                                                     dplyr::mutate(pair1=paste0(paste("intermediates","fastq",run_accession,sep="/"),"_1.fastq.gz")) %>%
                                                     dplyr::mutate(pair2=paste0(paste("intermediates","fastq",run_accession,sep="/"),"_2.fastq.gz"))->
  sra_metadata

sample.names <- sra_metadata$experiment_title

filtFs <- file.path("intermediates/filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("intermediates/filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

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

