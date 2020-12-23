from utils import metautils
from snakemake.utils import R
import os
project=metautils.srpMeta('SRP141397')
project.process()

#use S3 bucket with fastqs ready to go
useAWS=True

configfile: "config.yaml"

#set this to where your conda env containing python2 and qiime is located
qiime1env=config['qiime1env']

# Sequence data was processed using DADA2 [42]. Reads were trimmed from 251 bases to 240. Dereplication, error modeling, denoising, pair merging, and chimera removal were
# performed using default parameters. Taxonomic assignments were generated by comparison to the Silva reference database [43,44,45]. Samples with fewer than 100 reads were
# excluded from downstream analysis. Beta diversity analysis (unweighted UniFrac and weighted UniFrac) was performed using the beta_diversity.py script from Qiime 1.9.1 [46].
# Jaccard, PCoA, and PERMANOVA analysis were performed using functions vegdist, pcoa, and adonis from the R package “vegan” [47].

# To go from SRA do this first
# pysradb download  --out-dir ./raw -p SRP141397
rule fetchOrExtractFastq:
    output: project.get16SFiles(flatten=True,compress=True)
    threads: 8
    run:
        if useAWS:
          shell("aws s3 cp --recursive --no-sign-request s3://placenta/intermediates/fastq intermediates/fastq")
        else:
          shell("pysradb download -y -t {threads} --out-dir ./raw -p SRP141397")
          shell("""parallel-fastq-dump --threads 4 --outdir intermediates/fastq --split-files --tmpdir /tmp --gzip -s `find raw/ -name "*sra"`""")

rule getFiltered:
    input: project.get16SFiles(flatten=True,compress=True)
    output: project.getFiltered16SFiles()
    shell: """
           echo "source('dataFilterTrim.R')" |  R --quiet --no-save
           """

rule fastq2fasta:
    input: "intermediates/filtered/{file}.fastq.gz"
    output: "intermediates/fasta/{file}.fasta.gz"
    shell:
        """
        gunzip -c {input} | fastq_to_fasta -z -o {output}
        """

#named and sequencially numbered
rule fasta2qiimedef:
    input: pair1="intermediates/fasta/{code}_16S_F_filt.fasta.gz", pair2="intermediates/fasta/{code}_16S_R_filt.fasta.gz"
    output: "intermediates/qiime_ready/{code}_16S.fa"
    shell: """
           gunzip -c {input.pair1} {input.pair2} | awk '/^>/{{print ">{wildcards.code}_" ++i; next}}{{print}}' > {output}
           """

rule getFasta:
    input: project.getFastaFiles()

rule getQiimeReady:
    input: project.getQiimeReady()

# concatenate all the sequences into one file, keeping only samples with at least 100 sequences
rule combine:
    input: project.getQiimeReady()
    output: "intermediates/seqs.fna"
    shell: """
           for f in intermediates/qiime_ready/*.fa; do cat $f | perl -e 'while(<>){{$lines.=$_;/>/ && $cnt++;}}if($cnt>=100){{print $lines;}}' >> {output}; done
           """

##### qiime
useOpen=False
if useOpen:
  qiimeoutdir="intermediates/OTU_SILVA_OPEN"
  otutree="intermediates/OTU_SILVA_OPEN/rep_set.tre"
  taxonomy="intermediates/OTU_SILVA_OPEN/uclust_assigned_taxonomy/rep_set_tax_assignments.txt"
  biom="intermediates/OTU_SILVA_OPEN/otu_table_mc2_w_tax.biom"
  program="pick_open_reference_otus.py"
  phylumfile="otu_table_mc2_w_tax_L2.txt"
else:
  qiimeoutdir="intermediates/OTU_SILVA_CLOSED"
  otutree="SILVA_119_QIIME_release/97_FastTree_trees/Silva_119_rep_set97_aligned_16S_only_pfiltered.tre"
  taxonomy="SILVA_119_QIIME_release/taxonomy/97/taxonomy_97_7_levels.txt"
  biom="intermediates/OTU_SILVA_CLOSED/otu_table.biom"
  program="pick_closed_reference_otus.py"
  phylumfile="otu_table_L2.txt"

# 2 hours
# choice of UCLUST or usearch61
# rule OTUopen:
#     input: "intermediates/seqs.fna"
#     output: biom, otutree
#     shell: """
#            export PATH={qiime1env}/bin/:$PATH
#            {qiime1env}/bin/pick_open_reference_otus.py   -i {input} -r SILVA_119_QIIME_release/rep_set/97/Silva_119_rep_set97.fna -o {qiimeoutdir}  -s 0.1 -m usearch61 -p otu_SILVA_settings.txt
#            """

#35 minutes
rule makebiom:
    input: "intermediates/seqs.fna"
    output: biom
    shell: """
           export PATH={qiime1env}/bin/:$PATH
           {qiime1env}/bin/pick_closed_reference_otus.py -f -i {input} -r SILVA_119_QIIME_release/rep_set/97/Silva_119_rep_set97.fna -o {qiimeoutdir} -p otu_SILVA_settings.txt -t {taxonomy}
           """

rule humanreadableotutable:
    input: biom
    output: qiimeoutdir+'/otu_table.txt'
    shell: """
           export PATH={qiime1env}/bin/:$PATH
           {qiime1env}/bin/biom convert -i {input} -o {output} --to-tsv --header-key taxonomy
           """

rule phylumtable:
    input: biom
    output: qiimeoutdir+"/phylum/"+phylumfile
    shell: 
      """
      {qiime1env}/bin/summarize_taxa.py -i {input} -L 2 -o {qiimeoutdir}/phylum/
      """

#conda install  r-optparse bioconductor-metagenomeseq  r-biom r-plyr r-RJSONIO bioconductor-rhdf5 bioconductor-biomformat
#cd dependencies && git clone git@github.com:joey711/biom.git
rule normalize:
    input: biom
    output: qiimeoutdir+"/CSS_normalized_otu_table.biom"
    shell: """
           export PATH={qiime1env}/bin/:$PATH
           {qiime1env}/bin/normalize_table.py -i {input} -a CSS -o {output}
           """
           
#generate distance matrices
#we can also try this with non-normalized?
rule diversity:
    input: otu="intermediates/OTU_SILVA_CLOSED/CSS_normalized_otu_table.biom", repset="intermediates/OTU_SILVA_CLOSED/rep_set.tre"
    params: beta="intermediates/OTU_SILVA_CLOSED/beta_div"
    output: "intermediates/OTU_SILVA_CLOSED/beta_div/bray_curtis_CSS_normalized_otu_table.txt", "intermediates/OTU_SILVA_CLOSED/beta_div/unweighted_unifrac_CSS_normalized_otu_table.txt"
    shell: """
           {qiime1env}/bin/beta_diversity.py -i {input.otu} -m bray_curtis,unweighted_unifrac -t {input.repset} -o {params.beta}
           """
           
#generate the Principle Coordinate plot (PCoA) decompositions of the distance matrices
rule qiimepca:
    input: bray=qiimeoutdir+"/beta_div/bray_curtis_CSS_normalized_otu_table.txt",unifrac=qiimeoutdir+"/beta_div/beta_div_coords_unifrac.txt"
    output: bray=qiimeoutdir+"/beta_div/beta_div_coords_bray.txt", unifrac=qiimeoutdir+"/beta_div/beta_div_coords_unifrac.txt"
    shell:
            """
            {qiime1env}/bin/principal_coordinates.py -i {input.bray} -o {output.bray}
            {qiime1env}/bin/principal_coordinates.py -i {input.unifrac} -o {output.unifrac}
            """

rule phyloseqpca:
    input: biom=biom,tree = otutree
    output: "phyloseq.html"
    shell:
        """
        echo "biom<-'{input.biom}';tree<-'{input.tree}';rmarkdown::render('phyloseq.Rmd')" |  R --quiet --no-save
        """

#mamba install r-vegan
#get pvalues for 
rule comparisons:
  input: otu="intermediates/OTU_SILVA_CLOSED/beta_div/unweighted_unifrac_CSS_normalized_otu_table.txt",mapfile="mapfile.txt"
  output: term="intermediates/OTU_SILVA_CLOSED/beta_div/adonis_term",type="intermediates/OTU_SILVA_CLOSED/beta_div/adonis_type"
  shell: """
         {qiime1env}/bin/compare_categories.py --method adonis -i {input.otu} -m {input.mapfile} -c Case_Control -o intermediates/OTU_SILVA_CLOSED/beta_div/adonis_term -n 999
         {qiime1env}/bin/compare_categories.py --method adonis -i {input.otu} -m {input.mapfile} -c Type         -o intermediates/OTU_SILVA_CLOSED/beta_div/adonis_type -n 999
         {qiime1env}/bin/compare_categories.py --method adonis -i {input.otu} -m {input.mapfile} -c Delivery     -o intermediates/OTU_SILVA_CLOSED/beta_div/adonis_delivery -n 999
       """

#emacs /home/jnl47/miniconda3/envs/qiime1/lib/python2.7/site-packages/qiime-1.9.1-py2.7.egg/qiime/make_2d_plots.py
#change axisbg to facecolor to work with matplotlib==2.2.5
rule twodplots:
  input: bray="intermediates/OTU_SILVA_CLOSED/beta_div/beta_div_coords_bray.txt", unifrac="intermediates/OTU_SILVA_CLOSED/beta_div/beta_div_coords_unifrac.txt", mapfile="mapfile.txt"
  output: bray="intermediates/OTU_SILVA_CLOSED/beta_div/brayPCoA/beta_div_coords_bray_2D_PCoA_plots.html", unifrac="intermediates/OTU_SILVA_CLOSED/beta_div/unifracPCoA/beta_div_coords_unifrac_2D_PCoA_plots.html"
  shell:
    """
    {qiime1env}/bin/make_2d_plots.py -i {input.unifrac} -m {input.mapfile} -o intermediates/OTU_SILVA_CLOSED/beta_div/unifracPCoA
    {qiime1env}/bin/make_2d_plots.py -i {input.bray} -m {input.mapfile} -o intermediates/OTU_SILVA_CLOSED/beta_div/brayPCoA
    """

# Metagenomic
# The shotgun metagenomic sequence reads were processed using the Sunbeam pipeline [32].
# Reads were quality controlled with Trimmomatic [48], to remove low-quality bases and adapter sequences.
# Host reads were then identified using BWA [49] and removed.
# The remaining reads were taxonomically classified using Kraken [33] with a custom database built on all bacterial, fungal, archaeal, and viral genomes in RefSeq release 79 [50].
# To further filter our data of problematic, low-complexity reads, we used Komplexity to remove low-complexity reads that remained after host decontamination [32].
# As a final filtering step, we excluded reads that were classified by Kraken as Chordata, Arthropoda, and Apicomplexa. Sequencing reads were aligned to target genomes using BWA-MEM [49],
# and genome coverage was calculated using BEDtools [51].
import pandas as pd
df=pd.read_csv("metadata/SRP141397.metadata",sep="\t")
mgxdf=df[df['library_strategy'].str.contains("WGS")]
mgxdf=mgxdf.assign(pair1='intermediates/mgx/'+mgxdf['run_accession']+'_1.fastq.gz',pair2='intermediates/mgx/'+mgxdf['run_accession']+'_2.fastq.gz')
mgx=mgxdf['pair1'].values.tolist()+mgxdf['pair2'].values.tolist()

#isolate the metagenomics files (from the 16s)
rule movemgx:
    input: "intermediates/fastq/{file}"
    output: "intermediates/mgx/{file}"
    shell:
      """
      mkdir -p intermediates/mgx/
      mv {input} {output}
      """

rule mgx:
    input: mgx

# rule filterMetaReads:
#     input: project.getMetaReads()
#     output: project.getFilteredMetaReads()
#     shell: """
#            echo "Komplexity"
#            """
# 
# rule trimMetaReads:
#     input: project.getMetaReads()
#     
# rule classifyMetaReads:
#     input: project.getFilteredMetaReads()
#     output: "kraken"
#     shell: """
#           echo "kraken""
#           """

