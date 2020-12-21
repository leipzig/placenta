from utils import metautils
from snakemake.utils import R
import os
project=metautils.srpMeta('SRP141397')
project.process()

#use S3 bucket with fastqs ready to go
useAWS=True

# Sequence data was processed using DADA2 [42]. Reads were trimmed from 251 bases to 240. Dereplication, error modeling, denoising, pair merging, and chimera removal were
# performed using default parameters. Taxonomic assignments were generated by comparison to the Silva reference database [43,44,45]. Samples with fewer than 100 reads were
# excluded from downstream analysis. Beta diversity analysis (unweighted UniFrac and weighted UniFrac) was performed using the beta_diversity.py script from Qiime 1.9.1 [46].
# Jaccard, PCoA, and PERMANOVA analysis were performed using functions vegdist, pcoa, and adonis from the R package “vegan” [47].

# To go from SRA do this first
# pysradb download  --out-dir ./raw -p SRP141397
rule sra2fastq:
    output: project.get16SFiles(flatten=True,compress=True)
    run:
        """
        if useAWS:
          parallel-fastq-dump --threads 4 --outdir intermediates/fastq --split-files --tmpdir /tmp --gzip -s `find raw/ -name "*sra"`
        else:
          pysradb download  --out-dir ./raw -p SRP141397
          aws s3 cp --recursive s3://placenta/intermediates/fastq intermediates/fastq
        """

rule fastq2fasta:
    input: "intermediates/filtered/{file}.fastq.gz"
    output: "intermediates/fasta/{file}.fasta.gz"
    shell:
        """
        gunzip -c {input} | fastq_to_fasta -z -o {output}
        """

rule fasta2qiimedef:
    input: pair1="intermediates/fasta/{code}_16S_F_filt.fasta.gz", pair2="intermediates/fasta/{code}_16S_R_filt.fasta.gz"
    output: "intermediates/qiime_ready/{code}_16S.fa"
    shell: """
           gunzip -c {input.pair1} {input.pair2} | awk '/^>/{{print ">{wildcards.code}_" ++i; next}}{{print}}' > {output}
           """


rule getFiltered:
    input: project.get16SFiles(flatten=True,compress=True)
    output: project.getFiltered16SFiles()
    shell: """
           R("source('dataFilterTrim.R'"))
           """

rule getFasta:
    input: project.getFastaFiles()

rule getQiimeReady:
    input: project.getQiimeReady()

rule combine:
    input: project.getQiimeReady()
    output: "seqs.fna"
    shell: "cat intermediates/qiime_ready/*fa > intermediates/seqs.fna"


# 2 hours
rule OTU:
    input: "seqs.fna"
    output: "intermediates/OTU_SILVA/otu_table_mc2_w_tax.biom"
    shell: """
           pick_open_reference_otus.py -i intermediates/seqs.fna -r SILVA_119_QIIME_release/rep_set/97/Silva_119_rep_set97.fna -o OTU_SILVA  -s 0.1 -m usearch61 -p params.txt
           """

rule humanreadableotutable:
    input: "intermediates/OTU_SILVA/otu_table_mc2_w_tax.biom"
    output: "intermediates/OTU_SILVA/otu_table.txt"
    shell: """
           biom convert -i {input} -o {output} --to-tsv --header-key taxonomy
           """
rule phylumtable:
    input: "intermediates/OTU_SILVA/otu_table_mc2_w_tax.biom"
    output: "intermediates/OTU_SILVA/phylum/otu_table_mc2_w_tax_L2.txt"
    shell: 
      """
      summarize_taxa.py -i {input} -L 2 -o intermediates/OTU_SILVA/phylum/
      """

#conda install  r-optparse bioconductor-metagenomeseq  r-biom r-plyr r-RJSONIO bioconductor-rhdf5 bioconductor-biomformat
#cd dependencies && git clone git@github.com:joey711/biom.git
rule normalize:
    input: "intermediates/OTU_SILVA/otu_table_mc2_w_tax.biom"
    output: "intermediates/OTU_SILVA/CSS_normalized_otu_table.biom"
    shell: """
           normalize_table.py -i {input} -a CSS -o {output}
           """
           
#generate distance matrices
#we can also try this with non-normalized?
rule diversity:
    input: otu="intermediates/OTU_SILVA/CSS_normalized_otu_table.biom", repset="intermediates/OTU_SILVA/rep_set.tre"
    params: beta="intermediates/OTU_SILVA/beta_div"
    output: "intermediates/OTU_SILVA/beta_div/bray_curtis_CSS_normalized_otu_table.txt", "intermediates/OTU_SILVA/beta_div/unweighted_unifrac_CSS_normalized_otu_table.txt"
    shell: """
           beta_diversity.py -i {input.otu} -m bray_curtis,unweighted_unifrac -t {input.repset} -o {params.beta}
           """
           
#generate the Principle Coordinate plot (PCoA) decompositions of the distance matrices
rule pca:
    input: bray="intermediates/OTU_SILVA/beta_div/bray_curtis_CSS_normalized_otu_table.txt",unifrac="intermediates/OTU_SILVA/beta_div/beta_div_coords_unifrac.txt"
    output: bray="intermediates/OTU_SILVA/beta_div/beta_div_coords_bray.txt", unifrac="intermediates/OTU_SILVA/beta_div/beta_div_coords_unifrac.txt"
    shell:
            """
            principal_coordinates.py -i {input.bray} -o {output.bray}
            principal_coordinates.py -i {input.unifrac} -o {output.unifrac}
            """


#mamba install r-vegan
#get pvalues for 
rule comparisons:
  input: otu="intermediates/OTU_SILVA/beta_div/unweighted_unifrac_CSS_normalized_otu_table.txt",mapfile="mapfile.txt"
  output: term="intermediates/OTU_SILVA/beta_div/adonis_term",type="intermediates/OTU_SILVA/beta_div/adonis_type"
  shell: """
         compare_categories.py --method adonis -i {input.otu} -m {input.mapfile} -c Case_Control -o intermediates/OTU_SILVA/beta_div/adonis_term -n 999
         compare_categories.py --method adonis -i {input.otu} -m {input.mapfile} -c Type         -o intermediates/OTU_SILVA/beta_div/adonis_type -n 999
         compare_categories.py --method adonis -i {input.otu} -m {input.mapfile} -c Delivery     -o intermediates/OTU_SILVA/beta_div/adonis_delivery -n 999
       """

#emacs /home/jnl47/miniconda3/envs/qiime1/lib/python2.7/site-packages/qiime-1.9.1-py2.7.egg/qiime/make_2d_plots.py
#change axisbg to facecolor to work with matplotlib==2.2.5
rule twodplots:
  input: bray="intermediates/OTU_SILVA/beta_div/beta_div_coords_bray.txt", unifrac="intermediates/OTU_SILVA/beta_div/beta_div_coords_unifrac.txt", mapfile="mapfile.txt"
  output: bray="intermediates/OTU_SILVA/beta_div/brayPCoA/beta_div_coords_bray_2D_PCoA_plots.html", unifrac="intermediates/OTU_SILVA/beta_div/unifracPCoA/beta_div_coords_unifrac_2D_PCoA_plots.html"
  shell:
    """
    make_2d_plots.py -i {input.unifrac} -m {input.mapfile} -o intermediates/OTU_SILVA/beta_div/unifracPCoA
    make_2d_plots.py -i {input.bray} -m {input.mapfile} -o intermediates/OTU_SILVA/beta_div/brayPCoA
    """

# Metagenomic
# The shotgun metagenomic sequence reads were processed using the Sunbeam pipeline [32].
# Reads were quality controlled with Trimmomatic [48], to remove low-quality bases and adapter sequences.
# Host reads were then identified using BWA [49] and removed.
# The remaining reads were taxonomically classified using Kraken [33] with a custom database built on all bacterial, fungal, archaeal, and viral genomes in RefSeq release 79 [50].
# To further filter our data of problematic, low-complexity reads, we used Komplexity to remove low-complexity reads that remained after host decontamination [32].
# As a final filtering step, we excluded reads that were classified by Kraken as Chordata, Arthropoda, and Apicomplexa. Sequencing reads were aligned to target genomes using BWA-MEM [49],
# and genome coverage was calculated using BEDtools [51].


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

