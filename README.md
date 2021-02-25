# Reproducing Leiby et al "Lack of detection of a human placenta microbiome in samples from preterm and term deliveries"
https://doi.org/10.1186/s40168-018-0575-4


Quickstart:
```
docker run -it leipzig/placenta /bin/bash
snakemake --cores all comparisons        #generate comparisons by delivery type, term, and 
```

This analysis has been refactored in order to remain as faithful to the original implementation as described in Leiby et al. I have tried to make logical or default choices when detail was limited in the methods section.

### Reviewers are free to alter whatever they choose in order to assess the validity of this study. Some suggestions include, but are not limited to...
Possible decision points to alter in a test of robustness:
- dataFilterTrim.R contains a number of dada2 parameters which can affect downstream results
- try a more SILVA freeze (more recent that 119)

Possible major refactors:
- upgrade the workflow to use Qiime2
- use Kraken2 for the entire 16s analysis

# Data availability
Raw and intermediate data is available here at s3://placenta

├── intermediates
│   ├── OTU_SILVA
│   │   ├── CSS_normalized_otu_table.biom
│   │   ├── beta_div
│   │   │   ├── beta_div_coords.txt
│   │   │   ├── beta_div_coords_bray.txt
│   │   │   ├── beta_div_coords_unifrac.txt
│   │   │   ├── brayPCoA
│   │   │   │   ├── beta_div_coords_bray_2D_PCoA_plots.html
│   │   │   │   ├── js
│   │   │   │   │   └── overlib.js
│   │   │   │   └── ...
│   │   │   ├── bray_curtis_CSS_normalized_otu_table.txt
│   │   │   ├── unifracPCoA
│   │   │   │   ├── beta_div_coords_unifrac_2D_PCoA_plots.html
│   │   │   │   ├── js
│   │   │   │   │   └── overlib.js
│   │   │   │   └── ...
│   │   │   └── unweighted_unifrac_CSS_normalized_otu_table.txt
│   │   ├── final_otu_map.txt
│   │   ├── final_otu_map_mc2.txt
│   │   ├── index.html
│   │   ├── log_20201218171054.txt
│   │   ├── new_refseqs.fna
│   │   ├── otu_table.txt
│   │   ├── otu_table_mc2.biom
│   │   ├── otu_table_mc2_w_tax.biom
│   │   ├── otu_table_mc2_w_tax_no_pynast_failures.biom
│   │   ├── phylum
│   │   │   ├── otu_table_mc2_w_tax_L2.biom
│   │   │   └── otu_table_mc2_w_tax_L2.txt
│   │   ├── phylum_charts
│   │   │   ├── css
│   │   │   │   └── qiime_style.css
│   │   │   ├── js
│   │   │   │   └── overlib.js
│   │   │   └── raw_data
│   │   │       └── otu_table_mc2_w_tax_L2.txt
│   │   ├── pynast_aligned_seqs
│   │   │   ├── rep_set_aligned.fasta
│   │   │   ├── rep_set_aligned_pfiltered.fasta
│   │   │   ├── rep_set_failures.fasta
│   │   │   └── rep_set_log.txt
│   │   ├── rep_set.fna
│   │   ├── rep_set.tre
│   │   ├── step1_otus
│   │   │   ├── failures.fasta
│   │   │   ├── seqs_clusters.uc
│   │   │   ├── seqs_failures.txt
│   │   │   ├── seqs_otus.log
│   │   │   ├── seqs_otus.txt
│   │   │   └── step1_rep_set.fna
│   │   ├── step2_otus
│   │   │   ├── step2_rep_set.fna
│   │   │   ├── subsampled_failures.fasta
│   │   │   ├── subsampled_failures_clusters.uc
│   │   │   ├── subsampled_failures_otus.log
│   │   │   └── subsampled_failures_otus.txt
│   │   ├── step3_otus
│   │   │   ├── failures_clusters.uc
│   │   │   ├── failures_failures.fasta
│   │   │   ├── failures_failures.txt
│   │   │   ├── failures_otus.log
│   │   │   └── failures_otus.txt
│   │   ├── step4_otus
│   │   │   ├── failures_failures_clusters.uc
│   │   │   ├── failures_failures_otus.log
│   │   │   ├── failures_failures_otus.txt
│   │   │   └── step4_rep_set.fna
│   │   └── uclust_assigned_taxonomy
│   │       ├── rep_set_tax_assignments.log
│   │       └── rep_set_tax_assignments.txt
│   ├── fasta
│   │   ├── AS.B1_16S_F_filt.fasta.gz
...
│   ├── fastq
│   │   ├── SRR7048659_1.fastq.gz
...
│   ├── filtered
│   │   ├── AS.B1_16S_F_filt.fastq.gz
│   │   ├── AS.B1_16S_R_filt.fastq.gz
...
│   └── sunbeam_databases
│       └── krakendb
├── raw
│   └── SRP141397
│       ├── SRX3980107
│       │   └── SRR7049033.sra
...
└── sunbeam_databases
    └── krakendb
        ├── accmap_file.tmp
        ├── archaea
        │   ├── assembly_summary.txt
        │   ├── library.fna
        │   ├── library.fna.dusted
        │   ├── manifest.txt
        │   ├── prelim_map.txt
        │   └── rsync.err
        ├── bacteria
        │   ├── assembly_summary.txt
        │   ├── library.fna
        │   ├── library.fna.DEF0dea2
        │   ├── library.fna.dusted
        │   ├── manifest.txt
        │   ├── prelim_map.txt
        │   └── rsync.err
        ├── database.idx
        ├── database.jdb
        ├── database.jdb.big
        ├── database.kdb
        ├── database_0
        ├── database_1
        ├── database_10
        ├── database_11
        ├── database_12
        ├── database_13
        ├── database_14
        ├── database_15
        ├── database_2
        ├── database_3
        ├── database_4
        ├── database_5
        ├── database_6
        ├── database_7
        ├── database_8
        ├── database_9
        ├── human
        │   ├── assembly_summary.txt
        │   ├── library.fna
        │   ├── library.fna.dusted
        │   ├── manifest.txt
        │   ├── prelim_map.txt
        │   └── rsync.err
        ├── lca.complete
        ├── library
        │   └── added
        │       ├── TkNxTaeaEP.fna
        │       ├── UMUWy2OumF.fna
        │       ├── gmu7yUuoWB.fna
        │       ├── h1c1tKkdtz.fna
        │       ├── ntTrAq8P8N.fna
        │       └── prelim_map.txt
        ├── plasmid
        │   ├── library.fna
        │   ├── library.fna.dusted
        │   ├── manifest.txt
        │   └── prelim_map.txt
        ├── seqid2taxid.map
        ├── taxonomy
        │   ├── accmap.dlflag
        │   ├── citations.dmp
        │   ├── delnodes.dmp
        │   ├── division.dmp
        │   ├── gc.prt
        │   ├── gencode.dmp
        │   ├── merged.dmp
        │   ├── names.dmp
        │   ├── nodes.dmp
        │   ├── nucl_gb.accession2taxid
        │   ├── nucl_wgs.accession2taxid
        │   ├── prelim_map.txt
        │   ├── readme.txt
        │   ├── taxdump.dlflag
        │   ├── taxdump.tar.gz
        │   └── taxdump.untarflag
        ├── unmapped.txt
        └── viral
            ├── assembly_summary.txt
            ├── library.fna
            ├── library.fna.dusted
            ├── manifest.txt
            ├── prelim_map.txt
            └── rsync.err
