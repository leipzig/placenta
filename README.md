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

