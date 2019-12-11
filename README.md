# Reproducing Leiby et al "Lack of detection of a human placenta microbiome in samples from preterm and term deliveries"
https://doi.org/10.1186/s40168-018-0575-4

Miniconda
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
#restart shell
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Fetch sequence data and metadata
```
conda create -n sra pysradb
conda activate sra
mkdir metadata
wget -c https://starbuck1.s3.amazonaws.com/sradb/SRAmetadb.sqlite.gz && gunzip SRAmetadb.sqlite.gz
#if you dont use saveto its space delimited
pysradb metadata --db ./SRAmetadb.sqlite SRP141397 --detailed --expand --saveto metadata/SRP141397.metadata;  done
#this is reentrant - very cool
mkdir -p raw
pysradb download --db ./SRAmetadb.sqlite --out-dir ./raw -p SRP141397
```

Fetch supplemental tables
```
curl https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-018-0575-4/MediaObjects/40168_2018_575_MOESM1_ESM.xls > metadata/table1.xls
curl https://static-content.springer.com/esm/art%3A10.1186%2Fs40168-018-0575-4/MediaObjects/40168_2018_575_MOESM2_ESM.pdf > metadata/table2.pdf
```

Make qiime environment
```
conda create -n qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c conda-forge -c bioconda -c anaconda
conda activate qiime1
```

DADA2 for 16s
```
#this conflicts too much with qiime1 to share an environment
conda create -n dada bioconductor-dada2
source activate dada
```

https://benjjneb.github.io/dada2/tutorial.html