Fetch sequence data and metadata
```
mkdir metadata
conda install pysradb
wget -c https://starbuck1.s3.amazonaws.com/sradb/SRAmetadb.sqlite.gz && gunzip SRAmetadb.sqlite.gz
#if you dont use saveto its space delimited
pysradb metadata --db ./SRAmetadb.sqlite SRP141397 --detailed --expand --saveto metadata/SRP141397.metadata;  done
#this is reentrant - very cool
mkdir -p raw
pysradb download --db ./SRAmetadb.sqlite --out-dir ./raw -p SRP141397;  done
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

