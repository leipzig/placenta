FROM rocker/rstudio:4.0.3

MAINTAINER Jeremy Leipzig <leipzig@gmail.com>

#https://github.com/opencontainers/image-spec/blob/master/annotations.md
LABEL org.opencontainers.image.created="2020-12-19T19:39:14UTC"
LABEL org.opencontainers.image.title = "Leiby Placenta Reproduction"
LABEL org.opencontainers.image.description = "This generates the Leiby placenta paper analysis reverse engineered from https://doi.org/10.1186/s40168-018-0575-4"
LABEL org.opencontainers.image.url = "https://github.com/leipzig/placenta"

RUN apt-get update && apt-get install -y build-essential 

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git awscli && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

#https://stackoverflow.com/a/62674910/264696
SHELL ["/bin/bash", "-c"]
RUN conda init bash && conda create -y --name placenta python=3.7
SHELL ["conda", "run", "-n", "placenta", "/bin/bash", "-c"]

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

#main environment
RUN conda install -y mamba
RUN mamba install -y r-base==4.0.3 snakemake rpy2 pysradb bioconductor-dada2 r-dplyr parallel-fastq-dump fastx_toolkit
SHELL ["/bin/bash","-c"]

#qiime1 environment requires Python2.7
RUN conda init bash && conda create -y --name qiime1env python=2.7
SHELL ["conda", "run", "-n", "qiime1env", "/bin/bash", "-c"]
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --add channels anaconda
RUN conda install -y mamba
#RUN mamba install -y gxx_linux-64
RUN mamba install -y matplotlib=1.4.3 mock nose
RUN mamba install -y r-optparse
RUN mamba install -y bioconductor-metagenomeseq
RUN mamba install -y r-plyr r-RJSONIO
RUN mamba install -y bioconductor-rhdf5 bioconductor-biomformat
RUN mamba install -y numpy
RUN mamba install -y pandas==0.24.2
RUN mamba install -y scikit-bio==0.2.3
RUN mamba install -y cogent==1.5.3
RUN pip install --upgrade cython
RUN pip install biom-format==2.1.4
#this had trouble from conda so we install it from source later
#RUN mamba install -y r-biom


WORKDIR "/placenta"
#https://github.com/biocore/qiime/archive/1.9.1.tar.gz
COPY 1.9.1.tar.gz .
RUN gunzip 1.9.1.tar.gz && tar -xvf 1.9.1.tar
# make qiime work with matplotlib 1.4.3
RUN perl -p -i -e  's/axisbg/facecolor/' qiime-1.9.1/tests/test_make_2d_plots.py 
RUN perl -p -i -e  's/axisbg/facecolor/' qiime-1.9.1/scripts/make_2d_plots.py 
RUN perl -p -i -e  's/axisbg/facecolor/' qiime-1.9.1/qiime/make_2d_plots.py


COPY Snakefile .
COPY config.yaml.template config.yaml
COPY runDada.R .
COPY dadaFilterTrim.R .
COPY utils utils
COPY metadata metadata
COPY dependencies dependencies
RUN R CMD INSTALL dependencies/biom

RUN cd qiime-1.9.1 && cp ../dependencies/uclustq1.2.22_i86linux64 ./scripts/uclust && /opt/conda/envs/qiime1env/bin/python setup.py install
RUN rm 1.9.1.tar

RUN chown -R root:staff /placenta
RUN chmod 775 /placenta
RUN echo "setwd('/placenta')" > /home/rstudio/.Rprofile

ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD ["conda", "run", "-n", "placenta", "/bin/bash", "-c","snakemake --cores all"]
