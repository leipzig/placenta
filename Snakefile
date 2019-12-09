from utils import metautils
import os
project=metautils.srpMeta('SRP141397')

rule sra2fastq:
    input: "{file}.sra"
    params: dirname=lambda wildcards: os.path.dirname(wildcards.file+".sra")
    output: "{file}_1.fastq", "{file}_2.fastq"
    shell:
        """
        fastq-dump --split-files {input} --outdir {params.dirname}
        """

rule get16s:
    input: project.get16SFiles()