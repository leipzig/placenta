from utils import metautils
import os
project=metautils.srpMeta('SRP141397')
project.process()

rule sra2fastq:
    input: "{file}.sra"
    params: dirname=lambda wildcards: os.path.dirname(wildcards.file+".sra")
    output: "{file}_1.fastq", "{file}_2.fastq"
    shell:
        """
        fastq-dump --split-files {input} --outdir {params.dirname}
        """

rule fastq2fasta:
    input: "filtered/{file}.fastq.gz"
    output: "fasta/{file}.fasta.gz"
    shell:
        """
        gunzip -c {input} | fastq_to_fasta -z -o {output}
        """

rule get16s:
    input: project.get16SFiles(flatten=True,compress=True)

rule getFiltered:
    input: project.getFiltered16SFiles()

rule getFasta:
    input: project.getFastaFiles()
