rule bam2bed: #bedtools/2.30.0
    input: 
        SAMPLES_DIR + "/{sample}/align/{prefix}.toGenome.sorted.bam"
    output:
        ANALYSIS_DIR + "/{sample}/align/{prefix}.toGenome.sorted.bed"
    resources:
        mem_mb=1024*30
    envmodules:
        "bedtools/2.30.0"
    shell: 
        """
        bedtools bamtobed -bed12 -i {input} > {output}
        """

rule flair_correct:
    input:
        bed = ANALYSIS_DIR + "/{sample}/align/{prefix}.toGenome.sorted.bed",
        genome = DATA_DIR + "/" + ASSEMBLY + "/" + "genome/genome.fa",
        gtf= DATA_DIR + "/" + ASSEMBLY + "/genes.gtf",
    output: #_all_corrected.bed
        ANALYSIS_DIR + "/{sample}/flair/{prefix}_all_corrected.bed"
    params:
        out_prefix = ANALYSIS_DIR + "/{sample}/flair/{prefix}",
    resources:
        mem_mb = 1024 * 200,
        runtime = 60 * 12,
        lscratch=30
    threads: 16
    envmodules: #note that when FLAIR is loaded as a module, it must be run as flair.py, not flair on the cmd line
       "flair/1.6.1"
    shell:
        """
        flair.py correct --threads {threads} \
        --query {input.bed} --gtf {input.gtf} --genome {input.genome} \
        --output {params.out_prefix}
        """

def input_reads_flair_collapse(sample, fastq_prefix):
    s = sample

    if s.is_unstranded():
        return os.path.join(
            SAMPLES_DIR, s.name, "fastq", fastq_prefix + ".pychopped.fastq.gz")

    return os.path.join(
        SAMPLES_DIR, s.name, "fastq", fastq_prefix + ".fastq.gz")

rule flair_collapse:
    input:
        bed=ANALYSIS_DIR + "/{sample}/flair/{prefix}_all_corrected.bed",
        reads = lambda ws: input_reads_flair_collapse(samples[ws.sample], ws.prefix),
        genome = DATA_DIR + "/" + ASSEMBLY + "/" + "genome/genome.fa",
    output: 
        ANALYSIS_DIR + "/{sample}/flair/{prefix}_.isoforms.bed" #naming convention of '_.' is needed for flair, since output of flair command will have suffix '_.isoforms.bed'
    params: 
        out_prefix = ANALYSIS_DIR + "/{sample}/flair/{prefix}_",
        reads = SAMPLES_DIR + "/{sample}/fastq/{prefix}.pychopped.fastq.gz"
    resources:
        mem_mb = 1024 * 120,
        runtime = 60*12,
        lscratch=60
    threads: 40
    envmodules: 
       "flair/1.6.1"
    shell:
        """
        flair.py collapse --threads {threads} \
        --query {input.bed} --genome {input.genome} \
        --reads {input.reads} --output {params.out_prefix}
        """
# rule flair_quantify:
#     input:
#         ANALYSIS_DIR + "/{sample}/flair/{prefix}_isoforms.fa"
#     threads: 30
#     resources:
#         mem_mb = 1024 * 150,
#         runtime = 60*12,
#         lscratch=60
#     envmodules: 
#        "flair/1.6.1"

