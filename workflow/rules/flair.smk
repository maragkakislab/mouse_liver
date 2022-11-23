rule flair_bam2bed:
    input:
        SAMPLES_DIR + "/{sample}/align/reads.toGenome.sorted.bam"
    output:
        FLAIR_RES + "/{sample}/reads.toGenome.sorted.bed"
    resources:
        mem_mb = 30*1024
    envmodules:
        "flair/1.6.1"
    shell:
        """
        bam2Bed12 -i {input} > {output}
        """


rule flair_correct:
    input:
        bed = FLAIR_RES + "/{sample}/reads.toGenome.sorted.bed",
        genome = GENOME_FILE,
        gtf = GTF_FILE,
    output:
        FLAIR_RES + "/{sample}/reads_all_corrected.bed"
    params:
        out_prefix = FLAIR_RES + "/{sample}/reads",
    resources:
        mem_mb = 200*1024,
        runtime = 24*60,
        lscratch = 30
    threads: 16
    envmodules: # NOTE: when using module, use flair.py instead of flair
       "flair/1.6.1"
    shell:
        """
        flair.py correct \
            --threads {threads} \
            --query {input.bed} \
            --gtf {input.gtf} \
            --genome {input.genome} \
            --output {params.out_prefix}
        """


rule flair_concatenate_bed_files:
    input:
        bed = expand(FLAIR_RES + "/{s}/reads_all_corrected.bed", s = samples.keys())
    output:
        FLAIR_RES + "/all/reads_all_corrected.bed"
    resources:
        mem_mb = 20*1024,
        runtime = 12*60,
        lscratch = 30
    threads: 2
    envmodules:
       "bedops/2.4.41"
    shell:
        """
        bedops -u {input.bed} > {output}
        """


def input_reads_for_flair_collapse(samples):
    files = []
    for s in samples.values():
        if s.is_unstranded():
            files.append(os.path.join(
                SAMPLES_DIR, s.name, "fastq", "reads.pychopped.fastq.gz"))

        files.append(os.path.join(
            SAMPLES_DIR, s.name, "fastq", "reads.fastq.gz"))
    return files


rule flair_collapse:
    input:
        bed = FLAIR_RES + "/all/reads_all_corrected.bed",
        reads = lambda ws: input_reads_for_flair_collapse(samples),
        genome = GENOME_FILE,
        gtf = GTF_FILE,
    output:
        bed = FLAIR_RES + "/all/reads.isoforms.bed",
        fa = FLAIR_RES + "/all/reads.isoforms.fa",
        #gtf =  FLAIR_RES + "/all/reads.isoforms.gtf",
    params:
        out_prefix = FLAIR_RES + "/all/reads",
    resources:
        mem_mb = 120*1024,
        runtime = 12*60,
        lscratch = 60
    threads: 40
    envmodules:
       "flair/1.6.1"
    shell:
        """
        flair.py collapse \
            --threads {threads} \
            --query {input.bed} \
            --genome {input.genome} \
            --reads {input.reads} \
            --output {params.out_prefix}
        """


rule flair_make_metadata:
    output:
        FLAIR_RES + "/{sgroup}/metadata.tsv"
    threads: 1
    run:
        with open(output[0], 'w') as out:
            meta_list = FLAIR_METADATA[wildcards.sgroup]
            for l in meta_list:
                out.write("\t".join(l) + '\n')


rule flair_quantify:
    input:
        fa =  FLAIR_RES + "/all/reads.isoforms.fa",
        meta =  FLAIR_RES + "/all/metadata.tsv"
    output:
        tsv = FLAIR_RES + "/all/reads.flair.quantify"
    params:
        out_prefix = FLAIR_RES + "/all/reads.flair.quantify",
    resources:
        mem_mb = 120*1024,
        runtime = 12*60,
        lscratch = 60
    threads: 40
    envmodules:
       "R/4.2.2",
       "flair/1.6.1"
    shell:
        """
        flair.py quantify \
            --reads_manifest {input.meta} \
            --isoforms {input.fa} \
            --threads {threads} \
            --temp_dir temp/ \
            --output {params.out_prefix}
        """


rule flair_diffexp:
    input:
        tsv = FLAIR_RES + "/all/reads.flair.quantify"
    output:
        sentinel = FLAIR_RES + "/all/reads.flair.diffexp"
    params:
        out_prefix = FLAIR_RES + "/all/reads.flair.diffexp",
    resources:
        mem_mb = 120*1024,
        runtime = 24*60,
        lscratch = 60
    threads: 40
    envmodules:
       "R/4.2.2",
       "flair/1.6.1"
    shell:
        """
        flair.py diffexp \
            --counts_matrix {input.tsv} \
            --threads {threads} \
            --out_dir {params.out_prefix}
        touch {output.sentinel}
        """


rule flair_diffsplice:
    input:
        tsv = FLAIR_RES + "/all/reads.flair.quantify",
        bed = FLAIR_RES + "/all/reads.isoforms.bed"
    output:
        FLAIR_RES + "/all/reads.flair.diffsplice.alt3.events.quant.tsv",
        FLAIR_RES + "/all/reads.flair.diffsplice.alt5.events.quant.tsv",
        FLAIR_RES + "/all/reads.flair.diffsplice.es.events.quant.tsv",
        FLAIR_RES + "/all/reads.flair.diffsplice.ir.events.quant.tsv",
    params:
        out_prefix = FLAIR_RES + "/all/reads.flair.diffexp",
    resources:
        mem_mb = 120*1024,
        runtime = 12*60,
        lscratch = 60
    threads: 40
    envmodules:
       "flair/1.6.1"
    shell:
        """
        flair.py diffexp \
            --counts_matrix {input.tsv} \
            --threads {threads} \
            --test \
            --out_dir {params.out_prefix}
        """
