#!/bin/bash

module load snakemake || exit 1

source myconda
mamba activate base
snakemake -pr --keep-going --retries 3 --latency-wait 120 --use-conda --use-envmodules -s workflow/Snakefile --profile workflow/snakemake_profile #correction made
