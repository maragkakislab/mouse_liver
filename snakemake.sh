#!/bin/bash

module load snakemake || exit 1

snakemake -pr  --keep-going --latency-wait 120 -s workflow/Snakefile --profile workflow/snakemake_profile