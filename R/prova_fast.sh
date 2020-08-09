#!/bin/bash

#PBS -l select=1:ncpus=8:mem=16g
#PBS -S /bin/bash
#PBS -M stefano.depretis@iit.it
#PBS -m e

cd $PBS_O_WORKDIR # Replacement for -cwd

singularity exec -B /hpcnfs /hpcnfs/data/BA/MP/sdepreti/vm/download-align-count-from-sra.sif \
 fastq-dump ERR173280 -Z | rsem-calculate-expression -p 8 --no-bam-output - GRCh38_transcriptome ERR173280_fast