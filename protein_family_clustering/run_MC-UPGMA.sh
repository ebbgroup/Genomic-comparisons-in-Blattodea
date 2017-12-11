#!/bin/bash
#$ -N MCUPGMA
#$ -o /global/public/home/l_krem02/isoptera_annotation/family_clustering/mc-upgma-polistes/mcupgma_1.0.0/21_arthro_with_polistes/blast_allvsall_mcupgma_qsub.out
#$ -e /global/public/home/l_krem02/isoptera_annotation/family_clustering/mc-upgma-polistes/mcupgma_1.0.0/21_arthro_with_polistes/blast_allvsall_mcupgma_qsub.err
#$ -pe smp 80
#$ -V
#$ -S /bin/bash
#$ -wd /global/public/home/l_krem02/isoptera_annotation/family_clustering/mc-upgma-polistes/mcupgma_1.0.0/21_arthro_with_polistes/
#$ -l hostname=ebbsrv09

# the "#$" lines above are configurations of sun grid engine for execution on a computational cluster.

# x.edges.gz is the gzipped edges input file required for MC-UPGMA. This file is derived from BLAST all vs all results of all species
# proteomes, see the MC-UPGMA manual for details-

mcupgma_1.0.0/scripts/cluster.pl x.edges.gz -max_singleton 320670 -x 100.0 -M 800000000 -H 70 -K 69 -j 69 -r 0 -iterations 10000 -sleep 1 -tree x.mcupgma_tree
