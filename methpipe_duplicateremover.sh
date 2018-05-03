#!/bin/bash

# usage: ./scriptname PDxxxxa (sample name)

PD=$1


DUPREMOVE_DIRNAME="/scratch-shared/arjenbr/methpipe/2_dupremove"
SORTED_DIRNAME="/scratch-shared/arjenbr/methpipe/1_merged_sorted"

mkdir -p ${DUPREMOVE_DIRNAME}

/home/arjenbr/rmapbs/methpipe-v2.03/bin/duplicate-remover -S ${DUPREMOVE_DIRNAME}/${PD}_dupremove.txt -o ${DUPREMOVE_DIRNAME}/${PD}_dupremove.mr ${SORTED_DIRNAME}/${PD}.mr_sorted


