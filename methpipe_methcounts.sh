#!/bin/bash

# usage: ./scriptname PDxxxxa (sample name)

PD=$1


METHCOUNTS_DIRNAME="/scratch-shared/arjenbr/methpipe/3_methcounts"
DUPREMOVE_DIRNAME="/scratch-shared/arjenbr/methpipe/2_dupremove"
GENOME="/home/arjenbr/rmapbs/hg19_noUn_noRandom_lamba.tar.gz"

workdir=$TMPDIR/rmapbs
mkdir -p $workdir
cd $workdir
cp $GENOME .
tar -xf `basename $GENOME`
GENOME="rmap"
mkdir -p ${METHCOUNTS_DIRNAME}

/home/arjenbr/rmapbs/methpipe-v2.03/bin/methcounts -v -M 40 -c ${GENOME} -S ${METHCOUNTS_DIRNAME}/${PD}_CpG_methcounts_stat.txt -o ${METHCOUNTS_DIRNAME}/${PD}_CpG_methcounts.bed ${DUPREMOVE_DIRNAME}/${PD}_dupremove.mr 2> ${METHCOUNTS_DIRNAME}/${PD}_CpG_methcounts.log

rm -rf $TMPDIR/*
