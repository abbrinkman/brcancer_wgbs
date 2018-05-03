#!/bin/bash

# this script is for PE mapping of WGBS reads on Cartesius, the Dutch supercomputer cluster

FASTQgz1=$1 # give full path, quoted
FASTQgz2=$2 # give full path, quoted

RMAPBSpe="/home/arjenbr/rmapbs/methpipe-3.0.0/bin/rmapbs-pe"
GENOME="/home/arjenbr/rmapbs/hg19_noUn_noRandom_lamba.tar.gz" # hg19, without "Un" and "random" chromosomes, but including the genome of phage lambda (WGBS spikein for calculation of conversion efficiency)
PREFIXMAPPED="/scratch-shared/arjenbr/mr" # output path
FRAGLENGTH=600
ADAPTER="AGATCGGAAGAGC"
CORES=24

OUTPUT=`basename $FASTQgz1 .fastq.gz`

mkdir -p $PREFIXMAPPED
workdir=$TMPDIR/rmapbs
mkdir -p $workdir
cd $workdir
mkdir read1 read2
cp $FASTQgz1 read1
cp $FASTQgz2 read2
cp $GENOME .
tar -xf `basename $GENOME`

GENOME="rmap"

cd $workdir/read1
FASTQgz=`basename $FASTQgz1`
LENGTH=`zcat $FASTQgz | head -n2 |tail -n1 |wc -m`
LENGTH=$((LENGTH-1)); 

MISMATCHES=$(($LENGTH/10))
lines=`zcat $FASTQgz |wc -l`
Seqs=`echo $lines |awk '{print $1/4}'`
SeqsPerChunk=`echo $Seqs $CORES | awk '{print int(($1/$2)+1)}'`
LinesPerChunk=`echo $SeqsPerChunk | awk '{print $1*4}'`

zcat $FASTQgz |split -l $LinesPerChunk 
cd $workdir/read2
FASTQgz=`basename $FASTQgz2`
zcat $FASTQgz |split -l $LinesPerChunk

ls x* >$workdir/inputlist
cd $workdir

# Read all the files (from a text file, 1 per line) into an array.
IFS=$'\n' read -r -d '' -a files < inputlist


for f in ${files[@]} ; do
  $RMAPBSpe -C $ADAPTER -m $MISMATCHES -c $GENOME -L $FRAGLENGTH -o $f.mr read1/$f read2/$f > $f.out 2>&1 &
done

wait

cat xa*.mr |gzip > $OUTPUT.mr.gz
cat *.out > $OUTPUT.rmapbs.log
cp $OUTPUT.mr.gz $PREFIXMAPPED/
cp $OUTPUT.rmapbs.log $PREFIXMAPPED/

wait 

rm -rf ${TMPDIR}/*

# submit using the following:
# sbatch -t 24:00:00 --job-name=[FASTQgz1] --output=[OUTPUT].slurm.out [scriptname.h] [FASTQgz1] [FASTQgz2]

