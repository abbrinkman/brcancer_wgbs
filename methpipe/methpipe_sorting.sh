#!/bin/bash

# usage: ./scriptname PDxxxxa (sample name)

PD=$1
DIR_MAPPED="/scratch-shared/arjenbr/mr"
DIR_METHPIPE="/scratch-shared/arjenbr/methpipe"
DIRNAME="1_merged_sorted"
SUBDIR_1=${DIR_METHPIPE}/${DIRNAME}
mkdir -p ${SUBDIR_1}
SORTBF="60G"

FILES=${DIR_MAPPED}/${PD}_*.mr.gz
echo ${PD} > ${SUBDIR_1}/${PD}.log
echo "========" >> ${SUBDIR_1}/${PD}.log
for i in ${FILES[@]} ; do
  echo $i >> ${SUBDIR_1}/${PD}.log
done
#merge, take out faulty mappings, and sort
str="zcat ${FILES} | grep -vP '^.{3,}\t[0-9]{10,}\t' | sort -S ${SORTBF} -k 1,1 -k 2,2g -k 6,6 -k 3,3g > ${SUBDIR_1}/${PD}.mr_sorted 2>> ${SUBDIR_1}/${PD}.log" 
eval "${str}"
