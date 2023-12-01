#!/bin/bash
#Written by Joran Vanhollebeke & Isolde Cornelis
#Date: 

clear -x

#Input is asked
echo "Enter your input folder: (vb: /home/genomics/sderycke/Minion_run/basecalled_folder)"
read -e INPUTFOLDER
echo

cd $INPUTFOLDER
echo $INPUTFOLDER
echo

#Input is asked
echo "Enter your RUNNAME: (vb: MiFish_UE-S)"
#read -e RUN_NAME
echo

RUN_NAME=MiFish_UE-S
echo $RUN_NAME

# Define the variable GZFILES
GZFILES=$(find $INPUTFOLDER/${RUN_NAME}*/processed-reads/*/trimmed-R*/ -type f -name "*.gz")

echo Contains >> log.txt
echo $GZFILES >> log.txt

#Unzip all files
gunzip $GZFILES


echo $PWD
FILES=$(cat 20220330_Barcodefile_Samplelocation_tag.txt | awk '{ print $1}')

mkdir -p $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/sense/trimmed-R1/
mkdir -p $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/antisense/trimmed-R1/
mkdir -p $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/sense/trimmed-R2/
mkdir -p $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/antisense/trimmed-R2/

#Loop over all files and do all the commands
for f in $FILES
do 

BASENAME_f=$(basename $f)
#Define the variable SAMPLE who contains the basename where the extension is removed (-1.fastq.gz)

cp $INPUTFOLDER/${RUN_NAME}*/processed-reads/sense/trimmed-R1/${f}_S*.R1.fastq $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/sense/trimmed-R1
gzip $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/sense/trimmed-R1/*.R1.fastq

cp $INPUTFOLDER/${RUN_NAME}*/processed-reads/sense/trimmed-R2/${f}_S*.R2.fastq $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/sense/trimmed-R2
gzip $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/sense/trimmed-R2/*.R2.fastq

cp $INPUTFOLDER/${RUN_NAME}*/processed-reads/antisense/trimmed-R1/${f}_S*.R1.fastq $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/antisense/trimmed-R1
gzip $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/antisense/trimmed-R1/*.R1.fastq

cp $INPUTFOLDER/${RUN_NAME}*/processed-reads/antisense/trimmed-R2/${f}_S*.R2.fastq $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/antisense/trimmed-R2
gzip $INPUTFOLDER/${RUN_NAME}_concatenated/processed-reads/antisense/trimmed-R2/*.R2.fastq

done
  
 echo "finished"
