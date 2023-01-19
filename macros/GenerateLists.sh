#!/bin/bash

#-Set path to the femtoDst files
INPUT_DIR=$1
#-Set the maximum number of files in 1 list
NUM_DIVIDES=$2
#-Name pattern for output lists
NAME_PATTERN=$3

#-Example usage
# . GenerateLists.sh /mnt/pool/rhic/2/nigmatkulov/femtoDst/auau/200gev/12135/ 100

CURRENT_DIR=${PWD}
OUTPUT_DIR="../lists"

#TOTAL_NUM_FILES=`ls $INPUT_DIR/*.femtoDst.root | wc -l`
TOTAL_NUM_FILES=`find $INPUT_DIR -type f -name '*.femtoDst.root' | wc -l`
echo "Total number of DST files: ${TOTAL_NUM_FILES}"

MULT=$((TOTAL_NUM_FILES / NUM_DIVIDES))
RESIDUE=$((TOTAL_NUM_FILES % NUM_DIVIDES))

echo "Mult: $MULT, Residue: $RESIDUE"

for i in `seq 1 $MULT`
do
  OUTPUT_FILE=${NAME_PATTERN}${i}.list
  echo $OUTPUT_FILE
  find $INPUT_DIR -type f -name '*.femtoDst.root' | sort -n | head -n$((i*NUM_DIVIDES)) | tail -n$NUM_DIVIDES &> $OUTPUT_FILE
done

#-Generate last filelist for residual files
OUTPUT_FILE=${NAME_PATTERN}$((i+1)).list
echo $OUTPUT_FILE
find $INPUT_DIR -type f -name '*.femtoDst.root' | sort -n |  tail -n$RESIDUE &> $OUTPUT_FILE
