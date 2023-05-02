#!/bin/bash

#-Set path to the femtoDst files
INPUT_DIR=$1
#-Name pattern for output lists
OUT_DIR=$2

#-Example usage
# . MakeQaHadd.sh /mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/OUT/auau_200gev_run2/ ./runlists/

#/mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/OUT/auau_200gev_run2/388592/qa/qa_auau_200gev_run2_388592_1_run_12126101_part0.root
export START_SCRIPT=/mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/macros/start_qa_hadd.sh.template

CURRENT_DIR=${PWD}
export MAX_FILES=10

mkdir -p $OUT_DIR/lists
mkdir -p $OUT_DIR/qas

export ipart=0
export nlines=0

# Generate runlists of corresponding runid
find $INPUT_DIR -type f -name 'qa_*.root' | while read -r line
do
  export OUTPUT_FILE=${OUT_DIR}/lists/qa_part${ipart}.list
  if [[ ! -f $OUTPUT_FILE ]]
  then
    touch $OUTPUT_FILE
  fi
  if (( $nlines >= $MAX_FILES ))
  then
    ((ipart++))
    nlines=0
  fi
  if (( $nlines == 0 ))
  then
    echo "Processing part ${ipart} with line: ${line} -> ${OUTPUT_FILE}"
  fi
  export OUTPUT_FILE=${OUT_DIR}/lists/qa_part${ipart}.list
  echo $line &>> $OUTPUT_FILE
  ((nlines++))
done

ls ${OUT_DIR}/lists/*.list &> ${OUT_DIR}/run.list

export njobs=`cat ${OUT_DIR}/run.list | wc -l`

sed\
        -e 's!@NJOBS@!'$njobs'!g'\
        -e 's!@FILEDIR@!'$OUT_DIR'!g'\
        $START_SCRIPT > ${OUT_DIR}/start_qa_hadd.sh

#ipart=0
#find $OUT_DIR/lists -type f -name '*.list' | while read -r line
#do
#  echo "Merge file list ${line} to ${OUT_DIR}/files/qa_part${ipart}.root"
#  hadd ${OUT_DIR}/qas/qa_part${ipart}.root @${line}
#  ((ipart++))
#done
