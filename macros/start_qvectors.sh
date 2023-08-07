#!/bin/bash

#
#SBATCH -D /mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/TMP
#SBATCH -J run_makeQvectors_femtodst
#SBATCH -p compute
#SBATCH --time=18:30:00
#SBATCH -a 1-100
#
#SBATCH -o /mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/TMP/slurm_%A_%a.out
#SBATCH -e /mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/TMP/slurm_%A_%a.err
#

export SKIPED_TASKS=$1

export programm_name=newpid0

export JOB_ID=${SLURM_ARRAY_JOB_ID}
export TASK_ID=${SLURM_ARRAY_TASK_ID}

if (( $SKIPED_TASKS > 0 ))
then
  READ_ID=$(expr $TASK_ID + $SKIPED_TASKS)
else
  READ_ID=$TASK_ID
fi

# Set up main software via modules
source /mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/env.sh

# Main directory and starting directory to return to after this script is done
export START_DIR=${PWD}
export MAIN_DIR=/mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst

# File with correction calibration info (qa.root)
#export QA_FILELIST=/mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/macros/lists/runlists_qa_auau_200gev.list
#export ORIG_QA_FILE=`sed "${TASK_ID}q;d" $QA_FILELIST`
#export ORIG_QA_FILE=/mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/OUT/auau_200gev_qvec_run1/qa_1.root

# File with PID info (pidfit.root)
export PID_FILE=/mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/OUT/auau_200gev_aux_run1/pidfit_auau_200gev.root

# File list (of filelists) for UrQMD mcpico data at 5 GeV
export FILELIST=/mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/macros/runlists/runlists_femtodst_auau_200gev.list
export SHORTNAME1=`basename $FILELIST`
export SHORTNAME11=${SHORTNAME1%.list}
export SHORTNAME12=${SHORTNAME11#runlists_femtodst_}
export LABEL1=${SHORTNAME12}_qvec
#export LABEL2=flow_ch
export LABEL=${LABEL1}_${programm_name} #${LABEL1}_${LABEL2}

# Config file for Qn measurements
export CONVERT_EXE=${MAIN_DIR}/convertFemto.C
export MACRO_EXE=${MAIN_DIR}/makeQvectors.C

# Setting up main paths/filenames for output
export IN_FILE=`sed "${READ_ID}q;d" $FILELIST`
export line1=`basename $IN_FILE`
export line2=${line1#runlist_femtodst_runid_}
export line3=${line2%%.list}
export runid=$line3
export OUT_DIR=${MAIN_DIR}/OUT/${LABEL}
export OUT=${OUT_DIR}/${JOB_ID}
export OUT_LOG=${OUT}/log
export OUT_FILEDIR=${OUT}/files
export OUT_QADIR=${OUT}/qa
export QA_FILE=${OUT_QADIR}/qa_${LABEL}_${JOB_ID}_${TASK_ID}_run_${runid}.root
export OUT_FILE=${OUT_FILEDIR}/qn_${LABEL}_${JOB_ID}_${TASK_ID}_run_${runid}.root
export LOG=${OUT_LOG}/JOB_${JOB_ID}_${TASK_ID}_run_${runid}.log

export TMP_DIR=${MAIN_DIR}/TMP/TMP_${JOB_ID}_${TASK_ID}

# Creating output directories and log file
mkdir -p $TMP_DIR
mkdir -p $OUT_LOG
mkdir -p $OUT_QADIR
mkdir -p $OUT_FILEDIR
touch $LOG

# Main process
echo "Job Id:  ${JOB_ID}" &>> $LOG
echo "Task Id: ${TASK_ID}" &>> $LOG
echo "Read Id: ${READ_ID}" &>> $LOG
echo "Run Id:  ${runid}" &>> $LOG
echo "INFILE:  ${IN_FILE}" &>> $LOG
echo "OUTFILE: ${OUT_FILE}" &>> $LOG
echo "ORIG_QA: ${ORIG_QA_FILE}" &>> $LOG
echo "QA:      ${QA_FILE}" &>> $LOG
echo "Reader:  ${CONVERT_EXE}" &>> $LOG
echo "CONFIG:  ${MACRO_EXE}" &>> $LOG
echo "Config ROOT-macro contents:" &>> $LOG
echo "------------------------------------------------------------------------------" &>> $LOG
echo "" &>> $LOG
cat $MACRO_EXE &>> $LOG
echo "" &>> $LOG
echo "------------------------------------------------------------------------------" &>> $LOG
echo "" &>> $LOG


cd $MAIN_DIR
source ${MAIN_DIR}/env.sh &>> $LOG
#rsync -vuzh $MACRO_EXE                 ${TMP_DIR}/makeQvectors.C &>> $LOG
#rsync -vuzh ${MAIN_DIR}/makeQvectors.h ${TMP_DIR}/makeQvectors.h &>> $LOG
#rsync -vuzh ${MAIN_DIR}/utils.h        ${TMP_DIR}/utils.h &>> $LOG
root -l -b -q $CONVERT_EXE'("'${IN_FILE}'","'${TMP_DIR}/stardata'","'${PID_FILE}'")' &>> $LOG
cd $TMP_DIR

if [ ! -f "$ORIG_QA_FILE" ]
then
	export ORIG_QA_FILE=qa.root
else
  #ln -s $ORIG_QA_FILE ${TMP_DIR}/qa_orig.root &>> $LOG
  rsync -vuzh $ORIG_QA_FILE ${TMP_DIR}/qa.root &>> $LOG
fi

#root -l -b -q ${TMP_DIR}/makeQvectors.C'("'${TMP_DIR}/stardata.root'","'${ORIG_QA_FILE}'","'qn.root'")' &>> $LOG
root -l -b -q $MACRO_EXE'("'${TMP_DIR}/stardata.list'","'${TMP_DIR}/qa.root'","'qn.root'")' &>> $LOG

mv -v ${TMP_DIR}/qa.root ${QA_FILE} &>> $LOG
mv -v ${TMP_DIR}/qn.root ${OUT_FILE} &>> $LOG
rm -rfv $TMP_DIR &>> $LOG

cd $START_DIR
echo "Job is finished" &>> $LOG
