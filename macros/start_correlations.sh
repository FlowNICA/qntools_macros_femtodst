#!/bin/bash

#
#SBATCH -D /mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/TMP
#SBATCH -J run_QnTools_correlations
#SBATCH -p compute
#SBATCH --time=08:30:00
#SBATCH -a 1-100
#
#SBATCH -o /mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/TMP/slurm_%A_%a.out
#SBATCH -e /mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/TMP/slurm_%A_%a.err
#

export JOB_ID=${SLURM_ARRAY_JOB_ID}
export TASK_ID=${SLURM_ARRAY_TASK_ID}

# Set up main software via modules
source /mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/env.sh

# Main directory and starting directory to return to after this script is done
export START_DIR=${PWD}
export MAIN_DIR=/mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst

# File list (of filelists) for UrQMD mcpico data at 5 GeV
export FILELIST=/mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/macros/lists/runlists_qn_auau_200gev.list
export SHORTNAME1=`basename $FILELIST`
export SHORTNAME11=${SHORTNAME1%.list}
export SHORTNAME12=${SHORTNAME11#runlists_qn_}
export LABEL1=${SHORTNAME12}
#export LABEL2=flow_ch
export LABEL=${LABEL1} #${LABEL1}_${LABEL2}

# Config file for Qn measurements
export MACRO_EXE=${MAIN_DIR}/correlate.C

# Setting up main paths/filenames for output
export IN_FILE=`sed "${TASK_ID}q;d" $FILELIST`
export OUT_DIR=${MAIN_DIR}/OUT/${LABEL}
export OUT=${OUT_DIR}/${JOB_ID}
export OUT_LOG=${OUT}/log
export OUT_FILEDIR=${OUT}/correlations
export OUT_FILE=${OUT_FILEDIR}/qn_${LABEL}_${JOB_ID}_${TASK_ID}.root
export LOG=${OUT_LOG}/JOB_${JOB_ID}_${TASK_ID}.log

#export TMP_DIR=${MAIN_DIR}/TMP/TMP_${JOB_ID}_${TASK_ID}

# Creating output directories and log file
#mkdir -p $TMP_DIR
mkdir -p $OUT_LOG
mkdir -p $OUT_FILEDIR
touch $LOG

# Main process
echo "INFILE:  ${IN_FILE}" &>> $LOG
echo "OUTFILE: ${OUT_FILE}" &>> $LOG
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
root -l -b -q $MACRO_EXE'("'${IN_FILE}'","'${OUT_FILE}'")' &>> $LOG

#rm -rfv $TMP_DIR &>> $LOG

cd $START_DIR
echo "Job is finished" &>> $LOG
