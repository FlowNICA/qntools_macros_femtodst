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

export programm_name=run1

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

#!/bin/bash

#
#SBATCH -D /mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/TMP
#SBATCH -J run_aux_femtodst
#SBATCH -p compute
#SBATCH --time=18:30:00
#SBATCH -a 1-100
#
#SBATCH -o /mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/TMP/slurm_%A_%a.out
#SBATCH -e /mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/TMP/slurm_%A_%a.err
#

export SKIPED_TASKS=$1

export programm_name=run0

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

# File list (of filelists) for UrQMD mcpico data at 5 GeV
export FILELIST=/mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/macros/runlists/runlists_femtodst_auau_200gev.list
export SHORTNAME1=`basename $FILELIST`
export SHORTNAME11=${SHORTNAME1%.list}
export SHORTNAME12=${SHORTNAME11#runlists_femtodst_}
export LABEL1=${SHORTNAME12}_aux
#export LABEL2=flow_ch
export LABEL=${LABEL1}_${programm_name} #${LABEL1}_${LABEL2}

# Config file for Qn measurements
export MACRO_EXE=${MAIN_DIR}/RunDataAux.C

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
export OUT_FILE=${OUT_FILEDIR}/aux_${LABEL}_${JOB_ID}_${TASK_ID}_run_${runid}.root
export LOG=${OUT_LOG}/JOB_${JOB_ID}_${TASK_ID}_run_${runid}.log

export TMP_DIR=${MAIN_DIR}/TMP/TMP_${JOB_ID}_${TASK_ID}

# Creating output directories and log file
mkdir -p $TMP_DIR
mkdir -p $OUT_LOG
mkdir -p $OUT_FILEDIR
touch $LOG

# Main process
echo "Job Id:  ${JOB_ID}" &>> $LOG
echo "Task Id: ${TASK_ID}" &>> $LOG
echo "Read Id: ${READ_ID}" &>> $LOG
echo "Run Id:  ${runid}" &>> $LOG
echo "INFILE:  ${IN_FILE}" &>> $LOG
echo "OUTFILE: ${OUT_FILE}" &>> $LOG
echo "MACRO:   ${MACRO_EXE}" &>> $LOG
echo "Aux ROOT-macro contents:" &>> $LOG
echo "------------------------------------------------------------------------------" &>> $LOG
echo "" &>> $LOG
cat $MACRO_EXE &>> $LOG
echo "" &>> $LOG
echo "------------------------------------------------------------------------------" &>> $LOG
echo "" &>> $LOG


cd $MAIN_DIR
source ${MAIN_DIR}/env.sh &>> $LOG
root -l -b -q $MACRO_EXE'("'${IN_FILE}'","'${OUT_FILE}'")' &>> $LOG

rm -rfv $TMP_DIR &>> $LOG

cd $START_DIR
echo "Job is finished" &>> $LOG