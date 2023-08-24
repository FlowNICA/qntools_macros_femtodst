#!/bin/bash

#
#SBATCH -D /mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/TMP
#SBATCH -J run_PID
#SBATCH -p compute
#SBATCH --time=08:30:00
#SBATCH --cpus-per-task=20
#SBATCH -a 1-1
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

export INPUT_FILE=/mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/OUT/auau_200gev_aux_run0/aux_auau_200gev.root
export OUTPUT_FILE=/mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/OUT/auau_200gev_aux_run0/pidfit_auau_200gev.root
export LOG=/mnt/pool/nica/7/parfenovpeter/Soft/qntools_macros_femtodst/OUT/auau_200gev_aux_run0/pidfit_auau_200gev.log

# Config file for Qn measurements
export MACRO_EXE=${MAIN_DIR}/RunPidFits.C

export TMP_DIR=${MAIN_DIR}/TMP/TMP_${JOB_ID}_${TASK_ID}

# Creating output directories and log file
mkdir -p $TMP_DIR
touch $LOG

echo "" &>> $LOG
echo "Starting PID procedure..." &>> $LOG
echo "JOB ID: ${JOB_ID}" &>> $LOG
echo "" &>> $LOG

cp -v $MACRO_EXE ${TMP_DIR}/RunPidFits.C &>> $LOG
root -l -b -q ${TMP_DIR}/RunPidFits.C++'("'${INPUT_FILE}'","'${OUTPUT_FILE}'",true)' &>> $LOG

rm -rfv $TMP_DIR &>> $LOG

cd $START_DIR
echo "Job is finished" &>> $LOG
