
#!/bin/bash

#
#$ -wd /scratch2/$USER/TMP/
#$ -cwd
#$ -N run_makeQvectors
# -q all.q
# -l h=(ncx112|ncx115|ncx117|ncx121|ncx124|ncx12[6-7]|ncx130|ncx132|ncx134|ncx136|ncx138|ncx141|ncx144|ncx150|ncx152|ncx159|ncx16[0-9]|ncx17[0-2]|ncx17[4-6]|ncx18[0-1]|ncx18[4-5]|ncx20[1-3]|ncx20[5-8]|ncx21[1-8]|ncx22[4-8]|ncx23[2-8])
#$ -l h_rt=01:30:00
#$ -l s_rt=01:30:00
#$ -t 1-100
#$ -o /dev/null
#$ -e /dev/null
#

# Set up main software via modules
source /cvmfs/nica.jinr.ru/sw/os/login.sh
module add ROOT/v6.26.10-1

# Main directory and starting directory to return to after this script is done
export START_DIR=${PWD}
export MAIN_DIR=/scratch2/parfenov/Soft/qntools_macros_mcpico

# File list (of filelists) for UrQMD mcpico data at 5 GeV
export FILELIST=/scratch2/parfenov/Soft/QCumulant/lists/urqmd_auau_6gev_mcpico/runlists_urqmd_auau_6gev_mcpico_full.list
export SHORTNAME1=`basename $FILELIST`
export SHORTNAME11=${SHORTNAME1%.list}
export SHORTNAME12=${SHORTNAME11#runlists_}
export LABEL1=${SHORTNAME12}
export LABEL=${LABEL1}

# File with correction calibration info (qa.root)
export ORIG_QA_FILE= #

# Config file for Qn measurements
export MACRO_EXE=${MAIN_DIR}/makeQvectors.C

# Setting up main paths/filenames for output
export IN_FILE=`sed "${SGE_TASK_ID}q;d" $FILELIST`
export OUT_DIR=${MAIN_DIR}/OUT/${LABEL}
export OUT=${OUT_DIR}/${JOB_ID}
export OUT_LOG=${OUT}/log
export OUT_FILEDIR=${OUT}/files
export OUT_QADIR=${OUT}/qa
export QA_FILE=${OUT_QADIR}/qa_${LABEL}_${JOB_ID}_${SGE_TASK_ID}.root
export OUT_FILE=${OUT_FILEDIR}/qn_${LABEL}_${JOB_ID}_${SGE_TASK_ID}.root
export LOG=${OUT_LOG}/JOB_${JOB_ID}_${SGE_TASK_ID}.log

export TMP_DIR=${MAIN_DIR}/TMP/TMP_${JOB_ID}_${SGE_TASK_ID}

# Creating output directories and log file
mkdir -p $TMP_DIR
mkdir -p $OUT_LOG
mkdir -p $OUT_QADIR
mkdir -p $OUT_FILEDIR
touch $LOG

# Main process
echo "INFILE:  ${IN_FILE}" &>> $LOG
echo "OUTFILE: ${OUT_FILE}" &>> $LOG
echo "ORIG_QA: ${ORIG_QA_FILE}" &>> $LOG
echo "QA:      ${QA_FILE}" &>> $LOG
echo "CONFIG:  ${MACRO_EXE}" &>> $LOG
echo "Config ROOT-macro contents:" &>> $LOG
echo "------------------------------------------------------------------------------" &>> $LOG
echo "" &>> $LOG
cat $MACRO_EXE &>> $LOG
echo "" &>> $LOG
echo "------------------------------------------------------------------------------" &>> $LOG

if [ -f "$ORIG_QA_FILE" ]
then
  ln -s $ORIG_QA_FILE ${TMP_DIR}/qa.root &>> $LOG
fi

cd $MAIN_DIR
source ${MAIN_DIR}/env.sh &>> $LOG
root -l -b -q $MACRO_EXE'("'${IN_FILE}'","'${TMP_DIR}/qa.root'","'${OUT_FILE}'")' &>> $LOG

rm -rf $TMP_DIR &>> $LOG

cd $START_DIR
echo "Job is finished" &>> $LOG
