#!/bin/bash

#-Set path to the femtoDst files
INPUT_DIR=$1
#-Name pattern for output lists
NAME_PATTERN=$2

#-Example usage
# . GenerateLists.sh /mnt/pool/rhic/2/nigmatkulov/femtoDst/auau/200gev/ ./runlists/runlist_

#/mnt/pool/nica/5/rhic2/nigmatkulov/femtoDst/auau/200gev/1212/st_physics_12126101_raw_1010001.femtoDst.root

CURRENT_DIR=${PWD}
export MAX_FILES=30

mkdir -p `dirname $NAME_PATTERN`

export ipart=0
export nlines=0
export previousrun=0

# Generate runlists of corresponding runid
find $INPUT_DIR -type f -name '*.femtoDst.root' | sort -Vt _ -k3,3  | while read -r line
do
  export line1=`basename $line`
  export line2=${line1#st_physics_}
  export line3=${line2%%_raw*}
  export runid=$line3
  if (( $previousrun == 0 ))
  then
    previousrun=$runid
  fi
  if (( $previousrun != $runid ))
  then
    ipart=0
    nlines=0
  fi
  export OUTPUT_FILE=${NAME_PATTERN}femtodst_runid_${runid}_part${ipart}.list
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
    echo "Processing runid: ${runid} -> ${OUTPUT_FILE}"
  fi
  export OUTPUT_FILE=${NAME_PATTERN}femtodst_runid_${runid}_part${ipart}.list
  #echo "Processing runid: ${runid} -> ${OUTPUT_FILE}"
  echo $line &>> $OUTPUT_FILE
  previousrun=$runid
  ((nlines++))
done
