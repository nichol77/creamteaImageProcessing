#!/bin/bash

#if [ "$IND_VAL" = "" ]
#then
#   echo "IND_VAL BASE_NAME and DIR_NAME must be set " 1>&2
#   exit 1
#fi

source /home/wbrigg/.bash_profile
cd /home/wbrigg/imageProcessing/trunk


#for startVal in 5001 6001 7001 8001 9001; do
#let startVal=$IND_VAL-1
#let startVal=100*$startVal
#let startVal=$startVal+1
#let endVal=$startVal+99
#    echo $startVal $endVal
#OUTFILE=${DIR_NAME}/pca/pca_${BASE_NAME}_million_$IND_VAL.root
#echo $OUTFILE
#echo "time ./makePCAFile $DIR_NAME $BASE_NAME 100 $OUTFILE $startVal > log_${BASE_NAME}.txt 2>&1 "
#time ./makePCAFile $DIR_NAME $BASE_NAME 100 $OUTFILE $startVal 
#chmod a+r $OUTFILE 


./GenerateLambdaMap

if [ "$muon_num" = "" ]
then
   ./GenerateLambdaMap
   exit 1
fi
./GenerateLambdaMap $muon_num

