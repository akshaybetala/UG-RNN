#!/bin/bash

echo "GENERATE A MOLFILE DATASET FROM A CANONICAL SMILE DATASET"
echo "ATHOR ALESSANDRO LUSCI"
#parse input
if [ $# -lt 1 ] 
then
	echo "USAGE: generateDataset.sh datasetDirectory -r"
	exit
fi 
############

#generate dataset
ARGV=( $@ )
ARGC=${#ARGV[@]}
echo "GENERATING DATASET..."
reduceGraph=0
for i in `seq 0 $#`
do
	j=`expr $i + 1`
	if [ "${ARGV[$i]}" = "-r" ]
	then
		reduceGraph=`expr $reduceGraph + 1`
	fi	 
done
if [ $reduceGraph -eq 0 ]
	then
		CMD=`python generateDataset.py -d $1 2>&1`
		if [ $? -ne 0 ]
		then
			echo "ERROR EXITING NOW"
			echo $CMD
			exit
		fi
	fi
if [ $reduceGraph -eq 1 ]
	then
		echo "REDUCE GRAPH ON..."
		CMD=`python generateDataset.py -d $1 -r yes 2>&1`
		if [ $? -ne 0 ]
		then
			echo "ERROR EXITING NOW"
			echo $CMD
			exit
		fi
	fi

#delete carriage returns
echo "DELETING MSDOS CARRIAGE RETURNS..."
sed "s/\r$//" $1"/"$1"Temp" > $1"/"$1
rm $1"/"$1"Temp"
echo "...DONE"
