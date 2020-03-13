#!/bin/bash
fileName=$(basename "$1")
baseName=${1/"$fileName"}
echo "$fileName" >> $2
#echo "$fileName" >> $3
#echo "$baseName" >> $2
