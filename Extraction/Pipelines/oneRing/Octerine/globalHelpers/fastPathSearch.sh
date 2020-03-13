#!/bin/bash
setNumber=$(echo $1 | cut -d'*' -f1)
pathToMatch=$(echo $1 | cut -d'*' -f2)
fileName=$(basename"$2")
baseName=${2/"$fileName"}
#echo searching:$baseName
#echo matching:$pathToMatch
if [ "$baseName" == "$pathToMatch" ];then
  echo $2 >> $setNumber.txt
fi
