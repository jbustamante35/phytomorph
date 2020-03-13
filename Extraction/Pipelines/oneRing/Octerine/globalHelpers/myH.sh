#!/bin/bash
v=$(head -c$1 "$2")
hash=$(echo -n $v | md5sum )
hash=${hash/-/$2}
#echo $hash
#echo $3
echo $hash >> $3
