#!/bin/bash
fn=$(echo $1 | cut -d. -f3) 
xargs -a $1 -P 1 -I % -d '\n' ./fastOrderintoSets.sh % ./pathNames.txt$fn
