#!/bin/bash
# $1 pattern file
# $2 file to scan
xargs -P 1 -I % -d '\n' -a $1 ./superScanner.sh % $2 
