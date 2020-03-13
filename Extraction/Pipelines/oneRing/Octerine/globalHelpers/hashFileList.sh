#!/bin/bash
xargs -P 12 -I % -d '\n' -a $1 ./myH.sh 512 % $2
