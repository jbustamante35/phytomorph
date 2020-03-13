#!/bin/bash
xargs -P 10 -I % -d '\n' -a $2 ./fastPathSearch.sh $1 %
