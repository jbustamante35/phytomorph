#!/bin/bash
sep="_"
echo "Please enter your institution."
read institution
institution="${institution// /}"
echo "Please enter the lab's principal investigator."
read pi
pi="${pi// /}"
echo "Please enter the building name."
read buildingName
buildingName="${buildingName// /}"
echo "Please enter the room number."
read roomNumber
roomNumber="${roomNumber// /}"
RANDOM_ID=$(printf "%05d" $RANDOM)
uniqueID=$institution$sep$pi$sep$buildingName$sep$roomNumber$sep$RANDOM_ID
sed -i "s/uwtest0001/$uniqueID/" /etc/condor/config.d/20-Name
echo "Please email Jason Patton at jpatton@cs.wisc.edu for an authentication token."
echo "Please include the unique computer name in the email."
echo "Unique computer name:" $uniqueID

