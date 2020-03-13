#!/bin/bash
imageFile=$1
newWidth=$2
tmpFileLocation=$3
toDiv=$4
# query the width of the QR
W=$(identify -ping -format '%w' $imageFile)
# query the height of the QR
H=$(identify -ping -format '%h' $imageFile)
# compute side pad width
hW=$(($newWidth-$W))
hW=$((hW / 2))
# make sidePad
convert -size $hW'x'$H canvas:white $tmpFileLocation'tmpPadBottom.png'
# make divider line
convert -size '1x'$H canvas:black $tmpFileLocation'tmpDiv.png'
# attach the divider
#convert $tmpFileLocation'tmpPadBottom.png' $tmpFileLocation'tmpDiv.png' $imageFile $tmpFileLocation'tmpDiv.png' $tmpFileLocation'tmpPadBottom.png' +append $imageFile
convert $tmpFileLocation'tmpPadBottom.png' $imageFile $tmpFileLocation'tmpPadBottom.png' +append $imageFile
