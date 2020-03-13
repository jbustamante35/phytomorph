#!/bin/bash
tileString=$1
tileCount=$2
tmpFileLocation=$3
fullWidth=$4
fullHeight=$5

# build the QR tile
#echo "QR tile string is:"
#echo $tileString
qrencode $tileString -o $tmpFileLocation'tile'$tileCount.png

# query the width of the QR
W=$(identify -ping -format '%w' $tmpFileLocation'tile'$tileCount.png)

# query the height of the QR
H=$(identify -ping -format '%h' $tmpFileLocation'tile'$tileCount.png)

# compute bottom pad height
bottomPadHeight=$(($fullHeight-$H))

# make bottomPad
convert -size $W'x'$bottomPadHeight canvas:white $tmpFileLocation'tmpPadBottom.png'

# make divider line
convert -size $W'x1' canvas:black $tmpFileLocation'tmpDiv.png'

# attach the divider
#convert $tmpFileLocation'tile'$tileCount.png $tmpFileLocation'tmpDiv.png' -append $tmpFileLocation'tile'$tileCount.png

# attach the bottomPad
convert $tmpFileLocation'tile'$tileCount.png $tmpFileLocation'tmpPadBottom.png' -append $tmpFileLocation'tile'$tileCount.png

./padWidth.sh $tmpFileLocation'tile'$tileCount.png $fullWidth $tmpFileLocation 1

# query the width of the QR
W=$(identify -ping -format '%w' $tmpFileLocation'tile'$tileCount.png)

# query the height of the QR
H=$(identify -ping -format '%h' $tmpFileLocation'tile'$tileCount.png)

# build the label
convert -size $W'x'$textHeight -background white -fill black -pointsize 8 -font Helvetica label:"Step$tileCount" -gravity center $tmpFileLocation'tileLabel.png'

# attach the label
convert $tmpFileLocation'tile'$tileCount.png $tmpFileLocation'tileLabel.png' -append $tmpFileLocation'tile'$tileCount.png

#echo $tmpFileLocation'tile'$tileCount.png
