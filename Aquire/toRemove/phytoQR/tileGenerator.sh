#!/bin/bash
baseFolder=$HOME'/imageStore/'
kvSep='_'
anotherTile='y'
tileCount=1
folderKey='fr'
fileKey='fe'
tileString=''
fullHeight=110
fullWidth=200
tmpFileLocation=$HOME/phytoMorphTK/tmpData/
mkdir -p $tmpFileLocation

# make spacer image for tile strip
blankImage=$tmpFileLocation'white.png'

padSize='100x100'

# make pad image
convert -size $padSize canvas:white $blankImage

echo "Which USB port are you generating a tile array for?"
read USBport

# init command 
initCommand='{command_sessionStart}{device_'$USBport'}'

# echo the string to the strip file
strFile=$tmpFileLocation'tileStripText.txt'
echo $initCommand>>$strFile

# build init tile 
buildTile.sh $initCommand $tileCount $tmpFileLocation $fullWidth $fullHeight
tmpTileName=$tmpFileLocation'tile'$tileCount.png

# add padding image
convert  $tmpTileName $blankImage +append $tmpFileLocation'tileArray.png'

((tileCount++))

# save the end stip for the final close command
endCommand='{command_sessionEnd}{device_'$USBport'}'


while [ "$anotherTile" = "y" ]
do
	fileData=''
	folderData=''
	addTile='y'
	while [ "$addTile" = "y" ]
	do
		echo "Do you want to add file or folder meta data to the QR tile?"
		echo "[1] - folder"
		echo "[2] - file"
		read metaDataChoice
		case $metaDataChoice in
			1)
			echo "Please enter the folder meta data."
			read metaData
			folderData=$folderData$metaData'/'
			;;
			2)
			echo "Please enter the file meta_data key."
			read metaDataKey
			echo "Please enter the file meta_data value."
			read metaDataValue
			fileData=$fileData'{'$metaDataKey$kvSep$metaDataValue'}'
			;;
		esac
		echo "Current folder data is:"
		echo $baseFolder$folderData
		echo "Current file data is:"
		echo $fileData
		echo "Do you want to add more meta-data to this tile? (y/n)"
		read addTile
	done
	
	if [ "$folderData" != '' ]
	then
		echo "Here"
		tileString=$tileString'{'$folderKey$kvSep$folderData'}'
	fi
	
	if [ "$fileData" != '' ]
	then
		tileString=$tileString'{'$fileKey$kvSep$fileData'}'
	fi

	# echo the string to the strip file
	echo $tileString>>$strFile

	# build tile
	#buildTile.sh $tileString $tileCount $tmpFileLocation $fullWidth $fullHeight

	# query user for another
	echo "Do you want to add another tile? (y/n)"
	read anotherTile

	# add tile the command strip
	#convert  $tmpFileLocation'tileArray.png' $tmpFileLocation'tile'$tileCount.png  +append $tmpFileLocation'tileArray.png'

	# add padding image
	#convert  $tmpFileLocation'tileArray.png' $blankImage +append $tmpFileLocation'tileArray.png'

	# remove the if statement due to the 
	#if [ $tileCount -eq 1 ]
	#then
	# 	cp $tmpFileLocation'tile'$tileCount'.png' $tmpFileLocation'tileArray.png'
	#else
	#	convert  $tmpFileLocation'tileArray.png' $tmpFileLocation'tile'$tileCount.png  +append $tmpFileLocation'tileArray.png'
	#fi

	((tileCount++))
done

# echo the string to the strip file
echo $endCommand>>$strFile

# build end tile 
#buildTile.sh $endCommand $tileCount $tmpFileLocation $fullWidth $fullHeight

#tmpTileName=$tmpFileLocation'tile'$tileCount.png

# attach end
#convert  $tmpFileLocation'tileArray.png' $tmpTileName +append $tmpFileLocation'tileArray.png'
#gthumb $tmpFileLocation'tileArray.png'

#cp $tmpFileLocation'tileArray.png' $HOME/'tileStrip_usbPort'$USBport.png

# loop and build from file
tileCount=1
fileList=()
while read tileString; do
	# current file name
	tmpTileName=$tmpFileLocation'tile'$tileCount.png
	# add the tmp image file to the list
	fileList+=( $tmpTileName )
	# build tile
	buildTile.sh $tileString $tileCount $tmpFileLocation $fullWidth $fullHeight
	# increment the tile count	
	((tileCount++))
done <$strFile

# build the tile
convert  ${fileList[@]} +append $tmpFileLocation'tileArray.png'

# remove the tmp data
rm -r $tmpFileLocation


#echo $tmpFileLocation
#for fileName in $tmpFileLocation*.png
#do
#	echo $fileName
#done
