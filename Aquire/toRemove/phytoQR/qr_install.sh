sudo apt-get update
sudo apt-get install imagemagick
sudo apt-get install qrencode
sudo apt-get install libqrencode3
baseFolder = $HOME'/imageStore/'
kvSep='_'
while [ "$addTile" = "y"]
do
	fileData=''
	folderData=''
	echo "Do you want to add file or folder meta data to the QR tile?"
	echo "[1] - folder"
	echo "[2] - file"
	read metaDataChoice
	case $metaDataChoice in
		1)
		echo "Please enter the folder meta data."
		read metaData
		folderData=$folderData'/'$metaData
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
done
qrencode "test1" -o tile1.png
qrencode "test2" -o tile2.png
convert tile1.png tile2.png -append fullTile.png



