#!/bin/bash
# zip up the phtyoPhoto folder, back up the old copy, and place the new copy for download
packageName=$1
D=$(date +"%S-%M-%H-%m-%d-%y")
zip -r -j /mnt/scratch1/phytomorph_dev/Aquire/phyto$packageName.zip /mnt/scratch1/phytomorph_dev/Aquire/phyto$packageName
imkdir -p /iplant/home/nmiller/publicData/phyto$packageName/
imv /iplant/home/nmiller/publicData/phyto$packageName/phyto$packageName.zip /iplant/home/nmiller/publicData/phyto$packageName/phyto$packageName_$D.zip
iput -V -f /mnt/scratch1/phytomorph_dev/Aquire/phyto$packageName.zip /iplant/home/nmiller/publicData/phyto$packageName/
