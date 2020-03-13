function [] = updatePackage()
#!/bin/bash
# zip up the phtyoPhoto folder, back up the old copy, and place the new copy for download
packageName=$1

basePath = '/mnt/scratch1/phytomorph_dev/Deploy/';


CMD = 'date +"%S-%M-%H-%m-%d-%y"';
[~,date] = system(CMD);
date = strtrim(date);



zip -r -j /mnt/scratch1/phytomorph_dev/Aquire/phyto$packageName.zip /mnt/scratch1/phytomorph_dev/Aquire/phyto$packageName


imkdir -p /iplant/home/nmiller/publicData/phyto$packageName/
icp /iplant/home/nmiller/publicData/phyto$packageName/phyto$packageName.zip /iplant/home/nmiller/publicData/phyto$packageName/phyto$packageName_$D.zip
iput -V -f /mnt/scratch1/phytomorph_dev/Aquire/phyto$packageName.zip /iplant/home/nmiller/publicData/phyto$packageName/