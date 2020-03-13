#!/bin/bash
# zip up the phtyoStreams folder, back up the old copy, and place the new copy for download
D=$(date +"%S-%M-%H-%m-%d-%y")
zip -r -j /mnt/scratch1/phytomorph_dev/Aquire/phytoStreams.zip /mnt/scratch1/phytomorph_dev/Aquire/phytoStreams
imv /iplant/home/nmiller/publicData/phytoStreams/phytoStreams.zip /iplant/home/nmiller/publicData/phytoStreams/phytoStreams_$D.zip
iput -f /mnt/scratch1/phytomorph_dev/Aquire/phytoStreams.zip /iplant/home/nmiller/publicData/phytoStreams/
