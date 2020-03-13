#! /bin/bash
D=$(date +"%S-%M-%H-%m-%d-%y")
imkdir -p /iplant/home/nmiller/publicData/phytoTKinstalls/

echo "**********************************************************"
echo "Backing-up and Updating case-based install script"
echo "**********************************************************"
# update the case-based install script
imv /iplant/home/nmiller/publicData/phytoTKinstalls/install_phytoTK.sh /iplant/home/nmiller/publicData/phytoTKinstalls/install_phytoTK_$D.sh 
iput -f /mnt/scratch1/phytomorph_dev/Aquire/install_phytoTK.sh /iplant/home/nmiller/publicData/phytoTKinstalls/
echo "**********************************************************"

echo "**********************************************************"
echo "Backing-up and Updating master Phytomorph install script"
echo "**********************************************************"
# update the case-based install script
imv /iplant/home/nmiller/publicData/phytoTKinstalls/installPhytomorph.sh /iplant/home/nmiller/publicData/phytoTKinstalls/installPhytomorph_$D.sh 
iput -f /mnt/scratch1/phytomorph_dev/Aquire/installPhytomorph.sh /iplant/home/nmiller/publicData/phytoTKinstalls/
echo "**********************************************************"


