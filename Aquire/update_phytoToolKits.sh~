#!/bin/bash
case $1 in
	Streams)
		./helper_update.sh Streams
	;;
	Photo)
		./helper_update.sh Photo
	;;	
	all)
		./update_phytoToolKits.sh install
		./update_phytoToolKits.sh Photo
		./update_phytoToolKits.sh Streams
	;;
	install)
		D=$(date +"%S-%M-%H-%m-%d-%y")
		imkdir -p /iplant/home/nmiller/publicData/phytoTKinstalls/
		echo "**********************************************************"
		echo "Backing-up and Updating case-based install script"
		echo "**********************************************************"
		# update the case-based install script
		icp /iplant/home/nmiller/publicData/phytoTKinstalls/install_phytoMorph_toolKits.sh /iplant/home/nmiller/publicData/phytoTKinstalls/install_phytoMorph_toolKits_$D.sh 
		iput -f /mnt/scratch1/phytomorph_dev/Aquire/install_phytoMorph_toolKits.sh /iplant/home/nmiller/publicData/phytoTKinstalls/
		echo "**********************************************************"
	
	;;
	Meta)
		./helper_update.sh Meta
	;;
	main)
		D=$(date +"%S-%M-%H-%m-%d-%y")
		imkdir -p /iplant/home/nmiller/publicData/phytoTKinstalls/
		echo "**********************************************************"
		echo "Backing-up and Updating case-based install script"
		echo "**********************************************************"
		# update the case-based install script
		icp /iplant/home/nmiller/publicData/phytoTKinstalls/phytomorph /iplant/home/nmiller/publicData/phytomorph_$D.sh 
		iput -f /mnt/scratch1/phytomorph_dev/Aquire/phytomorph /iplant/home/nmiller/publicData/phytoTKinstalls/
		echo "**********************************************************"
	;;
	*)
		./update_phytoToolKits.sh main
		./update_phytoToolKits.sh install
		./update_phytoToolKits.sh Photo
		./update_phytoToolKits.sh Streams
		./helper_update.sh Meta
	;;
esac
