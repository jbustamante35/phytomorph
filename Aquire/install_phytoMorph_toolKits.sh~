#!/bin/bash

if [ $# -eq 0 ]
	then
		echo "Install Menu."
		echo "Please select one."
		echo "[1] phytoPhoto tool kit."
		echo "[2] phytoStreams tool kit."
		echo "[3] iCommands."
		echo "[4] Matlab MCR 2017b."
		echo "[5] sqlite database."
		echo "[6] gphoto2."
		echo "[7] Condor File Transfer."
		echo "[8] Install Condor Key."
		echo "[9] plantProfiler"
		echo "[10] JARs for Java" 
		echo "[11] phytoMeta"
		echo "[12] dcraw"
		echo "[13] phytoDAQ"

		read selection

		case $selection in
			1)
				installProgram="phytoPhoto"	
			;;
			2)
				installProgram="phytoStreams"
			;;
			3)
				installProgram="icommands"
			;;
			4)
				installProgram="matlabMCR2017b"
			;;
			5)
				installProgram="sqlite3"
			;;
			6)
				installProgram="gphoto2"
			;;
			7)
				installProgram="condorFileTransfer"
			;;
			8)
				installProgram="copyToken"
			;;
			9)
				installProgram="plantProfiler"
			;;
			10)
				installProgram='JARs'
			;;
			11)
				installProgram='phytoMeta'
			;;
			12)
				installProgram='dcraw'
			;;
			13)
				installProgram='phytoDAQ'
			;;

		esac
else
	installProgram=$1
fi

# make install location
mkdir -p $HOME/phytoMorphTK/

# is registered and first time run?
echo '0'> $HOME/phytoMorphTK/reg.txt
echo '0'> $HOME/phytoMorphTK/freg.txt

# make location for objects
PHYTO_OBJECTS=$HOME'/phytoMorphTK/objects/'
mkdir -p $PHYTO_OBJECTS
mkdir $PHYTO_OBJECTS'/.recycle/'
export "PHYTO_OBJECTS="$PHYTO_OBJECTS
echo "export PHYTO_OBJECTS="$PHYTO_OBJECTS"     #phytomorph - phytomorph objects" >> ~/.bashrc

# make base install directory
mkdir -p ~/phytoMorphTK/baseInstall

# make variable for phytomorph software reg key
PHYTO_REG=$HOME'/phytoMorphTK/baseInstall/reg.msg'

# make variable for phytomorph dictionary
PHYTO_DICT=$HOME'/phytoMorphTK/baseInstall/pmd.json'
mkdir -p $PHYTO_DICT
export "PHYTO_DICT="$PHYTO_DICT
echo "export PHYTO_DICT="$PHYTO_DICT"     #phytomorph - phytmorph dictionary" >> ~/.bashrc

# make temp location
PHYTO_TMP=$HOME"/phytoMorphTK/tmpFiles"
mkdir -p $PHYTO_TMP
export "PHYTO_TMP="$PHYTO_TMP
echo "export PHYTO_TMP="$PHYTO_TMP"     #phytomorph - tmpLocation" >> ~/.bashrc

# make files folder
PHYTO_FILES=$HOME"/phytoMorphTK/files"
mkdir -p $PHYTO_FILES
export "PHYTO_FILES="$PHYTO_FILES
echo "export PHYTO_FILES="$PHYTO_FILES"     #phytomorph - tmpLocation" >> ~/.bashrc

# set the image store base
PHYTO_IMAGE_BASE=$HOME"/imageStore/"
mkdir -p $PHYTO_IMAGE_BASE
echo "export PHYTO_IMAGE_BASE="$PHYTO_IMAGE_BASE"     #phytomorph - image_base" >> ~/.bashrc
export "PHYTO_IMAGE_BASE="$PHYTO_IMAGE_BASE

# set the signature base
PHYTO_SIGNATURE_BASE=$HOME"/imageStore/.signatures/"
mkdir -p $PHYTO_SIGNATURE_BASE
echo "export PHYTO_SIGNATURE_BASE="$PHYTO_SIGNATURE_BASE"     #phytomorph - signature_base" >> ~/.bashrc
export "PHYTO_SIGNATURE_BASE="$PHYTO_SIGNATURE_BASE


# export the current path for this install
PATH='~/phytoMorphTK/:'$PATH
export "PATH="$PATH

# copy the install script into the install base
cp install_phytoMorph_toolKits.sh ~/phytoMorphTK/

# install main script
wget -O ~/phytoMorphTK/phytomorph https://de.cyverse.org/dl/d/20E09C21-361E-4F95-8457-50C864A5CFC7/phytomorph

# run update and autoremove
sudo apt-get -q autoremove
sudo apt-get -q update

case $installProgram in
	phytoPhoto)
		echo "***************************************************************************"
		echo "*			Start Installing phytoPhoto				*"
		echo "***************************************************************************"
		./install_phytoMorph_toolKits.sh gphoto2
		# make the directory for phytoStreams
		mkdir -p ~/phytoMorphTK/phytoPhoto/
		# make the directory for the install file
		mkdir -p ~/phytoMorphTK/installFiles/
		# make the directory for the databases(s)
		mkdir -p ~/phytoMorphTK/dataBaseFiles/
		# get the current zip containing the phytoStreams toolkit
		wget -c -O ~/phytoMorphTK/installFiles/phytoPhoto.zip https://de.cyverse.org/dl/d/5CB594AD-DC65-4E73-A36F-B7FBB76AFBED/phytoPhoto.zip
		# unzip the phytoStreams toolkit
		unzip -o ~/phytoMorphTK/installFiles/phytoPhoto.zip -d ~/phytoMorphTK/phytoPhoto/
		# export the current path for this install
		export PATH=~/phytoMorphTK/phytoPhoto/:$PATH
		# export future shells to include the phytoStreams toolkit
		echo 'export PATH=~/phytoMorphTK/phytoPhoto/:$PATH        #phytomorph' >> ~/.bashrc
		configurePhoto
		#########################
		# call to import ENV vars
		OLDIFS=$IFS 
		IFS=$'\n' 
		for line in $(cat ~/phytoMorphTK/tmpFiles/envVar)
		do 
			eval $line
			export $line
			var=$(echo -n $line | cut -d= -f1)
			value=$(echo -n $line | cut -d= -f2)
			echo '---------------------'
			echo 'export var:'$var
			echo 'export value:'$value
			echo 'line:'$line
			echo '---------------------'
		done
		IFS=$OLDIFS
		#########################
		# call create streams
		echo "PRE-source:"$PHYTOPORTS_DB_FILE
		. ~/.bashrc
		echo "POST-source:"$PHYTOPORTS_DB_FILE
		createPortsDB
		echo "***************************************************************************"
		echo "*			End Installing phytoPhoto				*"
		echo "***************************************************************************"
		;;
	phytoStreams)
		echo "***************************************************************************"
		echo "*			Start Installing phytoStreams				*"
		echo "***************************************************************************"
		./install_phytoMorph_toolKits.sh sqlite3
		# make the directory for phytoStreams
		mkdir -p ~/phytoMorphTK/phytoStreams/
		# make the directory for the install file
		mkdir -p ~/phytoMorphTK/installFiles/
		# make the directory for the databases(s)
		mkdir -p ~/phytoMorphTK/dataBaseFiles/
		# get the current zip containing the phytoStreams toolkit
		wget -c -O ~/phytoMorphTK/installFiles/phytoStreams.zip https://de.cyverse.org/dl/d/129CBA87-46C6-4D67-90C7-CC7331CEF0BF/phytoStreams.zip
		# unzip the phytoStreams toolkit
		unzip -o ~/phytoMorphTK/installFiles/phytoStreams.zip -d ~/phytoMorphTK/phytoStreams/
		# export the current path for this install
		export "PATH=~/phytoMorphTK/phytoStreams/:"$PATH
		# export future shells to include the phytoStreams toolkit
		echo 'export PATH=~/phytoMorphTK/phytoStreams/:$PATH    #phytomorph' >> ~/.bashrc
		configureStreams
		#########################
		# call to import ENV vars
		OLDIFS=$IFS 
		IFS=$'\n' 
		for line in $(cat ~/phytoMorphTK/tmpFiles/envVar)
		do 
			eval $line
			export $line
			var=$(echo -n $line | cut -d= -f1)
			value=$(echo -n $line | cut -d= -f2)
			echo '---------------------'
			echo 'export var:'$var
			echo 'export value:'$value
			echo '---------------------'
			export $var
		done
		IFS=$OLDIFS
		#########################		
		# call create streams
		createStreamsDB
		echo "***************************************************************************"
		echo "*			End Installing phytoStreams				*"
		echo "***************************************************************************"
		;;
	icommands)
		echo "***************************************************************************"
		echo "*			Start Installing icommands				*"
		echo "***************************************************************************"
		# download the icommands DEB package
		wget -O ~/irods_install.deb https://files.renci.org/pub/irods/releases/4.1.10/ubuntu14/irods-icommands-4.1.10-ubuntu14-x86_64.deb
		# run the install
		sudo apt-get -y  -q install ~/irods_install.deb
		echo "***************************************************************************"
		echo "*			End Installing icommands				*"
		echo "***************************************************************************"
		;;
	matlabMCR2017b)
		echo "***************************************************************************"
		echo "*			Start Installing Matlab MCR 2017b			*"
		echo "***************************************************************************"
		# make the folder for the MCR 2017b
		mkdir -p ~/phytoMorphTK/installFiles/matlabMCR2017b/
		# download MCR from matlab
		wget ~/phytoMorphTK/installFiles/matlabMCR2017b.zip https://www.mathworks.com/supportfiles/downloads/R2017b/deployment_files/R2017b/installers/glnxa64/MCR_R2017b_glnxa64_installer.zip
		# unzip the zip file
		unzip -q MCR_R2017b_glnxa64_installer.zip -d ~/phytoMorphTK/installFiles/matlabMCR2017b/
		# run the installer
 		sudo ~/phytoMorphTK/installFiles/matlabMCR2017b/install -agreeToLicense yes -mode silent
		# remove the unzipped
		rm -r -f ~/phytoMorphTK/installFiles/matlabMCR2017b/
		# move files
		sudo mv /usr/local/MATLAB/MATLAB_Runtime/v93/bin/glnxa64/libcrypto.so.1.0.0 /usr/local/MATLAB/MATLAB_Runtime/v93/bin/glnxa64/libcrypto.so.1.0.0_BK
		sudo mv /usr/local/MATLAB/MATLAB_Runtime/v93//bin/glnxa64/libssl.so.1.0.0 /usr/local/MATLAB/MATLAB_Runtime/v93//bin/glnxa64/libssl.so.1.0.0_BK
		echo "***************************************************************************"
		echo "*			End Installing Matlab MCR 2017b			*"
		echo "***************************************************************************"
		;;
	sqlite3)
		echo "***************************************************************************"
		echo "*			Start Installing sqlite					*"
		echo "***************************************************************************"
		sudo apt-get -f -y -q install sqlite3
		echo "***************************************************************************"
		echo "*			End Installing sqlite					*"
		echo "***************************************************************************"
		;;
	gphoto2)
		echo "***************************************************************************"
		echo "*			Start Installing gphoto					*"
		echo "***************************************************************************"
		sudo apt-get -f -y -q install gphoto2
		echo "***************************************************************************"
		echo "*			End Installing gphoto					*"
		echo "***************************************************************************"
		;;
	condorFileTransfer)
		echo "***************************************************************************"
		echo "*			Start Installing Condor File Transfer			*"
		echo "***************************************************************************"
		export UBUNTU_CODENAME=$(grep UBUNTU_CODENAME /etc/os-release | cut -d'=' -f2)
		# grab a key file - signed by condor team
		wget https://research.cs.wisc.edu/htcondor/ubuntu/HTCondor-Release.gpg.key
		# add get to apt
		sudo apt-key add HTCondor-Release.gpg.key
		#  place condor repo into the computer repo
		sudo sh -c "echo deb http://research.cs.wisc.edu/htcondor/ubuntu/8.9/$UBUNTU_CODENAME $UBUNTU_CODENAME contrib >> /etc/apt/sources.list"
		sudo sh -c "echo deb-src http://research.cs.wisc.edu/htcondor/ubuntu/8.9/$UBUNTU_CODENAME $UBUNTU_CODENAME contrib >> /etc/apt/sources.list"
		# install condor and depend
		sudo apt-get -y -q install libglobus-gss-assist3 htcondor
		# get config files
		wget https://de.cyverse.org/dl/d/FA7254D5-2EAA-4A11-914A-5538C5CD5977/htcondor_file_transfer_config.zip
		sudo unzip htcondor_file_transfer_config.zip -d /etc/condor/config.d/
		sudo mkdir /etc/condor/{tokens.d,passwords.d}

		###########################################################################
		###########################################################################
		#sep="_"
		#echo "Please enter your institution."
		#read institution
		#institution="${institution// /}"
		#echo "Please enter the lab's principal investigator."
		#read pi
		#pi="${pi// /}"
		#echo "Please enter the building name."
		#read buildingName
		#buildingName="${buildingName// /}"
		#echo "Please enter the room number."
		#read roomNumber
		#roomNumber="${roomNumber// /}"
		#RANDOM_ID=$(printf "%05d" $RANDOM)
		#uniqueID=$institution$sep$pi$sep$buildingName$sep$roomNumber$sep$RANDOM_ID
		#sed -i "s/uwtest0001/$uniqueID/" /etc/condor/config.d/20-Name
		#sed -i "s/TEST/$RANDOM_ID/" /etc/condor/config.d/40-CCB
		#echo "Please email Jason Patton at jpatton@cs.wisc.edu for an authentication token."
		#echo "Please include the unique computer name in the email."
		#echo "Unique computer name:" $uniqueID
		###########################################################################
		###########################################################################


		# enable the condor service
		sudo systemctl enable condor.service
		echo "***************************************************************************"
		echo "*			End Installing Condor File Transfer                     *"
		echo "***************************************************************************"
		;;
	copyToken)
		echo "***************************************************************************"
		echo "*			Start Copy-Token Condor File Transfer                   *"
		echo "***************************************************************************"
		sudo cp ~/Downloads/STARTD_*.key /etc/condor/tokens.d/
		sudo systemctl restart condor.service
		echo "***************************************************************************"
		echo "*			End Copy-Token Condor File Transfer                     *"
		echo "***************************************************************************"
		;;
	plantProfiler)
		echo "***************************************************************************"
		echo "*			Start Install Plant Profiler                            *"
		echo "***************************************************************************"
		# install JARs
		./install_phytoMorph_toolKits.sh JARs3
		# install dcraw
		./install_phytoMorph_toolKits.sh dcraw3
		# install sqlite
		./install_phytoMorph_toolKits.sh sqlite3
		# install the photo tool kit
		./install_phytoMorph_toolKits.sh phytoPhoto
		# install the data stream tool kit
		./install_phytoMorph_toolKits.sh phytoStreams
		# install the icommands
		./install_phytoMorph_toolKits.sh icommands
		# install MCR
		./install_phytoMorph_toolKits.sh matlabMCR2017b
		mkdir -p ~/phytoMorphTK/plantProfiler/
		wget -O ~/phytoMorphTK/plantProfiler/plantProfiler https://de.cyverse.org/dl/d/28B9AEBE-2730-40DA-97D7-7820E707CCCB/plantProfiler
		wget -O ~/phytoMorphTK/plantProfiler/run_plantProfiler https://de.cyverse.org/dl/d/25C38603-11C2-4CA6-9926-64F878E53220/run_plantProfiler.sh
		chmod +x ~/phytoMorphTK/plantProfiler/run_plantProfiler
		chmod +x ~/phytoMorphTK/plantProfiler/plantProfiler

		#####################################
		# add pathing
		PATH='~/phytoMorphTK/:'~/phytoMorphTK/$packageName/
		export "PATH="$PATH
		echo 'export PATH=~/phytoMorphTK/'$packageName'/:$PATH' >> ~/.bashrc
		#####################################

		;;
	JARs)
		mkdir -p ~/phytoMorphTK/JARs/
		wget -O ~/phytoMorphTK/JARs/core-3.2.1.jar https://de.cyverse.org/dl/d/187AA822-2A3D-468C-BA30-69389C641247/core-3.2.1.jar	
		wget -O ~/phytoMorphTK/JARs/javase-3.2.1.jar https://de.cyverse.org/dl/d/D8436F43-220C-40F2-B9DE-F0CADE39C796/javase-3.2.1.jar
		echo "***************************************************************************"
		echo "*			End Install Plant Profiler                              *"
		echo "***************************************************************************"
		;;
	phytoMeta)
		echo "***************************************************************************"
		echo "*			Start Install phytoMeta                                 *"
		echo "***************************************************************************"
		sudo apt-get update
		sudo apt-get -y  -q install imagemagick
		sudo apt-get -y  -q install qrencode
		sudo apt-get -y  -q install libqrencode3
		packageName=phytoMeta
		# make the directory for phytoStreams
		mkdir -p ~/phytoMorphTK/$packageName/
		# make the directory for the install file
		mkdir -p ~/phytoMorphTK/installFiles/
		# get the current zip containing the phytoStreams toolkit
		wget -c -O ~/phytoMorphTK/installFiles/$packageName.zip https://de.cyverse.org/dl/d/1411C0CD-A6B3-41DC-B691-23413415B360/phytoMeta.zip
		# unzip the phytoStreams toolkit
		unzip -o ~/phytoMorphTK/installFiles/$packageName.zip -d ~/phytoMorphTK/$packageName/
		# export the current path for this install
		PATH=~/phytoMorphTK/$packageName/:$PATH
		export "PATH=$PATH"
		# export future shells to include the phytoQR toolkit
		echo 'export PATH=~/phytoMorphTK/'$packageName'/:$PATH' >> ~/.bashrc
		#configureQR
		echo "***************************************************************************"
		echo "*			End Install phytoMeta                                   *"
		echo "***************************************************************************"
		;;
	dcraw)
		echo "***************************************************************************"
		echo "*			Start Install dcraw                                     *"
		echo "***************************************************************************"
		sudo  -q -y apt-get install dcraw
		echo "***************************************************************************"
		echo "*			End Install dcraw					*"
		echo "***************************************************************************"
	;;
	phytoDAQ)
		echo "***************************************************************************"
		echo "*			Start Install phytoDAQ                                  *"
		echo "***************************************************************************"
		# install the icommands
		./install_phytoMorph_toolKits.sh phytoMeta
		# install JARs
		./install_phytoMorph_toolKits.sh JARs3
		# install dcraw
		./install_phytoMorph_toolKits.sh dcraw3
		# install sqlite
		./install_phytoMorph_toolKits.sh sqlite3
		# install the photo tool kit
		./install_phytoMorph_toolKits.sh phytoPhoto
		# install the data stream tool kit
		./install_phytoMorph_toolKits.sh phytoStreams
		# install the icommands
		./install_phytoMorph_toolKits.sh icommands
		# install MCR
		./install_phytoMorph_toolKits.sh matlabMCR2017b
		# install phytoPhoto
		./install_phytoMorph_toolKits.sh phytoPhoto
		#####################################
		# init the package name
		packageName='deviceBank'
		mkdir -p ~/phytoMorphTK/$packageName/
		wget -O ~/phytoMorphTK/$packageName/$packageName https://de.cyverse.org/dl/d/F1295C02-3887-4FE1-881C-219396D50786/deviceBank
		wget -O ~/phytoMorphTK/$packageName/run_$packageName https://de.cyverse.org/dl/d/6B2FA4A1-4D3F-4E60-BA47-3B71B85CA9A6/run_deviceBank.sh
		chmod +x ~/phytoMorphTK/$packageName/run_$packageName
		chmod +x ~/phytoMorphTK/$packageName/$packageName
		#####################################
		# add pathing
		PATH='~/phytoMorphTK/:'~/phytoMorphTK/$packageName/
		export "PATH="$PATH
		echo 'export PATH=~/phytoMorphTK/'$packageName'/:$PATH' >> ~/.bashrc
		#####################################
		echo "***************************************************************************"
		echo "*			End Install phytoDAQ                                    *"
		echo "***************************************************************************"
		;;
	all)
		echo "***************************************************************************"
		echo "*			Start Install All                                       *"
		echo "***************************************************************************"
		./install_phytoMorph_toolKits.sh phytoPhoto
		./install_phytoMorph_toolKits.sh phytoStreams
		./install_phytoMorph_toolKits.sh icommands
		./install_phytoMorph_toolKits.sh matlabMCR2017b
		./install_phytoMorph_toolKits.sh sqlite3
		./install_phytoMorph_toolKits.sh gphoto2
		./install_phytoMorph_toolKits.sh condorFileTransfer
		./install_phytoMorph_toolKits.sh copyToken
		./install_phytoMorph_toolKits.sh plantProfiler
		./install_phytoMorph_toolKits.sh JARs
		./install_phytoMorph_toolKits.sh phytoMeta
		./install_phytoMorph_toolKits.sh dcraw
		./install_phytoMorph_toolKits.sh phytoDAQ
		echo "***************************************************************************"
		echo "*			End Install All					*"
		echo "***************************************************************************"
	;;	
	*)
		echo "You did not make a proper selection."
		;;
esac

#resource
. ~/.bashrc
