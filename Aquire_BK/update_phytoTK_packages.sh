#!/bin/bash
case $1 in
	phytoStreams)
		./update_phytoStreams.sh
	;;
	installScripts)
		./update_installScripts.sh
	;;
	all)
		./update_phytoTK_packages.sh phytoStreams
		./update_phytoTK_packages.sh installScripts
	;;
	*)

	;;
esac
