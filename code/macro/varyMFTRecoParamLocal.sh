#!/bin/bash
param="$1"

macroDir="/home/cuipengyao/MFT/code/macro"
runningDir="/home/cuipengyao/MFT/code/run"


newfolder="RCut_0.0100"
if [ -z "$param" ]; then
	options="MFTTracking.LTFclsRCut=0.0100;MFTTracking.ROADclsRCut=0.04"
	[ -e dpl-config.json ] && rm dpl-config.json
	[ -e mftclusters.root ] && rm mftclusters.root
	[ -e mfttracks.root ] && rm mfttracks.root
	[ -e o2mftrecoflow_configuration.ini ] && rm o2mftrecoflow_configuration.ini
	[ -e MFTEffCheck.root ] && rm MFTEffCheck.root
	[ -e reco.log ] && rm reco.log
	echo "----> Running with default parameters"
	cd $runningDir
    eval o2-mft-reco-workflow -b --configKeyValues \"$options\" > reco.log
    mkdir $runningDir/$newfolder
	mv reco.log $runningDir/$newfolder
	echo Running efficiency script...
	root -b -q $macroDir/SimuMFTEff.C
	mv MFTEffCheck.root $runningDir/$newfolder
    echo Done!
	exit
else
	if [ "$param" = "PointsLTF" ]; then
		options = "MFTTracking.MinTrackStationsLTF="
		echo "Varying minimum number of points for a LTF track"
		parValues=( 3 4 5 6 7 8 9 10 )
	elif [ "$param" = "StationsLTF" ]; then
		options="MFTTracking.MinTrackPointsLTF="
		echo "Varying the minimum number of detector stations for a LTF track"
		parValues=( 2 3 4 5 6 )
	elif [ "$param" = "PointsCA" ]; then
		options="MFTTracking.MinTrackStationsCA="
		echo "Varying minimum number of points for a CA track"
		parValues=( 3 4 5 6 7 8 9 10 )
	elif [ "$param" = "StationsCA" ]; then
		options="MFTTracking.MinTrackPointsCA="
		echo "Varying the minimum number of detector stations for a CA track"
		parValues=( 2 3 4 5 6 )
	elif [ "$param" = "RCut" ]; then
		options="MFTTracking.LTFclsRCut="
		echo "Varying the maximum distance for a cluster to be attached to a seed line (LTF)"
		parValues=( 0.005 0.0075 0.0100 0.0125 0.0150 0.0200 )
	elif [ "$param" = "ROAD" ]; then
		options="MFTTracking.ROADclsRCut="
		echo  "Varying the maximum distance for a cluster to be attached to a seed line (CA road)"
		parValues=( 0.025 0.030 0.035 0.040 0.045 0.050 0.055 )
	elif [ "$param" = "RBins" ]; then
		options="MFTTracking.RBins="
		echo "Varying the number of bins in r-direction"
		parValues=( 35 40 45 50 55 60 65 )
	elif [ "$param" = "phi" ]; then
		options="MFTTracking.PhiBins="
		echo "Varying the number of bins in phi-direction"
		parValues=( 50 55 60 65 70 75 )
	elif [ "$param" = "LTFseed2BinWin" ]; then
		options="MFTTracking.LTFseed2BinWin="
		echo "Varying the RPhi search window bin width for the second point of a seed (LTF and CA)"
		parValues=( 1 2 3 4 5 )
	elif [ "$param" = "LTFinterBinWin" ]; then
		options="MFTTracking.LTFinterBinWin="
		echo "Varying the RPhi search window bin width for the intermediate points"
		parValues=( 1 2 3 4 5 )
	else
		echo Wrong parameter, exiting.
		exit
	fi
fi


for par in "${parValues[@]}"; do
	echo "------>Running with $param = $par"
	cd $runningDir
	echo Removing previous output
	[ -e dpl-config.json ] && rm  dpl-config.json
	[ -e mftclusters.root ] && rm  mftclusters.root
	[ -e mfttracks.root ] && rm  mfttracks.root
	[ -e o2mftrecoflow_configuration.ini ] && rm  o2mftrecoflow_configuration.ini
	[ -e MFTEffCheck.root ] && rm  MFTEffCheck.root 
	#rm reco*.log
	echo Running reconstruction workflow...
	eval o2-mft-reco-workflow -b --configKeyValues \"$options$par\"  > reco\_$param\_$par\.log
	echo Making new directory and copying files..
	mkdir -p $runningDir/$param\_$par
    mv reco\_$param\_$par\.log $runningDir/$param\_$par
	echo Running efficiency script...
	root -b -q $macroDir/SimuMFTEff.C
	mv MFTEffCheck.root $runningDir/$param\_$par
	echo Done!
done

#exit

#root -q -b $macroDir/SimuMFTEff.C $options
