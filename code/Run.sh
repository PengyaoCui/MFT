#!/bin/bash

##Default valaue
##Generate events and simulate detector
#o2-sim -m MFT -e TGeant3 -g fwmugen -n 10000
#
##Convert simulated hits to detector digits(This will generate mftdigits.root. To skip digitization of a detector, use --skipDet)
#o2-sim-digitizer-workflow -b --skipDet TPC,ITS,TOF,FT0,EMC,HMP,ZDC,TRD,MCH,MID,FDD,PHS,FV0,CPV
#
##The reconstruction worflow requires the cluster dictionary. See how to generate it on the next section.(This will save the mft clusters and mft tracks on mftclusters.root and mfttracks.root, respectively.)
#ln -s ~/alice/MFTdictionary.bin
#o2-mft-reco-workflow -b

o2-sim -n 50000 -g pythia8 -m MFT --configKeyValues=Diamond.width[2]=6.;
o2-sim-digitizer-workflow -b
o2-mft-reco-workflow -b
root -l MFTFitterTrackerChecker.C

