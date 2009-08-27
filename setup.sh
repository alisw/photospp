#!/bin/sh

if ! [[ -d ${ROOTSYS} ]]
then
	echo ""
	echo "Set the ROOTSYS variable first."
	echo "(You may also want to check the ROOT installation)"
	echo ""
	return
fi

HOME=/home/ndavidson/
export HEPMCLOCATION=/usr/local
export PHOTOSLOCATION=${HOME}/subversion_dir/PHOTOS
export PYTHIALOCATION=${HOME}/pythia8108
export MCTESTERLOCATION=${HOME}/subversion_dir/MC-TESTER/trunk/MC-TESTER
export PYTHIA8DATA=${PYTHIALOCATION}/xmldoc
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${HEPMCLOCATION}/lib
echo ""
echo "Path setup completed."
echo ""
