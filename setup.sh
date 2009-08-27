#!/bin/sh

if ! [[ -d ${ROOTSYS} ]]
then
	echo ""
	echo "Set the ROOTSYS variable first."
	echo "(You may also want to check the ROOT installation)"
	echo ""
	return
fi

HERE=`pwd`
export HEPMCLOCATION=../HepMC
export PYTHIALOCATION=../PYTHIA
export TAUOLALOCATION=../TAUOLA
export MCTESTERLOCATION=../MC-TESTER
export PYTHIA8DATA=${PYTHIALOCATION}/xmldoc
export LD_LIBRARY_PATH=${ROOTSYS}/lib:${HEPMCLOCATION}/lib
echo ""
echo "Path setup completed."
echo ""
