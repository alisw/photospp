#!/bin/sh

if test ! -d ../src; then
  echo "WARNING! Script must be running from PHOTOS/platform directory!"
  return
fi

echo "This option will overwrite all configuration scripts/makefiles"
echo "and modify the configuration procedure to match LCG setup."
echo ""
echo "You will need autotools version 2.59 or higher."
echo ""
echo "Proceed? (Yes/No)"
read ANSWER

ANSWER=`echo $ANSWER | tr "[:upper:]" "[:lower:]"`

if test "$ANSWER" = "yes" || test "$ANSWER" = "y"; then
  echo "Renaming .c files"
  mv ../src/photosFortranInterfaces/PH_HEPEVT_Interface.c ../src/photosFortranInterfaces/PH_HEPEVT_Interface.cxx

  mv ../examples/photos_standalone_example.c ../examples/photos_standalone_example.cxx
  mv ../examples/single_photos_gun_example.c ../examples/single_photos_gun_example.cxx
  mv ../examples/photos_pythia_example.c ../examples/photos_pythia_example.cxx
  mv ../examples/tauola_photos_pythia_example.c ../examples/tauola_photos_pythia_example.cxx

  echo "Removing previous installation scripts"
  rm -rf ../config* ../make* ../Make*
  rm -rf ../src/make.inc ../src/*/Makefile ../photos-fortran/make*
  rm -rf ../examples/config* ../examples/make* ../examples/Make*

  echo "Copying and configuring new scripts"
  cp -rf LCGCONFIG/* ../.
  cd ..
  autoreconf --install --force
	echo "Done."
else
	echo "Aborted."
fi
