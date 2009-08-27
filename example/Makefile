include ../platform/make.inc
ROOTLIB  := $(shell $(ROOTSYS)/bin/root-config --glibs)
ROOTINC := $(shell $(ROOTSYS)/bin/root-config --cflags)
MC_TESTER_LIBS = $(MCTESTERLOCATION)/lib/libHEPEvent.so \
	$(MCTESTERLOCATION)/lib/libHepMCEvent.so \
	$(MCTESTERLOCATION)/lib/libMCTester.so

#MAIN=photos_pythia_example
MAIN=tauola_photos_pythia_example

COMMAND = $(MAIN).exe
COMMAND_OBJECTS = $(MAIN).o
COMMAND_C = $(MAIN).c	

HERE=`pwd`
GBLIB_LIB  = $(HERE)/../photos-fortran/glib.a
PHOTOS_LIB  = $(HERE)/../photos-fortran/photos.o
PHOTOS_CXX_INT_LIB= $(HERE)/../lib/libPhotosCxxInterface.so

#TAUOLA LIBRARIES
TAUOLA_LIB  = $(TAUOLALOCATION)/tauola-fortran/tauola/glib.a
TAUOLA_CXX_INT_LIB=$(TAUOLALOCATION)/lib/libTauolaCxxInterface.a

############
all:	$(COMMAND)

.f.o:
	$(F77) $(FFLAGS) -c $<

.cxx.o:
	$(CC) $(CFLAGS) -c $< 

$(COMMAND_OBJECTS): $(COMMAND_C) pythia_loc hepmc_loc
	$(CC) -g $(CFLAGS) -I$(PYTHIALOCATION)/include $(ROOTINC) \
	-I$(MCTESTERLOCATION)/include -I$(HEPMCLOCATION)/include \
	-I$(TAUOLALOCATION)/include -I../include/ -c $(COMMAND_C)

# Link procedure
$(COMMAND): tauola_loc pythia_loc hepmc_loc $(COMMAND_OBJECTS)
	$(LD) -o $@ $(COMMAND_OBJECTS) $(PHOTOS_LIB) $(PHOTOS_CXX_INT_LIB) \
	$(ROOTINC) -L$(HEPMCLOCATION)/lib -lHepMC -lHepMCfio \
	-L$(PYTHIALOCATION)/lib/archive -lpythia8 -llhapdfdummy \
	-lhepmcinterface $(ROOTLIB) $(MC_TESTER_LIBS) \
	$(TAUOLA_CXX_INT_LIB) $(TAUOLA_LIB) \
	$(TAUOLA_CXX_INT_LIB) $(PHOTOS_CXX_INT_LIB)
	@echo "#########################################################"	
	@echo "Please run the example by typing ./"$(COMMAND)
	@echo "#########################################################"

pythia_loc: #check pythia root is set
ifeq ($(PYTHIALOCATION), )
	@echo "##############################################"
	@echo "Please set PYTHIALOCATION to the "
	@echo "directory where your Pythia 8 include/ and "
	@echo "lib/ directories can be found "
	@echo "##############################################"
	@false
else
	@echo "Using PYTHIALOCATION="$(PYTHIALOCATION)
endif	

hepmc_loc: #check HepMC root is set
ifeq ($(HEPMCLOCATION), )
	@echo "##############################################"
	@echo "Please set HEPMCLOCATION to the "
	@echo "directory where your HepMC include/ and "
	@echo "lib/ directories can be found "
	@echo "##############################################"
	@false
else
	@echo "Using HEPMCLOCATION="$(HEPMCLOCATION)
endif		

tauola_loc: #check tauola root is set
ifeq ($(TAUOLALOCATION), )
	@echo "##############################################"
	@echo "Please set TAUOLALOCATION to the "
	@echo "directory where your TAUOLA include/ and "
	@echo "lib/ directories can be found "
	@echo "##############################################"
	@false
else
	@echo "Using TAUOLALOCATION ="$(TAUOLALOCATION)
endif

clean:                                                     
	rm -f *.o; rm -f *.a; rm -f *~; rm -f *.exe; rm -f \#*

