include $(PHOTOSLOCATION)/platform/make.inc
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

GBLIB_LIB  = $(PHOTOSLOCATION)/photos/glib.a
PHOTOS_LIB  = $(PHOTOSLOCATION)/photos/photos.o
PHOTOS_CXX_INT_LIB= $(PHOTOSLOCATION)/demo_cxx/lib/libPhotosCxxInterface.so

#TAUOLA LIBRARIES
TAUOLALOCATION = ../../../TAUOLA/trunk/TAUOLA
TAUOLA_LIB  = $(TAUOLALOCATION)/tauola/glib.a
TAUOLA_CXX_INT_LIB=$(TAUOLALOCATION)/demo_cxx/lib/libTauolaCxxInterface.a

############
all:	$(COMMAND)

.f.o:
	$(F77) $(FFLAGS) -c $<

.cxx.o:
	$(CC) $(CFLAGS) -c $< 

$(COMMAND_OBJECTS): $(COMMAND_C) pythia_loc hepmc_loc
	$(CC) $(CFLAGS) -I$(PYTHIALOCATION)/include $(ROOTINC) \
	-I$(MCTESTERLOCATION)/include -I$(HEPMCLOCATION)/include \
	-I$(TAUOLALOCATION)/demo_cxx/include -I../include/ -c $(COMMAND_C)

# Link procedure
$(COMMAND): photos_loc pythia_loc hepmc_loc $(COMMAND_OBJECTS)
	$(LD) -o $@ $(COMMAND_OBJECTS) $(PHOTOS_LIB) $(PHOTOS_CXX_INT_LIB) \
	$(ROOTINC) -L$(HEPMCLOCATION)/lib -lHepMC -lHepMCfio \
	-L$(PYTHIALOCATION)/lib/archive -lpythia8 -llhapdfdummy \
	-lhepmcinterface $(ROOTLIB) $(MC_TESTER_LIBS) \
	$(TAUOLA_CXX_INT_LIB) $(TAUOLA_LIB) \
	$(TAUOLA_CXX_INT_LIB) $(PHOTOS_CXX_INT_LIB)
	@echo "#########################################################"	
	@echo "Please run the example by typing ./"$(COMMAND)
	@echo "#########################################################"
# -lstdc++ 
#$(PHOTOS_CXX_INT_LIB) 
#@echo "#####################################" \


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

photos_loc: #check tauola root is set
ifeq ($(PHOTOSLOCATION), )
	@echo "##############################################"
	@echo "Please set PHOTOSLOCATION to the "
	@echo "directory where tauola.o can be found "
	@echo "##############################################"
	@false
else
	@echo "Using PHOTOSLOCATION ="$(PHOTOSLOCATION)
endif

clean:                                                     
	rm -f *.o; rm -f *.a; rm -f *~; rm -f *.exe; rm -f \#*
