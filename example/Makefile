include ../platform/make.inc
include make.inc

PHOTOS        :=  $(shell cd .. && pwd)
ROOTLIB       :=  $(shell $(ROOTSYS)/bin/root-config --glibs)
ROOTINC       :=  $(shell $(ROOTSYS)/bin/root-config --cflags)
TAUOLA_LIBS    = -L$(TAUOLALOCATION)/lib -lTauolaCxxInterface -lTauolaFortran
PYTHIA_LIBS    = -L$(PYTHIALOCATION)/lib/archive -lpythia8 -llhapdfdummy -lhepmcinterface
MC_TESTER_LIBS = -L$(MCTESTERLOCATION)/lib -lHEPEvent -lHepMCEvent -lMCTester
HEPMC_LIBS     = -L$(HEPMCLOCATION)/lib -lHepMC
PHOTOS_LIBS    = -L../lib -lPhotosCxxInterface

############
all:
	@make tauola_photos_pythia_example.exe
	@echo "#########################################################"
	@echo " Please run the example by typing:"
	@echo " ./tauola_photos_pythia_example.exe"
	@echo "#########################################################"


# Universal example (all paths presents)
%.o: %.c
	$(CC) $(CFLAGS) -I../include -I$(PYTHIALOCATION)/include \
	-I$(MCTESTERLOCATION)/include -I$(TAUOLALOCATION)/include \
	-I$(HEPMCLOCATION)/include $(ROOTINC) -c $<

# Universal example (all paths presents) - linking procedure
%.exe: %.o
	$(CC) -o $@ $(LDFLAGS) $< $(ROOTLIB) $(TAUOLA_LIBS) \
	$(PYTHIA_LIBS) $(MC_TESTER_LIBS) $(HEPMC_LIBS) $(PHOTOS_LIBS) \
	-Wl,-rpath,$(PREFIX)/lib \
	-Wl,-rpath,$(PHOTOS)/lib \
	-Wl,-rpath,$(ROOTSYS)/lib \
	-Wl,-rpath,$(HEPMCLOCATION)/lib \
	-Wl,-rpath,$(TAUOLALOCATION)/lib \
	-Wl,-rpath,$(PYTHIALOCATION)/lib \
	-Wl,-rpath,$(MCTESTERLOCATION)/lib

clean:
	rm -f *.o *.a *~ *.exe;

Clean: clean
	rm -rf config.log config.status autom4te.cache configure.paths.sh
	rm -f make.inc

make.inc:
	@echo ""
	@echo "Please execute ./configure first!"
	@echo ""
	@false

../platform/make.inc:
	@echo ""
	@echo "Please compile the interface first!"
	@echo ""
	@false