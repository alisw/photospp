include ../make.inc
include make.inc

PHOTOS        :=  $(shell cd .. && pwd)
ROOTLIB       :=  $(shell root-config --glibs)
ROOTINC       :=  $(shell root-config --cflags)
ROOTLIBPATH   :=  $(shell root-config --libdir)
TAUOLA_LIBS    = -L${TAUOLALOCATION}/lib -lTauolaCxxInterface -lTauolaFortran
PYTHIA_LIBS    = -L$(PYTHIALOCATION)/lib/archive -lpythia8 -llhapdfdummy -lhepmcinterface
MC_TESTER_LIBS = -L$(MCTESTERLOCATION)/lib -lHEPEvent -lHepMCEvent -lMCTester
HEPMC_LIBS     = -L$(HEPMCLOCATION)/lib -lHepMC
PHOTOS_LIBS    = -L../lib -lPhotosCxxInterface

all:
ifneq ($(PYTHIALOCATION), )
ifneq ($(MCTESTERLOCATION), )
	@make photos_pythia_example.exe
ifneq ($(TAUOLALOCATION), )
	@make tauola_photos_pythia_example.exe
endif

	@echo "##################################################################"
	@echo " Examples compiled successfully. Available examples:"
	@if test -e photos_pythia_example.exe; then \
	echo " ./photos_pythia_example.exe"; fi
	@if test -e tauola_photos_pythia_example.exe; then \
	echo " ./tauola_photos_pythia_example.exe"; fi
	@echo "##################################################################"

else
	@echo "MC-Tester not present."
	@echo "Both MC-Tester and Pythia are required to compile the examples."
	@exit 0;
endif
else
	@echo "Pythia not present."
	@echo "Both MC-Tester and Pythia are required to compile the examples."
	@exit 0;
endif




# Universal example (all paths presents)
%.o: %.c
	$(CC) $(CFLAGS) -I../include -I$(TAUOLALOCATION)/include -I$(PYTHIALOCATION)/include \
	-I$(MCTESTERLOCATION)/include -I$(HEPMCLOCATION)/include $(ROOTINC) -c $<

# Universal example (all paths presents) - linking procedure
%.exe: %.o
	$(CC) -o $@ $(LDFLAGS) $< $(ROOTLIB) \
	$(TAUOLA_LIBS) $(PYTHIA_LIBS) $(MC_TESTER_LIBS) $(HEPMC_LIBS) $(PHOTOS_LIBS) \
	-Wl,-rpath,$(PREFIX)/lib \
	-Wl,-rpath,$(PHOTOS)/lib \
	-Wl,-rpath,$(ROOTLIBPATH) \
	-Wl,-rpath,$(TAUOLALOCATION)/lib \
	-Wl,-rpath,$(HEPMCLOCATION)/lib \
	-Wl,-rpath,$(PYTHIALOCATION)/lib \
	-Wl,-rpath,$(MCTESTERLOCATION)/lib

# Version for pythia example - HepMC, Pythia8 and MC-Tester required
photos_pythia_example.exe: photos_pythia_example.c
	$(CC) $(CFLAGS) -I../include -I$(PYTHIALOCATION)/include \
	-I$(MCTESTERLOCATION)/include -I$(HEPMCLOCATION)/include $(ROOTINC) -c $<
	$(CC) -o $@ $(LDFLAGS) $(<:.c=.o) $(ROOTLIB) $(HEPMC_LIBS) $(PYTHIA_LIBS) \
	$(PYTHIA_LIBS) $(MC_TESTER_LIBS) $(HEPMC_LIBS) $(PHOTOS_LIBS) \
	-Wl,-rpath,$(PREFIX)/lib \
	-Wl,-rpath,$(PHOTOS)/lib \
	-Wl,-rpath,$(ROOTLIBPATH) \
	-Wl,-rpath,$(HEPMCLOCATION)/lib \
	-Wl,-rpath,$(PYTHIALOCATION)/lib \
	-Wl,-rpath,$(MCTESTERLOCATION)/lib

clean:
	rm -f *.o *.a *~ *.exe

Clean: clean
	rm -rf config.log config.status autom4te.cache configure.paths.sh
	rm -f make.inc

make.inc:
	@echo ""
	@echo "Please execute ./configure first!"
	@echo ""
	@false

../make.inc:
	@echo ""
	@echo "Please compile the interface first!"
	@echo ""
	@false
