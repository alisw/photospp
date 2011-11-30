exampledir = $(top_srcdir)/examples

### Standalone example ###
example_PROGRAMS = photos_standalone_example.exe

photos_standalone_example_exe_SOURCES  = photos_standalone_example.cxx

INCLUDES = -I$(prefix)/include   \
           -I$(HEPMC_DIR)/include

AM_LDFLAGS = -R $(prefix)/lib -R $(HEPMC_DIR)/lib

LDADD =                          \
$(FLIBS)                         \
../lib/libPhotosCxxInterface.so  \
-L$(HEPMC_DIR)/lib -lHepMC

### Photos gun example (require Pythia8) ###
if HAS_PYTHIA8

AM_LDFLAGS  += -R $(PYTHIA8_DIR)/lib/archive

INCLUDES += -I$(PYTHIA8_DIR)/include
LDADD    += -L$(PYTHIA8_DIR)/lib/archive -lpythia8 -llhapdfdummy -lhepmcinterface

example_PROGRAMS += single_photos_gun_example.exe  photosLCG_pythia_example.exe

single_photos_gun_example_exe_SOURCES  = single_photos_gun_example.cxx
photosLCG_pythia_example_exe_SOURCES   = photosLCG_pythia_example.cxx

### Pythia example (requires Pythia8 and MC-Tester) ###
if HAS_MCTESTER

AM_LDFLAGS  += -R $(MCTESTER_DIR)/lib

INCLUDES += -I$(MCTESTER_DIR)/include -I$(ROOTINC)
LDADD    += -L$(MCTESTER_DIR)/lib -lHEPEvent -lHepMCEvent -lMCTester $(ROOTLIBS)

example_PROGRAMS += photos_pythia_example.exe

photos_pythia_example_exe_SOURCES  = photos_pythia_example.cxx

if HAS_TAUOLA

AM_LDFLAGS  += -R $(TAUOLA_DIR)/lib

INCLUDES += -I$(TAUOLA_DIR)/include
LDADD    += -L$(TAUOLA_DIR)/lib -lTauolaFortran -lTauolaCxxInterface

example_PROGRAMS += tauola_photos_pythia_example.exe

tauola_photos_pythia_example_exe_SOURCES  = tauola_photos_pythia_example.cxx

endif # TAUOLA

endif # MC-TESTER

endif # PYTHIA8
