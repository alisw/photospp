exampledir = $(top_srcdir)/examples

noinst_LTLIBRARIES = libphotos_hepevt_example_interface.la
libphotos_hepevt_example_interface_la_SOURCES = photos_hepevt_example_interface.cxx

### HEPEVT example ###
example_PROGRAMS = photos_hepevt_example.exe

photos_hepevt_example_exe_SOURCES      = photos_hepevt_example.f
photos_hepevt_example_exe_FFLAGS       = -fno-automatic -fno-backslash -ffixed-line-length-132

AM_CPPFLAGS = -I$(prefix)/include
AM_LDFLAGS  = -R $(prefix)/lib

LDADD =                          \
$(FLIBS)                         \
-L$(prefix)/lib/ -lPhotospp -lPhotosppHEPEVT

photos_hepevt_example_exe_LDADD = -lphotos_hepevt_example_interface $(LDADD) -lstdc++


### HepMC3 standalone example
if HAS_HEPMC3

AM_LDFLAGS += -R $(HEPMC3_DIR)/lib
AM_CPPFLAGS += -I$(HEPMC3_DIR)/include
LDADD       += -L$(HEPMC3_DIR)/lib -lHepMC3 -lPhotosppHepMC3

example_PROGRAMS += photos_hepmc3_standalone_example.exe

photos_hepmc3_standalone_example_exe_SOURCES  = photos_hepmc3_standalone_example.cxx

endif # HepMC3

### All other examples require HepMC ###
if HAS_HEPMC

AM_LDFLAGS  += -R $(HEPMC_DIR)/lib
AM_CPPFLAGS += -I$(HEPMC_DIR)/include
LDADD       += -L$(HEPMC_DIR)/lib -lHepMC -lPhotosppHepMC

example_PROGRAMS += photos_standalone_example.exe

photos_standalone_example_exe_SOURCES  = photos_standalone_example.cxx

### Photos gun example (require Pythia8) ###
if HAS_PYTHIA8

AM_LDFLAGS  += -R $(PYTHIA8_DIR)/lib
AM_CPPFLAGS += -I$(PYTHIA8_DIR)/include
LDADD       += -L$(PYTHIA8_DIR)/lib/ -lpythia8

example_PROGRAMS += single_photos_gun_example.exe  photosLCG_pythia_example.exe

single_photos_gun_example_exe_SOURCES  = single_photos_gun_example.cxx
photosLCG_pythia_example_exe_SOURCES   = photosLCG_pythia_example.cxx

### Pythia example (requires Pythia8 and MC-Tester) ###
if HAS_MCTESTER

AM_LDFLAGS  += -R $(MCTESTER_DIR)/lib
AM_CPPFLAGS += -I$(MCTESTER_DIR)/include -I$(ROOTINC)
LDADD       += -L$(MCTESTER_DIR)/lib -lHEPEvent -lHepMCEvent -lMCTester $(ROOTLIBS)

example_PROGRAMS += photos_pythia_example.exe

photos_pythia_example_exe_SOURCES  = photos_pythia_example.cxx

if HAS_TAUOLA

AM_LDFLAGS  += -R $(TAUOLA_DIR)/lib
AM_CPPFLAGS += -I$(TAUOLA_DIR)/include
LDADD       += -L$(TAUOLA_DIR)/lib -lTauolaCxxInterface -lTauolaFortran -lTauolaHepMC

example_PROGRAMS += tauola_photos_pythia_example.exe

tauola_photos_pythia_example_exe_SOURCES  = tauola_photos_pythia_example.cxx

endif # TAUOLA

endif # MC-TESTER

endif # PYTHIA8

endif # HEPMC
