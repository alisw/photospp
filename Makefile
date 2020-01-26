include make.inc

#Name of libraries to create
LIB_PHOTOSPP_SO        = libPhotospp.so
LIB_PHOTOSPP_A         = libPhotospp.a
LIB_PHOTOSPP_HEPEVT_SO = libPhotosppHEPEVT.so
LIB_PHOTOSPP_HEPEVT_A  = libPhotosppHEPEVT.a
LIB_PHOTOSPP_HEPMC_SO  = libPhotosppHepMC.so
LIB_PHOTOSPP_HEPMC_A   = libPhotosppHepMC.a
LIB_PHOTOSPP_HEPMC3_SO  = libPhotosppHepMC3.so
LIB_PHOTOSPP_HEPMC3_A   = libPhotosppHepMC3.a


#directories containing source code
EVENT_RECORD_INTERFACE_DIR = eventRecordInterfaces
PHOTOS_C_INTERFACE_DIR     = photosCInterfaces
PHOTOS_C_DIR               = photos-C
UTILITIES_DIR              = utilities

LIB_PHOTOSPP_OBJECTS = src/$(PHOTOS_C_INTERFACE_DIR)/*.o \
                       src/$(PHOTOS_C_DIR)/*.o \
                       src/$(UTILITIES_DIR)/*.o

LIB_PHOTOSPP_HEPEVT_OBJECTS =  src/$(EVENT_RECORD_INTERFACE_DIR)/PhotosHEPEVTParticle.o  src/$(EVENT_RECORD_INTERFACE_DIR)/PhotosHEPEVTEvent.o 
LIB_PHOTOSPP_HEPMC_OBJECTS  =  src/$(EVENT_RECORD_INTERFACE_DIR)/PhotosHepMCParticle.o  src/$(EVENT_RECORD_INTERFACE_DIR)/PhotosHepMCEvent.o
LIB_PHOTOSPP_HEPMC3_OBJECTS  = src/$(EVENT_RECORD_INTERFACE_DIR)/PhotosHepMC3Particle.o src/$(EVENT_RECORD_INTERFACE_DIR)/PhotosHepMC3Event.o

##### Link objects to make library ######
all: include_dir lib_dir $(EVENT_RECORD_INTERFACE_DIR) $(FORTRAN_PHOTOS_INTERFACE_DIR) $(PHOTOS_C_INTERFACE_DIR) $(PHOTOS_C_DIR) $(UTILITIES_DIR)
	ar cr lib/$(LIB_PHOTOSPP_A) $(LIB_PHOTOSPP_OBJECTS)
	$(LD) $(LDFLAGS) $(SOFLAGS) -o lib/$(LIB_PHOTOSPP_SO).$(LIB_VER) $(LIB_PHOTOSPP_OBJECTS)
	ar cr lib/$(LIB_PHOTOSPP_HEPEVT_A) $(LIB_PHOTOSPP_HEPEVT_OBJECTS)
	$(LD) $(LDFLAGS) $(SOFLAGS) -o lib/$(LIB_PHOTOSPP_HEPEVT_SO).$(LIB_VER) $(LIB_PHOTOSPP_HEPEVT_OBJECTS)
	ln -sf $(LIB_PHOTOSPP_SO).$(LIB_VER) lib/$(LIB_PHOTOSPP_SO)
	ln -sf $(LIB_PHOTOSPP_HEPEVT_SO).$(LIB_VER) lib/$(LIB_PHOTOSPP_HEPEVT_SO)
ifneq ($(HEPMCLOCATION), )
	ar cr lib/$(LIB_PHOTOSPP_HEPMC_A) $(LIB_PHOTOSPP_HEPMC_OBJECTS)
	$(LD) $(LDFLAGS) $(SOFLAGS) -o lib/$(LIB_PHOTOSPP_HEPMC_SO).$(LIB_VER) $(LIB_PHOTOSPP_HEPMC_OBJECTS)
	ln -sf $(LIB_PHOTOSPP_HEPMC_SO).$(LIB_VER) lib/$(LIB_PHOTOSPP_HEPMC_SO)
endif
ifneq ($(HEPMC3LOCATION), )
	ar cr lib/$(LIB_PHOTOSPP_HEPMC3_A) $(LIB_PHOTOSPP_HEPMC3_OBJECTS)
	$(LD) $(LDFLAGS) $(SOFLAGS) -o lib/$(LIB_PHOTOSPP_HEPMC3_SO).$(LIB_VER) $(LIB_PHOTOSPP_HEPMC3_OBJECTS)
	ln -sf $(LIB_PHOTOSPP_HEPMC3_SO).$(LIB_VER) lib/$(LIB_PHOTOSPP_HEPMC3_SO)
endif

	@echo "##################################################################"	
	@echo " Photos C++ libraries created and moved to lib/ directory         "
	@echo "##################################################################"
	@echo ""
	@echo "##################################################################"	
	@echo " To run examples, cd examples/ directory and there './configure'  "
	@echo " and 'make' again. Examples require Pythia8, ROOT and MC-Tester   "
	@echo "  installed. For details see examples/README.                     "
	@echo "##################################################################"

include_dir:
	mkdir -p include/Photos

lib_dir:
	mkdir -p lib

####### Make object files ########
$(EVENT_RECORD_INTERFACE_DIR):
	+make -C src/$(EVENT_RECORD_INTERFACE_DIR)
	cp src/$(EVENT_RECORD_INTERFACE_DIR)/*.h include/Photos

$(PHOTOS_C_INTERFACE_DIR):
	+make -C src/$(PHOTOS_C_INTERFACE_DIR)
	cp src/$(PHOTOS_C_INTERFACE_DIR)/*.h include/Photos

$(UTILITIES_DIR):
	+make -C src/$(UTILITIES_DIR)
	cp src/$(UTILITIES_DIR)/*.h include/Photos

$(PHOTOS_C_DIR):
	+make -C src/$(PHOTOS_C_DIR)
	cp src/$(PHOTOS_C_DIR)/*.h include/Photos

install:
	mkdir -p $(DESTDIR)$(PREFIX)/include/Photos/
	install  -m 644 include/Photos/* $(DESTDIR)$(PREFIX)/include/Photos/
	mkdir -p $(DESTDIR)$(LIBDIR)/
	install lib/*.$(LIB_VER)  $(DESTDIR)$(LIBDIR)/
	install lib/*.a $(DESTDIR)$(LIBDIR)/
	cp -P lib/*.so  $(DESTDIR)$(LIBDIR)/
clean:
	make clean -C src/$(EVENT_RECORD_INTERFACE_DIR)
	make clean -C src/$(PHOTOS_C_INTERFACE_DIR)
	make clean -C src/$(PHOTOS_C_DIR)
	make clean -C src/$(UTILITIES_DIR)
	rm -f *~

Clean: clean
	rm -f lib/* include/Photos/*
	rm -rf config.log config.status autom4te.cache 
	rm -rf configure.paths.sh configure.paths.csh
	rm -f platform/make.inc make.inc
	rm -f examples/make.inc

make.inc:
	@echo ""
	@echo "Please execute ./configure first!"
	@echo ""
	@false

always:
	@true
