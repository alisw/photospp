include make.inc

LIB_VER = 1.0.0

#Name of libraries to create
LIB_PHOTOS_CXX_INT_SO = libPhotosCxxInterface.so
LIB_PHOTOS_CXX_INT_A  = libPhotosCxxInterface.a

PHOTOS_CXX_INT_OBJECTS = src/$(FORTRAN_PHOTOS_INTERFACE_DIR)/*.o \
                         src/$(EVENT_RECORD_INTERFACE_DIR)/*.o \
                         src/$(C_PHOTOS_INTERFACE_DIR)/*.o \
                         src/$(UTILITIES_DIR)/*.o \
                         photos-fortran/*.o 

#directories containing source code
EVENT_RECORD_INTERFACE_DIR   = eventRecordInterfaces
FORTRAN_PHOTOS_INTERFACE_DIR = photosFortranInterfaces
C_PHOTOS_INTERFACE_DIR       = photosCInterfaces
UTILITIES_DIR                = utilities

all: photos_fortran src

##### Link objects to make library ######
src: include_dir $(EVENT_RECORD_INTERFACE_DIR) $(FORTRAN_PHOTOS_INTERFACE_DIR) $(C_PHOTOS_INTERFACE_DIR) $(UTILITIES_DIR)
	ar cr lib/$(LIB_PHOTOS_CXX_INT_A) $(PHOTOS_CXX_INT_OBJECTS)
	$(LD) $(LDFLAGS) $(SOFLAGS) -o lib/$(LIB_PHOTOS_CXX_INT_SO).$(LIB_VER) $(PHOTOS_CXX_INT_OBJECTS)
	ln -sf $(LIB_PHOTOS_CXX_INT_SO).$(LIB_VER) lib/$(LIB_PHOTOS_CXX_INT_SO)
	@echo "##################################################################"	
	@echo " Photos C++ Interface library created and moved to lib/ directory "
	@echo "##################################################################"
	@echo ""
	@echo "##################################################################"	
	@echo " To run examples, cd examples/ directory and there './configure'  "
	@echo " and 'make' again. Examples require Pythia8, ROOT and MC-Tester   "
	@echo "  installed. For details see examples/README.                     "
	@echo "##################################################################"

include_dir:
	mkdir -p include/Photos

####### Make object files ########
$(FORTRAN_PHOTOS_INTERFACE_DIR):
	make -C src/$(FORTRAN_PHOTOS_INTERFACE_DIR)
	cp src/$(FORTRAN_PHOTOS_INTERFACE_DIR)/*.h include/Photos

$(EVENT_RECORD_INTERFACE_DIR):
	make -C src/$(EVENT_RECORD_INTERFACE_DIR)
	cp src/$(EVENT_RECORD_INTERFACE_DIR)/*.h include/Photos

$(C_PHOTOS_INTERFACE_DIR):
	make -C src/$(C_PHOTOS_INTERFACE_DIR)
	cp src/$(C_PHOTOS_INTERFACE_DIR)/*.h include/Photos

$(UTILITIES_DIR):
	make -C src/$(UTILITIES_DIR)
	cp src/$(UTILITIES_DIR)/*.h include/Photos

photos_fortran:
	make -C photos-fortran

install:
	mkdir -p $(PREFIX)/include/Photos
	cp include/Photos/* $(PREFIX)/include/Photos/.
	mkdir -p $(PREFIX)/lib
	cp lib/* $(PREFIX)/lib/.

clean_src:
	make clean -C src/$(EVENT_RECORD_INTERFACE_DIR)
	make clean -C src/$(FORTRAN_PHOTOS_INTERFACE_DIR)
	make clean -C src/$(C_PHOTOS_INTERFACE_DIR)
	make clean -C src/$(UTILITIES_DIR)

clean: clean_src
	make clean -C photos-fortran
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
	@echo "(Consider using 'source platform/afs.paths.sh' [or .csh] first)"
	@echo ""
	@false

always:
	@true
