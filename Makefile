include platform/make.inc

#Name of libraries to create
LIB_NAME_SO = lib/libPhotosCxxInterface.so
LIB_NAME_A = lib/libPhotosCxxInterface.a

all:
	make -C photos-fortran
	make -C src
	cp src/*.h include
	ar cr $(LIB_NAME_A) src/*.o photos-fortran/photos.o
	$(LD) $(LDFLAGS) $(SOFLAGS) -o $(LIB_NAME_SO) src/*.o \
	 photos-fortran/photos.o
	@echo "##################################################################"
	@echo "Photos C++ Interface library created and moved to lib/ directory"
	@echo "##################################################################"

install:
	mkdir -p $(PREFIX)/include
	cp include/* $(PREFIX)/include/.
	mkdir -p $(PREFIX)/lib
	cp lib/* $(PREFIX)/lib/.

clean:
	@rm -f *.o; rm -f *.a; rm -f *~
	make clean -C photos-fortran
	make clean -C src
	rm -f include/*.h
	rm -f lib/*.h

Clean: clean
	rm -f lib/* include/*
	rm -rf config.log config.status autom4te.cache configure.paths.sh
	rm -f platform/make.inc
	rm -f examples/make.inc


platform/make.inc:
	@echo ""
	@echo "Please execute ./configure first!"
	@echo ""
	@false

