include platform/make.inc

#Name of libraries to create
LIB_NAME_SO = lib/libPhotosCxxInterface.so
LIB_NAME_A = lib/libPhotosCxxInterface.a

all:
	make -C photos-fortran
	make -C src
	cp src/*.h include
	ar cr $(LIB_NAME_A) src/*.o
#	$(LD) $(LDFLAGS) $(SOFLAGS) -o $(LIB_NAME_SO) $(LIB_NAME_A)
	$(LD) $(SOFLAGS) -o $(LIB_NAME_SO) src/*.o
	@echo "##################################################################"
	@echo "Photos C++ Interface library created and moved to lib/ directory"
	@echo "##################################################################"

clean:
	@rm -f *.o; rm -f *.a; rm -f *~
	make clean -C photos-fortran
	make clean -C src
	rm -f include/*.h
	rm -f lib/*.h

