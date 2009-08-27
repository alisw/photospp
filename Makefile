#############

all:
	make -C photos
	make -C demo_cxx
	make -C demo_cxx/example

clean:
	@rm -f *.o; rm -f *.a; rm -f *~
	make clean -C photos
	make clean -C demo_cxx
	make clean -C demo_cxx/example

include $(PHOTOSLOCATION)/platform/make.inc

#Name of libraries to create
LIB_NAME_SO = lib/libPhotosCxxInterface.so
LIB_NAME_A = lib/libPhotosCxxInterface.a

##### Link objects to make library ######
all:
	make -C src
	cp src/*.h include
	ar cr $(LIB_NAME_A) src/*.o
#	$(LD) $(LDFLAGS) $(SOFLAGS) -o $(LIB_NAME_SO) $(LIB_NAME_A)
	$(LD) $(SOFLAGS) -o $(LIB_NAME_SO) src/*.o

	@echo "##################################################################"	
	@echo "Photos C++ Interface library created and moved to lib/ directory"
	@echo "##################################################################"


clean:                                                     
	rm -f *~
	make clean -C src

clobber: clean
	rm -f $(LIB_NAME_SO) $(LIB_NAME_A)
	rm include/*