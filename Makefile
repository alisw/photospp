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

