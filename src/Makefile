include ../platform/make.inc

SOURCES_F=Photos_make.f
OBJECTS_F=$(SOURCES_F:.f=.o)

SOURCES_C= Photos.cxx PhotosHepMCEvent.cxx PhotosHepMCParticle.cxx \
	PhotosParticle.cxx PH_HEPEVT_Interface.cxx Log.cxx
OBJECTS_C=$(SOURCES_C:.cxx=.o)

############
all:
	$(F77) $(FFLAGS) -c $(SOURCES_F)
	$(CC) $(CFLAGS) -I $(HEPMCLOCATION)/include -c $(SOURCES_C)

clean:
	rm -f *.o; rm -f *~; rm -f \#*	

