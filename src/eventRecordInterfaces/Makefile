include ../make.inc

SOURCES_CC=PhotosHepMCParticle.cxx PhotosHepMCEvent.cxx
OBJECTS_CC=$(SOURCES_CC:.cxx=.o)

%.o: %.cxx
	$(CC) $(CFLAGS) $(DEBUG) -I$(HEPMCLOCATION)/include -I. -I../utilities -I../photosFortranInterfaces -I../photosCInterfaces -c $< -o $@

all:   $(OBJECTS_CC)

clean:
	rm -f *.o *~;
