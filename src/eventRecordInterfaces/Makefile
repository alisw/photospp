include ../make.inc

OBJECTS_CC = PhotosHEPEVTEvent.o PhotosHEPEVTParticle.o

ifneq ($(HEPMCLOCATION), )
  OBJECTS_CC += PhotosHepMCParticle.o PhotosHepMCEvent.o
endif

ifneq ($(HEPMC3LOCATION), )
  OBJECTS_CC += PhotosHepMC3Particle.o PhotosHepMC3Event.o
endif

all:   $(OBJECTS_CC)

PhotosHEPEVTEvent.o: PhotosHEPEVTEvent.cxx
	$(CC) $(CFLAGS) $(DEBUG)  -I. -I../utilities -I../photos-C -I../photosCInterfaces -c $< -o $@


PhotosHEPEVTParticle.o: PhotosHEPEVTParticle.cxx
	$(CC) $(CFLAGS) $(DEBUG)  -I. -I../utilities -I../photos-C -I../photosCInterfaces -c $< -o $@



PhotosHepMCParticle.o: PhotosHepMCParticle.cxx
	$(CC) $(CFLAGS) $(DEBUG) -I$(HEPMCLOCATION)/include -I. -I../utilities -I../photos-C -I../photosCInterfaces -c $< -o $@

PhotosHepMCEvent.o: PhotosHepMCEvent.cxx
	$(CC) $(CFLAGS) $(DEBUG) -I$(HEPMCLOCATION)/include -I. -I../utilities -I../photos-C -I../photosCInterfaces -c $< -o $@


PhotosHepMC3Particle.o: PhotosHepMC3Particle.cxx
	$(CC) $(CFLAGS) $(DEBUG) -I$(HEPMC3LOCATION)/include -I. -I../utilities -I../photos-C -I../photosCInterfaces -c $< -o $@

PhotosHepMC3Event.o: PhotosHepMC3Event.cxx
	$(CC) $(CFLAGS) $(DEBUG) -I$(HEPMC3LOCATION)/include -I. -I../utilities -I../photos-C -I../photosCInterfaces -c $< -o $@

clean:
	rm -f *.o *~;
