include ../make.inc

SOURCES_F=Photos_make.f forW-ME.f forZ-ME.f
OBJECTS_F=$(SOURCES_F:.f=.o)

SOURCES_C=PH_HEPEVT_Interface.c
OBJECTS_C=$(SOURCES_C:.c=.o)

%.o: %.c
	$(CC) $(CFLAGS) $(DEBUG) -I. -I../photosCInterfaces -I../photosFortranInterfaces -I../utilities -c $< -o $@

%.o: %.f
	$(F77) $(FFLAGS) $(DEBUG) -I. -c $< -o $@

all: $(OBJECTS_C) $(OBJECTS_F)

clean:
	rm -f *.o *~;
