#include ../../make.inc
COMMAND = phopar.exe
COMMAND_OBJECTS = phopar.o photos.o pseudo_jetset.o pary.o pairs.o
############
all:	$(COMMAND)
FFLAGS  = -fno-automatic -g -ffixed-line-length-132 -fno-backslash
LDFLAGS = -fno-automatic -g -lstdc++
CC      = g++
CFLAGS  = -fPIC -Wall -Wno-write-strings
F77 = gfortran
LD  = $(F77)
###########
.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

.f.o:
	$(F77) $(FFLAGS) -c $<
# Link procedure
$(COMMAND): $(COMMAND_OBJECTS)   
	$(LD) -o $@ $(LDFLAGS) $(COMMAND_OBJECTS)  
run:
	./phopar.exe
clean:                                                     
	rm -f *.o *.a  *~ *.exe
Clean:	clean
	cd ../lib; rm -f *.o *.a  *~ *.exe
	cd lib;    rm -f *.o *.a  *~ *.exe
