F77      = gfortran
F77FLAGS = -fPIC -O3 -fno-automatic -fno-backslash -ffixed-line-length-132
FFLAGS   = $(F77FLAGS)

LD      = gfortran
LDFLAGS = -O3
SOFLAGS  = -shared

RANLIB   = ranlib
AR       = ar

CPP     = /lib/cpp 
CPPFLAGS= -C -P -traditional-cpp -O3
DIFF    = diff -ib

CC=g++
CFLAGS += -fPIC -O3 -Wall

.SUFFIXES: .c .cxx .o .f