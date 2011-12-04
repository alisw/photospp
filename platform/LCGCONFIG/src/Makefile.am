includedir=$(prefix)/include/Photos

lib_LTLIBRARIES = libPhotosCxxInterface.la

libPhotosCxxInterface_la_SOURCES =              \
eventRecordInterfaces/PhotosHepMCEvent.cxx      \
eventRecordInterfaces/PhotosHepMCParticle.cxx   \
photosCInterfaces/Photos.cxx                    \
photosCInterfaces/PhotosParticle.cxx            \
photosCInterfaces/PhotosBranch.cxx              \
photosCInterfaces/PhotosEvent.cxx               \
photosFortranInterfaces/PH_HEPEVT_Interface.cxx \
utilities/Log.cxx                               \
utilities/PhotosRandom.cxx

include_HEADERS =                             \
eventRecordInterfaces/PhotosHepMCParticle.h   \
eventRecordInterfaces/PhotosHepMCEvent.h      \
photosCInterfaces/Photos.h                    \
photosCInterfaces/PhotosParticle.h            \
photosCInterfaces/PhotosBranch.h              \
photosCInterfaces/PhotosEvent.h               \
photosFortranInterfaces/PH_HEPEVT_Interface.h \
photosFortranInterfaces/f_Init.h              \
utilities/Log.h                               \
utilities/PhotosRandom.h

INCLUDES = -I$(top_srcdir)/src/photosFortranInterfaces \
           -I$(top_srcdir)/src/photosCInterfaces \
           -I$(top_srcdir)/src/utilities \
           -I$(HEPMC_DIR)/include

libPhotosCxxInterface_la_FFLAGS = -fno-automatic -fno-backslash -ffixed-line-length-132
