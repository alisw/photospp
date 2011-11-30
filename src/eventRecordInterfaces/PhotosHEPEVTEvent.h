#ifndef _PhotosHEPEVTEvent_h_included_
#define _PhotosHEPEVTEvent_h_included_

/**
 * @class PhotosHEPEVTParticle
 *
 * @brief Single particle of HEPEVT event record
 *
 * This class implements the virtual methods of
 * PhotosEvent. In this way it provides an
 * interface between the generic PhotosEvent class
 * and information stored in HEPEVT event record.
 *
 * @author Tomasz Przedzinski
 * @date 24 November 2011
 */

#include <iostream>
#include "PhotosEvent.h"
#include "PhotosParticle.h"
#include "PhotosHEPEVTParticle.h"

// Uncomment this line to use interface to common block HEPEVT
//#define USE_HEPEVT_INTERFACE

#ifdef USE_HEPEVT_INTERFACE

// Change this value to match HEPEVT size
const int NMXHEP = 10000;

extern "C" struct {
  int   nevhep;            // serial number
  int   nhep;              // number of particles
  int   isthep[NMXHEP];    // status code
  int   idhep [NMXHEP];    // particle PDG ID
  int   jmohep[NMXHEP][2]; // parent particles
  int   jdahep[NMXHEP][2]; // childreen particles
  float phep  [NMXHEP][5]; // four-momentum, mass [GeV]
  float vhep  [NMXHEP][4]; // vertex [mm]
} hepevt_;

#endif

class PhotosHEPEVTParticle;

class PhotosHEPEVTEvent : public PhotosEvent {

 public:

  /** Default destructor */
  ~PhotosHEPEVTEvent();

  /** Default constructor */
  PhotosHEPEVTEvent();

  /** Add particle at the end of event record */
  void addParticle(PhotosHEPEVTParticle *p);

  /** Get particle at index 'i' */
  PhotosHEPEVTParticle *getParticle(int i);

  /** Remove particle.

      Simplest implementation. This function does not change
      indexes of any particles and does not change position
      of the particles but instead creates empty space
      at index 'i' */
  void removeParticle(int i);

  /** Get higher-most index of the particles in event */
  int getParticleCount();

	/** Get an unfiltered list of particles from the event record */
	virtual vector<PhotosParticle*> getParticleList();

  /** Print out list of particles in the event */
  void print();
  
  /** Remove all particles from the event */
  void clear();

#ifdef USE_HEPEVT_INTERFACE
  /** Fill PhotosHEPEVTEvent from HEPEVT common block */
  static void fill_event_from_HEPEVT(PhotosHEPEVTEvent *evt);
  
  /** Write to HEPEVT common block content of PhotosHEPEVTEvent */
  static void write_event_to_HEPEVT(PhotosHEPEVTEvent *evt);
#endif

 private:

  /** List of all particles */
  std::vector<PhotosHEPEVTParticle*> particle_list;
};

#endif

