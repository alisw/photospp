#ifndef _PhotosHepMCEvent_h_included_
#define _PhotosHepMCEvent_h_included_

/**
 * @class PhotosHepMCEvent
 *
 * @brief Interface to HepMC::GenEvent objects
 *
 * This class implements the virtual methods of
 * PhotosEvent. In this way it provides an
 * interface between the generic PhotosEvent class
 * and a HepMC::GenEvent object.
 *
 * @author Nadia Davidson
 * @date 17 June 2008
 */

#include <iostream>
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "PhotosEvent.h"
#include "PhotosParticle.h"
#include "PhotosHepMCParticle.h"

class PhotosHepMCEvent : public PhotosEvent{

 public:

  /** Constructor which keeps a pointer to the HepMC::GenEvent*/
  PhotosHepMCEvent(HepMC::GenEvent * event);

  /** Destructor */
  ~PhotosHepMCEvent();

  /** Returns the HepMC::GenEvent */
  HepMC::GenEvent * getEvent();

  std::vector<PhotosParticle*> getBranchPoints();

 private:
  /** The event */
  HepMC::GenEvent * m_event;

  /** filter for branching points **/
  bool passBranchPointFilter(HepMC::GenParticle * particle);

  /** Check if the branching point should be skipped by PHOTOS. */
  bool passSuppressionFilter(PhotosHepMCParticle *particle);
  /** Check if the particles' mother is on the Photos::supParticles list.
      If it is, it will also be skipped. */
  bool passSuppressConsecutive(PhotosHepMCParticle *particle);

  /** branch points which should be given to PHOTOS */
  std::vector<PhotosParticle*> m_branch_points;

};

#endif  

