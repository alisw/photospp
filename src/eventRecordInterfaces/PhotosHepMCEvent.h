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

#include "HepMC/GenEvent.h"
#include "PhotosEvent.h"
#include "PhotosParticle.h"

class PhotosHepMCEvent : public PhotosEvent
{
public:
	/** Constructor which keeps a pointer to the HepMC::GenEvent*/
	PhotosHepMCEvent(HepMC::GenEvent * event);

	/** Returns the HepMC::GenEvent */
	HepMC::GenEvent * getEvent();

	std::vector<PhotosParticle*> getParticleList();

	void print();

private:
	/** The event */
	HepMC::GenEvent * m_event;
};

#endif  
