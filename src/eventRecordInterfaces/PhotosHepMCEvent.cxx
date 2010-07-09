#include "PhotosHepMCParticle.h"
#include "PhotosHepMCEvent.h"
#include "Log.h"
using namespace std;

PhotosHepMCEvent::PhotosHepMCEvent(HepMC::GenEvent * event)
{
	m_event=event;
}

HepMC::GenEvent * PhotosHepMCEvent::getEvent()
{
	return m_event;
}

void PhotosHepMCEvent::print()
{
	if(!m_event) return;
	m_event->print();
}

vector<PhotosParticle*> PhotosHepMCEvent::getParticleList()
{
	vector<PhotosParticle*> list;
	
	HepMC::GenEvent::particle_const_iterator part_itr = m_event->particles_begin();
	for( ; part_itr!=m_event->particles_end(); part_itr++)
	{
		PhotosParticle *particle = new PhotosHepMCParticle(*part_itr);
		list.push_back(particle);
	}
	return list;
}
