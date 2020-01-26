#include <vector>
#include "PhotosHepMC3Particle.h"
#include "PhotosHepMC3Event.h"
#include "Log.h"

#include "HepMC3/Print.h"

namespace Photospp
{
PhotosHepMC3Event::PhotosHepMC3Event(GenEvent * event)
{
        m_event=event;
        for(auto p: m_event->particles() )
        {
                PhotosParticle *particle = new PhotosHepMC3Particle(p);
                particles.push_back(particle);
        }
}

PhotosHepMC3Event::~PhotosHepMC3Event()
{
        while(particles.size())
        {
                PhotosParticle *p = particles.back();
                particles.pop_back();
                if(p) delete p;
        }
}

GenEvent * PhotosHepMC3Event::getEvent()
{
        return m_event;
}

void PhotosHepMC3Event::print()
{
        if(!m_event) return;
        Print::listing(*m_event);
}

vector<PhotosParticle*> PhotosHepMC3Event::getParticleList()
{
        return particles;
}

} // namespace Photospp
