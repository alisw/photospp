#include <vector>
#include <list>
#include "PH_HEPEVT_Interface.h"
#include "PhotosParticle.h"
#include "PhotosBranch.h"
#include "PhotosEvent.h"
#include "Log.h"
using std::vector;
using std::list;
typedef Photos::Log Log;

PhotosEvent::~PhotosEvent()
{
	while(m_branch_points.size()!=0)
	{
		PhotosBranch *temp = m_branch_points.back();
		m_branch_points.pop_back();
		delete temp;
	}
}

void PhotosEvent::process()
{
	//print();
	vector<PhotosParticle *> particles = getParticleList();

	m_branch_points = createBranches(particles);

	for(int i=0;i<(int)m_branch_points.size();i++)
	{
		int index = PH_HEPEVT_Interface::set(m_branch_points.at(i));
		photos_make_(&index);
		PH_HEPEVT_Interface::get();
		m_branch_points.at(i)->checkMomentumConservation();
	}
	//print();
}

vector<PhotosBranch *> PhotosEvent::createBranches(vector<PhotosParticle *> particles)
{
	list<PhotosParticle *> list(particles.begin(),particles.end());
	vector<PhotosBranch *> branches;
	while(!list.empty())
	{
		PhotosParticle *particle = list.front();
		list.pop_front();
		if(!particle) continue;
		if(!passParticleFilter(particle))
		{
			delete particle;
			continue;
		}
		Log::Debug(2)<<"Passed particle filter"<<endl;
		PhotosBranch *branch = new PhotosBranch(particle);
		if(branch->isValid()) branches.push_back(branch);
		else
		{
			delete branch;
			continue;
		}
		Log::Debug(3)<<"Passed branch point filter"<<endl;

		//In case we don't have mid-particle
		if(!branch->getDecayingParticle())
		{
			vector<PhotosParticle *> mothers = branch->getMothers();
			for(int i=0;i<(int)mothers.size();i++)
			{
				PhotosParticle *m = mothers.at(i);
				if(m->getBarcode()==particle->getBarcode()) continue;
				std::list<PhotosParticle *>::iterator it;
				for(it=list.begin();it!=list.end();it++)
					if(m->getBarcode()==(*it)->getBarcode())
					{
						it = list.erase(it);
						break;
					}
			}
		}
	}
	return branches;
}

bool PhotosEvent::passParticleFilter(PhotosParticle *particle)
{
	//check that the particle decays
	if(particle->getStatus()!=PhotosParticle::DECAYED) return false;

	//check for self decays
	if(particle->getDaughters().size()==1) return false;

	//Check if the particle is suppressed via consecutive suppression of one of its mothers
	if(Photos::supParticles && !passSuppressConsecutive(particle)) return false;

	return true;
}

bool PhotosEvent::passSuppressConsecutive(PhotosParticle *particle)
{
	vector<PhotosParticle*> mothers = particle->getMothers();
	if(mothers.size()==0) return true;

	//Recursive check of particles' mothers
	for(int i=0;i<(int)mothers.size();i++)
		for(int j=0;j<(int)Photos::supParticles->size();j++)
		{
			PhotosParticle* mother = mothers.at(i);
			if(mother->getBarcode()==Photos::supParticles->at(j)->getBarcode()) return false;
			if(passSuppressConsecutive(mother)==false) return false;
		}
	return true;
}
