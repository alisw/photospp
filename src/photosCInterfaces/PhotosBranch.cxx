#include <vector>
#include <list>
#include "PH_HEPEVT_Interface.h"
#include "PhotosParticle.h"
#include "PhotosBranch.h"
#include "Photos.h"
#include "Log.h"
using std::vector;
using std::list;
using std::endl;
typedef Photos::Log Log;

PhotosBranch::PhotosBranch(PhotosParticle* p)
{
	daughters = p->getDaughters();

	//Suppress if somehow got stable particle
	if(daughters.size()==0) suppression=1;
	
	if(daughters.at(0)->getMothers().size()==1)
	{
		// Regular case - one mother
		Log::Debug(1)<<"Regular case."<<endl;
		particle  = p;
		mothers   = p->findProductionMothers();
	}
	else
	{
		// Advanced case - branch with multiple mothers - no mid-particle
		Log::Debug(1)<<"Advanced case."<<endl;
		particle  = NULL;
		mothers   = daughters.at(0)->getMothers();
	}
	suppression = checkSuppressionLevel();
	if(!suppression && !checkMomentumConservation())
	{
		suppression=1;
		Log::Warning()<<"Branching ignored due to 4-momentum non conservation"<<endl;
	}
}

void PhotosBranch::process()
{
	int index = PH_HEPEVT_Interface::set(this);
	photos_make_(&index);
	PH_HEPEVT_Interface::get();
	checkMomentumConservation();
}

vector<PhotosParticle *> PhotosBranch::getParticles()
{
	vector<PhotosParticle *> ret = mothers;
	if(particle) ret.push_back(particle);
	ret.insert(ret.end(),daughters.begin(),daughters.end());
	return ret;
}

bool PhotosBranch::checkMomentumConservation()
{
	if(particle)         return particle->checkMomentumConservation();
	if(mothers.size()>0) return mothers.at(0)->checkMomentumConservation();
	return true;
}

vector<PhotosBranch *> PhotosBranch::createBranches(vector<PhotosParticle *> particles)
{
	list<PhotosParticle *> list(particles.begin(),particles.end());
	vector<PhotosBranch *> branches;
	while(!list.empty())
	{
		PhotosParticle *particle = list.front();
		list.pop_front();
		if(!particle) continue;

		PhotosBranch *branch = new PhotosBranch(particle);
		int suppression = branch->getSuppressionStatus();
		if(!suppression) branches.push_back(branch);
		else
		{
			//If suppressing consecutive decays
			if(suppression==2)
			{
				PhotosParticle *p = branch->getDecayingParticle();
				if(!p)
					if(branch->getMothers().size()>0) p = branch->getMothers().at(0);
					else continue;
				vector<PhotosParticle *> tree = p->getDecayTree();
				//Remove all particles from the list - max O(n*m)
				std::list<PhotosParticle *>::iterator it;
				for(it=list.begin();it!=list.end();it++)
				{
					for(int i=0;i<(int)tree.size();i++)
					{
						if(tree.at(i)->getBarcode()==(*it)->getBarcode())
						{
							it = list.erase(it);
							break;
						}
					}
					//In case we deleted the last particle
					if(it==list.end()) break;
				}
			}
			delete branch;
			continue;
		}
		Log::Debug(3)<<"Passed branch point filter"<<endl;

		//In case we don't have mid-particle erase rest of the mothers from list
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

int PhotosBranch::checkSuppressionLevel()
{
	if(!Photos::supBremList) return 0;
	
	int motherID;
	if(particle) motherID = particle->getPdgID();
	else
	{
		if(mothers.size()==0) return 0;
		motherID = mothers.at(0)->getPdgID();
	}

	vector< vector<int>* > &sup = *Photos::supBremList;
	vector<int> dID;
	for(int j=0;j<(int)daughters.size();j++) dID.push_back(daughters[j]->getPdgID());

	//Check if the mother and list of daughters matches any of the declared suppress patterns
	for(int j=0; j<(int)sup.size();j++)
	{
		//Skip patterns that don't have our mother
		if(motherID!=(*sup[j])[0]) continue;

		//Compare decay daughters with suppression pattern - max O(n*m)
		vector<int> &pList = *sup[j];
		bool fullMatch=true;
		for(int k = 1; k<(int)pList.size()-1; k++)
		{
			bool oneMatch=false;
			for(int l=0;l<(int)dID.size(); l++)
				if(pList[k]==dID[l]) { oneMatch=true; break; }
			if(!oneMatch) { fullMatch=false; break; }
		}
		if(pList.size()<=2 || (fullMatch && dID.size()==pList.size()-2))
		{
			//Check if the matching pattern is set for consecutive suppression
			return (pList.back()==1) ? 2 : 1;
		}
	}

	//Lastly, if not suppressed (so we are sure status!=2)
	//check for decay with photon already present
	for(int i=0;i<(int)daughters.size();i++)
		if(daughters.at(i)->getPdgID()==22) return 1;

	return 0;
}

