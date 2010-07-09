#include <vector>
#include "PhotosParticle.h"
#include "PhotosBranch.h"
#include "Photos.h"
#include "Log.h"
using std::vector;
using std::endl;

PhotosBranch::PhotosBranch(PhotosParticle* p)
{
	if(p->getMothers().size()<=1)
	{
		// Regular case - one mother
		Log::Debug(1)<<"Regular case."<<endl;
		particle  = p->findLastSelf();
		mothers   = p->findProductionMothers();
		daughters = particle->getDaughters();
	}
	else
	{
		// Advanced case - branch without intermediate mother
		Log::Debug(1)<<"Advanced case."<<endl;
		particle  = NULL;
		mothers   = p->getMothers().at(0)->getDaughters();
		daughters = p->getDaughters();
	}
	valid = passBranchPointFilter();
	if(valid)
	{
		valid = checkMomentumConservation();
		if(!valid)
			Log::Warning()<<"Branching ignored due to 4-momentum non conservation"<<endl;
	}
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
	if(particle) return particle->checkMomentumConservation();
	return mothers.at(0)->checkMomentumConservation();
}

bool PhotosBranch::passBranchPointFilter()
{
	if(daughters.size()<2) return false;

	for(int i=0;i<(int)daughters.size();i++)
		if(daughters.at(i)->getPdgID()==22) return false;

	if(!Photos::supBremList) return true;

	int motherID = particle->getPdgID();
	vector< vector<int>* > &sup = *Photos::supBremList;

	vector<int> dID;
	for(int j=0;j<(int)daughters.size();j++) dID.push_back(daughters[j]->getPdgID());

	//Check if the mother and list of daughters matches any of the declared suppress patterns
	for(int j=0; j<(int)sup.size();j++)
	{
		if(motherID!=(*sup[j])[0]) continue;
		vector<int> &pList = *sup[j];
		bool fullMatch=true;
		for(int k = 1; k<(int)pList.size()-1; k++)
		{
			bool oneMatch=false;
			for(int l=0;l<(int)dID.size(); l++)
			{
				if(pList[k]==dID[l]) { oneMatch=true; break; }
			}
			if(!oneMatch) { fullMatch=false; break; }
		}
		if(pList.size()<=2 || (fullMatch && dID.size()==pList.size()-2))
		{
			//If the matching pattern is set for consecutive suppression - add particle to the list
			if(pList[pList.size()-1]==1)
			{
				if(!Photos::supParticles) Photos::supParticles = new vector< PhotosParticle* >();
				Photos::supParticles->push_back(particle);
			}
			return false;
		}
	}
	return true;
}
