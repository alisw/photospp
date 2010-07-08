#include <vector>
#include "PH_HEPEVT_Interface.h"
#include "PhotosParticle.h"
#include "PhotosEvent.h"
#include "Log.h"
using std::vector;

PhotosEvent::~PhotosEvent()
{
	while(m_branch_points.size()!=0)
	{
		PhotosParticle *temp = m_branch_points.back();
		m_branch_points.pop_back();
		delete temp;
	}
}

void PhotosEvent::process()
{
	//print();
	int index = 2;
	m_branch_points = getBranchPoints();
	
	vector<PhotosParticle *> branch_points = filterBranchPoints();

	for(int i=0;i<(int)branch_points.size();i++)
	{	
		PH_HEPEVT_Interface::set(branch_points.at(i));
		photos_make_(&index);
		PH_HEPEVT_Interface::get();
		branch_points.at(i)->checkMomentumConservation();
	}
	//print();
}

vector<PhotosParticle *> PhotosEvent::filterBranchPoints()
{
	vector<PhotosParticle *> ret;
	for(int i=0;i<(int)m_branch_points.size();i++)
	{
		PhotosParticle *p = m_branch_points.at(i);
		if(!passBranchPointFilter(p)) continue;
		if(!passSuppressionFilter(p)) continue;
		if(p->checkMomentumConservation())
			ret.push_back(p);
		else
			Log::Warning()<<"Branching ignored due to 4-momentum non conservation"<<std::endl;
	}
	return ret;
}

bool PhotosEvent::passBranchPointFilter(PhotosParticle *particle)
{
	//check that the particle decays
	if(particle->getStatus()!=PhotosParticle::DECAYED) return false;
	//check for self decays
	if(particle->getDaughters().size()==1) return false;
	return true;
}

bool PhotosEvent::passSuppressionFilter(PhotosParticle *particle)
{
	if(!Photos::supBremList) return true;

	//Check if the particle is suppressed via consecutive suppression of one of its mothers
	if(Photos::supParticles && !passSuppressConsecutive(particle)) return false;

	int motherID = particle->getPdgID();
	vector< vector<int>* > &sup = *Photos::supBremList;

	//Search for last self to get daughters list
	vector< PhotosParticle* > daughters = particle->getDaughters();
	for(int j=0;j<(int)daughters.size();j++)
	{
		if(daughters[j]->getPdgID()==motherID)
		{
			daughters = daughters[j]->getDaughters();
			j=-1;
			continue;
		}
	}

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
