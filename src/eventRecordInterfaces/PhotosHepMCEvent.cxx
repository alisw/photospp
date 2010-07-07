#include "PhotosHepMCEvent.h"
#include "Photos.h"
#include "Log.h"

using namespace std;

PhotosHepMCEvent::PhotosHepMCEvent(HepMC::GenEvent * event){
  m_event=event;
};

PhotosHepMCEvent::~PhotosHepMCEvent(){

  while(m_branch_points.size()!=0){
    PhotosParticle * temp = m_branch_points.back();
    m_branch_points.pop_back();
    delete temp;
  }
}

HepMC::GenEvent * PhotosHepMCEvent::getEvent(){
  return m_event;
}

std::vector<PhotosParticle*> PhotosHepMCEvent::getBranchPoints(){

  //if the branch point vector has not been set yet, then fill it
  if(m_branch_points.size()==0){
    
    HepMC::GenEvent::particle_const_iterator part_itr = m_event->particles_begin();
    //loop over all particle in the event looking for a branch to give to PHOTOS
    for( ; part_itr!=m_event->particles_end(); part_itr++)
    {
	if(!passBranchPointFilter(*part_itr)) continue;
	PhotosHepMCParticle * particle = new PhotosHepMCParticle(*part_itr);
	if(!passSuppressionFilter(particle)) continue;
	if(particle->checkMomentumConservation())
	  m_branch_points.push_back(particle);
	else
	  Log::Warning()<<"Branching being ignored due to 4-momentum non conservation"<<endl;

    }
  }
  return m_branch_points;
}

bool PhotosHepMCEvent::passBranchPointFilter(HepMC::GenParticle * particle){

  //check that the particle decays
  if(!particle->end_vertex() || particle->end_vertex()->particles_out_size()==0)
    return false;

  //check for self decays
  if(particle->end_vertex()->particles_out_size()==1)
    return false;


  //check for more complicated scenarios
      /**  HepMC::GenVertex::particles_out_const_iterator pcle_itr;
  pcle_itr = particle->end_vertex()->particles_out_const_begin();
  HepMC::GenVertex::particles_out_const_iterator pcle_itr_end;
  pcle_itr_end = particle->end_vertex()->particles_out_const_end();
  
  for(;pcle_itr!=pcle_itr_end; pcle_itr++){
    if()

    }**/


  return true;
      
}

bool PhotosHepMCEvent::passSuppressionFilter(PhotosHepMCParticle *particle)
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
				if(!Photos::supParticles) Photos::supParticles = new vector< HepMC::GenParticle* >();
				Photos::supParticles->push_back(particle->getHepMC());
			}
			return false;
		}
	}
	return true;
}

bool PhotosHepMCEvent::passSuppressConsecutive(PhotosHepMCParticle *particle)
{
	vector< PhotosParticle* > mothers = particle->getMothers();
	if(mothers.size()==0) return true;

	//Recursive check of particles' mothers

	for(int i=0;i<(int)mothers.size();i++)
		for(int j=0;j<(int)Photos::supParticles->size();j++)
		{
			PhotosHepMCParticle* hepMother = (PhotosHepMCParticle*)mothers[i];
			if(hepMother->getHepMC()==Photos::supParticles->at(j)) return false;
			if(passSuppressConsecutive(hepMother)==false) return false;
		}
	return true;
}
