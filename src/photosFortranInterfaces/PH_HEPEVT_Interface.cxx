#include <vector>
#include <cmath>
#include "PhotosBranch.h"
#include "PhotosParticle.h"
#include "PH_HEPEVT_Interface.h"
#include "Log.h"
using namespace std;

namespace Photospp
{

vector<PhotosParticle*> PH_HEPEVT_Interface::m_particle_list;
int PH_HEPEVT_Interface::ME_channel=0;
int PH_HEPEVT_Interface::decay_idx=0;

void PH_HEPEVT_Interface::clear(){

  m_particle_list.clear();

  ph_hepevt_.nevhep=0; 
  ph_hepevt_.nhep=0;
  

  /**  for(int i=0; i < NMXHEP; i++){

    ph_hepevt_.isthep[i]=0;
    ph_hepevt_.idhep[i]=0;
    
    for(int j=0; j<2; j++){
      ph_hepevt_.jmohep[i][j]=0;
      ph_hepevt_.jdahep[i][j]=0;
    }
    
    for(int j=0; j<5; j++)
      ph_hepevt_.phep[i][j]=0;
    
    for(int j=0; j<4; j++)
      ph_hepevt_.vhep[i][j]=0;
  
      ph_phoqed_.qedrad[i]=0;
  
      }**/
}

void PH_HEPEVT_Interface::add_particle(int i,PhotosParticle * particle,
				       int first_mother, int last_mother,
				       int first_daughter, int last_daughter){

  if(i>0)
    i--; //account for fortran indicies begining at 1
  else
    Log::Warning()<<"Index given to PH_HEPEVT_Interface::add_particle "
	          <<"is too low (it must be > 0)."<<endl;

  //add to our internal list of pointer/index pairs
  m_particle_list.push_back(particle);

  //now set the element of PH_HEPEVT
  ph_hepevt_.nevhep=0; //dummy
  ph_hepevt_.nhep=ph_hepevt_.nhep++;
  ph_hepevt_.isthep[i]=particle->getStatus();
  ph_hepevt_.idhep[i]=particle->getPdgID();

  ph_hepevt_.jmohep[i][0]=first_mother;
  ph_hepevt_.jmohep[i][1]=last_mother;

  ph_hepevt_.jdahep[i][0]=first_daughter;
  ph_hepevt_.jdahep[i][1]=last_daughter;

  ph_hepevt_.phep[i][0]=particle->getPx();
  ph_hepevt_.phep[i][1]=particle->getPy();
  ph_hepevt_.phep[i][2]=particle->getPz();
  ph_hepevt_.phep[i][3]=particle->getE();
  
  // if massFrom4Vector=true (default) - get sqrt(e^2-p^2)
  // otherwise - get mass from event record
  if(!Photos::massFrom4Vector) ph_hepevt_.phep[i][4]=particle->getMass();
  else
  {
    double mass  = 0.0;
    int    pdgid = abs(particle->getPdgID());

    // if 'forceMass' for this PDGID was used - overwrite mass
    if(Photos::forceMassList)
    {
      for(unsigned int i=0;i<Photos::forceMassList->size();i++)
      {
        if(pdgid == Photos::forceMassList->at(i)->first)
        {
          mass = Photos::forceMassList->at(i)->second;
          
          // when 'forceMassFromEventRecord' is used mass is -1.0
          // in this case - get mass from event record
          if(mass<0.0) mass = particle->getMass();

          break;
        }
      }
    }

    // if PDGID was not found on the list - calculate from 4-vector
    if(mass==0.0) mass = particle->getVirtuality();

    ph_hepevt_.phep[i][4] = mass;
  }

  ph_hepevt_.vhep[i][0]=0;
  ph_hepevt_.vhep[i][1]=0;
  ph_hepevt_.vhep[i][2]=0;
  ph_hepevt_.vhep[i][3]=0;

  ph_phoqed_.qedrad[i]=1;

}

int PH_HEPEVT_Interface::set(PhotosBranch *branch)
{
	PH_HEPEVT_Interface::clear();
	int idx=1;

	//get mothers
	vector<PhotosParticle *> mothers = branch->getMothers();
	int nmothers=mothers.size();

	//check if mid-particle exist
	decay_idx=0;
	PhotosParticle *decay_particle = branch->getDecayingParticle();
	if(decay_particle) decay_idx=nmothers+1;

	//get daughters
	vector<PhotosParticle *> daughters = branch->getDaughters();
	int ndaughters=daughters.size();

	for(int i=0;i<nmothers;i++)
	{
		if(decay_idx)
			add_particle(idx++,mothers.at(i),
			             0,0, //mothers
			             decay_idx,decay_idx); //daughters
		else
			add_particle(idx++,mothers.at(i),
			             0,0, //mothers
			             nmothers+1,nmothers+ndaughters); //daughters
	}

	if(decay_particle)
		add_particle(idx++,decay_particle,
		             1,nmothers, //mothers
		             nmothers+2,nmothers+1+ndaughters); //daughters

	for(int i=0;i<ndaughters;i++)
	{
		if(decay_idx)
			add_particle(idx++,daughters.at(i),
			             decay_idx,decay_idx, //mothers
			             0,0); //daughters
		else
			add_particle(idx++,daughters.at(i),
			             1,nmothers, //mothers
			             0,0); //daughters
	}
	//Log::RedirectOutput( phodmp_ , Log::Debug(1000) );
	Log::Debug(1000,false)<<"PH_HEPEVT returning: "<<( (decay_idx) ? decay_idx : 1 )<<" from "<<idx-1<<" particles."<<endl;
	return (decay_idx) ? decay_idx : 1;
}

void PH_HEPEVT_Interface::get(){

  int index = 0;

  //if no photons have been added to the event record, do nothing.
  if(ph_hepevt_.nhep == (int) m_particle_list.size())
    return;

  //phodmp_();

  int  particle_count  = m_particle_list.size();
  int  daughters_start = ph_hepevt_.jmohep[ph_hepevt_.nhep-1][0];
  int  photons         = ph_hepevt_.nhep - m_particle_list.size();
  bool isPhotonCreated = (photons>0);

  index = particle_count;

  // Add extra photons
  for(;photons>0; photons--, index++){
    
    if(ph_hepevt_.idhep[index]!=PhotosParticle::GAMMA)
      Log::Fatal("PH_HEPEVT_Interface::get(): Extra particle added to the PH_HEPEVT common block in not a photon!",6);
    
    //create a new particle
    PhotosParticle * new_photon;
    new_photon = m_particle_list.at(0)->createNewParticle(ph_hepevt_.idhep[index],
							  ph_hepevt_.isthep[index],
							  ph_hepevt_.phep[index][4],
							  ph_hepevt_.phep[index][0],
							  ph_hepevt_.phep[index][1],
							  ph_hepevt_.phep[index][2],
							  ph_hepevt_.phep[index][3]);
    
    //add into the event record
    //get mother particle of photon
    PhotosParticle * mother =  m_particle_list.at(ph_hepevt_.jmohep[index][0]-1);
    mother->addDaughter(new_photon);
    
  }

  //otherwise loop over particles which are already in the
  //event record and modify their 4 momentum
  //4.03.2012: Fix to prevent kinematical trap in vertex of simultaneous:
  //           z-collinear and non-conservation pf E,p for dauthters of grandmothers
  for(index=daughters_start; index < particle_count && index < (int) m_particle_list.size(); index++){

    PhotosParticle * particle = m_particle_list.at(index);

    if(ph_hepevt_.idhep[index]!=particle->getPdgID())
      Log::Fatal("PH_HEPEVT_Interface::get(): Something is wrong with the PH_HEPEVT common block",5);

    // If photons were added - for each daughter create a history entry
    if(isPhotonCreated && Photos::isCreateHistoryEntries)
    {
      particle->createHistoryEntry();
    }
    
    //check to see if this particle's 4-momentum has been modified
    bool   update=false;

    // don't update particle if difference lower than THRESHOLD * particle energy (default threshold = 10e-8)
    double threshold = NO_BOOST_THRESHOLD*ph_hepevt_.phep[index][3];
    if( fabs(ph_hepevt_.phep[index][0]-particle->getPx()) > threshold ||
        fabs(ph_hepevt_.phep[index][1]-particle->getPy()) > threshold ||
        fabs(ph_hepevt_.phep[index][2]-particle->getPz()) > threshold ||
        fabs(ph_hepevt_.phep[index][3]-particle->getE())  > threshold    ) update=true;

    if(update)
    {
      
      //modify this particle's momentum and it's daughters momentum
      //Steps 1., 2. and 3. must be executed in order.

      //1. boost the particles daughters into it's (old) rest frame
      particle->boostDaughtersToRestFrame(particle);

      //2. change this particles 4 momentum
      particle->setPx(ph_hepevt_.phep[index][0]);
      particle->setPy(ph_hepevt_.phep[index][1]);
      particle->setPz(ph_hepevt_.phep[index][2]);
      particle->setE(ph_hepevt_.phep[index][3]);

      //3. boost the particles daughters back into the lab frame
      particle->boostDaughtersFromRestFrame(particle);
    }

  }
  
}

void PH_HEPEVT_Interface::prepare()
{
	check_ME_channel();
}

void PH_HEPEVT_Interface::complete()
{

}

void PH_HEPEVT_Interface::check_ME_channel()
{
	ME_channel=0;

// Check mothers:

	if(decay_idx==2)                              return; // Only one mother present
	if(ph_hepevt_.idhep[0]*ph_hepevt_.idhep[1]>0) return; // Mothers have same sign

	Log::Debug(900)<<"ME_channel: Mothers PDG:  "<<ph_hepevt_.idhep[0]<<" "<<ph_hepevt_.idhep[1]<<endl;
	if(decay_idx)
		Log::Debug(900,false)<<"            Intermediate: "<<ph_hepevt_.idhep[decay_idx-1]<<endl;

	int              firstDaughter=3;
	if(decay_idx==0) firstDaughter=2; // if no intermediate particle - daughters start at idx 2

	// Are mothers in range +/- 1-6; +/- 11-16?
	int mother1 = abs(ph_hepevt_.idhep[0]);
	int mother2 = abs(ph_hepevt_.idhep[1]);
	if( mother1<1 || (mother1>6 && mother1<11) || mother1>16 ) return;
	if( mother2<1 || (mother2>6 && mother2<11) || mother2>16 ) return;

//Check daughters

	// Z: check for pairs 11 -11 ; 13 -13 ; 15 -15
	// -------------------------------------------
	int firstPDG =0;
	int secondPDG=0;
	for(int i=firstDaughter; i<ph_hepevt_.nhep;i++)
	{
		int pdg = abs(ph_hepevt_.idhep[i]);
		if(pdg==11 || pdg==13 || pdg==15)
		{
			if(firstPDG==0) firstPDG=ph_hepevt_.idhep[i];
			else
			{
				secondPDG=ph_hepevt_.idhep[i];
				// Just in case two pairs are genereted - verify that we have a pair with oposite signs
				if(firstPDG*secondPDG>0) secondPDG=0;
				break;
			}
		}
	}

	if( ME_channel==0 && firstPDG!=0 && secondPDG!=0 &&
	    firstPDG==-secondPDG ) ME_channel=1;

	// W: check for pairs 11 -12; -11 12; 13 -14; -13 14; 15 -16; -15 16
	// -----------------------------------------------------------------
	firstPDG =0;
	secondPDG=0;
	for(int i=firstDaughter; i<ph_hepevt_.nhep;i++)
	{
		int pdg = abs(ph_hepevt_.idhep[i]);
		if(pdg>=11 && pdg<=16)
		{
			if(firstPDG==0) firstPDG=ph_hepevt_.idhep[i];
			else
			{
				secondPDG=ph_hepevt_.idhep[i];
				// Just in case two pairs are genereted - verify that we have a pair with oposite signs
				if(firstPDG*secondPDG>0) secondPDG=0;
				break;
			}
		}
	}

	firstPDG =abs(firstPDG);
	secondPDG=abs(secondPDG);

	if(  ME_channel==0 && firstPDG!=0 && secondPDG!=0 &&
	   ( ( firstPDG==11 && secondPDG==12 ) || (firstPDG == 12 && secondPDG == 11) ||
	     ( firstPDG==13 && secondPDG==14 ) || (firstPDG == 14 && secondPDG == 13) ||
	     ( firstPDG==15 && secondPDG==16 ) || (firstPDG == 16 && secondPDG == 15)
	   )
	  ) ME_channel=2;

	Log::Debug(901)<<"ME_channel: Found ME_channel: "<<ME_channel<<endl;

// Check intermediate particle (if exists):

	// Verify that intermediate particle PDG matches ME_channel found
	if(ME_channel>0 && decay_idx)
	{
		int pdg=ph_hepevt_.idhep[decay_idx-1];

		if(ME_channel==1 && !(pdg==22 || pdg==23) ) ME_channel=0; //gamma/Z
		if(ME_channel==2 && !(pdg==24 || pdg==-24)) ME_channel=0; //W+/W-

		if(ME_channel==0)
			Log::Debug(901,false)<<"            but set to 0: wrong intermediate particle: "<<pdg<<endl;
	}

// Check flags

	switch(ME_channel)
	{
		case  0: break; // Ok - no channel found
		case  1: if(!Photos::meCorrectionWtForZ) ME_channel=0; break;
		case  2: if(!Photos::meCorrectionWtForW) ME_channel=0; break;
		default: Log::Error()<<"Unimplemented ME channel: "<<ME_channel<<endl; break;
	}
	Log::Debug(902)<<"ME_channel: Finally, after flag check, ME_channel is: "<<ME_channel<<endl;
}

// Call from fortran: 'CALL ME_CHANNEL(X)'
extern "C" void me_channel_(int *x)
{
	*x=PH_HEPEVT_Interface::ME_channel;
}
// Call from fortran: 'CALL ME_SCALAR(X)'
// transmits flag if NLO correction for scalar is on/off
extern "C" void me_scalar_(int *x)
{
  *x=(int) Photos::meCorrectionWtForScalar;
}

} // namespace Photospp
