#include <vector>
#include "PhotosBranch.h"
#include "PhotosParticle.h"
#include "PH_HEPEVT_Interface.h"
#include "Log.h"

using namespace std;

vector<PhotosParticle*> PH_HEPEVT_Interface::m_particle_list;

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
  ph_hepevt_.phep[i][4]=particle->getMass();

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
	int decay_idx=0;
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

	for(int i=0;i<nmothers;i++)
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

	//phodmp_();
	//printf("AND RETURNING: %i\n",(decay_idx) ? decay_idx : 1 );
	return (decay_idx) ? decay_idx : 1;
}

void PH_HEPEVT_Interface::get(){

  int index = 0;

  //if no photons have been added to the event record, do nothing.
  if(ph_hepevt_.nhep == (int) m_particle_list.size())
    return;

  //phodmp_();

  //otherwise loop over particles which are already in the
  //event record and modify their 4 momentum
  for(; index < ph_hepevt_.nhep && index < (int) m_particle_list.size(); index++){

    PhotosParticle * particle = m_particle_list.at(index);

    if(ph_hepevt_.idhep[index]!=particle->getPdgID())
      Log::Fatal("PH_HEPEVT_Interface::get(): Something is wrong with the PH_HEPEVT common block",5);

    //check to see if this particle's 4-momentum has been modified
    if(ph_hepevt_.phep[index][0]!=particle->getPx()||
       ph_hepevt_.phep[index][1]!=particle->getPy()||
       ph_hepevt_.phep[index][2]!=particle->getPz()){
      
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


  //Now add extra photons
  int photons = ph_hepevt_.nhep - m_particle_list.size();
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
  
}
