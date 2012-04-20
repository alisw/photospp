/**
 * Example of use of Photos C++ interface.
 * e+, e- -> tau + tau - HEPEVT events are constructed.
 * Taus are subsequently decayed via Photos.
 *
 * @author Tomasz Przedzinski
 * @date 24 November 2011
 */

#include "Photos/Photos.h"
#include "Photos/PhotosHEPEVTParticle.h"
#include "Photos/PhotosHEPEVTEvent.h"
using namespace std;
using namespace Photospp;

/** Create a simple e+ + e- -> Z -> tau+ tau- HEPEVT event **/
PhotosHEPEVTEvent* make_simple_tau_event(){

  PhotosHEPEVTEvent * evt = new PhotosHEPEVTEvent();

  const double amell = 0.0005111;
  const double amtau = 1.777;

  // Create some four vectors for the electrons
  double e_mass_sq   = amell*amell;
  double tau_mass_sq = amtau*amtau;

  double e1_pz = -2.0; //change these
  double e2_pz =  3.5; //as needed
  double e1_e  = sqrt(e1_pz*e1_pz + e_mass_sq);
  double e2_e  = sqrt(e2_pz*e2_pz + e_mass_sq);

  // Make PhotosParticles for boosting
  PhotosHEPEVTParticle *first_e      = new PhotosHEPEVTParticle(-11, 3, 0., 0., e1_pz,      e1_e,     amell, -1, -1,  2,  2);
  PhotosHEPEVTParticle *second_e     = new PhotosHEPEVTParticle( 11, 3, 0., 0., e2_pz,      e2_e,     amell, -1, -1,  2,  2);
  PhotosHEPEVTParticle *intermediate = new PhotosHEPEVTParticle( 23, 3, 0., 0., e1_pz+e2_pz,e1_e+e2_e,0.,     0,  1,  3,  4);
  PhotosHEPEVTParticle *first_tau    = new PhotosHEPEVTParticle(-15, 1, 0., 0., 0.,         0.,       amtau,  2, -1, -1, -1);
  PhotosHEPEVTParticle *second_tau   = new PhotosHEPEVTParticle( 15, 1, 0., 0., 0.,         0.,       amtau,  2, -1, -1, -1);

  // Order matters!
  evt->addParticle(first_e     );
  evt->addParticle(second_e    );
  evt->addParticle(intermediate);
  evt->addParticle(first_tau   );
  evt->addParticle(second_tau  );

  double tau_energy = 0.5*sqrt( (e1_e+e2_e)*(e1_e+e2_e) - (e1_pz+e2_pz)*(e1_pz+e2_pz) );

  first_tau->setE  (tau_energy);
  first_tau->setPx ((1.0/sqrt(2.0))*sqrt(tau_energy*tau_energy-tau_mass_sq));
  first_tau->setPy ((1.0/sqrt(2.0))*sqrt(tau_energy*tau_energy-tau_mass_sq));

  second_tau->setE (tau_energy);
  second_tau->setPx(-1*(1.0/sqrt(2.0))*sqrt(tau_energy*tau_energy-tau_mass_sq));
  second_tau->setPy(-1*(1.0/sqrt(2.0))*sqrt(tau_energy*tau_energy-tau_mass_sq));

  // Boost particles from rest frame
  first_tau ->isInRestFrame = true;
  second_tau->isInRestFrame = true;

  first_tau ->boostFromRestFrame(intermediate);
  second_tau->boostFromRestFrame(intermediate);

  return evt;
}

/** Example of using Photos to process event stored in HEPEVT event record */
int main(void){

  int NumberOfEvents = 100000;

  Photos::initialize();

	int photonAdded=0,twoAdded=0,moreAdded=0,evtCount=0;

  // Begin event loop. Generate event.
  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {

    if(iEvent%10000==0) cout<<"Event: "<<iEvent<<"\t("<<iEvent*(100./NumberOfEvents)<<"%)"<<endl;

    // Create simple event
    PhotosHEPEVTEvent * event = make_simple_tau_event();

    int buf = -event->getParticleCount();

    //cout << "BEFORE:"<<endl;
    //event->print();

    event->process();

    buf += event->getParticleCount();
		if     (buf==1) photonAdded++;
		else if(buf==2) twoAdded++;
		else if(buf>2)  moreAdded++;

    //cout << "AFTER:"<<endl;
    //event->print();

    evtCount++;

    //clean up
    delete event;
  }

	// Print results
	cout.precision(3);
	cout.setf(ios::fixed);
	cout<<endl;
	cout<<"Summary of processing simple events:    e+ e- -> Z -> tau+ tau-"<<endl;
	cout<<evtCount   <<"\tevents processed"<<endl;
	cout<<photonAdded<<"\ttimes one photon added to the event           \t("<<(photonAdded*100./evtCount)<<"%)"<<endl;
	cout<<twoAdded   <<"\ttimes two photons added to the event          \t("<<(twoAdded*100./evtCount)<<"%)"<<endl;
	cout<<moreAdded  <<"\ttimes more than two photons added to the event\t("<<(moreAdded*100./evtCount)<<"%)"<<endl<<endl;
	cout<<"(Contrary to results from MC-Tester, these values are technical and infrared unstable)"<<endl<<endl;
}

