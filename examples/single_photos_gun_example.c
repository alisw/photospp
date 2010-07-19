/**
 * Example of use of processParticle routine.
 * Pythia events are generated and photos used on first tau+ found.
 *
 * @author Tomasz Przedzinski
 * @date 17 July 2010
 */

//pythia header files
#include "Pythia.h"
#include "HepMCInterface.h"

//PHOTOS header files
#include "Photos.h"
#include "PhotosHepMCParticle.h"
#include "Log.h"
typedef Photos::Log Log; //We're using Photos version of Log class

using namespace std;
using namespace Pythia8;

int main(int argc,char **argv)
{
	// Initialisation of pythia
	HepMC::I_Pythia8 ToHepMC;
	Pythia pythia;
	Event& event = pythia.event;
	pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
	pythia.readString("23:onMode = off");
	pythia.readString("23:onIfAny = 15");
	pythia.readString("HadronLevel:Hadronize = off");
	pythia.readString("SpaceShower:QEDshower = off");
	pythia.readString("SpaceShower:QEDshowerByL = off");
	pythia.readString("SpaceShower:QEDshowerByQ = off");
	pythia.readString("PartonLevel:ISR = off");
	pythia.readString("PartonLevel:FSR = off");
	pythia.init( 11, -11, 92.);

	Photos::initialize();
	Photos::setInfraredCutOff(0.001/200);

	int NumberOfEvents = 10000;
	if(argc>1) NumberOfEvents=atoi(argv[1]);

	int photonAdded=0,twoAdded=0,moreAdded=0,tauCount=0;
	// Begin event loop. Generate event.
	for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent)
	{
		if(iEvent%(NumberOfEvents/10)==0) Log::Info()<<iEvent<<endl;
		if(!pythia.next()) continue;

		HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();
		HepMCEvt->use_units(HepMC::Units::GEV,HepMC::Units::MM);
		ToHepMC.fill_next_event(event, HepMCEvt);

		// Find tau
		HepMC::GenParticle *tau=0;
		for(HepMC::GenEvent::vertex_const_iterator i = HepMCEvt->vertices_begin();i!=HepMCEvt->vertices_end();i++)
		{
			for(HepMC::GenVertex::particles_in_const_iterator p=(*i)->particles_in_const_begin();p!=(*i)->particles_in_const_end(); p++)
			{
				if((*p)->pdg_id()==15) tau=*p;
				break;
			}
			if(tau) break;
		}
		if(tau)
		{
			tauCount++;
			int buf = -HepMCEvt->particles_size();

			// Call photos
			Photos::processParticle( new PhotosHepMCParticle(tau) );

			buf+=HepMCEvt->particles_size();
			if(buf==1)      photonAdded++;
			else if(buf==2) twoAdded++;
			else if(buf>2)  moreAdded++;
		}

		//clean up
		delete HepMCEvt;
	}
	pythia.statistics();
	cout.precision(2);
	cout.setf(ios::fixed);
	cout<<endl;
	if(tauCount==0)
	{
		cout<<"Something went wrong with pythia generation."<<endl;
		cout<<"No taus were processed."<<endl<<endl;
		return 0;
	}
	cout<<"Summary (single tau decay processing):"<<endl;
	cout<<tauCount   <<"\ttaus processed"<<endl;
	cout<<photonAdded<<"\ttimes one photon added to the decay           \t("<<(photonAdded*100./tauCount)<<"%)"<<endl;
	cout<<twoAdded   <<"\ttimes two photons added to the decay          \t("<<(twoAdded*100./tauCount)<<"%)"<<endl;
	cout<<moreAdded  <<"\ttimes more than two photons added to the decay\t("<<(moreAdded*100./tauCount)<<"%)"<<endl<<endl;
	cout<<"(Contrary to results from MC-Tester, these values are technical and infrared unstable)"<<endl<<endl;
	cout<<"To proccess different number of events use:"<<endl<<" ./single_photos_gun_example <number_of_events>"<<endl<<endl;
}

