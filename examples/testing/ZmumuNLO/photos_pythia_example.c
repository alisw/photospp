/**
 * Example of use of photos C++ interface. Pythia events are
 * generated first and photos used for FSR.
 *
 * @author Nadia Davidson
 * @date 6 July 2009
 */

//pythia header files
#include "Pythia.h"
#include "HepMCInterface.h"

//MC-TESTER header files
#include "Generate.h"
#include "HepMCEvent.H"
#include "Setup.H"

//PHOTOS header files
#include "Photos.h"
#include "PhotosHepMCEvent.h"
#include "Log.h"
typedef Photos::Log Log; //We're using Photos version of Log class

using namespace std;
using namespace Pythia8;

bool ShowersOn=true;
unsigned long NumberOfEvents = 10000;

// Finds X Y -> 6 -6 decay and converts it to 100 -> 6 -6, where 100 = X + Y
void fixForMctester(HepMC::GenEvent *evt)
{
	for(HepMC::GenEvent::particle_const_iterator p=evt->particles_begin();p!=evt->particles_end(); p++)
	if((*p)->pdg_id()==6)
	{
		HepMC::GenParticle *pt = *p;
		int id=(* pt->production_vertex()->particles_in_const_begin() )->pdg_id();
		if(id!=21 && id!=11 && id>5) continue;

		// Get first mother and add 2x second mother to it
		HepMC::GenParticle *X = (* pt->production_vertex()->particles_in_const_begin());
		HepMC::GenParticle *Y = (* ++(pt->production_vertex()->particles_in_const_begin()) );
		HepMC::FourVector fX = X->momentum();
		HepMC::FourVector fY = Y->momentum();
		HepMC::FourVector fXY(fX.px()+fY.px(),fX.py()+fY.py(),fX.pz()+fY.pz(),fX.e()+fY.e());
		X->set_momentum(fXY);
		// Unique ID for MC-Tester to analyze
		X->set_pdg_id(100);

		// Set 2nd mother as decayed and delete it from production vertex
		Y->set_status(1);
		(* Y->production_vertex()->particles_in_const_begin())->set_status(1);
		pt->production_vertex()->remove_particle(Y);
		return;
	}
}

int main(int argc,char **argv)
{
	// Initialisation of pythia
	HepMC::I_Pythia8 ToHepMC;
	Pythia pythia;
	Event& event = pythia.event;
	//pythia.settings.listAll(); exit(0);

	// Console input parameters
	// (set by examples located in 'testing' directory)
	bool topDecays =false;
	bool zeeDecays =false;
	bool zmuDecays =false;
	bool zNLO      =false;
	if(argc>4)
	{
		// Advanced options
		topDecays = (atoi(argv[4])==1);
		zeeDecays = (atoi(argv[4])==2);
		zmuDecays = (atoi(argv[4])==3);
	}
	if(argc>3) NumberOfEvents=atol(argv[3]);
	if(argc>2) ShowersOn=atoi(argv[2]);

	if(true)
	{
		//pythia.readString("HadronLevel:all = off");
		pythia.readString("HadronLevel:Hadronize = off");
		pythia.readString("SpaceShower:QEDshower = off");
		pythia.readString("SpaceShower:QEDshowerByL = off");
		pythia.readString("SpaceShower:QEDshowerByQ = off");
	}
	pythia.readString("PartonLevel:ISR = off");
	pythia.readString("PartonLevel:FSR = off");
	
	// This one produces gamma in 11->11 decays
	pythia.readString("PartonLevel:Remnants = off");

	/*
	pythia.readString("TimeShower:QCDshower = off");
	pythia.readString("TimeShower:QEDshowerByGamma = off");
	pythia.readString("TimeShower:QEDshowerByL = off");
	pythia.readString("TimeShower:QEDshowerByQ = off");
	pythia.readString("SpaceShower:QCDshower = off");
	pythia.readString("SpaceShower:QEDshowerByL = off");
	pythia.readString("SpaceShower:QEDshowerByQ = off");
	pythia.readString("ParticleDecays:FSRinDecays = off");
	pythia.readString("PartonLevel:FSRinProcess = off");
	pythia.readString("PartonLevel:FSRinResonances = off");
	pythia.readString("PartonLevel:MI = off");
	*/
	
	if(argc>1)  //pre-set configs
	{
		pythia.readFile(argv[1]);
		if(zeeDecays || topDecays) pythia.init( -2212, -2212, 14000.0); //p  p  collisions
		else if(zmuDecays)         pythia.init( 11, -11, 91.17);        //e+ e- collisions
		else                       pythia.init( 11, -11, 200.);         //e+ e- collisions
	}
	else        //default config
	{
		pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
		pythia.readString("23:onMode = off");
		pythia.readString("23:onIfAny = 13");
		pythia.init( 11, -11, 91.17);                           //e+ e- collisions
	}

	MC_Initialize();

	Photos::initialize();
	Photos::setDoubleBrem(false);
	Photos::setExponentiation(false);
	Photos::setMeCorrectionWtForZ(zNLO);
	// Zee and ttbar require higher maxWtInterference
	if(zeeDecays || topDecays) Photos::maxWtInterference(2.0);

	Photos::setInfraredCutOff(0.001);//91.187);
	Photos::maxWtInterference(2.0);
	if( zNLO) Photos::maxWtInterference(2.0);
	Log::SummaryAtExit();
	cout.setf(ios::fixed);

	// Begin event loop
	int fatal=0;
	for(unsigned long iEvent = 0; iEvent < NumberOfEvents; ++iEvent)
	{
		if(iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<"\t("<<iEvent*(100./NumberOfEvents)<<"%)"<<endl;
		if(!pythia.next())
		{
			--iEvent;
			++fatal;
			continue;
		}

		HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();
		ToHepMC.fill_next_event(event, HepMCEvt);
		//HepMCEvt->print();

		// If event contains gamma - abort
		for(HepMC::GenEvent::particle_const_iterator i = HepMCEvt->particles_begin();i!=HepMCEvt->particles_end();i++)
			if((*i)->pdg_id()==22) { HepMCEvt->print(); cout<<" Gamma in the event! Abort."<<endl; exit(-1); }

		//Log::LogPhlupa(1,3);

		// Call photos
		PhotosHepMCEvent evt(HepMCEvt);
		evt.process();

		//HepMCEvt->print();

		// We mess with the event so MC-Tester can work on it as in LC analysis case
		if(topDecays) fixForMctester(HepMCEvt);
		//HepMCEvt->print();

		// Call MC-Tester
		HepMCEvent temp_event(*HepMCEvt,false);
		MC_Analyze(&temp_event);

		//if(iEvent>=NumberOfEvents-5) HepMCEvt->print();

		// Clean up
		delete HepMCEvt;
	}
	pythia.statistics();
	MC_Finalize();
	cout<<" Events aborted by PYTHIA: "<<fatal<<endl;
}
