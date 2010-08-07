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
int NumberOfEvents = 10000;


// Finds X Y -> 6 -6 decay and converts it to 100 -> Y 6 -6, where 100 = X + 2*Y
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
		HepMC::FourVector fXY(fX.px()+fY.px()+fY.px(),fX.py()+fY.py()+fY.py(),fX.pz()+fY.pz()+fY.pz(),fX.e()+fY.e()+fY.e());
		X->set_momentum(fXY);
		// Unique ID for MC-Tester to analyze
		X->set_pdg_id(100);

		// Set 2nd mother as decayed and move it to end vertex
		Y->set_status(1);
		(* Y->production_vertex()->particles_in_const_begin())->set_status(1);
		pt->production_vertex()->remove_particle(Y);
		pt->production_vertex()->add_particle_out(Y);
		//evt->remove_vertex( Y->production_vertex() );
		return;
	}
}

int main(int argc,char **argv)
{
	// Initialisation of pythia
	HepMC::I_Pythia8 ToHepMC;
	Pythia pythia;
	Event& event = pythia.event;
	//pythia.settings.listAll();
	bool topDecays   =false;
	bool ppGeneration=false;
	if(argc>4)
	{
		// Advanced options
		topDecays    = (atoi(argv[4])==1);
		ppGeneration = (atoi(argv[4])==2);
	}
	if(argc>3) NumberOfEvents=atoi(argv[3]);
	if(argc>2) ShowersOn=atoi(argv[2]);
	if(!ShowersOn)
	{
		//pythia.readString("HadronLevel:all = off");
		pythia.readString("HadronLevel:Hadronize = off");
		pythia.readString("SpaceShower:QEDshower = off");
		pythia.readString("SpaceShower:QEDshowerByL = off");
		pythia.readString("SpaceShower:QEDshowerByQ = off");
	}
	pythia.readString("PartonLevel:ISR = on");
	pythia.readString("PartonLevel:FSR = off");
	if(argc>1)  //pre-set configs
	{
		pythia.readFile(argv[1]);
		if(ppGeneration)   pythia.init( -2212, -2212, 14000.0); //p  p  collisions
		else if(topDecays) pythia.init( -2212, -2212, 14000.0); //p  p  collisions
		else               pythia.init( 11, -11, 91.17);         //e+ e- collisions
	}
	else        //default config
	{
		pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
		pythia.readString("23:onMode = off");
		pythia.readString("23:onIfAny = 11");
		pythia.init( 11, -11, 91.17);                           //e+ e- collisions
	}

	MC_Initialize();

	Photos::initialize();
	Photos::setInfraredCutOff(0.001/200);//91.187);
	Log::SummaryAtExit();

	// Begin event loop
	for(int iEvent = 0; iEvent < NumberOfEvents; ++iEvent)
	{
		if(iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<"\t("<<(iEvent*100)/NumberOfEvents<<"%)"<<endl;
		if (!pythia.next()) continue;

		HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();
		HepMCEvt->use_units(HepMC::Units::GEV,HepMC::Units::MM);
		ToHepMC.fill_next_event(event, HepMCEvt);
		//HepMCEvt->print();

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
}
