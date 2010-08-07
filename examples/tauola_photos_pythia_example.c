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
#include "PhotosHepMCParticle.h"
#include "PhotosHepMCEvent.h"
#include "Log.h"
typedef Photos::Log Log; //We're using Photos version of Log class

//TAUOLA header files
#include "Tauola.h"
#include "TauolaHepMCEvent.h"

using namespace std;
using namespace Pythia8;

bool ShowersOn=true;
int NumberOfEvents = 10000;

int main(int argc,char **argv)
{
	HepMC::I_Pythia8 ToHepMC;

	// Initialisation of pythia
	Pythia pythia;
	Event& event = pythia.event;

	int tauolaDecayMode=0;
	if(argc>4) tauolaDecayMode=atoi(argv[4]);
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
	pythia.readString("PartonLevel:ISR = off");
	pythia.readString("PartonLevel:FSR = off");
	if(argc>1)  //pre-set configs
	{
		pythia.readFile(argv[1]);
		if(tauolaDecayMode==3) pythia.init( 11, -11, 91.17);
		else pythia.init( 11, -11, 500.);
		//3 for Ztautau, 4 for Htautau
		Tauola::setSameParticleDecayMode(tauolaDecayMode);
		Tauola::setOppositeParticleDecayMode(tauolaDecayMode);
	}
	else        //default config
	{
		/** pythia.readString("HiggsSM:ffbar2H = on");
		pythia.readString("25:onMode = off");
		pythia.readString("25:onIfAny = 15");
		pythia.readString("25:m0 = 120"); */

		/** pythia.readString("WeakDoubleBoson:ffbar2WW = on");
		pythia.readString("24:onMode = off");
		pythia.readString("24:onIfAny = 15");**/

		pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
		pythia.readString("23:onMode = off");
		pythia.readString("23:onIfAny = 15");
		pythia.particleData.readString("15:mayDecay = off"); //<- uncomment for pythia+tauola
		//pythia.init( -2212, -2212, 14000.0);     //proton proton collisions
		pythia.init( 11, -11, 500.);             //electron positron collisions
	}
	Tauola::initialise();

	Photos::initialize();
	Photos::setInfraredCutOff(0.01/200);//91.187);

	Log::SummaryAtExit();
	//Log::LogInfo(false) //To turn printing of last five events and pythia statistics off

	// Example setup - suppress processing of whole Z0 decay,
	// leaving only the Z0 -> tau+ tau- decay and whole branch starting
	// from tau- to be processed
	//Photos::suppressBremForBranch(0,23);
	//Photos::forceBremForDecay (2,23,15,-15);
	//Photos::forceBremForBranch(0,15);

	MC_Initialize();
	// Begin event loop. Generate event.
	for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent)
	{
		if(iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<"\t("<<(iEvent*100)/NumberOfEvents<<"%)"<<endl;
		if(!pythia.next()) continue;

		HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();
		HepMCEvt->use_units(HepMC::Units::GEV,HepMC::Units::MM);
		ToHepMC.fill_next_event(event, HepMCEvt);

		TauolaHepMCEvent * t_event = new TauolaHepMCEvent(HepMCEvt);
		t_event->decayTaus();

		//HepMCEvt->print();
		//Log::LogPhlupa(2,4);
		PhotosHepMCEvent evt(HepMCEvt);
		evt.process();
		//HepMCEvt->print();

		HepMCEvent temp_event(*HepMCEvt,false);
		MC_Analyze(&temp_event);
		if(iEvent>=NumberOfEvents-5)
		{
			Log::RedirectOutput(Log::Info());
			//pythia.event.list();
			HepMCEvt->print();
			Log::RevertOutput();
		}
		//clean up
		delete HepMCEvt;
		delete t_event; //<- uncomment for pythia+tauola
	}
	Log::RedirectOutput(Log::Info());
	pythia.statistics();
	Log::RevertOutput();
	MC_Finalize();
}
