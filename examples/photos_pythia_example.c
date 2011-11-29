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
#include "Photos/Photos.h"
#include "Photos/PhotosHepMCEvent.h"
#include "Photos/Log.h"

using namespace std;
using namespace Pythia8;

bool ShowersOn=true;
unsigned long NumberOfEvents = 10000;

int main(int argc,char **argv)
{
	// Initialization of pythia
	HepMC::I_Pythia8 ToHepMC;
	Pythia pythia;
	Event& event = pythia.event;
	//pythia.settings.listAll();

	pythia.readString("PartonLevel:ISR = on");
	pythia.readString("PartonLevel:FSR = off");

	pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
	pythia.readString("23:onMode = off");
	pythia.readString("23:onIfAny = 13");
	pythia.init( 11, -11, 91.187);                           //e+ e- collisions

	MC_Initialize();

	Photos::initialize();
	//Photos::setDoubleBrem(false);
	//Photos::setExponentiation(false);

	Photos::setInfraredCutOff(0.01/91.187); // 10MeV for scale to M_Z=91.187
	Photos::maxWtInterference(3.0);

	Photos::iniInfo();
	Log::SummaryAtExit();
	cout.setf(ios::fixed);

	// Begin event loop
	for(unsigned long iEvent = 0; iEvent < NumberOfEvents; ++iEvent)
	{
		if(iEvent%1000==0) Log::Info()<<"Event: "<<iEvent<<"\t("<<iEvent*(100./NumberOfEvents)<<"%)"<<endl;
		if (!pythia.next()) continue;

		// Convert event record to HepMC
		HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();
		ToHepMC.fill_next_event(event, HepMCEvt);
		//HepMCEvt->print();

		//Log::LogPhlupa(1,3);

		// Run PHOTOS on the event
		PhotosHepMCEvent evt(HepMCEvt);
		evt.process();

		//HepMCEvt->print();

		// Run MC-TESTER on the event
		HepMCEvent temp_event(*HepMCEvt,false);
		MC_Analyze(&temp_event);

		// Print out last 5 events
		if(iEvent>=NumberOfEvents-5) HepMCEvt->print();

		// Clean up
		delete HepMCEvt;
	}
	pythia.statistics();
	MC_Finalize();
}
