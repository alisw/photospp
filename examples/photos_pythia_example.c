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

using namespace std;
using namespace Pythia8; 

bool ShowersOn=true;
int NumberOfEvents = 10000;

int main(int argc,char **argv){
  HepMC::I_Pythia8 ToHepMC;
  // Initialisation of pythia
  Pythia pythia;
  Event& event = pythia.event;

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
    pythia.init( 11, -11, 200.);  //e+ e-
    //    pythia.init( -2212, -2212, 14000.0); //proton proton collisions
  }
  else        //default config
  {
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    pythia.readString("23:onMode = off"); 
    pythia.readString("23:onIfAny = 11");
    pythia.init( 11, -11, 100.);          //electron positron collisions
  }

  MC_Initialize();

  Photos::initialize();
  Photos::setInfraredCutOff(0.001/200);//91.187);

  Log::SummaryAtExit();

  // Begin event loop. Generate event.
  for (int iEvent = 0; iEvent < NumberOfEvents; ++iEvent) {
    if(iEvent%1000==0) cout<<iEvent<<endl;
    if (!pythia.next()) continue;

    HepMC::GenEvent * HepMCEvt = new HepMC::GenEvent();
    HepMCEvt->use_units(HepMC::Units::GEV,HepMC::Units::MM);
    ToHepMC.fill_next_event(event, HepMCEvt);

    //call photos
    //HepMCEvt->print();
    //Log::LogPhlupa(2,4);
	PhotosHepMCEvent evt(HepMCEvt);
	evt.process();
    //HepMCEvt->print();

    //call mc-tester
    HepMCEvent temp_event(*HepMCEvt,false);
    MC_Analyze(&temp_event);

    if(iEvent>=NumberOfEvents-5)
      { 
	HepMCEvt->print();
      }

    //clean up
    delete HepMCEvt;
  }

  pythia.statistics();
  MC_Finalize();


}

