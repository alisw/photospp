/**
 * Example of photos usage.
 * Events are loaded from pre-generated set featuring Z0 -> tau+ tau- decays
 * and processed by photos.
 *
 * @author Tomasz Przedzinski
 * @date 26 January 2020
 */

// HepMC3 header files
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/Print.h"

// PHOTOS header files
#include "Photos/Photos.h"
#include "Photos/PhotosHepMC3Event.h"
#include "Photos/Log.h"

using namespace std;
using namespace Photospp;

int EventsToCheck=20;

// elementary test of HepMC typically executed before
// detector simulation based on http://home.fnal.gov/~mrenna/HCPSS/HCPSShepmc.html
// similar test was performed in Fortran
// we perform it before and after Photos (for the first several events)
void checkMomentumConservationInEvent(HepMC3::GenEvent &evt)
{
	double px = 0.0;
	double py = 0.0;
	double pz = 0.0;
	double e = 0.0;

	for (auto p : evt.particles())
	{
		if (p->status() == 1)
		{
			HepMC3::FourVector m = p->momentum();
			px += m.px();
			py += m.py();
			pz += m.pz();
			e  += m.e();
		}
	}
	cout.precision(6);
	cout.setf(ios_base::scientific, ios_base::floatfield);
	cout << endl << "Vector Sum: " << px << " " << py << " " << pz << " " << e << endl;
}

int main()
{
	HepMC3::ReaderAscii input_file("photos_hepmc3_standalone_example.dat");

	Photos::initialize();
	Photos::setInfraredCutOff(0.001/200);

	int photonAdded = 0;
	int twoAdded = 0;
	int moreAdded = 0;
	int evtCount = 0;

	// Begin event loop. Generate event.
	while (!input_file.failed()) {

        HepMC3::GenEvent evt(Units::GEV, Units::MM);

        // Read event from input file
        input_file.read_event(evt);

        // If reading failed - exit loop
        if (input_file.failed()) {
			break;
		}

		evtCount++;

		int buf = -evt.particles().size();

		//cout << "BEFORE:"<<endl;
	    //Print::listing(evt);

		if (evtCount < EventsToCheck)
		{
			cout << endl;
			cout << "Momentum conservation chceck BEFORE/AFTER Photos" << endl;
			checkMomentumConservationInEvent(evt);
		}

		// Process by photos
		PhotosHepMC3Event photosEvent(&evt);
		photosEvent.process();

		if (evtCount < EventsToCheck)
		{
			checkMomentumConservationInEvent(evt);
		}

		buf += evt.particles().size();

		if (buf == 1) {
			photonAdded++;
		}     
		else if (buf == 2) {
			twoAdded++;
		} 
		else if (buf > 2) {
			moreAdded++;
		}

		//cout << "AFTER:" << endl;
	    //Print::listing(evt);
	}

    input_file.close();

	// Print results
	cout.precision(2);
	cout.setf(ios_base::fixed, ios_base::floatfield);
	cout << endl;

	if (evtCount == 0)
	{
		cout<<"Something went wrong with the input file: photos_standalone_example.dat"<<endl;
		cout<<"No events were processed."<<endl<<endl;
		return 0;
	}

	cout << "Summary (whole event processing):" << endl;
	cout << evtCount    << "\tevents processed" << endl;
	cout << photonAdded << "\ttimes one photon added to the event           \t(" << (photonAdded*100./evtCount) << "%)" << endl;
	cout << twoAdded    << "\ttimes two photons added to the event          \t(" << (twoAdded*100./evtCount) << "%)" << endl;
	cout << moreAdded   << "\ttimes more than two photons added to the event\t(" << (moreAdded*100./evtCount) << "%)" << endl << endl;
	cout << "(Contrary to results from MC-Tester, these values are technical and infrared unstable)" << endl << endl;
}
