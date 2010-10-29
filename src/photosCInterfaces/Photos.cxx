#include <stdarg.h>
#include <iostream>
#include <vector>

#include "PhotosEvent.h"
#include "Photos.h"
#include "Log.h"
using std::vector;
using std::cout;
using std::endl;
using std::ios_base;
typedef Photos::Log Log;

Photos Photos::_instance;

vector<vector<int>* >    *Photos::supBremList    = 0;
vector<vector<int>* >    *Photos::forceBremList  = 0;
bool Photos::isSuppressed=false;
bool Photos::massFrom4Vector=true;

Photos::Photos()
{
	setAlphaQED           (0.00729735039);
	setInfraredCutOff     (0.01);
	setInterference       (true);
	setDoubleBrem         (true);
	setQuatroBrem         (false);
	setExponentiation     (true);
	setTopProcessRadiation(true);
	setCorrectionWtForW   (true);
}

void Photos::initialize()
{
// Should return if already initialized?

/*******************************************************************************
  All the following parameter setters can be called after PHOINI.
  Initialization of kinematic correction against rounding errors.
  The set values will be used later if called with zero.
  Default parameter is 1 (no correction) optionally 2, 3, 4
  In case of exponentiation new version 5 is needed in most cases.
  Definition given here will be thus overwritten in such a case
  below in routine PHOCIN
*******************************************************************************/
	if(!phokey_.iexp) initializeKinematicCorrections(1);

	int buf=1;
	iphqrk_(&buf); // Blocks emission from quarks if buf=1 (default); enables if buf=2
	               // Physical treatment will be 3, option 2 is not realistic and for tests only
	buf=2;
	iphekl_(&buf); // Blocks emission in  pi0 to gamma e+ e- if parameter is >1 (enables otherwise)

// Initialise status counter for warning messages
	for(int i=0;i<10;i++) phosta_.status[i]=0;

	pholun_.phlun=6; // Logical output unit for printing error messages

// Define pi and 2*pi
	phpico_.pi   =3.14159265358979324;
	phpico_.twopi=6.28318530717958648;

// Further initialization done automatically
// see places with - VARIANT A - VARIANT B - all over to switch between options

//----------- SLOWER VARIANT A, but stable ------------------
//--- it is limiting choice for small XPHCUT in fixed orer
//--- modes of operation

// Best choice is if FINT=2**N where N+1 is maximal number
// of charged daughters
// see report on overweihted events
	if(phokey_.interf) maxWtInterference(2.0);
	else               maxWtInterference(1.0);

/*
----------- FASTER VARIANT B  ------------------
-- it is good for tests of fixed order and small XPHCUT
-- but is less promising for more complex cases of interference
-- sometimes fails because of that

	if(phokey_.interf) maxWtInterference(1.8);
	else               maxWtInterference(0.0);
------------END VARIANTS A B -----------------------
*/


//------------------------------------------------------------------------------
// Print PHOTOS header
//------------------------------------------------------------------------------
	int                coutPrec = cout.precision(2);
	ios_base::fmtflags flags    = cout.setf(ios_base::fixed);
	cout<<endl;
	cout<<"********************************************************************************"<<endl<<endl;

	cout<<"                            ========================="<<endl;
	cout<<"                              PHOTOS, Version:  "<<VER_MAJOR<<"."<<VER_MINOR<<endl;
	cout<<"                              Released at:  "<<DAT_DAY<<"/"<<DAT_MONTH<<"/"<<DAT_YEAR<<endl;
	cout<<"                            ========================="<<endl<<endl;

	cout<<"                     Photos QED corrections in Particle Decays"<<endl<<endl;

	cout<<"           Monte Carlo Program - by E. Barberio, B. van Eijk and Z. Was"<<endl;
	cout<<"           From version 2.09   - by P. Golonka and Z. Was"<<endl;
	cout<<"           From version 3.00   - by N. Davidson, T. Przedzinski and Z. Was"<<endl;

	cout<<"********************************************************************************"<<endl<<endl;

	cout<<"                    Internal input parameters: "<<endl<<endl;
	cout<<"                    INTERF= "<<phokey_.interf<<" ISEC= " <<phokey_.isec <<" ITRE= "<<phokey_.itre
	                     <<" IEXP= "  <<phokey_.iexp  <<" IFTOP= "<<phokey_.iftop<<" IFW= " <<phokey_.ifw <<endl;
	cout<<"                    ALPHA_QED= "<<phocop_.alpha<<" XPHCUT= "<<phocop_.xphcut<<endl<<endl;

	if(phokey_.interf) cout<<"                    Option with interference is active"<<endl;
	if(phokey_.isec)   cout<<"                    Option with double photons is active"<<endl;
	if(phokey_.itre)   cout<<"                    Option with triple/quatric photons is active"<<endl;
	if(phokey_.iexp)   cout<<"                    Option with exponentiation is active EPSEXP="<<phokey_.expeps<<endl;
	if(phokey_.iftop)  cout<<"                    Emision in t tbar production is active"<<endl;
	if(phokey_.ifw)    cout<<"                    Correction wt in decay of W is active"<<endl;

	cout<<endl<<"          WARNING:  /HEPEVT/ is not anymore used."<<endl<<endl;
/*
	cout<<endl<<"            WARNING (1): /HEPEVT/ is not anymore the standard common block"<<endl<<endl;

	cout<<"            PHOTOS expects /HEPEVT/ to have REAL*8 variables. To change to"<<endl;
	cout<<"            REAL*4 modify its declaration in subr. PHOTOS_GET PHOTOS_SET:"<<endl;
	cout<<"                 REAL*8  d_h_phep,  d_h_vhep"<<endl;
	cout<<"            WARNING (2): check dims. of /hepevt/ /phoqed/ /ph_hepevt/."<<endl;
	cout<<"            HERE:                     d_h_nmxhep=10000 and  NMXHEP=10000"<<endl<<endl;
*/
	cout<<"********************************************************************************"<<endl;
	// Revert output stream flags and precision
	cout.precision(coutPrec);
	cout.flags    (flags);

// Initialize Marsaglia and Zaman random number generator
	PhotosRandom::initialize();
}

void Photos::processParticle(PhotosParticle *p)
{
	PhotosBranch b(p);
	if(!b.getSuppressionStatus()) b.process();
}

void Photos::processBranch(PhotosParticle *p)
{
	vector<PhotosParticle *> particles = p->getDecayTree();
	vector<PhotosBranch *>   branches = PhotosBranch::createBranches(particles);
	for(int i=0;i<(int)branches.size();i++) branches.at(i)->process();
}

void Photos::suppressBremForDecay(int count, int motherID, ... )
{
	va_list arg;
	va_start(arg, motherID);
	vector<int> *v = new vector<int>();
	v->push_back(motherID);
	for(int i = 0;i<count;i++)
	{
		v->push_back(va_arg(arg,int));
	}
	va_end(arg);
	v->push_back(0);
	if(!supBremList) supBremList = new vector< vector<int>* >();
	supBremList->push_back(v);
}

void Photos::suppressBremForBranch(int count, int motherID, ... )
{
	va_list arg;
	va_start(arg, motherID);
	vector<int> *v = new vector<int>();
	v->push_back(motherID);
	for(int i = 0;i<count;i++)
	{
		v->push_back(va_arg(arg,int));
	}
	va_end(arg);
	v->push_back(1);
	if(!supBremList) supBremList = new vector< vector<int>* >();
	supBremList->push_back(v);
}

void Photos::forceBremForDecay(int count, int motherID, ... )
{
	va_list arg;
	va_start(arg, motherID);
	vector<int> *v = new vector<int>();
	v->push_back(motherID);
	for(int i = 0;i<count;i++)
	{
		v->push_back(va_arg(arg,int));
	}
	va_end(arg);
	v->push_back(0);
	if(!forceBremList) forceBremList = new vector< vector<int>* >();
	forceBremList->push_back(v);
}

void Photos::forceBremForBranch(int count, int motherID, ... )
{
	va_list arg;
	va_start(arg, motherID);
	vector<int> *v = new vector<int>();
	v->push_back(motherID);
	for(int i = 0;i<count;i++)
	{
		v->push_back(va_arg(arg,int));
	}
	va_end(arg);
	v->push_back(1);
	if(!forceBremList) forceBremList = new vector< vector<int>* >();
	forceBremList->push_back(v);
}

void Photos::setMeCorrectionWtForW(bool corr)
{
	Log::Info()<<"Photos::setMeCorrectionWtForW: option not implemented in PHOTOS 3.0"<<endl;
}

void Photos::setMeCorrectionWtForZ(bool corr)
{
	Log::Info()<<"Photos::setMeCorrectionWtForZ: option not implemented in PHOTOS 3.0"<<endl;
}

