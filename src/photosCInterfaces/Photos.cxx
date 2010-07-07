#include <stdarg.h>
#include <iostream>
#include <vector>

#include "Photos.h"
#include "PhotosHepMCEvent.h"
#include "PH_HEPEVT_Interface.h"

#include "Log.h"
using std::vector;

vector< vector<int>* > *Photos::supBremList = 0;
vector< HepMC::GenParticle* > *Photos::supParticles = 0;


void Photos::initialize()
{
	phoini_();
	
	setAlphaQED(0.00729735039);
	setExponentiation(true);
	setInterference(true);
	setTopProcessRadiation(true);
	// debug mode on: printing will be activated for adresses in range 0 to 1000 
	//  phlupy_.ipoinm=0;
	//  phlupy_.ipoin=1000;
}

void Photos::process(HepMC::GenEvent * event)
{
	//event->print();
	int index = 2;
	
	PhotosHepMCEvent ph_evt(event);
	
	std::vector<PhotosParticle *> branch_points;
	branch_points = ph_evt.getBranchPoints();
	
	for(int i=0; i < (int) branch_points.size(); i++)
	{	
		PH_HEPEVT_Interface::set(branch_points.at(i));
		photos_make_(&index);
		PH_HEPEVT_Interface::get();
		branch_points.at(i)->checkMomentumConservation();
	}

	//  event->print();
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
