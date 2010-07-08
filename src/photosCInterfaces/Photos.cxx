#include <stdarg.h>
#include <vector>

#include "Photos.h"
#include "Log.h"
using std::vector;

vector<vector<int>* >    *Photos::supBremList  = 0;
vector<PhotosParticle* > *Photos::supParticles = 0;

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
