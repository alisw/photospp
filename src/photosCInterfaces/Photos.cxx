#include <stdarg.h>
#include <vector>

#include "PhotosEvent.h"
#include "Photos.h"
#include "Log.h"
using std::vector;
typedef Photos::Log Log;

bool Photos::isSuppressed=false;

vector<vector<int>* >    *Photos::supBremList    = 0;
vector<vector<int>* >    *Photos::forceBremList  = 0;

void Photos::initialize()
{
	phoini_();

	setAlphaQED(0.00729735039);
	setExponentiation(true);
	setInterference(true);
	setTopProcessRadiation(true);
	//	maxWtInterference(3.0); 
	//For Zee, where weight may be above unity
        // It needs investigation.
	//Log::LogPhlupa(0,1000);
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
