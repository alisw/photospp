#ifndef _Photos_h_included_
#define _Photos_h_included_

/**
 * @class Photos
 *
 * @brief Controls the configuration and initialisation of Photos.
 *
 * This is the main configuration class for Photos C++ Interface.
 * It is also used for invoking methods for processing single particle or branch.
 *
 * @author Nadia Davidson
 * @date 16th June 2008
 */
#include <stdarg.h>
#include <vector>
#include "PhotosParticle.h"
#include "f_Init.h"
using std::vector;

class PhotosParticle;

class Photos
{
public:
	/** Logging and memory leak tracking class */
	class Log;

	/** Initalise Photos with the parameters previously set via the
	   setter methods */
	static void initialize();

	/** Process decay of single particle */
	static void processParticle(PhotosParticle *p);
	/** Process decay of whole decay branch starting from given particle */
	static void processBranch(PhotosParticle *p);

	/** Suppress processing of a single decay */
	static void suppressBremForDecay (int count, int motherID, ... );
	/** Suppress processing of whole decay branch */
	static void suppressBremForBranch(int count, int motherID, ... );

	/** Suppress all processing.
	    Only forced decays will be processed. */
	static void suppressAll()                      { isSuppressed=true; }

	/** Force processing of a single decay */
	static void forceBremForDecay (int count, int motherID, ... );
	/** Force processing of a whole decay branch */
	static void forceBremForBranch(int count, int motherID, ... );

	//
	//ROUTINES BELLOW NEED COMMENTS FOR DOXYGEN!!!
	//

	//Wrappers for the PHOTOS configuration variables
	static void setInfraredCutOff(double cut_off)  { phocop_.xphcut=cut_off; }
	static void setAlphaQED(double alpha)          { phocop_.alpha=alpha; }
	static void setInterference(bool interference) { phokey_.interf=(int)interference; }
	static void setDoubleBrem(bool doub)           { phokey_.isec=(int)doub; }
	static void setHigherBrem(bool higherBrem)     { phokey_.itre=(int)higherBrem; }
	static void setExponentiation(bool expo)
	{
		phokey_.iexp = (int) expo;
		if(expo)
		{
			setDoubleBrem(false);
			setHigherBrem(false);
			setInfraredCutOff(0.0000001);
			initializeKinematicCorrections(5);
			phokey_.expeps=0.0001;
		}
	};

	static void setTopProcessRadiation(bool top)         { phokey_.iftop=(int)top; }
	static void initializeKinematicCorrections(int flag) { phcork_(&flag); }

	/** Is in suppressed mode */
	static bool isSuppressed;

	/** List of suppressed decays */
	static vector<vector<int>* >    *supBremList;

	/** List of forced decays */
	static vector<vector<int>* >    *forceBremList;
};

#endif

