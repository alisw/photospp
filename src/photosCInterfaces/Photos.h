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
#include "PhotosRandom.h"
#include "f_Init.h"
using std::vector;

class PhotosParticle;

class Photos
{
public:
	static const int VER_MAJOR=3, VER_MINOR=0;
	static const int DAT_DAY  =12,DAT_MONTH=8,DAT_YEAR=10;
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

	/** Suppress all processing. Only forced decays will be processed. */
	static void suppressAll()                      { isSuppressed=true; }

	/** Force processing of a single decay */
	static void forceBremForDecay (int count, int motherID, ... );

	/** Force processing of a whole decay branch */
	static void forceBremForBranch(int count, int motherID, ... );

public:
	/** Seed for RANMAR used by fortran part of the Photos */
	static void setSeed(int iseed1, int iseed2)    { PhotosRandom::setSeed(iseed1,iseed2); }

	/** Maximum interference weight */
	static void maxWtInterference(double interference) { phokey_.fint=interference; }

	/** Minimal energy (in units of decaying particle mass) for photons to be explicitly generated */
	static void setInfraredCutOff(double cut_off)  { phocop_.xphcut=cut_off; }

	/** Coupling constant alpha QED */
	static void setAlphaQED(double alpha)          { phocop_.alpha=alpha; }

	/** Key for interference, matrix element weight */
	static void setInterference(bool interference) { phokey_.interf=(int)interference; }

	/** Set double bremsstrahlung generation */
	static void setDoubleBrem(bool doub)           { phokey_.isec=(int)doub; }

	/** Set bremsstrahlung generation up to multiplicity of 4 */
	static void setQuatroBrem(bool quatroBrem)     { phokey_.itre=(int)quatroBrem; }

	/* Key for effects of initial state charge (in leptonic W decays) */
	static void setCorrectionWtForW(bool corr) { phokey_.ifw=(int)corr; }

	/** Set exponentiation mode */
	static void setExponentiation(bool expo)
	{
		phokey_.iexp = (int) expo;
		if(expo)
		{
			setDoubleBrem(false);
			setQuatroBrem(false);
			setInfraredCutOff(0.0000001);
			initializeKinematicCorrections(5);
			phokey_.expeps=0.0001;
		}
	};

	/** Switch for complete effects of matrix element (in leptonic W decays) */
	static void setMeCorrectionWtForW(bool corr);

	/** Switch for complete effects of matrix element (in leptonic Z decays) */
	static void setMeCorrectionWtForZ(bool corr);

	/** Set photon emission in top pair production in quark (gluon) pair annihilation */
	static void setTopProcessRadiation(bool top)         { phokey_.iftop=(int)top; }

	/** Initialize kinematic corrections */
	static void initializeKinematicCorrections(int flag) { phcork_(&flag); }

public:
	/** Is in suppressed mode */
	static bool isSuppressed;

	/** List of suppressed decays */
	static vector<vector<int>* >    *supBremList;

	/** List of forced decays */
	static vector<vector<int>* >    *forceBremList;

public:
	/** Get instance of Photos */
	Photos& getInstance() { return _instance; }
private:
	/* Singleton: only one instance allowed.
	   Constructor sets default values of PHOTOS parameters */
	 Photos();
	~Photos() {}
	Photos(const Photos&);
	Photos& operator=(const Photos&);
	static Photos _instance;
};

#endif

