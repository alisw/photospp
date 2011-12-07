#ifndef _Photos_h_included_
#define _Photos_h_included_

/**
 * @class Photos
 *
 * @brief Controls the configuration and initialization of Photos.
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
	static const int VER_MAJOR=3, VER_MINOR=3;
	static const int DAT_DAY  =7,DAT_MONTH=12,DAT_YEAR=11;
public:
	/** Logging and memory leak tracking class */
	class Log;

	/** Initalize Photos with the parameters previously set via the
	   setter methods */
	static void initialize();

	/** Prints info on  Photos initialization (reinitialization)
	   status */
	static void iniInfo();

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

	/* Key for partial effects of  matrix element (in leptonic W decays) */
	static void setCorrectionWtForW(bool corr) { phokey_.ifw=(int)corr; }

	/** Set exponentiation mode */
	static void setExponentiation(bool expo);

	/** Switch for complete effects of matrix element (in  scalar  to 2 scalars decays) */
	static void setMeCorrectionWtForScalar(bool corr);

	/** Switch for complete effects of matrix element (in leptonic W decays) */
	static void setMeCorrectionWtForW(bool corr);

	/** Switch for complete effects of matrix element (in leptonic Z decays) */
	static void setMeCorrectionWtForZ(bool corr);

	/** Set photon emission in top pair production in quark (gluon) pair annihilation */
	static void setTopProcessRadiation(bool top)         { phokey_.iftop=(int)top; }

	/** Initialize kinematic corrections */
	static void initializeKinematicCorrections(int flag) { phcork_(&flag); }

	/** Force mass value to be sqrt(e^2-p^2) for all particle momenta
	    taken from event record. May be important for numerical stability.
		May lead to faulty results due to rounding errors for
		hiper-relativistic electron, for example. */
	static void forceMassFrom4Vector(bool flag) { massFrom4Vector=flag; }
	
	/** set energy momentum conservation threshold */
	static void setMomentumConservationThreshold(double threshold){momentum_conservation_threshold=threshold; }

public:
	/** Is in suppressed mode */
	static bool isSuppressed;

	/** Is mass from 4-vector or from event record */
	static bool massFrom4Vector;
	
	/** List of suppressed decays */
	static vector<vector<int>* >    *supBremList;

	/** List of forced decays */
	static vector<vector<int>* >    *forceBremList;

 	/** Threshold for momentum conservation check */
	static double momentum_conservation_threshold;

	/** Flag for complete effects of matrix element (in scalars decays) */
	static bool meCorrectionWtForScalar;

	/** Flag for complete effects of matrix element (in leptonic Z decays) */
	static bool meCorrectionWtForZ;
	
	/** Flag for complete effects of matrix element (in leptonic W decays) */
	static bool meCorrectionWtForW;
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

