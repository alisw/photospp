#ifndef _PhotosBranch_h_included_
#define _PhotosBranch_h_included_

/**
 * @class PhotosBranch
 *
 * @brief Single branching point
 *
 * Contains information about daughters and mothers of a single branch.
 * Each branch will be converted to HEPEVT and processed by photos.
 *
 * @author Tomasz Przedzinski
 * @date 8 July 2010
 */

#include <vector>
#include "PhotosParticle.h"
using std::vector;

class PhotosBranch
{
public:
	PhotosBranch(PhotosParticle* p);
	vector<PhotosParticle *> getMothers()          { return mothers;   }
	vector<PhotosParticle *> getDaughters()        { return daughters; }
	PhotosParticle*          getDecayingParticle() { return particle;  }
	vector<PhotosParticle *> getParticles();
	bool isValid() { return valid; }
	bool checkMomentumConservation();
private:
	/** Filter for branches. Checks if branching is suppressed by PHOTOS. */
	bool passBranchPointFilter();
	bool valid;
	PhotosParticle          *particle;
	vector<PhotosParticle *> mothers;
	vector<PhotosParticle *> daughters;
};

#endif
