#ifndef _PhotosEvent_h_included_
#define _PhotosEvent_h_included_

/**
 * @class PhotosEvent
 *
 * @brief Abstract base class for containing the event information.
 *
 * PhotosEvent contains virtual methods, which need to be implemented
 * by the appropriate interface class to the event record. An object of
 * PhotosEvent type should be created by the user and processed
 * via the process() method.
 *
 * This class is responsible for finding branching points, and invoking
 * photos on them.
 *
 * @author Nadia Davidson
 * @date 16 June 2008
 */
#include <vector>
#include "PhotosParticle.h"
using std::vector;

class PhotosEvent
{
public:
	virtual std::vector<PhotosParticle*> getBranchPoints() = 0;
	virtual void print() = 0;

	void process();
	~PhotosEvent(); 
private:
	vector<PhotosParticle *> filterBranchPoints();

	/** filter for branching points **/
	bool passBranchPointFilter(PhotosParticle *particle);

	/** Check if the branching point should be skipped by PHOTOS. */
	bool passSuppressionFilter(PhotosParticle *particle);

	/** Check if the particles' mother is on the Photos::supParticles list.
	  If it is, it will also be skipped. */
	bool passSuppressConsecutive(PhotosParticle *particle);

	/** branch points which should be given to PHOTOS */
	vector<PhotosParticle*> m_branch_points;
};

#endif  
