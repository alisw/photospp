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
#include "PhotosBranch.h"
#include "PhotosParticle.h"
using std::vector;

class PhotosEvent
{
public:
	virtual std::vector<PhotosParticle*> getParticleList() = 0;
	virtual void print() = 0;

	void process();
	~PhotosEvent();
private:
	/** Creates branches from particles list removing the list in the process */
	vector<PhotosBranch *> createBranches(vector<PhotosParticle *> particles);

	/** Filter suppressed and invalid particles. **/
	bool passParticleFilter(PhotosParticle *particle);

	/** Check if the particles' mother is on the Photos::supParticles list.
	  If it is, it will also be skipped. */
	bool passSuppressConsecutive(PhotosParticle *particle);

	/** branch points which should be given to PHOTOS */
	vector<PhotosBranch *> m_branch_points;
};

#endif
