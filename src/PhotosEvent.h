#ifndef _PhotosEvent_h_included_
#define _PhotosEvent_h_included_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "PhotosParticle.h"

/**
 * @class PhotosEvent
 *
 * @brief Abstract base class for containing the event information.
 *
 * PhotosEvent contains virtual methods, which need to be implemented
 * by the appropriate interface class to the event record. Currently only
 * PhotosHepMCEvent does this. An object of PhotosEvent type should be
 * created by the user and can be decayed via the decayTaus() method.
 *
 * This class is responsible for finding taus, (or tau and 
 * it's neutrino) and creating PhotosParticlePairs out of them.
 *
 * @author Nadia Davidson
 * @date 16 June 2008
 */
class PhotosEvent{

 public:

   virtual std::vector<PhotosParticle*> getBranchPoints() = 0;

   virtual ~PhotosEvent(){};

 private:    

};

#endif  

