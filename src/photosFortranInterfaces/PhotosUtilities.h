#ifndef _PhotosUtilities_h_included_
#define _PhotosUtilities_h_included_

/**
 * @class PhotosUtilities
 *
 * @brief Support functions
 *
 * Functions for boosting, rotation, ...
 *
 * @author Tomasz Przedzinski, Zbigniew Was
 * @date 29 June 2013
 */

namespace Photospp
{

class PhotosUtilities
{
public:
	/** PHOton radiation in decays BOost routine '3' */
	static void PHOBO3(double ANGLE,double PVEC[4]);
};

} // namespace Photospp
#endif

