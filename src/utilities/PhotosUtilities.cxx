#include "PhotosUtilities.h"
#include <cstdlib>
#include <cstdio>

namespace Photospp
{

namespace PhotosUtilities
{

//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays BOost routine '3'
//
//    Purpose:  Boost  vector PVEC  along z-axis where ANGLE = EXP(ETA),
//              ETA is the hyperbolic velocity.
//
//    Input Parameters:  ANGLE, PVEC
//
//    Output Parameter:  PVEC
//
//    Author(s):  S. Jadach                       Created at:  01/01/89
//                B. van Eijk                     Last Update: 12/06/13
//
//----------------------------------------------------------------------
void PHOBO3(double ANGLE,double PVEC[4]){
  int j=1;  // convention of indices of Riemann space must be preserved.
  double QPL,QMI;
  QPL=(PVEC[4-j]+PVEC[3-j])*ANGLE;
  QMI=(PVEC[4-j]-PVEC[3-j])/ANGLE;
  PVEC[3-j]=(QPL-QMI)/2.0;
  PVEC[4-j]=(QPL+QMI)/2.0;
}

} // namespace PhotosUtilities
	
} // namespace Photospp

