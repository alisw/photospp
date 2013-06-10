#include "Photos.h"
#include <cmath>
using namespace Photospp;

/** Logical function used deep inside algorithm to check if emitted
    particles are to emit. For mother it blocks the vertex, 
    but for daughters individually: bad sisters will not prevent electron to emit.
    top quark has further exception method. */
bool F(int m, int i)
{
  return Photos::IPHQRK_setQarknoEmission(0,i) && (i<= 41 || i>100)
     && i != 21 
     && i != 2101 && i !=3101 && i !=3201 
     && i != 1103 && i !=2103 && i !=2203 
     && i != 3103 && i !=3203 && i !=3303;
}

//----------------------------------------------------------------------
//
//    PHOEPS:   PHOeps vector product (normalized to unity)
//
//    Purpose:  calculates vector product, then normalizes its length.
//              used to generate orthogonal vectors, i.e. to
//              generate polarimetric vectors for photons.
//
//    Input Parameters:  VEC1,VEC2 - input 4-vectors
//                                          
//    Output Parameters: EPS - normalized 4-vector, orthogonal to
//                             VEC1 and VEC2
//
//    Author(s):  Z. Was, P.Golonka               Created at:  19/01/05
//                                                Last Update: 10/06/13
//
//----------------------------------------------------------------------

void PHOEPS(double vec1[4], double vec2[4], double eps[4]){
  double xn;
  int j=1;  // convention of indices of Riemann space must be preserved.

  eps[1-j]=vec1[2-j]*vec2[3-j] - vec1[3-j]*vec2[2-j];
  eps[2-j]=vec1[3-j]*vec2[1-j] - vec1[1-j]*vec2[3-j];      
  eps[3-j]=vec1[1-j]*vec2[2-j] - vec1[2-j]*vec2[1-j];
  eps[4-j]=0.0;
      
  xn=sqrt( eps[1-j]*eps[1-j] + eps[2-j]*eps[2-j] + eps[3-j]*eps[3-j]);
      
  eps[1-j]=eps[1-j]/xn;
  eps[2-j]=eps[2-j]/xn;
  eps[3-j]=eps[3-j]/xn;

}
