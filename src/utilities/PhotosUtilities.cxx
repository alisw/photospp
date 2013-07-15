#include "PhotosUtilities.h"
#include <cstdlib>
#include <cstdio>

namespace Photospp
{

namespace PhotosUtilities
{

//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays calculation of TRIangle fie
//
//    Purpose:  Calculation of triangle function for phase space.
//
//    Input Parameters:  A, B, C (Virtual) particle masses.
//
//    Output Parameter:  Function value =
//                       SQRT(LAMBDA(A**2,B**2,C**2))/(2*A)
//
//    Author(s):  B. van Eijk                     Created at:  15/11/89
//                                                Last Update: 12/06/13
//
//----------------------------------------------------------------------
double PHOTRI(double A,double B,double C){
  double DA,DB,DC,DAPB,DAMB,DTRIAN;
  DA=A;
  DB=B;
  DC=C;
  DAPB=DA+DB;
  DAMB=DA-DB;
  DTRIAN=sqrt((DAMB-DC)*(DAPB+DC)*(DAMB+DC)*(DAPB-DC));
  return DTRIAN/(DA+DA);
}
//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays calculation of ANgle '1'
//
//    Purpose:  Calculate angle from X and Y
//
//    Input Parameters:  X, Y
//
//    Output Parameter:  Function value
//
//    Author(s):  S. Jadach                       Created at:  01/01/89
//                B. van Eijk                     Last Update: 12/06/13
//
//----------------------------------------------------------------------
double PHOAN1(double X,double Y){

  double phoan1 = 0.0;
  
  // we may want to use phpico_.pi phpico_.twopi defined in Photos::initialize()
  static double PI=3.14159265358979324, TWOPI=6.28318530717958648;
 
  if (fabs(Y)<fabs(X)){
    phoan1=atan(fabs(Y/X));
    if (X<0.0) phoan1=PI-phoan1;
  }
  else phoan1=acos(X/sqrt(X*X+Y*Y));
  //
  if (Y<0.0) phoan1=TWOPI-phoan1;
  return phoan1;
 
}

//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays calculation of ANgle '2'
//
//    Purpose:  Calculate angle from X and Y
//
//    Input Parameters:  X, Y
//
//    Output Parameter:  Function value
//
//    Author(s):  S. Jadach                       Created at:  01/01/89
//                B. van Eijk                     Last Update: 12/06/13
//
//----------------------------------------------------------------------
double PHOAN2(double X,double Y){

  double phoan2 = 0.0;

  // we may want to use phpico_.pi phpico_.twopi defined in Photos::initialize()
  static double PI=3.14159265358979324, TWOPI=6.28318530717958648;

  if (fabs(Y)<fabs(X)){
    phoan2=atan(fabs(Y/X));
    if (X<0.0) phoan2=PI-phoan2;
  }
  else phoan2=acos(X/sqrt(X*X+Y*Y));
  return phoan2;
}

//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays ROtation routine '2'
//
//    Purpose:  Rotate  x and z components  of vector PVEC  around angle
//              'ANGLE'.
//
//    Input Parameters:  ANGLE, PVEC
//
//    Output Parameter:  PVEC
//
//    Author(s):  S. Jadach                       Created at:  01/01/89
//                B. van Eijk                     Last Update: 12/06/13
//
//----------------------------------------------------------------------
void PHORO2(double ANGLE,double PVEC[4]){
  int j=1;  // convention of indices of Riemann space must be preserved.

  double CS,SN;
  CS= cos(ANGLE)*PVEC[1-j]+sin(ANGLE)*PVEC[3-j];
  SN=-sin(ANGLE)*PVEC[1-j]+cos(ANGLE)*PVEC[3-j];
  PVEC[1-j]=CS;
  PVEC[3-j]=SN;
}

//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays ROtation routine '3'
//
//    Purpose:  Rotate  x and y components  of vector PVEC  around angle
//              'ANGLE'.
//
//    Input Parameters:  ANGLE, PVEC
//
//    Output Parameter:  PVEC
//
//    Author(s):  S. Jadach     RO                 Created at:  01/01/89
//                B. van Eijk                     Last Update: 12/06/13
//
//----------------------------------------------------------------------
void PHORO3(double ANGLE,double PVEC[4]){
  int j=1;  // convention of indices of Riemann space must be preserved.
  double CS,SN;
  CS=cos(ANGLE)*PVEC[1-j]-sin(ANGLE)*PVEC[2-j];
  SN=sin(ANGLE)*PVEC[1-j]+cos(ANGLE)*PVEC[2-j];
  PVEC[1-j]=CS;
  PVEC[2-j]=SN;
}

//----------------------------------------------------------------------
//
//
//    PHOB:     PHotosBoost
//
//    Purpose:  Boosts VEC to (MODE=1)  rest frame of PBOOS1;  
//              or back (MODE=1)
//
//    Input Parameters:   MODE,PBOOS1,VEC
//
//    Output Parameters:  VEC
//
//    Author(s):                                  Created at:  08/12/05
//                Z. Was                          Last Update: 13/06/13
//
//----------------------------------------------------------------------

void PHOB(int MODE,double PBOOS1[4],double vec[4]){
  double BET1[3],GAM1,PB;
  static int j0=1;
  int I,J;


  PB=sqrt(PBOOS1[4-j0]*PBOOS1[4-j0]-PBOOS1[3-j0]*PBOOS1[3-j0]-PBOOS1[2-j0]*PBOOS1[2-j0]-PBOOS1[1-j0]*PBOOS1[1-j0]);
  for( J=1; J<4;J++){
    if (MODE==1) BET1[J-j0]=-PBOOS1[J-j0]/PB;
    else BET1[J-j0]= PBOOS1[J-j0]/PB;
  }

  GAM1=PBOOS1[4-j0]/PB;

  //--
  //--   Boost vector 

  PB=BET1[1-j0]*vec[1-j0]+BET1[2-j0]*vec[2-j0]+BET1[3-j0]*vec[3-j0];

  for( J=1; J<4;J++) vec[J-j0]=vec[J-j0]+BET1[J-j0]*(vec[4-j0]+PB/(GAM1+1.0));
  vec[4-j0]=GAM1*vec[4-j0]+PB;
  //--
}


//     *******************************
// Boost along arbitrary axis (as implemented by Ronald Kleiss).
// The method is described in book of Bjorken and Drell
// p boosted into r  from actual frame to rest frame of q
// forth (mode = 1) or back (mode = -1).
// q must be a timelike, p may be arbitrary.
void bostdq(int mode,double qq[4],double pp[4],double r[4]){
  double q[4],p[4],amq,fac;
  static int i=1;
  int k;

  for(k=1;k<=4;k++){
    p[k-i]=pp[k-i];
    q[k-i]=qq[k-i];
  }
  amq =sqrt(q[4-i]*q[4-i]-q[1-i]*q[1-i]-q[2-i]*q[2-i]-q[3-i]*q[3-i]);

  if    (mode == -1){
    r[4-i] = (p[1-i]*q[1-i]+p[2-i]*q[2-i]+p[3-i]*q[3-i]+p[4-i]*q[4-i])/amq;
    fac  = (r[4-i]+p[4-i])/(q[4-i]+amq);
  }
  else if(mode ==  1){
    r[4-i] =(-p[1-i]*q[1-i]-p[2-i]*q[2-i]-p[3-i]*q[3-i]+p[4-i]*q[4-i])/amq;
    fac  =-(r[4-i]+p[4-i])/(q[4-i]+amq);
  }
  else{
    cout << " ++++++++ wrong mode in boostdq " << endl;
    exit(0);
  }
  r[1-i]=p[1-i]+fac*q[1-i];
  r[2-i]=p[2-i]+fac*q[2-i];
  r[3-i]=p[3-i]+fac*q[3-i];
}


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

