#include <cmath>
#include <stdio.h>
#include <stdlib.h>  




// after investigations PHORO3 of PhotosUtilities.cxx will be used instead
// but it must be checked first if it works

extern "C" void rotod3_(double *pANGLE,double PVEC[4],double QVEC[4]){
  const double &ANGLE = *pANGLE;
 
  int j=1;  // convention of indices of Riemann space must be preserved.
  double CS,SN;
  //  printf ("%5.2f\n",cos(ANGLE));
  CS=cos(ANGLE)*PVEC[1-j]-sin(ANGLE)*PVEC[2-j];
  SN=sin(ANGLE)*PVEC[1-j]+cos(ANGLE)*PVEC[2-j];

  QVEC[1-j]=CS;
  QVEC[2-j]=SN;
  QVEC[3-j]=PVEC[3-j];
  QVEC[4-j]=PVEC[4-j];
}



extern "C" void   rotod2_(double *PH1,double PVEC[4],double QVEC[4]){
  const double &PHI = *PH1;
  double RVEC[4];
  int j=1;  // convention of indices of Riemann space must be preserved.
  double CS,SN;

  CS=cos(PHI);
  SN=sin(PHI);

  RVEC[1-j]=PVEC[1-j];
  RVEC[2-j]=PVEC[2-j];
  RVEC[3-j]=PVEC[3-j];
  RVEC[4-j]=PVEC[4-j];

  QVEC[1-j]= CS*RVEC[1-j]+SN*RVEC[3-j];
  QVEC[2-j]=RVEC[2-j];
  QVEC[3-j]=-SN*RVEC[1-j]+CS*RVEC[3-j];
  QVEC[4-j]=RVEC[4-j];
  //   printf ("%15.12f %15.12f %15.12f %15.12f \n",QVEC[0],QVEC[1],QVEC[2],QVEC[3]);
  // exit(-1);
}

