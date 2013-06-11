#include "Photos.h"
#include <cmath>
using std::max;
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

void fill_val(int beg, int end, double* array, double value); // Forward declaration

//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation  in decays function for SPIn determina-
//              tion
//
//    Purpose:  Calculate  the spin  of particle  with  code IDHEP.  The
//              code  of the particle  is  defined  by the Particle Data
//              Group in Phys. Lett. B204 (1988) 1.
//
//    Input Parameter:   IDHEP
//
//    Output Parameter:  Funtion  value = spin  of  particle  with  code
//                       IDHEP
//
//    Author(s):  E. Barberio and B. van Eijk     Created at:  29/11/89
//                                                Last update: 10/06/13
//
//----------------------------------------------------------------------
double PHOSPI(int idhep){
  static double SPIN[100] = { 0 };
  static int j=0;  
  //--
  //--   Array 'SPIN' contains the spin  of  the first 100 particles accor-
  //--   ding to the PDG particle code...

  if(j==0) // initialization
    {   
      j=1;
      fill_val(0 ,  8, SPIN, 0.5);
      fill_val(8 ,  9, SPIN, 1.0);
      fill_val(9 , 10, SPIN, 0.0);
      fill_val(10, 18, SPIN, 0.5);
      fill_val(18, 20, SPIN, 0.0);
      fill_val(20, 24, SPIN, 1.0);
      fill_val(24,100, SPIN, 0.0);
    }

  int idabs=abs(idhep);
  //--
  //--   Spin of quark, lepton, boson etc....
  if (idabs-1<100) return SPIN[idabs-1];

  //--   ...other particles, however...
  double xx=((idabs % 10)-1.0)/2.0;
  //--
  //--   ...K_short and K_long are special !!
  xx=max(xx,0.0);
  return xx;
}

//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays CHArge determination
//
//    Purpose:  Calculate the charge  of particle  with code IDHEP.  The
//              code  of the  particle  is  defined by the Particle Data
//              Group in Phys. Lett. B204 (1988) 1.
//
//    Input Parameter:   IDHEP
//
//    Output Parameter:  Funtion value = charge  of  particle  with code
//                       IDHEP
//
//    Author(s):  E. Barberio and B. van Eijk     Created at:  29/11/89
//                                                Last update: 11/06/13
//
//----------------------------------------------------------------------
double PHOCHA(int idhep){
  static double CHARGE[101] = { 0 };
  static int j=0;  
  //--
  //--   Array 'SPIN' contains the spin  of  the first 100 particles accor-
  //--   ding to the PDG particle code...

  if(j==0) // initialization
    {   
      j=1;
      fill_val(0 ,  1, CHARGE, 0.0         );
      fill_val(1 ,  2, CHARGE,-0.3333333333);
      fill_val(2 ,  3, CHARGE, 0.6666666667);
      fill_val(3 ,  4, CHARGE,-0.3333333333);
      fill_val(4 ,  5, CHARGE, 0.6666666667);
      fill_val(5 ,  6, CHARGE,-0.3333333333);
      fill_val(6 ,  7, CHARGE, 0.6666666667);
      fill_val(7 ,  8, CHARGE,-0.3333333333);
      fill_val(8 ,  9, CHARGE, 0.6666666667);
      fill_val(9 , 11, CHARGE, 0.0         );
      fill_val(11 ,12, CHARGE,-1.0         );
      fill_val(12 ,13, CHARGE, 0.0         );
      fill_val(13 ,14, CHARGE,-1.0         );
      fill_val(14, 15, CHARGE, 0.0         );
      fill_val(15 ,16, CHARGE,-1.0         );
      fill_val(16, 17, CHARGE, 0.0         );
      fill_val(17 ,18, CHARGE,-1.0         );
      fill_val(18, 24, CHARGE, 0.0         );
      fill_val(24, 25, CHARGE, 1.0         );
      fill_val(25, 37, CHARGE, 0.0         );
      fill_val(37, 38, CHARGE, 1.0         );
      fill_val(38,101, CHARGE, 0.0         );
    }

  int idabs=abs(idhep);
  double phoch=0.0;

  //--
  //--   Charge of quark, lepton, boson etc....
  if (idabs<=100) phoch=CHARGE[idabs];
  else {
    int Q3= idabs/1000 % 10;
    int Q2= idabs/100  % 10;
    int Q1= idabs/10   % 10;
    if (Q3==0){
      //--
      //-- ...meson...
      if(Q2 % 2==0) phoch=CHARGE[Q2]-CHARGE[Q1];
      else          phoch=CHARGE[Q1]-CHARGE[Q2];
    }
    else{
      //--
      //--   ...diquarks or baryon.
      phoch=CHARGE[Q1]+CHARGE[Q2]+CHARGE[Q3];
    }
  }
  //--
  //--   Find the sign of the charge...
  if (idhep<0.0) phoch=-phoch;
  if (phoch*phoch<0.000001) phoch=0.0;
  
  return phoch;
}




void fill_val(int beg, int end, double* array, double value) 
{
  for (int i = beg; i < end; i++)
    array[i] = value;
}

