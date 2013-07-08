#include "Photos.h"
#include <cmath>
#include <iostream>
#include "f_Init.h"
#include "PH_HEPEVT_Interface.h"
#include "PhotosUtilities.h"
using std::cout;
using std::endl;
using std::max;
using namespace Photospp;
using namespace PhotosUtilities;
// Declaration of structs defined in f_Init.h
struct PHOSTA phosta_;
struct PHLUPY phlupy_;
struct TOFROM tofrom_;
struct PHNUM  phnum_;
struct PHOLUN pholun_;
struct PHOREST phorest_;

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

//void fill_val(int beg, int end, double* array, double value); // Forward declaration

void fill_val(int beg, int end, double* array, double value) 
{
  for (int i = beg; i < end; i++)
    array[i] = value;
}

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



// --- can be used with  VARIANT A. For B use  PHINT1 or 2 --------------
//----------------------------------------------------------------------
//
//    PHINT:   PHotos universal INTerference correction weight
//
//    Purpose:  calculates correction weight as expressed by
//               formula (17) from CPC 79 (1994), 291. 
//
//    Input Parameters:  Common /PHOEVT/, with photon added.
//                                          
//    Output Parameters: correction weight
//
//    Author(s):  Z. Was, P.Golonka               Created at:  19/01/05
//                                                Last Update: 23/06/13
//
//----------------------------------------------------------------------

double PHINT(int IDUM){

  double PHINT2;
  double EPS1[4],EPS2[4],PH[4],PL[4];
  static int i=1;
  int K,L;
  //      DOUBLE PRECISION EMU,MCHREN,BETA,phophs_.costhg,MPASQR,XPH, XC1, XC2
  double  XNUM1,XNUM2,XDENO,XC1,XC2;

  //      REAL*8 PHOCHA
  //--

  //       Calculate polarimetric vector: ph, eps1, eps2 are orthogonal

  for( K=1;K<=4;K++){
    PH[K-i]= phoevt_.phep[phoevt_.nhep-i][K-i];
    EPS2[K-i]=1.0;
  }


  PHOEPS(PH,EPS2,EPS1);
  PHOEPS(PH,EPS1,EPS2);
    
 
  XNUM1=0.0;
  XNUM2=0.0;
  XDENO=0.0;

  for( K=phoevt_.jdahep[1-i][1-i]; K<=phoevt_.nhep-1;K++){  //! or jdahep[1-i][2-i]
      
    // momenta of charged particle in PL

    for( L=1;L<=4;L++) PL[L-i]=phoevt_.phep[K-i][L-i]; 

    // scalar products: epsilon*p/k*p

    XC1 = - PHOCHA(phoevt_.idhep[K-i]) * 
         ( PL[1-i]*EPS1[1-i] + PL[2-i]*EPS1[2-i] + PL[3-i]*EPS1[3-i] ) / 
	 ( PH[4-i]*PL[4-i]   - PH[1-i]*PL[1-i]   - PH[2-i]*PL[2-i] - PH[3-i]*PL[3-i] );
     
    XC2 = - PHOCHA(phoevt_.idhep[K-i]) * 
         ( PL[1-i]*EPS2[1-i] + PL[2-i]*EPS2[2-i] + PL[3-i]*EPS2[3-i] ) / 
	 ( PH[4-i]*PL[4-i]   - PH[1-i]*PL[1-i]   - PH[2-i]*PL[2-i] - PH[3-i]*PL[3-i] );
	

    // accumulate the currents
    XNUM1  = XNUM1+XC1;
    XNUM2  = XNUM2+XC2;

    XDENO = XDENO + XC1*XC1 + XC2*XC2;
  }

  PHINT2=(XNUM1*XNUM1 + XNUM2*XNUM2) / XDENO;
  return (XNUM1*XNUM1 + XNUM2*XNUM2) / XDENO;

}



//----------------------------------------------------------------------
//
//    PHINT:   PHotos INTerference (Old version kept for tests only.
//
//    Purpose:  Calculates interference between emission of photons from
//              different possible chaged daughters stored in
//              the  HEP common /PHOEVT/.  
//
//    Input Parameter:    commons /PHOEVT/ /PHOMOM/ /PHOPHS/
//    
//
//    Output Parameters:  
//                        
//
//    Author(s):  Z. Was,                         Created at:  10/08/93
//                                                Last Update: 15/03/99
//
//----------------------------------------------------------------------

double PHINT1(int IDUM){
# define pho phoevt_
  double PHINT;

  /*
      DOUBLE PRECISION phomom_.mchsqr,phomom_.mnesqr
      REAL*8 PNEUTR
      COMMON/PHOMOM/phomom_.mchsqr,phomom_.mnesqr,PNEUTR(5)
      DOUBLE PRECISION phophs_.costhg,SINTHG
      REAL*8 XPHMAX,phophs_.xphoto
      COMMON/PHOPHS/XPHMAX,phophs_.xphoto,phophs_.costhg,SINTHG

  */
  double MPASQR,XX,BETA;
  bool IFINT;
  int K,IDENT; 
  static int i=1;
  //
  for(K=pho.jdahep[1-i][2-i]; K>=pho.jdahep[1-i][1-i];K--){
    if(pho.idhep[K-i]!=22){
      IDENT=K;
      break;
    }
  }

  // check if there is a photon
  IFINT= pho.nhep>IDENT;
  // check if it is two body + gammas reaction
  IFINT= IFINT && (IDENT-pho.jdahep[1-i][1-i])==1;
  // check if two body was particle antiparticle
  IFINT= IFINT && pho.idhep[pho.jdahep[1-i][1-i]-i] == -pho.idhep[IDENT-i];
  // check if particles were charged
  IFINT= IFINT && PHOCHA(pho.idhep[IDENT-i]) != 0;
  // calculates interference weight contribution
  if(IFINT){
    MPASQR = pho.phep[1-i][5-i]*pho.phep[1-i][5-i];
    XX=4.0*phomom_.mchsqr/MPASQR*(1.0-phophs_.xphoto)/(1.0-phophs_.xphoto+(phomom_.mchsqr-phomom_.mnesqr)/MPASQR)/(1.0-phophs_.xphoto+(phomom_.mchsqr-phomom_.mnesqr)/MPASQR);
    BETA=sqrt(1.0-XX);
    PHINT  = 2.0/(1.0+phophs_.costhg*phophs_.costhg*BETA*BETA);
  }
  else{
    PHINT  = 1.0;
  }

  return  PHINT;
}


//----------------------------------------------------------------------
//
//    PHINT:   PHotos INTerference
//
//    Purpose:  Calculates interference between emission of photons from
//              different possible chaged daughters stored in
//              the  HEP common /PHOEVT/. 
//
//    Input Parameter:    commons /PHOEVT/ /PHOMOM/ /PHOPHS/
//    
//
//    Output Parameters:  
//                        
//
//    Author(s):  Z. Was,                         Created at:  10/08/93
//                                                Last Update: 
//
//----------------------------------------------------------------------

double PHINT2(int IDUM){
# define pho phoevt_

  /*
      DOUBLE PRECISION phomom_.mchsqr,phomom_.mnesqr
      REAL*8 PNEUTR
      COMMON/PHOMOM/phomom_.mchsqr,phomom_.mnesqr,PNEUTR(5)
      DOUBLE PRECISION phophs_.costhg,SINTHG
      REAL*8 XPHMAX,phophs_.xphoto
      COMMON/PHOPHS/XPHMAX,phophs_.xphoto,phophs_.costhg,SINTHG
  */
  double MPASQR,XX,BETA,pq1[4],pq2[4],pphot[4];
  double SS,PP2,PP,E1,E2,q1,q2,costhe,PHINT;
  bool IFINT;
  int K,k,IDENT; 
  static int i=1;
  //
  for(K=pho.jdahep[1-i][2-i]; K>=pho.jdahep[1-i][1-i];K--){
    if(pho.idhep[K-i]!=22){
      IDENT=K;
      break;
    }
  }

  // check if there is a photon
  IFINT= pho.nhep>IDENT;
  // check if it is two body + gammas reaction
  IFINT= IFINT&&(IDENT-pho.jdahep[1-i][1-i])==1;
  // check if two body was particle antiparticle (we improve on it !
  //      IFINT= IFINT.AND.pho.idhep(JDAPHO(1,1)).EQ.-pho.idhep(IDENT)
  // check if particles were charged
  IFINT= IFINT&&fabs(PHOCHA(pho.idhep[IDENT-i]))>0.01;
  // check if they have both charge
  IFINT= IFINT&&fabs(PHOCHA(pho.idhep[pho.jdahep[1-i][1-i]-i]))>0.01;
  // calculates interference weight contribution
  if(IFINT){
    MPASQR = pho.phep[1-i][5-i]*pho.phep[1-i][5-i];
    XX=4.0*phomom_.mchsqr/MPASQR*(1.0-phophs_.xphoto)/pow(1.-phophs_.xphoto+(phomom_.mchsqr-phomom_.mnesqr)/MPASQR,2);
    BETA=sqrt(1.0-XX);
    PHINT  = 2.0/(1.0+phophs_.costhg*phophs_.costhg*BETA*BETA);
    SS =MPASQR*(1.0-phophs_.xphoto);
    PP2=((SS-phomom_.mchsqr-phomom_.mnesqr)*(SS-phomom_.mchsqr-phomom_.mnesqr)-4*phomom_.mchsqr*phomom_.mnesqr)/SS/4;
    PP =sqrt(PP2);
    E1 =sqrt(PP2+phomom_.mchsqr);
    E2 =sqrt(PP2+phomom_.mnesqr);
    PHINT= (E1+E2)*(E1+E2)/((E2+phophs_.costhg*PP)*(E2+phophs_.costhg*PP)+(E1-phophs_.costhg*PP)*(E1-phophs_.costhg*PP));
    // return PHINT;
    //
    q1=PHOCHA(pho.idhep[pho.jdahep[1-i][1-i]-i]);
    q2=PHOCHA(pho.idhep[IDENT-i]);
    for( k=1;k<=4;k++){
      pq1[k-i]=pho.phep[pho.jdahep[1-i][1-i]-i][k-i];
      pq2[k-i]=pho.phep[pho.jdahep[1-i][1-i]+1-i][k-i];
      pphot[k-i]=pho.phep[pho.nhep-i][k-i];
    }
    costhe=(pphot[1-i]*pq1[1-i]+pphot[2-i]*pq1[2-i]+pphot[3-i]*pq1[3-i]);
    costhe=costhe/sqrt(pq1[1-i]*pq1[1-i]+pq1[2-i]*pq1[2-i]+pq1[3-i]*pq1[3-i]);
    costhe=costhe/sqrt(pphot[1-i]*pphot[1-i]+pphot[2-i]*pphot[2-i]+pphot[3-i]*pphot[3-i]);
    //
    // --- this IF checks whether JDAPHO(1,1) was MCH or MNE. 
    // --- phophs_.costhg angle (and in-generation variables) may be better choice 
    // --- than costhe. note that in the formulae below amplitudes were 
    // --- multiplied by (E2+phophs_.costhg*PP)*(E1-phophs_.costhg*PP). 
    if(phophs_.costhg*costhe>0){

      PHINT= pow(q1*(E2+phophs_.costhg*PP)-q2*(E1-phophs_.costhg*PP),2)/(q1*q1*(E2+phophs_.costhg*PP)*(E2+phophs_.costhg*PP)+q2*q2*(E1-phophs_.costhg*PP)*(E1-phophs_.costhg*PP));
    }
    else{

      PHINT= pow(q1*(E1-phophs_.costhg*PP)-q2*(E2+phophs_.costhg*PP),2)/(q1*q1*(E1-phophs_.costhg*PP)*(E1-phophs_.costhg*PP)+q2*q2*(E2+phophs_.costhg*PP)*(E2+phophs_.costhg*PP));
    }
  }
  else{
    PHINT  = 1.0;
  }
  return PHINT;
}


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
//    Author(s):  S. Jadach                       Created at:  01/01/89
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
//    PHOTOS:   PHOton radiation in decays ERRror handling
//
//    Purpose:  Inform user  about (fatal) errors and warnings generated
//              by either the user or the program.
//
//    Input Parameters:   IMES, TEXT, DATA
//
//    Output Parameters:  None
//
//    Author(s):  B. van Eijk                     Created at:  29/11/89
//                                                Last Update: 18/06/13
//
//----------------------------------------------------------------------
void PHOERR(int IMES,char *TEXT,double DATA){

  static int IERROR=0;
  double  SDATA;
  static int PHOMES=10;
  static int i=1;
  char star80[81]= "********************************************************************************";

  if (IMES<=PHOMES) phosta_.status[IMES-i]=phosta_.status[IMES-i]+1;
// 
//    Count number of non-fatal errors...
  if ((IMES ==  6) && (phosta_.status[IMES-i]>=2)) return;
  if ((IMES == 10) && (phosta_.status[IMES-i]>=2)) return;
  SDATA=DATA;
  //  int PHLUN=(int)pholun_.phlun;
  bool IFSTOP=phosta_.ifstop;
  FILE *PHLUN = stdout;
  int furthA=0;
  fprintf(PHLUN,"%s\n",star80);
  fprintf(PHLUN,"*\n");  //9120
  //      GOTO (10,20,30,40,50,60,70,80,90,100),IMES

  switch(IMES){
  case 1:
    fprintf(PHLUN,"* %s: Too many charged Particles, NCHARG = %6i\n", TEXT,(int)SDATA);   //I6
    furthA= 110;
    break;
  case 2:
    fprintf(PHLUN,"* %s: Too much Bremsstrahlung required, PRSOFT = %15.6f\n", TEXT,SDATA);//F15.6
    furthA= 110;
    break;
  case 3:
    fprintf(PHLUN,"* %s: Combined Weight is exceeding 1., Weight = %15.6f\n", TEXT,SDATA);   //F15.6
    furthA= 110;
    break;
  case 4:
    fprintf(PHLUN,"* %s: Error in Rescaling charged and neutral Vectors\n", TEXT);
    furthA= 110;
    break;
  case 5:
    fprintf(PHLUN,"* %s: Non matching charged Particle Pointer, NCHARG = %5i\n", TEXT,(int)SDATA);  //I5
    furthA= 110;
    break;
  case 6:
    fprintf(PHLUN,"* %s: Do you really work with a Particle of Spin: %4.1f\n", TEXT,SDATA);   //F4.1
    furthA= 130;
    break;
  case 7:
    fprintf(PHLUN,"* %s: Stack Length exceeded, NSTACK = %5i\n", TEXT,(int)(SDATA));//I5
    furthA= 110;
    break;
  case 8:
    fprintf(PHLUN,"* %s: Random Number Generator Seed(1) out of Range: %8i\n", TEXT,(int)SDATA);//I8
    furthA= 110;
    break;
  case 9:
    fprintf(PHLUN,"* %s: Random Number Generator Seed(2) out of Range: %8i\n", TEXT,(int)SDATA);//I8
    furthA= 110;
    break;
  case 10:
    fprintf(PHLUN,"* %s: Available Phase Space below Cut-off: %15.6f GeV/c^2\n", TEXT,SDATA);//F15.6
    furthA= 130;
    break;
  default:
    fprintf(PHLUN,"* Funny Error Message: %4i ! What to do ?\n", IMES);//I4
    furthA= 120;
    break;
  }

 switch(furthA){
 case 110:
   fprintf(PHLUN,"* Fatal Error Message, I stop this Run !\n");
   fprintf(PHLUN,"*\n"); //9120
   fprintf(PHLUN,"%s\n",star80);
   if (IFSTOP){ 
     exit(0);
   }
   else{
     fprintf(PHLUN,"*\n"); //9120
     fprintf(PHLUN,"%s\n",star80);
     break;
   }      
 case 120:
   IERROR=IERROR+1;
   if (IERROR>=10){
     fprintf(PHLUN,"* 10 Error Messages generated, I stop this Run !\n");
     fprintf(PHLUN,"*\n");//9120
     fprintf(PHLUN,"%s\n",star80);
     if (IFSTOP){
       exit(0);
     }
     else{
       fprintf(PHLUN,"*\n"); //9120
       fprintf(PHLUN,"%s\n",star80);
       break;
     }
   }  
 case 130:
  fprintf(PHLUN,"*\n");  //9120
  fprintf(PHLUN,"%s\n",star80);
  break;
 }
 return;


 //9120 FORMAT(1H ,'*',T81,'*')
 // 9140 FORMAT(1H ,'* Fatal Error Message, I stop this Run !',T81,'*')
 // 9150 FORMAT(1H ,'* 10 Error Messages generated, I stop this Run !',T81,
 //     &'*')
}


//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays run summary REPort
//
//    Purpose:  Inform user about success and/or restrictions of PHOTOS
//              encountered during execution.
//
//    Input Parameters:   Common /PHOSTA/
//
//    Output Parameters:  None
//
//    Author(s):  B. van Eijk                     Created at:  10/01/92
//                                                Last Update: 18/06/13
//
//----------------------------------------------------------------------
void PHOREP(){
  static int PHOMES=10;
  int I;
  bool ERROR=false;
  //  int PHLUN=(int)pholun_.phlun;
  char star80[81]= "********************************************************************************";
  char X26[27] = "                          ";
  char EQ25[26]= "=========================";
  char X30[31] = "                              ";
  char X22[23] = "                      ";
  char X23[24 ]= "                       ";
  char X16[17] = "                ";
  FILE *PHLUN = stdout;
  fprintf(PHLUN," \n");
  fprintf(PHLUN,"%s\n",star80);
  fprintf(PHLUN,"*\n");
  fprintf(PHLUN,"* %s %s\n",X26,EQ25);
  fprintf(PHLUN,"* %s PHOTOS Run Summary\n",X30);
  fprintf(PHLUN,"* %s %s\n",X26,EQ25);
  fprintf(PHLUN,"*\n");
  for(I=1;I<=PHOMES;I++){

    if (phosta_.status[I-1] == 0) break;
    if ((I == 6)|| (I == 10)){
      fprintf(PHLUN,"* %s Warning # %2i  occured %6i times\n",X22, I,phosta_.status[I-1]); // I2 I6 
    }
    else{
      ERROR=true;
      fprintf(PHLUN,"* %s Error # %2i occured %6i  times\n",X23, I,phosta_.status[I-1]);// I2 I6
    }	      
  }

  if (!ERROR) fprintf(PHLUN,"* %s PHOTOS Execution has successfully terminated\n",X16);
  fprintf(PHLUN,"*\n");
  fprintf(PHLUN,"%s\n",star80);
  return;

//      RETURN
// 9000 FORMAT(1H1)
// 9010 FORMAT(1H ,80('*'))
// 9020 FORMAT(1H ,'*',T81,'*')
// 9030 FORMAT(1H ,'*',26X,25('='),T81,'*')
// 9040 FORMAT(1H ,'*',30X,'PHOTOS Run Summary',T81,'*')
// 9050 FORMAT(1H ,'*',22X,'Warning #',I2,' occured',I6,' times',T81,'*')
// 9060 FORMAT(1H ,'*',23X,'Error #',I2,' occured',I6,' times',T81,'*')
// 9070 FORMAT(1H ,'*',16X,'PHOTOS Execution has successfully terminated',
//     &T81,'*')
}

//*****************************************************************
//*****************************************************************
//*****************************************************************
// beginning of the class of methods reading from  PH_HEPEVT
//*****************************************************************
//*****************************************************************
//*****************************************************************


void GETIDEIDF(int *IDE,int *IDF){
  // this method provides information on particles ID-s to be used
  // by  the Z-ME class in calculation of A_FB
  static int i=1;  
  *IDE=ph_hepevt_.idhep[1-i];
  *IDF=ph_hepevt_.idhep[4-i];
  if(abs(ph_hepevt_.idhep[4-i])==abs(ph_hepevt_.idhep[3-i])) *IDF=ph_hepevt_.idhep[3-i];
}


//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in decays event DuMP routine
//
//    Purpose:  Print event record.
//
//    Input Parameters:   Common /PH_HEPEVT/
//
//    Output Parameters:  None
//
//    Author(s):  B. van Eijk                     Created at:  05/06/90
//                                                Last Update: 20/06/13
//
//----------------------------------------------------------------------
void PHODMP(){

  double  SUMVEC[5];
  int I,J;
  static int i=1;
  char star80[81]= "********************************************************************************";
  char eq80[81]  = "================================================================================";
  char X26[27] = "                          ";
  char EQ25[26]= "=========================";
  char X30[31] = "                              ";
  char X29[30] = "                             ";
  char X22[23] = "                      ";
  char X23[24 ]= "                       ";
  char X16[17] = "                ";
  char X1[2] = " ";
  char X2[3] = "  ";
  char X3[4] = "   ";
  char X4[5] = "    ";
  char X5[6] = "     ";
  char X6[7] = "      ";
  char X7[8] = "       ";
  char X9[10]= "         ";
  FILE *PHLUN = stdout;

  for(I=0;I<5;I++)  SUMVEC[I]=0.0;
  //--
  //--   Print event number...
  fprintf(PHLUN,"%s",eq80);
  fprintf(PHLUN,"%s Event No.: %10i\n",X29,ph_hepevt_.nevhep);
  fprintf(PHLUN,"%s Particle Parameters\n",X6);
  fprintf(PHLUN,"%s Nr %s Type %s Parent(s) %s Daughter(s) %s Px %s Py %s Pz %s E %s Inv. M.\n",X1,X3,X3,X2,X6,X7,X7,X7,X4);
  for(I=1;I<=ph_hepevt_.nhep;I++){ 
    //--
    //--   For 'stable particle' calculate vector momentum sum
    if (ph_hepevt_.jdahep[I-i][1-i]==0){
      for(J=1; J<=4;J++){
	SUMVEC[J-i]=SUMVEC[J-i]+ph_hepevt_.phep[I-i][J-i];
      }
      if (ph_hepevt_.jmohep[I-i][2-i]==0){
	fprintf(PHLUN,"%4i %7i %s %4i %s Stable %9.2f %9.2f %9.2f %9.2f %9.2f\n" ,  I,ph_hepevt_.idhep[I-i],X3,ph_hepevt_.jmohep[I-i][1-i],X7,ph_hepevt_.phep[I-i][1-i],ph_hepevt_.phep[I-i][2-i],ph_hepevt_.phep[I-i][3-i],ph_hepevt_.phep[I-i][4-i],ph_hepevt_.phep[I-i][5-i]);
      }
      else{
	fprintf(PHLUN,"%4i %7i %4i - %4i %s Stable %9.2f %9.2f %9.2f %9.2f %9.2f\n",I,ph_hepevt_.idhep[I-i],ph_hepevt_.jmohep[I-i][1-i],ph_hepevt_.jmohep[I-i][2-i], X4,ph_hepevt_.phep[I-i][1-i],ph_hepevt_.phep[I-i][2-i],ph_hepevt_.phep[I-i][3-i],ph_hepevt_.phep[I-i][4-i],ph_hepevt_.phep[I-i][5-i]);
      }
    }
    else{
      if(ph_hepevt_.jmohep[I-i][2-i]==0){
	fprintf(PHLUN,"%4i %7i %s %4i %s %4i - %4i %9.2f %9.2f %9.2f %9.2f %9.2f\n" ,  I,ph_hepevt_.idhep[I-i],X3,ph_hepevt_.jmohep[I-i][1-i],X2,ph_hepevt_.jdahep[I-i][1-i],ph_hepevt_.jdahep[I-i][2-i],ph_hepevt_.phep[I-i][1-i],ph_hepevt_.phep[I-i][2-i],ph_hepevt_.phep[I-i][3-i],ph_hepevt_.phep[I-i][4-i],ph_hepevt_.phep[I-i][5-i]);
      }
      else{
	fprintf(PHLUN,"%4i %7i %4i - %4i %4i - %4i %9.2f %9.2f %9.2f %9.2f %9.2f\n",  I,ph_hepevt_.idhep[I-i],ph_hepevt_.jmohep[I-i][1-i],ph_hepevt_.jmohep[I-i][2-i],ph_hepevt_.jdahep[I-i][1-i],ph_hepevt_.jdahep[I-i][2-i],ph_hepevt_.phep[I-i][1-i],ph_hepevt_.phep[I-i][2-i],ph_hepevt_.phep[I-i][3-i],ph_hepevt_.phep[I-i][4-i],ph_hepevt_.phep[I-i][5-i]);
      }
    }
  }
  SUMVEC[5]=sqrt(SUMVEC[4-i]*SUMVEC[4-i]-SUMVEC[1-i]*SUMVEC[1-i]-SUMVEC[2-i]*SUMVEC[2-i]-SUMVEC[3-i]*SUMVEC[3-i]);
  fprintf(PHLUN,"%s  Vector Sum: %9.2f %9.2f %9.2f %9.2f %9.2f\n",X23,SUMVEC[1-i],SUMVEC[2-i],SUMVEC[3-i],SUMVEC[4-i],SUMVEC[5-i]);




// 9030 FORMAT(1H ,I4,I7,3X,I4,9X,'Stable',2X,5F9.2)
//"%4i %7i %s  %4i %s Stable %s  %9.2f %9.2f %9.2f %9.2f %9.2f "  X3,9X,X2

  // 9050 FORMAT(1H ,I4,I7,3X,I4,6X,I4,' - ',I4,5F9.2)
  //"%4i %7i %s  %4i %s %4i  -  %4i  %9.2f %9.2f %9.2f %9.2f %9.2f "  X3,X6

 


  //"%4i %7i %4i  -  %4i %s Stable %s  %9.2f %9.2f %9.2f %9.2f %9.2f "  X5,X2


 //9060 FORMAT(1H ,I4,I7,I4,' - ',I4,2X,I4,' - ',I4,5F9.2)
  //"%4i %7i %4i  -  %4i %s %4i -   %4i %9.2f %9.2f %9.2f %9.2f %9.2f "  X2,
}



//----------------------------------------------------------------------
//
//    PHLUPAB:   debugging tool
//
//    Purpose:  NONE, eventually may printout content of the 
//              /PH_HEPEVT/ common
//
//    Input Parameters:   Common /PH_HEPEVT/ and /PHNUM/ 
//                        latter may have number of the event. 
//
//    Output Parameters:  None
//
//    Author(s):  Z. Was                          Created at:  30/05/93
//                                                Last Update: 20/06/13
//
//----------------------------------------------------------------------

void PHLUPAB(int IPOINT){
  char name[12] = "/PH_HEPEVT/";
  int I,J;
  static int IPOIN0=-5;
  static int i=1;
  double  SUM[5];
  FILE *PHLUN = stdout;

  if (IPOIN0<0){
    IPOIN0=400000; //  ! maximal no-print point
    phlupy_.ipoin =IPOIN0;
    phlupy_.ipoinm=400001; // ! minimal no-print point
  }
  
  if (IPOINT<=phlupy_.ipoinm||IPOINT>=phlupy_.ipoin ) return;
  if ((int)phnum_.iev<1000){
    for(I=1; I<=5;I++) SUM[I-i]=0.0;
     
    fprintf(PHLUN,"EVENT NR= %i WE ARE TESTING %s at IPOINT=%i \n",(int)phnum_.iev,name,IPOINT);
    fprintf(PHLUN,"  ID      p_x      p_y      p_z      E        m        ID-MO_DA1 ID-MO_DA2\n");
    I=1;
    fprintf(PHLUN,"%4i %14.9f %14.9f %14.9f %14.9f %14.9f %9i %9i\n", ph_hepevt_.idhep[I-i],ph_hepevt_.phep[1-i][I-i],ph_hepevt_.phep[2-i][I-i],ph_hepevt_.phep[3-i][I-i],ph_hepevt_.phep[4-i][I-i],ph_hepevt_.phep[5-i][I-i],ph_hepevt_.jdahep[1-i][I-i],ph_hepevt_.jdahep[2-i][I-i]);
    I=2;
    fprintf(PHLUN,"%4i %14.9f %14.9f %14.9f %14.9f %14.9f %9i %9i\n", ph_hepevt_.idhep[I-i],ph_hepevt_.phep[1-i][I-i],ph_hepevt_.phep[2-i][I-i],ph_hepevt_.phep[3-i][I-i],ph_hepevt_.phep[4-i][I-i],ph_hepevt_.phep[5-i][I-i],ph_hepevt_.jdahep[1-i][I-i],ph_hepevt_.jdahep[2-i][I-i]);
    fprintf(PHLUN," \n");
    for(I=3;I<=ph_hepevt_.nhep;I++){
      fprintf(PHLUN,"%4i %14.9f %14.9f %14.9f %14.9f %14.9f %9i %9i\n", ph_hepevt_.idhep[I-i],ph_hepevt_.phep[1-i][I-i],ph_hepevt_.phep[2-i][I-i],ph_hepevt_.phep[3-i][I-i],ph_hepevt_.phep[4-i][I-i],ph_hepevt_.phep[5-i][I-i],ph_hepevt_.jmohep[1-i][I-i],ph_hepevt_.jmohep[2-i][I-i]);
      for(J=1;J<=4;J++) SUM[J-i]=SUM[J-i]+ph_hepevt_.phep[J-i][I-i];
    }


    SUM[5]=sqrt(fabs(SUM[4-i]*SUM[4-i]-SUM[1-i]*SUM[1-i]-SUM[2-i]*SUM[2-i]-SUM[3-i]*SUM[3-i]));
    fprintf(PHLUN," SUM %14.9f %14.9f %14.9f %14.9f %14.9f\n",SUM[1-i],SUM[2-i],SUM[3-i],SUM[4-i],SUM[5-i]);

  }


	// 10   FORMAT(1X,'  ID      ','p_x      ','p_y      ','p_z      ',
	//$                   'E        ','m        ',
	//$                   'ID-MO_DA1','ID-MO DA2' )
  // 20   FORMAT(1X,I4,5(F14.9),2I9)
  //"%i4 %14.9f %14.9f %14.9f %14.9f %i9 i9"
	// 30   FORMAT(1X,' SUM',5(F14.9))
}









//----------------------------------------------------------------------
//
//    PHLUPA:   debugging tool
//
//    Purpose:  NONE, eventually may printout content of the 
//              /PHOEVT/ common
//
//    Input Parameters:   Common /PHOEVT/ and /PHNUM/ 
//                        latter may have number of the event. 
//
//    Output Parameters:  None
//
//    Author(s):  Z. Was                          Created at:  30/05/93
//                                                Last Update: 21/06/13
//
//----------------------------------------------------------------------

void PHLUPA(int IPOINT){
  char name[9] = "/PHOEVT/";
  int I,J;
  static int IPOIN0=-5;
  static int i=1;
  double  SUM[5];
  FILE *PHLUN = stdout;

  if (IPOIN0<0){
    IPOIN0=400000; //  ! maximal no-print point
    phlupy_.ipoin =IPOIN0;
    phlupy_.ipoinm=400001; // ! minimal no-print point
  }
  
  if (IPOINT<=phlupy_.ipoinm||IPOINT>=phlupy_.ipoin ) return;
  if ((int)phnum_.iev<1000){
    for(I=1; I<=5;I++) SUM[I-i]=0.0;
     
    fprintf(PHLUN,"EVENT NR= %i WE ARE TESTING %s at IPOINT=%i \n",(int)phnum_.iev,name,IPOINT);
    fprintf(PHLUN,"  ID      p_x      p_y      p_z      E        m        ID-MO_DA1 ID-MO_DA2\n");
    I=1;
    fprintf(PHLUN,"%4i %14.9f %14.9f %14.9f %14.9f %14.9f %9i %9i\n", phoevt_.idhep[I-i],phoevt_.phep[1-i][I-i],phoevt_.phep[2-i][I-i],phoevt_.phep[3-i][I-i],phoevt_.phep[4-i][I-i],phoevt_.phep[5-i][I-i],phoevt_.jdahep[1-i][I-i],phoevt_.jdahep[2-i][I-i]);
    I=2;
    fprintf(PHLUN,"%4i %14.9f %14.9f %14.9f %14.9f %14.9f %9i %9i\n", phoevt_.idhep[I-i],phoevt_.phep[1-i][I-i],phoevt_.phep[2-i][I-i],phoevt_.phep[3-i][I-i],phoevt_.phep[4-i][I-i],phoevt_.phep[5-i][I-i],phoevt_.jdahep[1-i][I-i],phoevt_.jdahep[2-i][I-i]);
    fprintf(PHLUN," \n");
    for(I=3;I<=phoevt_.nhep;I++){
      fprintf(PHLUN,"%4i %14.9f %14.9f %14.9f %14.9f %14.9f %9i %9i\n", phoevt_.idhep[I-i],phoevt_.phep[1-i][I-i],phoevt_.phep[2-i][I-i],phoevt_.phep[3-i][I-i],phoevt_.phep[4-i][I-i],phoevt_.phep[5-i][I-i],phoevt_.jmohep[1-i][I-i],phoevt_.jmohep[2-i][I-i]);
      for(J=1;J<=4;J++) SUM[J-i]=SUM[J-i]+phoevt_.phep[J-i][I-i];
    }
  

    SUM[5]=sqrt(fabs(SUM[4-i]*SUM[4-i]-SUM[1-i]*SUM[1-i]-SUM[2-i]*SUM[2-i]-SUM[3-i]*SUM[3-i]));
    fprintf(PHLUN," SUM %14.9f %14.9f %14.9f %14.9f %14.9f\n",SUM[1-i],SUM[2-i],SUM[3-i],SUM[4-i],SUM[5-i]);

  }


	// 10   FORMAT(1X,'  ID      ','p_x      ','p_y      ','p_z      ',
	//$                   'E        ','m        ',
	//$                   'ID-MO_DA1','ID-MO DA2' )
  // 20   FORMAT(1X,I4,5(F14.9),2I9)
  //"%4i %14.9f %14.9f %14.9f %14.9f %9i %9i"
	// 30   FORMAT(1X,' SUM',5(F14.9))
}


void PHOtoRF(){
# define hep ph_hepevt_

  //      COMMON /PH_TOFROM/ QQ[4],XM,th1,fi1
  double PP[4],RR[4];

  int K,L;
  static int i=1;

  for(K=1;K<=4;K++){
    tofrom_.QQ[K-i]=0.0;
  }
  for( L=hep.jdahep[hep.jmohep[hep.nhep-i][1-i]-i][1-i];L<=hep.jdahep[hep.jmohep[hep.nhep-i][1-i]-i][2-i];L++){
    for(K=1;K<=4;K++){
      tofrom_.QQ[K-i]=tofrom_.QQ[K-i]+hep.phep[L-i][K-i];
    }
  }
  tofrom_.XM =tofrom_.QQ[4-i]*tofrom_.QQ[4-i]-tofrom_.QQ[3-i]*tofrom_.QQ[3-i]-tofrom_.QQ[2-i]*tofrom_.QQ[2-i]-tofrom_.QQ[1-i]*tofrom_.QQ[1-i];
  if(tofrom_.XM>0.0) tofrom_.XM=sqrt(tofrom_.XM);
  if(tofrom_.XM<=0.0) return;

  for(L=1;L<=hep.nhep;L++){
    for(K=1;K<=4;K++){       
      PP[K-i]=hep.phep[L-i][K-i];
    }
    bostdq(1,tofrom_.QQ,PP,RR);
    for(K=1;K<=4;K++){     
      hep.phep[L-i][K-i]=RR[K-i];
    }
  }

  tofrom_.fi1=0.0;
  tofrom_.th1=0.0;
  if(fabs(hep.phep[1-i][1-i])+fabs(hep.phep[1-i][2-i])>0.0) tofrom_.fi1=PHOAN1(hep.phep[1-i][1-i],hep.phep[1-i][2-i]);
  if(fabs(hep.phep[1-i][1-i])+fabs(hep.phep[1-i][2-i])+fabs(hep.phep[1-i][3-i])>0.0)  
    tofrom_.th1=PHOAN2(hep.phep[1-i][3-i],sqrt(hep.phep[1-i][1-i]*hep.phep[1-i][1-i]+hep.phep[1-i][2-i]*hep.phep[1-i][2-i]));

  for(L=1;L<=hep.nhep;L++){ 
    for(K=1;K<=4;K++){       
      RR[K-i]=hep.phep[L-i][K-i];
    }
     
    PHORO3(-tofrom_.fi1,RR);
    PHORO2(-tofrom_.th1,RR);
    for(K=1;K<=4;K++){     
      hep.phep[L-i][K-i]=RR[K-i];
    }
  }
  
  return;
}

void PHOtoLAB(){
# define hep ph_hepevt_
  //  //      REAL*8 QQ(4),XM,th1,fi1
  //     COMMON /PH_TOFROM/ QQ,XM,th1,fi1
  double PP[4],RR[4];
  int K,L;
  static int i=1;
  
  if(tofrom_.XM<=0.0) return;


  for(L=1;L<=hep.nhep;L++){
    for(K=1;K<=4;K++){
      PP[K-i]=hep.phep[L-i][K-i];
    }

    PHORO2( tofrom_.th1,PP);
    PHORO3( tofrom_.fi1,PP);
    bostdq(-1,tofrom_.QQ,PP,RR);

    for(K=1;K<=4;K++){
      hep.phep[L-i][K-i]=RR[K-i];
    }
  }
  return;
}





//             2) GENERAL INTERFACE:
//                                      PHOTOS_GET
//                                      PHOTOS_MAKE


//   COMMONS:
//   NAME     USED IN SECT. # OF OC//     Comment
//   PHOQED   1) 2)            3      Flags whether emisson to be gen. 
//   PHOLUN   1) 4)            6      Output device number
//   PHOCOP   1) 3)            4      photon coupling & min energy
//   PHPICO   1) 3) 4)         5      PI & 2*PI
//   PHSEED   1) 4)            3      RN seed 
//   PHOSTA   1) 4)            3      Status information
//   PHOKEY   1) 2) 3)         7      Keys for nonstandard application
//   PHOVER   1)               1      Version info for outside
//   HEPEVT   2)               2      PDG common
//   PH_HEPEVT2)               8      PDG common internal
//   PHOEVT   2) 3)           10      PDG branch
//   PHOIF    2) 3)            2      emission flags for PDG branch 
//   PHOMOM   3)               5      param of char-neutr system
//   PHOPHS   3)               5      photon momentum parameters
//   PHOPRO   3)               4      var. for photon rep. (in branch)
//   PHOCMS   2)               3      parameters of boost to branch CMS
//   PHNUM    4)               1      event number from outside         
//----------------------------------------------------------------------


//----------------------------------------------------------------------
//
//    PHOTOS_MAKE:   General search routine
//
//    Purpose:  Search through the /PH_HEPEVT/ standard HEP common, sta-
//              rting from  the IPPAR-th  particle.  Whenevr  branching 
//              point is found routine PHTYPE(IP) is called.
//              Finally if calls on PHTYPE(IP) modified entries, common
//               /PH_HEPEVT/ is ordered.
//
//    Input Parameter:    IPPAR:  Pointer   to   decaying  particle  in
//                                /PH_HEPEVT/ and the common itself,
//
//    Output Parameters:  Common  /PH_HEPEVT/, either with or without 
//                                new particles added.
//
//    Author(s):  Z. Was, B. van Eijk             Created at:  26/11/89
//                                                Last Update: 30/08/93
//
//----------------------------------------------------------------------

void PHOTOS_MAKE_C(int IPARR){
  static int i=1;
  int IPPAR,I,J,NLAST,MOTHER;

  //--
  PHLUPAB(3);

  //      write(*,*) 'at poczatek'
  //       PHODMP();
  IPPAR=abs(IPARR);
  //--   Store pointers for cascade treatement...
  NLAST=hep.nhep;


  //--
  //--   Check decay multiplicity and minimum of correctness..
  if ((hep.jdahep[IPPAR-i][1-i]==0)||(hep.jmohep[hep.jdahep[IPPAR-i][1-i]-i][1-i]!=IPPAR)) return;

  PHOtoRF();

  //      write(*,*) 'at przygotowany'
  //       PHODMP();

  //--
  //-- single branch mode 
  //-- IPPAR is original position where the program was called

  //-- let-s do generation
  phtype_(&IPPAR);


  //--   rearrange  /PH_HEPEVT/  for added particles.
  //--   at present this may be not needed as information 
  //--   is set at HepMC level.
  if (hep.nhep>NLAST){
    for(I=NLAST+1;I<=hep.nhep;I++){
      //--
      //--   Photon mother and vertex...
      MOTHER=hep.jmohep[I-i][1-i];
      hep.jdahep[MOTHER-i][2-i]=I;
      for( J=1;J<=4;J++){
        hep.vhep[I-i][J-i]=hep.vhep[I-1-i][J-i];
      }
    }
  }
  //      write(*,*) 'at po dzialaniu '
  //      PHODMP();
  PHOtoLAB();
  //      write(*,*) 'at koniec'
  //      PHODMP();
  return;
}




//----------------------------------------------------------------------
//
//    PHCORK: corrects kinmatics of subbranch needed if host program
//            produces events with the shaky momentum conservation
//
//    Input Parameters:   Common /PHOEVT/, MODCOR
//                        MODCOR >0 type of action
//                               =1 no action
//                               =2 corrects energy from mass
//                               =3 corrects mass from energy
//                               =4 corrects energy from mass for 
//                                  particles up to .4 GeV mass, 
//                                  for heavier ones corrects mass,
//                               =5 most complete correct also of mother
//                                  often necessary for exponentiation.
//                               =0 execution mode 
//
//    Output Parameters:  corrected /PHOEVT/
//
//    Author(s):  P.Golonka, Z. Was               Created at:  01/02/99
//						  Modified  :  07/07/13
//----------------------------------------------------------------------

void PHCORK(int MODCOR){
      
  double M,P2,PX,PY,PZ,E,EN,XMS;
  int I,K;
  FILE *PHLUN = stdout;


  static int MODOP=0;
  static int IPRINT=0;
  static double  MCUT=0.4;
  static int i=1;

  if(MODCOR !=0){
    //       INITIALIZATION
    MODOP=MODCOR;

    fprintf(PHLUN,"Message from PHCORK(MODCOR):: initialization\n");
    if(MODOP==1) fprintf(PHLUN,"MODOP=1 -- no corrections on event: DEFAULT\n");
    else if(MODOP==2) fprintf(PHLUN,"MODOP=2 -- corrects Energy from mass\n");
    else if(MODOP==3) fprintf(PHLUN,"MODOP=3 -- corrects mass from Energy\n");
    else if(MODOP==4){
      fprintf(PHLUN,"MODOP=4 -- corrects Energy from mass to Mcut\n");
      fprintf(PHLUN,"           and mass from  energy above  Mcut\n");
      fprintf(PHLUN," Mcut=%6.3f GeV",MCUT);
    }
    else if(MODOP==5) fprintf(PHLUN,"MODOP=5 -- corrects Energy from mass+flow\n");

    else{
      fprintf(PHLUN,"PHCORK wrong MODCOR=%4i\n",MODCOR);
      exit(0);
    }
    return;
  }

  if(MODOP==0&&MODCOR==0){
    fprintf(PHLUN,"PHCORK lack of initialization\n");
    exit(0);
  }

  // execution mode
  // ==============
  // ============== 

     
  PX=0.0;
  PY=0.0;
  PZ=0.0;
  E =0.0;

  if    (MODOP==1){
    //     -----------------------
    //       In this case we do nothing
    return;
  }
  else if(MODOP==2){
    //     -----------------------
    //      lets loop thru all daughters and correct their energies 
    //      according to E^2=p^2+m^2

    for( I=3;I<=pho.nhep;I++){
         
      PX=PX+pho.phep[I-i][1-i];
      PY=PY+pho.phep[I-i][2-i];
      PZ=PZ+pho.phep[I-i][3-i];

      P2=pho.phep[I-i][1-i]*pho.phep[I-i][1-i]+pho.phep[I-i][2-i]*pho.phep[I-i][2-i]+pho.phep[I-i][3-i]*pho.phep[I-i][3-i];

      EN=sqrt( pho.phep[I-i][5-i]*pho.phep[I-i][5-i] + P2);
         
      if (IPRINT==1)fprintf(PHLUN,"CORRECTING ENERGY OF %6i: %14.9f => %14.9f\n",I,pho.phep[I-i][4-i],EN);

      pho.phep[I-i][4-i]=EN;
      E = E+pho.phep[I-i][4-i];

    }
  }

  else if (MODOP==5){
    //     -----------------------
    //C      lets loop thru all daughters and correct their energies 
    //C      according to E^2=p^2+m^2

    for( I=3;I<=pho.nhep;I++){
      PX=PX+pho.phep[I-i][1-i];
      PY=PY+pho.phep[I-i][2-i];
      PZ=PZ+pho.phep[I-i][3-i];

      P2=pho.phep[I-i][1-i]*pho.phep[I-i][1-i]+pho.phep[I-i][2-i]*pho.phep[I-i][2-i]+pho.phep[I-i][3-i]*pho.phep[I-i][3-i];

      EN=sqrt( pho.phep[I-i][5-i]*pho.phep[I-i][5-i] + P2);
         
      if (IPRINT==1)fprintf(PHLUN,"CORRECTING ENERGY OF %6i: %14.9f => %14.9f\n",I,pho.phep[I-i][4-i],EN);

      pho.phep[I-i][4-i]=EN;
      E = E+pho.phep[I-i][4-i];

    }
    for( K=1;K<=4;K++){
      pho.phep[1-i][K-i]=0.0;
      for( I=3;I<=pho.nhep;I++){
        pho.phep[1-i][K-i]=pho.phep[1-i][K-i]+pho.phep[I-i][K-i];
      }
    }
    XMS=sqrt(pho.phep[1-i][4-i]*pho.phep[1-i][4-i]-pho.phep[1-i][3-i]*pho.phep[1-i][3-i]-pho.phep[1-i][2-i]*pho.phep[1-i][2-i]-pho.phep[1-i][1-i]*pho.phep[1-i][1-i]);
    pho.phep[1-i][5-i]=XMS;
  }
  else if(MODOP==3){
    //     -----------------------

    //      lets loop thru all daughters and correct their masses 
    //     according to E^2=p^2+m^2

    for (I=3;I<=pho.nhep;I++){
         
      PX=PX+pho.phep[I-i][1-i];
      PY=PY+pho.phep[I-i][2-i];
      PZ=PZ+pho.phep[I-i][3-i];
      E = E+pho.phep[I-i][4-i];

      P2=pho.phep[I-i][1-i]*pho.phep[I-i][1-i]+pho.phep[I-i][2-i]*pho.phep[I-i][2-i]+pho.phep[I-i][3-i]*pho.phep[I-i][3-i];

      M=sqrt(fabs( pho.phep[I-i][4-i]*pho.phep[I-i][4-i] - P2));

      if (IPRINT==1) fprintf(PHLUN,"CORRECTING MASS OF %6i: %14.9f => %14.9f\n",I,pho.phep[I-i][5-i],M);

      pho.phep[I-i][5-i]=M;

    }
      
  }
 else if(MODOP==4){
   //     -----------------------
            
   //      lets loop thru all daughters and correct their masses 
   //      or energies according to E^2=p^2+m^2

    for (I=3;I<=pho.nhep;I++){
         
      PX=PX+pho.phep[I-i][1-i];
      PY=PY+pho.phep[I-i][2-i];
      PZ=PZ+pho.phep[I-i][3-i];
      P2=pho.phep[I-i][1-i]*pho.phep[I-i][1-i]+pho.phep[I-i][2-i]*pho.phep[I-i][2-i]+pho.phep[I-i][3-i]*pho.phep[I-i][3-i];
      M=sqrt(fabs( pho.phep[I-i][4-i]*pho.phep[I-i][4-i] - P2));


      if(M>MCUT){
	if(IPRINT==1) fprintf(PHLUN,"CORRECTING MASS OF %6i: %14.9f => %14.9f\n",I,pho.phep[I-i][5-i],M);
	pho.phep[I-i][5-i]=M;
	E = E+pho.phep[I-i][4-i];
      }
      else{

      EN=sqrt( pho.phep[I-i][5-i]*pho.phep[I-i][5-i] + P2);
      if(IPRINT==1) fprintf(PHLUN,"CORRECTING ENERGY OF %6i: %14.9f =>% 14.9f\n",I ,pho.phep[I-i][4-i],EN);

      pho.phep[I-i][4-i]=EN;
      E = E+pho.phep[I-i][4-i];
      }
         

    }
 }

  //     -----

  if(IPRINT==1){
    fprintf(PHLUN,"CORRECTING MOTHER");
    fprintf(PHLUN,"PX:%14.9f =>%14.9f",pho.phep[1-i][1-i],PX-pho.phep[2-i][1-i]);
    fprintf(PHLUN,"PY:%14.9f =>%14.9f",pho.phep[1-i][2-i],PY-pho.phep[2-i][2-i]);
    fprintf(PHLUN,"PZ:%14.9f =>%14.9f",pho.phep[1-i][3-i],PZ-pho.phep[2-i][3-i]);
    fprintf(PHLUN," E:%14.9f =>%14.9f",pho.phep[1-i][4-i], E-pho.phep[2-i][4-i]);
  }

  pho.phep[1-i][1-i]=PX-pho.phep[2-i][1-i];
  pho.phep[1-i][2-i]=PY-pho.phep[2-i][2-i];
  pho.phep[1-i][3-i]=PZ-pho.phep[2-i][3-i];
  pho.phep[1-i][4-i]=E -pho.phep[2-i][4-i];


  P2=pho.phep[1-i][1-i]*pho.phep[1-i][1-i]+pho.phep[1-i][2-i]*pho.phep[1-i][2-i]+pho.phep[1-i][3-i]*pho.phep[1-i][3-i];
  if(pho.phep[1-i][4-i]*pho.phep[1-i][4-i]>P2){
    M=sqrt(pho.phep[1-i][4-i]*pho.phep[1-i][4-i] - P2 );
    if(IPRINT==1)fprintf(PHLUN," M: %14.9f => %14.9f\n",pho.phep[1-i][5-i],M);
    pho.phep[1-i][5-i]=M;
  }

  PHLUPA(25);

}






//----------------------------------------------------------------------
//
//    PHOTOS:   PHOton radiation in  decays DOing of KINematics
//
//    Purpose:  Starting  from   the  charged  particle energy/momentum,
//              PNEUTR, photon  energy  fraction and photon  angle  with
//              respect  to  the axis formed by charged particle energy/
//              momentum  vector  and PNEUTR, scale the energy/momentum,
//              keeping the original direction of the neutral system  in
//              the lab. frame untouched.
//
//    Input Parameters:   IP:      Pointer  to   decaying  particle   in
//                                 /PHOEVT/  and   the   common   itself
//                        NCHARB:  pointer to the charged radiating
//                                 daughter in /PHOEVT/.
//                        NEUDAU:  pointer to the first neutral daughter
//    Output Parameters:  Common /PHOEVT/, with photon added.
//
//    Author(s):  Z. Was, B. van Eijk             Created at:  26/11/89
//                                                Last Update: 27/05/93
//
//----------------------------------------------------------------------

void PHODO(int IP,int NCHARB,int NEUDAU){
  static int i=1;
  double QNEW,QOLD,EPHOTO,PMAVIR;
  double GNEUT,DATA;
  double PHORANC,CCOSTH,SSINTH,PVEC[4],PARNE;
  double TH3,FI3,TH4,FI4,FI5,ANGLE;
  int I,J,FIRST,LAST;

  //--
  EPHOTO=phophs_.xphoto*pho.phep[IP-i][5-i]/2.0;
  PMAVIR=sqrt(pho.phep[IP-i][5-i]*(pho.phep[IP-i][5-i]-2.0*EPHOTO));
  //--
  //--   Reconstruct  kinematics  of  charged particle  and  neutral system
  phorest_.fi1=PHOAN1(phomom_.pneutr[1-i],phomom_.pneutr[2-i]);
  //--
  //--   Choose axis along  z of  PNEUTR, calculate  angle  between x and y
  //--   components  and z  and x-y plane and  perform Lorentz transform...
  phorest_.th1=PHOAN2(phomom_.pneutr[3-i],sqrt(phomom_.pneutr[1-i]*phomom_.pneutr[1-i]+phomom_.pneutr[2-i]*phomom_.pneutr[2-i]));
  PHORO3(-phorest_.fi1,phomom_.pneutr);
  PHORO2(-phorest_.th1,phomom_.pneutr);
  //--
  //--   Take  away  photon energy from charged particle and PNEUTR !  Thus
  //--   the onshell charged particle  decays into virtual charged particle
  //--   and photon.  The virtual charged  particle mass becomes:
  //--   SQRT(pho.phep[5,IP)*(pho.phep[5,IP)-2*EPHOTO)).  Construct  new PNEUTR mo-
  //--   mentum in the rest frame of the parent:
  //--   1) Scaling parameters...
  QNEW=PHOTRI(PMAVIR,phomom_.pneutr[5-i],pho.phep[NCHARB-i][5-i]);
  QOLD=phomom_.pneutr[3-i];
  GNEUT=(QNEW*QNEW+QOLD*QOLD+phomom_.mnesqr)/(QNEW*QOLD+sqrt((QNEW*QNEW+phomom_.mnesqr)*(QOLD*QOLD+phomom_.mnesqr)));
  if(GNEUT<1.0){
    DATA=0.0;
    PHOERR(4,"PHOKIN",DATA);
  }
  PARNE=GNEUT-sqrt(max(GNEUT*GNEUT-1.0,0.0));
  //--
  //--   2) ...reductive boost...
  PHOBO3(PARNE,phomom_.pneutr);
  //--
  //--   ...calculate photon energy in the reduced system...
  pho.nhep=pho.nhep+1;
  pho.isthep[pho.nhep-i]=1;
  pho.idhep[pho.nhep-i] =22;
  //--   Photon mother and daughter pointers !
  pho.jmohep[pho.nhep-i][1]=IP;
  pho.jmohep[pho.nhep-i][2]=0;
  pho.jdahep[pho.nhep-i][1]=0;
  pho.jdahep[pho.nhep-i][2]=0;
  pho.phep[pho.nhep-i][4-i]=EPHOTO*pho.phep[IP-i][5-i]/PMAVIR;
  //--
  //--   ...and photon momenta
  CCOSTH=-phophs_.costhg;
  SSINTH=phophs_.sinthg;
  TH3=PHOAN2(CCOSTH,SSINTH);
  FI3=phpico_.twopi*Photos::randomDouble();
  pho.phep[pho.nhep-i][1-i]=pho.phep[pho.nhep-i][4-i]*phophs_.sinthg*cos(FI3);
  pho.phep[pho.nhep-i][2-i]=pho.phep[pho.nhep-i][4-i]*phophs_.sinthg*sin(FI3);
  //--
  //--   Minus sign because axis opposite direction of charged particle !
  pho.phep[pho.nhep-i][3-i]=-pho.phep[pho.nhep-i][4-i]*phophs_.costhg;
  pho.phep[pho.nhep-i][5-i]=0.0;
  //--
  //--   Rotate in order to get photon along z-axis
  PHORO3(-FI3,phomom_.pneutr);
  PHORO3(-FI3,pho.phep[pho.nhep-i]);
  PHORO2(-TH3,phomom_.pneutr);
  PHORO2(-TH3,pho.phep[pho.nhep-i]);
  ANGLE=EPHOTO/pho.phep[pho.nhep-i][4-i];
  //--
  //--   Boost to the rest frame of decaying particle
  PHOBO3(ANGLE,phomom_.pneutr);
  PHOBO3(ANGLE,pho.phep[pho.nhep-i]);
  //--
  //--   Back in the parent rest frame but PNEUTR not yet oriented !
  FI4=PHOAN1(phomom_.pneutr[1-i],phomom_.pneutr[2-i]);
  TH4=PHOAN2(phomom_.pneutr[3-i],sqrt(phomom_.pneutr[1-i]*phomom_.pneutr[1-i]+phomom_.pneutr[2-i]*phomom_.pneutr[2-i]));
  PHORO3(FI4,phomom_.pneutr);
  PHORO3(FI4,pho.phep[pho.nhep-i]);
  //--
  for(I=2; I<=4;I++) PVEC[I-i]=0.0;
  PVEC[1-i]=1.0;

  PHORO3(-FI3,PVEC);
  PHORO2(-TH3,PVEC);
  PHOBO3(ANGLE,PVEC);
  PHORO3(FI4,PVEC);
  PHORO2(-TH4,phomom_.pneutr);
  PHORO2(-TH4,pho.phep[pho.nhep-i]);
  PHORO2(-TH4,PVEC);
  FI5=PHOAN1(PVEC[1-i],PVEC[2-i]);
  //--
  //--   Charged particle restores original direction
  PHORO3(-FI5,phomom_.pneutr);
  PHORO3(-FI5,pho.phep[pho.nhep-i]);
  PHORO2(phorest_.th1,phomom_.pneutr);
  PHORO2(phorest_.th1,pho.phep[pho.nhep-i]);
  PHORO3(phorest_.fi1,phomom_.pneutr);
  PHORO3(phorest_.fi1,pho.phep[pho.nhep-i]);
  //--   See whether neutral system has multiplicity larger than 1...

  if((pho.jdahep[IP-i][2-i]-pho.jdahep[IP-i][1-i])>1){
    //--   Find pointers to components of 'neutral' system
    //--
    FIRST=NEUDAU;
    LAST=pho.jdahep[IP-i][2-i];
    for(I=FIRST;I<=LAST;I++){
      if(I!=NCHARB && ( pho.jmohep[I-i][1-i]==IP)){
	//--
	//--   Reconstruct kinematics...
	PHORO3(-phorest_.fi1,pho.phep[I-i]);
	PHORO2(-phorest_.th1,pho.phep[I-i]);
	//--
	//--   ...reductive boost
	PHOBO3(PARNE,pho.phep[I-i]);
	//--
	//--   Rotate in order to get photon along z-axis
	PHORO3(-FI3,pho.phep[I-i]);
	PHORO2(-TH3,pho.phep[I-i]);
	//--
	//--   Boost to the rest frame of decaying particle
	PHOBO3(ANGLE,pho.phep[I-i]);
	//--
	//--   Back in the parent rest-frame but PNEUTR not yet oriented.
	PHORO3(FI4,pho.phep[I-i]);
	PHORO2(-TH4,pho.phep[I-i]);
	//--
	//--   Charged particle restores original direction
	PHORO3(-FI5,pho.phep[I-i]);
	PHORO2(phorest_.th1,pho.phep[I-i]);
	PHORO3(phorest_.fi1,pho.phep[I-i]);
      }
    }
  }
  else{
    //--
    //   ...only one 'neutral' particle in addition to photon!
    for(J=1;J<=4;J++) pho.phep[NEUDAU-i][J-i]=phomom_.pneutr[J-i];
  }
  //--
  //--   All 'neutrals' treated, fill /PHOEVT/ for charged particle...
  for (J=1;J<=3;J++) pho.phep[NCHARB-i][J-i]=-(pho.phep[pho.nhep-i][J-i]+phomom_.pneutr[J-i]);
                     pho.phep[NCHARB-i][4-i]=pho.phep[IP-i][5-i]-(pho.phep[pho.nhep-i][4-i]+phomom_.pneutr[4-i]);
  //--
}

