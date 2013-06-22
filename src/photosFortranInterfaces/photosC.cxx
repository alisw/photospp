#include "Photos.h"
#include <cmath>
#include <iostream>
#include "f_Init.h"
#include "PH_HEPEVT_Interface.h"
using std::cout;
using std::endl;
using std::max;
using namespace Photospp;

// Declaration of structs defined in f_Init.h
struct PHOSTA phosta_;
struct PHLUPY phlupy_;

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
  for(I=0;I<ph_hepevt_.nhep;I++){ 
    //--
    //--   For 'stable particle' calculate vector momentum sum
    if (ph_hepevt_.jdahep[1-i][I-i]==0){
      for(J=1; J<=4;J++){
	SUMVEC[J-i]=SUMVEC[J-i]+ph_hepevt_.phep[J-i][I-i];
      }
      if (ph_hepevt_.jmohep[2-i][I-i]==0){
	fprintf(PHLUN,"%4i %7i %s %4i %s Stable %9.2f %9.2f %9.2f %9.2f %9.2f\n" ,  I,ph_hepevt_.idhep[I-i],X3,ph_hepevt_.jmohep[1-i][I-i],X7,ph_hepevt_.phep[1-i][I-i],ph_hepevt_.phep[2-i][I-i],ph_hepevt_.phep[3-i][I-i],ph_hepevt_.phep[4-i][I-i],ph_hepevt_.phep[5-i][I-i]);
      }
      else{
	fprintf(PHLUN,"%4i %7i %4i - %4i %s Stable %9.2f %9.2f %9.2f %9.2f %9.2f\n",I,ph_hepevt_.idhep[I-i],ph_hepevt_.jmohep[1-i][I-i],ph_hepevt_.jmohep[2-i][I-i], X4,ph_hepevt_.phep[1-i][I-i],ph_hepevt_.phep[2-i][I-i],ph_hepevt_.phep[3-i][I-i],ph_hepevt_.phep[4-i][I-i],ph_hepevt_.phep[5-i][I-i]);
      }
    }
    else{
      if(ph_hepevt_.jmohep[2-i][I-i]==0){
	fprintf(PHLUN,"%4i %7i %s %4i %s %4i - %4i %9.2f %9.2f %9.2f %9.2f %9.2f\n" ,  I,ph_hepevt_.idhep[I-i],X3,ph_hepevt_.jmohep[1-i][I-i], X6,ph_hepevt_.jdahep[1-i][I-i],ph_hepevt_.jdahep[2-i][I-i],ph_hepevt_.phep[1-i][I-i],ph_hepevt_.phep[2-i][I-i],ph_hepevt_.phep[3-i][I-i],ph_hepevt_.phep[4-i][I-i],ph_hepevt_.phep[5-i][I-i]);
      }
      else{
	fprintf(PHLUN,"%4i %7i %4i - %4i %4i - %4i %9.2f %9.2f %9.2f %9.2f %9.2f\n",  I,ph_hepevt_.idhep[I-i],ph_hepevt_.jmohep[1-i][I-i],ph_hepevt_.jmohep[2-i][I-i],ph_hepevt_.jdahep[1-i][I-i],ph_hepevt_.jdahep[2-i][I-i],ph_hepevt_.phep[1-i][I-i],ph_hepevt_.phep[2-i][I-i],ph_hepevt_.phep[3-i][I-i],ph_hepevt_.phep[4-i][I-i],ph_hepevt_.phep[5-i][I-i]);
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

  int I,J;
  int &IPOIN = phlupy_.ipoin;
  int &IPOINM= phlupy_.ipoinm;
  static int IPOIN0=-5;
  static int i=1;
  int IOUT;
  double  SUM[5];
  FILE *PHLUN = stdout;

  if (IPOIN0<0){
    IPOIN0=400000; //  ! maximal no-print point
    IPOIN =IPOIN0;
    IPOINM=400001; // ! minimal no-print point
  }
  if (IPOINT<=IPOINM||IPOINT>=IPOIN ) return;
  IOUT=56;
  if ((int)phnum_.iev<1000){
    for(I=1; I<=5;I++) SUM[I-i]=0.0;
     
    fprintf(PHLUN,"EVENT NR= %i WE ARE TESTING /PH_HEPEVT/ at IPOINT=%i \n",(int)phnum_.iev,IPOINT);
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

  int I,J;
  int &IPOIN = phlupy_.ipoin;
  int &IPOINM= phlupy_.ipoinm;
  static int IPOIN0=-5;
  static int i=1;
  int IOUT;
  double  SUM[5];
  FILE *PHLUN = stdout;

  if (IPOIN0<0){
    IPOIN0=400000; //  ! maximal no-print point
    IPOIN =IPOIN0;
    IPOINM=400001; // ! minimal no-print point
  }
  
  if (IPOINT<=IPOINM||IPOINT>=IPOIN ) return;
  IOUT=56;
  if ((int)phnum_.iev<1000){
    for(I=1; I<=5;I++) SUM[I-i]=0.0;
     
    fprintf(PHLUN,"EVENT NR= %i WE ARE TESTING /PHOEVT/ at IPOINT=%i \n",(int)phnum_.iev,IPOINT);
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











