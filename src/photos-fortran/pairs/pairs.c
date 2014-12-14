#include <cmath>
#include <stdio.h>
#include <stdlib.h>  
extern "C" double angfi_(double *pX,double *pY){                                               
  const double &X = *pX;
  const double &Y = *pY;
  double THE;
  const double PI=3.141592653589793238462643;

  if(X==0.0 && Y==0.0) return 0.0;           
   if(fabs(Y) < fabs(X)){                                         
    THE=atan(fabs(Y/X));                                              
    if(X<=0.0) THE=PI-THE;}                                         
  else{                                                              
    THE=acos(X/sqrt(X*X +Y*Y));                                     
  }
   if(Y<0.0) THE=2*PI-THE;
   return THE;
}

extern "C" double angxy_(double *pX,double *pY){                                               
  const double &X = *pX;
  const double &Y = *pY;
  double THE;
  const double PI=3.141592653589793238462643;

  if(X==0.0 && Y==0.0) return 0.0;

  if(fabs(Y) < fabs(X)){                                         
    THE=atan(fabs(Y/X));                                              
    if(X<=0.0) THE=PI-THE;}                                         
  else{                                                              
    THE=acos(X/sqrt(X*X +Y*Y));                                     
  }                                       
  return THE;                                                         
}

extern "C" void  bostd3_(double *pEXE,double PVEC[4],double QVEC[4]){
  // ----------------------------------------------------------------------
  // BOOST ALONG Z AXIS, EXE=EXP(ETA), ETA= HIPERBOLIC VELOCITY.
  //
  //     USED BY : KORALZ RADKOR
  /// ----------------------------------------------------------------------
  const double EXE = *pEXE;
  int j=1;  // convention of indices of Riemann space must be preserved.
  double RPL,RMI,QPL,QMI;
  double RVEC[4];


  RVEC[1-j]=PVEC[1-j];
  RVEC[2-j]=PVEC[2-j];
  RVEC[3-j]=PVEC[3-j];
  RVEC[4-j]=PVEC[4-j];
  RPL=RVEC[4-j]+RVEC[3-j];
  RMI=RVEC[4-j]-RVEC[3-j];
  QPL=RPL*EXE;
  QMI=RMI/EXE;
  QVEC[1-j]=RVEC[1-j];
  QVEC[2-j]=RVEC[2-j];
  QVEC[3-j]=(QPL-QMI)/2;
  QVEC[4-j]=(QPL+QMI)/2;
}

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


extern "C" void   lortra_(int *pKEY,double *pPRM,double PNEUTR[4],double PNU[4],double PAA[4],double PP[4],double PE[4]){
  // --------------------------------------------------------------------- 
  // THIS ROUTINE PERFORMS LORENTZ TRANSFORMATION ON MANY 4-VECTORS        
  // KEY   =1    BOOST    ALONG   3RD AXIS                                 
  //       =2    ROTATION AROUND 2ND AXIS                                  
  //       =3    ROTATION AROUND 3RD AXIS                                 
  // PRM         TRANSFORMATION PARAMETER - ANGLE OR EXP(HIPERANGLE).     
  //                                                                      
  //    called by : RADCOR                                                
  // --------------------------------------------------------------------- 
  const int &KEY = *pKEY;
  double PRM = *pPRM;                                                                       
  if(KEY==1){                                           
    bostd3_(&PRM,PNEUTR,PNEUTR);                                 
    bostd3_(&PRM,PNU ,PNU );                                     
    bostd3_(&PRM,PAA ,PAA );                                     
    bostd3_(&PRM,PE ,PE );                                     
    bostd3_(&PRM,PP ,PP );
  }                                     
  else if (KEY==2){                                           
    rotod2_(&PRM,PNEUTR,PNEUTR);                                 
    rotod2_(&PRM,PNU ,PNU );                                     
    rotod2_(&PRM,PAA ,PAA );                                    
    rotod2_(&PRM,PE  ,PE  );                                     
    rotod2_(&PRM,PP  ,PP );
  }                                    
  else if (KEY==3){                                          
    rotod3_(&PRM,PNEUTR,PNEUTR);                                
    rotod3_(&PRM,PNU ,PNU );                                     
    rotod3_(&PRM,PAA ,PAA );                                    
    rotod3_(&PRM,PE  ,PE  );                                     
    rotod3_(&PRM,PP  ,PP );
  }                                    
  else{                                                            
    printf (" STOP IN LOTRA. WRONG KEYTRA"); 
    exit(-1);
  }
}

extern "C" double amast_(double VEC[4]){
  int i=1;  // convention of indices of Riemann space must be preserved
  double ama= VEC[4-i]*VEC[4-i]-VEC[1-i]*VEC[1-i]-VEC[2-i]*VEC[2-i]-VEC[3-i]*VEC[3-i];
  ama=sqrt(fabs(ama));
  return ama;
} 

extern "C" void spaj_(int *pKUDA,double P2[4],double Q2[4],double PP[4],double PE[4]){    
  //     **********************     
  // THIS PRINTS OUT FOUR MOMENTA OF PHOTONS 
  // ON OUTPUT UNIT NOUT
  const int &KUDA = *pKUDA;

  double SUM[4];
  const int KLUCZ=0;
  if (KLUCZ==0) return;

  printf (" %10i =====================SPAJ==================== \n", KUDA);
  printf (" P2 %18.13f %18.13f %18.13f %18.13f \n",P2[0],P2[1],P2[2],P2[3]);
  printf (" Q2 %18.13f %18.13f %18.13f %18.13f \n",Q2[0],Q2[1],Q2[2],Q2[3]);
  printf (" PE %18.13f %18.13f %18.13f %18.13f \n",PE[0],PE[1],PE[2],PE[3]);
  printf (" PP %18.13f %18.13f %18.13f %18.13f \n",PP[0],PP[1],PP[2],PP[3]);
 
  for( int k=0;k<=3;k++) SUM[k]=P2[k]+Q2[k]+PE[k]+PP[k];

  printf ("SUM %18.13f %18.13f %18.13f %18.13f \n",SUM[0],SUM[1],SUM[2],SUM[3]);
}
extern "C" {
 extern struct { 
   double fi0; // FI0
   double fi1; // FI1
   double fi2; // FI2
   double fi3; // FI3
   double fi4; // FI4
   double fi5; // FI5
   double th0; // TH0
   double th1; // TH1
   double th3; // TH3
   double th4; // TH4
   double parneu; // PARNEU
   double parch; // PARCH
   double bpar; // BPAR
   double bsta; // BSTA
   double bstb; // BSTB
  } parkin_;
}

extern "C" void partra_(int *pIBRAN,double PHOT[4]){

   const int &IBRAN = *pIBRAN;

   // particularily stupid way of introdicing minus for argument of the function.
   double FI0 = -parkin_.fi0;
   double FI1 = -parkin_.fi1;
   double FI2 = -parkin_.fi2;
   double FI3 = -parkin_.fi3;
   double FI5 = -parkin_.fi5;
   double TH0 = -parkin_.th0;
   double TH1 = -parkin_.th1;
   double TH3 = -parkin_.th3;
   double TH4 = -parkin_.th4;

   rotod3_(        &FI0,PHOT,PHOT); 
   rotod2_(        &TH0,PHOT,PHOT);
   bostd3_(&parkin_.bsta,PHOT,PHOT);
   rotod3_(        &FI1,PHOT,PHOT);
   rotod2_(        &TH1,PHOT,PHOT);
   rotod3_(        &FI2,PHOT,PHOT);

   if(IBRAN==-1){
     bostd3_(&parkin_.parneu,PHOT,PHOT);
   }
   else{
     bostd3_(&parkin_.parch,PHOT,PHOT);
   }

   rotod3_(        &FI3,PHOT,PHOT);
   rotod2_(        &TH3,PHOT,PHOT);
   bostd3_(&parkin_.bpar,PHOT,PHOT);
   rotod3_(&parkin_.fi4,PHOT,PHOT);
   rotod2_(        &TH4,PHOT,PHOT);
   rotod3_(        &FI5,PHOT,PHOT);
   rotod3_(&parkin_.fi2,PHOT,PHOT);
   rotod2_(&parkin_.th1,PHOT,PHOT); 
   rotod3_(&parkin_.fi1,PHOT,PHOT);
   bostd3_(&parkin_.bstb,PHOT,PHOT);
   rotod2_(&parkin_.th0,PHOT,PHOT);
   rotod3_(&parkin_.fi0,PHOT,PHOT);
   
}


