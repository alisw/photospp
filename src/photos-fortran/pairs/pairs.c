#include <cmath>
#include <stdio.h>
#include <stdlib.h>  
//
  inline double xlam(double A,double B,double C){return sqrt((A-B-C)*(A-B-C)-4.0*B*C);}

  inline double max(double a, double b)
  {
  return (a > b) ? a : b;
  }
   // 
extern "C" void varran_( double RRR[], int *N);

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
extern "C" double angfi(double X,double Y){                                               
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

extern "C" double angxy(double X,double Y){                                               
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
extern "C" void  bostd3(double EXE,double PVEC[4],double QVEC[4]){
  // ----------------------------------------------------------------------
  // BOOST ALONG Z AXIS, EXE=EXP(ETA), ETA= HIPERBOLIC VELOCITY.
  //
  //     USED BY : KORALZ RADKOR
  /// ----------------------------------------------------------------------
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


extern "C" void rotod3(double ANGLE,double PVEC[4],double QVEC[4]){

 
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



extern "C" void   rotod2(double PHI,double PVEC[4],double QVEC[4]){

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

extern "C" void   lortra(int KEY,double PRM,double PNEUTR[4],double PNU[4],double PAA[4],double PP[4],double PE[4]){
  // --------------------------------------------------------------------- 
  // THIS ROUTINE PERFORMS LORENTZ TRANSFORMATION ON MANY 4-VECTORS        
  // KEY   =1    BOOST    ALONG   3RD AXIS                                 
  //       =2    ROTATION AROUND 2ND AXIS                                  
  //       =3    ROTATION AROUND 3RD AXIS                                 
  // PRM         TRANSFORMATION PARAMETER - ANGLE OR EXP(HIPERANGLE).     
  //                                                                      
  //    called by : RADCOR                                                
  // --------------------------------------------------------------------- 
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


extern "C" void spaj(int KUDA,double P2[4],double Q2[4],double PP[4],double PE[4]){    
  //     **********************     
  // THIS PRINTS OUT FOUR MOMENTA OF PHOTONS 
  // ON OUTPUT UNIT NOUT

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


extern "C" void trypar_(bool *pJESLI,double *pSTRENG,double PA[4],double PB[4],double PE[4],double PP[4]){       
  bool &JESLI = *pJESLI;
  const double &STRENG = *pSTRENG;
  //      COMMON  /PARKIN/ 
  double &FI0=parkin_.fi0;
  double &FI1=parkin_.fi1;
  double &FI2=parkin_.fi2;
  double &FI3=parkin_.fi3;
  double &FI4=parkin_.fi4;
  double &FI5=parkin_.fi5;
  double &TH0=parkin_.th0;
  double &TH1=parkin_.th1;
  double &TH3=parkin_.th3;
  double &TH4=parkin_.th4;
  double &PARNEU=parkin_.parneu;
  double &PARCH=parkin_.parch;
  double &BPAR=parkin_.bpar;
  double &BSTA=parkin_.bsta;
  double &BSTB=parkin_.bstb;

  double  PNEUTR[4],PAA[4],PHOT[4],PSUM[4];
  double VEC[4];                                                
  double RRR[8]; 
  bool JESLIK; 
  const double PI=3.141592653589793238462643;     
  const double AMEL=.511e-3;
  const double XK0 =  1.0e-3;                                 
  const double ALFINV= 137.01;
  const int j=1;  // convention of indices of Riemann space must be preserved.
   
  PA[4-j]=max(PA[4-j],sqrt(PA[1-j]*PA[1-j]+PA[2-j]*PA[2-j]+PA[3-j]*PA[3-j]));
  PB[4-j]=max(PB[4-j],sqrt(PB[1-j]*PB[1-j]+PB[2-j]*PB[2-j]+PB[3-j]*PB[3-j]));
  // 4-MOMENTUM OF THE NEUTRAL SYSTEM                                 

  for( int k=0;k<=3;k++){
    PE[k]    =0.0;
    PP[k]    =0.0;
    PSUM[k]  =PA[k]+PB[k];
    PAA[k]   =PA[k];
    PNEUTR[k]=PB[k];
  }
  if((PAA[4-j]+PNEUTR[4-j])< 0.01 ){
    JESLI=false;                                      
    return;
  }

  // MASSES OF THE NEUTRAL AND CHARGED SYSTEMS AND OVERALL MASS
  // FIRST WE HAVE TO GO TO THE RESTFRAME TO GET RID OF INSTABILITIES 
  // FROM BHLUMI OR ANYTHING ELSE            
  // THIRD AXIS ALONG PNEUTR                                             
  double X1 = PSUM[1-j];                                                  
  double X2 = PSUM[2-j]; 
  FI0  =angfi(X1,X2);           
  X1 = PSUM[3-j];                                                 
  X2 = sqrt(PSUM[1-j]*PSUM[1-j]+PSUM[2-j]*PSUM[2-j]);                            
  TH0  =angxy(X1,X2) ;
  spaj(-2,PNEUTR,PAA,PP,PE);
  lortra(3,-FI0,PNEUTR,VEC,PAA,PP,PE);
  lortra(2,-TH0,PNEUTR,VEC,PAA,PP,PE);
  rotod3(-FI0,PSUM,PSUM);
  rotod2(-TH0,PSUM,PSUM);
  BSTA=(PSUM[4-j]-PSUM[3-j])/sqrt(PSUM[4-j]*PSUM[4-j]-PSUM[3-j]*PSUM[3-j]);
  BSTB=(PSUM[4-j]+PSUM[3-j])/sqrt(PSUM[4-j]*PSUM[4-j]-PSUM[3-j]*PSUM[3-j]);
  lortra(1,BSTA,PNEUTR,VEC,PAA,PP,PE);
  spaj(-1,PNEUTR,PAA,PP,PE);                                  
  double AMNE=amast_(PNEUTR);                                              
  double AMCH=amast_(PAA);                                                
  if(AMCH<0.0) AMCH=AMEL;                                   
  if (AMNE<0.0) AMNE=0.0;
  double AMTO =PAA[4-j]+PNEUTR[4-j];
  int osm=8;
  varran_(RRR,&osm);

  double PRHARD;
  PRHARD= (1.0/PI/ALFINV)*(1.0/PI/ALFINV)* (2.0*log(AMTO/AMEL/2.0)) * 
          log(1.0/XK0) * log(AMTO*AMTO/2.0/AMEL/AMEL);

  // this just enforces hard pairs to be generated 'always'
  // this is for the sake of tests only.
  //  PRHARD=0.99;  
  //

  double XMP=2.0*AMEL*exp(RRR[1-j]*log(AMTO/2.0/AMEL)); 
  double XP =AMTO*XK0*exp(RRR[2-j]*log(1.0/XK0));    
  double C1 =1.0-4.0*AMEL*AMEL/AMTO/AMTO*exp(RRR[3-j]*log(AMTO*AMTO/2.0/AMEL/AMEL));
  double FIX1=2.0*PI*RRR[4-j];
  double C2  =1.0-2.0*RRR[5-j]; 
  double FIX2=2.0*PI*RRR[6-j];
 
  JESLI=(RRR[7-j]<PRHARD)       &&  
        (XMP<(AMTO-AMNE-AMCH))  &&  
        (XP >XMP)               &&  
        (XP <((AMTO*AMTO+XMP*XMP-(AMCH+AMNE)*(AMCH+AMNE))/2.0/AMTO));

  // histograming .......................
  JESLIK=     (XP <((AMTO*AMTO+XMP*XMP-(AMCH+AMNE)*(AMCH+AMNE))/2.0/AMTO));
  double WTA=0.0;
  WTA=WTA*5;
  if(JESLIK) WTA=1.0;
  //      GMONIT( 0,101   ,WTA,1D0,0D0)
  JESLIK= (XMP<(AMTO-AMNE-AMCH));
  WTA=0.0;
  if(JESLIK) WTA=1.0;
  //      GMONIT( 0,102   ,WTA,1D0,0D0)
  JESLIK= (XMP<(AMTO-AMNE-AMCH))&& (XP >XMP);

  WTA=0.0;
  if(JESLIK) WTA=1.0;
  //      GMONIT( 0,103   ,WTA,1D0,0D0)
  JESLIK=
         (XMP<(AMTO-AMNE-AMCH))&&
         (XP >XMP)             &&
         (XP <((AMTO*AMTO+XMP*XMP-(AMCH+AMNE)*(AMCH+AMNE))/2.0/AMTO));
  WTA=0.0;
  if (JESLIK) WTA=1.0;
  //      GMONIT( 0,104   ,WTA,1D0,0D0)
  // end of histograming ................  

  if (!JESLI) return;

  // ... jacobians weights etc. 

  double XMK2=AMTO*AMTO+XMP*XMP-2.0*XP*AMTO;
  double XPMAX=(AMTO*AMTO+XMP*XMP-(AMCH+AMNE)*(AMCH+AMNE))/2.0/AMTO;
  double YOT3=(1-C1+4.0*AMEL*AMEL/AMTO/AMTO)/(1-C1+XMP*XMP/AMTO/AMTO);
  double YOT2=(1-C1)*(1+C1)/2.0/(1-C1+XMP*XMP/AMTO/AMTO)*
               xlam(AMTO*AMTO,XMK2,XMP*XMP)/(AMTO*AMTO+XMP*XMP-XMK2);

  double YOT1=xlam(1.0,AMEL*AMEL/XMP/XMP, AMEL*AMEL/XMP/XMP)*
              ((1-C1+XMP*XMP/AMTO/AMTO)/(1-C1+XMP*XMP/XP/XP))*
              ((1-C1+XMP*XMP/AMTO/AMTO)/(1-C1+XMP*XMP/XP/XP))*
              (1-XP/XPMAX+0.5*(XP/XPMAX)*(XP/XPMAX))*2.0/3.0;

// note that the factor 2/3 in YOT1 above should be replaced by the 
// appropriate A-P kernel for gamma splitting to e+e- !!!!!!!
// the part of the weight below, should have average 1, but fluctuates 
// wildly. This cancelation is important ingredient of the leading logs.
//      YOT1=YOT1*
//     $     (1D0-XP/(0.5D0*AMTO))/(1D0-XP/XPMAX)*
//     $     XLAM(1D0,AMCH**2/XMK2,AMNEU**2/XMK2)/
//     $     XLAM(1D0,AMCH**2/AMTO**2,AMNEU**2/AMTO**2)

  double WT=YOT1*YOT2*YOT3;
  //C histograming .......................
  //      GMONIT( 0,105   ,WT  ,1D0,0D0) 
  //      GMONIT( 0,106   ,YOT1,1D0,0D0) 
  //      GMONIT( 0,107   ,YOT2,1D0,0D0) 
  //      GMONIT( 0,108   ,YOT3,1D0,0D0)
  //      GMONIT( 0,109   ,YOT4,1D0,0D0)
  // end of histograming ................ 

  if (RRR[8-j]>WT){
    JESLI=false;
    return;
  }


  //                                                                     
  //                                                                     
  // FRAGMENTATION COMES !!                                              
  //                                                                     
  // THIRD AXIS ALONG PNEUTR                                             
  X1 = PNEUTR[1-j];                                                  
  X2 = PNEUTR[2-j];                                                 
  FI1  =angfi(X1,X2);                                             
  X1 = PNEUTR[3-j];                                                 
  X2 = sqrt(PNEUTR[1-j]*PNEUTR[1-j]+PNEUTR[2-j]*PNEUTR[2-j]) ;                           
  TH1  =angxy(X1,X2); 
  spaj(0,PNEUTR,PAA,PP,PE);                                             
  lortra(3,-FI1,PNEUTR,VEC,PAA,PP,PE);
  lortra(2,-TH1,PNEUTR,VEC,PAA,PP,PE);
  VEC[4-j]=0.0;                                                  
  VEC[3-j]=0.0;                                                 
  VEC[2-j]=0.0;                                                 
  VEC[1-j]=1.0;                                                  
  FI2=angfi(VEC[1-j],VEC[2-j]);                                             
  lortra(3,-FI2,PNEUTR,VEC,PAA,PP,PE);
  spaj(1,PNEUTR,PAA,PP,PE);
                                        
  // STEALING FROM PAA AND PNEUTR ENERGY FOR THE pair
  // ====================================================  
  // NEW MOMENTUM OF PAA AND PNEUTR (IN THEIR VERY REST FRAME)
  // 1) PARAMETERS..... 
  double AMCH2=AMCH*AMCH;
  double AMNE2=AMNE*AMNE;
  double AMTOST=XMK2;
  double QNEW=xlam(AMTOST,AMNE2,AMCH2)/sqrt(AMTOST)/2.0;
  double QOLD=PNEUTR[3-j];
  double GCHAR=(QNEW*QNEW+QOLD*QOLD+AMCH*AMCH)/
               (QNEW*QOLD+sqrt((QNEW*QNEW+AMCH*AMCH)*(QOLD*QOLD+AMCH*AMCH)));
  double GNEU= (QNEW*QNEW+QOLD*QOLD+AMNE*AMNE)/
               (QNEW*QOLD+sqrt((QNEW*QNEW+AMNE*AMNE)*(QOLD*QOLD+AMNE*AMNE)));

  //      GCHAR=(QOLD**2-QNEW**2)/(
  //     &       QNEW*SQRT(QOLD**2+AMCH2)+QOLD*SQRT(QNEW**2+AMCH2)
  //     &                        )
  //      GCHAR=SQRT(1D0+GCHAR**2) 
  //      GNEU=(QOLD**2-QNEW**2)/(
  //     &       QNEW*SQRT(QOLD**2+AMNE2)+QOLD*SQRT(QNEW**2+AMNE2)
  //     &                        )
  //      GNEU=SQRT(1D0+GNEU**2) 
  if(GNEU<1.||GCHAR<1.){
    printf(" TRYPAR GBOOST LT 1., LIMIT OF PHASE SPACE %18.13f %18.13f %18.13f %18.13f %18.13f %18.13f %18.13f %18.13f \n" 
             ,GNEU,GCHAR,QNEW,QOLD,AMTO,AMTOST,AMNE,AMCH);
    double XK=STRENG,XKM=0,XK0DEC=0,AXK=0;
    printf(" %18.13f %18.13f %18.13f %18.13f ",XK,XKM,XK0DEC,AXK);
    return;
  }
  PARCH =GCHAR+sqrt(GCHAR*GCHAR-1.0);
  PARNEU=GNEU -sqrt(GNEU*GNEU -1.0);

  // 2) REDUCTIEV BOOSTS
  bostd3(PARNEU,VEC ,VEC );
  bostd3(PARNEU,PNEUTR,PNEUTR);
  bostd3(PARCH,PAA ,PAA );
  spaj(2,PNEUTR,PAA,PP,PE);
                                             
  // TIME FOR THE PHOTON that is electron pair
  double PMOD=xlam(XMP*XMP,AMEL*AMEL,AMEL*AMEL)/XMP/2.0;
  double S2=sqrt(1.0-C2*C2);
  PP[4-j]=XMP/2.0;
  PP[3-j]=PMOD*C2;
  PP[2-j]=PMOD*S2*cos(FIX2);
  PP[1-j]=PMOD*S2*sin(FIX2);
  PE[4-j]= PP[4-j];
  PE[3-j]=-PP[3-j];
  PE[2-j]=-PP[2-j];
  PE[1-j]=-PP[1-j];
  // PHOTON ENERGY and momentum IN THE REDUCED SYSTEM
  double PENE=(AMTO*AMTO-XMP*XMP-XMK2)/2.0/sqrt(XMK2);
  double PPED=sqrt(PENE*PENE-XMP*XMP);
  FI3=FIX1;
  double COSTHG=C1;
  double SINTHG=sqrt(1.0-C1*C1);
  X1 = -COSTHG;
  X2 =  SINTHG;
  TH3  =angxy(X1,X2);
  PHOT[1-j]=PMOD*SINTHG*cos(FI3); 
  PHOT[2-j]=PMOD*SINTHG*sin(FI3);
  // MINUS BECAUSE AXIS OPPOSITE TO PAA
  PHOT[3-j]=-PMOD*COSTHG;
  PHOT[4-j]=PENE;
  // ROTATE TO PUT PHOTON ALONG THIRD AXIS
  X1 = PHOT[1-j];
  X2 = PHOT[2-j]; 
  lortra(3,-FI3,PNEUTR,VEC,PAA,PP,PE);
  rotod3(-FI3,PHOT,PHOT);
  lortra(2,-TH3,PNEUTR,VEC,PAA,PP,PE);
  rotod2(-TH3,PHOT,PHOT);
  spaj(21,PNEUTR,PAA,PP,PE);                                             
  // ... now get the pair !
  double PAIRB=PENE/XMP+PPED/XMP;
  bostd3(PAIRB,PE,PE);  
  bostd3(PAIRB,PP,PP);   
  spaj(3,PNEUTR,PAA,PP,PE); 
  double GAMM=(PNEUTR[4-j]+PAA[4-j]+PP[4-j]+PE[4-j])/AMTO;
  BPAR=GAMM-sqrt(GAMM*GAMM-1.0);
  lortra(1, BPAR,PNEUTR,VEC,PAA,PP,PE);
  bostd3( BPAR,PHOT,PHOT);
  spaj(4,PNEUTR,PAA,PP,PE);                                             
  // BACK IN THE TAU REST FRAME BUT PNEUTR NOT YET ORIENTED.
  X1 = PNEUTR[1-j];
  X2 = PNEUTR[2-j];
  FI4  =angfi(X1,X2);
  X1 = PNEUTR[3-j];
  X2 = sqrt(PNEUTR[1-j]*PNEUTR[1-j]+PNEUTR[2-j]*PNEUTR[2-j]);
  TH4  =angxy(X1,X2);
  lortra(3, FI4,PNEUTR,VEC,PAA,PP,PE);
  rotod3( FI4,PHOT,PHOT);
  lortra(2,-TH4,PNEUTR,VEC,PAA,PP,PE);
  rotod2(-TH4,PHOT,PHOT);
  X1 = VEC[1-j];
  X2 = VEC[2-j];
  FI5=angfi(X1,X2);
  lortra(3,-FI5,PNEUTR,VEC,PAA,PP,PE);
  rotod3(-FI5,PHOT,PHOT);
  // PAA RESTORES ORIGINAL DIRECTION 
  lortra(3, FI2,PNEUTR,VEC,PAA,PP,PE);
  rotod3( FI2,PHOT,PHOT);
  lortra(2, TH1,PNEUTR,VEC,PAA,PP,PE);
  rotod2( TH1,PHOT,PHOT);
  lortra(3, FI1,PNEUTR,VEC,PAA,PP,PE);
  rotod3( FI1,PHOT,PHOT);
  spaj(10,PNEUTR,PAA,PP,PE);
  lortra(1,BSTB,PNEUTR,VEC,PAA,PP,PE);
  lortra(2,TH0,PNEUTR,VEC,PAA,PP,PE);
  lortra(3,FI0,PNEUTR,VEC,PAA,PP,PE);
  spaj(11,PNEUTR,PAA,PP,PE);
}  

