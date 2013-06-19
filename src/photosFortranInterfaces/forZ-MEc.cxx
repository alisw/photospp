#include "Photos.h"
#include <cmath>
#include <iostream>
using std::max;
using std::cout;
using std::endl;
using namespace Photospp;

// from photosC.cxx
extern void GETIDEIDF(int *IDE, int *IDF);

// ----------------------------------------------------------------------
// PROVIDES ELECTRIC CHARGE AND WEAK IZOSPIN OF A FAMILY FERMION
// IDFERM=1,2,3,4 DENOTES NEUTRINO, LEPTON, UP AND DOWN QUARK
// NEGATIVE IDFERM=-1,-2,-3,-4, DENOTES ANTIPARTICLE
// IHELIC=+1,-1 DENOTES RIGHT AND LEFT HANDEDNES ( CHIRALITY)
// SIZO3 IS THIRD PROJECTION OF WEAK IZOSPIN (PLUS MINUS HALF)
// AND CHARGE IS ELECTRIC CHARGE IN UNITS OF ELECTRON CHARGE
// KOLOR IS A QCD COLOUR, 1 FOR LEPTON, 3 FOR QUARKS
//
//     called by : EVENTE, EVENTM, FUNTIH, .....
// ----------------------------------------------------------------------

void GIVIZO(int IDFERM,int IHELIC,double *SIZO3,double *CHARGE,int *KOLOR) {
  //
  int IH, IDTYPE, IC, LEPQUA, IUPDOW; 
  if (IDFERM==0 || abs(IDFERM)>4 || abs(IHELIC)!=1){
    cout << "STOP IN GIVIZO: WRONG PARAMS" << endl;
    exit(0);
   }

  IH  =IHELIC;
  IDTYPE =abs(IDFERM);
  IC  =IDFERM/IDTYPE;
  LEPQUA=(int)(IDTYPE*0.4999999);
  IUPDOW=IDTYPE-2*LEPQUA-1;
  *CHARGE  =(-IUPDOW+2.0/3.0*LEPQUA)*IC;
  *SIZO3   =0.25*(IC-IH)*(1-2*IUPDOW);
  *KOLOR=1+2*LEPQUA;
  //** NOTE THAT CONVENTIONALY Z0 COUPLING IS
  //** XOUPZ=(SIZO3-CHARGE*SWSQ)/SQRT(SWSQ*(1-SWSQ))
  return;
}


////////////////////////////////////////////////////////////////////////////
///                                                                       //
/// This routine provides unsophisticated Born differential cross section //
/// at the crude x-section level, with Z and gamma s-chanel exchange.     //
///////////////////////////////////////////////////////////////////////////
double PHBORNM(double svar,double costhe,double T3e,double qe,double T3f,double qf,int NCf){

  double   s,t,Sw2,MZ,MZ2,GammZ,MW,MW2,AlfInv,GFermi;
  double   sum,deno,Ve,Ae,thresh;
  double   xe,yf,xf,ye,ff0,ff1,amx2,amfin,Vf,Af;
  double   ReChiZ,SqChiZ,RaZ,RaW,ReChiW,SqChiW;
  double   Born, BornS;
  int  KeyZet,HadMin,KFbeam;
  int  i,ke,KFfin,kf,IsGenerated,iKF;
  int  KeyWidFix;
  // we may want to use phpico_.pi phpico_.twopi defined in Photos::initialize()
  static double PI=3.14159265358979324, TWOPI=6.28318530717958648;
 
  AlfInv= 137.0359895;
  GFermi=1.16639e-5;

  //--------------------------------------------------------------------
  s = svar;
  //------------------------------
  //     EW paratemetrs taken from BornV
  MZ=91.187;
  GammZ=2.50072032;
  Sw2=.22276773;
  //------------------------------
  // Z and gamma couplings to beams (electrons)
  // Z and gamma couplings to final fermions
  // Loop over all flavours defined in m_xpar(400+i)


  //------ incoming fermion
  Ve=  2*T3e -4*qe*Sw2;
  Ae=  2*T3e;
  //------ final fermion couplings
  amfin = 0.000511; //  m_xpar(kf+6)
  Vf =  2*T3f -4*qf*Sw2;
  Af =  2*T3f;
  if(abs(costhe) > 1.0){
    cout << "+++++STOP in PHBORN: costhe>0 =" << costhe << endl;
    exit(0);
  }
  MZ2  = MZ*MZ;
  RaZ  = (GFermi *MZ2 *AlfInv  )/( sqrt(2.0) *8.0 *PI); //
  RaZ  = 1/(16.0*Sw2*(1.0-Sw2));
  KeyWidFix = 1;       // fixed width
  KeyWidFix = 0;       // variable width
  if( KeyWidFix == 0 ){
    ReChiZ=(s-MZ2)*s/((s-MZ2)*(s-MZ2)+(GammZ*s/MZ)*(GammZ*s/MZ)) *RaZ;     // variable width
    SqChiZ=      s*s/((s-MZ2)*(s-MZ2)+(GammZ*s/MZ)*(GammZ*s/MZ)) *RaZ*RaZ; // variable width
  }
  else{
      ReChiZ=(s-MZ2)*s/((s-MZ2)*(s-MZ2)+(GammZ*MZ)*(GammZ*MZ)) *RaZ;     // fixed width
      SqChiZ=      s*s/((s-MZ2)*(s-MZ2)+(GammZ*MZ)*(GammZ*MZ)) *RaZ*RaZ; // fixed width
  }
  xe= Ve*Ve +Ae*Ae;
  xf= Vf*Vf +Af*Af;
  ye= 2*Ve*Ae;
  yf= 2*Vf*Af;
  ff0= qe*qe*qf*qf +2*ReChiZ*qe*qf*Ve*Vf +SqChiZ*xe*xf;
  ff1=             +2*ReChiZ*qe*qf*Ae*Af +SqChiZ*ye*yf;
  Born    = (1.0+ costhe*costhe)*ff0 +2.0*costhe*ff1;
  // Colour factor
  Born = NCf*Born;
  // Crude method of correcting threshold, cos(theta) depencence incorrect!!!
  if(    svar <  4.0*amfin*amfin){ 
    thresh=0.0;
  }
  else if(svar < 16.0*amfin*amfin){
    amx2=4.0*amfin*amfin/svar;
    thresh=sqrt(1.0-amx2)*(1.0+amx2/2.0);
  }
  else{
    thresh=1.0;
  }

  Born= Born*thresh;
  return Born;
}


// ----------------------------------------------------------------------
// THIS ROUTINE CALCULATES  BORN ASYMMETRY.
// IT EXPLOITS THE FACT THAT BORN X. SECTION = A + B*C + D*C**2
//
//     called by : EVENTM
// ----------------------------------------------------------------------
//
double AFBCALC(double SVAR,int IDEE,int IDFF){
  int KOLOR,KOLOR1;
  double T3e,qe,T3f,qf,A,B;
  GIVIZO(IDEE,-1,&T3e,&qe,&KOLOR);
  GIVIZO(IDFF,-1,&T3f,&qf,&KOLOR1);

  A=PHBORNM(SVAR,0.5,T3e,qe,T3f,qf,KOLOR*KOLOR1);
  B=PHBORNM(SVAR,-0.5,T3e,qe,T3f,qf,KOLOR*KOLOR1);
  return (A-B)/(A+B)*5.0/2.0 *3.0/8.0;
}


int GETIDEE(int IDE){

  int IDEE;
  if((IDE==11)       || (IDE== 13) || (IDE== 15)){
    IDEE=2;
  }
  else if((IDE==-11) || (IDE==-13) || (IDE==-15)){
    IDEE=-2;
  }
  else if((IDE== 12) || (IDE== 14) || (IDE== 16)){
    IDEE=1;
  }
  else if((IDE==-12) || (IDE==-14) || (IDE==-16)){
    IDEE=-1;
  }
  else if((IDE==  1) || (IDE==  3) || (IDE==  5)){
    IDEE=4;
  }
  else if((IDE== -1) || (IDE== -3) || (IDE== -5)){
    IDEE=-4;
  }
  else if((IDE==  2) || (IDE==  4) || (IDE==  6)){
    IDEE=3;
  }
  else if((IDE==- 2) || (IDE== -4) || (IDE== -6)){
    IDEE=-3;
  }

  return IDEE;
}


//----------------------------------------------------------------------
//
//    PHASYZ:   PHotosASYmmetry of Z
//
//    Purpose:  Calculates born level asymmetry for Z
//              between distributions (1-C)**2 and (1+C)**2
//              At present dummy, requrires effective Z and gamma 
//              Couplings and also spin polarization states
//              For initial and final states.
//              To be correct this function need to be tuned
//              to host generator. Axes orientation polarisation
//              conventions etc etc. 
//              Modularity of PHOTOS would break. 
//
//    Input Parameters:   SVAR
//
//    Output Parameters:  Function value
//
//    Author(s):  Z. Was                          Created at:  10/12/05
//                                                Last Update: 19/06/13
//
//----------------------------------------------------------------------
double PHASYZ(double SVAR){


  double AFB;
  int IDE,IDF,IDEE,IDFF;
  GETIDEIDF(&IDE,&IDF);
  IDEE=abs(GETIDEE(IDE));
  IDFF=abs(GETIDEE(IDF));
  AFB= -AFBCALC(SVAR,IDEE,IDFF);
  //      AFB=0
  return 4.0/3.0*AFB;
  //      write(*,*) 'IDE=',IDE,'  IDF=',IDF,'  SVAR=',SVAR,'AFB=',AFB
}

