#include "Photos.h"
#include <cmath>
#include <iostream>
using std::max;
using std::cout;
using std::endl;
using namespace Photospp;
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
