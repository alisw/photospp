#include "f_Init.h"
#include "PH_HEPEVT_Interface.h"
#include "PhotosUtilities.h"
using namespace Photospp;
using PhotosUtilities::PHOCHA;

// not in any include file (yet)
extern "C" void partra_(int *pIBRAN,double PHOT[4]);
extern "C" void trypar_(bool *pJESLI,double *pSTRENG,double PA[4],double PB[4],double PE[4],double PP[4]);

/*----------------------------------------------------------------------

      PHOTOS:   Photon radiation in decays

      Purpose:  e+e- pairs  are  generated  in
                the decay of the IPPAR-th particle in the HEP-like
                common /PHOEVT/.  Radiation takes place from one
                of the charged daughters of the decaying particle IPPAR



      Input Parameter:    IPPAR:  Pointer   to   decaying  particle  in
                                  /PHOEVT/ and the common itself,
                                  NHEP0 length of the /HEPEVT/ entry
                                  before starting any activity on this
                                  IPPAR decay.
      Output Parameters:  Common  /HEPEVT/, either  with  or  without a
                                  e+e-(s) added.


      Author(s):  Z. Was,                         Created at:  01/06/93
                                                  Last Update:

  ----------------------------------------------------------------------*/
void PHOPAR(int IPARR,int NHEP0) {
  double MINMAS,STRENG,PCHAR[4],PNEU[4],PELE[4],PPOZ[4],BUF[4];
  float  MASSUM;
  int    IP,IPPAR,NLAST;
  bool   BOOST,JESLI;

  IPPAR = IPARR;
  // Store pointers for cascade treatment...
  IP    = IPPAR - 1;
  NLAST = pho.nhep;

  // Check decay multiplicity..
  PHOIN(IPPAR,&BOOST,&NHEP0);
  PHOCHK(hep.jdahep[IP][0]);
  PHLUPA(100);
  if(pho.jdahep[IP][0] == 0) return;
  if(pho.jdahep[IP][0] == pho.jdahep[IP][1]) return;

  // Loop over charged daughters
  for(int I=pho.jdahep[IP][0]-1; I <= pho.jdahep[IP][1]-1; ++I) {

    // Skip this particle if it has no charge
    if( PHOCHA(pho.idhep[I]) == 0 ) continue;

    // Set  3-vectors
    for(int J = 0; J < 3; ++J) {
      PCHAR[J] = pho.phep[I][J];
      PNEU [J] =-pho.phep[I][J];
    }

    // Set energy
    PNEU[3]  = pho.phep[IP][4] - pho.phep[I][3];
    PCHAR[3] = pho.phep[I][3];

    STRENG   = 1.0;

    //here we attempt generating pair from PCHAR. One of the charged
    //decay products; that is why algorithm works in a loop.
    //PNEU is four vector of all decay products except PCHAR
    //we do not care on rare cases when two pairs could be generated
    //we assume it is negligibly rare and fourth order in alpha anyway
    //TRYPAR should take as an input electron mass.
    //then it can be used for muons.

    trypar_(&JESLI,&STRENG,PCHAR,PNEU,PELE,PPOZ);

    //emitted pair four momenta are stored in PELE PPOZ
    //then JESLI=.true.

    // If JESLI = true, we modify old particles of the vertex
    if (JESLI) {

      // we have to correct 4-momenta
      // of all decay products
      // we use PARTRA for that
      // PELE PPOZ are in right frame
      for(int J = pho.jdahep[IP][0]-1; J<pho.jdahep[IP][1]-1; ++J) {
        for(int K = 0; K<4; ++K) {
          BUF[K] = pho.phep[I][K];
        }
        if (J == I) {
          int _TEMP = 1;
          partra_(&_TEMP,BUF);
        } else {
          int _TEMP = -1;
          partra_(&_TEMP,BUF);
        }
        for(int K = 0; K<4; ++K) {
          pho.phep[I][K] = BUF[K];
        }
      }

      PHLUPA(1011);

      // electron: adding to vertex
      pho.nhep = pho.nhep+1;
      pho.isthep[pho.nhep] = 1;
      pho.idhep [pho.nhep] = 11;
      pho.jmohep[pho.nhep][0] = IP;
      pho.jmohep[pho.nhep][1] = 0;
      pho.jdahep[pho.nhep][0] = 0;
      pho.jdahep[pho.nhep][1] = 0;

      for(int K = 1; K<4; ++K) {
        pho.phep[pho.nhep][K] = PELE[K];
      }

      pho.phep[pho.nhep][4] = 0.000511;

      // positron: adding
      pho.nhep = pho.nhep+1;
      pho.isthep[pho.nhep] = 1;
      pho.idhep [pho.nhep] =-11;
      pho.jmohep[pho.nhep][0] = IP;
      pho.jmohep[pho.nhep][1] = 0;
      pho.jdahep[pho.nhep][0] = 0;
      pho.jdahep[pho.nhep][1] = 0;

      for(int K = 1; K<4; ++K) {
        pho.phep[pho.nhep][K] = PPOZ[K];
      }

      pho.phep[pho.nhep][4] = 0.000511;

      // write in
      PHLUPA(1012);
      PHOOUT(IPPAR, BOOST, NHEP0);
      PHOIN (IPPAR,&BOOST,&NHEP0);
      PHLUPA(1013);
    } // end of if (JESLI)
  } // end of loop over charged particles
}
