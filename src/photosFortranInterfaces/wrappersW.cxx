#include <complex>
using std::complex;


extern double WDecayEikonalSqrKS_1ph(double p3[4],double p1[4],double p2[4],double k[4]);

extern "C" double wdecayeikonalsqrks_1ph_(double p3[4],double p1[4],double p2[4],double k[4])
{
  return WDecayEikonalSqrKS_1ph(p3,p1,p2,k);
}


extern double WDecayBornAmpSqrKS_1ph(double p3[4],double p1[4],double p2[4]);

extern "C" double wdecaybornampsqrks_1ph_(double p3[4],double p1[4],double p2[4])
{
  return WDecayBornAmpSqrKS_1ph(p3,p1,p2);
}


extern double WDecayAmplitudeSqrKS_1ph(double p3[4],double p1[4],double p2[4],double k[4]);

extern "C" double wdecayamplitudesqrks_1ph_(double p3[4],double p1[4],double p2[4],double k[4])
{
  return WDecayAmplitudeSqrKS_1ph(p3,p1,p2,k);
}



extern double SANC_WT(double PW[4],double PNE[4],double PMU[4],double PPHOT[4],double B_PW[4],double B_PNE[4],double B_PMU[4]);

extern "C" double sanc_wt_(double PW[4],double PNE[4],double PMU[4],double PPHOT[4],double B_PW[4],double B_PNE[4],double B_PMU[4])
{
  return SANC_WT(PW,PNE,PMU,PPHOT,B_PW,B_PNE,B_PMU);
}


void SANC_INIT1(double QB0,double QF20,double MF10,double MF20,double MB0);

extern "C" void sanc_init1_(double *QB0,double *QF20,double *MF10,double *MF20,double *MB0)
{
  return SANC_INIT1(*QB0,*QF20,*MF10,*MF20,*MB0);
}


void SANC_INIT(double ALPHA,int PHLUN);

extern "C" void sanc_init_(double *ALPHA,int *PHLUN)
{
  return SANC_INIT(*ALPHA,*PHLUN);
}


void PHOBWnlo(double WT);

extern "C" void phobwnlo_(double *WT)
{
  return PHOBWnlo(*WT);
}

