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
