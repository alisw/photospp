#include <complex>
using std::complex;

extern complex<double> WDecayEikonalKS_1ph(double p3[4],double p1[4],double p2[4],double k[4],int s);

extern "C" complex<double> wdecayeikonalks_1ph_(double p3[4],double p1[4],double p2[4],double k[4],int *s)
{
  return WDecayEikonalKS_1ph(p3,p1,p2,k,*s);
}
