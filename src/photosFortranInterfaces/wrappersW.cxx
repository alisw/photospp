#include <complex>
using std::complex;

extern complex<double> InProd_zero(double p1[4],int l1,double p2[4],int l2);

extern "C" complex<double> inprod_zero_(double p1[4],int l1,double p2[4],int l2)
{
  return InProd_zero(p1,l1,p2,l2);
}

double InSqrt(double p[4],double q[4]);

extern "C" double insqrt_(double p[4],double q[4])
{
  return InSqrt(p,q);
}


extern complex<double> BsFactor(int s,double k[4],double p[4],double m);

extern "C" complex<double> bsfactor_(int s,double k[4],double p[4],double m)
{
  return BsFactor(s,k,p,m);
}

extern complex<double> WDecayEikonalKS_1ph(double p3[4],double p1[4],double p2[4],double k[4],int s);

extern "C" complex<double> wdecayeikonalks_1ph_(double p3[4],double p1[4],double p2[4],double k[4],int *s)
{
  return WDecayEikonalKS_1ph(p3,p1,p2,k,*s);
}
