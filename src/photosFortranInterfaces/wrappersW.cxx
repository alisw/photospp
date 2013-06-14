#include <complex>
using std::complex;

extern complex<double> InProd_zero(double p1[4],int l1,double p2[4],int l2);

extern "C" complex<double> inprod_zero_(double p1[4],int *l1,double p2[4],int *l2)
{
  return InProd_zero(p1,*l1,p2,*l2);
}

double InSqrt(double p[4],double q[4]);

extern "C" double insqrt_(double p[4],double q[4])
{
  return InSqrt(p,q);
}


extern double InProd_mass(double p1[4],double m1,int l1,double p2[4],double m2,int l2);

extern "C" double insqrtinprod_mass_(double p1[4],double *m1,int *l1,double p2[4],double *m2,int *l2)
{
  return InProd_mass(p1,*m1,*l1,p2,*m2,*l2);
}

extern complex<double> BsFactor(int s,double k[4],double p[4],double m);

extern "C" complex<double> bsfactor_(int *s,double k[4],double p[4],double *m)
{
  return BsFactor(*s,k,p,*m);
}

extern complex<double> WDecayEikonalKS_1ph(double p3[4],double p1[4],double p2[4],double k[4],int s);

extern "C" complex<double> wdecayeikonalks_1ph_(double p3[4],double p1[4],double p2[4],double k[4],int *s)
{
  return WDecayEikonalKS_1ph(p3,p1,p2,k,*s);
}

extern complex<double> SoftFactor(int s,double k[4],double p1[4],double m1,double p2[4],double m2,double Gmass2);

extern "C" complex<double> softfactor_(int *s,double k[4],double p1[4],double *m1,double p2[4],double *m2,double *Gmass2)
{
  return SoftFactor(*s,k,p1,*m1,p2,*m2,*Gmass2);
}

extern complex<double> TrMatrix_zero(double p1[4],double m1,int l1,double k[4],int s,double p2[4],double m2,int l2);

extern "C" complex<double> trmatrix_zero_(double p1[4],double *m1,int *l1,double k[4],int *s,double p2[4],double *m2,int *l2)
{
  return TrMatrix_zero(p1,*m1,*l1,k,*s,p2,*m2,*l2);
}

extern complex<double> TrMatrix_mass(double p1[4],double m1,int l1,double k[4],double m,int s,double p2[4],double m2,int l2);

extern "C" complex<double> trmatrix_mass_(double p1[4],double *m1,int *l1,double k[4],double *m,int *s,double p2[4],double *m2,int *l2)
{
  return TrMatrix_mass(p1,*m1,*l1,k,*m,*s,p2,*m2,*l2);
}

extern complex<double> WDecayBornAmpKS_1ph(double p3[4],int l3,double p1[4],int l1,double p2[4],int l2);

extern "C" complex<double> wdecaybornampks_1ph_(double p3[4],int *l3,double p1[4],int *l1,double p2[4],int *l2)
{
  return WDecayBornAmpKS_1ph(p3,*l3,p1,*l1,p2,*l2);
}

extern complex<double> WDecayAmplitudeKS_1ph(double p3[4],int l3,double p1[4],int l1,double p2[4],int l2,double k[4],int s);

extern "C" complex<double> wdecayamplitudeks_1ph_(double p3[4],int *l3,double p1[4],int *l1,double p2[4],int *l2,double k[4],int *s)
{
  return WDecayAmplitudeKS_1ph(p3,*l3,p1,*l1,p2,*l2,k,*s);
}


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
