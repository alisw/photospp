#include "Photos.h"
#include "PhotosUtilities.h"
using namespace Photospp;
using namespace PhotosUtilities;

extern bool F(int m, int i);

extern "C" bool f_(int *IDABS)
{
  return F(0,*IDABS);
}

extern void PHOEPS(double vec1[4], double vec2[4], double eps[4]);

extern "C" void phoeps_(double VEC1[4], double VEC2[4], double EPS[4])
{
  PHOEPS(VEC1,VEC2,EPS);
}

extern double PHOSPI(int idabs);

extern "C" double phospi_(int *IDABS)
{
  return PHOSPI(*IDABS);
}

extern double PHOCHA(int idhep);

extern "C" double phocha_(int *IDHEP)
{
  return PHOCHA(*IDHEP);
}

extern double PHINT(int idhep);

extern "C" double phint_(int *IDHEP)
{
  return PHINT(*IDHEP);
}
extern double PHINT1(int idhep);

extern "C" double phint1_(int *IDHEP)
{
  return PHINT1(*IDHEP);
}
extern double PHINT2(int idhep);

extern "C" double phint2_(int *IDHEP)
{
  return PHINT2(*IDHEP);
}

extern double PHOTRI(double A,double B,double C);

extern "C" double photri_(double *A,double *B,double *C)
{
  return PHOTRI(*A,*B,*C);
}

extern double PHOAN1(double X,double Y);

extern "C" double phoan1_(double *X,double *Y)
{
  return PHOAN1(*X,*Y);
}

extern double PHOAN2(double X,double Y);

extern "C" double phoan2_(double *X,double *Y)
{
  return PHOAN2(*X,*Y);
}

extern "C" void phobo3_(double *ANGLE,double PVEC[4])
{
  PHOBO3(*ANGLE,PVEC);
}

extern double PHORO3(double ANGLE,double PVEC[4]);

extern "C" double phoro3_(double *ANGLE,double PVEC[4])
{
  return PHORO3(*ANGLE,PVEC);
}

extern double PHORO2(double ANGLE,double PVEC[4]);

extern "C" double phoro2_(double *ANGLE,double PVEC[4])
{
  return PHORO2(*ANGLE,PVEC);
}

extern void PHOB(int MODE,double PBOOS1[4],double VEC[4]);

extern "C" void phob_(int *MODE, double PBOOS1[4], double VEC[4])
{
  return PHOB(*MODE,PBOOS1,VEC);
}

extern void bostdq(int mode,double qq[4],double pp[4],double r[4]);

extern "C" void bostdq_(int *mode,double qq[4],double pp[4],double r[4])
{
  return bostdq(*mode,qq,pp,r);
}


extern void PHOERR(int IMES,char *TEXT,double DATA);

extern "C" void phoerr_(int *IMES,char *TEXT,double *DATA)
{
  return PHOERR(*IMES,TEXT,*DATA);
}


extern void PHOREP();

extern "C" void phorep_()
{
  return PHOREP();
}



extern void GETIDEIDF(int *IDE, int *IDF);

extern "C" void getideidf_(int *IDE, int *IDF)
{
  return GETIDEIDF(IDE,IDF);
}


extern void PHODMP();

extern "C" void phodmp_()
{
  return PHODMP();
}


extern void PHLUPAB(int IPOINT);

extern "C" void phlupab_(int *IPOINT)
{
  return PHLUPAB(*IPOINT);
}

extern void PHLUPA(int IPOINT);

extern "C" void phlupa_(int *IPOINT)
{
  return PHLUPA(*IPOINT);
}


extern void PHOtoRF();

extern "C" void photorf_()
{
  return PHOtoRF();
}

extern void PHOtoLAB();

extern "C" void photolab_()
{
  return PHOtoLAB();
}
extern void PHOTOS_MAKE_C(int IPPAR);

extern "C" void photos_make_c_(int *IPPAR)
{
  return PHOTOS_MAKE_C(*IPPAR);
}
extern void PHCORK(int MODCOR);

extern "C" void phcork_(int *MODCOR)
{
  return PHCORK(*MODCOR);
}

extern void PHODO(int IP,int NCHARB,int NEUDAU);

extern "C" void phodo_(int *IP,int *NCHARB,int *NEUDAU)
{
  PHODO(*IP,*NCHARB,*NEUDAU);
}

extern void PHOBW(double *WT);

extern "C" void phobw_(double *WT)
{
  return PHOBW(WT);
}


extern double PHOFAC(int MODE);

extern "C" double phofac_(int *MODE)
{
  return PHOFAC(*MODE);
}


extern double PHOCORN(double MPASQR,double MCHREN,int ME);

extern "C" double phocorn_(double *MPASQR,double *MCHREN, int *ME)
{
  return PHOCORN(*MPASQR,*MCHREN,*ME);
}

extern double PHOCOR(double MPASQR,double MCHREN,int ME);

extern "C" double phocor_(double *MPASQR,double *MCHREN, int *ME)
{
  return PHOCOR(*MPASQR,*MCHREN,*ME);
}


extern void PHOTWO(int MODE);

extern "C" void photwo_(int *MODE)
{
  return PHOTWO(*MODE);
}

extern void PHOIN(int IP, bool *BOOST, int nhep0);

extern "C" void phoin_(int *IP, bool *BOOST, int *nhep0)
{
  return PHOIN(*IP,BOOST,*nhep0);
}


extern void PHOOUT(int IP, bool BOOST, int nhep0);

extern "C" void phoout_(int *IP, bool *BOOST, int *nhep0)
{
  return PHOOUT(*IP,*BOOST,*nhep0);
}
