#include "Photos.h"
using namespace Photospp;

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

extern double PHOBO3(double ANGLE,double PVEC[4]);

extern "C" double phobo3_(double *ANGLE,double PVEC[4])
{
  return PHOBO3(*ANGLE,PVEC);
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

extern double PHOB(int MODE,double PBOOS1[4],double VEC[4]);

extern "C" double phob_(int *MODE, double PBOOS1[4], double VEC[4])
{
  return PHOB(*MODE,PBOOS1,VEC);
}