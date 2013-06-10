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
