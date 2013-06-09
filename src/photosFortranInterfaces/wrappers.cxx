#include "Photos.h"
using namespace Photospp;

extern bool F(int m, int i);

extern "C" bool f_(int *IDABS)
{
  return F(0,*IDABS);
}
