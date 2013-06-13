extern double GIVIZO(int idferm,int ihelic,double *sizo3,double *charge,int *kolor);

extern "C" double givizo_(int *IDFERM,int *IHELIC,double *SIZO3,double *CHARGE,int *KOLOR)
{
  return GIVIZO(*IDFERM,*IHELIC,SIZO3,CHARGE,KOLOR);
}
