extern double GIVIZO(int idferm,int ihelic,double *sizo3,double *charge,int *kolor);

extern "C" double givizo_(int *IDFERM,int *IHELIC,double *SIZO3,double *CHARGE,int *KOLOR)
{
  return GIVIZO(*IDFERM,*IHELIC,SIZO3,CHARGE,KOLOR);
}
extern double PHBORNM(double svar,double costhe,double T3e,double qe,double T3f,double qf,int Ncf);

extern "C" double phbornm_(double *svar,double *costhe,double *T3e,double *qe,double *T3f,double *qf,int *Ncf)
{
  return PHBORNM(*svar,*costhe,*T3e,*qe,*T3f,*qf,*Ncf);
}
