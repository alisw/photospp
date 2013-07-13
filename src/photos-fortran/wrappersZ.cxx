namespace Photospp
{

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

extern double AFBCALC(double SVAR,int IDEE,int IDFF);
extern "C" double afbcalc_(double *SVAR,int *IDEE,int *IDFF)
{
  return AFBCALC(*SVAR,*IDEE,*IDFF);
}

extern int GETIDEE(int IDE);
extern "C" int getidee_(int *IDE)
{
  return GETIDEE(*IDE);
}

extern double PHASYZ(double SVAR);
extern "C" double phasyz_(double *SVAR)
{
  return PHASYZ(*SVAR);
}

extern double Zphwtnlo(double svar,double xk, int IDHEP3, int IREP, double qp[4], double qm[4],double ph[4], double pp[4], double pm[4], double COSTHG, double BETA, double th1);

extern "C" double zphwtnlo_(double *svar,double *xk,int *IDHEP3,int *IREP,double qp[4],double qm[4],double ph[4],double pp[4], double pm[4],double *COSTHG, double *BETA, double *th1)
{
  return Zphwtnlo(*svar, *xk, *IDHEP3, *IREP,  qp,qm,ph,pp,pm, *COSTHG, *BETA, *th1);
}
extern double phwtnlo(double xdumm);
extern "C" double phwtnlo_(double *xdumm)
{
  return phwtnlo(*xdumm);
}

} // namespace Photospp

