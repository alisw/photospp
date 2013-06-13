extern "C" struct{
  // COMMON /Kleiss_Stirling/spV,bet
  double spV,bet
} kleiss_stirling_;
extern "C" struct{
  // COMMON /mc_parameters/pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN
  double pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN
} mc_parameters_;


extern double complex WDecayEikonalKS_1ph(double p3[4],double p1[4],double p2[4],int s);

extern "C" double complex WDecayEikonalKS_1ph_(double p3[4],double p1[4],double p2[4],int s)
{
  return WDecayEikonalKS_1ph(p3,p1,p2,s);
}
