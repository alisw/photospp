#include<complex>
using std::complex;

extern "C" struct{
  // COMMON /Kleiss_Stirling/spV,bet
  double spV,bet;
} kleiss_stirling_;

extern "C" struct{
  // COMMON /mc_parameters/pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN
  double pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN;
} mc_parameters_;

//=====================================================================                 
//     
// Eikonal factor of decay W->l_1+l_2+\gamma in terms of K&S objects !
// 
//   EikFactor = q1*eps.p1/k.p1 + q2*eps.p2/k.p2 - q3*eps.p3/k.p3
//
//   indices 1,2 are for charged decay products
//   index 3 is for W
//   
//   q - charge
//    
//======================================================================

extern "C" complex<double> bsfactor_(int *s, double *k, double *p, double *m);

complex<double> WDecayEikonalKS_1ph(double p3[4],double p1[4],double p2[4],double k[4],int s){

  double scalProd1,scalProd2,scalProd3;
  complex<double> wdecayeikonalks_1ph,BSoft1,BSoft2;  
  
  // COMMON /mc_paraneters/pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN    
  double &pi     = mc_parameters_.pi;
  double &alphaI = mc_parameters_.alphaI;
  double &qb     = mc_parameters_.qb;
  double &mf1    = mc_parameters_.mf1;
  double &mf2    = mc_parameters_.mf2;
  double &qf1    = mc_parameters_.qf1;
  double &qf2    = mc_parameters_.qf2;

  scalProd1 = p1[0]*k[0]-p1[1]*k[1]-p1[2]*k[2]-p1[3]*k[3];
  scalProd2 = p2[0]*k[0]-p2[1]*k[1]-p2[2]*k[2]-p2[3]*k[3];
  scalProd3 = p3[0]*k[0]-p3[1]*k[1]-p3[2]*k[2]-p3[3]*k[3];

  BSoft1  = bsfactor_(&s,k,p1,&mf1);
  BSoft2  = bsfactor_(&s,k,p2,&mf2);
 
  wdecayeikonalks_1ph = 
        sqrt(pi/alphaI)*(-(qf1/scalProd1+qb/scalProd3)*BSoft1   
                         +(qf2/scalProd2-qb/scalProd3)*BSoft2);

  return wdecayeikonalks_1ph;
}
