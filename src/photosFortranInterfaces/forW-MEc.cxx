#include<complex>
#include<iostream>
using std::cout;
using std::endl;
using std::complex;

extern "C" struct{
  // COMMON /Kleiss_Stirling/spV,bet
  double spV,bet;
} kleiss_stirling_;

extern "C" struct{
  // COMMON /mc_parameters/pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN
  double pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN;
} mc_parameters_;





//////////////////////////////////////////////////////////////////
//         small s_{+,-}(p1,p2) for massless case:              //
//                 p1^2 = p2^2 = 0                              // 
//                                                              //
//     k0(0) = 1.d0                                             //
//     k0(1) = 1.d0                                             //
//     k0(2) = 0.d0  Kleisse_Stirling k0 points to X-axis       // 
//     k0(3) = 0.d0                                             //
//                                                              //
//////////////////////////////////////////////////////////////////
complex<double> InProd_zero(double p1[4],int l1,double p2[4],int l2){

  int i; 
  double  forSqrt1,forSqrt2,sqrt1,sqrt2;
  complex<double>    Dcmplx;
  static complex<double>    i_= complex<double>(0.0,1.0);
  bool           equal;



  equal = true;  
  for (int i = 0; i < 4; i++){
 
    if (p1[i]!=p2[i])  equal = equal && false ;
  }               


  if ( (l1==l2) || equal ) return complex<double>(0.0,0.0);

 
  else if ( (l1==+1) && (l2==-1) ){

    forSqrt1 = (p1[0]-p1[1])/(p2[0]-p2[1]);
    forSqrt2 = 1.0/forSqrt1;
    sqrt1    = sqrt(forSqrt2);
    sqrt2    = sqrt(forSqrt1);

    return (p1[2]+i_*p1[3])*sqrt1 -
	   (p2[2]+i_*p2[3])*sqrt2 ;
  }
  else if ( (l1==-1) && (l2==+1) ){

    forSqrt1 = (p1[0]-p1[1])/(p2[0]-p2[1]);
    forSqrt2 = 1.0/forSqrt1;
    sqrt1    = sqrt(forSqrt2);
    sqrt2    = sqrt(forSqrt1);

    return (p2[2]-i_*p2[3])*sqrt2 -
           (p1[2]-i_*p1[3])*sqrt1 ;
  }
  else{
                 

    cout << " "<<endl;             
    cout << " ERROR IN InProd_zero:"<<endl;
    cout << "   WRONG VALUES FOR l1,l2: l1,l2 = -1,+1 "<<endl;
    cout << " "  <<endl;           
    exit(0);
  }
}

double InSqrt(double p[4],double q[4]){
            
  return sqrt( (p[0]-p[1]) / (q[0]-q[1]) );
}
    


/////////////////////////////////////////////////////////////////////
//                                                                 //
//  this is small B_{s}(k,p) function when TrMartix is diaagonal!! //
//                                                                 //
/////////////////////////////////////////////////////////////////////
  complex<double>  BsFactor(int s,double k[4],double p[4],double m){
    double forSqrt1,sqrt1;
    complex<double>  inPr1;
  
  // COMMON /mc_paraneters/pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN 
  // temporary solution these will be global variables of the class   
  double pi     = mc_parameters_.pi;
  double alphaI = mc_parameters_.alphaI;
  double qb     = mc_parameters_.qb;
  double mf1    = mc_parameters_.mf1;
  double mf2    = mc_parameters_.mf2;
  double qf1    = mc_parameters_.qf1;
  double qf2    = mc_parameters_.qf2;
  // end of temporary solution

  if ( s==1 ){ 

    inPr1    = InProd_zero(k,+1,p,-1);
    forSqrt1 = (p[0]-p[1])/(k[0]-k[1]);
    sqrt1    = sqrt(2.0*forSqrt1);  
    //BsFactor = 
    return inPr1*sqrt1;
  }

  else if ( s==-1 ){

    inPr1    = InProd_zero(k,-1,p,+1);
    forSqrt1 = (p[0]-p[1])/(k[0]-k[1]);
    sqrt1    = sqrt(2.0*forSqrt1); 
    //BsFactor = 
    return inPr1*sqrt1;
  }
  else{

    cout << " "<<endl;             
    cout << " ERROR IN BsFactor: "<<endl;
    cout << "       WRONG VALUES FOR s : s = -1,+1"<<endl;
    cout << " "  <<endl;           
    exit(0);
  }
}




//====================================================================== 
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
complex<double> WDecayEikonalKS_1ph(double p3[4],double p1[4],double p2[4],double k[4],int s){

  double scalProd1,scalProd2,scalProd3;
  complex<double> wdecayeikonalks_1ph,BSoft1,BSoft2;  
  
  // COMMON /mc_paraneters/pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN 
  // temporary solution these will be global variables of the class   
  double pi     = mc_parameters_.pi;
  double alphaI = mc_parameters_.alphaI;
  double qb     = mc_parameters_.qb;
  double mf1    = mc_parameters_.mf1;
  double mf2    = mc_parameters_.mf2;
  double qf1    = mc_parameters_.qf1;
  double qf2    = mc_parameters_.qf2;
  // end of temporary solution

  scalProd1 = p1[0]*k[0]-p1[1]*k[1]-p1[2]*k[2]-p1[3]*k[3];
  scalProd2 = p2[0]*k[0]-p2[1]*k[1]-p2[2]*k[2]-p2[3]*k[3];
  scalProd3 = p3[0]*k[0]-p3[1]*k[1]-p3[2]*k[2]-p3[3]*k[3];


  BSoft1  = BsFactor(s,k,p1,mf1);
  BSoft2  = BsFactor(s,k,p2,mf2);
 
  //WDecayEikonalKS_1ph =   
   return sqrt(pi/alphaI)*(-(qf1/scalProd1+qb/scalProd3)*BSoft1   
                           +(qf2/scalProd2-qb/scalProd3)*BSoft2);

}
