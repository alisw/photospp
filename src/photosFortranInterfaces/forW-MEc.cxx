#include<complex>
#include<iostream>
using std::cout;
using std::endl;
using std::complex;

extern "C" struct{
  // COMMON /Kleiss_Stirling/spV,bet
  double spV[4],bet[4];
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
    
//////////////////////////////////////////////////////////////////
//                                                              //
//  Inner product for massive spinors: Ub(p1,m1,l1)*U(p2,m2,l2) //
//                                                              //
//////////////////////////////////////////////////////////////////

complex<double> InProd_mass(double p1[4],double m1,int l1,double p2[4],double m2,int l2){
  double sqrt1,sqrt2,forSqrt1;


  if ((l1==+1)&&(l2==+1)) {               
    forSqrt1    = (p1[0]-p1[1])/(p2[0]-p2[1]);
    sqrt1       = sqrt(forSqrt1);
    sqrt2       = 1.0/sqrt1;
    return complex<double>(m1*sqrt2+m2*sqrt1,0.0);
  }
  else if  ((l1==+1)&&(l2==-1))                              
    return InProd_zero(p1,+1,p2,-1);

  else if ((l1==-1)&&(l2==+1))                         
    return  InProd_zero(p1,-1,p2,+1);               

  else if ((l1==-1)&&(l2==-1)){                             
    forSqrt1    = (p1[0]-p1[1])/(p2[0]-p2[1]);
    sqrt1       = sqrt(forSqrt1);
    sqrt2       = 1.0/sqrt1;
    return complex<double>(m1*sqrt2+m2*sqrt1,0.0);
  }
  else {        
    cout <<" " <<endl;            
    cout <<" ERROR IN InProd_mass.."<<endl;
    cout <<"       WRONG VALUES FOR l1,l2"<<endl;
    cout <<" " <<endl;            
    exit(0);
  }
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

//======================================================================
//
//       Gauge invariant soft factor for decay!!
//       Gmass2 -- photon mass square       
// 
//======================================================================
complex<double>  SoftFactor(int s,double k[4],double p1[4],double m1,double p2[4],double m2,double Gmass2){

  double ScalProd1,ScalProd2;
  complex<double>  BsFactor2,BsFactor1;
           

  ScalProd1 = k[0]*p1[0]-k[1]*p1[1]-k[2]*p1[2]-k[3]*p1[3];
  ScalProd2 = k[0]*p2[0]-k[1]*p2[1]-k[2]*p2[2]-k[3]*p2[3];
          
  BsFactor1 = BsFactor(s,k,p1,m1);
  BsFactor2 = BsFactor(s,k,p2,m2);

  return + BsFactor2/2.0/(ScalProd2-Gmass2)
	 - BsFactor1/2.0/(ScalProd1-Gmass2);
}

//############################################################################# 
//                                                                            #
//                         \ eps(k,0,s)                                       # 
//                         /                                                  #   
//                        _\                                                  # 
//                         /\                                                 #
//                         \                                                  #
//                         /                                                  #
//           ---<----------\-------------<---                                 #
//       Ub(p1,m1,l1)                  U(p2,m2,l2)                            #
//                                                                            #
//                                                                            #
//             definition of arbitrary light-like vector beta!!               #
//                                                                            #
//              bet[0] = 1.d0                                                 #
//              bet[1] = 1.d0                                                 #
//              bet[2] = 0.d0      <==> bet == k0  expression becomes easy!!  #
//              bet[3] = 0.d0                                                 #
//#############################################################################

complex<double> TrMatrix_zero(double p1[4],double m1,int l1,double k[4],int s,double p2[4],double m2,int l2){

  double forSqrt1,forSqrt2;
  double p1_1[4],p2_1[4];
  double sqrt1,sqrt2,scalProd1,scalProd2;
  complex<double>   inPr1,inPr2,inPr3;
  bool          equal;

  equal = true;    
  for (int i = 0; i < 4; i++) 
    if (p1[i] != p2[i])  equal = equal&&false;

                    

  if ( (m1==m2)&&(equal) ){
    //..          
    //..             when:  p1=p2=p <=> m1=m2 TrMatrix_zero is diagonal
    //..               
    if ( (l1==+1)&&(l2==+1) ){ 

      inPr1    = InProd_zero(k,+s,p1,-s);
      forSqrt1 = (p1[0]-p1[1])/(k[0]-k[1]); 
      sqrt1    = sqrt(2.0*forSqrt1);

      return sqrt1*inPr1;
    }  
 
    else if ( (l1==+1)&&(l2==-1) ){                

      return complex<double>(0.0,0.0);}
                     

    else if ( (l1==-1)&&(l2==+1) ){               

      return complex<double>(0.0,0.0);
    } 

    else if ( (l1==-1)&&(l2==-1) ){                

      inPr1    = InProd_zero(k,+s,p1,-s);
      forSqrt1 = (p1[0]-p1[1])/(k[0]-k[1]); 
      sqrt1    = sqrt(2.0*forSqrt1);

      return sqrt1*inPr1;
    }  
          
    else{ 
        
      cout << ""  <<endl;           
      cout << " ERROR IN  TrMatrix_zero: " <<endl;
      cout << "       WRONG VALUES FOR l1,l2,s" <<endl; 
      cout <<  "" <<endl;             
      exit(0);

    }       

  }

  if ( (l1==+1)&&(l2==+1)&&(s==+1) ){

    inPr1    = InProd_zero(k,+1,p1,-1);
    forSqrt1 = (p2[0]-p2[1])/(k[0]-k[1]);
    sqrt1    = sqrt(2.0*forSqrt1);                   
 
    return sqrt1*inPr1;
  }
  else if ( (l1==+1)&&(l2==-1)&&(s==+1) ) {
 
    return complex<double>(0.0,0.0);
  }

  else if( (l1==-1)&&(l2==+1)&&(s==+1) ){
  
    forSqrt1 = (p1[0]-p1[1])/(p2[0]-p2[1]);             
    forSqrt2 = 1.0/forSqrt1;
    sqrt1    = sqrt(2.0*forSqrt1);                   
    sqrt2    = sqrt(2.0*forSqrt2);                   
                     
    return complex<double>(m2*sqrt1-m1*sqrt2,0.0);
  }
  else if ( (l1==-1)&&(l2==-1)&&(s==+1) ){ 

    inPr1    = InProd_zero(k,+1,p2,-1);
    forSqrt1 = (p1[0]-p1[1])/(k[0]-k[1]);
    sqrt1    = sqrt(2.0*forSqrt1);                   
  
    return inPr1*sqrt1;
  }

  else if ( (l1==+1)&&(l2==+1)&&(s==-1) ){
 
    inPr1    = -InProd_zero(k,-1,p2,+1);
    forSqrt1 = (p1[0]-p1[1])/(k[0]-k[1]);
    sqrt1    = sqrt(2.0*forSqrt1);                   
 
    return   -sqrt1*inPr1;
  }

  else if ( (l1==+1)&&(l2==-1)&&(s==-1) ){ 
           
    forSqrt1 = (p1[0]-p1[1])/(p2[0]-p2[1]);     
    forSqrt2 = 1.0/forSqrt1;
    sqrt1    = sqrt(2.0*forSqrt1);                   
    sqrt2    = sqrt(2.0*forSqrt2);                   
                     
    return complex<double>(m2*sqrt1-m1*sqrt2,0.0);
  }

  else if ( (l1==-1)&&(l2==+1)&&(s==-1) ){ 

    return complex<double>(0.0,0.0);
  }

  else if( (l1==-1)&&(l2==-1)&&(s==-1) ){ 

    inPr1    = -InProd_zero(k,-1,p1,+1);
    forSqrt1 = (p2[0]-p2[1])/(k[0]-k[1]);
    sqrt1    = sqrt(2.0*forSqrt1);                   
  
    return -inPr1*sqrt1;
  }
  else {     

    cout << "" << endl;
    cout << " ERROR IN TrMatrix_zero: " << endl;
    cout << "    WRONG VALUES FOR l1,l2,s" << endl;
    cout << "" << endl;             
    exit(0);
  }

}



////////////////////////////////////////////////////////////////
//          transition matrix for massive boson               //
//                                                            // 
//                                                            //
//                         \ eps(k,m,s)                       //
//                         /                                  // 
//                        _\                                  //
//                         /\ k                               // 
//                         \                                  //
//             <-- p1      /         <-- p2                   //                       
//           ---<----------\----------<---                    //
//       Ub(p1,m1,l1)                  U(p2,m2,l2)            //
//                                                            // 
////////////////////////////////////////////////////////////////                         
complex<double> TrMatrix_mass(double p1[4],double m1,int l1,double k[4],double m,int s,double p2[4],double m2,int l2){


  double forSqrt1,forSqrt2;
  double k_1[4],k_2[4];
  double forSqrt3,forSqrt4,sqrt3,sqrt1,sqrt2,sqrt4;
  complex<double>   inPr1,inPr2,inPr3,inPr4;

  //           double  spV[4],bet[4]
  //           double  pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af
  //           COMMON /Kleiss_Stirling/spV,bet
  //           COMMON /mc_parameters/pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN      
  // COMMON /mc_paraneters/pi,sw,cw,alphaI,qb,mb,mf1,mf2,qf1,qf2,vf,af,mcLUN 
  // temporary solution these will be global variables of the class   
  double pi     = mc_parameters_.pi;
  double alphaI = mc_parameters_.alphaI;
  double qb     = mc_parameters_.qb;
  double mf1    = mc_parameters_.mf1;
  double mf2    = mc_parameters_.mf2;
  double vf     = mc_parameters_.vf;
  double af     = mc_parameters_.af;
  double *spV = kleiss_stirling_.spV;
  // end of temporary solution

  for (int i = 0; i < 4; i++) {
    k_1[i] = 1.0/2.0*(k[i] - m*spV[i]);
    k_2[i] = 1.0/2.0*(k[i] + m*spV[i]);                                
  }

  if ( (l1==+1)&&(l2==+1)&&(s==0) ){ 
                
    inPr1 = InProd_zero(p1,+1,k_2,-1);
    inPr2 = InProd_zero(p2,-1,k_2,+1);
    inPr3 = InProd_zero(p1,+1,k_1,-1);
    inPr4 = InProd_zero(p2,-1,k_1,+1);
    sqrt1 = sqrt(p1[0]-p1[1]);
    sqrt2 = sqrt(p2[0]-p2[1]);
    sqrt3 = m1*m2/sqrt1/sqrt2;

              return                 
                            (inPr1*inPr2-inPr3*inPr4)*(vf+af)/m 
		+ (k_1[0]-k_2[0]-k_1[1]+k_2[1])*sqrt3*(vf-af)/m; 
  }       
                 
  else if ( (l1==+1)&&(l2==-1)&&(s==0) ){

    inPr1 = InProd_zero(p1,+1,k_1,-1);
    inPr2 = InProd_zero(p1,+1,k_2,-1);
    inPr3 = InProd_zero(p2,+1,k_2,-1);
    inPr4 = InProd_zero(p2,+1,k_1,-1);

    forSqrt1 = (k_1[0]-k_1[1])/(p2[0]-p2[1]);
    forSqrt2 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);
    forSqrt3 = (k_2[0]-k_2[1])/(p1[0]-p1[1]);
    forSqrt4 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);
    sqrt1 = sqrt(forSqrt1);
    sqrt2 = sqrt(forSqrt2);
    sqrt3 = sqrt(forSqrt3);
    sqrt4 = sqrt(forSqrt4);     

              return 
                  (inPr1*sqrt1 - inPr2*sqrt2)*(vf+af)*m2/m
		+ (inPr3*sqrt3 - inPr4*sqrt4)*(vf-af)*m1/m;
  }
  else if ( (l1==-1)&&(l2==+1)&&(s==0) ){ 

    inPr1 = InProd_zero(p1,-1,k_1,+1);
    inPr2 = InProd_zero(p1,-1,k_2,+1);
    inPr3 = InProd_zero(p2,-1,k_2,+1);
    inPr4 = InProd_zero(p2,-1,k_1,+1);

    forSqrt1 = (k_1[0]-k_1[1])/(p2[0]-p2[1]);
    forSqrt2 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);
    forSqrt3 = (k_2[0]-k_2[1])/(p1[0]-p1[1]);
    forSqrt4 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);
    sqrt1 = sqrt(forSqrt1);
    sqrt2 = sqrt(forSqrt2);
    sqrt3 = sqrt(forSqrt3);
    sqrt4 = sqrt(forSqrt4);     
        
              return 
                  (inPr1*sqrt1 - inPr2*sqrt2)*(vf-af)*m2/m
		+ (inPr3*sqrt3 - inPr4*sqrt4)*(vf+af)*m1/m;
  }
  else if  ( (l1==-1)&&(l2==-1)&&(s==0) ){ 

    inPr1 = InProd_zero(p2,+1,k_2,-1);
    inPr2 = InProd_zero(p1,-1,k_2,+1);
    inPr3 = InProd_zero(p2,+1,k_1,-1);
    inPr4 = InProd_zero(p1,-1,k_1,+1);
    sqrt1 = sqrt(p1[0]-p1[1]);
    sqrt2 = sqrt(p2[0]-p2[1]);
    sqrt3 = m1*m2/sqrt1/sqrt2;

             return                    
                         (inPr1*inPr2 - inPr3*inPr4)*(vf-af)/m  
	       + (k_1[0]-k_2[0]-k_1[1]+k_2[1])*sqrt3*(vf+af)/m;
  }
  else if ( (l1==+1)&&(l2==+1)&&(s==+1) ){ 

    inPr1 = InProd_zero(p1,+1,k_1,-1);
    inPr2 = InProd_zero(k_2,-1,p2,+1);
    inPr3 = inPr1*inPr2;

    forSqrt1 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);                       
    forSqrt2 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);  
    sqrt1 = sqrt(forSqrt1);                   
    sqrt2 = sqrt(forSqrt2);                   
    sqrt3 = m1*m2*sqrt1*sqrt2;

             return
	       sqrt(2.0)/m*(inPr3*(vf+af)+sqrt3*(vf-af));
  }

  else if ( (l1==+1)&&(l2==-1)&&(s==+1) ){

    inPr1 = InProd_zero(p1,+1,k_1,-1);
    inPr2 = InProd_zero(p2,+1,k_1,-1); 

    forSqrt1 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);                      
    forSqrt2 = (k_2[0]-k_2[1])/(p1[0]-p1[1]);                       
    sqrt1 = m2*sqrt(forSqrt1);                   
    sqrt2 = m1*sqrt(forSqrt2);                                     
                     
              return
                      sqrt(2.0)/m*( + inPr1*sqrt1*(vf+af)
                                    - inPr2*sqrt2*(vf-af)
				  );
  }
  else if  ( (l1==-1)&&(l2==+1)&&(s==+1) ){

    inPr1 = InProd_zero(k_2,-1,p2,+1);
    inPr2 = InProd_zero(k_2,-1,p1,+1);

    forSqrt1 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);                       
    forSqrt2 = (k_1[0]-k_1[1])/(p2[0]-p2[1]);                       
    sqrt1 = m1*sqrt(forSqrt1);                   
    sqrt2 = m2*sqrt(forSqrt2);                                     
                     
              return
                      sqrt(2.0)/m*( + inPr1*sqrt1*(vf+af)
                                    - inPr2*sqrt2*(vf-af)
				  );
  }
  else if ( (l1==-1)&&(l2==-1)&&(s==+1) ){ 

    inPr1 = InProd_zero(p2,+1,k_1,-1);
    inPr2 = InProd_zero(k_2,-1,p1,+1);
    inPr3 = inPr1*inPr2;

    forSqrt1 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);                       
    forSqrt2 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);  
    sqrt1 = sqrt(forSqrt1);                  
    sqrt2 = sqrt(forSqrt2);                   
    sqrt3 = m1*m2*sqrt1*sqrt2;

              return 
		sqrt(2.0)/m*(inPr3*(vf-af)+sqrt3*(vf+af));
  }

  else if ( (l1==+1)&&(l2==+1)&&(s==-1) ){ 

    inPr1 = InProd_zero(p2,-1,k_1,+1);
    inPr2 = InProd_zero(k_2,+1,p1,-1);
    inPr3 = inPr1*inPr2;

    forSqrt1 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);                       
    forSqrt2 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);  
    sqrt1 = sqrt(forSqrt1);                   
    sqrt2 = sqrt(forSqrt2);                   
    sqrt3 = m1*m2*sqrt1*sqrt2;

             return               
	       sqrt(2.0)/m*(inPr3*(vf+af)+sqrt3*(vf-af));
  }
  else if ( (l1==+1)&&(l2==-1)&&(s==-1) ){ 

    inPr1 = InProd_zero(k_2,+1,p2,-1);
    inPr2 = InProd_zero(k_2,+1,p1,-1);

    forSqrt1 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);                       
    forSqrt2 = (k_1[0]-k_1[1])/(p2[0]-p2[1]);                       
    sqrt1 = m1*sqrt(forSqrt1);                   
    sqrt2 = m2*sqrt(forSqrt2);                                     
                     
              return
                      sqrt(2.0)/m*(+ inPr1*sqrt1*(vf-af)
                                   - inPr2*sqrt2*(vf+af)
				  );
  }
  else if ( (l1==-1)&&(l2==+1)&&(s==-1) ){

    inPr1 = InProd_zero(p1,-1,k_1,+1);
    inPr2 = InProd_zero(p2,-1,k_1,+1);

    forSqrt1 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);                       
    forSqrt2 = (k_2[0]-k_2[1])/(p1[0]-p1[1]);                       
    sqrt1 = m2*sqrt(forSqrt1);                   
    sqrt2 = m1*sqrt(forSqrt2);                                     
                     
              return
                      sqrt(2.0)/m*(+ inPr1*sqrt1*(vf-af)
                                   - inPr2*sqrt2*(vf+af) 
				  );
  }
  else if ( (l1==-1)&&(l2==-1)&&(s==-1) ){ 

    inPr1 = InProd_zero(p1,-1,k_1,+1);
    inPr2 = InProd_zero(k_2,+1,p2,-1);
    inPr3 = inPr1*inPr2;

    forSqrt1 = (k_1[0]-k_1[1])/(p1[0]-p1[1]);                       
    forSqrt2 = (k_2[0]-k_2[1])/(p2[0]-p2[1]);  
    sqrt1 = sqrt(forSqrt1);                   
    sqrt2 = sqrt(forSqrt2);                   
    sqrt3 = m1*m2*sqrt1*sqrt2;

             return 
	       sqrt(2.0)/m*(inPr3*(vf-af)+sqrt3*(vf+af));
  }

  else{ 

    cout << " "<< endl;             
    cout << " TrMatrix_mass: Wrong values for l1,l2,s:"<< endl;
    cout << "          l1,l2 = -1,+1; s = -1,0,1 "<< endl;
    cout << " "<< endl;             
    exit(0);

  } 
         
}
