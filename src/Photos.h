#ifndef _Photos_h_included_
#define _Photos_h_included_

/** 
 * @class Photos
 *
 * @brief Controls the configuration, initialisation of Photos.
 *
 *
 * @author Nadia Davidson
 * @date 16th June 2008
 */

#include <stdarg.h>
#include <vector>
#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "Photos_make.h"

extern "C" {
  
  extern struct {
    double alpha;
    double xphcut;
  } phocop_;

  extern struct {
    double fsec;
    double fint;
    double expeps;
    int interf;
    int isec;
    int itre;
    int iexp;
    int iftop;
  } phokey_;

  extern struct {
    int iseed[2];
    int i97;
    int j97;
    double uran[97];
    double cran;
    double cdran;
    double cmran;
  } phseed_;

  extern "C" void phcork_(int * modcor);
  
}

class Photos{

 public:

   /** Initalise Photos with the parameters previously set via the
       setter methods */
   static void initialize();
   static void process(HepMC::GenEvent * event);
   static void suppressBremForDecay (int count, int motherID, ... );
   static void suppressBremForBranch(int count, int motherID, ... );

   static void setDebugPrintRange(int from, int to)  { phlupy_.ipoinm=from; phlupy_.ipoin=to; }

   //Wrappers for the PHOTOS configuration variables
   static void setInfraredCutOff(double cut_off){
     phocop_.xphcut=cut_off;};

   static void setAlphaQED(double alpha){
     phocop_.alpha=alpha;};
  
   static void setInterference(bool interference){
     phokey_.interf =(int) interference;}; 

   static void setDoubleBrem(bool doub){
     phokey_.isec =(int) doub;}; 

   static void setHigherBrem(bool higherBrem){
     phokey_.itre =(int) higherBrem;};

   static void setExponentiation(bool expo){
     phokey_.iexp = (int) expo; 
     if(expo){
       setDoubleBrem(false);
       setHigherBrem(false);       
       setInfraredCutOff(0.0000001);
       initializeKinematicCorrections(5);
       phokey_.expeps=0.0001;
     }
   };

   static void setTopProcessRadiation(bool top){
     phokey_.iftop = (int) top; };


   static void initializeKinematicCorrections(int flag){
     phcork_(&flag);};

   static std::vector< std::vector<int>* > *supBremList;
   static std::vector< HepMC::GenParticle* > *supParticles;
};

#endif  

