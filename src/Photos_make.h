#ifndef _f_Photos_make_included_
#define _f_Photos_make_included_

#include <vector>
#include <iostream>

extern "C" {

  //debug mode on if ipoin <  1 and ipoinm > 1
  extern struct {    
    int ipoin;
    int ipoinm;
  } phlupy_;  


  //extern void dexay_(int *state, double pol[4]);
  extern void phoini_(); //PHOTOS initialisation
  extern void photos_(int * id);
  extern void photos_get_();
  extern void photos_set_();
  extern void photos_make_(int * id);


}

#endif
