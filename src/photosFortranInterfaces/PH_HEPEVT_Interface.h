#ifndef _PH_HEPEVT_Interface_included_
#define _PH_HEPEVT_Interface_included_

#include <vector>
#include "PhotosBranch.h"
#include "PhotosParticle.h"

const static int NMXHEP = 10000;

extern "C"
{
	/** Definition of the PH_HEPEVT common block */
	extern struct
	{
		int    nevhep;
		int    nhep;
		int    isthep[NMXHEP];
		int    idhep[NMXHEP];
		int    jmohep[NMXHEP][2];
		int    jdahep[NMXHEP][2];
		double phep[NMXHEP][5];
		double vhep[NMXHEP][4];
		//      COMMON/PH_HEPEVT/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
		//  &JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
	} ph_hepevt_;

	extern struct
	{ 
		int qedrad[NMXHEP];
	} ph_phoqed_;

	/** Fortran debug routine */
	extern "C" void phodmp_();
}

class PH_HEPEVT_Interface
{
public:
	/** Convert PhotosBranch to HEPEVT */
	static int  set(PhotosBranch* branch);
	/** Update event record with data from HEPEVT */
	static void get();
	/** Clear HEPEVT */
	static void clear();
private:
	/** Add single particle to HEPEVT */
	static void add_particle(int i, PhotosParticle * particle,
	                         int first_mother, int last_mother,
	                         int first_daughter, int last_daughter);

	/** List of particles added to HEPEVT */
	static std::vector<PhotosParticle*> m_particle_list;
};

#endif
