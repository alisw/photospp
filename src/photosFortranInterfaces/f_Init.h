#ifndef _f_Init_included_
#define _f_Init_included_

namespace Photospp
{

const static int NMXHEP = 10000;
const static double PI    = 3.14159265358979324;
const static double TWOPI = 6.28318530717958648;


// extern "C"
//{

	/** Definition of the PHOEVT common block */
	extern struct HEPEVT
	{
		int    nevhep;
		int    nhep;
		int    isthep[NMXHEP];
		int    idhep[NMXHEP];
		int    jmohep[NMXHEP][2];
		int    jdahep[NMXHEP][2];
		double phep[NMXHEP][5];
		double vhep[NMXHEP][4];
	        int    qedrad[NMXHEP];  // to be bool once compatibility with F77 removed       
	        const static int nmxhep=NMXHEP;
		//      NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
		//  JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
	        //   int qedrad[NMXHEP]  was an add up 
                //   for  HEPEVT in F77 times. Separate common PH_PHOQED
	        //   also phoif_.chkif[NMXPHO] was add up for PHOEVT
                //   now it is pho.qedrad
	} hep,pho;
        //ph_hepevt_,phoevt_;



	extern struct PHOCOP
	{
		double alpha;
		double xphcut;
	} phocop_;

	struct
	{
		double iev;
	} phnum;



	extern struct PHOKEY
	{
		double fsec;
		double fint;
		double expeps;
		int interf;
		int isec;
		int itre;
		int iexp;
		int iftop;
		int ifw;
	} phokey_;

	struct
	{
		int status[10];
	        int ifstop;
	} phosta;

	struct
	{
		int phlun;
	} pholun;

	struct
	{
		double xphmax;
		double xphoto;
		double costhg;
		double sinthg;

	} phophs;
	struct
	{
		double QQ[4];
		double XM;
		double th1;
		double fi1;

	} tofrom;

	struct
	{
		double probh;
		double corwt;
		double xf;
		int irep;
	} phopro;

	struct
	{
		double fi3;
		double fi1;
		double th1;
        	int irep;     //    provably line to be removed

	} phorest;

	struct
	{
		double beta;
		double wt1;
		double wt2;
		double wt3;

	} phwt;
	struct
	{
		double phocorwt3;
		double phocorwt2;
		double phocorwt1;

	} phocorwt;

	struct
	{
		double mchsqr;
		double mnesqr;
		double pneutr[5];
	} phomom;
	struct
	{
		double bet[3];
		double gam;
	} phocms;

	struct
	{
	        const static int NX = 10;
		double pro[NX];
		int nchan;
	        int expini;    // bool
	} phoexp;

	//debug mode on if ipoin <  1 and ipoinm > 1
	struct
	{
		int ipoin;
		int ipoinm;
	} phlupy;

	/** Initialize kinematic corrections */
	void PHCORK(int modcor);

	/** Single branch processing */
	void PHOTOS_MAKE_C(int id);

	/* Central management routine. Defines what action
	   will be performed at point ID. */
	void PHTYPE(int ID);
// }

} // namespace Photospp
#endif
