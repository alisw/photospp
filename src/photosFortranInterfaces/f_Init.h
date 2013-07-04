#ifndef _f_Init_included_
#define _f_Init_included_
const static int NMXHEP = 10000;
extern "C"
{

	/** Definition of the PHOEVT common block */
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
		//      NEVPHO,NPHO,ISTPHO(NMXPHO),IDPHO(NMXPHO),
		//  JMOPHO(2,NMXPHO),JDAPHO(2,NMXPHO),PPHO(5,NMXPHO),VPHO(4,NMXPHO)
	} phoevt_;

	extern struct
	{
		double alpha;
		double xphcut;
	} phocop_;

	extern struct
	{
		double iev;
	} phnum_;



	extern struct
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

	extern struct PHOSTA
	{
		int status[10];
	        int ifstop;
	} phosta_;

	extern struct
	{
		double pi;
		double twopi;
	} phpico_;

	extern struct
	{
		int phlun;
	} pholun_;

	extern struct
	{
		double xphmax;
		double xphoto;
		double costhg;
		double sinthg;

	} phophs_;

	extern struct
	{
		double probh;
		double corwt;
		double xf;
		int irep;
	} phopro_;

	extern struct
	{
		double fi3;
		double fi1;
		double th1;
		int irep;

	} phorest_;

	extern struct
	{
		double beta;
		double wt1;
		double wt2;
		double wt3;

	} phwt_;

	extern struct
	{
		double mchsqr;
		double mnesqr;
		double pneutr[5];
	} phomom_;

	//debug mode on if ipoin <  1 and ipoinm > 1
	extern struct PHLUPY
	{
		int ipoin;
		int ipoinm;
	} phlupy_;

	/** Initialize kinematic corrections */
	void phcork_(int * modcor);

	/** Single branch processing */
	void photos_make_c_(int * id);

	/* PHOINI subroutines */
	int  iphqrk_(int *i);
	int  iphekl_(int *i);

	/* Printout of error messages */
	void phoerr_(int *imes,char *text,double *data);

}

#endif
