#ifndef _f_Init_included_
#define _f_Init_included_

extern "C"
{

	extern struct
	{
		double alpha;
		double xphcut;
	} phocop_;

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
	} phokey_;

	extern struct
	{
		int iseed[2];
		int i97;
		int j97;
		double uran[97];
		double cran;
		double cdran;
		double cmran;
	} phseed_;

	//debug mode on if ipoin <  1 and ipoinm > 1
	extern struct
	{    
		int ipoin;
		int ipoinm;
	} phlupy_;  

	extern void phcork_(int * modcor);

	//extern void dexay_(int *state, double pol[4]);
	extern void phoini_(); //PHOTOS initialisation
	extern void photos_make_(int * id);
}

#endif
