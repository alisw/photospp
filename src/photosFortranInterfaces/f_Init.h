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

	/** Initialize kinematic corrections */
	extern void phcork_(int * modcor);

	/** PHOTOS initialization */
	extern void phoini_();
	/** Single branch processing */
	extern void photos_make_c_(int * id);
}

#endif
