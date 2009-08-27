#include "Log.h"
#include <fstream>
using std::streambuf;
using std::stringstream;
using std::ostream;
using std::cout;
using std::cerr;
using std::endl;

streambuf   *Log::bCout=cout.rdbuf(),*Log::bCerr=cerr.rdbuf();
ostream     *Log::out=&cout;
stringstream Log::buf;
int  Log::dCount =0,Log::dRangeS =65535,Log::dRangeE =65534;
int  Log::faCount=0,Log::faRangeS=65535,Log::faRangeE=65534;
int  Log::iCount =0,Log::wCount =0,Log::eCount =0,Log::asCount=0, Log::asFailedCount=0;
bool Log::iAction=1,Log::wAction=1,Log::eAction=1,Log::asAction=1,Log::rAction=1;

ostream& Log::Debug(unsigned short int code, bool count)
{
	if(count) ++dCount;
	if(code>=dRangeS && code<=dRangeE ) return *out<<"DEBUG("<<code<<"):\t";
	return buf.seekp(0);
}


ostream& Log::Info(bool count)
{
	if(count) ++iCount;
	if(iAction) return *out<<"INFO:   \t";
	return buf.seekp(0);
}


ostream& Log::Warning(bool count)
{
	if(count) ++wCount;
	if(wAction) return *out<<"WARNING:\t";
	return buf.seekp(0);
}


ostream& Log::Error(bool count)
{
	if(count) ++eCount;
	if(eAction) return *out<<"ERROR:  \t";
	buf.seekp(0);
	return buf;
}

void Log::Assert(bool check, char *text)
{
	++asCount;
	if(check) return;
	++asFailedCount;
	if(text==NULL)	*out<<"ASSERT:\t\tAssertion failed. "<<endl;
	else *out<<"ASSERT:\t\tAssertion failed: "<<text<<endl;
	if(asAction) exit(-1);
}

void Log::Fatal(char *text,unsigned short code)
{
	++faCount;
	if(text==NULL) *out<<"TERM:\t\tTerminated by a call to Log::Exit();"<<endl;
	else *out<<"TERM:\t\t"<<text<<endl;
	if(code<faRangeS || code>faRangeE) exit(-1);
}

void Log::RedirectOutput(void (*func)(), ostream& where)
{

	if(!rAction) { func(); return; }
	cout.rdbuf(where.rdbuf());
	cerr.rdbuf(where.rdbuf());
	where<<endl;
	func();
	cout.rdbuf(bCout);
	cerr.rdbuf(bCerr);
}

void Log::RedirectOutput(ostream& where)
{
	if(!rAction) return;
	cout.rdbuf(where.rdbuf());
	cerr.rdbuf(where.rdbuf());
	where<<endl;
}

void Log::Summary()
{
	*out<<"-------------------------------- Log Summary ---------------------------------"<<endl;
	*out<<" Debug:   \t";
	if(dRangeS>dRangeE) *out<<"(OFF)";
	*out<<"\t\t"<<dCount<<"\t";
	if(dRangeS<=dRangeE) *out<<"Debug range: "<<dRangeS<<" - "<<dRangeE;
	*out<<endl;
	*out<<" Info:    \t";
	if(!iAction) *out<<"(OFF)";
	*out<<"\t\t"<<iCount<<"\t"<<endl;
	*out<<" Warnings:\t";
	if(!wAction) *out<<"(OFF)";
	*out<<"\t\t"<<wCount<<"\t"<<endl;
	*out<<" Errors:  \t";
	if(!eAction) *out<<"(OFF)";
	*out<<"\t\t"<<eCount<<"\t"<<endl;
	if(asCount || !asAction || faRangeS<faRangeE) cout<<"-----------------------------------"<<endl;
	if(asCount>0) *out<<" Asserts:\t\t\t"<<asCount<<endl;
	if(!asAction) *out<<" Failed asserts ignored:\t"<<asFailedCount<<endl;
	if(faRangeS<=faRangeE) *out<<" Fatal errors ignored:  \t"<<faCount<<endl;
	*out<<"------------------------------------------------------------------------------"<<endl;
}
