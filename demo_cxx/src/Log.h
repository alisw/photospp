#ifndef __LOG_CLASS_HEADER__
#define __LOG_CLASS_HEADER__
#include <iostream>
#include <sstream>
using std::stringstream;
using std::streambuf;
using std::ostream;

class Log
{
public:
	/** Shows the summary of all messages. */
	static void Summary();

	/** Shows the summary at the end of the program. */
	static void SummaryAtExit()              { atexit(Summary);      }

	/** Four logging entries. Usage:
            Log::Info()<<"Logging some info: "<<8<<" > "<<7.9<<endl;
	    Use Log::Info(false) if You don't want the message to be counted.*/
	static ostream& Debug(unsigned short int code=0, bool count=true);
	static ostream& Info(bool count=true);
	static ostream& Warning(bool count=true);
	static ostream& Error(bool count=true);

	/** Turns off or on particular types of messages
	    By default, only debugging messages are turned off. */
	static void LogInfo   (bool flag=true)  { iAction=flag;         }
	static void LogWarning(bool flag=true)  { wAction=flag;         }
	static void LogError  (bool flag=true)  { eAction=flag;         }

	static void LogAll    (bool flag=true)  { iAction=wAction=eAction=flag; dRangeS=0; dRangeE=65535; }

        /** Sets the range of debug codes that will be printed.
            By default, the debug messages are turned off. */
	static void LogDebug(unsigned short s=0,unsigned short e=65535)         { dRangeS=s; dRangeE=e;   }

	/** Asserts logical value. If the assertion fails, the default message or 'text'
            will be printed and the program will terminate.
            Program termination can be suppressed by Log::IgnoreFailedAsserts(); */
	static void Assert(bool check, char *text=NULL);

	/** Terminates the program with added default message or 'text'.
            It can be suppressed by Log::IgnoreFatal(); */
	static void Fatal(char *text, unsigned short int code=0);
	static void Fatal(unsigned short int code=0)                            { Fatal(NULL,code);       }

	/** Redirects output to log. Redirection can be done for a block of code
	    or for one function only. Redirection can be turned off by using
	    Log::IgnoreRedirection(); If the target is one of the log streams
	    (for example): Log::RedirectOutput( someFunction, Log::Info() );
	    You can turn the function's messages off by turning the apropriate
	    log entries off. The redirected code will still be executed,
	    only messages are redirected. */
	static void RedirectOutput(void (*func)(), ostream& where=*out);
	static void RedirectOutput(ostream& where=*out);
	/** WARNING! If You're redirecting more than one function, do not forget
	    to use RevertOutput() afterwards. */
	static void RevertOutput()                      { std::cout.rdbuf(bCout); std::cerr.rdbuf(bCerr); }

	/** Do not exit when Log::Assert() check is false.
	    The number of failed asserts will be listed in the summary. */
	static void IgnoreFailedAssert(bool flag=true)                           { asAction=!flag;        }

	/** Ignores redirections of functions' output.
	    The function will still be called in a normal way. */
	static void IgnoreRedirection(bool flag=true)                            { rAction=!flag;         }

	/** Do not exit when Log::Fatal() with the code within the provided range is called.
            The number of ignored fatal errors will be listed in the summary. */
	static void IgnoreFatal(unsigned short s=0,unsigned short e=65535) { faRangeS=s; faRangeE=e; }

	/** Change the output of the logged messages.
	    Log::SetOutput(cerr);                    //changes the output to cerr
	    Log::SetOutput(new ofstream("log.txt")); //changes the output to a file "log.txt" */
	static void SetOutput(ostream *newOut)                                    { out=newOut;           }
	static void SetOutput(ostream &newOut)                                    { out=&newOut;          }
protected:
	static streambuf *bCout,*bCerr;
	static ostream *out;
	static stringstream buf;
	static int  dCount,dRangeS,dRangeE,faCount,faRangeS,faRangeE;
	static int  iCount, wCount, eCount, asCount, asFailedCount;
	static bool iAction,wAction,eAction,asAction,rAction;
};
#endif
