/** 

 @mainpage C++ Interface to PHOTOS 
 @brief Description of PHOTOS Interface in C++

 @authors Nadia Davidson, Tomasz Przedzinski, Zbigniew Was

 @section download  release

  will be given here

 @section developement Developement version

 The source code and documentation are updated daily from the repository. The following files are provided for download of developement version:
 - <a href="../Photos_interface_design.pdf">Photos_interface_design.pdf</a>  full software documentation.
 - <a href="../PHOTOS.daily_temp.tar.gz">PHOTOS.daily_temp.tar.gz</a> tarball and its <a href="../svn_info_photos.txt">revision info</a> SVN tag, tarball creation date/time, etc.   The tar file contain the c++ interface along with 
 the source code for PHOTOS (see <a href="http://wasm.web.cern.ch/wasm/goodies.html">here
 </a>, version from Oct 11 2005).

 @section intro Introduction/Status
 At present (Jul 8) tar ball  is still not  complete but
 has functionality  of FORTRAN predecessor for all cases except t bar t  production.
 Program is under tests and extension (main work direction for July).
The purpose is to extend functionality to all features present in FORTRAN.
Introduce better control of options and finally get better user
controllable steering.
 
 @section setup Requirements

 For compilation, and to run simple example, the interface requires:
 - <a href="http://lcgapp.cern.ch/project/simu/HepMC/">HepMC v2.04</a> or later.

 For multitude of examples, one need to install libraries necessary to
generate physics events stored in HepMC and to monitor such events:
 - <a href="http://home.thep.lu.se/~torbjorn/Pythia.html">PYTHIA 8.1</a> or later. PYTHIA must be compiled with HepMC 2 so that the PYTHIA library hepmcinterface exists.
 - <a href="http://www.ph.unimelb.edu.au/~ndavidson/tauola/doxygen/index.html">TAUOLA C++ Interface v1.0</a> or later.
 - <a href="http://mc-tester.web.cern.ch/MC-TESTER/">MC-TESTER v1.24</a> or later. Do not forget to compile the additional HepMC library libHepMCEvent as well.
 - <a href="http://root.cern.ch/drupal/">ROOT v5.18</a> or later

 @section compile Configuration and Compilation

 In order to compile the PHOTOS C++ interface:
 - Execute './configure' with additional command line options:
   - '--with-HepMC=\<path\> ' provides the path to HepMC installation directory. One can set HEPMCLOCATION variable instead of using this directive. This path is required for interface to compile.
   - '--prefix=\<path\>' provides the installation path. The 'include' and 'lib' directories will be copied there if 'make install' is executed later. If none has been provided, the default directory for installation is '/usr/local'.
 - Execute 'make'
 - Optionally, execute 'make install' to copy files to the directory provided during configuration.

 The PHOTOS C++ interface will be compiled and the '/lib' and '/include' directories will contain the appropriate library and include files.

 In order to compile the examples, enter 'examples' directory, and:
 - execute './configure' to determine which examples can be compiled. Additional paths can be provided as command line options:
   - '--with-Pythia8=\<path\>' provides the path to Pythia8 installation directory. One can set PYTHIALOCATION variable instead of using this directive. This path is required for all examples and tests.
   - '--with-MC-Tester=\<path\>' provides the path to MC-Tester installation directory (the libHepMCEvent must be compiled as well, check MC-Tester documentation for more details). One can set MCTESTERLOCATION variable instead of using this directive. This path is required for all additional examples and tests. It is assumed that using this option also implies that ROOT has already been installed (since it's required by MC-TESTER). The location of its binaries should be listed in PATH variable.
   - '--with-Tauola=\<path\>'  provides the path to TAUOLA C++ interface installation directory. One can set TAUOLALOCATION variable instead of using this directive. This path is required for additional examples.
 - execute 'make'

 Note that for examples working with PYTHIA 8.1, PYTHIA8DATA global variable must be set (refer to instructions provided during configuration).
 Similar, for examples in examples/testing directory to work, MCTESTERLOCATION global variable must be set.
 If neither PYTHIA nor MC-TESTER are available, only  simple example can be 
used. The '/examples' directory will contain the compiled example files.

 @section testing Testing

 In order to run  more elaborate and physics interesting tests both PYTHIA 
and MC-TESTER must be installed. In some cases TAUOLA C++ will be needed too.
 - Compile PHOTOS C++ interface as well as examples.
 - Check that appropriate system variables are set: normally set by script
 configure.paths.sh [.csh] (configuation step is mentioning this script).
 - Enter /examples/testing directory. Modify test.inc if needed.
 - Enter selected directory and execute 'make'.

 The appropriate .root files as well as .pdf files generated by MC-TESTER will be created inside the choosen directory. You can execute 'make clobber' to clean the directory. You can also execute 'make' inside 'PHOTOS/examples/testing' directory to run all aviable tests one after another.



 <hr>
 Last update 9 July 2010.
 Expect more information to be added at a later date (here in the doxygen 
web page and in software documentation).

*/
