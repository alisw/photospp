/** 

 @mainpage C++ Interface to PHOTOS 
 @brief Description of PHOTOS Interface in C++

 @authors Nadia Davidson, Tomasz Przedzinski, Zbigniew Was

 
 
 @section download1 New release

 The source code and documentation for release 3.3. The following files are provided for download:
 - <a href="http://hibiscus.if.uj.edu.pl/~przedzinski/PHOTOS.3.3/Photos_interface_design.3.3.pdf">Photos_interface_design.pdf</a> full software documentation.
 - <a href="http://hibiscus.if.uj.edu.pl/~przedzinski/PHOTOS.3.3/PHOTOS.3.3.tar.gz">PHOTOS 3.3 source code </a> tarball
   and its <a href="http://hibiscus.if.uj.edu.pl/~przedzinski/PHOTOS.3.3/svn_info_photos.3.3.txt">revision info</a> SVN tag, tarball creation date/time, etc.
   For updates with respect to release 3.0 see <a href="http://hibiscus.if.uj.edu.pl/~przedzinski/PHOTOS.3.3/changelog.3.3.txt">changelog.txt</a>

 @section developement Developement version

 The source code and documentation are updated daily from the repository. The following files are provided for download of the developement version:
 - <a href="../Photos_interface_design.pdf">Photos_interface_design.pdf</a>  full software documentation.
 - <a href="../PHOTOS.daily_temp.tar.gz">PHOTOS.daily_temp.tar.gz</a> tarball and its <a href="../svn_info_photos.txt">revision info</a> SVN tag, tarball creation date/time, etc.   The tar file contains the c++ interface along with parts of
 the source code for PHOTOS (see <a href="http://wasm.web.cern.ch/wasm/goodies.html">here
 </a>, version from Oct 11 2005). For updates  with respect to release 3.0 see <a href="../changelog.txt">changelog.txt</a>.
 - NEW: NLO in decays of Z, W and B (K) --> 2 scalar, <a href="http://annapurna.ifj.edu.pl/~wasm/phNLO.htm">photos NLO tests</a>.

 @section download Older releases

 The source code and documentation for release 3.0. The following files are provided for download:
 - <a href="http://arxiv.org/abs/1011.0937">arXiv:1011.0937</a> full software documentation.
 - <a href="../PHOTOS.3.0.tar.gz">PHOTOS 3.0 source code </a> tarball.

 Note that LCG/Genser
 <a href="http://sftweb.cern.ch/generators/">Generator 
 Services Subproject </a> distributes compiled, platform adopted  
 tar balls of our programs.

 @section intro Introduction/Status
 The tar-ball has the functionality of its FORTRAN predecessor for all cases.
 An extensive number of tests performed in X.2011 and XI.2011 has been collected on the webpage of <a href="http://annapurna.ifj.edu.pl/~wasm/phNLO.htm">photos NLO tests</a>.
 In particular, tests with SANC went at sub 0.01\% level
  
 @section setup Requirements

 For compilation, and to run the simple example, the interface requires:
 - <a href="http://lcgapp.cern.ch/project/simu/HepMC/">HepMC v2.04</a> or later.

 For a multitude of examples, one needs to install the libraries necessary to
generate physics events stored in HepMC and to monitor such events:
 - <a href="http://home.thep.lu.se/~torbjorn/Pythia.html">PYTHIA 8.1</a> or later. PYTHIA must be compiled with HepMC 2 so that the PYTHIA library hepmcinterface exists.
 - <a href="http://www.ph.unimelb.edu.au/~ndavidson/tauola/doxygen/index.html">TAUOLA C++ Interface v1.0</a> or later.
 - <a href="http://mc-tester.web.cern.ch/MC-TESTER/">MC-TESTER v1.24</a> or later. Do not forget to compile the additional HepMC library libHepMCEvent as well.
 - <a href="http://root.cern.ch/drupal/">ROOT v5.18</a> or later

 @section compile Configuration and Compilation

 In order to compile the PHOTOS C++ interface:
 - Execute './configure' with additional command line options:
   - '--with-hepm=\<path\> ' provides the path to the HepMC installation directory. One can set the HEPMCLOCATION variable instead of using this directive. This path is required for the interface to compile. To compile without HepMC use '--without-hepmc'.
   - '--prefix=\<path\>' provides the installation path. The 'include' and 'lib' directories will be copied there if 'make install' is executed later. If none has been provided, the default directory for installation is '/usr/local'.
 - Execute 'make'
 - Optionally, execute 'make install' to copy files to the directory provided during configuration.

 The PHOTOS C++ interface will be compiled and the '/lib' and '/include' directories will contain the appropriate library and include files.

 In order to compile the examples, enter 'examples' directory, and:
 - execute './configure' to determine which examples can be compiled. Additional paths can be provided as command line options:
   - '--with-pythia8=\<path\>' provides the path to the Pythia8 installation directory. One can set the PYTHIALOCATION variable instead of using this directive. This path is required for all examples and tests.
   - '--with-mc-tester=\<path\>' provides the path to the MC-Tester installation directory (the libHepMCEvent must be compiled as well, check the MC-Tester documentation for more details). One can set the MCTESTERLOCATION variable instead of using this directive. This path is required for all additional examples and tests. It is assumed that using this option also implies that ROOT has already been installed (since it's required by MC-TESTER). The location of its binaries should be listed in the PATH variable.
   - '--with-tauola=\<path\>'  provides the path to the TAUOLA C++ interface installation directory. One can set the TAUOLALOCATION variable instead of using this directive. This path is required for additional examples.
 - execute 'make'

 Note that for examples working with PYTHIA 8.1, the PYTHIA8DATA global variable must be set (refer to the instructions provided during configuration).
 Similarly, for examples in the examples/testing directory to work, the MCTESTERLOCATION global variable must be set.
 If neither PYTHIA nor MC-TESTER are available, only the simple example can be 
used. The '/examples' directory will contain the compiled example files.

 @section testing Testing

 In order to run  more elaborate and physics interesting tests both PYTHIA 
and MC-TESTER must be installed. In some cases TAUOLA C++ will be needed too.
 - Compile the PHOTOS C++ interface as well as examples.
 - Check that the appropriate system variables are set: normally set by the script
 configure.paths.sh [.csh] (the configuation step mentions this script).
 - Enter the /examples/testing directory. Modify test.inc if needed.
 - Enter the selected directory and execute 'make'.

 The appropriate .root files as well as .pdf files generated by MC-TESTER will be created inside the choosen directory. You can execute 'make clobber' to clean the directory. You can also execute 'make' inside the 'PHOTOS/examples/testing' directory to run all available tests one after another.



 <hr>
 Last update 07 December 2011.

*/
