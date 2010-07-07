/** 

 @mainpage C++ Interface to PHOTOS 
 @brief Description of PHOTOS Interface in C++

 @authors Nadia Davidson, Tomasz Przedzinski, Zbigniew Was

 @section download  release

  will be given here

 @section developement Developement version

 The source code and documentation are updated daily from the repository as well. The following files are provided for download of developement version:
 - <a href="../Photos_interface_design.pdf">Photos_interface_design.pdf</a>  full software documentation.
 - <a href="../PHOTOS.daily_temp.tar.gz">PHOTOS.daily_temp.tar.gz</a> tarball and its <a href="../svn_info_photos.txt">revision info</a> SVN tag, tarball creation date/time, etc.   The tar file contain the c++ interface along with 
 the source code for PHOTOS (see <a href="http://wasm.web.cern.ch/wasm/goodies.html">here
 </a>, version from Oct 11 2005).


 @section intro Introduction/Status
 At present (Jul 7) tar ball  is still not  complete, but these days
we work on it most of our time. 

Further content of doxygen web page is not updated since Feb 2 2010.
 

 @section setup Setup

 In order to run the interface and examples correctly you will need:
 - <a href="http://root.cern.ch/drupal/">ROOT v5.24</a>or later and <a href="http://lcgapp.cern.ch/project/simu/HepMC/">HepMC v2</a> or later
 - To run the PYTHIA example, you need <a href="http://home.thep.lu.se/~torbjorn/Pythia.html">PYTHIA 8.1</a> installed. PYTHIA 8.1 must be compiled with HepMC 2 so that the PYTHIA library hepmcinterface exists.
 - To run the TAUOLA+PYTHIA example, you need <a href="http://www.ph.unimelb.edu.au/~ndavidson/tauola/doxygen/index.html">TAUOLA C++ interface</a>.
 - You will also need <a href="http://mc-tester.web.cern.ch/MC-TESTER/">MC-TESTER</a> for the examples. Do not forget to type 'make libHepMCEvent' after compilation of MC-TESTER is done.
 - after downloading and uncompressing the interface, modify 'PHOTOS/setup.sh' to give the location of PYTHIA, HEPMC, TAUOLA and MC-TESTER.
   - PYTHIALOCATION should be the path to the base of the /include and /lib directories for PYTHIA 8.1
   - PYTHIA8DATA is the path to the directory containing PYTHIA xml documents (generally it should be "$(PYTHIA_INSTALL_LOCATION)/xmldoc").
   - HEPMCLOCATION should be the path to the base of the /include and /lib directories for HepMC
   - TAUOLALOCATION should be the path to the base of the /include and /lib directories for TAUOLA C++ interface
   - MCTESTERLOCATION should be the path to MC-TESTER (useful for validating the new interface)

 @section compile Compilation

 In order to compile the PHOTOS C++ interface:
 - modify setup.sh - set appropriate HEPMCLOCATION. If you plan to run examples, set PYTHIALOCATION, TAUOLALOCATION and MCTESTERLOCATION as well.
 - execute 'source setup.sh'
 - check if 'platform/make.inc' points to correct version of your platform and change it if necessary.
 - execute 'make'

 The '/lib' and '/include' directories will contain the appropriate library and include files.

 In order to compile the examples:
 - Compile PHOTOS C++ interface
 - Verify that you have both PYTHIA and MC-TESTER installed and compiled.
 - If you haven't done it yet, modify setup.sh setting all four variables and execute 'source setup.sh' again.
 - Enter 'PHOTOS/example' directory and execute 'make'
 - An example of PHOTOS C++ interface compiled along with TAUOLA C++ interface is provided as well. It can be compiled by changing the "MAIN" variable in the Makefile and executing 'make'.

 The '/example' directory will contain the compiled example files.

 @section testing Testing

 In order to run some more specific tests:
 - Compile PHOTOS C++ interface as well as examples. For all test to work you will need both aviable example files so modify 'Makefile' inside PHOTOS/example' directory, switching MAIN variable to compile both files. After both files are compiled enter 'PHOTOS/example/testing' directory.
 - modify test.inc if needed.
 - enter choosen directory and execute 'make'.

 The appropriate .root files as well as .pdf files generated by MC-TESTER will be created inside the choosen directory. You can execute 'make clobber' to clean the directory. You can also execute 'make' inside 'PHOTOS/example/testing' directory to run all aviable tests one after another.

 <hr>
 Expect more information to be added at a later date.

*/
