INSTALLATION PROCEDURE FOR CODE_SATURNE
=======================================
=======================================

For more information about the different modules and external libraries
necessary or compliant with Code_Saturne, refer to the COMPATIBILITY file
in the Kernel distribution (ncs-x.y.z).

Section I   gives information on the automatic installer
Section II  gives information for manual installation
Section III gives information for using the code


I) AUTOMATIC INSTALL
====================
The Install directory contains a python script for automatic
installation of the Code_Saturne elements and associated routines.
In most cases, it will be enough. In case of problems, switch to
section II for element by element install.
These scripts are given in the hope that they will be useful, but
WITHOUT ANY WARRANTY.

The script can download every packages needed by the code to run
properly. If this behaviour is not wanted, set the "download" variable
to "no" in the setup script.

It is possible to specify her/his own Python executable through the "python"
variable, as well as BLAS libraries through the "blas" variable, or Metis and
Scotch through similar variables.

SYRTHES installation path (for coupling with the thermal code SYRTHES) will be
provided with the "syrthes" variable.

Lastly, the possibility is given to compile Code_Saturne with debugging symbols
("debug" variable), to disable the Graphical User Interface ("disable_gui"
variable), and to specify the language (between English and French).

On some architectures and for some elements (MED and FVM for instance)
it is preferable if the "make" command refers to the GNU "make". Otherwise
some problems can occur (problem with libtool, need to copy the sources in
the build directory, ...)

* install_saturne.py:
  This python script will install the different elements of Code_Saturne and
  associated libraries. Due to dependencies between the different modules, the
  order of install should be the following:
  - libxml2 (it is advised to use the distrib own package)
  - swig (it is advised to use the distrib own package)
  - MPI
  - CGNS
  - HDF5
  - MED
  - BFT (Code_Saturne Base Functions and Types library)
  - FVM (Code_Saturne Finite Volume Mesh library)
  - MEI (Code_Saturne Mathematical Expressions Interpreter library)
  - ECS (Code_Saturne Preprocessor)
  - NCS (Code_Saturne Kernel and Graphical User Interface)

  The following packages cannot be installed
  - Zlib
  - BLAS
  - Metis
  - Scotch
  - PyQT
  - Python


  The install script uses the "setup" file to determine which library to
  install or to use. For each element, there are four options:

  - to not use the element (for optional libraries like MPI)
     In this case, specify "no" in the "Usage" and "Install" columns. The other
     elements will be installed in accordance. The "Path" column is not used.

  - to automatically detect some element (especially useful for libxml2 or swig)
     In this case, specify "auto" in the "Usage". The other elements will be
      installed in accordance. The "Path" and "Install" column are not used.

  - to use a pre-installed library in a non standard path
     In this case, specify "yes" in the Usage column and "no" in the Install
     column. The "Path" column should contain the location of the library
     (up to the name of the library itself).

  - to install and use a library
     In this case, specify "yes" in the "Usage" and "Install" columns. The
     script will download the library and install it default install directory.
     If download has been set to "no", package archive are looked for at the
     same location than the installation script (the right number and archive
     name are needed, accordingly to what is prescribed in the script).
     After each element has been installed, the "setup" file is modified, the
     column "Install" of the concerned element is set to "no" and the "Path"
     column is filled so that the element is not installed a second time if
     the script is relaunched (if there was a problem with a later element).

   Before using the "install_saturne.py" script, the C and Fortran compilers
   to be used can be specified next to the CompC and CompF keywords.
   An optional MPI wrapper compiler (for the C language) can be specified to
   be used for FVM and the Kernel installation.

   If the "use_arch" variable is set to "yes", then the "arch" keyword refers
   to the architecture of the machine. Leaving it blank will make it
   automatically detected with the "uname" command."arch" should be specified
   if you want different implementations on the same architecture
   (for instance Linux_LAM and Linux_MPICH).

   Commonly used compilers for different architectures :
   SunOS : compC = cc -Xa
           compF = f90
   
   IRIX64 : compC = cc -64
            compF = f90 -64



    
II) MANUAL INSTALL
==================
If the automatic install script fails, you should install each module
individually. Except for the modules developped by EDF (Kernel, Preprocessor,
BFT, FVM and MEI), you should refer to the installation procedure provided by
the distributor of the module.

Herebelow, $version stands for the current version of Code_Saturne and $prefix
for the directory where you want to install Code_Saturne.


  Installing BFT
  --------------
  - create a "build" directory (usually bft-x.y.z.build)

  - from within the build directory, run the configure command:

      ../bft-x.y.z/configure --prefix=$prefix/cs-$version

    Depending on the elements wanted, additional options can be added
    to the configure:
      --with-zlib=... for Zlib support
      CC=...          to specify the compiler if necessary

  - run the "make" commands:
      make
      make install
      make clean

    For further information, refer to the INSTALL file in the BFT directory.


  Installing FVM
  --------------
  - create a "build" directory (usually fvm-x.y.z.build)

  - from within the build directory, run the configure command:

      ../fvm-x.y.z/configure --prefix=$prefix/cs-$version \
         --with-bft=$prefix/cs-$version

    Depending on the elements wanted, additional options can be added
    to the configure:
      --with-cgns=... for CGNS support
      --with-hdf5=... for HDF5 (compulsory for MED)
      --with-med=...  for MED
      --with-mpi=...  for MPI
      CC=...          to specify the compiler
                         (especially if mpicc should be used, to get
                          the proper links to MPI libraries)

  - run the "make" commands:
      make
      make install
      make clean

    For further information, refer to the INSTALL file in the FVM directory.


  Installing MEI
  --------------
  - create a "build" directory (usually mei-x.y.z.build)

  - from within the build directory, run the configure command:

      ../mei-x.y.z/configure --prefix=$prefix/cs-$version \
         --with-bft=$prefix/cs-$version

    Depending on the elements wanted, additional options can be added
    to the configure:
      --with-python-exec=...  for specific Python executable
      --with-swig-exec=...    for SWIG (Python bindings)

  - run the "make" commands:
      make
      make install
      make clean

    For further information, refer to the INSTALL file in the MEI directory.


  Installing the Preprocessor (ecs)
  ---------------------------------
  - create a "build" directory (usually ecs-x.y.z.build)

  - from within the build directory, run the configure command:

      ../ecs-x.y.z/configure --prefix=$prefix/cs-$version \
         --with-bft=$prefix/cs-$version

    Depending on the elements wanted, additional options can be added
    to the configure:
      --with-adf=...    for ADF support (compulsory for CCM if no CGNS)
      --with-ccm=...    for CCM support
      --with-cgns=...   for CGNS support
      --with-hdf5=...   for HDF5 (compulsory for MED)
      --with-med=...    for MED
      --with-metis=...  for Metis optimised domain partitioning
                            (strongly advised for parallel computing)
      --with-scotch=... for Scotch optimised domain partitioning
                           (alternative for Metis)
      CC=...            to specify the compiler if necessary

  - run the "make" commands:
      make
      make install
      make clean

    For further information, refer to the INSTALL file in the ECS directory.


  Installing the Kernel (ncs)
  ---------------------------

  - create a "build" directory (usually ncs-x.y.z.build)

  - from within the build directory, run the configure command:

      ../ncs-x.y.z/configure --prefix=$prefix/cs-$version \
         --with-bft=$prefix/cs-$version \
         --with-fvm=$prefix/cs-$version \
         --with-mei=$prefix/cs-$version \
         --with-prepro=$prefix/cs-$version

    Depending on the elements wanted, additional options can be added
    to the configure:
      --with-libxml2=...      for Libxml2 support (compulsory for the Interface)
      --with-blas=...         for BLAS
      --with-mpi=...          for MPI
      --with-syrthes=...      for SYRTHES coupling
      --with-python-exec=...  for specific Python executable
      --with-pyqt4-exec=...   for PyQt4 developper tools
      CC=...           to specify the compiler if necessary
                         (especially if mpicc should be used, to get
                          the proper links to MPI libraries)
      FC=...           to specify the Fortran compiler if necessary

  - run the "make" commands:
      make
      make install
      make clean


  - documentation can be compiled (if LaTeX is available) and installed
    by running:
      make pdf
      make install-pdf

    The compiled libraries will be put in $prefix/cs-$version
    


III) Before using Code_Saturne
==============================
Each user of Code_Saturne must set her/his PATH accordingly with Code_Saturne
installation before using the code. The easiest way is to put the following
lines in each of the users ".profile" (depending on the shell).

cspath=$prefix/cs-$version/bin
#(adjust path to your system)
if [ -d $cspath ] ; then  
  export PATH=$cspath:$PATH
fi

After changing the user ".profile", it is advised to logout and login,
so that there is no mix-up in the PATH variable.

For more information refer to the Code_Saturne documentation, available
through the "code_saturne info -g refcard" and "code_saturnes info -g user"
commands.

Code_Saturne support: saturne-support@edf.fr
