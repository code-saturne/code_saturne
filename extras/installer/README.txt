INSTALLATION PROCEDURE FOR CODE_SATURNE
=======================================
=======================================

For more information about the different modules and external libraries
necessary or compliant with Code_Saturne, refer to the COMPATIBILITY file
in the Kernel distribution (code_saturne-x.y.z).

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

The script can download every package needed by the code to run
properly. If this behaviour is not wanted, set the "download" variable
to "no" in the setup script.

It is possible to specify a Python interpreter through the "python"
variable, as well as Metis and Scotch through similar variables.

Lastly, the possibility is given to compile Code_Saturne with debugging symbols
("debug" variable), to disable the Graphical User Interface ("disable_gui"
variable), and to specify the language (between English and French).

On some architectures and for some elements (MED and SCOTCH for instance)
it is preferable if the "make" command is a recent enough version of GNU "make".
Otherwise some problems can occur.

* install_saturne.py:
  This python script will install the different elements of Code_Saturne and
  associated libraries. Due to dependencies between the different modules, the
  order of install should be the following:
  - libxml2 (it is advised to use the distrib's own package)
  - MPI
  - CGNS
  - HDF5
  - MED

  The following packages cannot be installed
  - Zlib
  - Metis
  - Scotch
  - PyQT
  - Python


  The install script uses the "setup" file to determine which library to
  install or to use. For each element, there are four options:

  - do not use the element (for optional libraries like MPI)
     In this case, specify "no" in the "Usage" and "Install" columns. The other
     elements will be installed in accordance. The "Path" column is not used.

  - automatically detect some element (especially useful for libxml2)
     In this case, specify "auto" in the "Usage". The other elements will be
      installed in accordance. The "Path" and "Install" column are not used.

  - use a pre-installed library in a non standard path
     In this case, specify "yes" in the Usage column and "no" in the Install
     column. The "Path" column should contain the location of the library
     (up to the name of the library itself).

  - install and use a library
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
   be used Code_Saturne installation.

   If the "use_arch" variable is set to "yes", then the "arch" keyword refers
   to the architecture of the machine. Leaving it blank will make it
   automatically detected with the "uname" command."arch" should be specified
   if you want different implementations on the same architecture
   (for instance Linux_OMPI and Linux_MPICH).


II) MANUAL INSTALL
==================
If the automatic install script fails, you should install each module
individually. For all modules, except Code_Saturne, you should refer
to the installation procedure provided by the distributor of the module.

Herebelow, $version stands for the current version of Code_Saturne and $prefix
for the directory where you want to install Code_Saturne.

  Installing Code_Saturne
  -----------------------

  - create a "build" directory (usually code_saturne-x.y.z.build)

  - from within the build directory, run the configure command:

    ../code_saturne-x.y.z/configure --prefix=$prefix/code_saturne-$version \

      # Pre- and post-processing *optional* support

      --with-zlib=...       for Zlib support
      --with-hdf5=...       for HDF5 support (compulsory for MED)
      --with-med=...        for MED support
      --with-cgns=...       for CGNS support
      --with-adf=...        for ADF support (compulsory for CCM if no CGNS)
      --with-ccm=...        for CCM support

      # Partitioning *optional* support (strongly advised for parallel computing)

      --with-scotch=...     for Pt-SCOTCH / SCOTCH support
      --with-metis=...      for ParMETIS / METIS (an alternative to SCOTCH)

      # Graphical interface and scripts support

      --with-libxml2=...    for Libxml2 support
      PYTHON=...            for specific Python executable
      PYUIC4=...            for PyQt4 developer tools pyuic4
      PYRCC4=...            for PyQt4 developer tools pyrcc4

      # Run-time options

      --with-mpi=...        for MPI (parallel computing)

      CC=...   to specify the compiler if necessary (especially if mpicc
                  should be used, to get the proper links to MPI libraries)
      FC=...   to specify the Fortran compiler if necessary


  - run the "make" commands:
      make
      make install
      make clean

    The compiled libraries will be put in $prefix/cs-$version


III) Before using Code_Saturne
==============================
For some systems (such as when using a batch system or coupling with SYRTHES,
a post-install step may be required). In this case, copy
"$prefix/code_saturne-$version/etc/code_saturne.cfg.template" to
"$prefix/code_saturne-$version/etc/code_saturne.cfg" and adapt the file to
your needs.

Each user of Code_Saturne may set her/his PATH or define an alias accordingly
with the Code_Saturne installation before using the code.
The easiest way is to add the following
line in the user's ".profile" or ".alias" file(depending on the shell).

alias code_saturne="$prefix/code_saturne-$version/bin/code_saturne"

For more information please refer to the Code_Saturne documentation, available
through the "code_saturne info -g refcard" and "code_saturnes info -g user"
commands.

Code_Saturne support: saturne-support@edf.fr
