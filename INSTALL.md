<!--
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
-->

code_saturne installation documentation
=======================================

The code_saturne tool may be installed either directly through its GNU Autotools
based scripts (the traditional `configure`, `make`, `make install` sequence),
or using a 2-pass semi-automatic install script (`install_saturne.py`),
which generates an initial setup file when run a first time, and builds and
installs code_saturne and some optional libraries based on the edited setup
when run a second time.

The semi-automatic and regular installation modes may be combined: after an
initial run using the `install_saturne.py` script, the lower-level autotools
may be used to fine-tune the installation.

Preparing for build
===================

Source tree obtained through a source code repository
-----------------------------------------------------

For users obtaining the code was obtained through a Git repository, an additional
step is required:
```
$ cd code_saturne
$ ./sbin/bootstrap
$ cd ..
```
In this case, additional tools need to be available:
• GNU Autotools: Autoconf, Automake, Libtool (2.2 or 2.4).
• PdfLaTeX
• Doxygen (1.8.7 or more recent). The path to Doxygen can be specified during
the configure phase with `configure DOXYGEN=PATH TO DOXYGEN`.

These tools are not necessary for builds from tarballs; they are called when
building the tarball (using `make dist`), so as to reduce the number of
prerequisites for regular users, while developers building the code from a
repository can be expected to need a more complete development environment.
Also, to build and install the documentation when building the code from a
repository instead of a tarball, the following stages are required:
```
$ make doc
$ make install-doc
```

Source tree obtained through a tarball
--------------------------------------

When downloading a `.tar.gz` file, all that should be necessary is to unpack it:
```
tar xvzf code_saturne.tar.gz
```

Source, build, and install directories
--------------------------------------

3 types of directories are considered when installing the code:
* The *source* directory is the directory in which the source tree is cloned
  or unpacked
* The *build* directory is the directory in which the configuration and compilation
  will take place
* The *install* directory is the directry where the usable code will be installed.
  If unspecified, it will usually be `/usr/local` or `usr`.

In the following sections, we will assume a dedicated directory structure is used,
for example: `/home/user/code_Saturne/src` for source trees,
`/home/user/code_Saturne/src/build` for build trees, and `/home/user/code_Saturne/arch`
for installation (so as to keep directories both separate and easy to locate
relative to each other). Any directory structure may be chosen, as long as it is
easily comprehensible.

It is strongly recommended to build the code in a separate directory from the source
(the semi-automatic installer enforces this).
This also allows multiple builds, for example, building both an optimized and a
debugging version. In this case, choose a consistent naming scheme, using an additional
level of sub-directories, for example:

```
$ mkdir saturne_build
$ cd saturne_build
$ mkdir prod
$ cd prod
```

Also, as installing to a default system directory requires administrator-level write
permissions, and is not well adapted to multiple builds, this should be reserved
for Linux package maintainers.
In all other cases, it is strongly recommended to choose an installation path in a
user or project directory.
Installing code_saturne *does not require root privileges* if done correctly.

Semi-automatic installation
===========================

The `install_saturne.py` Python script present in the root directory of the
code_saturne source tree allows for semi-automatic installation of
code_saturne elements and some associated optional libraries. In most cases, it
will be sufficient to obtain a usable build. In case of problems, switching to
a regular install is recommended. This script is given in the hope that they it will
useful, but is simply a wrapper to the actual build system.

`install_saturne.py` works in 2 or more passes.

* on the first run, the script attemts to detect some tools, and generates
  a first `setup` configuration file;
* The `setup` file should be checked and modified according to need;
* on subsequent runs, libraries are downloaded and built, and code_saturne
  is configured or built, based on the `setup` configuration
* The script can be edited run multiple times, so that steps finished in
  prior runs do not need to be executed again if they succeeded.

The script can download the main packages needed by the code to run properly.
If these packages are not needed, or already present, set the *download* variable
to *no* in the `setup` script. If your system already provides adequate versions of
these tools, it is usually simpler to use the system-provided ones.

Lastly, the possibility is given to compile code_saturne with debugging symbols
(*debug* variable), and to disable the Graphical User Interface (*disable_gui*
variable).

Due to dependencies between the different modules, the order of install of
the optional modules should be the following:

  - HDF5
  - CGNS
  - MED
  - PT-Scotch
  - ParMETIS

The following packages are not handled by the installer so must be installed
first:

  - C, C++, and Fortran compilers
  - Python
  - MPI (optional, strongly recommended)
  - PyQt (optional, required for the GUI)
  - Zlib (optional)

Using the `setup` file
----------------------

As already mentioned, this file is generated the first time the script is
run in a given directory.

The C, Fortran, and optional C++ compilers to be used can be specified next
to the CompC and CompF keywords. The Python interpreter and MPI C and
C++ compiler wrappers can also be specified.

If the *use_arch* variable is set to *yes*, then the *arch* keyword refers
to the architecture of the machine. If left blank, the name will
be based on the `uname` command. When multiple builds are needed,
(such as simply using a production and debug build), *arch* should be
 specified for at least one of the builds.

The last section of the file relates to external libraries which can be
downloaded, installed, and used by code_saturne.

For each element, there are four options:

  - do not use the element (for optional libraries like CGNS)
     In this case, specify *no* in the *Use* and *Install* columns. The other
     elements will be installed in accordance. The *Path* column is not used.

  - automatically detect some element
     In this case, specify *auto* in the *Use* column. The other elements
     will be installed in accordance. The *Path* and *Install* column are
     not used.

  - use a pre-installed library in a non standard path
     In this case, specify *yes* in the *Use* column and *no* in the *Install*
     column. The *Path* column should contain the location of the library
     (up to the name of the library itself).

  - install and use a library
     In this case, specify *yes* in the *Use* and *Install* columns. The
     script will download the library and install it default install directory.
     If download has been set to *no*, package archive are looked for at the
     same location than the installation script (the right number and archive
     name are needed, accordingly to what is prescribed in the script).
     After each element has been installed, the `setup` file is modified, the
     column *Install* of the concerned element is set to *no* and the *Path*
     column is filled so that the element is not installed a second time if
     the script is relaunched (if there was a problem with a later element).

Note that when `install_saturne.py` is run, it generates a build directory, in
which the build may be modified (re-running configure, possibly adapting
the command logged at the beginning of the config.status file) and resumed.

So after an initial install is done using this script, it is possible to
switch to the regular installation to modify options and rebuild the code
starting from a first installation scheme.

Recommended post-install steps
==============================

Each user of _code_saturne_ may set her/his PATH or define an alias accordingly
with the _code_saturne_ installation, to avoid needing to type the full
install path name when using the code.
The easiest solution is to define a alias in the user's
`$HOME/.bashrc`, `$HOME/.alias`, or equivalent file (depending on the shell).

`alias code_saturne="<install-prefix>/bin/code_saturne"`

Note that when multiple versions of the code are installed side by side, using
a different alias for each will allow using them simultaneously, with no risk
of confusion. Aliases will be usable in a user shell but not in most scripts.
With no alias, using the absolute path is always an option.

Using the bash interpreter, automatic code completion help may also be set up,
using:

'source <install_prefix>/etc/bash_completion.d/code_saturne`

This may also be defined in a file such as `$HOME/.bashrc` or `$HOME/.bash_profile`
so as to be usable across sessions.

Global configuration file
-------------------------

For some systems (such as when using a batch system or coupling with SYRTHES),
additional post-install step may be required. In this case, copy
Copy or rename the `<install-prefix>/etc/code_saturne.cfg.template` to \\
`<install-prefix>/etc/code_saturne.cfg`,
and uncomment and define the applicable sections to adapt the file to
your needs.

This is useful for example for:

* Define which `batch` system is used. The name of the batch system should match
  one of the templates in `<install-prefix>/share/code_saturne/batch`,
  but an absolute path (with a file ending in the given batch system name) can
  be used to define a batch template taylored to a given system.

* Define `compute_versions` using the relative paths of alternate builds, so
  as to be able to use them from a "main" build. All specified builds are then
  available from the GUI and run script, which is especiall y useful to switch
  from a production to a debug build. In this case the secondary builds do not
  need to contain the full fromt-end (GUI, documentation, ...).

* Assign a name to the current build or system, so that the run configurations
  defined for each case (number of processors, batch template, ...) are saved
  specfically for that system. When using different systems on a same network,
  this allows keeping track of configurations specific to each system
  while moving or copying cases from one system to another.

* Defining the path to a SYRTHES installation to enable conjugate heat transfer.

* All default MPI execution commands and options may be overriden using the
  `mpi` section. Note that only the options defined in this section
  are overridden; defaults continue to apply for all others.

For more information on run configurations, please refer to the Doxygen
documentation.

Regular installation
====================

The installation scripts of code_saturne are based on the GNU Autotools,
(Autoconf, Automake, and Libtool), so it should be familiar for many
administrators. A few remarks are given here:

* As with most software with modern build systems, it is recommended
  to build the code in a separate directory from the sources. This
  allows multiple builds (for example production and debug), and is
  considered good practice. Building directly in the source tree is
  not regularly tested, and is not guaranteed to work, in addition
  to "polluting" the source directory with build files.
* By default, optional libraries which may be used by code_saturne are
  enabled automatically if detected in default search paths
  (i.e. `/usr/` and `/usr/local`. To find libraries
  associated with a package installed in an alternate path,
  a `--with-<package>=...` option to the `configure` script
  must given. To disable the use of a library which would be
  detected automatically, a matching `--without-<package>` option
  must be passed to `configure` instead.
* Most third-party libraries usable by code_saturne are considered optional,
  and are simply not used if not detected, but the libraries needed by
  the GUI are considered mandatory, unless the `--disable-gui`
  or `--disable-frontend` option is explicitly used.

When the prerequisites are available, and a build directory
created, building and installing code_saturne may be as simple as running:

```
$ ../../code_saturne/configure
$ make
$ make install
```

The following sections give more details on code_saturne's recommended
third-party libraries, configuration recommendations, troubleshooting,
and post-installation options.
We always assume that we are in build directory separate from sources.
In different examples, we assume that third-party libraries used by
code_saturne are either available as part of the base system (i.e. as
packages in a Linux distribution), as Environment Modules, or are
installed under a separate path.

Configuration
-------------

code_saturne uses a build system based on the GNU Autotools, which includes
its own documentation.

To obtain the full list of available configuration options,
run: `configure --help`.

Note that for all options starting with `--enable-`,
there is a matching option with `--disable-`. Similarly,
for every `--with-`, `--without-` is also possible.

Select configuration options, then run `configure`. For example,
if the code's source tree is in /home/user/code_Saturne/src/code_saturne:

```
$ /home/user/code_Saturne/src/code_saturne/configure \
--prefix=/home/user/code_saturne/<version>/arch/prod \
--with-med=/home/user/opt/med-4.1 \
CC=/home/user/opt/mpich-3.3/bin/mpicc FC=gfortran
```

Most available prerequisites are auto-detected, so to install the
code to the default `/usr/local` sub-directory, a command such as:

`$ ../../code_saturne-\verscs/configure`

should be sufficient.

### Defining compiler flags and associated environment variables

As usual when using an Autoconf-based `configure` script, some environment
variables may be used. `configure --help` will provide the list of
recognized variables. `CC`, `CXX`, and `FC` allow selecting the C; C++,
and Fortran compiler respectively (possibly using an MPI compiler wrapper).

Compiler options are usually defined automatically, based on
detection of the compiler (and depending on whether `--enable-debug`
was used). This is handled in a `config/cs_auto_flags.sh`
and `libple/config/ple_auto_flags.sh` scripts.
These files are sourced when running `configure`, so any modification to it
will be effective as soon as `configure` is run. When installing on an exotic
machine, or with a new compiler, adapting this file is useful (and providing
feedback to the code_saturne development team will enable support of a broader
range of compilers and systems in the future.

The usual `CPPFLAGS`, `CFLAGS`, `FCCFLAGS`, `LDFLAGS`, and `LIBS` environment
variables may also be used, and flags provided by the user are appended to
those automatic defined automatically. To completely disable automatic setting
of flags, the `--disable-auto-flags` option may be used.

Compile and install
-------------------

Once the code is configured, it may be compiled and installed;
for example, to compile the code (using 4 parallel threads),
then install it:

`$ make -j 4 && make install`

To compile the documentation, add:

`$ make doc && make install-doc`

To clean the build directory, keeping the configuration,
use `make clean`;
To uninstall an installed build, use `make uninstall`.
To clear all configuration info, use `make distclean`
(`make uninstall` will not work after this).

Troubleshooting
---------------

If `configure` fails and reports an error, the message should
be sufficiently clear in most case to understand the cause of the
error and fix it. Do not forget that for libraries installed using
packages, the development versions of those packages are also
necessary, so if configure fails to detect a package which you
believe is installed, check the matching development package.

Also, whether it succeeds or fails, `configure` generates
a file named `config.log`, which contains details on tests
run by the script, and is very useful to troubleshoot
configuration problems. When `configure` fails due to a given
third-party library, details on tests relative to that library
are found in the `config.log` file. The interesting information
is usually in the middle of the file, so you will need to search
for strings related to the library to find the test that failed
and detailed reasons for its failure.

Installing to a system directory
--------------------------------

When installing to a system directory in the default library search
path, such as `/usr` or `/usr/local`, some Linux systems may require
running `ldconfig` as root or sudoer for the code to work correctly.

Compilers and interpreters
==========================

For a minimal build of code_saturne on a Linux or Posix system, the requirements are:
* A C compiler, conforming at least to the C11 standard.
* A Fortran compiler, conforming at least to the Fortran 95 standard
  and supporting the ISO_C_BINDING Fortran 2003 module.
* A Python interpreter, with Python version 3.4 or above.

For parallel runs, an MPI library is also necessary (MPI-2 or MPI-3 conforming).
To build and use the GUI, PyQt 4 or 5 (which in turn requires Qt 4 or 5 and SIP)
are required. Other libraries may be used for additional mesh format options,
as well as to improve performance. A list of those libraries
and their role is given in a dedicated [section](@ref cs_install_list_ext_lib)).

For some external libraries, such as MED, MEDCoupling, and ParaView Catalyst,
a C++11 compliant C++ compiler is also required.

The SALOME platform V9 and above may require Python 3.6 at least,
and a matching version should be used when building with SALOME support.

In practice, the code is known to build and function properly at least with the
GNU compilers 4.4 and above (up to 9.x at this date), Intel compilers 11 and
above (up to 2020 at this date), and Clang (tested with 3.7 or above).

Note also that while code_saturne makes heavy use of Python, this is for scripts and
for the GUI only; The solver only uses compiled code, so we could for example use
a 32-bit version of Python with 64-bit code_saturne libraries and executables.
Also, the version of Python used by ParaView/Catalyst may be independent
from the one used for building code_saturne.

MPI compiler wrappers
---------------------

MPI environments generally provide compiler wrappers, usually with names such as
`mpicc` for C, `mpicxx` for C++, and `mpif90` for Fortran 90. Wrappers conforming
to the MPI standard recommendations should provide a `-show` option, to show
which flags are added to the compiler so as to enable MPI. Using wrappers is
fine as long as several third-party tools do not provide their own wrappers,
in which case either a priority must be established. For example, using HDF5's
`h5pcc` compiler wrapper includes the options used by `mpicc` when building HDF5
with parallel IO, in addition to HDF5's own flags, so it could be used instead of
`mpicc`. On the contrary, when using a serial build of HDF5 for a parallel
build of code_saturne, the `h5cc` and `mpicc` wrappers
contain different flags, so they are in conflict.

Also, some MPI compiler wrappers may include optimization options used to build
MPI, which may be different from those we wish to use that were passed.

To avoid issues with MPI wrappers, it is possible to select an MPI library using
the `--with-mpi` option to `configure`. For finer control,
`--with-mpi-include` and `--with-mpi-lib` may be defined separately.

Still, this may not work in all cases, as a fixed list of libraries
is tested for, so using MPI compiler wrappers remains the simplest and safest
solution. Simply use a `CC=mpicc` or similar option instead
of `--with-mpi`.

There is no need to use an `FC=mpif90` or equivalent option:
in code_saturne, MPI is never called directly from Fortran code,
so Fortran MPI bindings are not necessary nor recommended.

Interaction with system environment modules
===========================================

On some systems, especially compute clusters, _Environment Modules_
allow the administrators to provide multiple versions of many scientific
libraries, as well us compilers or MPI libraries, using the `module` command.
More details on different environment module systems may be found at
[http://modules.sourceforge.net](http://modules.sourceforge.net)
or [https://github.com/TACC/Lmod](https://github.com/TACC/Lmod).

The code_saturne `configure` script checks for modules loaded with the
`module` command, and records the list of loaded modules. Whenever running
that build of code_saturne, the modules detected at installation time will
be used, rather than those defined by default in the user's environment.
This allows using versions of code_saturne built with different modules safely
and easily, even if the user may be experimenting with other modules for
various purposes. This is especially useful when different builds use different
compiler or MPI versions.

Given this, it is recommended that when configuring and installing code_saturne,
only the modules necessary for that build of for profiling or debugging be
loaded. Note that as code_saturne uses the module environment detected and
runtime instead of the user's current module settings, debuggers requiring
a specific module may not work directly under a standard run script if they were
not loaded when installing the code.

The detection of environment modules may be disabled using the
`--without-modules` option, or the use of a specified (colon-separated) list
of modules may be forced using the `--with-modules=` option.

Pre-loading an environment
--------------------------

If installing and running code_saturne requires sourcing a given environment
or loading environment modules, the `--with-shell-env` option allows defining
the path for a file to source, or if no path is given, loading default modules.

By default, the main `code_saturne` command is a Python script. When sourcing an
environment, a launcher shell script is run first, loads the required environment,
then calls Python with the `code_saturne.py` script.

This is useful mainly when the Python version used is not installed in a
default system path, and must be loaded using an environment module, or
otherwise sourcing a specific `LD_LIBRARY_PATH` environment. Otherwise, the
main script may fail due to "library not found" issues.

Optional third-party libraries
==============================

In addition to an MPI library (for parallelism) and PyQt library (for the GUI),
other libraries may be used for additional mesh format options,
as well as to improve performance. A list of those libraries
and their role is detailed in this section.

Installing third-party libraries for code_saturne
-------------------------------------------------

Third-Party libraries usable with code_saturne may be installed in several
ways:

* On many Linux systems, most of optional third-party libraries
  are available through the distribution's package manager.

  Note that distributions usually split libraries or tools into runtime
  and development packages, and that although some packages are installed
  by default on many systems, this is generally not the case for the
  associated development headers. Development packages usually have the
  same name as the matching runtime package, with a `-dev` postfix added.
  Names might also differ slightly. For example, on a Debian system, the main
  package  for Open MPI is `openmpi-bin`, but `libopenmpi-dev` must also be
  installed for the code_saturne build to be able to use the former.

* If not otherwise available, third-party software may be compiled and installed
  by an administrator or a user. An administrator will choose where software
  may be installed, but for a user without administrator privileges or write
  access to `usr/local`, installation to a user account is often the only option.
  None of the third-party libraries usable by code_saturne require administrator
  privileges, so they may all be installed normally in a user account, provided the
  user has sufficient expertise to install them. This is usually not complicated
  (provided one reads the installation instructions, and is prepared to read error
  messages if something goes wrong), but even for an experienced user or
  administrator, compiling and installing 5 or 6 libraries as a prerequisite
  significantly increases the time and effort required to install code_saturne.

  Compiling and installing third-party software may still be necessary when no
  matching packages or Environment Modules are available, or when a more recent
  version or a build with different options is desired.

* When code_saturne is configured to use the SALOME platform, some libraries
  included in that platform may be used directly; this is described in the
  [Installation with the SALOME platform](#Installation-with-the-SALOME-platform)
  section.

List of third-party libraries usable by code_saturne
----------------------------------------------------

The list of third-party software usable with code_saturne is provided here:

* [PyQt](https://riverbankcomputing.com/software/pyqt/intro)
  version 4 or 5 is required by the code_saturne GUI. PyQt in turn requires
  Qt (4 or 5), Python, and SIP. Without this library, the GUI may not be
  built, although XML files generated with another install of code_saturne
  may be used.

  If desired, Using Qt for Python (PySide2) instead of PyQt should
  require a relatively small porting effort, as most of the preparatory work
  has been done. The development team should be contacted in this case.

* [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
  is necessary for MED, and may also be used by CGNS.

* [CGNS](http://cgns.github.io/) is necessary to read or write mesh and
  visualization files using the CGNS format, available as an export format with
  many third-party meshing tools. CGNS version 3.1 or above is required.

* [MED](https://old.salome-platform.org/user-section/about/med) is necessary to
  read or write mesh and visualization files using the MED format, mainly used by
  the SALOME platform.

* libCCMIO is necessary to read or write mesh and visualization files
  generated or readable by STAR-CCM+ using its native format.
  As recent versions of STAR-CCM+ also handle CGNS output, this is deprecated,
  and CGNS should be preferred.

* [Scotch or PT-Scotch](https://gitlab.inria.fr/scotch/scotch) may be used
  to optimize mesh partitioning. Depending on the mesh, parallel computations
  with meshes partitioned with these libraries may be from 10% to 50% faster
  than using the built-in space-filling curve based partitioning.

  As Scotch and PT-Scotch are part of the same package and use symbols with
  the same names, only one of the 2 may be used. If both are detected,
  PT-Scotch is used. Versions 6.0 and above are supported.

* [METIS](https://github.com/KarypisLab/METIS) or
  [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview)
  are alternative mesh partitioning libraries. These libraries have a
  separate source tree, but some of their functions have identical names,
  so only one of the 2 may be used. If both are available, ParMETIS will
  be used. Partitioning quality is similar to that obtained with Scotch or
  PT-Scotch.

  Though broadly available, the ParMETIS license is quite restrictive,
  so PT-Scotch may be preferred (code_saturne may be built with both METIS
  and Scotch libraries). METIS uses the Apache 2 licence since March 2013,
  but it seems that the ParMETIS licence has not been updated so far.
  METIS 5.0 or above and ParMETIS 4.0 or above are supported.

* [ParaView Catalyst](https://www.paraview.org/in-situ) or full ParaView
  may be used for co-visualization or in-situ visualization.
  This requires ParaView 4.2 or above, though 5.6 or above is recommended.
  Note that ParaView must be built with MPI support for Catalyst to
  be usable.

* [Melissa](https://melissa-sa.github.io) may be used for in-situ
  statistical analysis and post-processing of ensemble runs.

* `EOS` 1.2 or above may be used for thermodynamic properties of fluids.
   it is not currently free, so usually available only to users at EDF,
   CEA, or organisms participating in projects with those entities.

* [freesteam](http://freesteam.sourceforge.net) is a free software thermodynamic
  properties library, implementing the IAPWS-IF97 steam tables, from the
  [International Association for the Properties of Water and Steam (IAPWS)](http://www.iapws.org).
  Version 2.0 or above may be used.

* [CoolProp](http://www.coolprop.org) is a quite recent open source library,
  which provides pure and pseudo-pure fluid equations of state and transport properties
  for 122 components (as of version 6.3), mixture properties using high-accuracy Helmholtz
  energy formulations (or cubic EOS), correlations of properties
  of incompressible fluids and brines, fast IAPWS-IF97 (Industrial Formulation)
  for Water/Steam, and cubic equations of state (SRK, PR).
  Its validation is based at least in part on comparisons with NIST REFPROP.

* [SSH-aerosol](https://sshaerosol.wordpress.com/) model represents the
  physico chemical transformation undergone by aerosols in the troposphere.
  It is developed by Cerea (joint EDF/ENPC lab) and Ineris, and replaces
  the SIze REsolved Aerosol Model (SIREAM) available in code_saturne versions
  3.2 to 6.0.

* BLAS (Basic Linear Algebra Subroutines) may be used by the `cs_blas_test` unit
  test to compare the cost of operations such as vector sums and dot products with
  those provided by the code and compiler. If no third-party BLAS is provided,
  code_saturne reverts to its own implementation of BLAS routines, so no
  functionality is lost here. Optimized BLAS libraries such as Atlas or MKL
  may be very fast for BLAS3 (dense matrix/matrix operations),
  but the advantage is usually much less significant for BLAS 1 (vector/vector)
  operations, which are almost the only ones code_saturne has the opportunity of
  using. code_saturne uses its own dot product implementation (using a superblock
  algorithm, for better precision), and *y <- ax+y* operations, so external
  BLAS1 are not used for computation, but only for unit testing (so as
  to be able to compare performance of built-in BLAS with external BLAS).

  The only exception is the Intel MKL BLAS, which may also be used for sparse
  matrix-vector products, and is linked with the solver when available.
  It is then used at runtime on a case per case basis based on initial
  performance timings.

  Note that in some cases, threaded BLAS routines might oversubscribe
  processor cores in some MPI calculations, depending on the way both
  code_saturne and the BLAS were configured and interact, and this can actually
  lead to lower performance.
  Use of BLAS libraries is thus useful as a unit benchmarking feature,
  but has no influence on full calculations.

* [PETSc](https://www.mcs.anl.gov/petsc/) (Portable, Extensible Toolkit
  for Scientific Computation), consists  of a variety of libraries,
  which may be used by code_saturne for the resolution of linear equation systems.
  In addition to providing many solver options, it may be used as a bridge
  to other major solver libraries.
  Currently, only PETSc using install structures based using the
  `--prefix` configure option are correctly detected by the
  code_saturne build scripts.

* [HYPRE](http://www.llnl.gov/casc/hypre/) (high performance preconditioners)
  is a library of high performance preconditioners and solvers featuring
  multigrid methods for the solution of large, sparse linear systems of equations
  on massively parallel computers.

* [AmgX](https://developer.nvidia.com/amgx/) is a high performance
  multigrid and preconditioned iterative method library for NVIDIA GPUs.
  It includes a flexible solver composition system that allows
  a user to easily construct complex nested solvers and preconditioners.

* The [SYRTHES](https://www.edf.fr/en/the-edf-group/world-s-largest-power-company/activities/research-and-development/scientific-communities/simulation-softwares?logiciel=10818)
  code may be used for conjugate heat transfer.

For developers, the GNU Autotools (Autoconf, Automake, Libtool)
will be necessary. To build the documentation, pdfLaTeX,
and Doxygen are needed, and dot (from Graphviz) recommended.

PLE (Parallel Location and Exchange) Library
--------------------------------------------

By default, the PLE library is built as a subset of code_saturne, but a
standalone version may be configured and built, using the `libple/configure`
script from the code_saturne source tree, instead of the top-level `configure`
script. In this case, code_saturne may then be configured to use an existing
install of PLE using the `--with-ple` configure option.

Special notes on some third-party tools and libraries
-----------------------------------------------------

### Python and PyQt

The GUI is written in PyQt (Python bindings for Qt), so Qt (version 4 or 5)
and the matching Python bindings must be available. On most modern
Linux distributions, this is available through the package manager,
which is by far the preferred solution.

On systems on which both PyQt4 and Pyqt5 are available, PyQt5 will be selected
by default, but the selection may be forced by defining
`QT_SELECT=4` or `QT_SELECT=5`.

When running on a system which does
not provide these libraries, there are several alternatives:

* build code_saturne without the GUI. XML files produced with the GUI are
  still usable, so if an install of code_saturne with the GUI is available
  on an other machine, the XML files may be copied on the current machine.
  This is not an optimal solution, but in the case where users have a mix of
  desktop or virtual machines with full Linux desktop including PyQt installed,
  and a compute cluster with a more limited or older system, this may avoid
  requiring a build of Qt and PyQt on the cluster if users find this too daunting.

* Install a local Python interpreter, and add Qt5 bindings to this interpreter.

  [Python](http://www.python.org) and [Qt](https://www.qt.io) must be downloaded
  and installed first, in any order. The installation instructions of both of these
  tools are quite clear, and though the installation of these large packages
  (especially Qt) may be a lengthy process in terms of compilation time, but is
  well automated and usually devoid of nasty surprises.

  Once Python is installed, the [SIP](http://riverbankcomputing.co.uk/software/sip/intro)
  bindings generator must also be installed. This is a small package, and
  configuring it simply requires running `python configure.py` in its source
  tree, using the Python interpreter just installed.

  Finally, the [PyQt](http://riverbankcomputing.co.uk/software/pyqt/intro)
  bindings, in a manner similar to SIP.

  When this is finished, the local Python interpreter contains the PyQt bindings, and
  may be used by code_saturne's `configure` script by passing
  `PYTHON=<path_to_python_executable`.

* add Python Qt bindings as a Python extension module for an existing
  Python installation. This is a more elegant solution than the previous
  one, and avoids requiring rebuilding Python, but if the user does not
  have administrator privileges, the extensions will be placed in a
  directory that is not on the default Python extension search path, and
  that must be added to the `PYTHONPATH` environment variable.
  This works fine, but for all users using this build of code_saturne, the
  `PYTHONPATH` environment variable will need to be set.

  The process is similar to the previous one, but SIP and PyQt installation
  requires a few additional configuration options in this case. See the SIP
  and PyQt reference guides for detailed instructions, especially the
  *Building a Private Copy of the SIP Module` section of the SIP guide*.

### Scotch and PT-Scotch

Note that both Scotch and PT-Scotch may be built from the same source
tree, and installed together with no name conflicts.

For better performance, PT-Scotch may be built to use threads with concurrent
MPI calls. This requires initializing MPI with `MPI_Init_thread`
with `MPI_THREAD_MULTIPLE` (instead of the more restrictive
`MPI_THREAD_SERIALIZED`, `MPI_THREAD_FUNNELED`, or
`MPI_THREAD_SINGLE`, or simply using `MPI_Init`).
As code_saturne does not currently require thread models in which different
threads may call MPI functions simultaneously, and the use of `MPI_THREAD_MULTIPLE`
may carry a performance penalty, we prefer to sacrifice some of PT-Scotch's *
performance by requiring that it be compiled without the `-DSCOTCH_PTHREAD` flag.
This is not detected at compilation time, but with recent MPI libraries,
PT-Scotch will complain at run time if it notices that the MPI thread safety
level in insufficient.

Detailed build instructions, including troubleshooting instructions, are given
in the source tree's `INSTALL.txt` file. In case of trouble, note especially the
explanation relative to the `dummysizes` executable, which is run to determine the
sizes of structures. On machines with different front-end and compute node
architectures, it may be necessary to start the build process, let it fail,
run this executable manually using `mpiexec`}, then pursue the build process.

Note that PT-Scotch may be installed by code_saturne's semi-automatic
installer, which chooses settings that should work on most machines.

### MED

MED can be built using either CMake or the GNU Autotools.
The Autotools installation of MED is simple on most machines,
but a few remarks may be useful for specific cases.

Note that up to MED 3.3.1, HDF5 1.8 was required, while MED 4.x
uses  HDF5 1.10. It does not accept HDF5 1.12 yet.

MED has a C API, is written in a mix of C and C++ code, and provides both
a C (`libmedC`) and an Fortran API (`libmed`) by default (i.e. unless
the `--disable-fortran` configure option is used. code_saturne only requires
the C API.

MED requires a C++ runtime library, which is usually transparent when shared
libraries are used. When built with static libraries only, this is not sufficient,
so when testing for a MED library, the code_saturne `configure` script also tries
linking with a C++ compiler if linking with a C compiler fails. This must be the
same compiler that was used for MED, to ensure the runtime matches. The choice of
this C++ compiler may be defined passing the standard `CXX` variable to `configure`.

Also, when building MED in a cross-compiling situation, `--med-int=int` or
`--med-int=int64_t` (depending on whether 32 or 64 bit ids should be used) should
be passed to its `configure` script to avoid a run-time test.

### libCCMIO

Different versions of this library may use different build systems, and use
different names for library directories, so using both the `--with-ccm=`
or `--with-ccm-include=` and `--with-ccm-lib=` options to `configure` is
usually necessary. Also, the include directory should be the toplevel library,
as header files are searched under a `libccmio` subdirectory.
This is made necessary by libCCMIO version 2.6.1, in which this is hard-coded in
headers including other headers. In more recent versions such as 2.06.023, this
is not the case anymore, and an `include` subdirectory is present, but
it does not contain the `libccmioversion.h` file, which is
found only under the `libccmio` subdirectory, but is required
by code_saturne to handle differences between versions, so that source
directory is preferred to the installation `include`.

A libCCMIO distribution usually contains precompiled binaries, but recompiling
the library is recommended. Note that at least for version 2.06.023, the
build will fail building dump utilities, due to the `-l adf` option
being placed too early in the link command. To avoid this, add `LDLIBS=-ladf`
to the makefile command, for example:

`make -f Makefile.linux SHARED=1 DEBUG=1 LDLIBS=-ladf`

`SHARED=1` and `DEBUG=1` may be used to force shared library or debug
builds respectively, but some issues have been observed on some systems
using non-debug builds, so adding that option is recommended.

Finally, if using libCCMIO 2.6.1, remove the `libcgns*` files from the libCCMIO
libraries directory if also building code_saturne with CGNS support, as those
libraries are not required for CCMIO, and are are an obsolete version of
CGNS, which may interfere with the version used by code_saturne.

Note that libCCMIO uses a modified version of CGNS's ADF library, which may not
be compatible with that of CGNS. When building with shared libraries, the reader
for libCCMIO uses a plugin architecture to load the library dynamically.
For a static build with both libCCMIO and CGNS support, reading ADF-based
CGNS files may fail. To work around this issue, CGNS files may be converted
to HDF5 using the `adf2hdf` utility (from the CGNS tools). By default,
CGNS post-processing output files use HDF5, so this issue is rare on output.

### freesteam

This library's build instructions mention bindings with ascend, but those are not
necessary in the context of code_saturne, so building without them is simplest.
Its build system is based on scons, and builds on relatively recent systems should
be straightforward.

### CoolProp

This library's build system is based on CMake, and building it is straightforward,
though some versions seem to have build issues (the 5.1.0 release is missing a file,
while 5.1.1 release builds fine). CoolProp uses submodules which are downloaded
using `git clone https://github.com/CoolProp/CoolProp.git --recursive`,
but may be missing when downloading a zip file.

Its user documentation is good, but its installation documentation not so much,
so recommendations are provided here

To download and prepare CoolProp for build, using an out-of-tree build
(so as to avoid polluting the source tree with cache files), the
following commands are recommended:

```
$ git clone https://github.com/CoolProp/CoolProp.git --recursive
$ cd CoolProp
$ git checkout release
$ cd ..
$ mkdir CoolProp_build
$ cd CoolProp_build
```

Then configure, build, and install, run:

```
$ cmake \
-DCOOLPROP_INSTALL_PREFIX=${INSTALL_PATH} \
-DCOOLPROP_SHARED_LIBRARY=ON \
${COOLPROP_SRC_PATH}
```

Followed by:

```
$ make
$ make install
$ make clean
```

CoolProp's installer only installs one C wrapper header, not the
C++ headers required for lower-level access, so the following commands
must also be run:

```
$ cp -rp ${COOLPROP_SRC_PATH}/include ${INSTALL_PATH}
$ rm -f ${INSTALL_PATH}/CoolPropLib.h
```

Alternatively, to copy less files and avoid changing the structure provided
by CoolProp:

```
$ cp -r ${COOLPROP_SRC_PATH}/include ${INSTALL_PATH}
$ cp -r ${COOLPROP_SRC_PATH}/externals/fmtlib/fmt ${INSTALL_PATH}/include/
```

To install CoolProp's Python bindings (used by the GUI when available),
the straigthforward method is to go into the CoolProp source directory,
into the wrappers/Python subdirectory, then run:

```
$ export PYTHONPATH=${COOLPROP_INSTALL_PREFIX}/lib64/${python_version}/site-packages:$PYTHONPATH
$ python setup.py install --prefix=${COOLPROP_INSTALL_PREFIX}
```

Although this is not really an out-of-tree build, the Python
setup also cleans the directory.

### Paraview or Catalyst

By default, this library is built with a GUI, but it may also be be
built using OSMesa for offscreen rendering. The build documentation
on the ParaView website and Wiki details this. On a workstation,
a regular build of ParaViw with MPI support may be sufficient.

For a compute cluster or server in which code_saturne will run
outside a graphical (X11) environment, the recommended solution is to build
or use a standard ParaView build for interactive visualization, and to use
its Catalyst/Define Exports option to generate Python co-processing scripts.
A second build, using OSMesa (or EGL), may be used for in-situ visualization.
This is the Version code_saturne will be linked to. A recommended cmake
command for this build contains:

```
$ cmake \
-DCMAKE_INSTALL_PREFIX=${INSTALL_PATH}_osmesa \
-DPARAVIEW_USE_QT=OFF \
-DPARAVIEW_USE_MPI=ON \
-DPARAVIEW_USE_PYTHON=ON \
-DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
-DOSMESA_INCLUDE_DIR=${MESA_INSTALL_PREFIX}/include \
-DOSMESA_LIBRARY=${MESA_INSTALL_PREFIX}/lib/libOSMesa.so \
-DVTK_OPENGL_HAS_OSMESA=ON \
-DVTK_USE_X=OFF \
${PARAVIEW_SRC_PATH}
```

More info may also be found on the
[ParaView Wiki](http://www.paraview.org/Wiki/ParaView/ParaView_And_Mesa_3D).

Note that when ParaView uses libraries which are in non-standard locations,
it may be necessary to specify those locations in the CMake prefix path
for ParaView detection by code_saturne. Actually, the option
passed to `--with-paraview` when running code_saturne's `configure` step is
the `CMAKE_PREFIX_PATH`, so if multiple directories need to be included,
an enquoted and semicolon-separated path may be used, for example:

```
--with-catalyst="/home/user/opt/paraview-5.8;/home/user/opt/ospray2"
```

Also, if the detection of Catalyst fails due to incorrect detection
of the TBB library, the `TBB_INCLUDE_DIR` environment variable may
be set to pass the correct path to the configuration scripts.

On some systems, loading the Catalyst module as a plug-in (which is the
default) seems to interfere with the detection of required OpenGL2 features
or extensions required by ParaView 5.2 and above. In this case, Catalyst
support may be linked in the standard manner by using the
`--disable-catalyst-as-plugin` configuration option.
A less extreme option is to use the `--enable-dlopen-rtld-global`
option, which changes the system options with which libraries are loaded
(possibly impacting all plugins). This seems to be sufficient with
OSMesa 17.x versions. Using the `LD_PRELOAD` environment variable
at runtime to preload the OSMesa library also avoids the issue.

On at least one system (with RHEL 8.3 and the gcc 8.3.1 compiler),
a crash at finalization of Catalyst has been observed. If this is the case,
setting the `CS_PV_CP_DELETE_CRASH_WORKAROUND` environment variable to 1
should avoid calling the offending code.

### Coupling with SYRTHES

Coupling with SYRTHES requires defining the path to a SYRTHES installation
at the post-install stage (in the `code_saturne.cfg` file).

Both code_saturne and SYRTHES must use the same MPI library, and must use
the same major version of the PLE library from
code_saturne. If SYRTHES is installed after code_saturne, the simplest solution
is to configure its build to use the PLE library from the existing
code_saturne install.

Specific build types
====================

Debug builds
------------

It may be useful to install debug builds alongside production builds of
code_saturne, especially when user-defined functions are used and the risk of
crashes due to user programming error is high.
Running the code using a debug build is significantly slower, but more
information may be available in the case of a crash, helping understand
and fix the problem faster.

Here, having a consistent and practical naming scheme is useful.
For a side-by-side debug build for the initial example, we simply replace `prod` by
`dbg` in the `--prefix` option, and add:
`--enable-debug` to the configure command:

```
$ cd ..
$ mkdir dbg
$ cd dbg
$ ../../code_saturne/configure \
--prefix=/home/user/code_saturne/<version>/arch/dbg \
--with-med}=/home/user/opt/med-4.1 \
--enable-debug \
CC=/home/user/opt/mpich-3.3/bin/mpicc FC=gfortran
```

Shared or static builds
-----------------------

By default, on most architectures, code_saturne will be built with shared
libraries. Shared libraries may be disabled (in which case static libraries
are automatically enabled) by adding  `--disable-shared` to the options
passed to `configure`.
On some systems, the build may default to static libraries instead.

It is possible to build both shared and static libraries by not adding
`--disable-static` to the `configure` options, but the
executables will be linked with the shared version of the libraries,
so this is rarely useful (the build process is also slower in this case, as
each file is compiled twice).

In some cases, a shared build may fail due to some dependencies
on static-only libraries. In this case, `--disable-shared`
will be necessary. Disabling shared libraries is also necessary
to avoid issues with linking user functions on Mac OS-X systems.

In any case, be careful if you switch from one option to the other: as
linking will be done with shared libraries by default, a build
with static libraries only will not completely overwrite a build using
shared libraries, so uninstalling the previous build first
is recommended.

Relocatable builds
------------------

By default, a build of code_saturne is not movable, as not only
are library paths hard-coded using `rpath` type info,
but the code's scripts also contain absolute paths.

To ensure a build is movable, pass the `--enable-relocatable` option
to `configure`.

*Movable builds assume a standard directory hierarchy*, so when running
`configure`, the `--prefix` option may be used, but fine tuning
of installation directories using options such as `--bindir`,
`--libdir`, or `--docdir` must not be used (these options are useful
to install to strict directory hierarchies, such as when packaging the
code for a Linux distribution, in which case making the build relocatable
would be nonsense anyways, so this is not an issue.

**Remark**: for relocatable builds using ParaView/Catalyst, a `CATALYST_ROOT_DIR`
environment variable may be used to specify the Catalyst location
in case that was moved also.

### For package maintainers

In the special case of packaging the code, which may require both
 fine-grained control of the installation directories and the possibility to
support options such as `dpkg`'s `--instdir`, it is assumed the packager has
sufficient knowledge to update both `rpath` information and paths in scripts
in the executables and python package directories of a non-relocatable build,
and that the packaging mechanism includes the necessary tools and scripts to
enable this.

Remarks for large meshes
------------------------

If code_saturne is to be run on large meshes, several precautions regarding
its configuration and that of third-party software must be taken.

in addition to local connectivity arrays, code_saturne uses global element ids
for some operations, such as reading and writing meshes and restart files,
parallel interface element matching, and post-processing output.
For a hexahedral mesh with *N* cells, the number of faces is about *3N*
(6 faces per cell, shared by 2 cells each). With 4 cells per face, the
*face to vertices* array is of size of the order of *4 * 3N*, so global ids
used in that array's index will reach *2<sup>31</sup>* for a mesh in the range
of *2<sup>31</sup> / 12* (approximately 178 million cells).
In practice, we have encountered a limit with slightly smaller meshes.

Above 150 million hexahedral cells or so, it is thus imperative to configure
the build to use 64-bit global element ids. This is the default.
Local indexes use the default `int` size. To slightly decrease memory
consumption if meshes of this size are never expected (for example on a workstation
or a small cluster), the `--disable-long-gnum` option may be used.

Recent versions of some third-party libraries may also optionally use 64-bit ids,
independently of each other or of code_saturne. This is the case for the
Scotch, METIS, MED and CGNS libraries. In the case of graph-based partitioning,
only global cell ids are used, so 64-bit ids should not in theory be necessary
for meshes under 2 billion cells. In a similar vein, for post-processing output
using nodal connectivity, 64-bit global ids should only be an imperative
when the number of cells or vertices approaches 2 billion.
Practical limits may be lower, if some intermediate internal counts
reach these limits earlier.

Partitioning a 158 million hexahedral mesh using serial METIS 5 or Scotch
on a front-end node with 128 Gb memory is possible, but partitioning the
same mesh on cluster nodes with "only" 24 Gb each may not, so using parallel
partitioning PT-Scotch or ParMETIS should be preferred.

Installation with the SALOME platform
-------------------------------------

To simplify the detection of some libraries provided with a [SALOME platform]
(http://www.salome-platform.org) installation, the `--with-salome` configuration
option may be used, so as to specify the directory of the SALOME installation
(note that this should be the main installation directory, not the default
application directory, also generated by SALOME's installers).

Specifying a SALOME directory is only a first step, and does not automatically
force the code_saturne `configure` script to find some libraries
which may be available in the SALOME distribution, such as HDF5,
MED, or CGNS. To indicate that the versions from SALOME should be used,
without needing to provide the full paths, the following configuration
options may be used for HDF5, CGNS, and MED respectively, as well
as for Catalyst when available in a given Salome platform variant.

 ```
--with-hdf5=salome
--with-cgns=salome
--with-med=salome
--with-catalyst=salome
```

When those libraries are available as system packages and not
inside the SALOME directory, these options might not be available.

As CGNS and MED file formats are portable, MED or CGNS files produced
by either code_saturne or SALOME remain interoperable. At the least,
files produced with a given version of CGNS or MED should be readable
with the same or a newer version of that library.

Unless a specific `--with-medcoupling` option is given, a compatible
MEDCoupling library is also searched for in the SALOME distribution.

Also note that for SALOME builds containing their own Python interpreter
and library, using that same interpreter for code_saturne may avoid some
issues with Python-based tools, but may then require sourcing the SALOME
environment or at least its Python-related `LD_LIBRARY_PATH` for the
main code_saturne script to be usable. The `--with-shell-ext` configure
option of code_saturne is useful here.

SALOME platform extensions
--------------------------

The SALOME platform extensions are now provided in a separate repository,
with a public mirror at
[https://github.com/code-saturne/salome_cfd_extensions](https://github.com/code-saturne/salome_cfd_extensions). They can be installed separately once
the main code_saturne install is completed (though in fact installation order
is not important).

These extensions contain the CFDSTUDY salome module (available by running
`salome_cfd` after install), which provide workbench integration and
links with the [OpenTURNS](https://openturns.github.io/www/index.html)
and the [PERSALYS](https://persalys.fr) graphical interface for sensitivity
studies.

To simply use code_saturne from a salome shell session,
these extensions are not necessary.

Note finally that SALOME expects a specific directory tree when loading modules,
so salome_cfd extensions might fail when installing with a specified
(i.e. non-default) `--datarootdir` path in the code_saturne `configure` options.

Configuration command examples
==============================

For the following advanced examples, Let us define environment variables
respectively reflecting the code_saturne source path, installation path,
and a path where optional libraries are installed:

```
$ SRC_PATH=/home/projects/code_saturne/<x.y>/src/code_saturne-<x.y.z>
$ INSTALL_PATH=/home/projects/code_saturne/<x.y>
$ CS_OPT=/home/projects/opt
```

Here, we have a general `/home/projects/code_saturne` directory for all
code_saturne directories, with a `x.y` directory for each (major.minor) version,
including `src` and `arch` subdirectories for sources and installations
respectively. Only the installation directories are needed to use the code,
but keeping the sources nearby may be practical for reference.

For an install on which multiple versions and architectures of the code should
be available, configure commands with all bells and whistles (except SALOME
support) for a build on a cluster named gaia}, using the Intel compilers
(made available through environment modules) may look like this:

```
module purge
module load intel_compilers/2019.0.045
module load open_mpi/gcc/4.0.1
${SRC_PATH}/configure \
--prefix=${INSTALL_PATH}/arch/gaia_ompi \
--with-blas=/opt/mkl-2019.0.045/mkl \
--with-hdf5=${CS_OPT}/hdf5-1.10/arch/gaia \
--with-med=${CS_OPT}/med-4.1/arch/gaia \
--with-cgns=${CS_OPT}/cgns-3.4/arch/gaia \
--with-ccm=${CS_OPT}/libccmio-2.06.23/arch/gaia \
--with-scotch=${CS_OPT}/scotch-6.0/arch/gaia_ompi \
--with-metis=${CS_OPT}/parmetis-4.0/arch/gaia_ompi \
--with-eos/${CS_OPT}/eos-1.2.0/arch/gaia_ompi \
CC=mpicc FC=ifort CXX=mpicxx
```

In the example above, we have appended the `_ompi` postfix to the architecture
name for libraries using MPI, in case we intend to install 2 builds, with
different MPI libraries (such as Open MPI and MPICH-based Intel MPI).
Note that optional libraries using MPI must also use the same MPI library.
This is the case for PT-Scotch or ParMETIS, but also HDF5, CGNS, and MED if
they are built with MPI-IO support. Similarly, C++ and Fortran libraries,
and even C libraries built with recent optimizing C compilers, may require
runtime libraries associated to that compiler, so if versions using different
compilers are to be installed, it is recommended to use a naming scheme which
reflects this. In this example, HDF5, CGNS and MED were built without MPI-IO
support, as code_saturne does not yet exploit MPI-IO for these libraries.

To avoid copying platform-independent data (such as the documentation)
from different builds multiple times, we may use the same
`--datarootdir` option for each build so as to install that
data to the same location for each build. In this case, we recommend
`--datarootdir=${INSTALL_PATH}/share` to add a `share` directory alongside
`arch`.

Cross-compiling
---------------

On machines with different front-end and compute node architectures,
cross-compiling may be necessary.
To install and run code_saturne, 2 builds are then required:

* a *front-end* build, based on front-end node's architecture. This is
  the build whose `code_saturne` command, GUI, and documentation
  will be used, and with which meshes may be imported (i.e. whose
  Preprocessor will be used). This build is not intended for calculations,
  though it could be used for mesh quality criteria checks.
  This build will thus usually not need  MPI.
* a *compute* build, cross-compiled to run on the compute nodes.
  This build does not need to include the GUI, documentation, or
  even the preprocessor module (see user documentation).

A debug variant of the compute build is also recommended, as always.
Providing a debug variant of the front-end build is not generally useful.

A post-install step (see Post-install setup) will allow the scripts of the
front-end build to access the compute build in a transparent manner, so it
will appear to the users that they are simply working with that build.

Depending on their role, optional third-party libraries should be installed
either for the front-end, for the compute nodes, or both:

* BLAS will be useful only for the compute nodes, and are generally
  always available on large compute facilities.
* Python and PyQt will run on the front-end node only.
* HDF5, MED, CGNSlib, and libCCMIO may be used by the Preprocessor on
  the front-end node to import meshes, and by the main solver on the
  compute nodes to output visualization meshes and fields.
* Scotch or METIS may be used by a front-end node build of the solver,
  as serial partitioning of large meshes requires a lot of memory.
* PT-Scotch or ParMETIS may be used by the main solver on the compute nodes.

Compiling for Cray X series
---------------------------

For Cray X series, when using the GNU compilers, installation should
be similar to that on standard clusters. Using The Cray compilers,
options such as in the following example are recommended:

```
\$ ${SRC_PATH}/configure \
--prefix=${INSTALL_PATH}/arch/xc30 \
--with-hdf5=${CS_OPT}/hdf5-1.10/arch/xc30 \
--with-med}=${CS_OPT}/med-4.1/arch/xc30 \
--with-cgns}=${CS_OPT}/cgns-3.4/arch/xc30 \
--with-scotch}=${CS_OPT}/scotch-6.0/arch/xc30 \
--disable-sockets \
--disable-shared \
--host=x86_64-unknown-linux-gnu \
CC=craycc CXX=crayc++ FC=crayftn
```

In case the automated environment modules handling causes issues,
adding the `--without-modules` option may be necessary.
In that case, caution must be exercised so that the user will load
the same modules as those used for installation. This is not an issue
if modules for code_saturne is also built, and the right dependencies
handled at that level.

Note that to build without OpenMP with the Cray compilers,
`CFLAGS=-h noomp` and `FCFLAGS=-h noomp` need to be added.

Caveats
=======

Moving an existing installation
-------------------------------

**Never move a non-relocatable installation** of code_saturne.

Using `LD_LIBRARY_PATH` or `LD_PRELOAD` may allow the executable to run
despite `rpath` info not being up-to-date, but in environments where
different library versions are available, there is a strong risk of not using
the correct library. In addition, the scripts will not work unless paths in
the installed scripts are updated.

To build a relocatable installation, see the Relocatable builds section.

If you are packaging the code and need both fine-grained control of
the installation directories, and the possibility to support
options such as `dpkg`'s `--instdir`, it is assumed
you have sufficient knowledge to update both `rpath` information
and paths in scripts in the executables and python package directories,
and that the packaging mechanism includes the necessary tools and
scripts to enable this.

In any other case, you should not even think about moving a non-relocatable
(i.e. default) instead of using a relocatable one.

If you need to test an installation in a test directory before
installing it in a production directory, use the
`make install DESTDIR=<test_prefix>` provided by the Autotools mechanism rather
than configuring an install for a test directory and then moving it to a
production directory. Another less elegant but safe solution is to configure the
build for installation to a test directory, and once it is tested,
re-configure the build for installation to the final production
directory, and rebuild and install.

Dynamic linking and path issues on some systems
-----------------------------------------------

On Linux and Unix-like systems, there are several ways for a library or executable
to find dynamic libraries, listed here in decreasing priority:

* the `LD_PRELOAD` environment variable explicitly lists
  libraries to be loaded with maximum priority, before the libraries
  otherwise specified (useful mainly for instrumentation
  and debugging, and should be avoided otherwise);
* the `RPATH` binary header of the dependent library or executable;
  (if both are present, the library has priority);
* the `LD_LIBRARY_PATH` environment variable;
* the `RUNPATH` binary header of executable;
* `/etc/ld.so.cache`};
* base library directories (`/lib` and /`/usr/lib`);

Note that changing the last two items usually require administrator privileges,
and we have encountered cases where installing to /`/usr/lib` was not
sufficient without updating `/etc/ld.so.cache`. We do not consider
`LD_PRELOAD` here, as it has other specific uses.

So basically, when using libraries in non-default paths, the remaining
options are between `RPATH` or `RUNPATH` binary headers, or the
`LD_LIBRARY_PATH` environment variable.

The major advantage of using binary headers is that the executable
can be run without needing to source a specific environment, which
is very useful, especially when running under MPI (where the propagation
of environment variables may depend on the MPI library and batch system's
configuration), or running under debuggers (where library paths would have
to be sourced first).

In addition, the `RPATH` binary header has priority over
`LD_LIBRARY_PATH`, allowing the installation to be "protected"
from settings in the user environment required by other tools.
this is why the code_saturne installation chooses this mode by default,
unless the `--enable-relocatable` option is passed to
`configure`.

Unfortunately, the ELF library spec indicates that the use of the
`DT_RPATH` entry (for `RPATH`) has been superseded by the
`DT_RUNPATH` (for `RUNPATH}). Most systems still use
`RPATH}, but some (such as SUSE and Gentoo) have defaulted to
`RUNPATH`, which provides no way of "protecting" an executable
or library from external settings.

Also, the `--enable-new-dtags` linker option allows replacing
`RPATH` with `RUNPATH}, so adding `-Wl,--enable-new-dtags`
to the `configure` options will do this.

The addition of `RUNPATH` to the ELF specifications may have
corrected the oversight of not being able to supersede an executable's
settings when needed (though `LD_PRELOAD` is usually sufficient
for debugging, but the bright minds who decided that it should
replace `RPATH` and not simply supplement it did not provide
a solution for the following scenario:

* code_saturne in installed, along with the MED and HDF libraries,
  on a system where `--enable-new-dtags` is enabled by default.
* Another code is installed, with its own (older) versions of MED and HDF
  libraries; this second code requires sourcing environment variables
  including `LD_LIBRARY_PATH` to work at all, so the user
  adds those libraries to his environment (or the admins add it to
  environment modules).
* code_saturne now uses the older libraries, and is not capable of reading
  files generated with more recent versions

The aforementioned scenario occurs with code_saturne and `rpath`, on some
machines, and could occur with code_saturne and some SALOME libraries, and there is
no way around it short of changing the installation logic of these other tools,
or using a cumbersome wrapper to launch code_saturne, which could still fail when
code_saturne needs to load a `rpath` or SALOME environment for coupled cases.
A wrapper would lead to its own problems, as for example Qt is needed
by the GUI but not the executable, so to avoid causing issues with
a debugger using its own version of Qt, separate sections would need
to be defined. None of those issues exist with `RPATH`.

To avoid most issues, the code_saturne scripts also update `LD_LIBRARY_PATH`
before calling executable modules, but you could be affected if running
them directly from the command line.
