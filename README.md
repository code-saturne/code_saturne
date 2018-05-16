General Information
===================

This directory contains the [_Code_Saturne_](https://code-saturne.org) CFD tool,
EDF's general purpose Computational Fluid Dynamics (CFD) software.

The basic capabilities of _Code_Saturne_ enable the handling of either
incompressible or expandable flows with or without heat transfer and
turbulence. Dedicated modules are available for specific physics such as
radiative heat transfer, combustion (gas, coal, heavy fuel oil, ...),
magneto-hydrodynamics, compressible flows, two-phase flows
(Euler-Lagrange approach with two-way coupling), or atmospheric flows.

Presentation
============

_Code_Saturne_ is portable on all Linux flavors and UNIX platforms tested so far.
It runs in parallel with MPI on distributed memory machines (Intel, Cray X series,
IBM Blue Gene, IBM Power, ...).
Developed since 1997 at EDF R&D, it is based on a co-located Finite Volume
approach that accepts meshes with any type of cell (tetrahedral, hexahedral,
prismatic, pyramidal, polyhedral...) and any type of grid structure
(unstructured, block structured, hybrid, conforming or with hanging nodes, ...).

Meshes may be imported using the CGNS, MED, GMSH, I-Deas, GAMBIT, or Simail
formats, and Post-processing output is available in EnSight, CGNS
and MED formats. In-situ postprocessing is available using the
[ParaView Catalyst](https://www.paraview.org/in-situ) and
[Melissa](https://melissa-sa.github.io) libraries.

_Code_Saturne_ can be coupled in parallel to EDF's thermal software
[SYRTHES](https://www.edf.fr/en/the-edf-group/world-s-largest-power-company/activities/research-and-development/scientific-communities/simulation-softwares?logiciel=10818)
(conjugate heat transfer). It can also produce output usable by EDF's structural
analysis software [code_aster](https://code-aster.org), in particular in the
[SALOME platform](https://www.salome-platform.org/). SYRTHES and
code_aster are developed by EDF and distributed under a GNU GPL licence.

Copying
=======

_Code_Saturne_ is distributed under the GNU General Public Licence, v2.
See the COPYING file for details.

Installation
============

Manual installation
-------------------

Detailed installation instructions are also available as a pdf file,
available on the _Code_Saturne_ web site, or as part of this package.

For generic instructions relative to the GNU autotools-based
installation, see also the file 'INSTALL'.

Semi-automatic installation
---------------------------

The Install directory contains a python script for automatic
installation of the _Code_Saturne_ elements and associated routines.
In most cases, it will be enough. In case of problems, switch to
section II for element by element install.
These scripts are given in the hope that they will be useful, but
WITHOUT ANY WARRANTY.

The script can download every package needed by the code to run
properly. If this behaviour is not wanted, set the "download" variable
to "no" in the setup script.

Lastly, the possibility is given to compile _Code_Saturne_ with debugging symbols
("debug" variable), to disable the Graphical User Interface ("disable_gui"
variable), and to specify the language (between English and French).

* install_saturne.py:
  This python script will install the different elements of _Code_Saturne_ and
  associated libraries. Due to dependencies between the different modules, the
  order of install should be the following:

  - libxml2 (it is advised to use the distribution's own package instead)
  - HDF5
  - CGNS
  - MED
  - PT-Scotch
  - ParMETIS

  The following packages are not handled by the installer so must be installed
  first:

  - Zlib (optional)
  - BLAS (optional)
  - PyQT (optional, required for the GUI)
  - C, C++, and Fortran compilers
  - Python
  - MPI (optional)

  The install script uses a "setup" file to determine which library to
  install or to use. In not already present, This file is generated the
  first time the script is run in a given directory.
  For each element, there are four options:

  - do not use the element (for optional libraries like CGNS)
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

   Before using the "install_saturne.py" script, the C, Fortran, and optional
   C++ compilers to be used can be specified next to the CompC and CompF keywords.
   The Python interpreter may also be specified.
   MPI compiler wrappers (for the C and C++ languages) can also be specified.

   If the "use_arch" variable is set to "yes", then the "arch" keyword refers
   to the architecture of the machine. Leaving it blank will make it
   automatically detected with the "uname" command."arch" should be specified
   if you want different implementations on the same architecture
   (for instance Linux_OMPI and Linux_MPICH).

Post-install setup
------------------

For some systems (such as when using a batch system or coupling with SYRTHES,
a post-install step may be required). In this case, copy
`"$prefix/code_saturne-$version/etc/code_saturne.cfg.template"` to
`"$prefix/code_saturne-$version/etc/code_saturne.cfg"` and adapt the file to
your needs.

Each user of _Code_Saturne_ may set her/his PATH or define an alias accordingly
with the _Code_Saturne_ installation before using the code.
The easiest way is to add the following
line in the user's ".profile" or ".alias" file (depending on the shell).

`alias code_saturne="$prefix/code_saturne-$version/bin/code_saturne"`

For more information please refer to the _Code_Saturne_ documentation, available
through the "code_saturne info -g refcard" and "code_saturnes info -g user"
commands.

Code_Saturne support: saturne-support@edf.fr
