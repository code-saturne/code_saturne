General Information
===================

This directory contains the salome_cfd extensions for the
code_saturne (https://code-saturne.org) CFD tool,
EDF's general purpose Computational Fluid Dynamics (CFD) software.

These extension provide better integration with the SALOME plarform
(https://www.salome-platform.org).

Presentation
============

The salome_cfd extensions provide increased integration of code_saturne
in the Salome platform, especially:

* the CFD_STUDY module for integration of the code_saturne GUI in the SALOME
  workbench; this also allows visualization of probe positions and boundary
  groups relative to the mesh;

* integration with OpenTURNS and the PERSALYS graphical interface for
  sensitivity studies.

Copying
=======

the salome_cfd extensions for code_saturne are distributed under the GNU
General Public Licence, v2. or higher. See the COPYING file for details.

Installation
============

Detailed installation instructions are also available as a pdf file,
available on the code_saturne web site, or as part of the
code_saturne package.

Installation may be done as part of the code_saturne installation or
as a post-install to extend an existing installation.

It is based on GNU autotools, so the classical
`configure && make && make install` paradigm may be used here.
It is strongly recommended to build ouside the source tree (i.e. run
`configure` from outside the source tree, in a dedicated build directory),
and in-tree builds are not supported.

To obtain available options, run `configure --help`

Code_Saturne support: saturne-support@edf.fr
