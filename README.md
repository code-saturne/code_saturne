General Information
===================

This directory contains the [code_caturne](https://code-saturne.org) CFD tool,
EDF's general purpose Computational Fluid Dynamics (CFD) software.

The basic capabilities of code_saturne enable the handling of either
incompressible or expandable flows with or without heat transfer and
turbulence. Dedicated modules are available for specific physics such as
radiative heat transfer, combustion (gas, coal, heavy fuel oil, ...),
magneto-hydrodynamics, compressible flows, two-phase flows
(Euler-Lagrange approach with two-way coupling), or atmospheric flows.

Presentation
============

code_saturne is portable on all Linux flavors and UNIX platforms tested so far.
It runs in parallel with MPI on distributed memory machines (Intel, Cray X series,
IBM Power, ...).
Developed since 1997 at EDF R&D, it is based on a co-located Finite Volume
approach that accepts meshes with any type of cell (tetrahedral, hexahedral,
prismatic, pyramidal, polyhedral...) and any type of grid structure
(unstructured, block structured, hybrid, conforming or with hanging nodes, ...).

Meshes may be imported using the CGNS, MED, GMSH, I-Deas, GAMBIT, or Simail
formats, and Post-processing output is available in EnSight, CGNS
and MED formats. In-situ postprocessing is available using the
[ParaView Catalyst](https://www.paraview.org/in-situ) and
[Melissa](https://melissa-sa.github.io) libraries.

code_saturne can be coupled in parallel to EDF's thermal software
[SYRTHES](https://www.edf.fr/en/the-edf-group/world-s-largest-power-company/activities/research-and-development/scientific-communities/simulation-softwares?logiciel=10818)
(conjugate heat transfer). It can also produce output usable by EDF's structural
analysis software [code_aster](https://code-aster.org), in particular in the
[SALOME platform](https://www.salome-platform.org/). SYRTHES and
code_aster are developed by EDF and distributed under a GNU GPL licence.
The atmospheric model can include chemistry modeling based on the
[SSH-aerosol](https://sshaerosol.wordpress.com/) library.

Copying
=======

code_saturne is distributed under the GNU General Public Licence, v2, or
(at your option) any later version.

See the COPYING file for details.

Installation
============

Installation from source files can be done either through a semi-automatic
`install_saturne.py` script, directly by using GNU-Autotools-based scripts,
or by a combination thereof.

Install instructions (from source) are provided in the
accompanying [INSTALL.md](INSTALL.md) file.

code_saturne support: saturne-support@edf.fr
