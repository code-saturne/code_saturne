/*============================================================================
 * Main Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of the PLE (Parallel Location and Exchange) library.

  Copyright (C) 2011  EDF

  The PLE library is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

  The PLE library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this library; if not, write to the Free Software Foundation, Inc.,
  51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*-----------------------------------------------------------------------------*/

/*!
  \mainpage PLE (Parallel Location and Exchange) documentation
 
  \section intro Introduction

  PLE is a libary designed to simplify coupling of distributed parallel
  computational codes. It is maintained as a part of Code_Saturne,
  EDF's general purpose Computational Fluid Dynamics (CFD) software,
  but it may also be used with other tools, and is distributed under
  a broader licence (LGPL instead of GPL).

  PLE provides support for 2 categories of tasks: synchronizing
  parallel codes at predifined points, and enabling parallel mapping
  of points to meshes, and transfer of variables using this mapping.

  \section install Installation

  PLE may be configured and installed using the \c configure shell script
  and \c make. Please read the \c INSTALL file in the toplevel source
  directory if you are not familiar with configuration scripts generated
  through GNU autoconf and automake.

  \section FVM

  The "Finite Volume Mesh" library is intended to provide mesh and associated
  fields I/O and manipulation services for unstructure Finite Volume codes
  or other tools with similar requirements.

  \subsection PLE Coupling API

  The ple_coupling_...() functions allow identifying applications and
  defining MPI communicators necessary to the ple_locator_...()
  functions, as well as providing each of a set of coupled codes
  with info on the other code's time steps, convergence status, and
  other synchronization data at predifined points (usually once per
  time step).

  \subsection PLE Locator subset

  The ple_locator_...() functions allow mapping points to a mesh
  in parallel, given serial functions providing this functionnality
  for the associated data structures, then exchanging variables using
  this mapping.

*/
