/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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
*/

/*-----------------------------------------------------------------------------*/

/*!
  \page remarks Remarks
 
  \section note_1 Note 1: ghost cells - (halos)
  A cell (real cell) is an elementary mesh element of the spatial
  discretisation of the calculation domain. The mesh is made of \ref ncel cells.
  When using periodicity and parallelism, extra "ghost" cells
  (also called "halo" cells) are defined for temporary storage of some information
  (on a given process).
  The total number of real and ghost cells is \ref ncelet. 
   - Indeed, when periodicity is enabled, the cells with
  periodic faces do not have any real neighboring cell across these
  particular faces. Their neighboring cell is elsewhere in the calculation
  domain (its position is determined by the periodicity). In order to
  temporarily store the information coming from this ``distant''
  neighboring cell, a ghost cell (``halo'') is created. 
   - The same kind of problem exists in the case of a
  calculation on parallel machines: due to the decomposition of the
  calculation domain, some cells no longer have access to all
  their neighboring cells, some of them being assigned to another parallel domain.
  The creation of ghost cells allows to temporarily store the information
  coming from real neighboring cells treated by other MPI ranks.
  The variables are generally arrays of size \ref ncelet (number of real and
  fictitious cells). The calculations (loops) are made on \ref ncel cells (only
  the real cells, the fictitious cells are only used to store information).


  \section note_2 Note 2: internal faces
  An internal face is an inferface shared by two cells (real or ghost
  ones) of the mesh. A boundary face is a face which has only one real
  neighboring cell. In the case of periodic calculations, a periodic face
  is an internal face. In the case of parallel running calculations, the
  faces situated at the boundary of a partition may be internal faces or
  boundary faces (of the whole mesh);

  \section note_3 Note 3: faces-vertices connectivity
  The faces - vertices connectivity is accessed by
  means of four integer functions: \ref ipnfac and \ref nodfac for the
  internal faces, \ref ipnfbr and \ref nodfbr for the boundary faces.
  \ref nodfac accesses the list of all the vertices of all the internal faces;
  first the vertices of the first face, then the vertices of the second face,
  and so on.
  \ref ipnfac (size: \ref nfac+1) gives the position \ref ipnfac "ipnfac"(ifac)
  in \ref nodfac of the first node of each internal face \c ifac.
  Therefore, the reference numbers of all
  the vertices of the internal face \c ifac are:
  \ref nodfac "nodfac"(\ref ipnfac "ipnfac"(ifac)),
  \ref nodfac "nodfac"(\ref ipnfac "ipnfac"(ifac)+1), ...,
  \ref nodfac "nodfac"(\ref ipnfac "ipnfac"(ifac+1)-1).
  In order for this last formula to be valid even for \c ifac=nfac,
  \ref ipnfac is of size \ref nfac+1 and \ref ipnfac(nfac+1)
  is equal to \ref lndfac+1.
  For boundary faces, the array access functions \ref nodfbr and \ref ipnfbr
  are used in a similar fashion.

  \section note_4 Note 4: modules
  The user must not modify the existing modules in user source directories,
  as this may break the code in more or less subtle fashion.
  Developers modifying modules must recompile the whole code. As this may
  not be possible for users (and in any case breaks versioning and quality
  assurance), modules should be defined in \ref cs_user_modules.f90.

*/
