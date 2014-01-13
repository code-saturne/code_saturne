/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2014 EDF S.A.

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
  When using periodicity and parallelism, extra ``ghost'' cells
  (called ``halo'' cells) are defined for temporary storage of some information
  (on a given processor).
  The total number of real and ghost cells is \ref ncelet. 
   - Indeed, when periodicity is enabled, the cells with
  periodic faces do not have any real neighbouring cell across these
  particular faces. Their neighbouring cell is elsewhere in the calculation
  domain (its position is determined by the periodicity). In order to
  temporarily store the information coming from this ``distant''
  neighbouring cell, a ghost cell (``halo'') is created. 
   - The same kind of problem exists in the case of a
  calculation on parallel machines: due to the decomposition of the
  calculation domain, some cells no longer have access to all
  their neighbouring cells, some of them being treated by another processor. The
  creation of ghost cells allows to temporarily store the information
  coming from real neighbouring cells treated by other processors.
  The variables are generally arrays of size \ref ncelet (number of real and
  fictitious cells). The calculations (loops) are made on \ref ncel cells (only
  the real cells, the fictitious cells are only used to store information).


  \section note_2 Note 2: internal faces
  An internal face is an inferface shared by two cells (real or ghost
  ones) of the mesh. A boundary face is a face which has only one real
  neighbouring cell. In the case of periodic calculations, a periodic face
  is an internal face. In the case of parallel running calculations, the
  faces situated at the boundary of a partition may be internal faces or
  boundary faces (of the whole mesh);

  \section note_3 Note 3: faces-nodes connectivity
  The faces - nodes connectivity is stored by
  means of four integer arrays: \ref ipnfac and \ref nodfac for the
  internal faces, \ref ipnfbr and \ref nodfbr for the boundary faces.
  \ref nodfac (size \ref lndfac)
  contains the list of all the nodes of all the internal faces; first the nodes of
  the first face, then the nodes of the second face, and so on.
  \ref ipnfac (size: \ref nfac+1) gives the position \ref ipnfac(ifac)
  in \ref nodfac of the first node of each internal face \ref ifac.
  Therefore, the reference numbers of all
  the nodes of the internal face \ref ifac are: \ref nodfac(ipnfac(ifac)),
  \ref nodfac(ipnfac(ifac)+1), ..., \ref nodfac(ipnfac(ifac+1)-1).
  In order for this last formula to be valid even for \ref ifac=nfac,
  \ref ipnfac is of size \ref nfac+1 and \ref ipnfac(nfac+1)
  is equal to \ref lndfac+1.
  The composition of the arrays \ref nodfbr and \ref ipnfbr is similar.

  \section note_4 Note 4: modules
  The user will not modify the existing modules. This would require the
  recompilation of the complete version, operation which is not allowed in
  standard use.

*/
