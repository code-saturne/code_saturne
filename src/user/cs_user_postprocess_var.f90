!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!===============================================================================
! Purpose:
! -------

!> \file cs_user_postprocess_var.f90
!> \brief Output additional variables on a postprocessing mesh.
!>
!> Several "automatic" postprocessing meshes may be defined:
!> - The volume mesh (ipart=-1)
!> - The boundary mesh (ipart=-2)
!> - SYRTHES coupling surface (ipart < -2)
!>
!> Additional meshes (cells or faces) may also be defined through the GUI or
!> using the \ref cs_user_postprocess_meshes function from the
!> cs_user_postprocess.c file.
!>
!> This subroutine is called once for each post-processing mesh
!> (with a different value of 'ipart') for each time step at which output
!> on this mesh is active.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ipart         number of the post-processing mesh (< 0 or > 0)
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     nignor        ignored (set to 0, kept for argument list
!                               compatibility)
!> \param[in]     ncelps        number of cells in post-processing mesh
!> \param[in]     nfacps        number of interior faces in post-process. mesh
!> \param[in]     nfbrps        number of boundary faces in post-process. mesh
!> \param[in]     itypps        global presence flag (0 or 1) for cells (1),
!>                              interior faces (2), or boundary faces (3) in
!>                              post-processing mesh
!> \param[in]     lstcel        list of cells in post-processing mesh
!> \param[in]     lstfac        list of interior faces in post-processing mesh
!> \param[in]     lstfbr        list of boundary faces in post-processing mesh
!_______________________________________________________________________________

subroutine usvpst &
 ( ipart  ,                                                       &
   nvar   , nscal  , nignor ,                                     &
   ncelps , nfacps , nfbrps ,                                     &
   itypps ,                                                       &
   lstcel , lstfac , lstfbr )

!===============================================================================

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use pointe
use entsor
use optcal
use numvar
use parall
use period
use mesh
use field
use post
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          ipart
integer          nvar,   nscal , nignor
integer          ncelps, nfacps, nfbrps

integer          itypps(3)
integer          lstcel(ncelps), lstfac(nfacps), lstfbr(nfbrps)

!===============================================================================

return

end subroutine usvpst
