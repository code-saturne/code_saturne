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
!>
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

! Local variables

integer          ntindp, f_id
integer          iel, ifac, iloc, ivar
integer          idimt, ii , jj
logical          ientla, ivarpr
integer          imom1, imom2
double precision pnd
double precision rvoid(1)

double precision, dimension(:), allocatable :: scel, sfac, sfbr
double precision, dimension(:,:), allocatable :: vcel, vfac, vfbr
double precision, dimension(:), pointer :: bfpro_rom
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cvar_var

double precision, dimension(:), pointer :: cmom_1, cmom_2

double precision, dimension(:,:), pointer :: cvar_vel

integer          intpst
data             intpst /0/
save             intpst

!===============================================================================

!===============================================================================
! Increment call counter once per time step (possibly used in some tests)
!===============================================================================

if (ipart .eq. -1) then
  intpst = intpst + 1
endif

!===============================================================================
! 1. Handle variables to output
!    MUST BE FILLED IN by the user at indicated places
!===============================================================================

! The ipart argument matches a post-processing maehs id (using the EnSight
! vocabulary; the MED and CGNS equivalents are "mesh" and "base" respectively).
! The user will have defined post-processing meshes using the GUI or the
! cs_user_postprocess_meshes() function from the cs_user_postprocess.c
! file.

! This subroutine is called once for each post-processing mesh
! (with a different value of 'ipart') for each time step at which output
! on this mesh is active. For each mesh and for all variables we wish to
! post-process here, we must define certain parameters and pass them to
! the 'post_write_var' subroutine, which is in charge of the actual output.
! These parameters are:

! namevr <-- variable name
! idimt  <-- variable dimension
!            (1: scalar, 3: vector, 6: symmetric tensor, 9: tensor)
! ientla <-- when idimt >1, this flag specifies if the array containing the
!            variable values is interlaced when ientla = .true.
!            (x1, y1, z1, x2, y2, z2, x3, y3, z3...), or non-interlaced when
!            ientla = .false. (x1, x2, x3,...,y1, y2, y3,...,z1, z2, z3,...).
! ivarpr <-- specifies if the array containing the variable is defined on
!            the "parent" mesh or locally.
!            Even if the 'ipart' post-processing mesh contains all the
!            elements of its parent mesh, their numbering may be different,
!            especially when different element types are present.
!            A local array passed as an argument to 'post_write_var' is built
!            relative to the numbering of the 'ipart' post-processing mesh.
!            To post-process a variable contained for example in the 'user'
!            array, it should first be re-ordered, as shown here:
!              do iloc = 1, ncelps
!                iel = lstcel(iloc)
!                scel(iloc) = user(iel)
!              enddo
!            An alternative option is provided, to avoid unnecessary copies:
!            an array defined on the parent mesh, such our 'user' example,
!            may be passed directly to 'post_write_var', specifying that values
!            are defined on the parent mesh instead of the post-processing mesh,
!            by setting the 'ivarpr' argument of 'post_write_var' to .true..

! Note: be cautious with variable name lengths.

! We allow up to 32 characters here, but names may be truncted depending on the
! output format.

! The name length is not limited internally, so in case of 2 variables whoses
! names differ only after the truncation character, the corresponding names will
! both appear in the ".case" file; simply renaming one of the field descriptors
! in this text file will correct the output.

! Whitespace at the beginning or the end of a line is truncated automatically.
! Depending on the format used, prohibited characters (under EnSight, characters
! (  ) ] [ + - @           ! # * ^ $ / as well as white spaces and tabulations
! are automatically replaced by the _ character.

! Examples:

!   For post-processing mesh 2, we output the velocity, pressure, and prescribed
!   temperature at boundary faces (as well as 0 on possible interior faces)

!   For post-processing mesh 1, we output all the variables usually
!   post-processed, using a more compact coding.

!   Examples given here correspond to the meshes defined in
!   cs_user_postprocess.c

!===============================================================================
! Examples for initialization
!===============================================================================

!< [postprocess_var_ex_3]

!===============================================================================
! Examples of volume variables on the boundary mesh (ipart = -2)
!===============================================================================

if (ipart .eq. -2) then

  ! Output of the density at the boundary
  ! -------------------------------------

  idimt = 1        ! 1: scalar, 3: vector, 6/9: symm/non-symm tensor
  ientla = .true.  ! dimension 1 here, so no effect
  ivarpr = .true.  ! we use the bfpro_rom array defined on the parent mesh

  ! Output values; as we have no cell or interior face values, we can pass a
  ! trivial array for those.
  call field_get_val_s(ibrom, bfpro_rom)
  call post_write_var(ipart, 'Density at boundary', idimt, ientla, ivarpr,    &
                      ntcabs, ttcabs, rvoid, rvoid, bfpro_rom)

endif

!< [postprocess_var_ex_3]

!< [postprocess_var_ex_4]

!===============================================================================
! Examples of volume variables on user meshes 1 or 2
!===============================================================================

if (ipart.eq.1 .or. ipart.eq.2) then
  ! Map field arrays
  call field_get_val_v(ivarfl(iu), cvar_vel)

  ! Output of the velocity
  ! ----------------------

  ! Compute variable values on interior faces.
  ! In this example, we use a simple linear interpolation.
  ! For parallel calculations, if neighbors are used, they must be synchronized
  ! first. This also applies for periodicity.

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvin(cvar_vel)
  endif

  allocate(vfac(3,nfacps), vfbr(3,nfbrps))

  do iloc = 1, nfacps

    ifac = lstfac(iloc)
    ii = ifacel(1, ifac)
    jj = ifacel(2, ifac)
    pnd = pond(ifac)

    vfac(1,iloc) = pnd  * cvar_vel(1,ii) + (1.d0 - pnd) * cvar_vel(1,jj)
    vfac(2,iloc) = pnd  * cvar_vel(2,ii) + (1.d0 - pnd) * cvar_vel(2,jj)
    vfac(3,iloc) = pnd  * cvar_vel(3,ii) + (1.d0 - pnd) * cvar_vel(3,jj)

  enddo

  ! Compute variable values on boundary faces.
  ! In this example, we use a simple copy of the adjacent cell value.

  do iloc = 1, nfbrps

    ifac = lstfbr(iloc)
    ii = ifabor(ifac)

    vfbr(1,iloc) = cvar_vel(1,ii)
    vfbr(2,iloc) = cvar_vel(2,ii)
    vfbr(3,iloc) = cvar_vel(3,ii)

  enddo

  idimt = 3        ! 1: scalar, 3: vector, 6/9: symm/non-symm tensor
  ientla = .true.  ! interleaved
  ivarpr = .false. ! defined on the work array, not on the parent

  ! Output values; as we have no cell values, we can pass a
  ! trivial array for those.
  call post_write_var(ipart, 'Interpolated velocity', idimt, ientla, ivarpr,  &
                      ntcabs, ttcabs, rvoid, vfac, vfbr)

  deallocate(vfac, vfbr)

!< [postprocess_var_ex_4]

!< [postprocess_var_ex_5]

  ! Output of the pressure
  ! ----------------------

  ! Variable number
  ivar = ipr
  call field_get_val_s(ivarfl(ivar), cvar_var)

  ! Compute variable values on interior faces.
  ! In this example, we use a simple linear interpolation.
  ! For parallel calculations, if neighbors are used, they must be synchronized
  ! first. This also applies for periodicity.

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(cvar_var)
  endif

  allocate(sfac(nfacps), sfbr(nfbrps))

  do iloc = 1, nfacps

    ifac = lstfac(iloc)
    ii = ifacel(1, ifac)
    jj = ifacel(2, ifac)
    pnd = pond(ifac)

    sfac(iloc) =           pnd  * cvar_var(ii)  &
                 + (1.d0 - pnd) * cvar_var(jj)

  enddo

  ! Compute variable values on boundary faces.
  ! In this example, we use a simple copy of the adjacent cell value.

  do iloc = 1, nfbrps

    ifac = lstfbr(iloc)
    ii = ifabor(ifac)

    sfbr(iloc) = cvar_var(ii)

  enddo

  idimt = 1        ! 1: scalar, 3: vector, 6/9: symm/non-symm tensor
  ientla = .true.  ! dimension 1 here, so no effect
  ivarpr = .false. ! defined on the work array, not on the parent

  ! Output values; as we have no cell values, we can pass a
  ! trivial array for those.
  call post_write_var(ipart, 'Interpolated pressure', idimt, ientla, ivarpr,  &
                      ntcabs, ttcabs, rvoid, sfac, sfbr)

  deallocate(sfac, sfbr)

!< [postprocess_var_ex_5]

  ! The examples below illustrate how to output a same variable in different
  ! ways (interlaced or not, using an indirection or not).

!< [postprocess_var_ex_6]

  ! Output of the centers of gravity, interlaced
  ! --------------------------------

  if (intpst.eq.1) then

    allocate(vfac(3,nfacps), vfbr(3,nfbrps))

    do iloc = 1, nfacps

      ifac = lstfac(iloc)

      vfac(1,iloc) = cdgfac(1, ifac)
      vfac(2,iloc) = cdgfac(2, ifac)
      vfac(3,iloc) = cdgfac(3, ifac)

    enddo

    ! Compute variable values on boundary faces

    do iloc = 1, nfbrps

      ifac = lstfbr(iloc)

      vfbr(1, iloc) = cdgfbo(1, ifac)
      vfbr(2, iloc) = cdgfbo(2, ifac)
      vfbr(3, iloc) = cdgfbo(3, ifac)

    enddo

    ! We assign a negative time step and output this variable once only
    ! to avoid duplicating it at each output time (assuming a fixed mesh).
    ntindp = -1

    idimt = 3        ! 1: scalar, 3: vector, 6/9: symm/non-symm tensor
    ientla = .true.  ! interleaved
    ivarpr = .false. ! defined on the work array, not on the parent

    ! Output values; as we have no cell values, we can pass a
    ! trivial array for those.
    call post_write_var(ipart, 'face cog (interlaced)', idimt,               &
                        ientla, ivarpr,                                      &
                        ntindp, ttcabs, rvoid, vfac, vfbr)

    deallocate(vfac, vfbr)

  endif

!< [postprocess_var_ex_6]

!< [postprocess_var_ex_7]

  ! Output of the centers of gravity, non-interlaced, time independent
  ! --------------------------------

  if (intpst.eq.1) then

    allocate(vfac(nfacps, 3), vfbr(nfbrps, 3))

    do iloc = 1, nfacps

      ifac = lstfac(iloc)

      vfac(iloc,1) = cdgfac(1, ifac)
      vfac(iloc,2) = cdgfac(2, ifac)
      vfac(iloc,3) = cdgfac(3, ifac)

    enddo

    ! Compute variable values on boundary faces

    do iloc = 1, nfbrps

      ifac = lstfbr(iloc)

      vfbr(iloc,1) = cdgfbo(1, ifac)
      vfbr(iloc,2) = cdgfbo(2, ifac)
      vfbr(iloc,3) = cdgfbo(3, ifac)

    enddo

    ! We assign a negative time step and output this variable once only
    ! to avoid duplicating it at each output time (assuming a fixed mesh).
    ntindp = -1

    idimt = 3         ! 1: scalar, 3: vector, 6/9: symm/non-symm tensor
    ientla = .false.  ! not interleaved
    ivarpr = .false.  ! defined on the work array, not on the parent

    ! Output values; as we have no cell values, we can pass a
    ! trivial array for those.
    call post_write_var(ipart, 'face cog (non interlaced)', idimt,           &
                        ientla, ivarpr,                                      &
                        ntindp, ttcabs, rvoid, vfac, vfbr)

    deallocate(vfac, vfbr)

  endif

!< [postprocess_var_ex_7]

!< [postprocess_var_ex_8]

  ! Output of the centers of gravity, with indirection (parent-based)
  ! --------------------------------

  if (intpst.eq.1) then

    ! We assign a negative time step and output this variable once only
    ! to avoid duplicating it at each output time (assuming a fixed mesh).
    ntindp = -1

    idimt = 3        ! 1: scalar, 3: vector, 6/9: symm/non-symm tensor
    ientla = .true.  ! interleaved
    ivarpr = .true.  ! defined on the parent

    ! Output values; as we have no cell values, we can pass a
    ! trivial array for those.
    call post_write_var(ipart, 'face cog (parent)', idimt, ientla, ivarpr,   &
                        ntindp, ttcabs, rvoid, cdgfac, cdgfbo)

  endif

!< [postprocess_var_ex_8]

!< [postprocess_var_ex_9]

!===============================================================================
! Examples of volume variables on user meshes 3 or 4
!===============================================================================

else if (ipart.ge.3 .and. ipart.le.4) then

  ! Map field arrays
  call field_get_val_v(ivarfl(iu), cvar_vel)

  ! Output of the velocity
  ! ----------------------

  ! Compute variable values on interior faces.
  ! In this example, we use a simple linear interpolation.
  ! For parallel calculations, if neighbors are used, they must be synchronized
  ! first. This also applies for periodicity.

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvin(cvar_vel)
  endif

  allocate(vfac(3,nfacps), vfbr(3,nfbrps))

  do iloc = 1, nfacps

    ifac = lstfac(iloc)
    ii = ifacel(1, ifac)
    jj = ifacel(2, ifac)
    pnd = pond(ifac)

    vfac(1,iloc) =            pnd  * cvar_vel(1,ii)   &
                    + (1.d0 - pnd) * cvar_vel(1,jj)
    vfac(2,iloc) =            pnd  * cvar_vel(2,ii)   &
                    + (1.d0 - pnd) * cvar_vel(2,jj)
    vfac(3,iloc) =            pnd  * cvar_vel(3,ii)   &
                    + (1.d0 - pnd) * cvar_vel(3,jj)

  enddo

  ! Compute variable values on boundary faces.
  ! In this example, we use a simple copy of the adjacent cell value.

  do iloc = 1, nfbrps

    ifac = lstfbr(iloc)
    ii = ifabor(ifac)

    vfbr(1,iloc) = cvar_vel(1,ii)
    vfbr(2,iloc) = cvar_vel(2,ii)
    vfbr(3,iloc) = cvar_vel(3,ii)

  enddo

  idimt = 3         ! 1: scalar, 3: vector, 6/9: symm/non-symm tensor
  ientla = .true.   ! interleaved
  ivarpr = .false.  ! defined on the work array

  ! Output values; as we have no cell values, we can pass a
  ! trivial array for those.
  call post_write_var(ipart, 'Velocity', idimt, ientla, ivarpr,              &
                      ntcabs, ttcabs, rvoid, vfac, vfbr)

  deallocate(vfac, vfbr)

!< [postprocess_var_ex_9]

!< [postprocess_var_ex_10]

  ! Output of the pressure
  ! ----------------------

  ! Variable number
  ivar = ipr
  call field_get_val_s(ivarfl(ivar), cvar_var)

  ! Compute variable values on interior faces.
  ! In this example, we use a simple linear interpolation.
  ! For parallel calculations, if neighbors are used, they must be synchronized
  ! first. This also applies for periodicity.

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(cvar_var)
  endif

  allocate(sfac(nfacps), sfbr(nfbrps))

  do iloc = 1, nfacps

    ifac = lstfac(iloc)
    ii = ifacel(1, ifac)
    jj = ifacel(2, ifac)
    pnd = pond(ifac)

    sfac(iloc)  =           pnd  * cvar_var(ii)   &
                  + (1.d0 - pnd) * cvar_var(jj)

  enddo

  ! Compute variable values on boundary faces.
  ! In this example, we use a simple copy of the adjacent cell value.

  do iloc = 1, nfbrps

    ifac = lstfbr(iloc)
    ii = ifabor(ifac)

    sfbr(iloc) = cvar_var(ii)

  enddo

  idimt = 1         ! 1: scalar, 3: vector, 6/9: symm/non-symm tensor
  ientla = .true.   ! interleaved
  ivarpr = .false.  ! defined on the work array

  ! Output values; as we have no cell values, we can pass a
  ! trivial array for those.
  call post_write_var(ipart, 'Pressure', idimt, ientla, ivarpr,              &
                      ntcabs, ttcabs, rvoid, sfac, sfbr)

  deallocate(sfac, sfbr)

endif ! end of test on post-processing mesh number

!< [postprocess_var_ex_10]

return

end subroutine usvpst
