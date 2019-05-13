!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

!> \file clpsca.f90
!> \brief This subroutine clips the values of a given scalar or variance.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iscal         scalar number
!_______________________________________________________________________________

subroutine clpsca &
 ( iscal )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use dimens
use cstphy
use cstnum
use parall
use field
use cs_c_bindings
use mesh

!===============================================================================

implicit none

! Arguments

integer          iscal

! Local variables

integer          ivar, iel, iflid
integer          iclmax(1), iclmin(1), iiscav
integer          kscmin, kscmax, f_id
double precision vmin(1), vmax(1), vfmin, vfmax
double precision scmax, scmin
double precision scmaxp, scminp

double precision, dimension(:), pointer :: cvar_scal, cvar_scav

!===============================================================================

!===============================================================================
! 1. Initialisation
!===============================================================================

! --- Numero de variable de calcul et de post associe au scalaire traite
ivar   = isca(iscal)
iflid  = ivarfl(ivar)

! --- Numero du scalaire eventuel associe dans le cas fluctuation
iiscav = iscavr(iscal)

! Map field arrays
call field_get_val_s(ivarfl(ivar), cvar_scal)

! Key ids for clipping
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

!===============================================================================
! 2. Printings and clippings
!===============================================================================

! Variances are always clipped at positive values

! Compute min. and max. values
vmin(1) = cvar_scal(1)
vmax(1) = cvar_scal(1)
do iel = 1, ncel
  vmin(1) = min(vmin(1),cvar_scal(iel))
  vmax(1) = max(vmax(1),cvar_scal(iel))
enddo

if (iiscav.eq.0) then

  ! Clipping of non-variance scalars

  iclmax(1) = 0
  iclmin(1) = 0

  ! Get the min clipping
  call field_get_key_double(iflid, kscmin, scminp)
  call field_get_key_double(iflid, kscmax, scmaxp)

  if(scmaxp.gt.scminp)then
    do iel = 1, ncel
      if(cvar_scal(iel).gt.scmaxp)then
        iclmax(1) = iclmax(1) + 1
        cvar_scal(iel) = scmaxp
      endif
      if(cvar_scal(iel).lt.scminp)then
        iclmin(1) = iclmin(1) + 1
        cvar_scal(iel) = scminp
      endif
    enddo
  endif

else

  ! Clipping of variances

  f_id = ivarfl(isca(iiscav))
  call field_get_val_s(f_id, cvar_scav)

  iclmax(1) = 0
  iclmin(1) = 0

  ! Minimal clipping at minimum 0.
  if(iclvfl(iscal).eq.0) then
    do iel = 1, ncel
      if(cvar_scal(iel).lt.0.d0) then
        iclmin(1) = iclmin(1) + 1
        cvar_scal(iel) = 0.d0
      endif
    enddo

  ! Clipping based on associated scalar values (or 0 at min.)
  elseif(iclvfl(iscal).eq.1) then
    do iel = 1, ncel
      if(cvar_scal(iel).lt.0.d0) then
        iclmin(1) = iclmin(1) + 1
        cvar_scal(iel) = 0.d0
      endif
    enddo

    ! Get the min clipping
    call field_get_key_double(f_id, kscmin, scmin)
    call field_get_key_double(f_id, kscmax, scmax)

    do iel = 1, ncel
      vfmax = (cvar_scav(iel)-scmin)*(scmax-cvar_scav(iel))
      if(cvar_scal(iel).gt.vfmax) then
        iclmax(1) = iclmax(1) + 1
        cvar_scal(iel) = vfmax
      endif
    enddo

!   -- Clipping a partir des valeurs donnees par l'utilisateur
!        (ou 0 au min)
  elseif(iclvfl(iscal).eq.2) then
    vfmin = 0.d0
    ! Get the min clipping
    call field_get_key_double(iflid, kscmin, scminp)
    call field_get_key_double(iflid, kscmax, scmaxp)
    vfmin = max(scminp,vfmin)
    vfmax = scmaxp
    if(vfmax.gt.vfmin)then
      do iel = 1, ncel
        if(cvar_scal(iel).gt.vfmax)then
          iclmax(1) = iclmax(1) + 1
          cvar_scal(iel) = vfmax
        endif
        if(cvar_scal(iel).lt.vfmin)then
          iclmin(1) = iclmin(1) + 1
          cvar_scal(iel) = vfmin
        endif
      enddo
    endif
  endif

endif

call log_iteration_clipping_field(iflid, iclmin(1), iclmax(1), vmin, vmax, iclmin(1), iclmax(1))

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine
