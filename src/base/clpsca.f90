!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine clpsca &
!================

 ( ncelet , ncel   , iscal  , scandd , rtp    )

!===============================================================================
! FONCTION :
! ----------

! CLIPPING
!   POUR UN SCALAIRE OU VARIANCE

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! iscal            ! i  ! <-- ! scalar number                                  !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant        )          !
! scandd           ! tr ! <-- ! scalaire auquel est associe la                 !
! (ncelet)         !    !     !    variance traitee (si c'en est une)          !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

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

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel
integer          iscal

double precision rtp(ncelet,nflown:nvar)
double precision scandd(ncelet)

! Local variables

integer          ivar, iel, iflid
integer          iclmax, iclmin, iiscav
integer          kscmin, kscmax, f_id
double precision vmin(1), vmax(1), vfmin, vfmax
double precision scmax, scmin
double precision scmaxp, scminp

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! --- Numero de variable de calcul et de post associe au scalaire traite
ivar   = isca(iscal)
iflid  = ivarfl(ivar)

! --- Numero du scalaire eventuel associe dans le cas fluctuation
iiscav = iscavr(iscal)

! Key ids for clipping
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

!===============================================================================
! 2. IMPRESSIONS ET CLIPPINGS
!===============================================================================

!      IL Y A TOUJOURS CLIPPING DES VARIANCES A DES VALEURS POSITIVES

! --- Calcul du min et max
vmin(1) = rtp(1,ivar)
vmax(1) = rtp(1,ivar)
do iel = 1, ncel
  vmin(1) = min(vmin(1),rtp(iel,ivar))
  vmax(1) = max(vmax(1),rtp(iel,ivar))
enddo

if (iiscav.eq.0) then

  ! Clipping of non-variance scalars

  iclmax = 0
  iclmin = 0

  ! Get the min clipping
  call field_get_key_double(iflid, kscmin, scminp)
  call field_get_key_double(iflid, kscmax, scmaxp)

  if(scmaxp.gt.scminp)then
    do iel = 1, ncel
      if(rtp(iel,ivar).gt.scmaxp)then
        iclmax = iclmax + 1
        rtp(iel,ivar) = scmaxp
      endif
      if(rtp(iel,ivar).lt.scminp)then
        iclmin = iclmin + 1
        rtp(iel,ivar) = scminp
      endif
    enddo
  endif

else

  ! Clipping of variances

  f_id = ivarfl(isca(iiscav))

  iclmax = 0
  iclmin = 0

  ! -- Clipping minimal au minimum 0.
  if(iclvfl(iscal).eq.0) then
    do iel = 1, ncel
      if(rtp(iel,ivar).lt.0.d0) then
        iclmin = iclmin + 1
        rtp(iel,ivar) = 0.d0
      endif
    enddo

  ! -- Clipping a partir des valeurs du scalaire (ou 0 au min)
  elseif(iclvfl(iscal).eq.1) then
    do iel = 1, ncel
      if(rtp(iel,ivar).lt.0.d0) then
        iclmin = iclmin + 1
        rtp(iel,ivar) = 0.d0
      endif
    enddo

    ! Get the min clipping
    call field_get_key_double(f_id, kscmin, scmin)
    call field_get_key_double(f_id, kscmax, scmax)

    do iel = 1, ncel
      vfmax = (scandd(iel)-scmin)*(scmax-scandd(iel))
      if(rtp(iel,ivar).gt.vfmax) then
        iclmax = iclmax + 1
        rtp(iel,ivar) = vfmax
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
        if(rtp(iel,ivar).gt.vfmax)then
          iclmax = iclmax + 1
          rtp(iel,ivar) = vfmax
        endif
        if(rtp(iel,ivar).lt.vfmin)then
          iclmin = iclmin + 1
          rtp(iel,ivar) = vfmin
        endif
      enddo
    endif
  endif

endif

call log_iteration_clipping_field(iflid, iclmin, iclmax, vmin, vmax)

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine
