!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine visv2f &
!================

 ( rtp    , rtpa   , propce )

!===============================================================================
! FONCTION :
! --------

! CALCUL DE LA VISCOSITE TURBULENTE POUR
!          LE MODELE K-OMEGA V2F-BL


! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
use dimens, only: ndimfb
use numvar
use optcal
use cstphy
use entsor
use mesh
use field

!===============================================================================

implicit none

! Arguments

double precision rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)

! Local variables

integer          iel, iccocg, inc
integer          ipcvis, ipcvst

double precision s11, s22, s33
double precision dudy, dudz, dvdx, dvdz, dwdx, dwdy
double precision xk, xe, xrom, xnu
double precision ttke, ttmin, ttlim, tt

logical          ilved

double precision, allocatable, dimension(:) :: s2
double precision, dimension(:,:,:), allocatable :: gradv
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: crom

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)

! --- Memoire
allocate(s2(ncelet))

! --- Rang des variables dans PROPCE (prop. physiques au centre)
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)
call field_get_val_s(icrom, crom)

!===============================================================================
! 2.  CALCUL DES GRADIENTS DE VITESSE ET DE
!       S2 = 2* (S11**2+S22**2+S33**2+2*(S12**2+S13**2+S23**2)
!===============================================================================

! Allocate temporary arrays for gradients calculation
allocate(gradv(3,3,ncelet))

iccocg = 1
inc = 1

ilved = .false.

! WARNING: gradv(xyz, uvw, iel)
call grdvec &
!==========
( iu  , imrgra , inc    ,                               &
  nswrgr(iu) , imligr(iu) , iwarni(iu) ,                &
  epsrgr(iu) , climgr(iu) ,                             &
  ilved  ,                                              &
  rtpa(1,iu) ,  coefau , coefbu,                        &
  gradv  )

do iel = 1, ncel

  s11  = gradv(1,1,iel)
  s22  = gradv(2,2,iel)
  s33  = gradv(3,3,iel)
  dudy = gradv(2,1,iel)
  dudz = gradv(3,1,iel)
  dvdx = gradv(1,2,iel)
  dvdz = gradv(3,2,iel)
  dwdx = gradv(1,3,iel)
  dwdy = gradv(2,3,iel)

  s2(iel) = 2.d0*(s11**2 + s22**2 + s33**2)                   &
       + (dudy+dvdx)**2 + (dudz+dwdx)**2 + (dvdz+dwdy)**2
  s2(iel) = sqrt(max(s2(iel),1.d-10))

enddo

! Free memory
deallocate(gradv)

!===============================================================================
! 3.  CALCUL DE LA VISCOSITE
!===============================================================================

do iel = 1, ncel

  xk = rtp(iel,ik)
  xe = rtp(iel,iep)
  xrom = crom(iel)
  xnu = propce(iel,ipcvis)/xrom

  ttke = xk / xe
  ttmin = cpalct*sqrt(xnu/xe)
  ttlim = 0.6d0/rtp(iel,iphi)/sqrt(3.d0)/cpalmu/s2(iel)
  tt = min(ttlim,sqrt(ttke**2 + ttmin**2))

  propce(iel,ipcvst) = cpalmu*xrom*tt*rtp(iel,iphi)*rtp(iel,ik)

enddo

! Free memory
deallocate(s2)

!----
! FORMAT
!----


!----
! FIN
!----

return
end subroutine
