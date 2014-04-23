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

subroutine vissst &
!================

 ( rtpa   , propce )

!===============================================================================
! FONCTION :
! --------

! CALCUL DE LA VISCOSITE TURBULENTE POUR
!          LE MODELE K-OMEGA SST

! VISCT = ROM * A1 * K /MAX(A1*W ; SQRT(S2KW)*F2)
! AVEC S2KW =  2 * Sij.Sij
!       Sij = (DUi/Dxj + DUj/Dxi)/2

! ET F2 = TANH(ARG2**2)
! ARG2**2 = MAX(2*SQRT(K)/CMU/W/Y ; 500*NU/W/Y**2)

! DIVU EST CALCULE EN MEME TEMPS QUE S2KW POUR ETRE REUTILISE
! DANS TURBKW

! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at previous time step)                       !
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
use cstnum
use pointe, only: s2kw, divukw, ifapat, dispar
use numvar
use optcal
use dimens, only: nvar
use cstphy
use entsor
use mesh
use field
use field_operator

!===============================================================================

implicit none

! Arguments

double precision rtpa(ncelet,nflown:nvar)
double precision propce(ncelet,*)

! Local variables

integer          iel, inc
integer          ipcvis, ipcvst
integer          nswrgp, imligp, iwarnp
integer          ifacpt, iprev

double precision d1s3, d2s3
double precision epsrgp, climgp, extrap
double precision xk, xw, rom, xmu, xdist, xarg2, xf2

double precision, allocatable, dimension(:) :: w1
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

! --- Rang des variables dans PROPCE (prop. physiques au centre)
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)
call field_get_val_s(icrom, crom)

d1s3 = 1.d0/3.d0
d2s3 = 2.d0/3.d0

!===============================================================================
! 2. Compute the scalar s2kw rate SijSij and the trace of the velocity
!    gradient

!      (Sij^D) (Sij^D)  is stored in    s2kw (deviatoric s2kw tensor rate)
!      tr(Grad u)       is stored in    divukw
!===============================================================================


! Allocate temporary arrays for gradients calculation
allocate(gradv(3,3,ncelet))

inc = 1
iprev = 1

call field_gradient_vector(ivarfl(iu), iprev, imrgra, inc,    &
                           gradv)

! s2kw = Stain rate of the deviatoric part of the s2kw tensor
!      = 2 (Sij^D).(Sij^D)
! divukw   = trace of the velocity gradient
!          = dudx + dvdy + dwdz

do iel = 1, ncel

  s2kw(iel) = 2.d0                                                           &
    *( ( d2s3*gradv(1,1,iel) - d1s3*gradv(2,2,iel) - d1s3*gradv(3,3,iel))**2   &
     + (-d1s3*gradv(1,1,iel) + d2s3*gradv(2,2,iel) - d1s3*gradv(3,3,iel))**2   &
     + (-d1s3*gradv(1,1,iel) - d1s3*gradv(2,2,iel) + d2s3*gradv(3,3,iel))**2   &
     )                                                                         &
    + (gradv(2,1,iel) + gradv(1,2,iel))**2                                     &
    + (gradv(3,1,iel) + gradv(1,3,iel))**2                                     &
    + (gradv(3,2,iel) + gradv(2,3,iel))**2

  divukw(iel) = gradv(1,1,iel) + gradv(2,2,iel) + gradv(3,3,iel)

enddo

! Free memory
deallocate(gradv)

!===============================================================================
! 3.  CALCUL DE LA DISTANCE A LA PAROI
!===============================================================================

! Allocate a work array
allocate(w1(ncelet))

if(abs(icdpar).eq.2) then
  do iel = 1 , ncel
    ifacpt = ifapat(iel)
    if (ifacpt.gt.0) then
      w1(iel) =  (cdgfbo(1,ifacpt)-xyzcen(1,iel))**2           &
               + (cdgfbo(2,ifacpt)-xyzcen(2,iel))**2           &
               + (cdgfbo(3,ifacpt)-xyzcen(3,iel))**2
      w1(iel) = sqrt(w1(iel))
    else
      w1(iel) = grand
    endif
  enddo
else
  do iel = 1 , ncel
    w1(iel) =  max(dispar(iel),epzero)
  enddo
endif

!===============================================================================
! 4.  CALCUL DE LA VISCOSITE
!===============================================================================

do iel = 1, ncel

  xk = rtpa(iel,ik)
  xw = rtpa(iel,iomg)
  rom = crom(iel)
  xmu = propce(iel,ipcvis)
  xdist = w1(iel)
  xarg2 = max (                                                   &
       2.d0*sqrt(xk)/cmu/xw/xdist,                                &
       500.d0*xmu/rom/xw/xdist**2 )
  xf2 = tanh(xarg2**2)

  propce(iel,ipcvst) = rom*ckwa1*xk                               &
       /max( ckwa1*xw , sqrt(s2kw(iel))*xf2 )

enddo

! Free memory
deallocate(w1)

!-------
! Format
!-------

!----
! End
!----

return
end subroutine
