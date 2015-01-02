!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

subroutine diverv &
!================

 ( div    , ux     , vy     , wz     ,                            &
   coefax , coefay , coefaz ,                                     &
   coefbx , coefby , coefbz )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!        CALCUL DE LA DIVERGENCE D'UN VECTEUR

!   (On ne s'embete pas, on appelle 3 fois le gradient)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! div(ncelet)      ! tr ! --> ! divergence du vecteur                          !
! ux,uy,uz         ! tr ! --> ! composante du vecteur                          !
! (ncelet)         !    !     !                                                !
! coefax,...       ! tr ! ->  ! conditions aux limites pour les                !
! coefbz           !    !     ! faces de bord                                  !
! (nfabor)         !    !     !                                                !
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
use cstphy
use cstnum
use parall
use period
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

double precision div(ncelet)
double precision ux(ncelet) , vy(ncelet) , wz(ncelet)
double precision coefax(nfabor) , coefay(nfabor) , coefaz(nfabor)
double precision coefbx(nfabor) , coefby(nfabor) , coefbz(nfabor)

! Local variables

integer          ivar0
integer          iel
integer          inc, iccocg
integer          nswrgp, imligp, iwarnp

double precision epsrgp, climgp, extrap

double precision, allocatable, dimension(:,:) :: gradu, gradv, gradw

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate work arrays
allocate(gradu(ncelet,3), gradv(ncelet,3), gradw(ncelet,3))


! En periodique et parallele, echange avant calcul du gradient
if (irangp.ge.0.or.iperio.eq.1) then
  call synvec(ux, vy, wz)
  !==========
endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas Rij)
ivar0 = 0

inc = 1
iccocg = 1
nswrgp = 100
imligp = -1
iwarnp = 2
epsrgp = 1.d-8
climgp = 1.5d0
extrap = 0.d0

!===============================================================================
! 1. Calcul du gradient de UX DANS W1
!===============================================================================

call grdcel                                                       &
!==========
( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,           &
  iwarnp , nfecra , epsrgp , climgp , extrap ,                    &
  ux     , coefax , coefbx ,                                      &
  gradu  )

!===============================================================================
! 2. Calcul du gradient de VY DANS W2
!===============================================================================

call grdcel                                                       &
!==========
( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,           &
  iwarnp , nfecra , epsrgp , climgp , extrap ,                    &
  vy     , coefay , coefby ,                                      &
  gradv  )

!===============================================================================
! 3. Calcul du gradient de VZ DANS W3
!===============================================================================

call grdcel                                                       &
!==========
( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,           &
  iwarnp , nfecra , epsrgp , climgp , extrap ,                    &
  wz     , coefaz , coefbz ,                                      &
  gradw  )

!===============================================================================
! 4. Calcul de la divergence du vecteur (UX,VY,WZ)
!===============================================================================

do iel = 1,ncel
  div(iel) = gradu(iel,1) + gradv(iel,2) + gradw(iel,3)
enddo

! Free memory
deallocate(gradu, gradv, gradw)

!----
! FIN
!----

end subroutine
