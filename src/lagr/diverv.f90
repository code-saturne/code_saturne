!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine diverv &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   ia     ,                                                       &
   dt     ,                                                       &
   div    , ux     , vy     , wz     ,                            &
   coefax , coefay , coefaz ,                                     &
   coefbx , coefby , coefbz ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )

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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! div(ncelet)      ! tr ! --> ! divergence du vecteur                          !
! ux,uy,uz         ! tr ! --> ! composante du vecteur                          !
! (ncelet)         !    !     !                                                !
! coefax,...       ! tr ! ->  ! conditions aux limites pour les                !
! coefbz           !    !     ! faces de bord                                  !
! (nfabor)         !    !     !                                                !
! w1...9(ncelet    ! tr ! --- ! tableau de travail                             !
! ra(*)            ! ra ! --- ! main real work array                           !
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

integer          idbia0 , idbra0
integer          nvar   , nscal

integer          ia(*)

double precision dt(ncelet)
double precision div(ncelet)
double precision ux(ncelet) , vy(ncelet) , wz(ncelet)
double precision coefax(nfabor) , coefay(nfabor) , coefaz(nfabor)
double precision coefbx(nfabor) , coefby(nfabor) , coefbz(nfabor)
double precision w1(ncelet) , w2(ncelet) , w3(ncelet)
double precision w4(ncelet) , w5(ncelet) , w6(ncelet)
double precision w7(ncelet) , w8(ncelet) , w9(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ivar0
integer          iel
integer          inc, iccocg
integer          nswrgp, imligp, iwarnp
double precision epsrgp, climgp, extrap

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

! En periodique et parallele, echange avant calcul du gradient
if (irangp.ge.0.or.iperio.eq.1) then
  call synvec(ux, vy, wz)
  !==========
endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
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
  ia     ,                                                        &
  ux     , coefax , coefbx ,                                      &
  w1     , w4     , w5     ,                                      &
  w6     , w7     , w8     ,                                      &
  ra     )

!===============================================================================
! 2. Calcul du gradient de VY DANS W2
!===============================================================================

call grdcel                                                       &
!==========
( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,           &
  iwarnp , nfecra , epsrgp , climgp , extrap ,                    &
  ia     ,                                                        &
  vy     , coefay , coefby ,                                      &
  w4     , w2     , w5     ,                                      &
  w6     , w7     , w8     ,                                      &
  ra     )

!===============================================================================
! 3. Calcul du gradient de VZ DANS W3
!===============================================================================

call grdcel                                                       &
!==========
( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,           &
  iwarnp , nfecra , epsrgp , climgp , extrap ,                    &
  ia     ,                                                        &
  wz     , coefaz , coefbz ,                                      &
  w5     , w6     , w3     ,                                      &
  w7     , w8     , w9     ,                                      &
  ra     )

!===============================================================================
! 4. Calcul de la divergence du vecteur (UX,VY,WZ)
!===============================================================================

do iel = 1,ncel
  div(iel) = w1(iel) + w2(iel) + w3(iel)
enddo

!----
! FIN
!----

end subroutine
