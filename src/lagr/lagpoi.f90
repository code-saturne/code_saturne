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

subroutine lagpoi &
!================

 ( nbpmax , nvp    , nivep  ,                            &
   nvlsta ,                                              &
   itepa  ,                                              &
   ettp   , statis )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     RESOLUTION DE L'EQUATION DE POISSON POUR LES VITESSE MOYENNES
!                 DES PARTICULES
!       ET CORRECTION DES VITESSES INSTANTANNEES
!                 DES PARTICULES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! itepa            ! te ! --> ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! statis           ! tr ! <-- ! moyennes statistiques                          !
!(ncelet,nvlsta    !    !     !                                                !
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
use optcal
use entsor
use cstphy
use cstnum
use pointe
use parall
use period
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          nbpmax , nvp    , nivep
integer          nvlsta

integer          itepa(nbpmax,nivep)

double precision ettp(nbpmax,nvp)
double precision statis(ncelet,nvlsta)

! Local variables

integer          npt , iel , ifac
integer          ivar0
integer          inc, iccocg
integer          nswrgp , imligp , iwarnp

double precision epsrgp , climgp , extrap

double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: phil
double precision, allocatable, dimension(:) :: coefap, coefbp

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================


!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Allocate a temporary array
allocate(phil(ncelet))

do iel=1,ncel
  if ( statis(iel,ilpd) .gt. seuil ) then
    statis(iel,ilvx) = statis(iel,ilvx) / statis(iel,ilpd)
    statis(iel,ilvy) = statis(iel,ilvy) / statis(iel,ilpd)
    statis(iel,ilvz) = statis(iel,ilvz) / statis(iel,ilpd)
    statis(iel,ilfv) = statis(iel,ilfv) / ( dble(npst) * volume(iel) )
  else
    statis(iel,ilvx) = 0.d0
    statis(iel,ilvy) = 0.d0
    statis(iel,ilvz) = 0.d0
    statis(iel,ilfv) = 0.d0
  endif
enddo

call lageqp                                                       &
!==========
 ( statis(1,ilvx)  , statis(1,ilvy)  , statis(1,ilvz)  ,          &
   statis(1,ilfv)  ,                                              &
   phil   )

! Calcul du gradient du Correcteur PHI
! ====================================

! Allocate temporary arrays
allocate(coefap(nfabor))
allocate(coefbp(nfabor))

do ifac = 1, nfabor
  iel = ifabor(ifac)
  coefap(ifac) = phil(iel)
  coefbp(ifac) = zero
enddo

inc = 1
iccocg = 1
nswrgp = 100
imligp = -1
iwarnp = 2
epsrgp = 1.d-8
climgp = 1.5d0
extrap = 0.d0

! Allocate a work array
allocate(grad(ncelet,3))

! En periodique et parallele, echange avant calcul du gradient
if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(phil)
  !==========
endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
ivar0 = 0

call grdcel                                                       &
!==========
 ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   phil   , coefap , coefbp ,                                     &
   grad   )

! Free memory
deallocate(phil)
deallocate(coefap, coefbp)

! CORRECTION DES VITESSES MOYENNES ET RETOUR AU CUMUL

do iel = 1,ncel
  if ( statis(iel,ilpd) .gt. seuil ) then
    statis(iel,ilvx) = statis(iel,ilvx) - grad(iel,1)
    statis(iel,ilvy) = statis(iel,ilvy) - grad(iel,2)
    statis(iel,ilvz) = statis(iel,ilvz) - grad(iel,3)
  endif
enddo

do iel = 1,ncel
  if ( statis(iel,ilpd) .gt. seuil ) then
    statis(iel,ilvx) = statis(iel,ilvx)*statis(iel,ilpd)
    statis(iel,ilvy) = statis(iel,ilvy)*statis(iel,ilpd)
    statis(iel,ilvz) = statis(iel,ilvz)*statis(iel,ilpd)
    statis(iel,ilfv) = statis(iel,ilfv)*( dble(npst) * volume(iel) )
  endif
enddo

! CORRECTION DES VITESSES INSTANTANNES

do npt = 1,nbpart
  if (itepa(npt,jisor).gt.0) then
    iel = itepa(npt,jisor)
    ettp(npt,jup) = ettp(npt,jup) - grad(iel,1)
    ettp(npt,jvp) = ettp(npt,jvp) - grad(iel,2)
    ettp(npt,jwp) = ettp(npt,jwp) - grad(iel,3)
  endif
enddo

! Free memory
deallocate(grad)

!===============================================================================

!--------
! FORMATS
!--------

!----
! FIN
!----

end subroutine
