!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
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

subroutine perinu &
!================

 ( nvar   , nscal  ,                                              &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   dudxyz )

!===============================================================================
! FONCTION :
! --------

! PREPARATION DE LA PERIODICITE DE ROTATION POUR LA VITESSE

!  ON CALCULE ICI UNE ESTIMATION DU GRADIENT DE LA VITESSE, QUI
!    N'EST PAS UNE VARIABLE SCALAIRE, MAIS VECTORIELLE.
!  LE GRADIENT EST PRIS SANS RECONSTRUCTION (SINON, ON A BESOIN
!    DU GRADIENT, ET DONC DE LA PERIODICITE). IL SEMBLE POSSIBLE
!    D'UTILISER GRADMC.
!  LE GRADIENT EST ENSUITE STOCKE OU IL FAUT DANS UN TABLEAU
!    REPRESENTANT LE HALO PUIS
!    SOUMIS A ROTATION LA OU C'EST NECESSAIRE POUR ETRE PRET
!    A L'EMPLOI (CF PERING).


! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! dudxyz           ! tr ! <-- ! gradient de u aux cellules halo pour           !
!                  !    !     ! l'approche explicite en periodicite            !
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
use cstnum
use entsor
use pointe
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal


double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision dudxyz(ncelet-ncel,3,3)

! Local variables

integer          inc, iccocg, ipiph, nswrgp, imligp, iwarnp
integer          isou, isou1

double precision epsrgp, climgp,extrap

double precision, allocatable, dimension(:,:) :: grad

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Allocate a work array
allocate(grad(ncelet,3))


inc = 0
iccocg = 1

do isou = 1,3
  if(isou.eq.1) ipiph  = iu
  if(isou.eq.2) ipiph  = iv
  if(isou.eq.3) ipiph  = iw


! On ne reconstruit pas et on ne limite pas car
!   on ne connait pas a priori le gradient des voisins (justement,
!   on le calcule ici). En fait, a partir du second pas de temps relatif
!   on dispose des valeurs calculees precedemment, et on pourrait donc les
!   utiliser pour reconstruire : a tester. (le pas de temps relatif est compte
!   a partir du premier pas de temps du "run" courant, donc a partir de la
!   lecture du fichier suite eventuellement)

! Attention, on precise bien qu'a la sortie de grdcel ci dessous, le halo
!   des rotations contient :
!       - rien du tout au premier pas de temps relatif
!       - sinon, le gradient calcule au pas de temps precedent
! On fera donc attention a ne pas utiliser le halo ici dans le cas general.

!       NSWRGP = NSWRGR(IPIPH)
  nswrgp = 1
!        IMLIGP = IMLIGR(IPIPH)
  imligp = -1
  iwarnp = iwarni(ipiph)
  epsrgp = epsrgr(ipiph)
  climgp = climgr(ipiph)
  extrap = extrag(ipiph)

  call grdcel                                                     &
  !==========
 ( ipiph  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rtp(1,ipiph)  , coefa(1,ipiph) , coefb(1,ipiph) ,              &
   grad   )

  isou1 = isou
  call peinu1                                                     &
  !==========
  ( isou1  ,                                                      &
    dudxyz ,                                                      &
    grad(1,1) , grad(1,2) , grad(1,3) )

enddo

! --> ON FAIT TOURNER LE TENSEUR DUDXYZ PAR MANQUE DE TABLEAUX DE
!     TRAVAIL (ON A LE MEME PROBLEME POUR RIJ)

call peinu2 ( dudxyz )
!==========

! On a calcule les gradients dans DUDXYZ
iguper = 1

! Free memory
deallocate(grad)

!----
! FIN
!----

return
end subroutine
