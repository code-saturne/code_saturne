!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

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

subroutine perinr &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  , iphas  ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idevel , ituser , ia     ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   rdevel , rtuser ,                                              &
   drdxyz  ,                                                      &
   ra     )

!===============================================================================
! FONCTION :
! --------

! PREPARATION DE LA PERIODICITE DE ROTATION POUR LES TENSIONS DE
!  REYNOLDS

!  ON CALCULE ICI UNE ESTIMATION DU GRADIENT DE RIJ, QUI
!    N'EST PAS UNE VARIABLE SCALAIRE, MAIS TENSORIELLE.
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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1...8(ncelet    ! tr ! --- ! tableau de travail                             !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! drdxyz           ! tr ! <-- ! gradient de r aux cellules halo pour           !
!                  !    !     ! l'approche explicite en periodicite            !
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

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas  , iphas
integer          nideve , nrdeve , nituse , nrtuse

integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet),w7(ncelet),w8(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision drdxyz(ncelet-ncel,6,3,nphas)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          inc, iccocg,ipiph,nswrgp,imligp,iwarnp
integer          isou, isou1, iphydp

double precision epsrgp, climgp,extrap

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

idebia = idbia0
idebra = idbra0

inc = 0
iccocg = 1

do isou = 1,6
  if(isou.eq.1) ipiph  = ir11(iphas)
  if(isou.eq.2) ipiph  = ir22(iphas)
  if(isou.eq.3) ipiph  = ir33(iphas)
  if(isou.eq.4) ipiph  = ir12(iphas)
  if(isou.eq.5) ipiph  = ir13(iphas)
  if(isou.eq.6) ipiph  = ir23(iphas)


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

!        NSWRGP = NSWRGR(IPIPH)
  nswrgp = 1
!        IMLIGP = IMLIGR(IPIPH)
  imligp = -1
  iwarnp = iwarni(ipiph)
  epsrgp = epsrgr(ipiph)
  climgp = climgr(ipiph)
  extrap = extrag(ipiph)
  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ipiph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   idevel , ituser , ia     ,                                     &
   w1     , w1     , w1     ,                                     &
   rtp(1,ipiph)  , coefa(1,ipiph) , coefb(1,ipiph) ,              &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )

  isou1 = isou
  call peinr1                                                     &
  !==========
  ( isou1  , iphas   ,                                            &
    drdxyz ,                                                      &
    w1     , w2      , w3     )

enddo

! --> ON FAIT TOURNER LE TENSEUR DRDXYZ PAR MANQUE DE TABLEAUX DE
!     TRAVAIL (ON A LE MEME PROBLEME POUR U)

call peinr2                                                       &
!==========
  ( iphas  , drdxyz )

! On a calcule les gradients dans DRDXYZ
igrper = 1


!----
! FIN
!----

return
end subroutine
