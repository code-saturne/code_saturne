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

subroutine ppprcl &
!================

 ( nvar   , nscal  ,                                              &
   icodcl , izfppp ,                                              &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   rcodcl )

!===============================================================================
! FONCTION :
! --------

!    PREPARATION DU REMPLISSAGE DES CONDITIONS AUX LIMITES

!           AIGUILLAGE SPECIFIQUE AUX PHYSIQUES PARTICULIERES


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
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
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use cfpoin
use atincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvar)
integer          izfppp(nfabor)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ifac, izone, icha, iclapc
integer          ivar

!===============================================================================
!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================


! ---> Combustion gaz USEBUC
!      Flamme de diffusion : chimie 3 points

if ( ippmod(icod3p).ge.0 ) then

  do izone = 1, nozppm
    qimp(izone)   = zero
    iqimp(izone)  = 0
    ientox(izone) = 0
    ientfu(izone) = 0
  enddo

  tinoxy        = zero
  tinfue        = zero

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

! ---> Combustion gaz USEBUC
!      Flamme de premelange : modele EBU

elseif ( ippmod(icoebu).ge.0 ) then

  do izone = 1, nozppm
    iqimp(izone)  = 0
    qimp(izone)   = zero
    icalke(izone) = 0
    dh(izone)     = zero
    xintur(izone) = zero
    fment(izone)  = zero
    tkent(izone)  = zero
    ientgf(izone) = 0
    ientgb(izone) = 0
  enddo

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo


! ---> Combustion charbon pulverise USCPCL

elseif ( ippmod(icp3pl).ge.0 ) then

  do izone = 1, nozppm
    iqimp(izone)  = 0
    icalke(izone) = 0
    ientcp(izone) = 0
    ientat(izone) = 0
    dh(izone)     = zero
    xintur(izone) = zero
    qimpat(izone) = zero
    timpat(izone) = zero
    do icha = 1, ncharm
      qimpcp(izone,icha) = zero
      timpcp(izone,icha) = zero
      do iclapc = 1, ncpcmx
        distch(izone,icha,iclapc) = zero
      enddo
    enddo
  enddo

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

elseif ( ippmod(iccoal).ge.0 ) then

  do izone = 1, nozppm
    iqimp(izone)  = 0
    icalke(izone) = 0
    ientcp(izone) = 0
    ientat(izone) = 0
    dh(izone)     = zero
    xintur(izone) = zero
    qimpat(izone) = zero
    timpat(izone) = zero
    do icha = 1, ncharm
      qimpcp(izone,icha) = zero
      timpcp(izone,icha) = zero
      do iclapc = 1, ncpcmx
        distch(izone,icha,iclapc) = zero
      enddo
    enddo
  enddo

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

! ---> Combustion charbon pulverise couple Lagrangien USCPLC

elseif ( ippmod(icpl3c).ge.0 ) then

  do izone = 1, nozppm
    iqimp(izone)  = 0
    icalke(izone) = 0
    ientat(izone) = 0
    dh(izone)     = zero
    xintur(izone) = zero
    qimpat(izone) = zero
    timpat(izone) = zero
    do icha = 1, ncharm
      qimpcp(izone,icha) = zero
    enddo
  enddo

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

! ---> Combustion fuel  USFUCL

elseif ( ippmod(icfuel).ge.0 ) then

  do izone = 1, nozppm
    iqimp(izone)  = 0
    icalke(izone) = 0
    ientat(izone) = 0
    dh(izone)     = zero
    xintur(izone) = zero
    qimpat(izone) = zero
    timpat(izone) = zero
    ientfl(izone) = 0
    qimpfl(izone) = zero
    timpfl(izone) = zero
  enddo

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

! ---> Compressible

elseif ( ippmod(icompf).ge.0 ) then

!     Zones
  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

!     Marqueur d'utilisation de Rusanov au bord (0 = non)
!     Marqueur de flux conductif impose au bord (0 = non)
  do ifac = 1, nfabor
    ifbrus(ifac) = 0
    ifbet(ifac) = 0
  enddo

!     Flux de Rusanov au bord pour Qdm et E
  do ifac = 1, nfabor
    propfb(ifac,ipprob(ifbrhu)) = 0.d0
    propfb(ifac,ipprob(ifbrhv)) = 0.d0
    propfb(ifac,ipprob(ifbrhw)) = 0.d0
    propfb(ifac,ipprob(ifbene)) = 0.d0
  enddo

!     Initialisation des RCODCL(IFAC,.,1) à -RINFIN
!       pour savoir si l'utilisateur les a modifies (ils sont
!       initialises par defaut à 0)
  do ivar = 1, nvar
    do ifac = 1, nfabor
      rcodcl(ifac,ivar,1) =-rinfin
    enddo
  enddo

! ---> Version electrique
!      Effet Joule
!      Conduction ionique

elseif ( ippmod(ieljou).ge.1 .or.                                 &
     ippmod(ielarc).ge.1 .or.                                     &
     ippmod(ielion).ge.1       ) then

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

! ---> Version ecoulements atmospheriques

elseif ( ippmod(iatmos).ge.0  ) then

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo
  do izone = 1, nozppm
    iprofm(izone) = 0
  enddo

!     Initialisation des RCODCL(IFAC,.,1) à RINFIN
!       pour savoir si l'utilisateur les a modifies (ils sont
!       initialises par defaut à 0)
  do ivar = 1, nvar
    do ifac = 1, nfabor
      rcodcl(ifac,ivar,1) = rinfin
    enddo
  enddo

! ---> Version aerorefrigerants

elseif ( ippmod(iaeros).ge.0 ) then

  do ifac = 1, nfabor
    izfppp(ifac) = 0
  enddo

endif

!----
! FORMATS
!----


!----
! FIN
!----

return
end subroutine
