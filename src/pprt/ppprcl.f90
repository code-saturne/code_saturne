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

subroutine ppprcl &
!================

 ( nvar   ,                                                       &
   izfppp ,                                                       &
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
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
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
!                  !    !     !        cp*(viscls+visct/turb_schmidt)*gradt    !
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

integer          nvar

integer          izfppp(nfabor)

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
    icvfli(ifac) = 0
    ifbet(ifac) = 0
  enddo

! ---> Version electrique
!      Effet Joule
!      Conduction ionique

elseif ( ippmod(ieljou).ge.1 .or.                                 &
         ippmod(ielarc).ge.1       ) then

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
