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

subroutine coprop &
!================

 ( ipropp , ipppst )

!===============================================================================
!  FONCTION  :
!  ---------

!     INIT DES POSITIONS DES VARIABLES D'ETAT POUR
!              POUR LA COMBUSTION
!        FLAMME DE DIFFUSION ET DE PREMELANGE
!         (DANS VECTEURS PROPCE, PROPFB)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ipropp           ! e  ! <-- ! numero de la derniere propriete                !
!                  !    !     !  (les proprietes sont dans propce ou prpb)     !
! ipppst           ! e  ! <-- ! pointeur indiquant le rang de la               !
!                  !    !     !  derniere grandeur definie aux                 !
!                  !    !     !  cellules (rtp,propce...) pour le              !
!                  !    !     !  post traitement                               !
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
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use radiat
use ihmpre

!===============================================================================

implicit none

! Arguments

integer       ipropp, ipppst

! Local variables

integer       iprop, icg, idirac

!===============================================================================
!===============================================================================
! 1. POSITIONNEMENT DES PROPRIETES : PROPCE
!    Physique particuliere : Flamme de diffusion chimie 3 points
!===============================================================================

if ( ippmod(icod3p).ge.0 ) then

! ---> Definition des pointeurs relatifs aux variables d'etat


  iprop = ipropp

  iprop  = iprop + 1
  itemp  = iprop
  do icg = 1, ngazg
    iprop    = iprop + 1
    iym(icg) = iprop
  enddo
  if ( ippmod(icod3p).eq.1 .and. iirayo.gt.0 ) then
    iprop = iprop + 1
    ickabs= iprop
    iprop = iprop + 1
    it4m  = iprop
    iprop = iprop + 1
    it3m  = iprop
  endif


! ----  Nb de variables algebriques (ou d'etat)
!         propre a la physique particuliere NSALPP
!         total NSALTO

  nsalpp = iprop - ipropp
  nsalto = iprop

! ----  On renvoie IPROPP au cas ou d'autres proprietes devraient
!         etre numerotees ensuite

  ipropp = iprop


! ---> Positionnement dans le tableau PROPCE
!      et reperage du rang pour le post-traitement

  iprop                 = nproce

  iprop                 = iprop + 1
  ipproc(itemp)         = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  do icg = 1, ngazg
    iprop                 = iprop + 1
    ipproc(iym(icg))      = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst
  enddo

  if ( ippmod(icod3p).eq.1 .and. iirayo.gt.0 ) then

    iprop                 = iprop + 1
    ipproc(ickabs)        = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    iprop                 = iprop + 1
    ipproc(it4m)          = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    iprop                 = iprop + 1
    ipproc(it3m)          = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

  endif

  nproce = iprop

 endif


!===============================================================================
! 3. POSITIONNEMENT DES PROPRIETES : PROPCE
!    Physique particuliere : Flamme de premelange - Modele EBU
!===============================================================================

if ( ippmod(icoebu).ge.0 ) then

! ---> Definition des pointeurs relatifs aux variables d'etat

  iprop = ipropp

  iprop  = iprop + 1
  itemp  = iprop
  do icg = 1, ngazg
    iprop    = iprop + 1
    iym(icg) = iprop
  enddo
  if ( ( ippmod(icoebu).eq.1 .or. ippmod(icoebu).eq.3 )           &
       .and. ( iirayo.gt.0 )  ) then
    iprop = iprop + 1
    ickabs= iprop
    iprop = iprop + 1
    it4m  = iprop
    iprop = iprop + 1
    it3m  = iprop
  endif


! ----  Nb de variables algebriques (ou d'etat)
!         propre a la physique particuliere NSALPP
!         total NSALTO

  nsalpp = iprop - ipropp
  nsalto = iprop

! ----  On renvoie IPROPP au cas ou d'autres proprietes devraient
!         etre numerotees ensuite

  ipropp = iprop


! ---> Positionnement dans le tableau PROPCE
!      et reperage du rang pour le post-traitement

  iprop         = nproce

  iprop                 = iprop + 1
  ipproc(itemp)         = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  do icg = 1, ngazg
    iprop                 = iprop + 1
    ipproc(iym(icg))      = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst
  enddo

  if ( ( ippmod(icoebu).eq.1 .or. ippmod(icoebu).eq.3 )           &
       .and. ( iirayo.gt.0 )  ) then

    iprop                 = iprop + 1
    ipproc(ickabs)        = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    iprop                 = iprop + 1
    ipproc(it4m)          = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    iprop                 = iprop + 1
    ipproc(it3m)          = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

  endif

  nproce = iprop

 endif


!===============================================================================
! 5. POSITIONNEMENT DES PROPRIETES : PROPCE
!    Physique particuliere : Flamme de premelange - Modele LWC
!===============================================================================

if ( ippmod(icolwc).ge.0 ) then

! ---> Definition des pointeurs relatifs aux variables d'etat

  iprop = ipropp

  iprop  = iprop + 1
  itemp  = iprop

  iprop  = iprop + 1
  imam   = iprop

  iprop  = iprop + 1
  itsc   = iprop

  do icg = 1, ngazg
    iprop    = iprop + 1
    iym(icg) = iprop
  enddo

  do idirac = 1, ndirac
! masse volumique Locale
    iprop         = iprop + 1
    irhol(idirac) = iprop
! Temperature L    .
    iprop         = iprop + 1
    iteml(idirac) = iprop
! Fraction de Melange L.
    iprop         = iprop + 1
    ifmel(idirac) = iprop
! Fraction Massique L.
    iprop         = iprop + 1
    ifmal(idirac) = iprop
! Amplitude L.
    iprop         = iprop + 1
    iampl(idirac) = iprop
! Terme Source Chimique L.
    iprop         = iprop + 1
    itscl(idirac) = iprop
! MAsse Molaire L.
    iprop         = iprop + 1
    imaml(idirac) = iprop
  enddo

  if ( ( ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3             &
                             .or. ippmod(icolwc).eq.5 )           &
       .and. ( iirayo.gt.0 )  ) then
    iprop = iprop + 1
    ickabs= iprop
    iprop = iprop + 1
    it4m  = iprop
    iprop = iprop + 1
    it3m  = iprop
  endif

! ----  Nb de variables algebriques (ou d'etat)
!         propre a la physique particuliere NSALPP
!         total NSALTO

  nsalpp = iprop - ipropp
  nsalto = iprop

! ----  On renvoie IPROPP au cas ou d'autres proprietes devraient
!         etre numerotees ensuite

  ipropp = iprop


! ---> Positionnement dans le tableau PROPCE

  iprop         = nproce

  iprop         = iprop + 1
  ipproc(itemp) = iprop
  ipppst        = ipppst + 1
  ipppro(iprop) = ipppst

  iprop         = iprop + 1
  ipproc(imam)  = iprop
  ipppst        = ipppst + 1
  ipppro(iprop) = ipppst

  iprop         = iprop + 1
  ipproc(itsc)  = iprop
  ipppst        = ipppst + 1
  ipppro(iprop) = ipppst

  do icg = 1, ngazg
    iprop            = iprop + 1
    ipproc(iym(icg)) = iprop
    ipppst           = ipppst + 1
    ipppro(iprop)    = ipppst
  enddo

  do idirac = 1, ndirac
    iprop                 = iprop + 1
    ipproc(irhol(idirac)) = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    iprop                 = iprop + 1
    ipproc(iteml(idirac)) = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    iprop                 = iprop + 1
    ipproc(ifmel(idirac)) = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    iprop                 = iprop + 1
    ipproc(ifmal(idirac)) = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    iprop                 = iprop + 1
    ipproc(iampl(idirac)) = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    iprop                 = iprop + 1
    ipproc(itscl(idirac)) = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    iprop                 = iprop + 1
    ipproc(imaml(idirac)) = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

  enddo

  if ( ( ippmod(icolwc).eq.1 .or. ippmod(icolwc).eq.3             &
                             .or. ippmod(icolwc).eq.5 )           &
       .and. ( iirayo.gt.0 )  ) then

    iprop                 = iprop + 1
    ipproc(ickabs)        = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    iprop                 = iprop + 1
    ipproc(it4m)          = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    iprop                 = iprop + 1
    ipproc(it3m)          = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

  endif

  nproce = iprop

endif

!   - Interface Code_Saturne
!     ======================
!     Construction de l'indirection entre la numerotation du noyau et XML
if (iihmpr.eq.1) then
  call uicopr (nsalpp, ippmod, ipppro, ipproc,                       &
               icolwc, iirayo, itemp, imam, iym, ickabs, it4m, it3m, &
               itsc, irhol, iteml, ifmel, ifmal, iampl, itscl,       &
               imaml)

endif

return

end subroutine
