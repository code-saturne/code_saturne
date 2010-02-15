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

subroutine elprop &
!================

 ( ipropp , ipppst )

!===============================================================================
!  FONCTION  :
!  ---------

!     INIT DES POSITIONS DES VARIABLES D'ETAT POUR
!                LE MODULE ELECTRIQUE
!         (DANS VECTEURS PROPCE, PROPFA, PROPFB)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ipropp           ! e  ! ->  ! numero de la derniere case utlisee             !
!                  !    !     ! dans ipproc, ipprob, ipprof                    !
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

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "dimens.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "entsor.h"
include "cstnum.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "elincl.h"

!===============================================================================

! Arguments

integer       ipropp , ipppst

! Local variables

integer       iprop, idimve

!===============================================================================
!===============================================================================
! 1. DEFINITION DES POINTEURS
!===============================================================================

!     Pointeurs dans propce (ca n'implique pas qu'on ne calcule pas
!     les variables non definies ici)

iprop = ipropp

! ---> Temperature en K

iprop  = iprop + 1
itemp  = iprop

! ---> Puissance volumique dissipee par effet Joule W/m3

iprop  = iprop + 1
iefjou = iprop

! ---> Densite de courant electrique reelle A/m2

do idimve = 1, ndimve
  iprop        = iprop + 1
  idjr(idimve) = iprop
enddo

! Variables specifiques Effet Joule
! =================================

if ( ippmod(ieljou).eq.4 ) then

! ---> Densite de courant electrique imaginaire A/m2

  do idimve = 1, ndimve
    iprop        = iprop + 1
    idji(idimve) = iprop
  enddo

endif


! Variables specifiques Arc Electrique
! ====================================

if ( ippmod(ielarc).ge.1 ) then

! ---> Forces electromagnetiques de Laplace en N/m3

  do idimve = 1, ndimve
    iprop          = iprop + 1
    ilapla(idimve) = iprop
  enddo

! ---> Puissance volumique rayonnee W/m3
!      ou coefficient d'absorption

  if ( ixkabe .gt.0 ) then
    iprop = iprop + 1
    idrad = iprop
  endif
endif

! Variables specifiques Conduction Ionique
! ========================================

if ( ippmod(ielion).ge.1 ) then

! ---> Charge electrique volumique C/m3

  iprop  = iprop + 1
  iqelec = iprop

endif

! ----  Nb de variables algebriques (ou d'etat)
!         propre a la physique particuliere NSALPP
!         total NSALTO

nsalpp = iprop - ipropp
nsalto = iprop

! ----  On renvoie IPROPP au cas ou d'autres proprietes devraient
!         etre numerotees ensuite

ipropp = iprop

!===============================================================================
! 2. POSITIONNEMENT DES PROPRIETES : PROPCE, PROPFB, PROPFA
!===============================================================================

! ---> Positionnement dans le tableau PROPCE

iprop         = nproce

iprop         = iprop + 1
ipproc(itemp) = iprop
ipppst        = ipppst + 1
ipppro(iprop) = ipppst

iprop          = iprop + 1
ipproc(iefjou) = iprop
ipppst         = ipppst + 1
ipppro(iprop)  = ipppst

do idimve = 1, ndimve
  iprop                = iprop + 1
  ipproc(idjr(idimve)) = iprop
  ipppst               = ipppst + 1
  ipppro(iprop)        = ipppst
enddo

if ( ippmod(ieljou).eq.4 ) then

! ---> Densite de courant electrique imaginaire A/m2

  do idimve = 1, ndimve
    iprop                = iprop + 1
    ipproc(idji(idimve)) = iprop
    ipppst               = ipppst + 1
    ipppro(iprop)        = ipppst
  enddo

endif

if ( ippmod(ielarc).ge.1 ) then

  do idimve = 1, ndimve
    iprop                  = iprop + 1
    ipproc(ilapla(idimve)) = iprop
    ipppst                 = ipppst + 1
    ipppro(iprop)          = ipppst
  enddo

  if ( ixkabe .gt. 0 ) then
    iprop          = iprop + 1
    ipproc(idrad)  = iprop
    ipppst         = ipppst + 1
    ipppro(iprop)  = ipppst
  endif

endif

if ( ippmod(ielion).ge.1 ) then

! ---> Charge electrique volumique C/m3

  iprop          = iprop + 1
  ipproc(iqelec) = iprop
  ipppst         = ipppst + 1
  ipppro(iprop)  = ipppst

endif

nproce = iprop


!===============================================================================
! 3. POSITIONNEMENT DES PROPRIETES : PROPFB, PROPFA
!===============================================================================

! ---> Positionnement dans le tableau PROPFB
!      Au centre des faces de bord

  iprop = nprofb
  nprofb = iprop


! ---> Positionnement dans le tableau PROPFA
!      Au centre des faces internes (flux de masse)

 iprop = nprofa
 nprofa = iprop

return
end subroutine
