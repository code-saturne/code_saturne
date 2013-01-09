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

subroutine cpprop &
!================

 ( ipropp , ipppst )

!===============================================================================
!  FONCTION  :
!  ---------

! INIT DES POSITIONS DES VARIABLES D'ETAT SELON
!         COMBUSTION CHARBON PULVERISE
!   (DANS VECTEURS PROPCE, PROPFA, PROPFB)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ipropp           ! e  ! <-- ! numero de la derniere propriete                !
!                  !    !     !  (les proprietes sont dans propce,             !
!                  !    !     !   propfa ou prpfb)                             !
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
use ppcpfu
use coincl
use cpincl
use ppincl
use ihmpre

!===============================================================================

implicit none

! Arguments

integer       ipropp, ipppst

! Local variables

integer       iprop, ige , icla, iprop2

!===============================================================================

! ---> Definition des pointeurs relatifs aux variables d'etat

iprop = ipropp

!    Phase continue (melange gazeux)
iprop   = iprop + 1
itemp1  = iprop
iprop   = iprop + 1
irom1  = iprop
do ige = 1, (ngaze-2*ncharb)
! ---- Cf. definition de NGAZE dans cplecd.F
  iprop     = iprop + 1
  iym1(ige) = iprop
enddo
iprop = iprop + 1
immel = iprop

iprop2 = iprop

!   Phase dispersee (classes de particules)
do icla = 1, nclacp
  iprop        = iprop2 + icla
  itemp2(icla) = iprop
  iprop        = iprop2 + 1*nclacp + icla
  ix2(icla)    = iprop
  iprop        = iprop2 + 2*nclacp + icla
  irom2(icla)  = iprop
  iprop        = iprop2 + 3*nclacp + icla
  idiam2(icla) = iprop
  iprop        = iprop2 + 4*nclacp + icla
  igmdch(icla) = iprop
  iprop        = iprop2 + 5*nclacp + icla
  igmdv1(icla) = iprop
  iprop        = iprop2 + 6*nclacp + icla
  igmdv2(icla) = iprop
  iprop        = iprop2 + 7*nclacp + icla
  igmhet(icla) = iprop
  if ( ihtco2 .eq. 1 ) then
    iprop        = iprop2 + 8*nclacp + icla
    ighco2(icla) = iprop
    if ( ippmod(icp3pl) .eq. 1 ) then
      iprop        = iprop2 + 9*nclacp + icla
      igmsec(icla) = iprop
    endif
  else
    if ( ippmod(icp3pl) .eq. 1 ) then
      iprop        = iprop2 + 8*nclacp + icla
      igmsec(icla) = iprop
    endif
  endif
enddo

! ---- Nb de variables algebriques (ou d'etat)
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

!    Phase continue (melange gazeux)
iprop           = iprop + 1
ipproc(itemp1)  = iprop
ipppst          = ipppst + 1
ipppro(iprop)   = ipppst

iprop           = iprop + 1
ipproc(irom1)   = iprop
ipppst          = ipppst + 1
ipppro(iprop)   = ipppst

do ige = 1, (ngaze-2*ncharb)
! ---- Cf. definition de NGAZE dans cplecd.F
  iprop                 = iprop + 1
  ipproc(iym1(ige))     = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
enddo

iprop                 = iprop + 1
ipproc(immel)         = iprop
ipppst                = ipppst + 1
ipppro(iprop)         = ipppst

iprop2 = iprop

!   Phase dispersee (classes de particules)
do icla = 1, nclacp

  iprop                 = iprop2 + icla
  ipproc(itemp2(icla))  = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iprop2 + 1*nclacp + icla
  ipproc(ix2(icla))     = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iprop2 + 2*nclacp + icla
  ipproc(irom2(icla))   = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iprop2 + 3*nclacp + icla
  ipproc(idiam2(icla))  = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iprop2 + 4*nclacp + icla
  ipproc(igmdch(icla))  = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iprop2 + 5*nclacp + icla
  ipproc(igmdv1(icla))  = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iprop2 + 6*nclacp + icla
  ipproc(igmdv2(icla))  = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iprop2 + 7*nclacp + icla
  ipproc(igmhet(icla))  = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  if ( ihtco2 .eq. 1 ) then
    iprop                 = iprop2 + 8*nclacp + icla
    ipproc(ighco2(icla))  = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    if ( ippmod(icp3pl) .eq. 1 ) then
      iprop                 = iprop2 + 9*nclacp + icla
      ipproc(igmsec(icla))  = iprop
      ipppst                = ipppst + 1
      ipppro(iprop)         = ipppst
    endif
  else
    if ( ippmod(icp3pl) .eq. 1 ) then
      iprop                 = iprop2 + 8*nclacp + icla
      ipproc(igmsec(icla))  = iprop
      ipppst                = ipppst + 1
      ipppro(iprop)         = ipppst
    endif
  endif

enddo

nproce = iprop


! ---> Positionnement dans le tableau PROPFB
!      Au centre des faces de bord

iprop = nprofb
nprofb = iprop

! ---> Positionnement dans le tableau PROPFA
!      Au centre des faces internes (flux de masse)

iprop = nprofa
nprofa = iprop


!   - Interface Code_Saturne
!     ======================
!     Construction de l'indirection entre la numerotation du noyau et XML
if (iihmpr.eq.1) then
  call uicppr (nclacp, nsalpp, nsalto, ippmod, icp3pl, ipppro,    &
               ipproc, ihtco2, itemp1, irom1, iym1, immel,        &
               itemp2, ix2, irom2, idiam2, igmdch, igmdv1,        &
               igmdv2, igmhet, ighco2, igmsec)

endif

return
end subroutine
