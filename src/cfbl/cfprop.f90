!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine cfprop &
!================

 ( ipropp , ipppst )

!===============================================================================
!  FONCTION  :
!  ---------

!     INIT DES POSITIONS DES VARIABLES D'ETAT POUR
!              POUR LE COMPRESSIBLE SANS CHOC
!         (DANS VECTEURS PROPCE)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ipropp           ! e  ! <-- ! numero de la derniere propriete                !
!                  !    !     !  (les proprietes sont dans propce)             !
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
use ppincl

!===============================================================================

implicit none

! Arguments

integer       ipropp, ipppst

! Local variables

integer       iprop, ipp

!===============================================================================
!===============================================================================
! 1. POSITIONNEMENT DES PROPRIETES : PROPCE
!    Physique particuliere : Compressible sans choc
!===============================================================================

if ( ippmod(icompf).ge.0 ) then

! ---> Definition des pointeurs relatifs aux variables d'etat


  iprop = ipropp

!  Proprietes des phases : CV s'il est variable
  if(icv.ne.0) then
    iprop         = iprop + 1
    icv    = iprop
  endif

!  Proprietes des phases : Viscosite en volume
  if(iviscv.ne.0) then
    iprop         = iprop + 1
    iviscv = iprop
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

  iprop = nproce

  if(icv.gt.0) then
    iprop                 = iprop + 1
    ipproc(icv   ) = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst
  endif

  if(iviscv.gt.0) then
    iprop                 = iprop + 1
    ipproc(iviscv) = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst
  endif

  nproce = iprop

!===============================================================================
! 2. ENTREES SORTIES (entsor.h)
!===============================================================================

!     Comme pour les autres variables,
!       si l'on n'affecte pas les tableaux suivants,
!       les valeurs par defaut seront utilisees

!     nomprp( ) = nom de la propriete
!     ichrvr( ) = sortie chono (oui 1/non 0)
!     ilisvr( ) = suivi listing (oui 1/non 0)
!     ihisvr( ) = sortie historique (nombre de sondes et numeros)
!     si ihisvr(.,1)  = -1 sortie sur toutes les sondes

!     NB : Seuls les 8 premiers caracteres du nom seront repris dans le
!          listing le plus detaille

  !-->  chaleur specifique a volume constant
  if (icv   .gt.0) then
    ipp = ipppro(ipproc(icv))
    nomprp(ipproc(icv))   = 'Specific Heat Cst Vol'
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = 0
  endif

  !-->  viscosite laminaire
  if (iviscv.gt.0) then
    ipp = ipppro(ipproc(iviscv))
    nomprp(ipproc(iviscv))   = 'Volume Viscosity'
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = 0
  endif

endif


return
end subroutine
