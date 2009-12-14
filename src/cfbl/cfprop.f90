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

subroutine cfprop &
!================

 ( ipropp , ipppst )

!===============================================================================
!  FONCTION  :
!  ---------

!     INIT DES POSITIONS DES VARIABLES D'ETAT POUR
!              POUR LE COMPRESSIBLE SANS CHOC
!         (DANS VECTEURS PROPCE, PROPFA, PROPFB)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
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

implicit none

!===============================================================================
!     DONNEES EN COMMON
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

!===============================================================================

! Arguments

integer       ipropp, ipppst

! VARIABLES LOCALES

integer       iprop, ipp, iphas

!===============================================================================
!===============================================================================
! 1. POSITIONNEMENT DES PROPRIETES : PROPCE, PROPFA, PROPFB
!    Physique particuliere : Compressible sans choc
!===============================================================================

if ( ippmod(icompf).ge.0 ) then

! ---> Definition des pointeurs relatifs aux variables d'etat


  iprop = ipropp

!  Proprietes des phases : CV s'il est variable
  do iphas = 1, nphas
    if(icv(iphas).ne.0) then
      iprop         = iprop + 1
      icv   (iphas) = iprop
    endif
  enddo

!  Proprietes des phases : Viscosite en volume
  do iphas = 1, nphas
    if(iviscv(iphas).ne.0) then
      iprop         = iprop + 1
      iviscv(iphas) = iprop
    endif
  enddo

!   Flux de masse specifique pour la vitesse (si on en veut un)
  do iphas = 1, nphas
    if(iflmau(iphas).gt.0) then
      iprop         = iprop + 1
      ifluma(iu  (iphas)) = iprop
      ifluma(iv  (iphas)) = iprop
      ifluma(iw  (iphas)) = iprop
    endif
  enddo

!    Flux de Rusanov au bord pour Qdm et E
  do iphas = 1, nphas
    iprop         = iprop + 1
    ifbrhu(iphas) = iprop
    iprop         = iprop + 1
    ifbrhv(iphas) = iprop
    iprop         = iprop + 1
    ifbrhw(iphas) = iprop
    iprop         = iprop + 1
    ifbene(iphas) = iprop
  enddo


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

  do iphas = 1, nphas

    if(icv(iphas).gt.0) then
      iprop                 = iprop + 1
      ipproc(icv   (iphas)) = iprop
      ipppst                = ipppst + 1
      ipppro(iprop)         = ipppst
    endif

    if(iviscv(iphas).gt.0) then
      iprop                 = iprop + 1
      ipproc(iviscv(iphas)) = iprop
      ipppst                = ipppst + 1
      ipppro(iprop)         = ipppst
    endif

  enddo

  nproce = iprop


! ---> Positionnement dans le tableau PROPFB
!      Au centre des faces de bord

  iprop = nprofb

  do iphas = 1, nphas
    iprop                 = iprop + 1
    ipprob(ifbrhu(iphas)) = iprop
    iprop                 = iprop + 1
    ipprob(ifbrhv(iphas)) = iprop
    iprop                 = iprop + 1
    ipprob(ifbrhw(iphas)) = iprop
    iprop                 = iprop + 1
    ipprob(ifbene(iphas)) = iprop
  enddo

  nprofb = iprop


! ---> Positionnement dans le tableau PROPFA
!      Au centre des faces internes (flux de masse)

  iprop = nprofa

  do iphas = 1, nphas
    if(iflmau(iphas).gt.0) then
      iprop                     = iprop + 1
      ipprof(ifluma(iu(iphas))) = iprop
    endif
  enddo

  nprofa = iprop

!===============================================================================
! 2. ENTREES SORTIES (entsor.h)
!===============================================================================

!     Comme pour les autres variables,
!       si l'on n'affecte pas les tableaux suivants,
!       les valeurs par defaut seront utilisees

!     NOMVAR( ) = nom de la variable
!     ICHRVR( ) = sortie chono (oui 1/non 0)
!     ILISVR( ) = suivi listing (oui 1/non 0)
!     IHISVR( ) = sortie historique (nombre de sondes et numeros)
!     si IHISVR(.,1)  = -1 sortie sur toutes les sondes

!     NB : Seuls les 8 premiers caracteres du nom seront repris dans le
!          listing le plus detaille

!-->  chaleur specifique a volume constant
  if(icv   (iphas).gt.0) then
    ipp = ipppro(ipproc(icv   (iphas)))
    NOMVAR(IPP)   = 'Specific Heat Cst Vol'
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = 0
  endif

!-->  viscosite laminaire
  if(iviscv(iphas).gt.0) then
    ipp = ipppro(ipproc(iviscv(iphas)))
    NOMVAR(IPP)   = 'Volume Viscosity'
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = 0
  endif

endif


return
end subroutine
