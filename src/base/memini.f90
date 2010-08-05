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

subroutine memini &
!================

 ( longia , longra ,                                              &
   nideve , nrdeve , nituse , nrtuse )

!===============================================================================
!  FONCTION
!  --------

!         GESTION MEMOIRE INITIALE
!         MACROS TABLEAUX, DEVELOPPEUR, UTILISATEUR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! longia longra    ! e  ! <-- ! longueur de ia     ra                          !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
!__________________.____._____.________________________________________________.

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
use entsor
use ihmpre

!===============================================================================

implicit none

! Arguments

integer longia , longra
integer nideve , nrdeve , nituse , nrtuse

! Local variables

integer ii, iok
integer icoftu(16)

!===============================================================================

!===============================================================================
! 0. INITIALISATIONS
!===============================================================================

! --- Dimension des macros tableaux entier IA et reel RA

longia = 0
longra = 0

! --- Dimension des tableaux utilisateurs et developpeurs
!      supplementaires (entiers et reels) generalement un multiple
!      de NCEL, NFABOR ...

nideve = 0
nrdeve = 0

nituse = 0
nrtuse = 0

do ii = 1, 16
  icoftu(ii) = 0
enddo

!===============================================================================
! 1. MEMOIRE INITIALE POUR UTILISATEUR, DEVELOPPEUR
!===============================================================================

!   - Interface Code_Saturne
!     ======================

if(iihmpr.eq.1) then

  call uiusar(icoftu)
  !==========

  nituse = icoftu(1)*ncelet + icoftu(2)*nfac + icoftu(3)*nfabor   &
         + icoftu(4)
  nrtuse = icoftu(5)*ncelet + icoftu(6)*nfac + icoftu(7)*nfabor   &
         + icoftu(8)
  longia = icoftu( 9)*ncelet + icoftu(10)*nfac + icoftu(11)*nfabor&
         + icoftu(12)
  longra = icoftu(13)*ncelet + icoftu(14)*nfac + icoftu(15)*nfabor&
         + icoftu(16)

endif

!   - Sous-programme utilisateur
!     ==========================

call ustbtr                                                       &
!==========
 ( ncel   , ncelet , nfac   , nfabor , nnod  ,                    &
   longia , longra ,                                              &
   nideve , nituse , nrdeve , nrtuse )


!===============================================================================
! 2. VERIFICATIONS
!===============================================================================

iok = 0

if (longia.lt.0) then
  WRITE(NFECRA,1000)'LONGIA',LONGIA
  iok = iok + 1
endif
if (longra.lt.0) then
  WRITE(NFECRA,1000)'LONGRA',LONGIA
  iok = iok + 1
endif

if (nideve.lt.0) then
  WRITE(NFECRA,1000)'NIDEVE',NIDEVE
  iok = iok + 1
endif
if (nrdeve.lt.0) then
  WRITE(NFECRA,1000)'NRDEVE',NRDEVE
  iok = iok + 1
endif
if (nituse.lt.0) then
  WRITE(NFECRA,1000)'NITUSE',NITUSE
  iok = iok + 1
endif
if (nrtuse.lt.0) then
  WRITE(NFECRA,1000)'NRTUSE',NRTUSE
  iok = iok + 1
endif

if (iok.gt.0) call csexit(1)
              !==========

!===============================================================================
! 3. DIMENSIONNEMENT MEMOIRE PAR DEFAUT
!===============================================================================

if (longia.eq.0) then
  longia =  30 * max(ncelet, nfabor)
endif

if (longra.eq.0) then
  longra = 120 * max(ncelet, nfabor)
endif

!----
! FORMATS
!----


 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER POSITIF OU NUL               ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou ustbtr.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

return
end subroutine
