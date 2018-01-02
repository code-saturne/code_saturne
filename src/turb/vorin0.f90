!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

subroutine vorin0 &
 ( nfabor )

!===============================================================================
!  FONCTION  :
!  ---------

! INITIALISATION PAR DEFAUT DES PARAMETRES POUR LA METHODE
!   DES VORTEX AVANT DE PASSER LA MAIN A L'UTILISATEUR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nfabor           ! i  ! <-- ! number of boundary faces                       !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use vorinc

!===============================================================================

implicit none

! Arguments

integer          nfabor

! Local variables

integer          ii, jj, ifac

!===============================================================================


!===============================================================================
! 1. PARAMETRES DE LA METHODE
!===============================================================================

! ---> Nombre d'entrees

nnent = -999

! ---> Nombre de vortex

do ii = 1, nentmx
  nvort(ii) = -999
enddo

! ---> Defintion du cas

do ii = 1, nentmx
  icas(ii)   = - 999
enddo

! ---> Tableau de reperage et repere local

do ifac = 1, nfabor
  irepvo(ifac) = 0
enddo

do ii = 1, nentmx
  do jj = 1, 3
    dir1(jj,ii) = 0.d0
    dir2(jj,ii) = 0.d0
    cen(jj,ii)  = 0.d0
  enddo
enddo

! ---> Conditions aux limites et grandeurs caracteristiques

do ii = 1, nentmx
  do jj = 1, 4
    iclvor(jj,ii) = -999
  enddo
  lly(ii) = -999.d0
  llz(ii) = -999.d0
  lld(ii) = -999.d0
enddo

! ---> Parametres physiques

do ii = 1, nentmx
  itlivo(ii) = -999
  tlimvo(ii) = -999.d0
  isgmvo(ii) = -999
  xsgmvo(ii) = -999.d0
  idepvo(ii) = -999
  ud(ii)     = 0.d0
enddo

! ---> Donnees utilisateurs

do ii = 1, nentmx
  WRITE(FICVOR(II),'(1A6,I2.2)') 'vordat',II
  udebit(ii) = 0.d0
  kdebit(ii) = -999.d0
  edebit(ii) = -999.d0
enddo

return
end subroutine
