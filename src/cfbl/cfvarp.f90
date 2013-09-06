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

subroutine cfvarp
!================


!===============================================================================
!  FONCTION  :
!  ---------

!              INIT DES POSITIONS DES VARIABLES
!            POUR LE COMPRESSIBLE SANS CHOC SELON
! REMPLISSAGE DES PARAMETRES (DEJA DEFINIS) POUR LES SCALAIRES PP

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
use ihmpre

!===============================================================================

implicit none

! Local variables

integer          ii, iprop, iccfth, imodif
double precision dblpre(1)

!===============================================================================
!===============================================================================
! 1. DEFINITION DES POINTEURS
!===============================================================================


if ( ippmod(icompf).ge.0 ) then

  iprop =0

! ---- Energie totale
  iprop = iprop + 1
  ienerg = iscapp(iprop)
!     Alias pour les C.L.
  irunh = ienerg

! ---- Temperature (post)
  iprop = iprop + 1
  itempk = iscapp(iprop)

! ---- Viscosite dynamique de reference relative au scalaire ITEMPK
  ivisls(itempk) = 0
  visls0(itempk) = epzero

! ---- Viscosite dynamique de reference relative au scalaire IENERG
  ivisls(ienerg) = 0
  visls0(ienerg) = epzero

! ---- Initialisation par defaut de la viscosite en volume (cste)
  iviscv = 0
  viscv0 = 0.d0


!===============================================================================
! 2. OPTIONS DE CALCUL
!===============================================================================

! --> Cv constant ou variable (par defaut : constant)
  icv = 0
  cv0 = 0.d0

  iccfth = -1
  imodif = 0
  ii     = 1
  dblpre(1) = 0.d0
  call cfther                                                     &
  !==========
 ( ii ,                                                           &
   iccfth , imodif  ,                                             &
   dblpre , dblpre , dblpre , dblpre , dblpre ,                   &
   dblpre , dblpre , dblpre , dblpre , dblpre , dblpre )

! --> Utilisation d'un flux de masse specifique pour la vitesse

!     ATTENTION   PAS ENCORE IMPLEMENTE
!========   LAISSER IFLMAU = 0

  iflmau = 0

!===============================================================================
! 3. ON REDONNE LA MAIN A L'UTILISATEUR
!===============================================================================
!   - Interface Code_Saturne
!     ======================
!     Construction de l'indirection entre la numerotation du noyau et XML
  if (iihmpr.eq.1) then
    call uicfsc(ienerg, itempk)
    call csvvva(iviscv)
  endif

endif

!--------
! FORMATS
!--------


return
end subroutine

