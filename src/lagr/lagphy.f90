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

subroutine lagphy &
!================

 ( ntersl , nvlsta , nvisbr ,                                     &
   dt     , rtp    , propce ,                                     &
   taup   , tlag   ,                                              &
   tempct , tsvar  , auxl   ,                                     &
   cpgd1  , cpgd2  , cpght  )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     INTEGRATION DES EDS CONCERNANT LES PHYSIQUES PARTICULIERES
!       LIEES AUX PARTICULES :

!         - Temperature du fluide vu par les particules,
!         - Temperature des particules,
!         - Diametre des particules
!         - Masse des particules
!         - Variables liees aux grains de charbon (Temp,MCH,MCK),
!         - Variables Utilisateur supplementaires.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! taup(nbpmax)     ! tr ! <-- ! temps caracteristique dynamique                !
! tlag(nbpmax)     ! tr ! <-- ! temps caracteristique fluide                   !
! tempct           ! tr ! <-- ! temps caracteristique thermique                !
!  (nbpmax,2)      !    !     !                                                !
! tsvar            ! tr ! <-- ! prediction 1er sous-pas pour la                !
! (nbpmax,nvp1)    !    !     !   variable courante, utilise pour la           !
!                  !    !     !   correction au 2eme sous-pas                  !
! auxl(nbpmax,3    ! tr ! --- ! tableau de travail lagrangien                  !
! cpgd1,cpgd2,     ! tr ! --> ! termes de devolatilisation 1 et 2 et           !
!  cpght(nbpmax    !    !     !   de combusion heterogene (charbon             !
!                  !    !     !   avec couplage retour thermique)              !
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
use cstphy
use cstnum
use optcal
use dimens, only: nvar
use entsor
use lagdim, only: nbpmax, nvp1
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          ntersl , nvlsta , nvisbr

double precision dt(ncelet) , rtp(ncelet,nflown:nvar)
double precision propce(ncelet,*)
double precision taup(nbpmax) , tlag(nbpmax,3) , tempct(nbpmax,2)
double precision tsvar(nbpmax,nvp1) , auxl(nbpmax,3)
double precision cpgd1(nbpmax) , cpgd2(nbpmax) , cpght(nbpmax)

! Local variables

!===============================================================================

!===============================================================================
! 0.  GESTION MEMOIRE
!===============================================================================


!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================


!===============================================================================
! 2. INTEGRATION DE LA TEMPERATURE FLUIDE VU PAR LES PARTICULES
!===============================================================================

if ( iphyla.eq.2 .or. (iphyla.eq.1 .and. itpvar.eq.1) ) then

  call lagitf                                                     &
  !==========
  ( rtp    , propce , tsvar  , auxl(1,1) )

endif

!===============================================================================
! 3. INTEGRATION DE LA TEMPERATURE DES PARTICULES
!===============================================================================

if ( iphyla.eq.1 .and. itpvar.eq.1 ) then

  call lagitp                                                     &
  !==========
  ( propce ,                                                      &
    tempct ,                                                      &
    tsvar  , auxl(1,1) , auxl(1,2)  )

endif

!===============================================================================
! 4. INTEGRATION DU DIAMETRE DES PARTICULES
!===============================================================================

if ( iphyla.eq.1 .and. idpvar.eq.1 ) then

  call lagidp                                                     &
  !==========
  ( tsvar  , auxl(1,1) , auxl(1,2)  )

endif

!===============================================================================
! 5. INTEGRATION DE LA MASSE DES PARTICULES
!===============================================================================

if (iphyla.eq.1 .and. impvar.eq.1) then

  call lagimp                                                     &
  !==========
  ( tsvar  , auxl(1,1) , auxl(1,2)  )

endif

!===============================================================================
! 6. INTEGRATION DES EQUATIONS DU CHARBON : HP, MCH, MCK
!===============================================================================

if (iphyla.eq.2) then

  call lagich                                                     &
  !==========
  ( propce , tempct , tsvar  , cpgd1  , cpgd2  ,                  &
    cpght  )

endif

!===============================================================================
! 7. INTEGRATION DES VARIABLES UTILISATEURS SUPPLEMENTAIRES
!===============================================================================

if (nvls.ge.1) then

  call uslaed                                                     &
  !==========
    ( ntersl , nvlsta , nvisbr ,                                  &
      dt     ,                                                    &
      taup   , tlag   , tempct , tsvar  ,                         &
      auxl(1,1) ,  auxl(1,2) ,  auxl(1,3) )

endif

!===============================================================================

!----
! FIN
!----

end subroutine
