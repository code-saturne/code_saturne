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

subroutine raydak &
!================

 ( ncel   , ncelet ,                                              &
   ck     , pco2   , ph2o   , fv     , temp   )

!===============================================================================
! FONCTION :
! --------

!   SOUS-PROGRAMME DU MODULE RAYONNEMENT :
!   --------------------------------------

!  CALCUL DES PROPRIETES RADIATIVES D'UN GAZ EN FONCTION DE LA
!  TEMPERATURE,DE LA COMPOSITION DES PRODUITS EN CO2, H2O ET SUIES
!  EN UTILISANT LES REGRESSIONS LINEAIRES ETABLIES PAR MODAK.

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncel             ! i  ! <-- ! number of cells                                !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ck (ncelet)      ! tr ! --> ! coefficient d'absorption du milieu             !
!                  !    !     ! (nul si transparent)                           !
! pco2(ncelet)     ! tr ! <-- ! pression partielle de co2                      !
! pco2(ncelet)     ! tr ! <-- ! pression partielle de h2o                      !
! fv  (ncelet)     ! tr ! <-- ! fraction volumique de suies                    !
! temp(ncelet)     ! tr ! <-- ! temperature                                    !
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
use entsor
use cstnum

!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel
double precision ck(ncelet), temp(ncelet), fv(ncelet)
double precision pco2(ncelet), ph2o(ncelet)

! Local variables

integer          iel
double precision alpha, path, te, ts, sootk, tmin, tmax

!===============================================================================
!===============================================================================
! CALCULS
!===============================================================================

! --- Longueur moyenne de penetration du rayonnement

path = 15.d0

tmax = 2000.d0
tmin = 300.d0

! ATTENTION : LES TEMPERATURES UTILISEES DANS MODAK SONT MISES EN KELVIN
! =========

do iel = 1, ncel

! --- Temperature du melange gazeux

  te = temp(iel)

! --- Temperature du corps noir

  ts = temp(iel)

! --- Limitation A TMAX = 2000 K et TMIN = 300 K

  if ( temp(iel).gt.tmax ) then
    ts = tmax
    te = tmax
  endif
  if ( temp(iel).lt.tmin ) then
    ts = tmin
    te = tmin
  endif

! --- Fraction volumique de suies

  sootk = 7.d0*fv(iel)/0.95d-6

! --- Calcul de l'absorptivite du fluide

  call absorb                                                     &
  !==========
  ( ts , te , path , sootk , pco2(iel) , ph2o(iel) , alpha )

! --- Test d'erreur

  if ( (1.d0-alpha).le.epzero ) then
    write(nfecra,1000) iel, alpha, pco2(iel), ph2o(iel),          &
                       sootk, te, path, fv(iel)
    call csexit(1)
  endif

! --- Calcul du coeffcient d'absorption

  ck(iel) = - log(1.d0-alpha)/path

enddo

!========
! FORMATS
!========

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ERREUR RAYDAK : CALCUL DE ABSORPTIVITE                  ',/,&
'@    =============                                           ',/,&
'@ IEL   = ', I10                                              ,/,&
'@ ALPHA = ', G15.7                                            ,/,&
'@ PCO2  = ', G15.7                                            ,/,&
'@ PH2O  = ', G15.7                                            ,/,&
'@ SOOTK = ', G15.7                                            ,/,&
'@ TE    = ', G15.7                                            ,/,&
'@ PATH  = ', G15.7                                            ,/,&
'@ FV    = ', G15.7                                            ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----
return

end subroutine
