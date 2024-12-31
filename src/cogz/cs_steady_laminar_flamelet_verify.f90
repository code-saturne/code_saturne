!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

subroutine cs_steady_laminar_flamelet_verify &
!================

 ( iok    )

!===============================================================================
!  FONCTION  :
!  ---------

! VERIFICATION DES PARAMETRES DE CALCUL
!   COMBUSTION GAZ : FLAMME DE DIFFUSION Steady laminar flamelet
!     APRES INTERVENTION UTILISATEUR
!       (COMMONS)
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
use coincl
use ppincl
use radiat

!===============================================================================

implicit none

! Arguments

integer          iok

! Local variables

!===============================================================================

!===============================================================================
! 1. OPTIONS DU CALCUL : TABLEAUX DE ppincl.h : formats 2000
!===============================================================================

! --> Coefficient de relaxation de la masse volumique

if( srrom.lt.0d0 .or. srrom.ge.1d0) then
  WRITE(NFECRA,2000)'SRROM ', SRROM
  iok = iok + 1
endif

!===============================================================================
! 2. Physical constants
!===============================================================================

! --> Masse volumique

if( ro0.lt.0d0) then
  WRITE(NFECRA,3000)'RO0   ', RO0
  iok = iok + 1
endif

! --> Fuel and oxydant reference entalpie

if (hinfue.le.-grand) then
  write(nfecra,3000)'hinfue', hinfue
  iok = iok + 1
endif
if (hinoxy.le.-grand) then
  write(nfecra,3000)'hinoxy', hinoxy
  iok = iok + 1
endif

if (ngazfl.gt.ngazgm - 1) then
  write(nfecra,3001)'ngazfl',  ngazgm, ngazfl
  iok = iok + 1
endif

if (    flamelet_zm  .eq.-1 .or. flamelet_zvar.eq.-1  &
    .or.flamelet_xr  .eq.-1 .or. flamelet_temp.eq.-1  &
    .or.flamelet_rho .eq.-1 .or. flamelet_vis .eq.-1  &
    .or.flamelet_dt  .eq.-1 ) then
    write(nfecra,3002) 'flamelet_zm'  , flamelet_zm   ,  &
                       'flamelet_zvar', flamelet_zvar ,  &
                       'flamelet_xr'  , flamelet_xr   ,  &
                       'flamelet_temp', flamelet_temp ,  &
                       'flamelet_rho' , flamelet_rho  ,  &
                       'flamelet_vis' , flamelet_vis  ,  &
                       'flamelet_dt'  , flamelet_dt
    iok = iok + 1
endif

if (ippmod(islfm).lt.2) then
  if (flamelet_ki.eq.-1) then
    write(nfecra,3003) 'flamelet_ki', flamelet_ki
    iok = iok + 1
  endif
else
  if (flamelet_c.eq.-1.or.flamelet_omg_c.eq.-1) then
    write(nfecra,3004) 'flamelet_c'    , flamelet_c    , &
                       'flamelet_omg_c', flamelet_omg_c
    iok = iok + 1
  endif
endif

!===============================================================================
! 3. Working array of coincl.h (Soot)
!===============================================================================

if (isoot.ge.1.and.iirayo.eq.0) then
  write(nfecra,4000) isoot,iirayo
  iok = iok + 1
endif

!===============================================================================
! 4. FORMATS VERIFICATION
!===============================================================================

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,                            ' DOIT ETRE UN REEL    ',/,&
'@    SUPERIEUR OU EGAL A ZERO ET INFERIEUR STRICTEMENT A 1   ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usd3p1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE RENGEIGNEE PAR L''UTILISATEUR          ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uppmod.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE INFERIEUR A', I8                        ,/,&
'@    IL VAUT ICI ', I8                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uppmod.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@  Les indices des variables dans le tableau de flammelettes ',/,&
'@  doivent etre renseignes par l''utilisateur,               ',/,&
'@  Ils valent ici:                                           ',/,&
'@   ',A15, 4x, i8                                             ,/,&
'@   ',A15, 4x, i8                                             ,/,&
'@   ',A15, 4x, i8                                             ,/,&
'@   ',A15, 4x, i8                                             ,/,&
'@   ',A15, 4x, i8                                             ,/,&
'@   ',A15, 4x, i8                                             ,/,&
'@   ',A15, 4x, i8                                             ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uppmod.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3003 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@  L''indice du scalar dissipation rate dans le tableau de   ',/,&
'@  flammelettes doit etre renseigne par l''utilisateur       ',/,&
'@  Il vaut ici:                                              ',/,&
'@   ',A15, 4x, i8                                             ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uppmod.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3004 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@  L''indice du progress variable et de la production du     ',/,&
'@  progress variable dans le tableau de flammelettes doivent ',/,&
'@  etre renseignes par l''utilisateur,                       ',/,&
'@  Ils valent ici:                                           ',/,&
'@   ',A15, 4x, i8                                             ,/,&
'@   ',A15, 4x, i8                                             ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uppmod.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ISOOT EST POSITIONNE A ',I8,'                           ',/,&
'@    SANS MODELE DE RAYONNEMENT (iirayo = ',i8,')            ',/,&
'@                                                            ',/,&
'@  Ce calcul sans interet ne sera pas execute.               ',/,&
'@                                                            ',/,&
'@  Verifier usppmo et usray1.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


!===============================================================================
! 6. SORTIE
!===============================================================================

return
end subroutine
