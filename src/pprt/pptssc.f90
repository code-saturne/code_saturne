!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

subroutine pptssc &
!================

 ( iscal  ,                                                       &
   smbrs  , rovsdt , tslagr )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE
! ON PRECISE LES TERMES SOURCES POUR UN SCALAIRE PP
!   SUR UN PAS DE TEMPS

! ATTENTION : LE TRAITEMENT DES TERMES SOURCES EST DIFFERENT
! ---------   DE CELUI DE USTSSC.F

! ON RESOUT ROVSDT*D(VAR) = SMBRS

! ROVSDT ET SMBRS CONTIENNENT DEJA D'EVENTUELS TERMES SOURCES
!  UTILISATEUR. IL FAUT DONC LES INCREMENTER ET PAS LES
!  ECRASER

! POUR DES QUESTIONS DE STABILITE, ON NE RAJOUTE DANS ROVSDT
!  QUE DES TERMES POSITIFS. IL N'Y A PAS DE CONTRAINTE POUR
!  SMBRS

! DANS LE CAS D'UN TERME SOURCE EN CEXP + CIMP*VAR ON DOIT
! ECRIRE :
!          SMBRS  = SMBRS  + CEXP + CIMP*VAR
!          ROVSDT = ROVSDT + MAX(-CIMP,ZERO)

! ON FOURNIT ICI ROVSDT ET SMBRS (ILS CONTIENNENT RHO*VOLUME)
!    SMBRS en kg variable/s :
!     ex : pour la vitesse            kg m/s2
!          pour les temperatures      kg degres/s
!          pour les enthalpies        Joules/s
!    ROVSDT en kg /s

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iscal            ! i  ! <-- ! scalar number                                  !
! smbrs(ncelet)    ! tr ! --> ! second membre explicite                        !
! rovsdt(ncelet    ! tr ! --> ! partie diagonale implicite                     !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!(ncelet,*)        !    !     !     lagrangien                                 !
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
use entsor
use optcal
use cstphy
use cstnum
use pointe, only: itypfb
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          iscal

double precision smbrs(ncelet), rovsdt(ncelet)
double precision tslagr(ncelet,*)

! Local variables

!===============================================================================

! Soot model

if (isoot.eq.1) then
  call sootsc(iscal, smbrs, rovsdt)
endif

! ---> Flamme de premelange : Modele EBU

if (ippmod(icoebu).ge.0) then
  call ebutss(iscal, smbrs, rovsdt)
endif

! ---> Flamme de premelange : Modele BML

!      IF ( IPPMOD(ICOBML).GE.0 )
!           CALL BMLTSS

! ---> Flamme de premelange : Modele LWC

if (ippmod(icolwc).ge.0) then
  call lwctss(iscal, smbrs, rovsdt)
endif

! ---> Flamme charbon pulverise

if ( ippmod(iccoal).ge.0 ) then
  call cs_coal_scast(iscal, smbrs, rovsdt)
endif

! ---> Flamme charbon pulverise couplee Transport Lagrangien
!      des particules de charbon

if ( ippmod(icpl3c).ge.0 .and. iilagr.eq.2 ) then
  call cpltss(iscal, itypfb, smbrs, rovsdt, tslagr)
endif

! ---> Flamme fuel

if ( ippmod(icfuel).ge.0 ) then
  call cs_fuel_scast(iscal, smbrs, rovsdt)
endif

! ---> Versions electriques :
!             Effet Joule
!             Arc Electrique
!             Conduction ionique

if (ippmod(ieljou).ge.1 .or.                                      &
    ippmod(ielarc).ge.1       ) then
  call eltssc(iscal, smbrs)
endif

! ---> Version atmospherique :

if (ippmod(iatmos).ge.0) then
  call attssc(iscal,smbrs)
endif


! ---> Version aerorefrigerant :

if (ippmod(iaeros).ge.0) then
  call cttssc(iscal, smbrs, rovsdt)
endif

!----
! End
!----

return

end subroutine
