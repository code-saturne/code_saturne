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

subroutine pptssc &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , rtp    , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iscal            ! i  ! <-- ! scalar number                                  !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
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
use pointe, only: itypfb, izfppp
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision smbrs(ncelet), rovsdt(ncelet)
double precision tslagr(ncelet,*)

! Local variables


!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================



!===============================================================================
! 2. AIGUILLAGE VERS LE MODELE ADEQUAT
!===============================================================================

! Soot model

if (isoot.eq.1) then
  call sootsc                                                     &
  !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp ,                                                       &
   dt     , rtpa   , rtp    , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt )
endif

! ---> Flamme de premelange : Modele EBU

if ( ippmod(icoebu).ge.0 ) then
  call ebutss                                                     &
  !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp ,                                                       &
   dt     , rtpa   , rtp    , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt )
endif


! ---> Flamme de premelange : Modele BML

!      IF ( IPPMOD(ICOBML).GE.0 )
!           CALL BMLTSS

! ---> Flamme de premelange : Modele LWC

if ( ippmod(icolwc).ge.0 ) then
  call lwctss                                                     &
  !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp ,                                                       &
   dt     , rtpa   , rtp    , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt )
endif

! ---> Flamme charbon pulverise

if ( ippmod(iccoal).ge.0 ) then
  call cs_coal_scast                                              &
   !================
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   itypfb ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp ,                                                       &
   dt     , rtpa   , rtp    , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt )
endif

! ---> Flamme charbon pulverise couplee Transport Lagrangien
!      des particules de charbon

if ( ippmod(icpl3c).ge.0 .and. iilagr.eq.2 ) then
  call cpltss                                                     &
   !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   itypfb ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp ,                                                       &
   dt     , rtpa   , rtp    , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt , tslagr )
endif

! ---> Flamme fuel

if ( ippmod(icfuel).ge.0 ) then
  call cs_fuel_scast                                              &
   !================
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   itypfb ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp ,                                                       &
   dt     , rtpa   , rtp    , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt )
endif

! ---> Versions electriques :
!             Effet Joule
!             Arc Electrique
!             Conduction ionique

if ( ippmod(ieljou).ge.1 .or.                                     &
     ippmod(ielarc).ge.1 .or.                                     &
     ippmod(ielion).ge.1       ) then
   call eltssc                                                    &
   !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   itypfb ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp ,                                                       &
   dt     , rtpa   , rtp    , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt )
endif

! ---> Version atmospherique :

if ( ippmod(iatmos).ge.0 ) then
   call attssc                                                    &
   !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   itypfb ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp ,                                                       &
   dt     , rtpa   , rtp    , propce , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt )
endif


! ---> Version aerorefrigerant :

if ( ippmod(iaeros).ge.0 ) then
   call cttssc                                                    &
   !==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   itypfb ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   izfppp ,                                                       &
   dt     , rtpa   , rtp    ,                                     &
   coefa  , coefb  , ckupdc , smacel ,                            &
   smbrs  , rovsdt )
endif

!----
! FIN
!----

return

end subroutine
