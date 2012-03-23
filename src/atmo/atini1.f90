!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine atini1
!================


!===============================================================================
!  FONCTION  :
!  ---------

!   INIT DES OPTIONS DES VARIABLES POUR LA VERSION ATMOSPHERIQUE
!      EN COMPLEMENT DE CE QUI A DEJA ETE FAIT DANS USIPSU

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
use ihmpre
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use ppincl

use atincl

!===============================================================================

implicit none

! Local variables

integer          ii, isc, jj, ipp

!===============================================================================


!===============================================================================
! 0. VERIFICATION ISCALT, ISCSTH pour IPPMOD(IATMOS) = 1 or 2
!===============================================================================
!     L'utilisateur ne doit pas y avoir touche.

if ( ippmod(iatmos).ge.1 ) then
  if(iscalt.ne.-1) then
    write(nfecra,1000)iscalt
    call csexit (1)
    !==========
  endif

  do ii = 1, nscapp
    if(iscsth(iscapp(ii)).ne.-10) then
    write(nfecra,1001)ii,iscapp(ii),iscapp(ii),iscsth(iscapp(ii))
     call csexit (1)
     !==========
    endif
  enddo
endif


!===============================================================================
! 1. INFORMATIONS GENERALES
!===============================================================================

! ---> Masse volumique et viscosite
irovar = 0
ivivar = 0

!===============================================================================
! 2. VARIABLES TRANSPORTEES pour IPPMOD(IATMOS) = 1 or 2
!===============================================================================

! 2.1  Dry atmosphere
! =====================

if ( ippmod(iatmos).eq.1 ) then

  iscsth(itempp) = 1
  iscalt = itempp
  scamin(itempp)   = 0.d0
  scamax(itempp)   = +grand

!  for the dry atmosphere case, non constant density
  irovar = 1

! --> Donnees physiques ou numeriques propres aux scalaires

  do isc = 1, nscapp

    jj = iscapp(isc)

    if (iscavr(jj).le.0) then
      visls0(jj) = viscl0
    endif

    blencv(isca(jj)) = 1.d0

  enddo

  ipp = ipprtp(isca(itempp))
  nomvar(IPP)  = 'PotTemp'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

endif

! 2.2  Humid atmosphere
! =====================

if ( ippmod(iatmos).eq.2 ) then

  iscsth(itempl) = 1
  iscsth(itotwt) = 0
  iscsth(intdrp) = 0

  iscalt = itempl

!  for the humid atmosphere case, non constant density
  irovar = 1


  ipp = ipprtp(isca(itempl))
  nomvar(IPP)  = 'LqPotTmp'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

  ipp = ipprtp(isca(itotwt))
  nomvar(IPP)  = 'TotWater'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

  ipp = ipprtp(isca(intdrp))
  nomvar(IPP)  = 'TotDrop'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

endif

!===============================================================================
! 3. VARIABLES D'ETAT pour IPPMOD(IATMOS) = 1 or 2
!===============================================================================

! 3.1  Dry or humid atmosphere
! =============================

if (ippmod(iatmos).eq.1 .or. ippmod(iatmos).eq.2) then

  ipp = ipppro(ipproc(itempc))
  nomvar(IPP)   = 'TempC'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

endif


!===============================================================================
! 4. One scale turbulent model for k-eps closure for IPPMOD(IATMOS) = 1 or 2
!===============================================================================

if (ippmod(iatmos).eq.1 .or. ippmod(iatmos).eq.2) then

  if (itytur.eq.2) then
    ideuch = 0
  endif

endif

!===============================================================================
! 5. Turbulent Schmidt number for atmospheric flows
!===============================================================================

if (nscal.gt.0) then
  do ii = 1, nscal
    sigmas(ii) = 0.7d0
  enddo
endif

!===============================================================================
! 6. ON DONNE LA MAIN A L'UTLISATEUR
!===============================================================================

!   - Interface Code_Saturne
!     ======================

if (iihmpr.eq.1) then

  call uiati1 (imeteo)
  !==========

endif

call usati1
!==========

!--------
! FORMATS
!--------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@                                                            ',/,&
'@  La valeur de ISCALT est renseignee automatiquement.       ',/,&
'@                                                            ',/,&
'@  L''utilisateur ne doit pas la renseigner dans usipsu, or  ',/,&
'@    elle a ete affectee comme suit :                        ',/,&
'@    ISCALT = ',I10                                           ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usipsu (dans cs_user_parameters.f90)             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@                                                            ',/,&
'@  Les valeurs de ISCSTH sont renseignees automatiquement.   ',/,&
'@                                                            ',/,&
'@  L''utilisateur ne doit pas les renseigner dans usipsu, or ',/,&
'@    pour le scalaire ',I10   ,' correspondant au scalaire   ',/,&
'@    physique particuliere ',I10   ,' on a                   ',/,&
'@    ISCSTH(',I10   ,') = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usipsu (cs_user_parameters.f90)                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

return
end subroutine
