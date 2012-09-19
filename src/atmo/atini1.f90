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
use atsoil

!===============================================================================

implicit none

! Local variables

integer          ii, isc, jj, ipp

!===============================================================================


!===============================================================================
! 0. VERIFICATION ISCALT, ISCSTH pour IPPMOD(IATMOS) = 1 or 2
!===============================================================================
!     L'utilisateur ne doit pas y avoir touche.

if (ippmod(iatmos).ge.1) then
  if(iscalt.ne.-1) then
    write(nfecra,1000)iscalt
    call csexit (1)
    !==========
  endif

  do ii = 1, nscapp
    if(iscsth(iscapp(ii)).ne.-10) then
    write(nfecra,1001)ii,iscapp(ii),iscsth(iscapp(ii))
     call csexit (1)
     !==========
    endif
  enddo
endif

if (ippmod(iatmos).ge.2) then
  if (itytur.ne.2) then
    write(nfecra, 1002)
    call csexit(1)
  endif
endif

if (ippmod(iatmos).le.1) then
  if (iatra1.eq.1.or.iatsoil.eq.1) then
    write(nfecra, 1003)
    call csexit(1)
  endif
endif

!===============================================================================
! 1. INFORMATIONS GENERALES
!===============================================================================

!--> constants used in the atmospheric physics module
!    (see definition in atincl.h):

ps = 1.0d5
rvsra = 1.608d0
cpvcpa = 1.866d0
clatev = 2.501d6
gammat = -6.5d-03
rvap = rvsra*rair

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
  scamin(itempp) = 0.d0
  scamax(itempp) = +grand

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

  ipp = ipprtp(isca(itempp))  ! NB : already defined by the user interface ?
  nomvar(ipp)  = 'PotTemp'
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

! --> Donnees physiques ou numeriques propres aux scalaires

  do isc = 1, nscapp

    jj = iscapp(isc)

    if (iscavr(jj).le.0) then
      visls0(jj) = viscl0
    endif

    blencv(isca(jj)) = 1.d0

  enddo

  ipp = ipprtp(isca(itempl))
  nomvar(ipp)  = 'LqPotTmp'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

  ipp = ipprtp(isca(itotwt))
  nomvar(ipp)  = 'TotWater'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

  ipp = ipprtp(isca(intdrp))
  nomvar(ipp)  = 'TotDrop'
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
  nomvar(ipp)   = 'RealTemp'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

endif

if (ippmod(iatmos).eq.2) then

  ipp = ipppro(ipproc(iliqwt))
  nomvar(ipp)   = 'LiqWater'
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
! 5. Turbulent Schmidt and Prandtl number for atmospheric flows
!===============================================================================

if (nscal.gt.0) then
  do ii = 1, nscal
    sigmas(ii) = 0.7d0
  enddo
endif

!===============================================================================
! 6. Other variables for atmosĥeric flows
!===============================================================================

! Space and time reference of the run:
! ------------------------------------

!  syear --> starting hour
!  squant --> starting quantile
!  shour --> starting hour (UTC)
!  smin --> starting minute
!  ssec --> starting second

! xlon --> longitude of the domain origin
! xlat --> latitude of the domain origin

syear = 1994
squant = 1
shour = 1
smin = 0
ssec = 0
xlon = 0.d0
xlat = 45.d0

! Option for the meteo profile computation
!ihpm   --> flag to compute the hydrostastic pressure by Laplace integration
!           in the meteo profiles
!       = 0 : bottom to top Laplace integration, based on P(sea level) (default)
!       = 1 : top to bottom Laplace integration based on P computed for
!            the standard atmosphere at z(nbmaxt)

ihpm = 0

! 1d radiative transfer model:
! ----------------------------

! iatra1 -->  flag for the use of the 1d atmo radiative model
! nfatr1 --> 1d radiative model pass frequency
! ivert  --> flag for the definition of the vertical grid
!            -1: no vertical grid
!            0 : automatic definition !!!!!!!!!MM 2do
!            1 : definition by the user in usatdv
! iqv0   --> flag for humidity above the domain (0 : no humidity; 1 : decreasing)

iatra1 = 0
nfatr1 = 1
ivert = 1  ! if iatra1=1 alors ivert=1
iqv0 = 0

! flag to use the soil model (if humid atmosphere)
iatsoil = 0
! Initial values for soil variables
tsini  = 20.d0   !TSINI  : Surface ground temperature
tprini = 20.d0   !TPRINI : Deep ground temperature
qvsini = 0.d0    !QVSINI : Ground humidity
tmer   = 20.d0   !Sea temperature

!  -------------------------------------------------------------------------------
!  Microphysics parameterization options
!  -------------------------------------------------------------------------------
!  --> Option for subgrid models
!  modsub = 0 : the simplest parameterization (for numerical verifications)
!  modsub = 1 : Bechtold et al. 1995 (Luc Musson-Genon)
!  modsub = 2 : Bouzereau et al. 2004
!  modsub = 3 : Cuijpers and Duynkerke 1993, Deardorff 1976, Sommeria and Deardorff 1977
 modsub = 0

!  --> Option for liquid water content distribution models
!  moddis = 1 : all or nothing
!  moddis = 2 : Gaussian distribution
moddis = 1

!  modnuc = 0 : without nucleation
!  modnuc = 1 : Pruppacher and Klett 1997
!  modnuc = 2 : Cohard et al. 1998,1999
!  modnuc = 3 : Abdul-Razzak et al. 1998,2000 NOT IMPLEMENTED YET
modnuc = 0

! sedimentation flag
modsedi = 0

! logaritmic standard deviation of the log-normal law of the droplet spectrum
! adimensional
sigc = 0.53 ! other referenced values are 0.28, 0.15

! clipping of the humid atmosphere variables :
if ( ippmod(iatmos).eq.2 ) then
  scamin(iscapp(1)) = 200.d0
  scamin(iscapp(2)) = 0.d0
  scamin(iscapp(3)) = 0.d0
endif

!===============================================================================
! 7. ON DONNE LA MAIN A L'UTLISATEUR
!===============================================================================

!   - Interface Code_Saturne
!     ======================

if (iihmpr.eq.1) then

  call uiati1 (imeteo, ficmet, len(ficmet))
  !==========

endif

! initmeteo --> use meteo profile for variables initialization
!               (0: not used; 1: used )
! NB : should eventually be set by interface

initmeteo = 1


call usati1
!==========

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

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
'@    pour le scalaire ',I10 ,'.                              ',/,&
'@    on a :                                                  ',/,&
'@    ISCSTH(',I10   ,') = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usipsu (cs_user_parameters.f90)                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@                                                            ',/,&
'@  Seul le modele de turbulence k-eps est disponible avec    ',/,&
'@   le module atmosphere humide (ippmod(iatmos) = 2).        ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usipsu (cs_user_parameters.f90)                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1003 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@                                                            ',/,&
'@  Les modeles de sol (iatsoil) et de rayonnement (iatra1)   ',/,&
'@   ne sont disponilbes qu''avec le module atmosphere        ',/,&
'@   humide (ippomod(iatmos) = 2).                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usipsu (cs_user_parameters.f90)                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                ATMOSPHERIC  MODULE                         ',/,&
'@                                                            ',/,&
'@  ISCALT IS SPECIFIED AUTOMATICALLY.                        ',/,&
'@  iscalt should not be specified in usipsu, here:           ',/,&
'@       ISCALT  = ', I10                                      ,/,&
'@  Computation CAN NOT run.                                  ',/,&
'@                                                            ',/,&
'@  Check the input data given through the User Interface     ',/,&
'@   or in cs_user_parameters.f90.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                ATMOSPHERIC  MODULE                         ',/,&
'@                                                            ',/,&
'@  ISCSTH IS SPECIFIED AUTOMATICALLY.                        ',/,&
'@  For the scalar ', I10 ,' iscsth  should not be specified  ',/,&
'@   in usipsu, here:                                         ',/,&
'@          ISCSTH(',I10   ,') = ',I10                         ,/,&
'@  Computation CAN NOT run.                                  ',/,&
'@                                                            ',/,&
'@  Check the input data given through the User Interface     ',/,&
'@   or in cs_user_parameters.f90.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1002 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                ATMOSPHERIC  MODULE                         ',/,&
'@                                                            ',/,&
'@  Only k-eps turbulence model is available with humid       ',/,&
'@   atmosphere module (ippmod(iatmos) = 2).                  ',/,&
'@  Computation CAN NOT run.                                  ',/,&
'@                                                            ',/,&
'@  Check the input data given through the User Interface     ',/,&
'@   or in cs_user_parameters.f90.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1003 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                ATMOSPHERIC  MODULE                         ',/,&
'@                                                            ',/,&
'@  Ground model (iatsoil) and radiative model (iatra1)       ',/,&
'@   are only available with humid atmosphere module          ',/,&
'@   (ippmod(iatmos) = 2).                                    ',/,&
'@  Computation CAN NOT run.                                  ',/,&
'@                                                            ',/,&
'@  Check the input data given through the User Interface     ',/,&
'@   or in cs_user_parameters.f90.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! End
!----

return
end subroutine atini1
