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
use atchem
use siream

!===============================================================================

implicit none

! Local variables

integer          ii, isc, jj, ipp
character*2                :: nbin
!===============================================================================


!===============================================================================
! 0. VERIFICATIONS
!===============================================================================

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
! ===================

if (ippmod(iatmos).eq.1) then

  itpscl = 1

  iscacp(iscalt) = 1
  scamin(iscalt) = 0.d0
  scamax(iscalt) = +grand

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

  ipp = ipprtp(isca(iscalt))  ! NB : already defined by the user interface ?
  nomvar(ipp)  = 'PotTemp'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

endif

! 2.2  Humid atmosphere
! =====================

if (ippmod(iatmos).eq.2) then

  iscacp(iscalt) = 1
  iscacp(itotwt) = 0
  iscacp(intdrp) = 0

  itpscl = 1

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

  ipp = ipprtp(isca(iscalt))
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

! Atmospheric gaseous chemistry
! Do not change this order
if (iaerosol.eq.1) ichemistry = 3
! if a chemical scheme is solved, a concentration profiles
! file must be used
if (ichemistry.ge.1) ifilechemistry = ichemistry
if (inogaseouschemistry.eq.1) ichemistry = 0

if (ifilechemistry.ge.1) then

  ! logical unit and name of the chemical profiles file
  impmec = 27
  ficmec = 'chemistry'

  ! Initialization of the chemical scheme
  ! quasi stationary equilibrium NOx scheme with 4 species and 5 reactions
  if (ifilechemistry.eq.1) then
    nrg = 5
    nespg = 4
    allocate(dmmk(nespg))
    allocate(chempoint(nespg))

    nomvar(ipprtp(isca(1)))='NO'
    nomvar(ipprtp(isca(2)))='NO2'
    nomvar(ipprtp(isca(3)))='O3'
    nomvar(ipprtp(isca(4)))='O3P'
    dmmk(1)=30.d-3 ! Molar mass NO
    dmmk(2)=46.d-3 ! Molar mass NO2
    dmmk(3)=48.d-3 ! Molar mass O3
    dmmk(4)=16.d-3 ! Molar mass O3P
    chempoint = (/ 4, 3, 2, 1 /)
  ! scheme with 20 species and 34 reactions
  else if (ifilechemistry.eq.2) then
    nrg = 34
    nespg = 20
    allocate(dmmk(nespg))
    allocate(chempoint(nespg))

    nomvar(ipprtp(isca(1)))='NO'
    nomvar(ipprtp(isca(2)))='NO2'
    nomvar(ipprtp(isca(3)))='O3'
    nomvar(ipprtp(isca(4)))='O3P'
    nomvar(ipprtp(isca(5)))='O1D'
    nomvar(ipprtp(isca(6)))='OH'
    nomvar(ipprtp(isca(7)))='HO2'
    nomvar(ipprtp(isca(8)))='H2O2'
    nomvar(ipprtp(isca(9)))='NO3'
    nomvar(ipprtp(isca(10)))='N2O5'
    nomvar(ipprtp(isca(11)))='HONO'
    nomvar(ipprtp(isca(12)))='HNO3'
    nomvar(ipprtp(isca(13)))='CO'
    nomvar(ipprtp(isca(14)))='HCHO'
    nomvar(ipprtp(isca(15)))='ALD2'
    nomvar(ipprtp(isca(16)))='C2O3'
    nomvar(ipprtp(isca(17)))='PAN'
    nomvar(ipprtp(isca(18)))='XO2'
    nomvar(ipprtp(isca(19)))='SO2'
    nomvar(ipprtp(isca(20)))='H2SO4'
    dmmk(1)=30.d-3     ! Molar mass NO
    dmmk(2)=46.d-3     ! Molar mass NO2
    dmmk(3)=48.d-3     ! Molar mass O3
    dmmk(4)=16.d-3     ! Molar mass O3P
    dmmk(5)=16.d-3     ! Molar mass O1D
    dmmk(6)=17.01d-3   ! Molar mass OH
    dmmk(7)=33.01d-3   ! Molar mass HO2
    dmmk(8)=34.01d-3   ! Molar mass H2O2
    dmmk(9)=62.01d-3   ! Molar mass NO3
    dmmk(10)=108.01d-3 ! Molar mass N2O5
    dmmk(11)=47.01d-3  ! Molar mass HONO
    dmmk(12)=63.01d-3  ! Molar mass HNO3
    dmmk(13)=28.01d-3  ! Molar mass CO
    dmmk(14)=30.03d-3  ! Molar mass HCHO
    dmmk(15)=44.05d-3  ! Molar mass ALD2
    dmmk(16)=75.04d-3  ! Molar mass C2O3
    dmmk(17)=121.05d-3 ! Molar mass PAN
    dmmk(18)=47.03d-3  ! Molar mass XO2
    dmmk(19)=64.06d-3  ! Molar mass SO2
    dmmk(20)=98.08d-3  ! Molar mass H2SO4
    chempoint = (/ 20, 19, 16, 17, 2, 15, 14, 3, 18, 7, 8, 9, 4, &
                   10, 1, 12, 11, 13, 5, 6 /)
  ! scheme CB05 with 52 species and 155 reactions
  else if (ifilechemistry.eq.3) then
    if (iaerosol.eq.1) then
      nrg = 162
      nespg = 65
    else
      nrg = 155
      nespg = 52
    endif
    allocate(dmmk(nespg))
    allocate(chempoint(nespg))

    nomvar(ipprtp(isca(1)))='NO'
    nomvar(ipprtp(isca(2)))='NO2'
    nomvar(ipprtp(isca(3)))='O3'
    nomvar(ipprtp(isca(4)))='O3P'
    nomvar(ipprtp(isca(5)))='O1D'
    nomvar(ipprtp(isca(6)))='OH'
    nomvar(ipprtp(isca(7)))='HO2'
    nomvar(ipprtp(isca(8)))='H2O2'
    nomvar(ipprtp(isca(9)))='NO3'
    nomvar(ipprtp(isca(10)))='N2O5'
    nomvar(ipprtp(isca(11)))='HONO'
    nomvar(ipprtp(isca(12)))='HNO3'
    nomvar(ipprtp(isca(13)))='HNO4'
    nomvar(ipprtp(isca(14)))='CO'
    nomvar(ipprtp(isca(15)))='HCHO'
    nomvar(ipprtp(isca(16)))='ALD2'
    nomvar(ipprtp(isca(17)))='C2O3'
    nomvar(ipprtp(isca(18)))='PAN'
    nomvar(ipprtp(isca(19)))='ALDX'
    nomvar(ipprtp(isca(20)))='CXO3'
    nomvar(ipprtp(isca(21)))='PANX'
    nomvar(ipprtp(isca(22)))='XO2'
    nomvar(ipprtp(isca(23)))='XO2N'
    nomvar(ipprtp(isca(24)))='NTR'
    nomvar(ipprtp(isca(25)))='ETOH'
    nomvar(ipprtp(isca(26)))='CH4'
    nomvar(ipprtp(isca(27)))='MEO2'
    nomvar(ipprtp(isca(28)))='MEOH'
    nomvar(ipprtp(isca(29)))='MEPX'
    nomvar(ipprtp(isca(30)))='FACD'
    nomvar(ipprtp(isca(31)))='ETHA'
    nomvar(ipprtp(isca(32)))='ROOH'
    nomvar(ipprtp(isca(33)))='AACD'
    nomvar(ipprtp(isca(34)))='PACD'
    nomvar(ipprtp(isca(35)))='PAR'
    nomvar(ipprtp(isca(36)))='ROR'
    nomvar(ipprtp(isca(37)))='ETH'
    nomvar(ipprtp(isca(38)))='OLE'
    nomvar(ipprtp(isca(39)))='IOLE'
    nomvar(ipprtp(isca(40)))='ISOP'
    nomvar(ipprtp(isca(41)))='ISPD'
    nomvar(ipprtp(isca(42)))='TERP'
    nomvar(ipprtp(isca(43)))='TOL'
    nomvar(ipprtp(isca(44)))='XYL'
    nomvar(ipprtp(isca(45)))='CRES'
    nomvar(ipprtp(isca(46)))='TO2'
    nomvar(ipprtp(isca(47)))='OPEN'
    nomvar(ipprtp(isca(48)))='CRO'
    nomvar(ipprtp(isca(49)))='MGLY'
    nomvar(ipprtp(isca(50)))='SO2'
    nomvar(ipprtp(isca(51)))='H2SO4'
    nomvar(ipprtp(isca(52)))='HCO3'
    dmmk(1)=30.d-3     ! Molar mass NO
    dmmk(2)=46.d-3     ! Molar mass NO2
    dmmk(3)=48.d-3     ! Molar mass O3
    dmmk(4)=16.d-3     ! Molar mass O3P
    dmmk(5)=16.d-3     ! Molar mass O1D
    dmmk(6)=17.01d-3   ! Molar mass OH
    dmmk(7)=33.01d-3   ! Molar mass HO2
    dmmk(8)=34.01d-3   ! Molar mass H2O2
    dmmk(9)=62.01d-3   ! Molar mass NO3
    dmmk(10)=108.01d-3 ! Molar mass N2O5
    dmmk(11)=47.01d-3  ! Molar mass HONO
    dmmk(12)=63.01d-3  ! Molar mass HNO3
    dmmk(13)=79.01d-3  ! Molar mass HNO4
    dmmk(14)=28.01d-3  ! Molar mass CO
    dmmk(15)=30.03d-3  ! Molar mass HCHO
    dmmk(16)=44.05d-3  ! Molar mass ALD2
    dmmk(17)=75.04d-3  ! Molar mass C2O3
    dmmk(18)=121.05d-3 ! Molar mass PAN
    dmmk(19)=43.04d-3  ! Molar mass ALDX
    dmmk(20)=74.04d-3  ! Molar mass CXO3
    dmmk(21)=120.04d-3 ! Molar mass PANX
    dmmk(22)=47.03d-3  ! Molar mass XO2
    dmmk(23)=47.03d-3  ! Molar mass XO2N
    dmmk(24)=77.04d-3  ! Molar mass NTR
    dmmk(25)=46.07d-3  ! Molar mass ETOH
    dmmk(26)=16.04d-3  ! Molar mass CH4
    dmmk(27)=47.03d-3  ! Molar mass MEO2
    dmmk(28)=32.04d-3  ! Molar mass MEOH
    dmmk(29)=48.04d-3  ! Molar mass MEPX
    dmmk(30)=46.03d-3  ! Molar mass FACD
    dmmk(31)=30.07d-3  ! Molar mass ETHA
    dmmk(32)=47.03d-3  ! Molar mass ROOH
    dmmk(33)=60.05d-3  ! Molar mass AACD
    dmmk(34)=76.05d-3  ! Molar mass PACD
    dmmk(35)=15.03d-3  ! Molar mass PAR
    dmmk(36)=16.d-3    ! Molar mass ROR
    dmmk(37)=28.05d-3  ! Molar mass ETH
    dmmk(38)=27.05d-3  ! Molar mass OLE
    dmmk(39)=56.11d-3  ! Molar mass IOLE
    dmmk(40)=68.18d-3  ! Molar mass ISOP
    dmmk(41)=70.09d-3  ! Molar mass ISPD
    dmmk(42)=136.24d-3 ! Molar mass TERP
    dmmk(43)=92.14d-3  ! Molar mass TOL
    dmmk(44)=106.16d-3 ! Molar mass XYL
    dmmk(45)=108.14d-3 ! Molar mass CRES
    dmmk(46)=141.15d-3 ! Molar mass TO2
    dmmk(47)=48.04d-3  ! Molar mass OPEN
    dmmk(48)=108.14d-3 ! Molar mass CRO
    dmmk(49)=72.06d-3  ! Molar mass MGLY
    dmmk(50)=64.06d-3  ! Molar mass SO2
    dmmk(51)=98.08d-3  ! Molar mass H2SO4
    dmmk(52)=63.03d-3  ! Molar mass HCO3
    if (iaerosol.eq.1) then
      nomvar(ipprtp(isca(53)))='HC8'
      nomvar(ipprtp(isca(54)))='API'
      nomvar(ipprtp(isca(55)))='LIM'
      nomvar(ipprtp(isca(56)))='CVARO1'
      nomvar(ipprtp(isca(57)))='CVARO2'
      nomvar(ipprtp(isca(58)))='CVALK1'
      nomvar(ipprtp(isca(59)))='CVOLE1'
      nomvar(ipprtp(isca(60)))='CVAPI1'
      nomvar(ipprtp(isca(61)))='CVAPI2'
      nomvar(ipprtp(isca(62)))='CVLIM1'
      nomvar(ipprtp(isca(63)))='CVLIM2'
      nomvar(ipprtp(isca(64)))='NH3'
      nomvar(ipprtp(isca(65)))='HCL'
      nomvar(ipprtp(isca(66)))='CVBIBMP'
      nomvar(ipprtp(isca(67)))='CVANCLP'
      nomvar(ipprtp(isca(68)))='CVBIISO1'
      nomvar(ipprtp(isca(69)))='CVBIISO2'
      dmmk(53)=114.0d0   ! Molar mass HC8
      dmmk(54)=136.0d0   ! Molar mass API
      dmmk(55)=136.0d0   ! Molar mass LIM
      dmmk(56)=150.0d0   ! Molar mass CVARO1
      dmmk(57)=150.0d0   ! Molar mass CVARO2
      dmmk(58)=140.0d0   ! Molar mass CVALK1
      dmmk(59)=140.0d0   ! Molar mass CVOLE1
      dmmk(60)=184.0d0   ! Molar mass CVAPI1
      dmmk(61)=184.0d0   ! Molar mass CVAPI2
      dmmk(62)=200.0d0   ! Molar mass CVLIM1
      dmmk(63)=200.0d0   ! Molar mass CVLIM2
      dmmk(64)=17.0d0    ! Molar mass NH3
      dmmk(65)=36.5d0    ! Molar mass HCL
    endif
    if (iaerosol.eq.1) then
      chempoint = (/ 64, 65, 59, 57, 3, 55, 61, 20, 56, 16, 23, 50,&
                     17, 51, 54, 60, 62, 13, 48, 58, 18, 63, 52, 46,&
                     4, 5, 53, 14, 22, 35, 6, 32, 49, 33, 47, 19, 34,&
                     36, 37, 44, 45, 39, 7, 8, 40, 15, 41, 43, 42, 9,&
                     10, 21, 11, 24, 27, 30, 31, 12, 38, 25, 26, 28, 29,&
                     1, 2 /)
    else
      chempoint = (/ 48, 52, 47, 43, 1, 42, 50, 17, 44, 9, 15, 38, 13, 37,&
                     41, 45, 51, 10, 35, 46, 14, 49, 39, 33, 2, 3, 40, 11,&
                     19, 20, 4, 21, 36, 22, 34, 16, 23, 24, 25, 31, 32, 26,&
                     5, 6, 27, 12, 28, 30, 29, 7, 8, 18 /)
    endif
 endif

endif

! Atmospheric aerosol chemistry
if (iaerosol.eq.1) then

  write (nfecra,*) "Atmospheric aerosol chemistry activated"

  ! logical unit and name of the chemical profiles file
  impmea=28
  ficmea='aerosols'

  ! 4 species not followed in gaseous scheme 3 but necessary to siream
  nespg_siream=nespg+4


  ! Verification
  if (ifilechemistry.ne.3) then
    write(nfecra,1004)
    call csexit (1)
  endif

  ! The user must declare at least as many user scalars as the number
  ! of chemical species nespg_siream+nbin_aer*nesp_aer
  ! The user can declare more user scalars than nespg. But only the first
  ! nespg will be considered for chemistry resolution
  if (nscaus.lt.nespg_siream+nesp_aer*nbin_aer+nbin_aer) then
    write(nfecra,1005) nscaus, nespg_siream+nesp_aer*nbin_aer+nbin_aer
    call csexit (1)
  endif

  ! Filling the names array
  esp_siream=(/'MD    ', 'BC    ', 'NA    ', 'H2SO4 ', 'NH3   ', 'HNO3  ',&
 'HCL   ','ARO1  ', 'ARO2  ', 'ALK1  ', 'OLE1  ', 'API1  ', 'API2  ',&
 'LIM1  ', 'LIM2  ', 'ANCLP ', 'BIISO1', 'BIISO2',&
 'BIBMP ', 'POA   ', 'H2O   '/)
  do ii = 1, nbin_aer
    write(nbin,"(i2)") ii
    do jj = 1, nesp_aer
      nomvar(ipprtp(isca(nespg_siream+ii+(jj-1)*nbin_aer)))=&
             trim(esp_siream(jj))//'_bin'//trim(adjustl(nbin))
    enddo
      nomvar(ipprtp(isca(nespg_siream+nesp_aer*nbin_aer+ii)))=&
            'Naero_bin'//trim(adjustl(nbin))
  enddo
endif

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

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

 1004 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    CHIMIE ATMOSPHERIQUE PARTICULAIRE DEMANDEE              ',/,&
'@                                                            ',/,&
'@  Lorsque le modele de chimie particulaire siream est utilise',/,&
'@   un schema gazeux complet(CB05) est automatiquement utilise',/,&
'@  L''utilisateur ne doit pas en specifier d''autre (ichemistry)',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usati1 (cs_user_parameters.f90)                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1005 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    CHIMIE ATMOSPHERIQUE PARTICULAIRE DEMANDEE              ',/,&
'@                                                            ',/,&
'@  Le nombre de scalaires utilisateurs declares doit etre    ',/,&
'@  superieur ou egal au nombre d''especes chimiques         ',/,&
'@                                                            ',/,&
'@   Nombre de scalaires utilisateurs declares : ',I10         ,/,&
'@   Nombre d''especes chimiques declarees : ',I10             ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
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
'@  ISCACP IS SPECIFIED AUTOMATICALLY.                        ',/,&
'@  For the scalar ', I10 ,' iscacp  should not be specified  ',/,&
'@   in usipsu, here:                                         ',/,&
'@          ISCACP(',I10   ,') = ',I10                         ,/,&
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

 1004 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ATMOSPHERIC AEROSOL CHEMISTRY MODULE                    ',/,&
'@                                                            ',/,&
'@  When aerosol chemistry model siream is used               ',/,&
'@   a full gaseous scheme (CB05) is automatically used       ',/,&
'@  The user cannot specify any other scheme (ichemistry)     ',/,&
'@  Computation CAN NOT run.                                  ',/,&
'@                                                            ',/,&
'@  Check the input data given through the User Interface     ',/,&
'@   or in cs_user_parameters.f90.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1005 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
     '@ @@  WARNING:   STOP WHILE READING INPUT DATA          ',/,&
'@    =========                                               ',/,&
'@    ATMOSPHERIC AEROSOL CHEMISTRY MODULE                    ',/,&
'@                                                            ',/,&
'@  The number of user scalars must be greater                ',/,&
'@  than the number of chemical species                       ',/,&
'@                                                            ',/,&
'@   Number of user scalars declared: ',I10                    ,/,&
'@   Number of chemical species declared : ',I10               ,/,&
'@                                                            ',/,&
'@  The computation will not be run                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! End
!----

return
end subroutine atini1
