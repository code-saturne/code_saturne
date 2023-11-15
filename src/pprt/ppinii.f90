!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!-------------------------------------------------------------------------------

!> \file ppinii.f90
!> \brief Default initialization of specific modules
!> (only non-map fortran common variables of modules)
!>
!------------------------------------------------------------------------------

subroutine ppinii

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use cstphy
use ppppar
use ppthch
use coincl
use cpincl
use cs_coal_incl
use ppincl
use ppcpfu
use atincl
use atimbr
use atchem
use atsoil
use field
use sshaerosol

!===============================================================================

implicit none

! Local variables

integer         icla, icha, igg, it, ir, ih, if, izone
integer         ige, iat
integer         idirac


!===============================================================================

! Mappings to C

call pp_models_init
call thch_models_init
call co_models_init
call cp_models_init
call ppcpfu_models_init

!===============================================================================
! 1. REMPLISSAGE INCLUDE ppincl.h
!                INCLUDE GENERAL PROPRE A LA PHYSIQUE PARTICULIERE
!===============================================================================

ihm = 0 ! enthalpy, common to many models

i_comb_drift = 0

!> --- Specific condensation modelling
!>     (icondb, icondv = -1 : not activated by default)
icondb = -1
icondv = -1

!> Wall condensation modelling
icondb_model  = 0 ! Wall condensation correlation (default: COPAIN)

! ---> Initialisation pour la combustion gaz
!       Variables transportees
ifm    = 0
ifp2m  = 0
ifsqm  = 0
ipvm   = 0
iygfm  = 0
icm    = 0
icp2m  = 0
ifpcpm = 0
iyfm   = 0
iyfp2m = 0
icoyfp = 0
!       Variables algebriques ou d'etat
itemp  = 0
irecvr = 0
itotki = 0
ihrr   = 0
ixr    = 0
iomgc  = 0
do igg = 1, ngazgm
  iym(igg) = 0
  ibym(igg) = -1
enddo
ickabs = 0
it2m   = 0
it4m   = 0
it3m   = 0
do idirac = 1, ndracm
  irhol (idirac) = 0
  iteml (idirac) = 0
  ifmel (idirac) = 0
  ifmal (idirac) = 0
  iampl (idirac) = 0
  itscl (idirac) = 0
  imaml (idirac) = 0
enddo

! ---> Initialisation for soot model
inpm = 0
ifsm = 0

! ---> Initialisation pour la combustion du charbon
!       Variables transportees
do icha = 1, ncharm
  if1m(icha) = 0
  if2m(icha) = 0
enddo
if4m   = 0
if5m   = 0
if6m   = 0
if7m   = 0
if8m   = 0
if9m   = 0
ifvp2m = 0
iyco2  = 0
do icla = 1, nclcpm
  ixck(icla)   = 0
  ixch(icla)   = 0
  inp(icla)    = 0
  ih2(icla)    = 0
  ixwt(icla)   = 0
  inagecp(icla) = 0
enddo
!
!       Variables algebriques ou d'etat
do ige = 1, ngazem
  iym1(ige) = 0
enddo
immel = 0
do icla = 1, nclcpm
  ix2(icla)    = 0
  itemp2(icla) = 0
  irom2(icla)  = 0
  idiam2(icla) = 0
  igmdch(icla) = 0
  igmdv1(icla) = 0
  igmdv2(icla) = 0
  igmhet(icla) = 0
  ighco2(icla) = 0
  ighh2o(icla) = 0
  igmsec(icla) = 0
enddo
do ige = 1, ngazem
  af3(ige) = 0.d0
  af4(ige) = 0.d0
  af5(ige) = 0.d0
  af6(ige) = 0.d0
  af7(ige) = 0.d0
  af8(ige) = 0.d0
  af9(ige) = 0.d0
enddo

! ---> Initialisation pour la combustion fuel
!       Variables transportees

do icla = 1, nclcpm
  ing(icla)   = 0
  iyfol(icla) = 0
  ihlf (icla) = 0
enddo
iyco2   = 0
iyhcn   = 0
iynh3   = 0
iyno    = 0
itaire  = 0

!       Variables algebriques ou d'etat

do ige = 1, ngazem
  iym1(ige) = 0
enddo

do icla=1,nclcpm
  igmeva(icla) = 0
  igmhtf(icla) = 0
enddo

ighcn1 = 0
ighcn2 = 0
ignoth = 0

ignh31 = 0
ignh32 = 0
ifhcnd = 0
ifhcnc = 0
ifnh3d = 0
ifnh3c = 0
ifnohc = 0
ifnonh = 0
ifnoch = 0
ifnoth = 0
icnohc = 0
icnonh = 0
ifhcnr = 0
icnorb = 0
igrb   = 0

ieqnox = 1
imdnox = 0
irb = 0

! Kinetic model
! solve transport equation of CO2 mass fraction by default
! (for coal or fuel combustion)
ieqco2 = 1

! ---> Coefficient de relation de la masse volumique
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
srrom = 0.95d0 ! GUI default

! Initialization for compressible module

! Standard compressible module scalars
ienerg = 0
itempk = 0
! Compressible homogeneous two-phase flow scalars
ifracv = 0
ifracm = 0
ifrace = 0

!===============================================================================
! 2. REMPLISSAGE INCLUDE ppthch.h
!                INCLUDE THERMOCHIMIE POUR LA PHYSIQUE PARTICULIERE
!===============================================================================


! ---> Initialisation Common / TCHPPI /

npo   = 0

! ---> Initialisation Common / TCHPPR /

do it = 1, npot
  th(it) = zero
enddo

do ir = 1, nrgazm
  fs(ir) = zero
enddo

! TODO : passer en c ?
do igg = 1, ngazgm
  do it = 1, npot
    ehgazg(igg,it) = zero
  enddo
  do ir = 1, nrgazm
    stoeg(igg,ir) = zero
  enddo
  wmolg(igg) = zero
  ckabsg(igg)= zero
enddo

ckabs1 = zero

do ige = 1, ngazem
  do it = 1, npot
    ehgaze(ige,it) = zero
  enddo
  wmole(ige) = zero
enddo

do iat = 1, natom
  wmolat(iat) = zero
enddo

xco2 = zero
xh2o = zero

!===============================================================================
! 3. REMPLISSAGE INCLUDE coincl.h
!                INCLUDE POUR LA PHYSIQUE PARTICULIERE RELATIF A
!                LA COMBUSTION GAZ
!===============================================================================

! ---> Modele de flamme de diffusion (chimie 3 points)

!if (     ippmod(icod3p).ge.0 .or. ippmod(islfm).ge.0          &
! .or. ippmod(icoebu).ge.0 .or. ippmod(icolwc).ge.0) then
!  call co_models_bc_map !necessaire pour mapper ientox et ientfu en C mais il faut call cs_boundary_conditions_create en amont qui alloue ces tabs, or il est appelé après cette routine. Ainsi je commente l'init de ientox sachant qu'il est déja initialisé en C dans cs_bc_create
  !endif

nmaxh = 0
nmaxf = 0
tinoxy = zero
tinfue = zero
!do izone = 1, nozppm
!  ientox(izone) = 0
!  ientfu(izone) = 0
!enddo
hinfue = -grand
hinoxy = -grand
hstoea = -grand
do ih = 1, nmaxhm
  hh(ih) = -grand
enddo
do if = 1, nmaxfm
  ff(if)= zero
  do ih = 1, nmaxhm
    tfh(if,ih) = zero
  enddo
enddo

! ---> Modele de la flamme de diffusion Steady laminar flamelet

nki     = -1
nxr     = -1
nzm     = -1
nzvar   = -1
nlibvar = -1
ngazfl  = -1

flamelet_zm    = -1
flamelet_zvar  = -1
flamelet_ki    = -1
flamelet_xr    = -1
flamelet_temp  = -1
flamelet_rho   = -1
flamelet_vis   = -1
flamelet_dt    = -1
flamelet_temp2 = -1
flamelet_hrr   = -1

flamelet_species(:)  = -1
flamelet_species_name(:)  = ' '


flamelet_c     = -1
flamelet_omg_c = -1

! ---> Modele de flamme de premelange (modele EBU et LWC)
!      On prend 300K pour temperature des gaz frais.
!        En suite de calcul c'est ecrase par le fichier suite
!        Cette valeur ne sert que lors du premier appel a ebuphy
!          (ensuite, la valeur imposee dans les CL prend le pas)

cebu  = 2.5d0
vref  = zero
lref  = zero
ta    = zero
tstar = zero
frmel = zero
tgf   = 300.d0
!do izone = 1, nozppm
!  ientgf(izone) = 0
!  ientgb(izone) = 0
!  fment(izone)  = zero
!  tkent(izone)  = zero
!enddo
hgf   = zero
tgbad = zero

fmin = zero
fmax = 1.d0
hmin = zero
hmax = zero
coeff1 = zero
coeff2 = zero
coeff3 = zero

!===============================================================================
! 4. REMPLISSAGE INCLUDE cpincl.h
!                INCLUDE POUR LA PHYSIQUE PARTICULIERE RELATIF A
!                LA COMBUSTION CP
!===============================================================================

! ---> Donnees relatives au charbon

do icha = 1, ncharm

  cch(icha)    = zero
  hch(icha)    = zero
  och(icha)    = zero
  sch(icha)    = zero
  nch(icha)    = zero

  alpha(icha)  = zero
  beta(icha)   = zero
  teta (icha)  = zero
  omega(icha)  = zero

  pcich(icha)  = zero
  rho0ch(icha) = zero
  thcdch(icha) = zero

  cck(icha)    = zero
  hck(icha)    = zero
  ock(icha)    = zero
  sck(icha)    = zero
  nck(icha)    = zero

  rhock(icha)  = zero
  gamma(icha)  = zero
  delta(icha)  = zero
  kappa(icha)  = zero
  zeta (icha)  = zero
  pcick(icha)  = zero

  cpashc(icha) = zero
  h0ashc(icha) = zero

  iy1ch(icha)  = 0
  y1ch(icha)   = zero
  a1ch(icha)   = zero
  e1ch(icha)   = zero
  crepn1(1,icha) = zero
  crepn1(2,icha) = zero

  iy2ch(icha)  = 0
  y2ch(icha)   = zero
  a2ch(icha)   = zero
  e2ch(icha)   = zero
  crepn2(1,icha) = zero
  crepn2(2,icha) = zero

  ahetch(icha) = zero
  ehetch(icha) = zero
  iochet(icha) = 0

  ahetc2(icha) = zero
  ehetc2(icha) = zero
  ioetc2(icha) = 0

  ahetwt(icha) = zero
  ehetwt(icha) = zero
  ioetwt(icha) = 0

  repnck(icha) = zero
  repnle(icha) = zero
  repnlo(icha) = zero

  ychxle(icha) = zero
  ychxlo(icha) = zero
  yhcnle(icha) = zero
  yhcnlo(icha) = zero
  ynh3le(icha) = zero
  ynh3lo(icha) = zero
  ycoch1(icha) = zero
  yhcnc1(icha) = zero
  ynoch1(icha) = zero
  ycoch2(icha) = zero
  yhcnc2(icha) = zero
  ynoch2(icha) = zero

  nnch(icha)   = zero
  nnckle(icha) = zero
  nhckle(icha) = zero
  ncckle(icha) = zero
  nncklo(icha) = zero
  nhcklo(icha) = zero
  nccklo(icha) = zero

  wchx1c(icha) = zero
  wchx2c(icha) = zero
enddo

wmchx1 = zero
wmchx2 = zero

! ---> Definition des Pointeurs du tableau TBMCR utilise dans cpphy1.F
!      et les sous-programmes appeles

  do icha = 1, ncharm
    if1mc(icha) = 0
    if2mc(icha) = 0
  enddo
  ix1mc   = 0
  ix2mc   = 0
  ichx1f1 = 0
  ichx2f2 = 0
  icof1   = 0
  icof2   = 0
  ih2of1  = 0
  ih2of2  = 0
  ih2sf1  = 0
  ih2sf2  = 0
  ihcnf1  = 0
  ihcnf2  = 0

! ---> Donnees relatives a la combustion des especes gazeuses

do icha = 1, ncharm
  ichx1c(icha) = 0
  ichx2c(icha) = 0
enddo
ichx1 = 0
ichx2 = 0
ico   = 0
ih2s  = 0
ih2   = 0
ihcn  = 0
io2   = 0
ico2  = 0
ih2o  = 0
iso2  = 0
inh3  = 0
in2   = 0

xsi   = 3.76d0

do icha = 1, ncharm
  chx1(icha) = zero
  chx2(icha) = zero
  a1(icha)   = zero
  a2(icha)   = zero
  b1(icha)   = zero
  b2(icha)   = zero
  c1(icha)   = zero
  c2(icha)   = zero
  d1(icha)   = zero
  d2(icha)   = zero
  e1(icha)   = zero
  e2(icha)   = zero
  f1(icha)   = zero
  f2(icha)   = zero
enddo

!===============================================================================
! 5. REMPLISSAGE INCLUDE fuelincl.h
!                INCLUDE POUR LA PHYSIQUE PARTICULIERE RELATIF A
!                LA COMBUSTION FUEL
!===============================================================================

! ---> Donnees relatives a la combustion des especes gazeuses

ico   = 0
io2   = 0
ico2  = 0
ih2o  = 0
in2   = 0

xsi   = 3.76d0

!===============================================================================
! 6. Global variables for atmospheric flows (module atincl.f90)
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
! iqv0   --> flag for humidity above the domain (0 : no humidity; 1 : decreasing)

iatra1 = 0
nfatr1 = 1
iqv0 = 0

!  -------------------------------------------------------------------------------
!  Microphysics parameterization options
!  -------------------------------------------------------------------------------

! logaritmic standard deviation of the log-normal law of the droplet spectrum
! adimensional
sigc = 0.53 ! other referenced values are 0.28, 0.15

!  -----------------------------------------------------------------------------
!  Atmospheric imbrication on large scale meteo (atimbr module)
!  -----------------------------------------------------------------------------

! activation flag
imbrication_flag = .false.
imbrication_verbose = .false.

! ------------------------------------------------------------------------------
! flags for activating the cressman interpolation for the boundary conditions
! ------------------------------------------------------------------------------

cressman_u = .false.
cressman_v = .false.
cressman_tke = .false.
cressman_eps = .false.
cressman_theta = .false.
cressman_qw = .false.
cressman_nc = .false.

! --------------------------------------------------------------
! numerical parameters for the cressman interpolation formulas
! --------------------------------------------------------------
horizontal_influence_radius = 8500.d0
vertical_influence_radius = 100.d0

! key id for optimal interpolation

call field_get_key_id("opt_interp_id", kopint)

! initmeteo --> use meteo profile for variables initialization
!               (0: not used; 1: used )
! NB : should eventually be set by interface and by zone (as BCs)

initmeteo = 1

! --> Initialisation for the chemistry models:

init_at_chem = 1

! --> Initialisation for the gaseous chemistry model:

ifilechemistry = 0
nbchim = 0
nbchmz = 0
nespgi = 0
dtchemmax = 10.d0

! --> Initialisation for the aerosol chemistry model:

do izone = 1, nozppm
  iprofa(izone) = 0
enddo

! Default values (climatic ones) for radiative transfer and
! aerosols
aod_o3_tot=0.20d0
aod_h2o_tot=0.10d0
gaero_o3=0.66d0
gaero_h2o=0.64d0
piaero_o3=0.84d0
piaero_h2o=0.84d0
black_carbon_frac=0.d0
zaero = 6000d0

return
end subroutine ppinii
