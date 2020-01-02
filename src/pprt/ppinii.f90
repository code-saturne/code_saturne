!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

subroutine ppinii
!================


!===============================================================================
!  FONCTION  :
!  ---------

! INITIALISATION PAR DEFAUT DES COMMONS PHYSIQUE PARTICULIERE

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
use cstnum
use cstphy
use ppppar
use ppthch
use coincl
use cpincl
use cs_coal_incl
use cs_fuel_incl
use ppincl
use ppcpfu
use atincl
use atchem
use siream

!===============================================================================

implicit none

! Local variables

integer         icla, icha, igg, it, ir, ih, if, izone
integer         isol, ige, iat , iii
integer         idirac


!===============================================================================

! Mappings to C

call pp_models_init
call thch_models_init
call co_models_init
call cp_models_init
call fuel_models_init

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

! ---> Initialisation pour la combustion gaz
!       Variables transportees
ifm    = 0
ifp2m  = 0
iygfm  = 0
icm    = 0
icp2m  = 0
ifpcpm = 0
iyfm   = 0
iyfp2m = 0
icoyfp = 0
!       Variables algebriques ou d'etat
itemp  = 0
do igg = 1, ngazgm
  iym(igg) = 0
  ibym(igg) = -1
enddo
ickabs = 0
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
do icha = 1, ncharb
  if1m(icha) = 0
  if2m(icha) = 0
enddo
if4m   = 0
if5m   = 0
if6m   = 0
if7m   = 0
if8m   = 0
if9m   = 0
if4p2m = 0
ifvp2m = 0
iyco2  = 0
iage = 0
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
itemp1 = 0
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
ifvap   = 0
if4p2m  = 0
iyco2   = 0
iyhcn   = 0
iynh3   = 0
iyno    = 0
itaire  = 0

!       Variables algebriques ou d'etat

itemp1 = 0
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
srrom = 1.d0

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
ngaze = 0
ngazg = 0
nato  = 0
nrgaz = 0

! ---> Initialisation Common / TCHPPR /

do it = 1, npot
  th(it) = zero
enddo

do ir = 1, nrgazm
  fs(ir) = zero
enddo

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

diftl0 = zero

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

nmaxh = 0
nmaxf = 0
tinoxy = zero
tinfue = zero
do izone = 1, nozppm
  ientox(izone) = 0
  ientfu(izone) = 0
enddo
hinfue = zero
hinoxy = zero
hstoea = zero
do ih = 1, nmaxhm
  hh(ih) = zero
enddo
do if = 1, nmaxfm
  ff(if)= zero
  do ih = 1, nmaxhm
    tfh(if,ih) = zero
  enddo
enddo

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
do izone = 1, nozppm
  ientgf(izone) = 0
  ientgb(izone) = 0
  qimp(izone)   = zero
  fment(izone)  = zero
  tkent(izone)  = zero
enddo
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

! ----   Par charbon
ncharb = 0
nsolid = 0

do icha = 1, ncharm

  nclpch(icha) = 0

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

  xashch(icha) = zero
  cpashc(icha) = zero
  h0ashc(icha) = zero

  xwatch(icha) = zero

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

  ich(icha)    = 0
  ick(icha)    = 0
  iash(icha)   = 0
  iwat(icha)   = 0

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

do isol = 1, nsolim
  do it = 1, npot
    ehsoli(isol,it) = zero
  enddo
  wmols(isol) = zero
enddo

! ----   Par classe

nclacp = 0
do icla = 1, nclcpm
  ichcor(icla) = 0
  diam20(icla) = zero
  dia2mn(icla) = zero
  rho20(icla)  = zero
  rho2mn(icla) = zero
  xmp0(icla)   = zero
  xmash(icla)  = zero
enddo

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


! ---> Donnees complementaires relatives au calcul de rho sur les facettes
!      d'entree

do izone = 1, nozppm
  ientat(izone) = 0
  ientcp(izone) = 0
  timpat(izone) = zero
  do icla = 1, nclcpm
    x20(izone,icla) = zero
  enddo
enddo

!===============================================================================
! 5. REMPLISSAGE INCLUDE fuelincl.h
!                INCLUDE POUR LA PHYSIQUE PARTICULIERE RELATIF A
!                LA COMBUSTION FUEL
!===============================================================================

! ---> Donnees relatives au fuel

nsolid = 0

! ----   Par classe

nclafu = 0

!    1 seul fioul

cfol    = zero
hfol    = zero
ofol     = zero
pcifol  = zero
rho0fl = zero

ckf    = zero
hkf    = zero
okf    = zero
rhokf  = zero
pcikf = zero

yfol   = zero
afol   = zero
efol   = zero

ahetfl = zero
ehetfl = zero

iofhet = 0

ifol   = 0
ikf   = 0

do isol = 1, nsolim
  do it = 1, npot
    ehsoli(isol,it) = zero
  enddo
  wmols(isol) = zero
enddo

! ---> Donnees relatives a la combustion des especes gazeuses

ifov  = 0
ico   = 0
io2   = 0
ico2  = 0
ih2o  = 0
in2   = 0

xsi   = 3.76d0
fvapmx =  2.d0*0.012d0 / (2.d0*0.028d0 + xsi*0.028d0)

fov  = zero
a  = zero
b  = zero


! ---> Donnees complementaires relatives au calcul de rho sur les facettes
!      d'entree

do izone = 1, nozppm
  ientat(izone)  = 0
  ientfl(izone) = 0
  timpat(izone)  = zero
  timpfl(izone) = zero
enddo

!===============================================================================
! 6. REMPLISSAGE INCLUDE atincl.h
!                INCLUDE POUR LA VERSION ATMOSPHERIQUE
!===============================================================================

!--> Initialisation for the meteo profile
imeteo = 0
nbmetd = 0
nbmett = 0
nbmetm = 0
nbmaxt = 0
do izone = 1, nozppm
  iprofm(izone) = 0
enddo
ihpm=0

! --> Initialisation for the 1d radiative model:

iatra1 = 0
nfatr1 = 1
nvert = 1
kvert = 20
ivert = 0
iqv0 = 0
iatsoil = 0

! --> Initialisation for the chemistry models:

init_at_chem = 1

! --> Initialisation for the gaseous chemistry model:
ichemistry = 0
ifilechemistry = 0
isepchemistry = 2
iphotolysis = 1
nbchim = 0
nbchmz = 0
nespgi = 0
dtchemmax = 10.d0
do izone = 1, nozppm
  iprofc(izone) = 0
enddo


! --> Initialisation for the aerosol chemistry model:
iaerosol = 0
inogaseouschemistry = 0
nespg_siream = 0
! Number of iterations for time splitting
ncycle_aer = 1
! Flag to activate or not coagulation (1=yes, 0=no)
icoag_siream = 1
! Flag to activate or not condensation/evaporation (1=yes, 0=no)
icond_siream = 1
! Flag to activate or not nucleation (1=yes, 0=no)
inucl_siream = 1
! Flag to consider or not kelvin effect (1=yes, 0=no)
ikelv_siream = 1
! Cutting bin between equilibrium (1 to icut_siream)
! and dynamic bins (icut_siream to nbin_aer)
icut_siream = nbin_aer
! Sulfate condensation computation method (1=dynamic, 0=equilibrium)
isulfcond_siream = 1
! Solver for dynamic condensation (1='etr', 2='ros2', 3='ebi')
kdslv_siream = 2
! Redistribution method of lagrangian bins
! (1=number conserving, 2=interpolation)
iredist_siream = 1
! Nucleation model (0=binary, 1=ternary)
itern_siream = 0
! Method used for estimation of wet diameter (0=isoropia, 1=gerber)
ithrm_siream = 1
! m
bin_bound_aer = (/0.01d-6, 0.0398d-6, 0.1585d-6, 0.6310d-6, 2.5119d-6, 10.d-6/)
fixed_density_aer = 1400.0d0 ! kg/m3
do iii = 1, nbin_aer
  density_aer(iii) = fixed_density_aer
enddo
do izone = 1, nozppm
  iprofa(izone) = 0
enddo

return
end subroutine ppinii
