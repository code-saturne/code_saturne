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

!-------------------------------------------------------------------------------

!> \file ppinii.f90
!> \brief Default initialization of specific modules
!> (only non-map fortran common variables of modules)
!>
!------------------------------------------------------------------------------

subroutine ppinii () &
  bind(C, name='cs_f_ppinii')

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

integer         igg, it, ir, ih, if
integer         idirac

!===============================================================================

! Mappings to C

call pp_models_init

!===============================================================================
! 1. REMPLISSAGE INCLUDE ppincl.h
!                INCLUDE GENERAL PROPRE A LA PHYSIQUE PARTICULIERE
!===============================================================================

ihm = 0 ! enthalpy, common to many models

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

! Initialization for compressible module

! Standard compressible module scalars
ienerg = 0
itempk = 0

!===============================================================================
! 2. REMPLISSAGE INCLUDE ppthch.h
!                INCLUDE THERMOCHIMIE POUR LA PHYSIQUE PARTICULIERE
!===============================================================================

! ---> Initialisation Common / TCHPPR /

! TODO : passer en c ?

do igg = 1, ngazgm
  do it = 1, npot
    ehgazg(igg,it) = zero
  enddo
  do ir = 1, nrgazm
    stoeg(igg,ir) = zero
  enddo
  ckabsg(igg)= zero
enddo

!===============================================================================
! 3. REMPLISSAGE INCLUDE coincl.h
!                INCLUDE POUR LA PHYSIQUE PARTICULIERE RELATIF A
!                LA COMBUSTION GAZ
!===============================================================================

! ---> Modele de flamme de diffusion (chimie 3 points)

nmaxh = 0
nmaxf = 0
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
hgf   = zero
tgbad = zero

!fmin = zero
!fmax = 1.d0
!hmin = zero
!hmax = zero
coeff1 = zero
coeff2 = zero
coeff3 = zero

!===============================================================================
! 6. Global variables for atmospheric flows (module atincl.f90)
!===============================================================================

! Space and time reference of the run:

! ------------------------------------
! 1d radiative transfer model:
! ----------------------------

! iatra1 -->  flag for the use of the 1d atmo radiative model

iatra1 = 0

call atmo_init_imbrication()

! key id for optimal interpolation

call field_get_key_id("opt_interp_id", kopint)

! initmeteo --> use meteo profile for variables initialization
!               (0: not used; 1: used )
! NB : should eventually be set by interface and by zone (as BCs)

initmeteo = 1

! --> Initialisation for the chemistry models:

init_at_chem = 1

! --> Initialisation for the gaseous chemistry model:

nbchim = 0
nbchmz = 0
nespgi = 0
dtchemmax = 10.d0

! --> Initialisation for the aerosol chemistry model:

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
