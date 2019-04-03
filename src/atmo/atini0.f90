!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
!> \file atini0.f90
!> \brief Initialisation of variable options for Code_Saturne atmospheric
!>     module in addition to what is done previously in iniini function
subroutine atini0

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
use atimbr
use siream
use field

!===============================================================================

implicit none

! Local variables

!===============================================================================

!===============================================================================
! 0. General info
!===============================================================================

!--> constants used in the atmospheric physics module
!    (see definition in atincl.f90):

ps = 1.0d5
cpvcpa = 1.866d0
gammat = -6.5d-03
rvap = rvsra*rair

!===============================================================================
! 1. Global variables for atmospheric flows
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
!  modsub = 3 : Cuijpers and Duynkerke 1993, Deardorff 1976, Sommeria
!               and Deardorff 1977
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

iatmst = 0

theo_interp = 0

! -------------------------------------
! flags for 1d radiative transfer model
! -------------------------------------

! no computation / storage of downward and upward infrared radiative fluxes
irdu = 0
! no computation / storage of downward and upward solar radiative fluxes
soldu = 0

!===============================================================================
! 2. ON DONNE LA MAIN A L'UTLISATEUR
!===============================================================================

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

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine atini0
