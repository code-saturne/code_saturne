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
!> \file atphyv.f90
!> \brief Functions that compute physical variables for each cell
!>  for the atmospheric module
!
!> \brief Initialise physical variables of the atmospheric module \n
!> Remarques :
!> This routine is called at the beginning of each time step
!>
!-------------------------------------------------------------------------------

subroutine atphyv

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum, only: pi
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use atincl
use field
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

integer          ivart, iel

double precision xvart, rhum, rscp, pp, zent
double precision lrhum, theta0
double precision qsl, deltaq
double precision yw_liq, qwt, tliq, dum

double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: cvar_vart, cvar_totwt
double precision, dimension(:), pointer :: cpro_tempc, cpro_liqwt
double precision, dimension(:), pointer :: cpro_beta
double precision, dimension(:), pointer :: cpro_met_p, cpro_met_rho

logical activate

!===============================================================================
! 0. INITIALISATIONS A CONSERVER
!===============================================================================

activate = .false.

! Initialize variables to avoid compiler warnings

ivart = -1

if (idilat.eq.0) then
  call field_get_val_s_by_name("thermal_expansion", cpro_beta)
endif

if (imeteo.ge.2) then
  call field_get_val_s_by_name('meteo_pressure', cpro_met_p)
  call field_get_val_s_by_name('meteo_density', cpro_met_rho)
endif

! This routine computes the density and the thermodynamic temperature.
! The computations may require the pressure profile which is here taken from
! the meteo file. If no meteo file is used, the user can
! give the laws for RHO and T in cs_user_physical_properties.f90

!===============================================================================

!   Positions des variables, coefficients
!   -------------------------------------

! --- Numero de variable thermique
!       (et de ses conditions limites)
!       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
!          IVART = ISCA(2)

if (iscalt.gt.0) then
  ivart = isca(iscalt)
else
  write(nfecra,9010) iscalt
  call csexit (1)
endif

! Density
call field_get_val_s(icrom, crom)

call field_get_val_s(itempc, cpro_tempc)
call field_get_val_s(ivarfl(ivart), cvar_vart)
if (ippmod(iatmos).ge.2) then
  call field_get_val_s(ivarfl(isca(iymw)), cvar_totwt)
  call field_get_val_s(iliqwt,cpro_liqwt)
endif

! From potential temperature, compute:
! - Temperature in Celsius
! - Density
! ----------------------

! Computes the perfect gas constants according to the physics

rhum = rair
rscp = rair/cp0
! Adiabatic (constant) potential temperature
theta0 = t0 * (p0/ps)**rscp
do iel = 1, ncel

  zent = xyzcen(3,iel)

  ! Reference pressure
  if (imeteo.eq.0) then
    call atmstd(zent,pp,dum,dum)
  else if (imeteo.eq.1) then
    ! Pressure profile from meteo file:
    call intprf &
       ( nbmett, nbmetm,                                            &
         ztmet , tmmet , phmet , zent, ttcabs, pp )
  else
    pp = cpro_met_p(iel)
  endif

  ! Potential temperature
  ! or liquid potential temperature for humid atmosphere
  xvart = cvar_vart(iel)

  ! (liquid) temperature
  ! law: T = theta * (p/ps) ** (Rair/Cp0)
  tliq = xvart*(pp/ps)**rscp

  ! Dry atmosphere, total water fraction
  qwt = 0.d0

  ! Humid atmosphere
  if (ippmod(iatmos).ge.2) then
    qwt = cvar_totwt(iel)
  endif

  ! Density in cell centers:
  ! ------------------------
  ! law: rho = P / ( R_mixture * T_mixture(K) )

  if (idilat.eq.0) then
    ! Boussinesq with respect to the adiabatic density
    if (imeteo.ge.2) then
      crom(iel) = cpro_met_rho(iel)
      ! "delta rho = - beta rho0 delta theta" gives
      ! "beta = 1 / theta"
      cpro_beta(iel) = 1.d0 / theta0
      ! Compute T in Celisus
      cpro_tempc(iel) = tliq - tkelvi
    ! Boussinesq with respect to rho0
    else
      crom(iel) = ro0
      ! "delta rho = - beta rho0 delta theta" gives
      ! "beta = 1 / theta"
      cpro_beta(iel) = 1.d0 / xvart
      ! Compute T in Celisus
      cpro_tempc(iel) = tliq - tkelvi
    endif
  else

    call cs_rho_humidair(qwt, tliq, pp, yw_liq, cpro_tempc(iel), crom(iel))
  endif

  ! Humid atmosphere
  if (ippmod(iatmos).ge.2) then
    cpro_liqwt(iel) = yw_liq
  endif

enddo

! Gaussian subgrid condensation scheme for humid atmosphere physics
if (ippmod(iatmos).ge.2.and. moddis.ge.2) then
  call gaussian()
endif

!===============================================================================
! FORMATS
!----

9010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    APPEL A csexit DANS LE SOUS PROGRAMME atphyv            ',/,&
'@                                                            ',/,&
'@    La variable dont dependent les proprietes physiques ne  ',/,&
'@      semble pas etre une variable de calcul.               ',/,&
'@    En effet, on cherche a utiliser la temperature alors que',/,&
'@      ISCALT = ',I10                                         ,/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de cs_user_physical_properties,',      /,&
'@     (et le test lors de la definition de IVART).'           ,/,&
'@    Verifier la definition des variables de calcul dans     ',/,&
'@      usipsu. Si un scalaire doit jouer le role de la       ',/,&
'@      temperature, verifier que ISCALT a ete renseigne.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

return
contains

! *******************************************************************
!> \brief Internal function -
!>   subgrid condensation scheme assuming a gaussian distribution for the
!> fluctuations of both qw and thetal.
!-------------------------------------------------------------------------------
subroutine gaussian()
double precision, dimension(:,:), allocatable :: dtlsd
double precision, dimension(:,:), allocatable :: dqsd
integer    inc

double precision a_const
double precision a_coeff
double precision alpha,al
double precision sig_flu ! standard deviation of qw'-alpha*theta'
double precision var_q_tl
double precision q1,qsup, rvap, rscp
double precision ek, ep

double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_nusa, cvar_omg
double precision, dimension(:,:), pointer :: cvar_rij
double precision, dimension(:), pointer :: nn, nebdia

rvap = rair*rvsra
rscp = rair/cp0

allocate(dtlsd(3,ncelet))
allocate(dqsd(3,ncelet))

inc = 1

! computation of grad(theta_l)
call field_gradient_scalar(ivarfl(isca(iscalt)), 1, inc, dtlsd)

! computation of grad(qw)
call field_gradient_scalar(ivarfl(isca(iymw)), 1, inc, dqsd)

if (itytur.eq.2) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)
elseif (itytur.eq.3) then
  call field_get_val_v(ivarfl(irij), cvar_rij)
  call field_get_val_s(ivarfl(iep), cvar_ep)
elseif (iturb.eq.50) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)
elseif (iturb.eq.60) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iomg), cvar_omg)
elseif (iturb.eq.70) then
  call field_get_val_s(ivarfl(inusa), cvar_nusa)
endif


call field_get_val_s_by_name("nebulosity_frac", nn)
call field_get_val_s_by_name("nebulosity_diag", nebdia)

! -------------------------------------------------------------
! Gradients are used for estimating standard deviations of the
! subgrid fluctuations.
! -------------------------------------------------------------

lrhum = rhum

a_const = 2.d0*cmu/2.3d0
do iel = 1, ncel
  if (itytur.eq.2) then
    a_coeff = a_const*cvar_k(iel)**3/cvar_ep(iel)**2 ! 2 cmu/c2 * k**3 / eps**2
  elseif (itytur.eq.3) then
    ek = 0.5d0*(cvar_rij(1,iel) + cvar_rij(2,iel) + cvar_rij(3,iel))
    a_coeff = a_const*ek**3/cvar_ep(iel)**2 ! 2 cmu/c2 * k**3 / eps**2
    call field_get_val_s(ivarfl(iep), cvar_ep)
  elseif (iturb.eq.50) then
    a_coeff = a_const*cvar_k(iel)**3/cvar_ep(iel)**2 ! 2 cmu/c2 * k**3 / eps**2
  elseif (iturb.eq.60) then
    ep = cvar_omg(iel)*cvar_k(iel)*cmu
    a_coeff = a_const*cvar_k(iel)**3/cvar_ep(iel)**2 ! 2 cmu/c2 * k**3 / eps**2
  elseif (iturb.eq.70) then
    ! using cvar_nusa(iel) = cmu*xkent**2/xeent
    ! FIXME: There is no good way to calculate tke and eps from nusa.
    ! For the moment we use tke**4/eps**2 instead of tke**3/eps**2
    ! Need to return WARNING that in case of Spalart-Allmaras we use bad assumpltion
    ! or RETURN error for this case.
    a_coeff = a_const*cvar_nusa(iel)**2/cmu**2 ! 2 cmu/c2 * k**3 / eps**2
  endif

  zent = xyzcen(3,iel)

  if (imeteo.eq.0) then
    call atmstd(zent,pp,dum,dum)
  else if (imeteo.eq.1) then
    ! Pressure profile from meteo file:
    call intprf(nbmett, nbmetm, ztmet , tmmet , phmet , zent, ttcabs, pp)
  else
    pp = cpro_met_p(iel)
  endif

  xvart = cvar_vart(iel) ! thermal scalar: liquid potential temperature
  tliq = xvart*(pp/ps)**rscp ! liquid temperature
  qsl = cs_air_yw_sat(tliq-tkelvi, pp) ! saturated vapor content
  alpha = (clatev*qsl/(rvap*tliq**2))*(pp/ps)**rscp

  var_q_tl = a_coeff * ( (dqsd(1,iel) - alpha * dtlsd(1,iel))**2  &
                       + (dqsd(2,iel) - alpha * dtlsd(2,iel))**2  &
                       + (dqsd(3,iel) - alpha * dtlsd(3,iel))**2)

  sig_flu = max(sqrt(var_q_tl), 1.d-30)

  qwt  = cvar_totwt(iel) ! total water content
  deltaq = qwt - qsl
  q1 = deltaq/sig_flu
  al = 1.d0/(1.d0 + qsl*clatev**2/(rair*rvsra*cp0*tliq**2))
  qsup = qsl/sig_flu

  nebdia(iel) = 0.5d0*(1.d0 + erf(q1/dsqrt(2.d0)))

  !FIXME MF : put in input of the global function...
  yw_liq = (sig_flu                                                               &
        /(1.d0 + qsl*clatev**2/(rvap*cp0*tliq**2)))                             &
        *(nebdia(iel)*q1 + exp(-q1**2/2.d0)/sqrt(2.d0*pi))
  yw_liq = max(yw_liq, 0.d0)
  nn(iel) = nebdia(iel) - (nebdia(iel)*q1                                       &
          + exp(-q1**2/2.d0)/sqrt(2.d0*pi))*exp(-q1**2/2.d0)/sqrt(2.d0*pi)

  ! go back to all or nothing
  if(qwt.lt.yw_liq)then

    nn(iel) = 0.d0

    ! deltaq set to 0 if unsaturted air parcel
    if (deltaq.le.0.d0) then ! unsaturated air parcel
      deltaq = 0.d0
      nebdia(iel) = 0.d0
    else ! saturated (ie. with liquid water) air parcel
      nebdia(iel) = 1.d0
    endif

    ! TODO input ?
    ! 0 if unsaturated air parcel
    yw_liq = deltaq / (1.d0 + qsl*clatev**2/(rvap*cp0*tliq**2))
  endif ! qwt.lt.yw_liq

  ! Celcius temperature of the air parcel
  cpro_tempc(iel) = tliq + (clatev/cp0)*yw_liq - tkelvi
  ! liquid water content
  cpro_liqwt(iel) = yw_liq
  !density
  lrhum = rair*(1.d0 - yw_liq + (rvsra - 1.d0)*(qwt - yw_liq))
  crom(iel) = pp/(lrhum*(tliq + (clatev/cp0)*yw_liq))

enddo

! when properly finished deallocate dtlsd
deallocate(dtlsd)
deallocate(dqsd)

end subroutine gaussian

end subroutine atphyv
