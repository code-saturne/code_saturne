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
use spefun
use field_operator

!===============================================================================

implicit none

! Local variables

integer          ivart, iel

double precision xvart, rhum, rscp, pp, zent
double precision lrhum, lrscp
double precision qsl, deltaq
double precision yw_liq, qwt, tliq, dum

double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: cvar_vart, cvar_totwt
double precision, dimension(:), pointer :: cpro_tempc, cpro_liqwt

logical activate

!===============================================================================
! 0. INITIALISATIONS A CONSERVER
!===============================================================================

activate = .false.

! Initialize variables to avoid compiler warnings

ivart = -1

! --- Initialisation memoire

! This routine computes the density and the thermodynamic temperature.
! The computations require the pressure profile which is here taken from
! the meteo file. If no meteo file is used, the user should
! give the laws for RHO and T in cs_user_physical_properties.f90

if (imeteo.eq.0) return

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

call field_get_val_s(itempc,cpro_tempc)
call field_get_val_s(ivarfl(ivart), cvar_vart)
if (ippmod(iatmos).ge.2) then
  call field_get_val_s(ivarfl(isca(iymw)), cvar_totwt)
  call field_get_val_s(iliqwt,cpro_liqwt)
endif

! From potential temperature, compute:
! - Temperature in Celsius
! - Density
! ----------------------

! Computes the perfect gaz constants according to the physics

rhum = rair
rscp = rair/cp0

lrscp = rair/cp0

do iel = 1, ncel

  zent = xyzcen(3,iel)

  ! Reference pressure
  if (imeteo.eq.0) then !FIXME useless...
    call atmstd(zent,pp,dum,dum)
  else
    ! Pressure profile from meteo file:
    call intprf &
    !==========
       ( nbmett, nbmetm,                                            &
         ztmet , tmmet , phmet , zent, ttcabs, pp )
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

  call cs_rho_humidair(qwt, tliq, pp, yw_liq, cpro_tempc(iel), crom(iel))

  ! Humid atmosphere
  if (ippmod(iatmos).ge.2) then
    cpro_liqwt(iel) = yw_liq
  endif

enddo

! Gaussian subgrid condensation scheme for humid atmosphere physics
if (ippmod(iatmos).ge.2.and. moddis.ge.2) then
  call gaussian()
endif

! User re-definition:
call usatph ()

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
integer    iccocg
integer    inc

double precision a_const
double precision a_coeff
double precision alpha,al
double precision sig_flu ! standard deviation of qw'-alpha*theta'
double precision var_tl,var_q,cov_tlq
double precision q1,qsup

double precision, dimension(:), pointer :: cvar_k, cvar_ep

! rvap = rair*rvsra

allocate(dtlsd(3,ncelet))
allocate(dqsd(3,ncelet))

iccocg = 1
inc = 1

! computation of grad(theta_l)
call field_gradient_scalar(ivarfl(isca(iscalt)), 1, imrgra, inc,    &
                           iccocg,                                  &
                           dtlsd)

! computation of grad(qw)
call field_gradient_scalar(ivarfl(isca(iymw)), 1, imrgra, inc,      &
                           iccocg,                                  &
                           dqsd)


call field_get_val_s(ivarfl(ik), cvar_k)
call field_get_val_s(ivarfl(iep), cvar_ep)

! -------------------------------------------------------------
! Gradients are used for estimating standard deviations of the
! subgrid fluctuations.
! -------------------------------------------------------------

lrhum = rhum

a_const = 2.d0*cmu/2.3d0
do iel = 1, ncel

  a_coeff = a_const*cvar_k(iel)**3/cvar_ep(iel)**2 ! 2 cmu/c2 * k**3 / eps**2
  var_tl= a_coeff*(dtlsd(1,iel)**2 + dtlsd(2,iel)**2 + dtlsd(3,iel)**2)
  var_q = a_coeff*( dqsd(1,iel)**2 + dqsd(2,iel)**2 + dqsd(3,iel)**2)
  cov_tlq = a_coeff*(  dtlsd(1,iel)*dqsd(1,iel)   &
                     + dtlsd(2,iel)*dqsd(2,iel)   &
                     + dtlsd(3,iel)*dqsd(3,iel))

  zent = xyzcen(3,iel)

  if (imeteo.eq.0) then
    call atmstd(zent,pp,dum,dum)
  else
    ! Pressure profile from meteo file:
    call intprf(nbmett, nbmetm, ztmet , tmmet , phmet , zent, ttcabs, pp)
  endif

  xvart = cvar_vart(iel) ! thermal scalar: liquid potential temperature
  tliq = xvart*(pp/ps)**rscp ! liquid temperature
  qwt  = cvar_totwt(iel) ! total water content
  qsl = cs_air_yw_sat(tliq-tkelvi, pp) ! saturated vapor content
  deltaq = qwt - qsl
  alpha = (clatev*qsl/(rvap*tliq**2))*(pp/ps)**rscp
  sig_flu = sqrt(var_q + alpha**2*var_tl - 2.d0*alpha*cov_tlq) !FIXME a^2 +b^2 -2ab = (a-b)^2, a = a_coeff dqsd; b = alpha dtlsd

  if (sig_flu.lt.1.d-30) sig_flu = 1.d-30
  q1 = deltaq/sig_flu
  al = 1.d0/(1.d0 + qsl*clatev**2/(rair*rvsra*cp0*tliq**2))
  qsup = qsl/sig_flu

  nebdia(iel) = 0.5d0*(1.d0 + ferf(q1/sqrt(2.d0)))

  !FIXME MF : put in input of the global function...
  yw_liq = (sig_flu                                                               &
        /(1.d0 + qsl*clatev**2/(rvap*cp0*tliq**2)))                             &
        *(nebdia(iel)*q1 + exp(-q1**2/2.d0)/sqrt(2.d0*pi))
  yw_liq = max(yw_liq, 1d-15)
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

  cpro_liqwt(iel) = yw_liq
  ! Celcius temperature of the air parcel
  cpro_tempc(iel) = tliq + (clatev/cp0)*yw_liq - tkelvi

  !FIXME back to the previous formulation
  call cs_rho_humidair(qwt, tliq, pp,   &
                       cpro_liqwt(iel), &
                       cpro_tempc(iel), &
                       crom(iel))

enddo

! when properly finished deallocate dtlsd
deallocate(dtlsd)
deallocate(dqsd)

end subroutine gaussian

end subroutine atphyv
