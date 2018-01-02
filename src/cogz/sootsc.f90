!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!===============================================================================
! Function:
! ---------

!> \file sootsc.f90
!>
!> \brief Specific physic subroutine: two equations soot model.
!>
!> This subroutine defines the source terms for the soot mass fraction
!> and the precursor number for soot model of Moss et al for one time step.
!
!  The equations read: \f$ rovsdt \delta a = smbrs \f$
!
!  \f$ rovsdt \f$ et \f$ smbrs \f$ could already contain source term
!  and don't have to be erased but incremented.
!
!  For stability sake, only positive terms should be add in \f$ rovsdt \f$.
!  There is no constrain for \f$ smbrs \f$.
!
!  For a source term written \f$ S_{exp} + S_{imp} a \f$, source terms are:
!           \f$ smbrs  = smbrs  + S_{exp} + S_{imp} a \f$
!           \f$ rovsdt = rovsdt + \max(-S_{imp},0) \f$
!
!  Here are set \f$ rovsdt \f$ and \f$ smbrs \f$ containning \f$ \rho \Omega \f$
!   - \f$ smbrs \f$ in \f$ kg_a.s^{-1} \f$ (ex: for velocity:
!     \f$ kg.m.s^{-2} \f$, for temperature: \f$ kg.°C.s^{-1} \f$,
!     for enthalpy: \f$ J.s^{-1} \f$)
!   - \f$ rovsdt \f$ en \f$ kg.s^{-1} \f$
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iscal         scalar index
!> \param[in,out] smbrs         explicit right hand side
!> \param[in,out] rovsdt        implicit terms
!_______________________________________________________________________________

subroutine sootsc &
 ( iscal  ,                                                       &
   smbrs  , rovsdt )

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
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          iscal

double precision smbrs(ncelet), rovsdt(ncelet)

! Local variables

character(len=80) :: chaine
integer          ivar, iel

double precision epsi
parameter       (epsi = 1.d-6)
double precision d1s3, d2s3, cexp, cimp
double precision zetan, zetas, rho, xfu, xm, temp, nn0
double precision ka, kb, kz, kt, chi, po2, wox
double precision aa, bb, cc, taa, tcc, caa, cbb, ccc, dd
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: cvar_scal, cvara_scal
double precision, dimension(:), pointer :: cvara_fsm, cvara_npm
double precision, dimension(:), pointer :: cpro_temp
double precision, dimension(:), pointer :: cpro_ym1, cpro_ym2, cpro_ym3

type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

ivar = isca(iscal)
call field_get_label(ivarfl(ivar), chaine)
call field_get_val_s(icrom, crom)

if (ivar.eq.isca(ifsm).or.ivar.eq.isca(inpm)) then
  call field_get_val_s(ivarfl(isca(iscal)), cvar_scal)
  call field_get_val_s(itemp, cpro_temp)
  call field_get_val_s(iym(1), cpro_ym1)
  call field_get_val_s(iym(2), cpro_ym2)
  call field_get_val_s(iym(3), cpro_ym3)
  call field_get_val_prev_s(ivarfl(isca(iscal)), cvara_scal)
  call field_get_val_prev_s(ivarfl(isca(ifsm)), cvara_fsm)
  call field_get_val_prev_s(ivarfl(isca(inpm)), cvara_npm)
endif

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

!===============================================================================
! 2. Writings
!===============================================================================

if (vcopt%iwarni.ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

!=======================================================================
! --- Moss et al.:
! zeta_s (isca(ifsm)) soot mass fraction zeta_s = (rho_s/rho).f_v
! zeta_n (isca(inpm)) precursor density  zeta_n = n / (rho.No)
!=======================================================================

if (ivar.eq.isca(ifsm).or.ivar.eq.isca(inpm)) then

  ! To be changed for other combustible !FIXME
  ! Methane CH4 (Syed, Stewart and Moss Symposium 1990)
  caa = 6.54d4 !m^3/kg^2.K^0.5.s
  cbb = 1.3d7 ! m^3.K^-1/2.s^-1
  ccc = 0.1d0 ! m^3.kg^-2/3.K^-1/2.s^-1
  taa = 46.1d3 ! K
  tcc = 12.6d3 ! K

  d1s3 = 1.d0/3.d0
  d2s3 = 2.d0/3.d0

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(cpro_temp(1))
    call synsca(cvar_scal)
  endif

  do iel = 1, ncel

    cexp = 0.d0
    cimp = 0.d0

    nn0 = 6.0223d23
    rho = crom(iel)                  ! Mixture density (kg/m3)
    temp = cpro_temp(iel) ! Temperature

    xm = 1.d0/ (  cpro_ym1(iel)/wmolg(1)             &
                + cpro_ym2(iel)/wmolg(2)             &
                + cpro_ym3(iel)/wmolg(3) )

    xfu = cpro_ym1(iel) * xm / wmolg(1) ! Fuel molar fraction

    ! --- rate of particule nucleation
    aa = caa * rho**2 * temp**0.5d0 * xfu * exp(-taa/temp)

    ! --- coagulation
    bb = cbb * temp**0.5d0

    ! --- surface growth of soot
    cc = ccc * rho * temp**0.5d0 * xfu * exp(-tcc/temp)

    po2 = cpro_ym2(iel)*xm/wmolg(2)*1.d0/4.76d0

    ! --- oxidation
    ka = 20.d0*exp(-15098.d0/temp)
    kb = 4.46d-3*exp(-7650.d0/temp)
    kt = 1.51d5*exp(-48817.d0/temp)
    kz = 21.3d0*exp(2063.d0/temp)

    chi = kb*po2/(kb*po2+kt)

    wox = 1.2d2*( (ka*po2*chi)/(1.d0+kz*po2) + kb*po2*(1.d0-chi) )

    dd = (36.d0*acos(-1.d0)/rosoot**2.d0)**d1s3

   ! -------------------------------------------------------------

    zetas = cvara_fsm(iel) ! fraction massique de suies (SU)
    zetan = cvara_npm(iel) ! densite de precurseurs (SU)

    if (ivar.eq.isca(ifsm)) then

      ! --- Surface growth : quadratic
      if (zetas.gt.epsi) cimp = volume(iel) *                     &
        (  nn0**d1s3 * rho * cc * zetas**(-d1s3) * zetan**d1s3    &
         - rho * dd *nn0**d1s3 *zetan**d1s3 *zetas**(-d1s3)*wox )
      cexp = volume(iel) * (  144.d0*aa )
    endif

    if (ivar.eq.isca(inpm)) then
      cimp = volume(iel) * ( - rho**2.d0 * bb * zetan )
      cexp = volume(iel) * ( aa )
    endif

    smbrs(iel)  = smbrs(iel)  + cexp + cimp*cvara_scal(iel)
    rovsdt(iel) = rovsdt(iel) + max(-cimp,0.d0)

  enddo

endif

!--------
! Formats
!--------

 1000 format(' TERMES SOURCES PHYSIQUE PARTICULIERE POUR LA VARIABLE '  &
       ,a8,/)

!----
! End
!----

return

end subroutine
