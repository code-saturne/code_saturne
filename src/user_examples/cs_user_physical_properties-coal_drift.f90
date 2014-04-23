!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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
! Purpose:
! -------

!> \file cs_user_physical_properties-coal_drift.f90
!> \brief Definition of physical variable laws for scalars with a drift.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     mbrom         indicator of filling of romb array
!> \param[in]     dt            time step (per cell)
!> \param[in]     rtp, rtpa     calculated variables at cell centers
!> \param[in]                    (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!_______________________________________________________________________________

subroutine usphyv &
 ( nvar   , nscal  ,                                              &
   mbrom  ,                                                       &
   dt     , rtp    , rtpa   ,                                     &
   propce )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use field
use ppincl
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          mbrom

double precision dt(ncelet), rtp(ncelet,nflown:nvar), rtpa(ncelet,nflown:nvar)
double precision propce(ncelet,*)

! Local variables

!< [loc_var_dec]
integer          ivart, iel, ifac
integer          ipcvis, ipccp
integer          ipcvsl, ivar
integer          f_id
integer          iclapc, idecal
integer          ll
integer          ipcte1, icla, iromf
integer          iscdri, keydri, iflid, nfld, keyccl

double precision xrtp
double precision aa1, bb1, cc1, dd1
double precision aa2, bb2, cc2, dd2
double precision aa3, bb3, cc3, dd3
double precision aa4, bb4, cc4, dd4
double precision aa5, bb5, cc5, dd5
double precision aa6, bb6, cc6, dd6
double precision aa7, bb7, cc7, dd7
double precision visco_O2, visco_CO, visco_H2, visco_N2
double precision visco_SO2, visco_NH3, visco_CO2
double precision g

character*80     fname

double precision, allocatable, dimension(:) :: x1, visco

double precision, dimension(:), pointer :: taup
double precision, dimension(:), pointer :: taupg
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! 0. Initializations to keep
!===============================================================================

!< [init]
ipcvis = ipproc(iviscl)
allocate(x1(ncelet), visco(ncelet))

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)

! Key id of the coal scalar class
call field_get_key_id("scalar_class", keyccl)

! Number of fields
call field_get_n_fields(nfld)
!< [init]

!===============================================================================

!   The following examples should be adapted by the user
!   ====================================================

!===============================================================================
!  Example:
!  =======
!===============================================================================

!  ===================================================================

!< [example_1]

g = sqrt(gx**2 + gy**2 + gz**2)

! Temperature
ipcte1 = ipproc(itemp1)

! Gas density
iromf = ipproc(irom1)

! First initialization
if (ntcabs.le.1) then
  do iel = 1, ncel
    visco(iel) = viscl0
    propce(iel,iromf) = ro0
    do icla = 1, nclacp
      propce(iel,ipproc(irom2(icla)))  = ro0
      propce(iel,ipproc(idiam2(icla))) = diam20(icla)
    enddo
  enddo
endif

!----------------Gas viscosity function of temperature--------------------------
!
!--------------------------1-O2 2-CO 3-H2 4-N2 5-SO2 6-NH3 7-CO2----------------

aa1 = 4.0495d-6
bb1 = 6.22d-8
cc1 = -2.3032d-11
dd1 = 4.4077d-15

aa2 = 9.9987d-6
bb2 = 5.1578d-8
cc2 = -1.8383d-11
dd2 = 3.33307d-15

aa3 = 2.894d-6
bb3 = 2.22508d-8
cc3 = -8.041d-12
dd3 = 1.4619d-15

aa4 = 4.3093d-6
bb4 = 5.0516d-8
cc4 = -1.7869d-11
dd4 = 3.2136d-15

aa5 = -1.9889d-6
bb5 = 5.365d-8
cc5 = -1.4286d-11
dd5 = 2.1639d-15

aa6 = -1.293d-6
bb6 = 4.1194d-8
cc6 = -1.772d-11
dd6 = 1.8699d-15

aa7 = 4.4822d-7
bb7 = 5.4327d-8
cc7 = -1.7581d-11
dd7 = 2.9979d-15

!-------------------------------------------------------------------------------
!      law                    mu   = a + b T + c T**2 + d T**3
!      so      propce(iel, ipcvis) = a +b*xrtp+c*xrtp**2 + d*xrtp**3
!-------------------------------------------------------------------------------

if (ntcabs.gt.1) then
  do iel = 1, ncel

    xrtp = propce(iel,ipcte1)
    visco_O2  = aa1 + xrtp*bb1 + cc1*xrtp**2 + dd1*xrtp**3
    visco_CO  = aa2 + xrtp*bb2 + cc2*xrtp**2 + dd2*xrtp**3
    visco_H2  = aa3 + xrtp*bb3 + cc3*xrtp**2 + dd3*xrtp**3
    visco_N2  = aa4 + xrtp*bb4 + cc4*xrtp**2 + dd4*xrtp**3
    visco_SO2 = aa5 + xrtp*bb5 + cc5*xrtp**2 + dd5*xrtp**3
    visco_NH3 = aa6 + xrtp*bb6 + cc6*xrtp**2 + dd6*xrtp**3
    visco_CO2 = aa7 + xrtp*bb7 + cc7*xrtp**2 + dd7*xrtp**3

    ! Viscosity of the mixing
    visco(iel) = ( propce(iel,ipproc(iym1(8))) * visco_O2                     &
                 + propce(iel,ipproc(iym1(3))) * visco_CO                     &
                 + propce(iel,ipproc(iym1(5))) * visco_H2                     &
                 + propce(iel,ipproc(iym1(12)))* visco_N2                     &
                 + propce(iel,ipproc(iym1(11)))* visco_SO2                    &
                 + propce(iel,ipproc(iym1(7))) * visco_NH3                    &
                 + propce(iel,ipproc(iym1(9))) * visco_CO2 )/                 &
                 ( propce(iel,ipproc(iym1(8))) + propce(iel,ipproc(iym1(3)))  &
                 + propce(iel,ipproc(iym1(5))) + propce(iel,ipproc(iym1(12))) &
                 + propce(iel,ipproc(iym1(11)))+ propce(iel,ipproc(iym1(7)))  &
                 + propce(iel,ipproc(iym1(9))))

  enddo
endif

! Compute x1 = 1 - sum x2
do iel = 1, ncel

  x1(iel) = 1.d0

  do icla = 1, nclacp
    x1(iel) = x1(iel) - propce(iel,ipproc(ix2(icla))) !FIXME is propce(iel,ipproc(ix2(icla))) initialized ?
  enddo

enddo

! All gas scalars have the same drift as if1m(1)
!-----------------------------------------------

do iflid = 0, nfld-1

  ! Index of the scalar class (<0 if the scalar belongs to the gas phase)
  call field_get_key_int(iflid, keyccl, icla)

  call field_get_key_int(iflid, keydri, iscdri)

  ! We only handle here one scalar with a drift per gas class
  if (icla.le.-1.and.btest(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)) then

    ! Position of variables, coefficients
    ! -----------------------------------

    ! Name of the drift scalar
    call field_get_name(iflid, fname)

    ! Index of the corresponding "relaxation time" (taupg) for the gas
    ! WARNING: for the gas, this tau might be negative
    call field_get_id('drift_tau_'//trim(fname), f_id)
    call field_get_val_s(f_id, taupg)

    ! Initialize to 0
    do iel = 1, ncel
      taupg(iel) = 0.d0
    enddo

  endif
enddo

! Loop over coal particle classes
! We only handle here coal class with a drift
!--------------------------------------------

do iflid = 0, nfld-1

  ! index of the coal particle class
  call field_get_key_int(iflid, keyccl, icla)

  call field_get_key_int(iflid, keydri, iscdri)

  ! We only handle here one scalar with a drift per particle class
  if (icla.ge.1.and.btest(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)) then

    ! Position of variables, coefficients
    ! -----------------------------------

    ! Name of the drift scalar
    call field_get_name(iflid, fname)

    ! Index of the corresponding relaxation time (taup)
    call field_get_id('drift_tau_'//trim(fname), f_id)
    call field_get_val_s(f_id, taup)

    ! Computation of the relaxation time of the particles
    ! the drift is therefore v_g = tau_p * g
    !----------------------------------------------------

    do iel = 1, ncel

      ! Simple model for Low Reynolds Numbers
      taup(iel) = x1(iel) * propce(iel,ipproc(irom2(icla)))          &
                          * propce(iel,ipproc(idiam2(icla)))**2      &
                          / (18.d0*visco(iel))

    enddo

    ! Drift for the gas:
    ! tau_pg = - Sum_i X2_i v_gi
    do iel = 1, ncel

      taupg(iel) = taupg(iel)                                       &
                 - ( taup(iel) * propce(iel,ipproc(ix2(icla))) )

    enddo

  endif
enddo

!< [example_1]

!Free memory
deallocate(x1, visco)

!===============================================================================

!----
! End
!----

return
end subroutine usphyv
