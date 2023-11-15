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

!===============================================================================
! Function:
! ---------

!> \file prehyd.f90
!>
!> \brief Compute an "a priori" hydrostatic pressure and its gradient associated
!> before the Navier Stokes equations
!> (prediction and correction steps \ref cs_solve_navier_stokes.c).
!>
!> This function computes a hydrostatic pressure \f$ P_{hydro} \f$ solving an
!> a priori simplified momentum equation:
!> \f[
!> \rho^n \dfrac{(\vect{u}^{hydro} - \vect{u}^n)}{\Delta t} =
!> \rho^n \vect{g}^n - \grad P_{hydro}
!> \f]
!> and using the mass equation as following:
!> \f[
!> \rho^n \divs \left( \delta \vect{u}_{hydro} \right) = 0
!> \f]
!> with: \f$ \delta \vect{u}_{hydro} = ( \vect{u}^{hydro} - \vect{u}^n) \f$
!>
!> finally, we resolve the simplified momentum equation below:
!> \f[
!> \divs \left( K \grad P_{hydro} \right) = \divs \left(\vect{g}\right)
!> \f]
!> with the diffusion coefficient (\f$ K \f$) defined as:
!> \f[
!>      K \equiv \dfrac{1}{\rho^n}
!> \f]
!> with a Neumann boundary condition on the hydrostatic pressure:
!> \f[
!>    D_\fib \left( K, \, P_{hydro} \right) =
!>    \vect{g} \cdot \vect{n}_\ib
!> \f]
!> (see the theory guide for more details on the boundary condition
!>  formulation).
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    grdphd         the a priori hydrostatic pressure gradient
!>                              \f$ \partial _x (P_{hydro}) \f$
!> \param[in]     iterns        Navier-Stokes iteration number
!_______________________________________________________________________________

subroutine prehyd  &
  (grdphd, iterns) &
  bind(C, name='cs_hydrostatic_pressure_prediction')

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

use paramx
use numvar
use entsor
use cstphy
use cstnum
use optcal
use albase
use parall
use period
use lagran
use cplsat
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

real(c_double) :: grdphd(ndim, ncelet)
integer(c_int), value :: iterns

! Local variables

integer          inc, isym  , f_id
integer          iel   , ifac
integer          nswrgp, iwarnp
integer          iflmas, iflmab
integer          ibsize
integer          iescap, f_id0, f_id_post
integer          imucpp
integer          iharmo
integer          icvflb, hyd_p_flag
integer          ivoid(1)

double precision hint, qimp
double precision normp

double precision rvoid(1)

real(kind=c_double), dimension(:,:), pointer :: pvoid2
real(kind=c_double), dimension(1,1), target :: rvoid2

double precision, allocatable, dimension(:) :: coefap, cofafp, coefbp, cofbfp

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: xinvro
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:) :: smbr, rovsdt
double precision, dimension(:), pointer :: imasfl, bmasfl, prhyd
double precision, dimension(:), pointer :: crom, post_phydr

type(var_cal_opt) :: vcopt
type(var_cal_opt), target   :: vcopt_loc
type(var_cal_opt), pointer  :: p_k_value
type(c_ptr)                 :: c_k_value

procedure() :: set_neumann_scalar, compute_td_pressure_perfect_gas

!===============================================================================
! 1. Initializations
!===============================================================================

! Map arrays
call field_get_val_s_by_name('hydrostatic_pressure_prd', prhyd)

! Allocate temporary arrays

! Boundary conditions for delta P
allocate(coefap(nfabor), cofafp(nfabor), coefbp(nfabor), cofbfp(nfabor))

! --- Physical properties
call field_get_val_s(icrom, crom)

call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

! --- Resolution options
isym  = 1
if (vcopt%iconv.gt.0) then
  isym  = 2
endif

! --- Matrix block size
ibsize = 1

! Void arrays
pvoid2 => rvoid2

!===============================================================================
! 2. Solving a diffusion equation with source term to obtain
!    the a priori hydrostatic pressure
!===============================================================================

! --- Allocate temporary arrays
allocate(dpvar(ncelet))
allocate(viscf(nfac), viscb(nfabor))
allocate(xinvro(ncelet))
allocate(smbr(ncelet), rovsdt(ncelet))

! --- Initialization of the variable to solve from the interior cells
do iel = 1, ncel
  xinvro(iel) = 1.d0/crom(iel)
  rovsdt(iel) = 0.d0
  smbr(iel)   = 0.d0
enddo

! --- Viscosity (k_t := 1/rho )
iharmo = 1
call viscfa (iharmo, xinvro, viscf, viscb)

! Neumann boundary condition for the pressure increment
!------------------------------------------------------

do ifac = 1, nfabor

  iel = ifabor(ifac)

  ! Prescribe the pressure gradient: kt.grd(Phyd)|_b = (g.n)|_b

  hint = 1.d0 /(crom(iel)*distb(ifac))
  qimp = - (gx*surfbo(1,ifac)                &
           +gy*surfbo(2,ifac)                &
           +gz*surfbo(3,ifac))/(surfbn(ifac))

  call set_neumann_scalar(coefap(ifac), cofafp(ifac),      &
                          coefbp(ifac), cofbfp(ifac),      &
                          qimp, hint)

enddo

!--------------------------------------------------------------------------
! Solve the diffusion equation

! By default, the hydrostatic pressure variable is resolved with 5 sweeps
! for the reconstruction gradient. Here we make the assumption that the
! mesh is orthogonal (any reconstruction gradient is done for the
! hydrostatic pressure variable)

! We do not yet use the multigrid to resolve the hydrostatic pressure
!--------------------------------------------------------------------------

f_id0  = -1
iescap = 0
imucpp = 0
! all boundary convective flux with upwind
icvflb = 0
normp = -1.d0

! Solve the diffusion equation

call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
vcopt_loc = vcopt

vcopt_loc%iconv = 0
vcopt_loc%istat = 0
vcopt_loc%icoupl = -1
vcopt_loc%ndircl = 0
vcopt_loc%idiff  = 1
vcopt_loc%idifft = -1
vcopt_loc%idften = ISOTROPIC_DIFFUSION
vcopt_loc%nswrsm = 1      ! no reconstruction gradient
                          ! (important for mesh with reconstruction)
vcopt_loc%imrgra = imrgra
vcopt_loc%iwgrec = 0      ! Warning, may be overwritten if a field
vcopt_loc%blend_st = 0    ! Warning, may be overwritten if a field

p_k_value => vcopt_loc
c_k_value = equation_param_from_vcopt(c_loc(p_k_value))

call cs_equation_iterative_solve_scalar      &
 ( idtvar , iterns ,                         &
   f_id0  , "Prhydro"       ,                &
   iescap , imucpp , normp  , c_k_value   ,  &
   prhyd  , prhyd  ,                         &
   coefap , coefbp , cofafp , cofbfp ,       &
   imasfl , bmasfl ,                         &
   viscf  , viscb  , viscf  , viscb  ,       &
   rvoid , rvoid , rvoid ,                   &
   icvflb , ivoid  ,                         &
   rovsdt , smbr  , prhyd        , dpvar  ,  &
   rvoid  , rvoid  )

! Save hydrostatic pressure for postprocessing if requested
call field_get_id_try("phydr", f_id_post)
if (f_id_post.ge.0) then
  call field_get_val_s(f_id_post,post_phydr)
  do iel=1, ncel
    post_phydr(iel) = prhyd(iel)
  enddo
endif

! Free memory
deallocate(dpvar)

inc    = 1
nswrgp = 1
f_id = -1

hyd_p_flag = 0

call gradient_weighted_s &
 ( f_id   , imrgra , inc    , nswrgp , vcopt_loc%imligr ,        &
   hyd_p_flag,                                                   &
   iwarnp , vcopt_loc%epsrgr , vcopt_loc%climgr , pvoid2 ,       &
   prhyd  , xinvro , coefap , coefbp ,                           &
   grdphd   )

!===============================================================================
! Free memory
!===============================================================================

deallocate(coefap, cofafp, coefbp, cofbfp)
deallocate(viscf, viscb)
deallocate(xinvro)
deallocate(smbr, rovsdt)

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine
