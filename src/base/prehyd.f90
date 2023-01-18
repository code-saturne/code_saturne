!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
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
!> (prediction and correction steps \ref navstv.f90).
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
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[out]    grdphd         the a priori hydrostatic pressure gradient
!>                              \f$ \partial _x (P_{hydro}) \f$
!> \param[in]     iterns        Navier-Stokes iteration number
!_______________________________________________________________________________

subroutine prehyd &
 ( grdphd  , iterns )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

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

double precision grdphd(ndim, ncelet)
integer          iterns

! Local variables

integer          iccocg, inc, isym  , f_id
integer          iel   , ifac
integer          nswrgp, imligp, iwarnp
integer          iflmas, iflmab
integer          idiffp, iconvp, ndircp
integer          ibsize
integer          iescap, ircflp, ischcp, isstpp, f_id0
integer          nswrsp
integer          imucpp, idftnp, iswdyp
integer          iharmo
integer          icvflb, hyd_p_flag
integer          ivoid(1)

double precision thetap
double precision epsrgp, climgp, extrap, epsilp
double precision hint, qimp, epsrsp, blencp, relaxp
double precision normp

double precision rvoid(1)

real(kind=c_double), dimension(:,:), pointer :: pvoid2
real(kind=c_double), dimension(1,1), target :: rvoid2
type(var_cal_opt), target   :: vcoptph
type(var_cal_opt), pointer  :: p_k_value
type(c_ptr)                 :: c_k_value

double precision, allocatable, dimension(:) :: coefap, cofafp, coefbp, cofbfp

double precision, allocatable, dimension(:) :: viscf, viscb
double precision, allocatable, dimension(:) :: xinvro
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:) :: smbr, rovsdt
double precision, dimension(:), pointer :: imasfl, bmasfl, prhyd
double precision, dimension(:), pointer :: crom

type(var_cal_opt) :: vcopt

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

call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
f_id0  = -1
iconvp = 0
idiffp = 1
ndircp = 0
nswrsp = 1           ! no reconstruction gradient
nswrgp = vcopt%nswrgr
imligp = vcopt%imligr
ircflp = vcopt%ircflu
ischcp = vcopt%ischcv
isstpp = vcopt%isstpc
iescap = 0
imucpp = 0
idftnp = ISOTROPIC_DIFFUSION
iswdyp = vcopt%iswdyn
iwarnp = vcopt%iwarni
blencp = vcopt%blencv
epsilp = vcopt%epsilo
epsrsp = vcopt%epsrsm
epsrgp = vcopt%epsrgr
climgp = vcopt%climgr
extrap = 0.d0
relaxp = vcopt%relaxv
thetap = vcopt%thetav
! all boundary convective flux with upwind
icvflb = 0
normp = -1.d0

! Solve the diffusion equation
p_k_value => vcoptph
c_k_value = c_loc(p_k_value)

vcoptph%iwarni = iwarnp

vcoptph%iconv  = iconvp
vcoptph%istat  = 0

vcoptph%icoupl = -1
vcoptph%ndircl = ndircp
vcoptph%idiff  = idiffp
vcoptph%idifft = -1
vcoptph%idften = idftnp
vcoptph%iswdyn = iswdyp
vcoptph%ischcv = ischcp
vcoptph%isstpc = isstpp
vcoptph%nswrgr = nswrgp
vcoptph%nswrsm = nswrsp ! Important for mesh with reconstruction
vcoptph%imrgra = imrgra
vcoptph%imligr = imligp
vcoptph%ircflu = ircflp
vcoptph%iwgrec = 0      ! Warning, may be overwritten if a field
vcoptph%thetav = thetap
vcoptph%blencv = blencp
vcoptph%blend_st = 0    ! Warning, may be overwritten if a field
vcoptph%epsilo = epsilp
vcoptph%epsrsm = epsrsp
vcoptph%epsrgr = epsrgp
vcoptph%climgr = climgp
vcoptph%relaxv = relaxp

call cs_equation_iterative_solve_scalar          &
 ( idtvar , iterns ,                             &
   f_id0    , "Prhydro"      ,                   &
   iescap , imucpp , normp  , c_k_value       ,  &
   prhyd       , prhyd      ,                    &
   coefap , coefbp , cofafp , cofbfp ,           &
   imasfl , bmasfl ,                             &
   viscf  , viscb  , viscf  , viscb  ,           &
   rvoid , rvoid , rvoid ,                       &
   icvflb , ivoid  ,                             &
   rovsdt , smbr  , prhyd        , dpvar  ,  &
   rvoid   , rvoid  )

! Free memory
deallocate(dpvar)

inc    = 1
iccocg = 1
nswrgp = 1
extrap = 0.d0
f_id = -1

hyd_p_flag = 0

call gradient_weighted_s &
 ( f_id   , imrgra , inc    , nswrgp , imligp ,                  &
   hyd_p_flag,                                                   &
   iwarnp , epsrgp , climgp , pvoid2 ,                           &
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
