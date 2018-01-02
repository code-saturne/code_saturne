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

!> \file richards.f90
!>
!> \brief This routine solves the Richards equation, then compute the new
!> velocities deducted from the gradients of hydraulic head and
!> from the permeability.
!> These velocities are used for post-processing, calculation of dispersion
!> coefficients, convergence criterion of Newton scheme... but not
!> for transport.
!> In order to ensure the exact conservation of mass, the mass fluxes are
!> computed following the procedure of the standard subroutine resopv
!> (See the <a href="../../theory.pdf#resopv"><b>resopv</b></a>
!> section of the theory guide for more informations).
!>
!> The hydraulic head is contained in the pressure array, while the velocities
!> are contained in the velocity arrays.
!>
!> The Richards equation for underground flows is :
!> d theta(h) / dt - div( K(h) grad(h) ) = 0 (or source(h) in very
!> particular cases).
!> This equation can also be written, denoting C(h) = d theta(h)/dh :
!> C(h) dh/dt - div( K(h) grad(h) ) = C(h) dh/dt - d theta(h) / dt (1)
!> The right hand term is close to zero and is not considered in ESTEL.
!> We consider it here for exact mass conservation.
!>
!> (1) is solved with a 'Newton scheme' handled by tridim. The structure used
!> for this Newton scheme is the same as the one used in the loop over
!> nterup for standard flows.
!> For going from time step n to time step n+1, we use sub-iterations
!> indexed by m :
!> C(h^n) (h^(n+1,m+1)-h^n)/detla t - div( K(h^(n+1,m)) grad(h^(n+1,m+1)) )
!>  = C(h^n) (h^(n+1,m)-h^n)/detla t - (theta(h^(n+1,m))-theta(h(n)))/delta t.
!> These sub-iterations, if they converge, converge to the solution i
!> of the problem.
!>
!> The Darcy velocity q is then computed thanks to the relation :
!>   q = -K(h) grad(h).
!>
!> This routine is essentially inspired from navstv, resopv and codits.
!>
!> Please refer to the <a href="../../theory.pdf#groundwater"><b>groundwater flows</b></a>
!> section of the theory guide for more theoretical informations.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     icvrge        indicator of convergence
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________


subroutine richards &
 (icvrge, dt)

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use entsor
use cstphy
use cstnum
use optcal
use pointe
use albase
use parall
use period
use mesh
use field
use field_operator
use cs_c_bindings
use darcy_module
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer icvrge
double precision, pointer, dimension(:)   :: dt

! Local variables

integer iccocg, inc   , iel, isou, init
integer nswrgp, imligp, iwarnp
integer imucpp, ircflp, isweep, isym, lchain
integer ndircp, niterf, nswmpr
integer iinvpe, iflmas, iflmab, iesize, idiffp, iconvp, ibsize
integer fid
integer iflid , iflwgr, f_dim, f_id0, iwgrp, iprev, iitsm

double precision epsrgp, climgp, extrap
double precision thetap, xdu, xdv, xdw, xnrmul
double precision relaxp, residu, ressol, epsilp

character(len=80) :: chaine

type(solving_info) sinfo
type(var_cal_opt) :: vcopt_p

double precision rvoid(1)

double precision, dimension(:,:), allocatable :: gradp
double precision, allocatable, dimension(:,:) :: xam, weighf, uvwk
double precision, allocatable, dimension(:) :: dpvar, dam, presa, rhs, w1, weighb
double precision, allocatable, dimension(:) :: rovsdt, rhs0, viscf, viscb

double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: coefa_p, coefb_p
double precision, dimension(:), pointer :: coefaf_p, coefbf_p
double precision, dimension(:), pointer :: cpro_permeability, cproa_capacity
double precision, dimension(:), pointer :: cproa_sat, cpro_sat
double precision, dimension(:,:), pointer :: cpro_permeability_6
double precision, dimension(:,:), pointer :: cvar_vel
double precision, dimension(:), pointer :: cvar_pr, cvara_pr
double precision, dimension(:), pointer :: cpro_prtot
double precision, dimension(:,:), pointer :: cpro_wgrec_v
double precision, dimension(:), pointer :: cpro_wgrec_s

!===============================================================================
! 0. Initialization
!===============================================================================

if (darcy_convergence_criterion.eq.0) then
  allocate(uvwk(1,ncelet))
else
  allocate(uvwk(3,ncelet))
endif
allocate(dpvar(ncelet))
allocate(presa(ncelet))
allocate(rhs(ncelet), rovsdt(ncelet), rhs0(ncelet))
allocate(viscf(nfac), viscb(ndimfb))

! Map field arrays
call field_get_val_v(ivarfl(iu), cvar_vel)

if (darcy_anisotropic_permeability.eq.0) then
  call field_get_val_s_by_name('permeability', cpro_permeability)
else
  call field_get_id('permeability', fid)
  call field_get_val_v(fid, cpro_permeability_6)
endif

call field_get_val_prev_s_by_name('capacity', cproa_capacity)
call field_get_val_prev_s_by_name('saturation', cproa_sat)
call field_get_val_s_by_name('saturation', cpro_sat)
call field_get_val_s(ivarfl(ipr), cvar_pr)
call field_get_val_prev_s(ivarfl(ipr), cvara_pr)

sinfo%nbivar = 0

! Preparation of convergence criterion for Newton scheme
if (nterup.gt.1) then
  !$omp parallel do
  do iel = 1,ncel
    if (darcy_convergence_criterion.eq.0) then
      uvwk(1,iel) = cvar_pr(iel)
    else
      uvwk(1,iel) = cvar_vel(1,iel)
      uvwk(2,iel) = cvar_vel(2,iel)
      uvwk(3,iel) = cvar_vel(3,iel)
    endif
  enddo
  xnrmul = 0.d0
  !$omp parallel do reduction(+:xnrmul)
  do iel = 1, ncel
    if (darcy_convergence_criterion.eq.0) then
      xnrmul = xnrmul +cvar_pr(iel)**2*volume(iel)
    else
      xnrmul = xnrmul +(cvar_vel(1,iel)**2 &
                      + cvar_vel(2,iel)**2 &
                      + cvar_vel(2,iel)**2)*volume(iel)
    endif
  enddo
  if (irangp.ge.0) then
    call parsom (xnrmul)
  endif
  xnrmu0 = sqrt(xnrmul)
endif

! Index of the field
iflid = ivarfl(ipr)

call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_p)

if (vcopt_p%iwgrec.eq.1) then
  ! Id weighting field for gradient
  call field_get_key_int(iflid, kwgrec, iflwgr)
  call field_get_dim(iflwgr, f_dim)
  if (f_dim.gt.1) then
    call field_get_val_v(iflwgr, cpro_wgrec_v)
  else
    call field_get_val_s(iflwgr, cpro_wgrec_s)
  endif
endif

! Preparation of writing
call field_get_name(ivarfl(ipr), chaine)
lchain = 16

f_id0 = -1
iwgrp = 0
if (vcopt_p%iwgrec.eq.1) iwgrp = 1

! --- Boundary conditions
call field_get_coefa_s(ivarfl(ipr), coefa_p)
call field_get_coefb_s(ivarfl(ipr), coefb_p)
call field_get_coefaf_s(ivarfl(ipr), coefaf_p)
call field_get_coefbf_s(ivarfl(ipr), coefbf_p)

! --- Physical fluxes
call field_get_key_int(ivarfl(ipr), kimasf, iflmas)
call field_get_key_int(ivarfl(ipr), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)


! --- Solving options. Copied from resopv.
isym  = 1
if (darcy_unsteady.eq.1) isym  = 2

! Matrix block size. Copied from reopv.
ibsize = 1
iesize = 1

allocate(dam(ncelet), xam(isym,nfac))

!===============================================================================
! 1. Building of the linear system to solve
!===============================================================================

! Unsteady term
if (darcy_unsteady.eq.1) then
  do iel = 1, ncel
    rovsdt(iel) = volume(iel)*cproa_capacity(iel)/dt(iel)
  enddo
else
  do iel = 1, ncel
    rovsdt(iel) = 0.d0
  enddo
endif

! We keep the capacity in explicit because the Newton scheme converges
! much better like this.

! Computation of diffusion coefficients at the centers of faces
if (darcy_anisotropic_permeability.eq.0) then
  call viscfa &
  ( imvisf ,                                  &
    cpro_permeability     ,                   &
    viscf  , viscb  )
  if (vcopt_p%iwgrec.eq.1) then
    ! Weighting for gradient
    do iel = 1, ncel
      cpro_wgrec_s(iel) = cpro_permeability(iel)
    enddo
    call synsca(cpro_wgrec_s)
  endif
else if (darcy_anisotropic_permeability.eq.1) then
  allocate(weighf(2,nfac))
  allocate(weighb(ndimfb))
  iwarnp = vcopt_p%iwarni
  call vitens &
  ( cpro_permeability_6 , iwarnp ,       &
    weighf , weighb ,                    &
    viscf  , viscb  )
  if (vcopt_p%iwgrec.eq.1) then
    ! Weighting for gradient
    do iel = 1, ncel
      do isou = 1, 6
        cpro_wgrec_v(isou,iel) = cpro_permeability_6(isou,iel)
      enddo
    enddo
    call syntis(cpro_wgrec_v)
  endif
endif

iconvp = 0 ! no convection
idiffp = 1 ! diffusion
ndircp = ndircl(ipr) ! no diagonal stepped aside
thetap = 1.d0 ! implicit scheme
imucpp = 0 ! do not multiply the convectiv term by anything

! Note that the computed matrix is always the same in the loop over iterns
! (loop of the Newton scheme, handled by tridim)
! Another option would be to store the arrays dam and xam in a global array
! rather than computing the matrix at each call to 'richards'.
call matrix &
( iconvp , idiffp , ndircp , isym   ,                            &
  thetap , imucpp ,                                              &
  coefb_p , coefbf_p     , rovsdt ,                              &
  imasfl , bmasfl , viscf  , viscb  ,                            &
  rvoid  , dam    , xam    )

!===============================================================================
! 2. Solving (Loop over the non-orthogonalities)
!===============================================================================

! For this part, see documentation on resopv.
! The way used to handle non-orthogonalities is the same.
! The difference is the instationnary term.

nswmpr = vcopt_p%nswrsm

! Initialization of dpvar for avoiding warnings
do iel = 1, ncel
  dpvar(iel) = 0.d0
enddo

! We compute the first residue (rhs)

relaxp = vcopt_p%relaxv ! relaxation for loop over isweep
isweep = 1 ! counter for non-orthogonalities
iccocg = 1 ! no calculation of cocg. What does it mean?
init = 1 ! it is an initialization of the residue
inc  = 1 ! 0 increment, 1 otherwise
nswrgp = vcopt_p%nswrgr ! number of sweeps for gradients reconstruction
imligp = vcopt_p%imligr
iwarnp = vcopt_p%iwarni
epsrgp = vcopt_p%epsrgr
climgp = vcopt_p%climgr
extrap = vcopt_p%extrag
ircflp = vcopt_p%ircflu

if (darcy_anisotropic_permeability.eq.0) then

  call itrgrp &
( ivarfl(ipr), init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
  iwarnp ,                                                                    &
  epsrgp , climgp , extrap ,                                                  &
  rvoid  ,                                                                    &
  cvar_pr   ,                                                                 &
  coefa_p  , coefb_p  ,                                                       &
  coefaf_p , coefbf_p ,                                                       &
  viscf  , viscb  ,                                                           &
  cpro_permeability,                                                          &
  rhs   )

else if (darcy_anisotropic_permeability.eq.1) then

  call itrgrv &
( ivarfl(ipr), init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp ,        &
  iphydr , iwarnp ,                                                                  &
  epsrgp , climgp , extrap ,                                                         &
  rvoid  ,                                                                           &
  cvar_pr   ,                                                                        &
  coefa_p  , coefb_p  ,                                                              &
  coefaf_p , coefbf_p ,                                                              &
  viscf  , viscb  ,                                                                  &
  cpro_permeability_6 ,                                                              &
  weighf , weighb ,                                                                  &
  rhs   )

endif

do iel = 1, ncel
  if (darcy_unsteady.eq.1) then
    ! We take care of exact mass conservation (as precise as the convergence
    ! of the Newton loop is good)
    rhs0(iel) = cproa_capacity(iel)*volume(iel)/dt(iel)                   &
       *cvar_pr(iel) + volume(iel)/dt(iel)                                &
       *(cproa_sat(iel) - cpro_sat(iel))
  else
    rhs0(iel) = 0.d0
  endif
enddo

if (ncetsm.gt.0) then
  !$omp parallel do private(iel) if(ncetsm > thr_n_min)
  do iitsm = 1, ncetsm
    iel = icetsm(iitsm)
    rhs0(iel) = rhs0(iel) + volume(iel)*smacel(iitsm,ipr)
  enddo
endif

do iel = 1, ncel
  if (darcy_unsteady.eq.1) then
    rhs(iel) = rhs0(iel) - rhs(iel)                                        &
             - volume(iel)/dt(iel)*cproa_capacity(iel)*cvar_pr(iel)
  else
    rhs(iel) = rhs0(iel) - rhs(iel)
  endif
enddo

residu = sqrt(cs_gdot(ncel, rhs, rhs))

sinfo%rnsmbr = residu

! We compute the normalisation residue, which is used as a stop criterion
! in the loop of non-orthogonalities. We have to "normalize" the problem,
! taking into account the boundary conditions.
! This part is inspired from codits (call to promav)
allocate(w1(ncelet))
iinvpe = 1 ! for processing communication before calculation
call promav(isym, ibsize, iesize, iinvpe, ivarfl(ipr), dam, xam, cvar_pr, w1)

!$omp parallel do
do iel = 1, ncel
  w1(iel) = w1(iel) + rhs(iel)
enddo
rnormp = sqrt(cs_gdot(ncel, w1, w1))

! Free memory
deallocate(w1)

! Writing
if (vcopt_p%iwarni.ge.2) then
  write(nfecra,1400)chaine(1:16), isweep, residu, relaxp
endif

! Loop for non-orthogonalities
! TODO: dynamic relaxation
nswmpr = vcopt_p%nswrsm
do while ( (isweep.le.nswmpr.and.residu.gt.vcopt_p%epsrsm*rnormp) &
            .or. (isweep.eq.1) )
            ! We pass at least once to ensure exactness of mass flux

  iwarnp = vcopt_p%iwarni
  epsilp = vcopt_p%epsilo
  iinvpe = 1
  ressol = residu

  call sles_solve_native(ivarfl(ipr), chaine,                         &
                         isym, ibsize, iesize, dam, xam, iinvpe,      &
                         epsilp, rnormp, niterf, ressol, rhs, dpvar)

  if (isweep.le.nswmpr.and.residu.gt.vcopt_p%epsrsm*rnormp) then
    do iel = 1, ncel
      presa(iel) = cvar_pr(iel)
      cvar_pr(iel) = presa(iel) + vcopt_p%relaxv*dpvar(iel)
    enddo

    ! If it is the last sweep, update with the total increment
  else
    do iel = 1, ncel
      presa(iel) = cvar_pr(iel)
      cvar_pr(iel) = presa(iel) + dpvar(iel)
    enddo
  endif

  iccocg = 1
  init = 1
  inc  = 1

  if (darcy_anisotropic_permeability.eq.0) then

    call itrgrp &
  ( ivarfl(ipr), init  , inc , imrgra , iccocg , nswrgp , imligp , iphydr ,  &
    iwarnp ,                                                                 &
    epsrgp , climgp , extrap ,                                               &
    rvoid  ,                                                                 &
    cvar_pr   ,                                                              &
    coefa_p  , coefb_p  ,                                                    &
    coefaf_p , coefbf_p ,                                                    &
    viscf  , viscb  ,                                                        &
    cpro_permeability,                                                       &
    rhs   )

  else if (darcy_anisotropic_permeability.eq.1) then

    call itrgrv &
    ( ivarfl(ipr), init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp ,        &
    iphydr , iwarnp ,                                                                    &
    epsrgp , climgp , extrap ,                                                           &
    rvoid  ,                                                                             &
    cvar_pr   ,                                                                          &
    coefa_p  , coefb_p  ,                                                                &
    coefaf_p , coefbf_p ,                                                                &
    viscf  , viscb  ,                                                                    &
    cpro_permeability_6 ,                                                                &
    weighf , weighb ,                                                                    &
    rhs   )

  endif

  do iel = 1, ncel
    if (darcy_unsteady.eq.1) then
      rhs(iel) = rhs0(iel) - rhs(iel)                                    &
               - volume(iel)/dt(iel)*cproa_capacity(iel)*cvar_pr(iel)
    else
      rhs(iel) = rhs0(iel) - rhs(iel)
    endif
  enddo

  ! --- Convergence test
  residu = sqrt(cs_gdot(ncel, rhs, rhs))

  ! Writing
  sinfo%nbivar = sinfo%nbivar + niterf

  if (vcopt_p%iwarni.ge.2) then
    write(nfecra,1400)chaine(1:16), isweep, residu, relaxp
  endif

  if (iwarnp.ge.3) then
    write(nfecra,1500) chaine(1:16), isweep, residu, rnormp, niterf
  endif

  isweep = isweep + 1

enddo
! --- Reconstruction loop (end)

if (abs(rnormp).gt.epzero) then
  sinfo%resvar = residu/rnormp
  sinfo%dervar = - sinfo%rnsmbr/rnormp !FIXME
else
  sinfo%resvar = 0.d0
  sinfo%dervar = 0.d0
endif

! Writing
if (vcopt_p%iwarni.ge.2) then
  if (isweep.gt.nswmpr) then
    write(nfecra,1600) chaine(1:16), nswmpr
  endif
endif

call sles_free_native(ivarfl(ipr), chaine)
deallocate(dam, xam, rhs)

!===============================================================================
! 3. Updating of mass fluxes
!===============================================================================

iccocg = 1
init = 1
inc = 1
nswrgp = vcopt_p%nswrgr
imligp = vcopt_p%imligr
iwarnp = vcopt_p%iwarni
epsrgp = vcopt_p%epsrgr
climgp = vcopt_p%climgr
extrap = vcopt_p%extrag

! We compute the new mass flux, taking care not to reconstruct
! the last increment of the lopp on isweep, to ensure an
! exact conservation of the mass.

if (darcy_anisotropic_permeability.eq.0) then

  call itrmas &
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr , &
   iwgrp  , iwarnp ,                                              &
   epsrgp , climgp , extrap ,                                     &
   rvoid  ,                                                       &
   presa  ,                                                       &
   coefa_p , coefb_p , coefaf_p , coefbf_p ,                      &
   viscf  , viscb  ,                                              &
   cpro_permeability,                                             &
   imasfl , bmasfl )

  ! The last increment is not reconstructed to fullfill exactly the continuity
  ! equation (see theory guide).
  iccocg = 0
  nswrgp = 0
  inc = 0
  init = 0

  call itrmas &
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr ,     &
   iwgrp  , iwarnp ,                                                           &
   epsrgp , climgp , extrap ,                                                  &
   rvoid  ,                                                                    &
   dpvar  ,                                                                    &
   coefa_p , coefb_p , coefaf_p , coefbf_p ,                                   &
   viscf  , viscb  ,                                                           &
   cpro_permeability,                                                          &
   imasfl , bmasfl )

else if (darcy_anisotropic_permeability.eq.1) then

  call itrmav &
 ( f_id0, init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
   iphydr , iwgrp  , iwarnp ,                                            &
   epsrgp , climgp , extrap ,                                            &
   rvoid  ,                                                              &
   presa  ,                                                              &
   coefa_p , coefb_p , coefaf_p , coefbf_p ,                             &
   viscf  , viscb  ,                                                     &
   cpro_permeability_6 ,                                                 &
   weighf , weighb ,                                                     &
   imasfl , bmasfl )

  ! The last increment is not reconstructed to fullfill exactly the continuity
  ! equation (see theory guide).
  iccocg = 0
  nswrgp = 0
  inc = 0
  init = 0
  ircflp = 0

  call itrmav &
 ( f_id0, init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp ,  &
   iphydr , iwgrp  , iwarnp ,                                             &
   epsrgp , climgp , extrap ,                                             &
   rvoid  ,                                                               &
   dpvar  ,                                                               &
   coefa_p , coefb_p , coefaf_p , coefbf_p ,                              &
   viscf  , viscb  ,                                                      &
   cpro_permeability_6 ,                                                  &
   weighf , weighb ,                                                      &
   imasfl , bmasfl )

endif

!===============================================================================
! 4.  Updating of the velocity field and the pressure head
!===============================================================================

! We compute the gradient of hydraulique head and multiply it by the
! cpro_permeability in order to get the new velocities. These velocities will
! only be used for post-processing, computation of the dispersion coefficients
! of scalars and possibly convergence criterion of the Newton scheme. The
! transport equation uses the mass fluxes at the center of faces, and not
! the velocities at the center of cells.

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(cvar_pr)
endif

iccocg = 1
inc = 1
iprev  = 0

allocate(gradp(3,ncelet))

! We use gradient_scalar instead of potential because iphydr is
! not taken into account in case of Darcy calculation.
call field_gradient_scalar(ivarfl(ipr), iprev, imrgra, inc,    &
                           iccocg, gradp)

! Computation of velocity
if (darcy_anisotropic_permeability.eq.1) then

  !$omp parallel do
  do iel = 1, ncel
    cvar_vel(1,iel) = - ( cpro_permeability_6(1,iel)*gradp(1,iel)          &
                         + cpro_permeability_6(4,iel)*gradp(2,iel)     &
                         + cpro_permeability_6(6,iel)*gradp(3,iel))
    cvar_vel(2,iel) = - ( cpro_permeability_6(4,iel)*gradp(1,iel)          &
                         + cpro_permeability_6(2,iel)*gradp(2,iel)     &
                         + cpro_permeability_6(5,iel)*gradp(3,iel))
    cvar_vel(3,iel) = - (cpro_permeability_6(6,iel)*gradp(1,iel)           &
                         + cpro_permeability_6(5,iel)*gradp(2,iel)     &
                         + cpro_permeability_6(3,iel)*gradp(3,iel))
  enddo

else

  !$omp parallel do
  do iel = 1, ncel
    cvar_vel(1,iel) = - cpro_permeability(iel)*gradp(1,iel)
    cvar_vel(2,iel) = - cpro_permeability(iel)*gradp(2,iel)
    cvar_vel(3,iel) = - cpro_permeability(iel)*gradp(3,iel)
  enddo

endif

! update pressure head (h = H - z) for post-processing
! Only used when gravity is taken into account
if (darcy_gravity.ge.1) then
  call field_get_val_s(iprtot, cpro_prtot)
  !$omp parallel do
  do iel = 1, ncel
    cpro_prtot(iel) = cvar_pr(iel) - xyzcen(1,iel)*darcy_gravity_x             &
                                   - xyzcen(2,iel)*darcy_gravity_y             &
                                   - xyzcen(3,iel)*darcy_gravity_z
  enddo

endif

!===============================================================================
! 5.  Checking of convergence criterion for the Newton scheme
!===============================================================================

! Computation of the convergence criterion for the Newton sheme. If
! the convergence is reached, we set icvrge=1, which will be used by tridim
! for stopping the loop.

if (nterup.gt.1) then
  icvrge = 1
  xnrmul = 0.d0
  !$omp parallel do reduction(+:xnrmul) private(xdu, xdv, xdw)
  do iel = 1,ncel
    if (darcy_convergence_criterion.eq.0) then
      xdu = cvar_pr(iel) - uvwk(1,iel)
      xnrmul = xnrmul + xdu**2*volume(iel)
    else
      xdu = cvar_vel(1,iel) - uvwk(1,iel)
      xdv = cvar_vel(2,iel) - uvwk(2,iel)
      xdw = cvar_vel(3,iel) - uvwk(3,iel)
      xnrmul = xnrmul + (xdu**2 + xdv**2 + xdw**2)*volume(iel)
    endif
  enddo
  if (irangp.ge.0) call parsom (xnrmul)
  xnrmu = sqrt(xnrmul)
  ! Indicator of convergence for the Newton scheme
  if (xnrmu.ge.epsup*xnrmu0) icvrge = 0
endif

!===============================================================================
! 6.  Finalization
!===============================================================================

! Save convergence info
call field_set_key_struct_solving_info(ivarfl(ipr), sinfo)

! Free memory
deallocate(gradp)
deallocate(uvwk)
deallocate(dpvar)
deallocate(presa)
deallocate(rovsdt)
if (allocated(weighf)) deallocate(weighf, weighb)
deallocate(viscf, viscb)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)
 1400 format(1X,A16,' : SWEEP = ',I5,' NORME SECOND MEMBRE = ',E14.6,  &
             ', RELAXP = ',E14.6)
 1500 format ( &
 1X,A16,' : Current reconstruction sweep = ',I5                     ,/,&
'           sweep residual = ',E12.5,', norm = ',E12.5              ,/,&
'           number of sweeps for solver = ',I5)
 1600 format( &
'@'                                                                 ,/,&
'@ @@ ATTENTION : ', A16,' RICHARDS'                                ,/,&
'@    ========='                                                    ,/,&
'@  Nombre d''iterations maximal ',I10   ,' atteint'                ,/,&
'@' )
#else
 1400 format(1X,A16,' : SWEEP = ',I5,' RIGHT HAND SIDE NORM = ',E14.6, &
             ', RELAXP = ',E14.6)
 1500 format ( &
 1X,A16,' : Current reconstruction sweep = ',I5                     ,/,&
'           sweep residual = ',E12.5,', norm = ',E12.5              ,/,&
'           number of sweeps for solver = ',I5)
 1600 format( &
'@'                                                                 ,/,&
'@ @@ WARNING: ', A16,' RICHARDS'                                   ,/,&
'@    ========'                                                     ,/,&
'@  Maximum number of iterations ',I10   ,' reached'                ,/,&
'@'                                                              )
#endif

!--------
! End
!--------

return
end subroutine richards
