!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

!> \file coditv.f90
!>
!> \brief This function solves an advection diffusion equation with source terms
!> for one time step for the vector variable \f$ \vect{a} \f$.
!>
!> The equation reads:
!>
!> \f[
!> \tens{f_s}^{imp}(\vect{a}^{n+1}-\vect{a}^n)
!> + \divv \left( \vect{a}^{n+1} \otimes \rho \vect {u}
!>              - \mu \gradt \vect{a}^{n+1}\right)
!> = \vect{Rhs}
!> \f]
!>
!> This equation is rewritten as:
!>
!> \f[
!> \tens{f_s}^{imp} \delta \vect{a}
!> + \divv \left( \delta \vect{a} \otimes \rho \vect{u}
!>              - \mu \gradt \delta \vect{a} \right)
!> = \vect{Rhs}^1
!> \f]
!>
!> where \f$ \delta \vect{a} = \vect{a}^{n+1} - \vect{a}^n\f$ and
!> \f$ \vect{Rhs}^1 = \vect{Rhs}
!> - \divv \left( \vect{a}^n \otimes \rho \vect{u}
!>              - \mu \gradt \vect{a}^n \right)\f$
!>
!>
!> It is in fact solved with the following iterative process:
!>
!> \f[
!> \tens{f_s}^{imp} \delta \vect{a}^k
!> + \divv \left( \delta \vect{a}^k \otimes \rho \vect{u}
!>              - \mu \gradt \delta \vect{a}^k \right)
!> = \vect{Rhs}^k
!> \f]
!>
!> where \f$ \vect{Rhs}^k = \vect{Rhs}
!> - \tens{f_s}^{imp} \left(\vect{a}^k-\vect{a}^n \right)
!> - \divv \left( \vect{a}^k \otimes \rho \vect{u}
!>              - \mu \gradt \vect{a}^k \right)\f$
!>
!> Be careful, it is forbidden to modify \f$ \tens{f_s}^{imp} \f$ here!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     idtvar        indicator of the temporal scheme
!> \param[in]     ivar          index of the current variable
!> \param[in]     iconvp        indicator
!>                               - 1 convection,
!>                               - 0 otherwise
!> \param[in]     idiffp        indicator
!>                               - 1 diffusion,
!>                               - 0 otherwise
!> \param[in]     ndircp        indicator (0 if the diagonal is stepped aside)
!> \param[in]     imrgra        indicateur
!>                               - 0 iterative gradient
!>                               - 1 least squares gradient
!> \param[in]     nswrsp        number of reconstruction sweeps for the
!>                               Right Hand Side
!> \param[in]     nswrgp        number of reconstruction sweeps for the
!>                               gradients
!> \param[in]     imligp        clipping gradient method
!>                               - < 0 no clipping
!>                               - = 0 thank to neighboring gradients
!>                               - = 1 thank to the mean gradient
!> \param[in]     ircflp        indicator
!>                               - 1 flux reconstruction,
!>                               - 0 otherwise
!> \param[in]     ivisep        indicator to take \f$ \divv
!>                               \left(\mu \gradt \transpose{\vect{a}} \right)
!>                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
!>                               - 1 take into account,
!>                               - 0 otherwise
!> \param[in]     ischcp        indicator
!>                               - 1 centered
!>                               - 0 2nd order
!> \param[in]     isstpp        indicator
!>                               - 1 without slope test
!>                               - 0 with slope test
!> \param[in]     iescap        compute the predictor indicator if 1
!> \param[in]     idftnp        indicator
!>                               - 1 the diffusivity is scalar
!>                               - 6 the diffusivity is a symmetric tensor
!> \param[in]     iswdyp        indicator
!>                               - 0 no dynamic relaxation
!>                               - 1 dynamic relaxation depending on
!>                                 \f$ \delta \vect{\varia}^k \f$
!>                               - 2 dynamic relaxation depending on
!>                                 \f$ \delta \vect{\varia}^k \f$  and
!>                                 \f$ \delta \vect{\varia}^{k-1} \f$
!> \param[in]     iwarnp        verbosity
!> \param[in]     blencp        fraction of upwinding
!> \param[in]     epsilp        precision pour resol iter
!> \param[in]     epsrsp        relative precision for the iterative process
!> \param[in]     epsrgp        relative precision for the gradient
!>                               reconstruction
!> \param[in]     climgp        clipping coefficient for the computation of
!>                               the gradient
!> \param[in]     relaxp        coefficient of relaxation
!> \param[in]     thetap        weighting coefficient for the theta-schema,
!>                               - thetap = 0: explicit scheme
!>                               - thetap = 0.5: time-centered
!>                               scheme (mix between Crank-Nicolson and
!>                               Adams-Bashforth)
!>                               - thetap = 1: implicit scheme
!> \param[in]     pvara         variable at the previous time step
!>                               \f$ \vect{a}^n \f$
!> \param[in]     pvark         variable at the previous sub-iteration
!>                               \f$ \vect{a}^k \f$.
!>                               If you sub-iter on Navier-Stokes, then
!>                               it allows to initialize by something else than
!>                               pvara (usually pvar=pvara)
!> \param[in]     coefav        boundary condition array for the variable
!>                               (explicit part)
!> \param[in]     coefbv        boundary condition array for the variable
!>                               (implicit part)
!> \param[in]     cofafv        boundary condition array for the diffusion
!>                               of the variable (Explicit part)
!> \param[in]     cofbfv        boundary condition array for the diffusion
!>                               of the variable (Implicit part)
!> \param[in]     flumas        mass flux at interior faces
!> \param[in]     flumab        mass flux at boundary faces
!> \param[in]     viscfm        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the matrix
!> \param[in]     viscbm        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at boundary faces for the matrix
!> \param[in]     viscfs        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the r.h.s.
!> \param[in]     viscbs        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at boundary faces for the r.h.s.
!> \param[in]     secvif        secondary viscosity at interior faces
!> \param[in]     secvib        secondary viscosity at boundary faces
!> \param[in]     icvflb        global indicator of boundary convection flux
!>                               - 0 upwind scheme at all boundary faces
!>                               - 1 imposed flux at some boundary faces
!> \param[in]     icvfli        boundary face indicator array of convection flux
!>                               - 0 upwind scheme
!>                               - 1 imposed flux
!> \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$
!> \param[in]     smbrp         Right hand side \f$ \vect{Rhs}^k \f$
!> \param[in,out] pvar          current variable
!> \param[out]    eswork        prediction-stage error estimator
!>                              (if iescap > 0)
!_______________________________________________________________________________

subroutine coditv &
 ( idtvar , ivar   , iconvp , idiffp , ndircp ,                   &
   imrgra , nswrsp , nswrgp , imligp , ircflp , ivisep ,          &
   ischcp , isstpp , iescap , idftnp , iswdyp ,                   &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp ,                   &
   relaxp , thetap ,                                              &
   pvara  , pvark  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab ,                                              &
   viscfm , viscbm , viscfs , viscbs , secvif , secvib ,          &
   icvflb , icvfli ,                                              &
   fimp   ,                                                       &
   smbrp  ,                                                       &
   pvar   ,                                                       &
   eswork )

!===============================================================================
! Module files
!===============================================================================

use dimens, only: ndimfb
use paramx
use numvar
use cstnum
use entsor
use parall
use period
use mesh
use field
use cs_f_interfaces
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          idtvar , ivar   , iconvp , idiffp , ndircp
integer          imrgra , nswrsp , nswrgp , imligp , ircflp
integer          ischcp , isstpp , iescap
integer          idftnp , iswdyp , icvflb
integer          iwarnp
integer          ivisep

integer          icvfli(nfabor)

double precision blencp , epsilp , epsrgp , climgp
double precision relaxp , thetap , epsrsp

double precision pvara(3,ncelet)
double precision pvark(3,ncelet)
double precision pvar(3,ncelet)
double precision coefav(3,ndimfb)
double precision cofafv(3,ndimfb)
double precision coefbv(3,3,ndimfb)
double precision cofbfv(3,3,ndimfb)
double precision flumas(nfac), flumab(nfabor)
double precision viscfm(*), viscbm(nfabor)
double precision viscfs(*), viscbs(nfabor)
double precision secvif(nfac), secvib(nfabor)
double precision fimp(3,3,ncelet)
double precision smbrp(3,ncelet)
double precision eswork(3,ncelet)

! Local variables

character(len=80) :: chaine
character(len=16) :: cnom
integer          f_id,isym,isqrt
integer          inc,isweep,niterf,iel,nswmod
integer          iinvpe
integer          idtva0
integer          isou , jsou
integer          ibsize, iesize
integer          lvar, insqrt

double precision residu, rnorm, ressol, rnorm2
double precision thetex
double precision alph, beta
double precision paxkrk, nadxk, paxm1rk, nadxkm1, paxm1ax

type(solving_info) sinfo

double precision, allocatable, dimension(:,:,:) :: dam
double precision, allocatable, dimension(:,:) :: xam
double precision, allocatable, dimension(:,:) :: dpvar, smbini, w1
double precision, allocatable, dimension(:,:) :: adxk, adxkm1, dpvarm1, rhs0

!===============================================================================

!===============================================================================
! 0.  Initialization
!===============================================================================

! Matrix block size
ibsize = 3
if (idftnp.eq.1) iesize = 1
if (idftnp.eq.6) iesize = 3

! Allocate temporary arrays
allocate(dam(3,3,ncelet))
allocate(dpvar(3,ncelet), smbini(3,ncelet))
if (iswdyp.ge.1) then
  allocate(adxk(3,ncelet), adxkm1(3,ncelet), dpvarm1(3,ncelet))
  allocate(rhs0(3,ncelet))
endif

! Name
if (ivar.gt.0) then
  f_id = ivarfl(ivar)
  call field_get_name(f_id, chaine)
  call field_get_key_struct_solving_info(f_id, sinfo)
else
  f_id = -1
  chaine = nomva0
endif
cnom= chaine(1:16)

! Symmetric matrix, except if advection
isym  = 1
if (iconvp.gt.0) isym  = 2

! be carefull here, xam is interleaved
if (iesize.eq.1) allocate(xam(isym,nfac))
if (iesize.eq.3) allocate(xam(3*3*isym,nfac))

! PRISE DE SQRT DANS PS
isqrt = 1

! iinvpe is useless in the vectorial framework
iinvpe = 0

!===============================================================================
! 1.  Building of the "simplified" matrix
!===============================================================================

if (iesize.eq.1) then

  call matrxv &
  !==========
   ( iconvp , idiffp , ndircp , isym   ,                            &
     thetap ,                                                       &
     coefbv , cofbfv , fimp   ,                                     &
     flumas , flumab , viscfm , viscbm ,                            &
     dam    , xam    )

elseif (iesize.eq.3) then

  call matrvv &
  !==========
   ( iconvp , idiffp , ndircp , isym   ,                            &
     thetap ,                                                       &
     coefbv , cofbfv , fimp   ,                                     &
     flumas , flumab , viscfm , viscbm ,                            &
     dam    , xam    )

endif

! For steady computations, the diagonal is relaxed
if (idtvar.lt.0) then
  !$omp parallel do private(isou, jsou)
  do iel = 1, ncel
    do isou = 1, 3
      do jsou = 1, 3
        dam(isou,jsou,iel) = dam(isou,jsou,iel)/relaxp
      enddo
    enddo
  enddo
endif

!===============================================================================
! 3. Iterative process to handle non orthogonlaities (starting from the second
! iteration).
!===============================================================================

! Application du theta schema

! On calcule le bilan explicite total
thetex = 1.d0 - thetap


! Si THETEX=0, ce n'est pas la peine d'en rajouter
if (abs(thetex).gt.epzero) then
  inc    = 1

  call bilscv &
  !==========
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , ivisep ,                   &
   iwarnp , idftnp ,                                              &
   blencp , epsrgp , climgp , relaxp , thetex ,                   &
   pvara  , pvara  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab , viscfs , viscbs , secvif , secvib ,          &
   icvflb , icvfli ,                                              &
   smbrp  )
endif

! Before looping, the RHS without reconstruction is stored in smbini

!$omp parallel do private(isou)
do iel = 1, ncel
  do isou = 1, 3
    smbini(isou,iel) = smbrp(isou,iel)
  enddo
enddo

! pvar is initialized on ncelet to avoid a synchronization

!$omp parallel do private(isou)
do iel = 1, ncelet
  do isou = 1, 3
    pvar(isou,iel) = pvark(isou,iel)
  enddo
enddo

! In the following, bilscv is called with inc=1,
! except for Weight Matrix (nswrsp=-1)
inc = 1
if (nswrsp.eq.-1) then
  nswrsp = 1
  inc = 0
endif

! ---> INCREMENTATION ET RECONSTRUCTION DU SECOND MEMBRE

!  On est entre avec un smb explicite base sur PVARA.
!     si on initialise avec PVAR avec autre chose que PVARA
!     on doit donc corriger SMBR (c'est le cas lorsqu'on itere sur navsto)

!$omp parallel do private(isou)
do iel = 1, ncel
  do isou = 1, 3
    smbrp(isou,iel) = 0.d0
  enddo
enddo

call bilscv &
!==========
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , ivisep ,                   &
   iwarnp , idftnp ,                                              &
   blencp , epsrgp , climgp , relaxp , thetap ,                   &
   pvar   , pvara  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab , viscfs , viscbs , secvif , secvib ,          &
   icvflb , icvfli ,                                              &
   smbrp  )

! Dynamic relaxation
if (iswdyp.ge.1) then

  !$omp parallel do private(isou)
  do iel = 1, ncel
    do isou = 1, 3
      rhs0(isou,iel) = smbrp(isou,iel)
      smbini(isou,iel) = smbini(isou,iel)                        &
                 -fimp(isou,1,iel)*(pvar(1,iel) - pvara(1,iel))  &
                 -fimp(isou,2,iel)*(pvar(2,iel) - pvara(2,iel))  &
                 -fimp(isou,3,iel)*(pvar(3,iel) - pvara(3,iel))
      smbrp(isou,iel) = smbrp(isou,iel) + smbini(isou,iel)

      adxkm1(isou,iel) = 0.d0
      adxk(isou,iel) = 0.d0
      dpvar(isou,iel) = 0.d0
    enddo
  enddo

  ! ||A.dx^0||^2 = 0
  nadxk = 0.d0

else

  !$omp parallel do private(isou)
  do iel = 1, ncel
    do isou = 1, 3
      smbini(isou,iel) = smbini(isou,iel)                        &
                 -fimp(isou,1,iel)*(pvar(1,iel) - pvara(1,iel))  &
                 -fimp(isou,2,iel)*(pvar(2,iel) - pvara(2,iel))  &
                 -fimp(isou,3,iel)*(pvar(3,iel) - pvara(3,iel))
      smbrp(isou,iel) = smbrp(isou,iel) + smbini(isou,iel)
    enddo
  enddo
endif

! --- Convergence test
call prodsc(3*ncel, isqrt, smbrp, smbrp, residu)

! ---> RESIDU DE NORMALISATION
!    (NORME C.L +TERMES SOURCES+ TERMES DE NON ORTHOGONALITE)

allocate(w1(3, ncelet))  ! Allocate a temporary array

call promav(isym, ibsize, iesize, iinvpe, dam, xam, pvar, w1)

!$omp parallel do private(isou)
do iel = 1, ncel
   do isou = 1, 3
      w1(isou,iel) = w1(isou,iel) + smbrp(isou,iel)
   enddo
enddo

call prodsc(3*ncel, isqrt, w1, w1, rnorm)

sinfo%rnsmbr = rnorm
rnorm2 = rnorm**2

deallocate(w1)  ! Free memory

! Warning: for Weight Matrix, one and only one sweep is done.
nswmod = max(nswrsp, 1)

isweep = 1

! Reconstruction loop (beginning)
!--------------------------------
sinfo%nbivar = 0

do while (isweep.le.nswmod.and.residu.gt.epsrsp*rnorm.or.isweep.eq.1)

  ! --- Solving on the increment dpvar

  ! Dynamic relaxation of the system
  if (iswdyp.ge.1) then

    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou = 1, 3
        dpvarm1(isou,iel) = dpvar(isou,iel)
        dpvar(isou,iel) = 0.d0
      enddo
    enddo
  else

    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou =1,3
        dpvar(isou,iel) = 0.d0
      enddo
    enddo
  endif

  ! iinvpe is useless in the vectorial framework
  iinvpe = 0

  ! Matrix block size
  ibsize = 3

  ! Solver residual
  ressol = residu

  call sles_solve_native(f_id, chaine,                                 &
                         isym, ibsize, iesize, dam, xam, iinvpe,       &
                         epsilp, rnorm, niterf, ressol, smbrp, dpvar)

  ! Dynamic relaxation of the system
  if (iswdyp.ge.1) then

    ! Computation of the variable relaxation coefficient
    lvar = 0

    !$omp parallel do private(isou)
    do iel = 1, ncelet
      do isou = 1, 3
        adxkm1(isou,iel) = adxk(isou,iel)
        adxk(isou,iel) = - rhs0(isou,iel)
      enddo
    enddo

    call bilscv &
    !==========
   ( idtvar , lvar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
     ischcp , isstpp , inc    , imrgra , ivisep ,                   &
     iwarnp , idftnp ,                                              &
     blencp , epsrgp , climgp , relaxp , thetap ,                   &
     dpvar  , dpvar  ,                                              &
     coefav , coefbv , cofafv , cofbfv ,                            &
     flumas , flumab , viscfs , viscbs , secvif , secvib ,          &
     icvflb , icvfli ,                                              &
     adxk   )

    insqrt = 0

    ! ||E.dx^(k-1)-E.0||^2
    nadxkm1 = nadxk

    ! ||E.dx^k-E.0||^2
    call prodsc(3*ncel , insqrt , adxk , adxk , nadxk)

    ! < E.dx^k-E.0; r^k >
    call prodsc(3*ncel , insqrt , smbrp , adxk , paxkrk)

    ! Relaxation with respect to dx^k and dx^(k-1)
    if (iswdyp.ge.2) then

      ! < E.dx^(k-1)-E.0; r^k >
      call prodsc(3*ncel , insqrt , smbrp , adxkm1 , paxm1rk)

      ! < E.dx^(k-1)-E.0; E.dx^k-E.0 >
      call prodsc(3*ncel , insqrt , adxk, adxkm1 , paxm1ax)

      if (nadxkm1.gt.1.d-30*rnorm2.and.                     &
          (nadxk*nadxkm1-paxm1ax**2).gt.1.d-30*rnorm2) then
        beta = (paxkrk*paxm1ax - nadxk*paxm1rk)/(nadxk*nadxkm1-paxm1ax**2)
      else
        beta = 0.d0
      endif

    else
      beta = 0.d0
      paxm1ax =1.d0
      paxm1rk = 0.d0
      paxm1ax = 0.d0
    endif

    ! The first sweep is not relaxed
    if (isweep.eq.1) then
      alph = 1.d0
      beta = 0.d0
    elseif (isweep.eq.2) then
      beta = 0.d0
      alph = -paxkrk/max(nadxk, 1.d-30*rnorm2)
    else
      alph = -(paxkrk + beta*paxm1ax)/max(nadxk, 1.d-30*rnorm2)
    endif

    ! Writing
    if (iwarnp.ge.3) then
      write(nfecra,1200) cnom, isweep, alph, beta, &
                         paxkrk, nadxk, paxm1rk, nadxkm1, paxm1ax
    endif

  endif

  ! --- Update the solution with the increment

  if (iswdyp.eq.0) then

    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou = 1, 3
         pvar(isou,iel) = pvar(isou,iel) + dpvar(isou,iel)
      enddo
    enddo
  elseif (iswdyp.eq.1) then

    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou = 1, 3
         pvar(isou,iel) = pvar(isou,iel) + alph*dpvar(isou,iel)
      enddo
    enddo
  elseif (iswdyp.ge.2) then

    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou = 1, 3
         pvar(isou,iel) = pvar(isou,iel) + alph*dpvar(isou,iel)  &
                        + beta*dpvarm1(isou,iel)
      enddo
    enddo
  endif

  ! ---> Handle parallelism and periodicity

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvin(pvar)
  endif

  ! --- Update the right hand and compute the new residual

  if (iswdyp.eq.0) then

    !$omp parallel do private(isou)
    do iel = 1, ncel
      ! smbini already contains unsteady terms and mass source terms
      ! of the RHS updated at each sweep
      do isou = 1, 3
        smbini(isou,iel) = smbini(isou,iel)                 &
                  - fimp(isou,1,iel)*dpvar(1,iel)           &
                  - fimp(isou,2,iel)*dpvar(2,iel)           &
                  - fimp(isou,3,iel)*dpvar(3,iel)
        smbrp(isou,iel) = smbini(isou,iel)
      enddo
    enddo

  elseif (iswdyp.eq.1) then

    !$omp parallel do private(isou)
    do iel = 1, ncel
      ! smbini already contains unsteady terms and mass source terms
      ! of the RHS updated at each sweep
      do isou = 1, 3
        smbini(isou,iel) = smbini(isou,iel)                 &
                  - fimp(isou,1,iel)*alph*dpvar(1,iel)      &
                  - fimp(isou,2,iel)*alph*dpvar(2,iel)      &
                  - fimp(isou,3,iel)*alph*dpvar(3,iel)
        smbrp(isou,iel) = smbini(isou,iel)
      enddo
    enddo

  elseif (iswdyp.eq.2) then

    !$omp parallel do private(isou)
    do iel = 1, ncel
      ! smbini already contains unsteady terms and mass source terms
      ! of the RHS updated at each sweep
      do isou = 1, 3
        smbini(isou,iel) = smbini(isou,iel)                                  &
                  - fimp(isou,1,iel)*(alph*dpvar(1,iel)+beta*dpvarm1(1,iel)) &
                  - fimp(isou,2,iel)*(alph*dpvar(2,iel)+beta*dpvarm1(2,iel)) &
                  - fimp(isou,3,iel)*(alph*dpvar(3,iel)+beta*dpvarm1(3,iel))
        smbrp(isou,iel) = smbini(isou,iel)
      enddo
    enddo
  endif

  call bilscv &
  !==========
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , ivisep ,                   &
   iwarnp , idftnp ,                                              &
   blencp , epsrgp , climgp , relaxp , thetap ,                   &
   pvar   , pvara  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab , viscfs , viscbs , secvif , secvib ,          &
   icvflb , icvfli ,                                              &
   smbrp  )

  ! --- Convergence test
  call prodsc(3*ncel, isqrt, smbrp, smbrp, residu)

  ! Writing
  sinfo%nbivar = sinfo%nbivar + niterf
  ! Writing
  if (iwarnp.ge.3) then
     write(nfecra,1000) cnom, isweep, residu, rnorm
  endif

  isweep = isweep + 1

enddo
! --- Reconstruction loop (end)

! Writing: convergence
if (abs(rnorm)/sqrt(3.d0).gt.epzero) then
  sinfo%resvar = residu/rnorm
else
  sinfo%resvar = 0.d0
endif

if (iwarnp.ge.1) then
  if (residu.le.epsrsp*rnorm) then
    write(nfecra,1000) cnom,isweep-1,residu,rnorm

! Writing: non-convergence
  else if (isweep.gt.nswmod) then
    write(nfecra,1100) cnom, nswmod
  endif
endif

! Save convergence info
if (ivar.gt.0) then
  call field_set_key_struct_solving_info(ivarfl(ivar), sinfo)
endif


!===============================================================================
! 4. After having computed the new value, an estimator is computed for the
! prediction step of the velocity.
!===============================================================================

if (iescap.gt.0) then

  ! ---> Computation of the estimator of the current component

  ! smbini already contains unsteady terms and mass source terms
  ! of the RHS updated at each sweep

  !$omp parallel do private(isou)
  do iel = 1, ncel
    do isou = 1, 3
      smbrp(isou,iel) = smbini(isou,iel) - fimp(isou,1,iel)*dpvar(1,iel) &
                                         - fimp(isou,2,iel)*dpvar(2,iel) &
                                         - fimp(isou,3,iel)*dpvar(3,iel)
    enddo
  enddo

  inc    = 1
  ! Without relaxation even for a statonnary computation
  idtva0 = 0

  call bilscv &
  !==========
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , ivisep ,                   &
   iwarnp , idftnp ,                                              &
   blencp , epsrgp , climgp , relaxp , thetap ,                   &
   pvar   , pvara  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab , viscfs , viscbs , secvif , secvib ,          &
   icvflb , icvfli ,                                              &
   smbrp   )

  ! Contribution of the current component to the L2 norm stored in eswork

  !$omp parallel do private(isou)
  do iel = 1, ncel
    do isou = 1, 3
      eswork(isou,iel) = (smbrp(isou,iel)/ volume(iel))**2
    enddo
  enddo

endif

!===============================================================================
! 5. Free solver setup
!===============================================================================

call sles_free_native(f_id, chaine)

! Free memory
deallocate(dam, xam)
deallocate(dpvar, smbini)
if (iswdyp.ge.1) deallocate(adxk, adxkm1, dpvarm1, rhs0)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format ( &
 1X,A16,' : CV-DIF-TS',I5,' IT - RES= ',E12.5,' NORME= ', E12.5)
 1100 format ( &
'@'                                                                 ,/,&
'@ @@ ATTENTION : ',A8 ,' CONVECTION-DIFFUSION-TERMES SOURCES'      ,/,&
'@    ========='                                                    ,/,&
'@  Nombre d''iterations maximal ',I10   ,' atteint'                ,/,&
'@' )
 1200 format ( &
 1X,A16,' Sweep: ',I5,' Dynamic relaxation: alpha = ',E12.5,' beta = ',E12.5,/,&
'    < dI^k  ; R^k > = ',E12.5,' ||dI^k  ||^2 = ',E12.5                     ,/,&
'    < dI^k-1; R^k > = ',E12.5,' ||dI^k-1||^2 = ',E12.5                     ,/,&
'   < dI^k-1; dI^k > = ',E12.5)

#else

 1000 format ( &
 1X,A16,' : CV-DIF-TS',I5,' IT - RES= ',E12.5,' NORM= ', E12.5)
 1100 format ( &
'@'                                                                 ,/,&
'@ @@ WARNING: ',A8 ,' CONVECTION-DIFFUSION-SOURCE TERMS'           ,/,&
'@    ========'                                                     ,/,&
'@  Maximum number of iterations ',I10   ,' reached'                ,/,&
'@' )
 1200 format ( &
 1X,A16,' Sweep: ',I5,' Dynamic relaxation: alpha = ',E12.5,' beta = ',E12.5,/,&
'    < dI^k  ; R^k > = ',E12.5,' ||dI^k  ||^2 = ',E12.5                     ,/,&
'    < dI^k-1; R^k > = ',E12.5,' ||dI^k-1||^2 = ',E12.5                     ,/,&
'   < dI^k-1; dI^k > = ',E12.5)

#endif

!----
! End
!----

end subroutine
