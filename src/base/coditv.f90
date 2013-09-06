!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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
!>                               - 0 sinon
!> \param[in]     idiffp        indicator
!>                               - 1 diffusion,
!>                               - 0 sinon
!> \param[in]     ireslp        indicator
!>                               - 0 conjugate gradient
!>                               - 1 jacobi
!>                               - 2 bi-cgstab
!> \param[in]     ndircp        indicator (0 if the diagonal is stepped aside)
!> \param[in]     nitmap        maximum number of iteration to solve
!>                               the iterative process
!> \param[in]     imrgra        indicateur
!>                               - 0 iterative gradient
!>                               - 1 least square gradient
!> \param[in]     nswrsp        number of reconstruction sweeps for the
!>                               Right Hand Side
!> \param[in]     nswrgp        number of reconstruction sweeps for the
!>                               gradients
!> \param[in]     imligp        clipping gradient method
!>                               - < 0 no clipping
!>                               - = 0 thank to neighbooring gradients
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
!>                               - 1 centred
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
!> \param[in]     imgrp         indicator
!>                               - 0 no multi-grid
!>                               - 1 otherwise
!> \param[in]     ncymxp        max. number of multigrid cycles
!> \param[in]     nitmfp        number of equivalent iterations on fine mesh
!> \param[in]     ippu          index of the variable for post-processing
!> \param[in]     ippv          index of the variable for post-processing
!> \param[in]     ippw          index of the variable for post-processing
!> \param[in]     iwarnp        verbosity
!> \param[in]     blencp        fraction of upwinding
!> \param[in]     epsilp        precision pour resol iter
!> \param[in]     epsrsp        relative precision for the iterative process
!> \param[in]     epsrgp        relative precision for the gradient
!>                               reconstruction
!> \param[in]     climgp        clipping coeffecient for the computation of
!>                               the gradient
!> \param[in]     extrap        coefficient for extrapolation of the gradient
!> \param[in]     relaxp        coefficient of relaxation
!> \param[in]     thetap        weightening coefficient for the theta-schema,
!>                               - thetap = 0: explicit scheme
!>                               - thetap = 0.5: time-centred
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
!>                               (Explicit part)
!> \param[in]     coefbv        boundary condition array for the variable
!>                               (Impplicit part)
!> \param[in]     cofafv        boundary condition array for the diffusion
!>                               of the variable (Explicit part)
!> \param[in]     cofbfv        boundary condition array for the diffusion
!>                               of the variable (Implicit part)
!> \param[in]     flumas        mass flux at interior faces
!> \param[in]     flumab        mass flux at boundary faces
!> \param[in]     viscfm        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the matrix
!> \param[in]     viscbm        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at border faces for the matrix
!> \param[in]     viscfs        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the r.h.s.
!> \param[in]     viscbs        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at border faces for the r.h.s.
!> \param[in]     icvflb
!> \param[in]     icvfbr
!> \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$
!> \param[in]     smbrp         Right hand side \f$ \vect{Rhs}^k \f$
!> \param[in,out] pvar          current variable
!> \param[out]    eswork        prediction-stage error estimator
!>                              (if iescap > 0)
!_______________________________________________________________________________

subroutine coditv &
 ( idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp , ivisep ,          &
   ischcp , isstpp , iescap , idftnp , iswdyp ,                   &
   imgrp  , ncymxp , nitmfp , ippu   , ippv   , ippw   , iwarnp , &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   pvara  , pvark  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab ,                                              &
   viscfm , viscbm , viscfs , viscbs , secvif , secvib ,          &
   icvflb , icvfbr ,                                              &
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
use mltgrd
use optcal, only: rlxp1
use mesh

!===============================================================================

implicit none

! Arguments

integer          idtvar , ivar   , iconvp , idiffp , ndircp
integer          nitmap
integer          imrgra , nswrsp , nswrgp , imligp , ircflp
integer          ischcp , isstpp , iescap , imgrp
integer          idftnp , iswdyp , icvflb
integer          ncymxp , nitmfp
integer          iwarnp
integer          ippu   , ippv   , ippw   , ivisep

integer          icvfbr(nfabor)

double precision blencp , epsilp , epsrgp , climgp , extrap
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

character*80     chaine
character*16     cnom(3)
integer          lchain
integer          isym,ireslp,ireslq,ipol,isqrt
integer          inc,isweep,niterf,iel,icycle,nswmod
integer          iinvpe
integer          idtva0
integer          nagmax, npstmg
integer          isou , jsou
integer          ibsize, iesize
integer          incp, insqrt

double precision residu, rnorm, ressol, rnorm2
double precision thetex
double precision alph, beta
double precision paxkrk, nadxk, paxm1rk, nadxkm1, paxm1ax

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
! be carefull here, xam is interleaved
if (iesize.eq.1) allocate(xam(2,nfac))
if (iesize.eq.3) allocate(xam(3*3*2,nfac))
allocate(dpvar(3,ncelet), smbini(3,ncelet))
if (iswdyp.ge.1) then
  allocate(adxk(3,ncelet), adxkm1(3,ncelet), dpvarm1(3,ncelet))
  allocate(rhs0(3,ncelet))
endif

! Names
chaine = nomvar(ippu)
cnom(1)= chaine(1:16)
chaine = nomvar(ippv)
cnom(2)= chaine(1:16)
chaine = nomvar(ippw)
cnom(3)= chaine(1:16)

! Symmetric matrix, except if advection
isym  = 1
if (iconvp.gt.0) isym  = 2

! METHODE DE RESOLUTION ET DEGRE DU PRECOND DE NEUMANN
!     0 SI CHOIX AUTOMATIQUE GRADCO OU BICGSTAB
!     0 SI CHOIX AUTOMATIQUE JACOBI
!     DONNE PAR IRESLP/1000 SI NON AUTOMATIQUE
if (ireslp.eq.-1) then
  ireslq = 0
  ipol   = 0
  if (iconvp.gt.0) then
    ireslq = 1
    ipol   = 0
  endif
else
  ireslq = mod(ireslp+10000,1000)
  ipol   = (ireslp-ireslq)/1000
endif

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
   ( ncelet , ncel   , nfac   , nfabor ,                            &
     iconvp , idiffp , ndircp , isym   , nfecra ,                   &
     thetap ,                                                       &
     ifacel , ifabor ,                                              &
     coefbv , cofbfv , fimp   ,                                     &
     flumas , flumab , viscfm , viscbm ,                            &
     dam    , xam    )

elseif (iesize.eq.3) then

  call matrvv &
  !==========
   ( ncelet , ncel   , nfac   , nfabor ,                            &
     iconvp , idiffp , ndircp , isym   , nfecra ,                   &
     thetap ,                                                       &
     ifacel , ifabor ,                                              &
     coefbv , cofbfv , fimp   ,                                     &
     flumas , flumab , viscfm , viscbm ,                            &
     dam    , xam    )

endif

! For stationary computations, the diagonal is relaxed
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
! 2. Preparation of the Algebraic Multigrid
!===============================================================================

if (imgrp.gt.0) then

  ! --- Building of the mesh hierarchy

  chaine = nomvar(ippu)
  iwarnp = iwarni(ivar)
  nagmax = nagmx0(ivar)
  npstmg = ncpmgr(ivar)
  lchain = 16

  call clmlga &
  !==========
 ( chaine(1:16) ,    lchain ,                                     &
   isym   , ibsize , iesize , nagmax , npstmg , iwarnp ,          &
   ngrmax , ncegrm ,                                              &
   rlxp1  ,                                                       &
   dam    , xam    )

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
   ippu   , ippv   , ippw   , iwarnp , idftnp ,                   &
   blencp , epsrgp , climgp , extrap , relaxp , thetex ,          &
   pvar   , pvara  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab , viscfs , viscbs , secvif , secvib ,          &
   icvflb , icvfbr ,                                              &
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
   ippu   , ippv   , ippw   , iwarnp , idftnp ,                   &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   pvar   , pvara  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab , viscfs , viscbs , secvif , secvib ,          &
   icvflb , icvfbr ,                                              &
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

rnsmbr(ippu) = rnorm
rnsmbr(ippv) = rnorm
rnsmbr(ippw) = rnorm
rnorm2 = rnorm**2

deallocate(w1)  ! Free memory

! Warning: for Weight Matrix, one and only one sweep is done.
nswmod = max(nswrsp, 1)

isweep = 1

! Reconstruction loop (beginning)
!--------------------------------
nbivar(ippu) = 0
nbivar(ippv) = 0
nbivar(ippw) = 0

do while (isweep.le.nswmod.and.residu.gt.epsrsp*rnorm.or.isweep.eq.1)

  ! --- Solving on the increment dpvar

  ! Dynamic relaxation of the system
  if (iswdyp.ge.1) then

    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou =1,3
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

  call invers &
  !==========
 ( cnom(1), isym   , ibsize , iesize ,                            &
   ipol   , ireslq , nitmap , imgrp  ,                            &
   ncymxp , nitmfp ,                                              &
   iwarnp , nfecra , niterf , icycle , iinvpe ,                   &
   epsilp , rnorm  , ressol ,                                     &
   dam    , xam    , smbrp  , dpvar  )

  ! Dynamic relaxation of the system
  if (iswdyp.ge.1) then

    ! Computation of the variable ralaxation coefficient

    !$omp parallel do private(isou)
    do iel = 1, ncelet
      do isou = 1, 3
        adxkm1(isou,iel) = adxk(isou,iel)
        adxk(isou,iel) = - rhs0(isou,iel)
      enddo
    enddo

    call bilscv &
    !==========
   ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
     ischcp , isstpp , incp   , imrgra , ivisep ,                   &
     ippu   , ippv   , ippw   , iwarnp , idftnp ,                   &
     blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
     dpvar  , dpvar  ,                                              &
     coefav , coefbv , cofafv , cofbfv ,                            &
     flumas , flumab , viscfs , viscbs , secvif , secvib ,          &
     icvflb , icvfbr ,                                              &
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
      ! smbini already contains instationnary terms and mass source terms
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
      ! smbini already contains instationnary terms and mass source terms
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
      ! smbini already contains instationnary terms and mass source terms
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
   ippu   , ippv   , ippw   , iwarnp , idftnp ,                   &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   pvar   , pvara  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab , viscfs , viscbs , secvif , secvib ,          &
   icvflb , icvfbr ,                                              &
   smbrp  )

  ! --- Convergence test
  call prodsc(3*ncel, isqrt, smbrp, smbrp, residu)

  ! Writing
  nbivar(ippu) = nbivar(ippu) + niterf
  nbivar(ippv) = nbivar(ippv) + niterf
  nbivar(ippw) = nbivar(ippw) + niterf
  ! Writing
  if (iwarnp.ge.3) then
     write(nfecra,1000) cnom(1), isweep, residu, rnorm
     write(nfecra,1000) cnom(2), isweep, residu, rnorm
     write(nfecra,1000) cnom(3), isweep, residu, rnorm
  endif

  isweep = isweep + 1

enddo
! --- Reconstruction loop (end)

! Writing: convergence
if (abs(rnorm)/sqrt(3.d0).gt.epzero) then
  resvar(ippu) = residu/rnorm
  resvar(ippv) = residu/rnorm
  resvar(ippw) = residu/rnorm
else
  resvar(ippu) = 0.d0
  resvar(ippv) = 0.d0
  resvar(ippw) = 0.d0
endif

if (iwarnp.ge.1) then
  if (residu.le.epsrsp*rnorm) then
    write(nfecra,1000) cnom(1),isweep-1,residu,rnorm
    write(nfecra,1000) cnom(2),isweep-1,residu,rnorm
    write(nfecra,1000) cnom(3),isweep-1,residu,rnorm
  endif
endif

! Writing: non-convergence
if (iwarnp.ge.1) then
  if (isweep.gt.nswmod) then
    write(nfecra,1100) cnom(1), nswmod
    write(nfecra,1100) cnom(2), nswmod
    write(nfecra,1100) cnom(3), nswmod
  endif
endif

!===============================================================================
! 4. After having computed the new value, an estimator is computed for the
! prediction step of the velocity.
!===============================================================================

if (iescap.gt.0) then

  ! ---> Computation of the estimator of the current component

  ! smbini already contains instationnary terms and mass source terms
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
   ippu   , ippv   , ippw   , iwarnp , idftnp ,                   &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   pvar   , pvara  ,                                              &
   coefav , coefbv , cofafv , cofbfv ,                            &
   flumas , flumab , viscfs , viscbs , secvif , secvib ,          &
   icvflb , icvfbr ,                                              &
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
! 5. Suppression of the mesh hierarchy
!===============================================================================

if (imgrp.gt.0) then
  chaine = nomvar(ippu)
  lchain = 16
  call dsmlga(chaine(1:16), lchain)
endif

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
