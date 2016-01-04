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

!> \file coditts.f90
!>
!> \brief This function solves an advection diffusion equation with source terms
!> for one time step for the symmetric tensor variable \f$ \tens{\variat} \f$.
!>
!> The equation reads:
!>
!> \f[
!> \tens{f_s}^{imp}(\tens{\variat}^{n+1}-\tens{\variat}^n)
!> + \divt \left( \tens{\variat}^{n+1} \otimes \rho \vect {u}
!>              - \mu \gradtt \tens{\variat}^{n+1}\right)
!> = \tens{Rhs}
!> \f]
!>
!> This equation is rewritten as:
!>
!> \f[
!> \tens{f_s}^{imp} \delta \tens{\variat}
!> + \divt \left( \delta \tens{\variat} \otimes \rho \vect{u}
!>              - \mu \gradtt \delta \tens{\variat} \right)
!> = \tens{Rhs}^1
!> \f]
!>
!> where \f$ \delta \tens{\variat} = \tens{\variat}^{n+1} - \tens{\variat}^n\f$ and
!> \f$ \tens{Rhs}^1 = \tens{Rhs}
!> - \divt \left( \tens{\variat}^n \otimes \rho \vect{u}
!>              - \mu \gradtt \tens{\variat}^n \right)\f$
!>
!>
!> It is in fact solved with the following iterative process:
!>
!> \f[
!> \tens{f_s}^{imp} \delta \tens{\variat}^k
!> + \divt \left( \delta \tens{\variat}^k \otimes \rho \vect{u}
!>              - \mu \gradtt \delta \tens{\variat}^k \right)
!> = \tens{Rhs}^k
!> \f]
!>
!> where \f$ \tens{Rhs}^k = \tens{Rhs}
!> - \tens{f_s}^{imp} \left(\tens{\variat}^k-\tens{\variat}^n \right)
!> - \divt \left( \tens{\variat}^k \otimes \rho \vect{u}
!>              - \mu \gradtt \tens{\variat}^k \right)\f$
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
!> \param[in]     ischcp        indicator
!>                               - 1 centered
!>                               - 0 2nd order
!> \param[in]     isstpp        indicator
!>                               - 1 without slope test
!>                               - 0 with slope test
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
!> \param[in]     coefats        boundary condition array for the variable
!>                               (Explicit part)
!> \param[in]     coefbts        boundary condition array for the variable
!>                               (Impplicit part)
!> \param[in]     cofafts        boundary condition array for the diffusion
!>                               of the variable (Explicit part)
!> \param[in]     cofbfts        boundary condition array for the diffusion
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
!> \param[in]     visccs        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
!> \param[in]     weighf        internal face weight between cells i j in case
!>                               of tensor diffusion
!> \param[in]     weighb        boundary face weight for cells i in case
!>                               of tensor diffusion
!> \param[in]     icvflb        global indicator of boundary convection flux
!>                               - 0 upwind scheme at all boundary faces
!>                               - 1 imposed flux at some boundary faces
!> \param[in]     icvfli        boundary face indicator array of convection flux
!>                               - 0 upwind scheme
!>                               - 1 imposed flux
!> \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$
!> \param[in]     smbrp         Right hand side \f$ \vect{Rhs}^k \f$
!> \param[in,out] pvar          current variable
!_______________________________________________________________________________

subroutine coditts &
 ( idtvar , ivar   , iconvp , idiffp , ndircp ,                   &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , idftnp , iswdyp ,                            &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp ,                   &
   relaxp , thetap ,                                              &
   pvara  , pvark  ,                                              &
   coefats , coefbts , cofafts , cofbfts ,                        &
   flumas , flumab ,                                              &
   viscfm , viscbm , viscfs , viscbs ,   visccs ,                 &
   weighf , weighb ,                                              &
   icvflb , icvfli ,                                              &
   fimp   ,                                                       &
   smbrp  ,                                                       &
   pvar   )

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
integer          ischcp , isstpp
integer          idftnp , iswdyp , icvflb
integer          iwarnp

integer          icvfli(nfabor)

double precision blencp , epsilp , epsrgp , climgp
double precision relaxp , thetap , epsrsp

double precision pvara(6,ncelet)
double precision pvark(6,ncelet)
double precision pvar(6,ncelet)
double precision coefats(6,ndimfb)
double precision cofafts(6,ndimfb)
double precision coefbts(6,6,ndimfb)
double precision cofbfts(6,6,ndimfb)
double precision flumas(nfac), flumab(nfabor)
double precision viscfm(*), viscbm(nfabor)
double precision viscfs(*), viscbs(nfabor)
double precision visccs(6,ncelet)
double precision weighf(2,nfac), weighb(nfabor)
double precision fimp(6,6,ncelet)
double precision smbrp(6,ncelet)

! Local variables

character(len=80) :: chaine
character(len=16) :: cnom
integer          f_id,isym
integer          inc,isweep,niterf,iel,nswmod
integer          iinvpe
integer          idtva0
integer          isou , jsou
integer          ibsize, iesize
integer          lvar, imasac

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
ibsize = 6
if (idftnp.eq.1) iesize = 1
if (idftnp.eq.6) iesize = 1

! Allocate temporary arrays
allocate(dam(6,6,ncelet))
allocate(dpvar(6,ncelet), smbini(6,ncelet))
if (iswdyp.ge.1) then
  allocate(adxk(6,ncelet), adxkm1(6,ncelet), dpvarm1(6,ncelet))
  allocate(rhs0(6,ncelet))
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
if (iesize.eq.6) allocate(xam(6*6*isym,nfac))

! iinvpe is useless in the vectorial framework
iinvpe = 0

!===============================================================================
! 1.  Building of the "simplified" matrix
!===============================================================================

if (iesize.eq.1) then

  call matrxts &
  !==========
   ( iconvp , idiffp , ndircp , isym   ,                            &
     thetap ,                                                       &
     coefbts , cofbfts , fimp   ,                                     &
     flumas , flumab , viscfm , viscbm ,                            &
     dam    , xam    )

elseif (iesize.eq.6) then

  call matrvts &
  !==========
   ( iconvp , idiffp , ndircp , isym   ,                            &
     thetap ,                                                       &
     coefbts , cofbfts , fimp   ,                                     &
     flumas , flumab , viscfm , viscbm ,                            &
     dam    , xam    )

endif

! For steady computations, the diagonal is relaxed
if (idtvar.lt.0) then
  !$omp parallel do private(isou, jsou)
  do iel = 1, ncel
    do isou = 1, 6
      do jsou = 1, 6
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

  ! The added convective scalar mass flux is:
  !      (thetex*Y_\face-imasac*Y_\celli)*mf.
  ! When building the explicit part of the rhs, one
  ! has to impose 0 on mass accumulation.
  imasac = 0
  call bilscts &
  !==========
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra ,                            &
   iwarnp , idftnp , imasac ,                                     &
   blencp , epsrgp , climgp , relaxp , thetex ,                   &
   pvara  , pvara  ,                                              &
   coefats , coefbts , cofafts , cofbfts ,                        &
   flumas , flumab , viscfs , viscbs ,   visccs ,                 &
   weighf , weighb ,                                              &
   icvflb , icvfli ,                                              &
   smbrp  )
endif

! Before looping, the RHS without reconstruction is stored in smbini

!$omp parallel do private(isou)
do iel = 1, ncel
  do isou = 1, 6
    smbini(isou,iel) = smbrp(isou,iel)
  enddo
enddo

! pvar is initialized on ncelet to avoid a synchronization

!$omp parallel do private(isou)
do iel = 1, ncelet
  do isou = 1, 6
    pvar(isou,iel) = pvark(isou,iel)
  enddo
enddo

! In the following, bilscts is called with inc=1,
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
  do isou = 1, 6
    smbrp(isou,iel) = 0.d0
  enddo
enddo

! The added convective scalar mass flux is:
!      (thetap*Y_\face-imasac*Y_\celli)*mf.
! When building the implicit part of the rhs, one
! has to impose 1 on mass accumulation.
imasac = 1

call bilscts &
!==========
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra ,                            &
   iwarnp , idftnp , imasac ,                                     &
   blencp , epsrgp , climgp , relaxp , thetap ,                   &
   pvar   , pvara  ,                                              &
   coefats , coefbts , cofafts , cofbfts ,                        &
   flumas , flumab , viscfs , viscbs ,  visccs ,                  &
   weighf , weighb ,                                              &
   icvflb , icvfli ,                                              &
   smbrp  )

! Dynamic relaxation
if (iswdyp.ge.1) then

  !$omp parallel do private(isou)
  do iel = 1, ncel
    do isou = 1, 6
      rhs0(isou,iel) = smbrp(isou,iel)
      smbini(isou,iel) = smbini(isou,iel)                        &
                 -fimp(isou,1,iel)*(pvar(1,iel) - pvara(1,iel))  &
                 -fimp(isou,2,iel)*(pvar(2,iel) - pvara(2,iel))  &
                 -fimp(isou,3,iel)*(pvar(3,iel) - pvara(3,iel))  &
                 -fimp(isou,4,iel)*(pvar(4,iel) - pvara(4,iel))  &
                 -fimp(isou,5,iel)*(pvar(5,iel) - pvara(5,iel))  &
                 -fimp(isou,6,iel)*(pvar(6,iel) - pvara(6,iel))
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
    do isou = 1, 6
      smbini(isou,iel) = smbini(isou,iel)                        &
                 -fimp(isou,1,iel)*(pvar(1,iel) - pvara(1,iel))  &
                 -fimp(isou,2,iel)*(pvar(2,iel) - pvara(2,iel))  &
                 -fimp(isou,3,iel)*(pvar(3,iel) - pvara(3,iel))  &
                 -fimp(isou,4,iel)*(pvar(4,iel) - pvara(4,iel))  &
                 -fimp(isou,5,iel)*(pvar(5,iel) - pvara(5,iel))  &
                 -fimp(isou,6,iel)*(pvar(6,iel) - pvara(6,iel))
      smbrp(isou,iel) = smbrp(isou,iel) + smbini(isou,iel)
    enddo
  enddo
endif

! --- Convergence test
residu = sqrt(cs_gdot(6*ncel, smbrp, smbrp))

! ---> RESIDU DE NORMALISATION
!    (NORME C.L +TERMES SOURCES+ TERMES DE NON ORTHOGONALITE)

allocate(w1(6, ncelet))  ! Allocate a temporary array

call promav(isym, ibsize, iesize, iinvpe, dam, xam, pvar, w1)

!$omp parallel do private(isou)
do iel = 1, ncel
   do isou = 1, 6
      w1(isou,iel) = w1(isou,iel) + smbrp(isou,iel)
   enddo
enddo

rnorm2 = cs_gdot(6*ncel, w1, w1)
rnorm = sqrt(rnorm2)

sinfo%rnsmbr = rnorm

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
      do isou = 1, 6
        dpvarm1(isou,iel) = dpvar(isou,iel)
        dpvar(isou,iel) = 0.d0
      enddo
    enddo
  else

    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou =1,6
        dpvar(isou,iel) = 0.d0
      enddo
    enddo
  endif

  ! iinvpe is useless in the vectorial framework
  iinvpe = 0

  ! Matrix block size
  ibsize = 6

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
      do isou = 1, 6
        adxkm1(isou,iel) = adxk(isou,iel)
        adxk(isou,iel) = - rhs0(isou,iel)
      enddo
    enddo

    call bilscts &
    !==========
   ( idtvar , lvar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
     ischcp , isstpp , inc    , imrgra ,                            &
     iwarnp , idftnp , imasac ,                                     &
     blencp , epsrgp , climgp , relaxp , thetap ,                   &
     dpvar  , dpvar  ,                                              &
     coefats , coefbts , cofafts , cofbfts ,                        &
     flumas , flumab , viscfs , viscbs ,  visccs ,                  &
     weighf , weighb ,                                              &
     icvflb , icvfli ,                                              &
     adxk   )

    ! ||E.dx^(k-1)-E.0||^2
    nadxkm1 = nadxk

    ! ||E.dx^k-E.0||^2
    nadxk = cs_gdot(6*ncel, adxk, adxk)

    ! < E.dx^k-E.0; r^k >
    paxkrk = cs_gdot(6*ncel, smbrp, adxk)

    ! Relaxation with respect to dx^k and dx^(k-1)
    if (iswdyp.ge.2) then

      ! < E.dx^(k-1)-E.0; r^k >
      paxm1rk = cs_gdot(6*ncel, smbrp, adxkm1)

      ! < E.dx^(k-1)-E.0; E.dx^k-E.0 >
      paxm1ax = cs_gdot(6*ncel, adxk, adxkm1)

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
      do isou = 1, 6
         pvar(isou,iel) = pvar(isou,iel) + dpvar(isou,iel)
      enddo
    enddo
  elseif (iswdyp.eq.1) then

    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou = 1, 6
         pvar(isou,iel) = pvar(isou,iel) + alph*dpvar(isou,iel)
      enddo
    enddo
  elseif (iswdyp.ge.2) then

    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou = 1, 6
         pvar(isou,iel) = pvar(isou,iel) + alph*dpvar(isou,iel)  &
                        + beta*dpvarm1(isou,iel)
      enddo
    enddo
  endif

  ! ---> Handle parallelism and periodicity

  if (irangp.ge.0.or.iperio.eq.1) then
    call syntis(pvar)
  endif

  ! --- Update the right hand and compute the new residual

  if (iswdyp.eq.0) then

    !$omp parallel do private(isou)
    do iel = 1, ncel
      ! smbini already contains unsteady terms and mass source terms
      ! of the RHS updated at each sweep
      do isou = 1, 6
        smbini(isou,iel) = smbini(isou,iel)                 &
                  - fimp(isou,1,iel)*dpvar(1,iel)           &
                  - fimp(isou,2,iel)*dpvar(2,iel)           &
                  - fimp(isou,3,iel)*dpvar(3,iel)           &
                  - fimp(isou,4,iel)*dpvar(4,iel)           &
                  - fimp(isou,5,iel)*dpvar(5,iel)           &
                  - fimp(isou,6,iel)*dpvar(6,iel)
        smbrp(isou,iel) = smbini(isou,iel)
      enddo
    enddo

  elseif (iswdyp.eq.1) then

    !$omp parallel do private(isou)
    do iel = 1, ncel
      ! smbini already contains unsteady terms and mass source terms
      ! of the RHS updated at each sweep
      do isou = 1, 6
        smbini(isou,iel) = smbini(isou,iel)                 &
                  - fimp(isou,1,iel)*alph*dpvar(1,iel)      &
                  - fimp(isou,2,iel)*alph*dpvar(2,iel)      &
                  - fimp(isou,3,iel)*alph*dpvar(3,iel)      &
                  - fimp(isou,4,iel)*alph*dpvar(4,iel)      &
                  - fimp(isou,5,iel)*alph*dpvar(5,iel)      &
                  - fimp(isou,6,iel)*alph*dpvar(6,iel)
        smbrp(isou,iel) = smbini(isou,iel)
      enddo
    enddo

  elseif (iswdyp.eq.2) then

    !$omp parallel do private(isou)
    do iel = 1, ncel
      ! smbini already contains unsteady terms and mass source terms
      ! of the RHS updated at each sweep
      do isou = 1, 6
        smbini(isou,iel) = smbini(isou,iel)                                  &
                  - fimp(isou,1,iel)*(alph*dpvar(1,iel)+beta*dpvarm1(1,iel)) &
                  - fimp(isou,2,iel)*(alph*dpvar(2,iel)+beta*dpvarm1(2,iel)) &
                  - fimp(isou,3,iel)*(alph*dpvar(3,iel)+beta*dpvarm1(3,iel)) &
                  - fimp(isou,4,iel)*(alph*dpvar(4,iel)+beta*dpvarm1(4,iel)) &
                  - fimp(isou,5,iel)*(alph*dpvar(5,iel)+beta*dpvarm1(5,iel)) &
                  - fimp(isou,6,iel)*(alph*dpvar(6,iel)+beta*dpvarm1(6,iel))
        smbrp(isou,iel) = smbini(isou,iel)
      enddo
    enddo
  endif

  ! The added convective scalar mass flux is:
  !      (thetex*Y_\face-imasac*Y_\celli)*mf.
  ! When building the implicit part of the rhs, one
  ! has to impose 1 on mass accumulation.
  imasac = 1

  call bilscts &
  !==========
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra ,                            &
   iwarnp , idftnp , imasac ,                                     &
   blencp , epsrgp , climgp , relaxp , thetap ,                   &
   pvar   , pvara  ,                                              &
   coefats , coefbts , cofafts , cofbfts ,                        &
   flumas , flumab , viscfs , viscbs ,   visccs ,                 &
   weighf , weighb ,                                              &
   icvflb , icvfli ,                                              &
   smbrp  )

  ! --- Convergence test
  residu = sqrt(cs_gdot(6*ncel, smbrp, smbrp))

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
if (abs(rnorm)/sqrt(6.d0).gt.epzero) then
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
! 4. Free solver setup
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
