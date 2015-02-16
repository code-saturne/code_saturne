!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

!> \file codits.f90
!>
!> \brief This function solves an advection diffusion equation with source terms
!> for one time step for the variable \f$ a \f$.
!>
!> The equation reads:
!>
!> \f[
!> f_s^{imp}(a^{n+1}-a^n)
!> + \divs \left( a^{n+1} \rho \vect{u} - \mu \grad a^{n+1} \right)
!> = Rhs
!> \f]
!>
!> This equation is rewritten as:
!>
!> \f[
!> f_s^{imp} \delta a
!> + \divs \left( \delta a \rho \vect{u} - \mu \grad \delta a \right)
!> = Rhs^1
!> \f]
!>
!> where \f$ \delta a = a^{n+1} - a^n\f$ and
!> \f$ Rhs^1 = Rhs - \divs( a^n \rho \vect{u} - \mu \grad a^n)\f$
!>
!>
!> It is in fact solved with the following iterative process:
!>
!> \f[
!> f_s^{imp} \delta a^k
!> + \divs \left(\delta a^k \rho \vect{u}-\mu\grad\delta a^k \right)
!> = Rhs^k
!> \f]
!>
!> where \f$Rhs^k=Rhs-f_s^{imp}(a^k-a^n)
!> - \divs \left( a^k\rho\vect{u}-\mu\grad a^k \right)\f$
!>
!> Be careful, it is forbidden to modify \f$ f_s^{imp} \f$ here!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     idtvar        indicateur du schema temporel
!> \param[in]     ivar          index of the current variable
!> \param[in]     iconvp        indicator
!>                               - 1 convection,
!>                               - 0 otherwise
!> \param[in]     idiffp        indicator
!>                               - 1 diffusion,
!>                               - 0 otherwise
!> \param[in]     ndircp        indicator (0 if the diagonal is stepped aside)
!> \param[in]     imrgra        indicator
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
!> \param[in]     iescap        compute the predictor indicator if 1
!> \param[in]     imucpp        indicator
!>                               - 0 do not multiply the convectiv term by Cp
!>                               - 1 do multiply the convectiv term by Cp
!> \param[in]     idftnp        indicator
!>                               - 0 the diffusivity is scalar
!>                               - 1 the diffusivity is a diagonal tensor
!>                               - 2 the diffusivity is a symmetric tensor
!> \param[in]     iswdyp        indicator
!>                               - 0 no dynamic relaxation
!>                               - 1 dynamic relaxation depending on
!>                                 \f$ \delta \varia^k \f$
!>                               - 2 dynamic relaxation depending on
!>                                 \f$ \delta \varia^k \f$  and
!>                                 \f$ \delta \varia^{k-1} \f$
!> \param[in]     iwarnp        verbosity
!> \param[in]     blencp        fraction of upwinding
!> \param[in]     epsilp        precision pour resol iter
!> \param[in]     epsrsp        relative precision for the iterative process
!> \param[in]     epsrgp        relative precision for the gradient
!>                               reconstruction
!> \param[in]     climgp        clipping coefficient for the computation of
!>                               the gradient
!> \param[in]     extrap        coefficient for extrapolation of the gradient
!> \param[in]     relaxp        coefficient of relaxation
!> \param[in]     thetap        weighting coefficient for the theta-schema,
!>                               - thetap = 0: explicit scheme
!>                               - thetap = 0.5: time-centered
!>                               scheme (mix between Crank-Nicolson and
!>                               Adams-Bashforth)
!>                               - thetap = 1: implicit scheme
!> \param[in]     pvara         variable at the previous time step
!>                               \f$ a^n \f$
!> \param[in]     pvark         variable at the previous sub-iteration
!>                               \f$ a^k \f$.
!>                               If you sub-iter on Navier-Stokes, then
!>                               it allows to initialize by something else than
!>                               pvara (usually pvar=pvara)
!> \param[in]     coefap        boundary condition array for the variable
!>                               (explicit part)
!> \param[in]     coefbp        boundary condition array for the variable
!>                               (implicit part)
!> \param[in]     cofafp        boundary condition array for the diffusion
!>                               of the variable (explicit part)
!> \param[in]     cofbfp        boundary condition array for the diffusion
!>                               of the variable (implicit part)
!> \param[in]     flumas        mass flux at interior faces
!> \param[in]     flumab        mass flux at boundary faces
!> \param[in]     viscfm        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the matrix
!> \param[in]     viscbm        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at boundary faces for the matrix
!> \param[in]     visccm        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
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
!> \param[in]     rovsdt        \f$ f_s^{imp} \f$
!> \param[in]     smbrp         Right hand side \f$ Rhs^k \f$
!> \param[in,out] pvar          current variable
!> \param[in,out] dpvar         last variable increment
!> \param[in]     xcpp          array of specific heat (Cp)
!> \param[out]    eswork        prediction-stage error estimator
!>                              (if iescap > 0)
!_______________________________________________________________________________

subroutine codits &
 ( idtvar , ivar   , iconvp , idiffp , ndircp ,                   &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   pvara  , pvark  ,                                              &
   coefap , coefbp , cofafp , cofbfp , flumas , flumab ,          &
   viscfm , viscbm , visccm , viscfs , viscbs , visccs ,          &
   weighf , weighb ,                                              &
   icvflb , icvfli ,                                              &
   rovsdt , smbrp  , pvar   , dpvar  ,                            &
   xcpp   , eswork )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

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
integer          iwarnp
integer          imucpp , idftnp , iswdyp , icvflb

integer          icvfli(nfabor)

double precision blencp , epsilp , epsrgp , climgp , extrap
double precision relaxp , thetap , epsrsp

double precision pvara(ncelet), pvark(ncelet)
double precision coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision flumas(nfac), flumab(nfabor)
double precision viscfm(nfac), viscbm(nfabor)
double precision visccm(ncelet)
double precision viscfs(nfac), viscbs(nfabor)
double precision visccs(ncelet)
double precision weighf(2,nfac), weighb(nfabor)
double precision rovsdt(ncelet), smbrp(ncelet)
double precision pvar(ncelet)
double precision dpvar(ncelet)
double precision eswork(ncelet)
double precision xcpp(ncelet)

! Local variables

character(len=80) :: chaine
character(len=16) :: cnom
integer          f_id,isym
integer          inc,isweep,niterf,iccocg,iel,nswmod
integer          itenso,iinvpe, iinvpp
integer          idtva0
integer          lvar
integer          ibsize, iesize

double precision residu, rnorm, ressol, rnorm2
double precision thetex, nadxkm1, nadxk, paxm1ax, paxm1rk, paxkrk, alph, beta

type(solving_info) sinfo

double precision, allocatable, dimension(:) :: dam
double precision, allocatable, dimension(:,:) :: xam
double precision, allocatable, dimension(:) :: smbini, w1, adxk, adxkm1, dpvarm1
double precision, allocatable, dimension(:) :: rhs0

!===============================================================================

!===============================================================================
! 0.  Initialization
!===============================================================================

! Allocate temporary arrays
allocate(dam(ncelet))
allocate(smbini(ncelet))
if (iswdyp.ge.1) then
  allocate(adxk(ncelet), adxkm1(ncelet), dpvarm1(ncelet))
  allocate(rhs0(ncelet))
endif

! Names
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

allocate(xam(isym,nfac))

! Matrix block size
ibsize = 1
iesize = 1

! PRISE EN COMPTE DE LA PERIODICITE

! Initialisation pour test avant promav
itenso = 0
iinvpe = 0

if (iperio.eq.1) then

!    Par defaut, toutes les periodicites seront traitees,
!      les variables etant assimilees a des scalaires (meme si ce sont
!      des composantes de vecteurs ou de tenseur)

  iinvpe = 1

  if (ivar.eq.ir11.or.ivar.eq.ir12.or.             &
      ivar.eq.ir13.or.ivar.eq.ir22.or.             &
      ivar.eq.ir23.or.ivar.eq.ir33) then

    !    Pour les tensions de Reynolds, et les tpucou
    !      seules seront echangees les informations sur les faces periodiques
    !      de translation ; on ne touche pas aux informations
    !      relatives aux faces de periodicite de rotation.
    itenso = 1

    !      Lors de la resolution par increments, on echangera egalement les
    !      informations relatives aux faces de periodicite de translation.
    !      Pour les faces de periodicite de rotation, l'increment sera
    !      annule en appelant syncmp au lieu de synsca (iinvpe=2).
    iinvpe = 2

  endif

endif

!===============================================================================
! 1.  Building of the "simplified" matrix
!===============================================================================

call matrix &
!==========
 ( iconvp , idiffp , ndircp , isym   ,                            &
   thetap , imucpp ,                                              &
   coefbp , cofbfp , rovsdt , flumas , flumab , viscfm , viscbm , &
   xcpp   , dam    , xam    )

! For steady computations, the diagonal is relaxed
if (idtvar.lt.0) then
  !$omp parallel do
  do iel = 1, ncel
    dam(iel) = dam(iel)/relaxp
  enddo
endif

!===============================================================================
! 2. Iterative process to handle non orthogonlaities (starting from the second
! iteration).
!===============================================================================

! Application du theta schema

! On calcule le bilan explicite total
thetex = 1.d0 - thetap

! Si THETEX=0, ce n'est pas la peine d'en rajouter
if (abs(thetex).gt.epzero) then
  inc    = 1
  iccocg = 1

  call bilsca &
  !==========
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   iwarnp , imucpp , idftnp ,                                     &
   blencp , epsrgp , climgp , extrap , relaxp , thetex ,          &
   pvara  , pvara  , coefap , coefbp , cofafp , cofbfp ,          &
   flumas , flumab , viscfs , viscbs , visccs , xcpp   ,          &
   weighf , weighb ,                                              &
   icvflb , icvfli ,                                              &
   smbrp  )
endif

! Before looping, the RHS without reconstruction is stored in smbini

!$omp parallel do
do iel = 1, ncel
  smbini(iel) = smbrp(iel)
enddo

! pvar is initialized on ncelet to avoid a synchronization

!$omp parallel do
do iel = 1, ncelet
  pvar(iel) = pvark(iel)
enddo

! In the following, bilsca is called with inc=1,
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

!$omp parallel do
do iel = 1, ncel
  smbrp(iel)  = 0.d0
enddo

iccocg = 1

call bilsca &
!==========
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   iwarnp , imucpp , idftnp ,                                     &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   pvar   , pvara  , coefap , coefbp , cofafp , cofbfp ,          &
   flumas , flumab , viscfs , viscbs , visccs , xcpp   ,          &
   weighf , weighb ,                                              &
   icvflb , icvfli ,                                              &
   smbrp  )

if (iswdyp.ge.1) then
  !$omp parallel do
  do iel = 1, ncel
    rhs0(iel) = smbrp(iel)
    smbini(iel) = smbini(iel) - rovsdt(iel)*(pvar(iel) - pvara(iel))
    smbrp(iel)  = smbrp(iel) + smbini(iel)

    adxkm1(iel) = 0.d0
    adxk(iel) = 0.d0
    dpvar(iel) = 0.d0
  enddo

  ! ||A.dx^0||^2 = 0
  nadxk = 0.d0

else

  !$omp parallel do
  do iel = 1, ncel
    smbini(iel) = smbini(iel) - rovsdt(iel)*(pvar(iel) - pvara(iel))
    smbrp(iel)  = smbrp(iel) + smbini(iel)
  enddo
endif

! --- Right hand side residual
residu = sqrt(cs_gdot(ncel,smbrp,smbrp))

! ---> RESIDU DE NORMALISATION
!    (NORME C.L +TERMES SOURCES+ TERMES DE NON ORTHOGONALITE)

!       Attention, lors de l'appel a promav, ici pour une variable qui
!         n'est pas en increments et qui est supposee initialisee
!         y compris dans le halo.
!         Pour les variables vitesse et les tensions de Reynolds
!         (IINVPE=2), il ne faudra donc pas annuler le halo
!         des periodicites de rotation, mais au contraire le laisser
!         inchange.
!         Pour les autres variables (scalaires) IINVPE=1 permettra de
!         tout echanger, meme si c'est superflu.

! Allocate a temporary array
allocate(w1(ncelet))

if (iinvpe.eq.2) then
  iinvpp = 3
else
  iinvpp = iinvpe
endif

call promav(isym,ibsize,iesize,iinvpp,dam,xam,pvar,w1)

!$omp parallel do
do iel = 1, ncel
  w1(iel) = w1(iel) + smbrp(iel)
enddo

rnorm2 = cs_gdot(ncel,w1,w1)
rnorm = sqrt(rnorm2)
sinfo%rnsmbr = rnorm

! Free memory
deallocate(w1)

! Warning: for Weight Matrix, one and only one sweep is done.
nswmod = max(nswrsp, 1)

! Reconstruction loop (beginning)
!--------------------------------
sinfo%nbivar = 0
isweep = 1

do while (isweep.le.nswmod.and.residu.gt.epsrsp*rnorm.or.isweep.eq.1)

  ! --- Solving on the increment dpvar

  if (iswdyp.ge.1) then
    !$omp parallel do
    do iel = 1, ncel
      dpvarm1(iel) = dpvar(iel)
      dpvar(iel) = 0.d0
    enddo
  else
    !$omp parallel do
    do iel = 1, ncel
      dpvar(iel) = 0.d0
    enddo
  endif

  ! Solver residual
  ressol = residu

  call sles_solve_native(f_id, chaine,                                &
                         isym, ibsize, iesize, dam, xam, iinvpe,      &
                         epsilp, rnorm, niterf, ressol, smbrp, dpvar)

  ! Dynamic relaxation of the system
  if (iswdyp.ge.1) then

    ! Computation of the variable relaxation coefficient
    lvar = 0

    !$omp parallel do
    do iel = 1, ncelet
      adxkm1(iel) = adxk(iel)
      adxk(iel) = - rhs0(iel)
    enddo

    call bilsca &
    !==========
   ( idtvar , lvar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
     ischcp , isstpp , inc    , imrgra , iccocg ,                   &
     iwarnp , imucpp , idftnp ,                                     &
     blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
     dpvar  , dpvar  , coefap , coefbp , cofafp , cofbfp ,          &
     flumas , flumab , viscfs , viscbs , visccs , xcpp   ,          &
     weighf , weighb ,                                              &
     icvflb , icvfli ,                                              &
     adxk   )

    ! ||E.dx^(k-1)-E.0||^2
    nadxkm1 = nadxk

    ! ||E.dx^k-E.0||^2
    nadxk = cs_gdot(ncel, adxk, adxk)

    ! < E.dx^k-E.0; r^k >
    paxkrk = cs_gdot(ncel, smbrp, adxk)

    ! Relaxation with respect to dx^k and dx^(k-1)
    if (iswdyp.ge.2) then

      ! < E.dx^(k-1)-E.0; r^k >
      paxm1rk = cs_gdot(ncel, smbrp, adxkm1)

      ! < E.dx^(k-1)-E.0; E.dx^k-E.0 >
      paxm1ax = cs_gdot(ncel, adxk, adxkm1)

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
    if (iwarnp.ge.2) then
      write(nfecra,1200) cnom, isweep, alph, beta, &
                         paxkrk, nadxk, paxm1rk, nadxkm1, paxm1ax
    endif

  endif

  ! --- Update the solution with the increment

  if (iswdyp.eq.0) then
    !$omp parallel do
    do iel = 1, ncel
      pvar(iel) = pvar(iel) + dpvar(iel)
    enddo
  elseif (iswdyp.eq.1) then
    if (alph.lt.0.d0) goto 100
    !$omp parallel do
    do iel = 1, ncel
      pvar(iel) = pvar(iel) + alph*dpvar(iel)
    enddo
  elseif (iswdyp.ge.2) then
    !$omp parallel do
    do iel = 1, ncel
      pvar(iel) = pvar(iel) + alph*dpvar(iel) + beta*dpvarm1(iel)
    enddo
  endif

  ! ---> Handle parallelism and periodicity
  !      (periodicity of rotation is not ensured here)
  if (irangp.ge.0 .or. iperio.eq.1) then
    if (itenso.eq.0) then
      call synsca (pvar)
    else if (itenso.eq.1) then
      call syncmp (pvar)
    endif
  endif

  ! --- Update the right hand side And compute the new residual

  iccocg = 0

  if (iswdyp.eq.0) then
    !$omp parallel do
    do iel = 1, ncel
      ! smbini already contains unsteady terms and mass source terms
      ! of the RHS updated at each sweep
      smbini(iel) = smbini(iel) - rovsdt(iel)*dpvar(iel)
      smbrp(iel)  = smbini(iel)
    enddo
  elseif (iswdyp.eq.1) then
    !$omp parallel do
    do iel = 1, ncel
      ! smbini already contains unsteady terms and mass source terms
      ! of the RHS updated at each sweep
      smbini(iel) = smbini(iel) - rovsdt(iel)*alph*dpvar(iel)
      smbrp(iel)  = smbini(iel)
    enddo
  elseif (iswdyp.ge.2) then
    !$omp parallel do
    do iel = 1, ncel
      ! smbini already contains unsteady terms and mass source terms
      ! of the RHS updated at each sweep
      smbini(iel) = smbini(iel)                                     &
                  - rovsdt(iel)*(alph*dpvar(iel)+beta*dpvarm1(iel))
      smbrp(iel)  = smbini(iel)
    enddo
  endif

  call bilsca &
  !==========
 ( idtvar , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   iwarnp , imucpp , idftnp ,                                     &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   pvar   , pvara  , coefap , coefbp , cofafp , cofbfp ,          &
   flumas , flumab , viscfs , viscbs , visccs , xcpp   ,          &
   weighf , weighb ,                                              &
   icvflb , icvfli ,                                              &
   smbrp  )

  ! --- Convergence test
  residu = sqrt(cs_gdot(ncel, smbrp, smbrp))

  ! Writing
  sinfo%nbivar = sinfo%nbivar + niterf

  ! Writing
  if (iwarnp.ge.2) then
    write(nfecra,1000) cnom, isweep, residu, rnorm
    write(nfecra,1010) cnom, isweep, niterf
  endif

  isweep = isweep + 1

enddo
! --- Reconstruction loop (end)

100 continue

! Writing: convergence
if (abs(rnorm).gt.epzero) then
  sinfo%resvar = residu/rnorm
else
  sinfo%resvar = 0.d0
endif

! Save convergence info
if (ivar.gt.0) then
  call field_set_key_struct_solving_info(ivarfl(ivar), sinfo)
endif

if (iwarnp.ge.1) then
  if (residu.le.epsrsp*rnorm) then
    write(nfecra,1000) cnom,isweep-1,residu,rnorm

  ! Writing: non-convergence
  else if (isweep.gt.nswmod) then
    write(nfecra,1100) cnom, nswmod
  endif
endif

!===============================================================================
! 3. After having computed the new value, an estimator is computed for the
! prediction step of the velocity.
!===============================================================================

if (iescap.gt.0) then

  ! ---> Computation of the estimator of the current component

  ! smbini already contains unsteady terms and mass source terms
  ! of the RHS updated at each sweep

  !$omp parallel do
  do iel = 1,ncel
    smbrp(iel) = smbini(iel) - rovsdt(iel)*dpvar(iel)
  enddo

  inc    = 1
  iccocg = 1
  ! Without relaxation even for a statonnary computation
  idtva0 = 0

  call bilsca &
  !==========
 ( idtva0 , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   iwarnp , imucpp , idftnp ,                                     &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   pvar   , pvara  , coefap , coefbp , cofafp , cofbfp ,          &
   flumas , flumab , viscfs , viscbs , visccs , xcpp   ,          &
   weighf , weighb ,                                              &
   icvflb , icvfli ,                                              &
   smbrp  )

  ! Contribution of the current component to the L2 norm stored in eswork

  !$omp parallel do
  do iel = 1,ncel
    eswork(iel) = (smbrp(iel) / volume(iel))**2
  enddo

endif

!===============================================================================
! 4. Free solver setup
!===============================================================================

call sles_free_native(f_id, chaine)

! Free memory
deallocate(dam, xam)
deallocate(smbini)
if (iswdyp.ge.1) deallocate(adxk, adxkm1, dpvarm1,rhs0)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format ( &
 1X,A16,' : CV-DIF-TS',I5,' IT - RES= ',E12.5,' NORME= ', E12.5)
 1010 format ( &
 1X,A16,' : Current reconstruction sweep = ',I5,' - Sweeps for solver = ',I5)
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
 1010 format ( &
 1X,A16,' : Current reconstruction sweep = ',I5,' - Sweeps for solver = ',I5)
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
