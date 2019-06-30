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

!===============================================================================
! Function:
! ---------

!> \file resopv.f90
!>
!> \brief This subroutine performs the pressure correction step of the Navier
!> Stokes equations for incompressible or slightly compressible flows for
!> the coupled velocity components solver.
!>
!> This function solves the following Poisson equation on the pressure:
!> \f[
!>     D \left( \Delta t, \delta p \right) =
!> \divs \left( \rho \vect{\widetilde{u}}\right)
!>     - \Gamma^n
!>     + \dfrac{\rho^n - \rho^{n-1}}{\Delta t}
!> \f]
!> The mass flux is then updated as follows:
!> \f[
!>  \dot{m}^{n+1}_\ij = \dot{m}^{n}_\ij
!>                    - \Delta t \grad_\fij \delta p \cdot \vect{S}_\ij
!> \f]
!>
!> Remarks:
!> - an iterative process is used to solve the Poisson equation.
!> - if the coefficient arak is set to 1, the the Rhie & Chow filter is
!>   activated.
!>
!> Please refer to the
!> <a href="../../theory.pdf#resopv"><b>resopv</b></a>
!> section of the theory guide for more informations.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     iterns        Navier-Stokes iteration number
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     nfbpcd        number of faces with condensation source term
!> \param[in]     ncmast        number of cells with condensation source terms
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     ifbpcd        index of faces with condensation source term
!> \param[in]     ltmast        index of cells with condensation source terms
!> \param[in]     isostd        indicator of standard outlet and index
!>                               of the reference outlet face
!> \param[in]     dt            time step (per cell)
!> \param[in]     vel           velocity
!> \param[in]     coefav        boundary condition array for the variable
!>                               (explicit part)
!> \param[in]     coefbv        boundary condition array for the variable
!>                               (implicit part)
!> \param[in]     coefa_dp      boundary conditions for the pressure increment
!> \param[in]     coefb_dp      boundary conditions for the pressure increment
!> \param[in]     smacel        variable value associated to the mass source
!>                               term (for ivar=ipr, smacel is the mass flux
!>                               \f$ \Gamma^n \f$)
!> \param[in]     spcond        variable value associated to the condensation
!>                               source term (for ivar=ipr, spcond is the flow rate
!>                               \f$ \Gamma_{s,cond}^n \f$)
!> \param[in]     svcond        variable value associated to the condensation
!>                              source term (for ivar=ipr, svcond is the flow rate
!>                              \f$ \Gamma_{v, cond}^n \f$)
!> \param[in]     frcxt         external forces making hydrostatic pressure
!> \param[in]     dfrcxt        variation of the external forces
!>                              composing the hydrostatic pressure
!> \param[in]     tpucou        non scalar time step in case of
!>                               velocity pressure coupling
!> \param[in]     viscf         visc*surface/dist aux faces internes
!> \param[in]     viscb         visc*surface/dist aux faces de bord
!> \param[in]     phi           potential to be solved (pressure increment)
!> \param[in]     tslagr        coupling term for the Lagrangian module
!_______________________________________________________________________________

subroutine resopv &
 ( nvar   , iterns , ncesmp , nfbpcd , ncmast ,                   &
   icetsm , ifbpcd , ltmast , isostd ,                            &
   dt     , vel    ,                                              &
   coefav , coefbv , coefa_dp        , coefb_dp ,                 &
   smacel , spcond , svcond ,                                     &
   frcxt  , dfrcxt , tpucou ,                                     &
   viscf  , viscb  ,                                              &
   phi    , tslagr )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use atincl
use paramx
use dimens, only: ndimfb
use numvar
use entsor
use cstphy
use cstnum
use optcal
use pointe, only: itypfb, b_head_loss, gamcav, dgdpca
use albase
use parall
use period
use ppincl, only: icondv
use lagran
use cplsat
use mesh
use field
use field_operator
use cavitation
use vof
use cs_f_interfaces
use cs_c_bindings
use cs_tagms, only:s_metal

!===============================================================================

implicit none

! Arguments

integer          nvar  , iterns
integer          ncesmp, nfbpcd, ncmast

integer          icetsm(ncesmp), ifbpcd(nfbpcd)
integer          ltmast(ncelet)
integer          isostd(nfabor+1)

double precision, dimension (1:ncelet), target :: dt
double precision smacel(ncesmp,nvar), spcond(nfbpcd,nvar)
double precision svcond(ncelet,nvar)
double precision frcxt(3,ncelet), dfrcxt(3,ncelet)
double precision, dimension (1:6,1:ncelet), target :: tpucou
double precision viscf(nfac), viscb(nfabor)
double precision phi(ncelet)
double precision tslagr(ncelet,*)
double precision coefav(3  ,nfabor)
double precision coefbv(3,3,nfabor)
double precision vel   (3  ,ncelet)
double precision coefa_dp(nfabor)
double precision coefb_dp(nfabor)

! Local variables

character(len=80) :: chaine
integer          lchain, iautof
integer          iccocg, inc   , iprev, init  , isym
integer          ii, jj, iel   , ifac  , ifac0 , iel0
integer          nswmpr
integer          isweep, niterf
integer          iflmb0, ifcsor
integer          nswrgp, imligp, iwarnp
integer          iflmas, iflmab
integer          idiffp, iconvp, ndircp
integer          indhyd
integer          itypfl
integer          isou  , ibsize, iesize
integer          imucpp, idftnp, iswdyp
integer          iescap, ircflp, ischcp, isstpp, ivar, f_id0
integer          nswrsp
integer          imvisp
integer          iflid, iflwgr, f_dim, imasac
integer          f_id
integer          icvflb
integer          ivoid(1)

double precision residu, phydr0, rnormp
double precision ardtsr, arsr  , thetap
double precision dtsrom
double precision epsrgp, climgp, extrap, epsilp
double precision drom  , dronm1, relaxp
double precision hint, qimp, qimpv(3), epsrsp, blencp
double precision ressol, rnorm2
double precision nadxkm1, nadxk, paxm1ax, paxm1rk, paxkrk, alph, beta
double precision visci(3,3), fikis, viscis, distfi
double precision cfl, kpdc, rho, pimp, bpmasf
double precision normp

type(solving_info) sinfo
type(var_cal_opt) :: vcopt_p, vcopt_u

double precision rvoid(1)

double precision, allocatable, dimension(:) :: dam, xam
double precision, allocatable, dimension(:) :: res, phia
double precision, dimension(:,:), allocatable :: gradp
double precision, allocatable, dimension(:) :: coefaf_dp, coefbf_dp
double precision, allocatable, dimension(:) :: coefap, coefbp, coefa_dp2
double precision, allocatable, dimension(:) :: coefa_rho, coefb_rho
double precision, allocatable, dimension(:) :: cofafp, cofbfp, coefaf_dp2
double precision, allocatable, dimension(:) :: rhs, rovsdt
double precision, allocatable, dimension(:) :: hydro_pres
double precision, allocatable, dimension(:) :: velflx, velflb, ddphi
double precision, allocatable, dimension(:,:) :: coefar, cofafr
double precision, allocatable, dimension(:,:,:) :: coefbr, cofbfr
double precision, allocatable, dimension(:) :: adxk, adxkm1, dphim1, rhs0
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb
double precision, allocatable, dimension(:,:) :: frchy, dfrchy
double precision, dimension(:), pointer :: coefa_p, coefb_p
double precision, dimension(:), pointer :: coefaf_p, coefbf_p
double precision, allocatable, dimension(:) :: iflux, bflux
double precision, allocatable, dimension(:) :: xunsro
double precision, allocatable, dimension(:), target :: xdtsro
double precision, allocatable, dimension(:), target  :: divu
double precision, allocatable, dimension(:,:), target :: tpusro
double precision, dimension(:), pointer :: viscap
double precision, dimension(:,:), pointer :: vitenp
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: brom, crom, croma, broma
double precision, dimension(:), pointer :: brom_eos, crom_eos
double precision, dimension(:), pointer :: cvar_pr
double precision, dimension(:), pointer :: cpro_divu
double precision, dimension(:), pointer :: cpro_wgrec_s
double precision, dimension(:,:), pointer :: cpro_wgrec_v
double precision, dimension(:), pointer :: c_estim_der
double precision, dimension(:), pointer :: cpro_tsrho
double precision, allocatable, dimension(:) :: surfbm, dphi
double precision, allocatable, dimension(:) :: ipro_visc, bpro_visc
double precision, allocatable, dimension(:) :: cpro_visc
double precision, allocatable, dimension(:,:) :: cpro_vitenp
double precision, allocatable, dimension(:,:) :: trav
double precision, dimension(:,:), pointer :: cpro_poro_div_duq
double precision, dimension(:), pointer :: cpro_rho_mass, bpro_rho_mass
double precision, dimension(:), allocatable, target :: cpro_rho_tc, bpro_rho_tc

!===============================================================================

!===============================================================================
! 1. Initialisations
!===============================================================================

! Initializations to avoid compiler warnings
rnorm2 = 0.d0
niterf = 0

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_p)

! Allocate temporary arrays
allocate(dam(ncelet), xam(nfac))
allocate(res(ncelet), phia(ncelet))
allocate(rhs(ncelet), rovsdt(ncelet))
allocate(iflux(nfac), bflux(ndimfb))
allocate(dphi(ncelet))
allocate(trav(3, ncelet))
iswdyp = vcopt_p%iswdyn
if (iswdyp.ge.1) allocate(adxk(ncelet), adxkm1(ncelet),   &
                          dphim1(ncelet), rhs0(ncelet))
if (icalhy.eq.1) allocate(frchy(ndim,ncelet),             &
                          dfrchy(ndim,ncelet), hydro_pres(ncelet))

! Diffusive flux Boundary conditions for delta P
allocate(coefaf_dp(ndimfb), coefbf_dp(ndimfb))

! Associate pointers to pressure diffusion coefficient
viscap => dt(:)
if (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then
  vitenp => tpucou(:,:)
endif

! Index of the field
iflid = ivarfl(ipr)
call field_get_key_struct_solving_info(iflid, sinfo)

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

call field_get_id_try("predicted_vel_divergence", f_id)
if (f_id.ge.0) then
  call field_get_val_s(f_id, cpro_divu)
else
  allocate(divu(ncelet))
  cpro_divu => divu
endif

! --- Writing
call field_get_name(ivarfl(ipr), chaine)
lchain = 16

f_id0 = -1

! --- Boundary conditions

call field_get_coefa_s(ivarfl(ipr), coefa_p)
call field_get_coefb_s(ivarfl(ipr), coefb_p)
call field_get_coefaf_s(ivarfl(ipr), coefaf_p)
call field_get_coefbf_s(ivarfl(ipr), coefbf_p)

! --- Physical quantities
call field_get_val_s(icrom, crom_eos)
if (icalhy.eq.1.or.idilat.gt.1.or.irovar.eq.1) then
  call field_get_val_prev_s(icrom, croma)
endif
call field_get_val_s(ibrom, brom_eos)
if (irovar.eq.1) then
  call field_get_val_prev_s(ibrom, broma)
endif

if (irovar.eq.1.and.idilat.gt.1) then
  call field_get_id("density_mass", f_id)
  call field_get_val_s(f_id, cpro_rho_mass)
  call field_get_id("boundary_density_mass", f_id)
  call field_get_val_s(f_id, bpro_rho_mass)

  ! Time interpolated density
  if (vcopt_u%thetav .lt. 1.d0 .and. iterns.gt.1) then
    allocate(cpro_rho_tc(ncelet))
    allocate(bpro_rho_tc(nfabor))

    do iel = 1, ncelet
      cpro_rho_tc(iel) = vcopt_u%thetav * cpro_rho_mass(iel) &
        + (1.d0 - vcopt_u%thetav) * croma(iel)
    enddo

    crom => cpro_rho_tc

    do ifac = 1, nfabor
      bpro_rho_tc(ifac) = vcopt_u%thetav * bpro_rho_mass(ifac) &
        + (1.d0 - vcopt_u%thetav) * broma(ifac)
    enddo

    brom => bpro_rho_tc

  else
    crom => cpro_rho_mass
    brom => bpro_rho_mass
  endif

! Weakly variable density algo. (idilat <=1) or constant density
else
  crom => crom_eos
  brom => brom_eos
endif

call field_get_val_s(ivarfl(ipr), cvar_pr)

call field_get_key_int(ivarfl(ipr), kimasf, iflmas)
call field_get_key_int(ivarfl(ipr), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

! --- Solving options
isym  = 1
if(vcopt_p%iconv.gt.0 ) then
  isym  = 2
endif

! Matrix block size
ibsize = 1
iesize = 1

! Initialization dedicated to VOF algo.
if (ivofmt.ge.0) then
  ! The pressure correction is done through the volumetric flux (that is
  ! the convective flux of the void fraction), not the mass flux
  call field_get_key_int(ivarfl(ivolf2), kimasf, iflmas)
  call field_get_key_int(ivarfl(ivolf2), kbmasf, iflmab)
  call field_get_val_s(iflmas, imasfl)
  call field_get_val_s(iflmab, bmasfl)
endif

! Calculation of dt/rho
if (ivofmt.ge.0.or.idilat.eq.4) then

  ! Allocate and initialize specific arrays
  allocate(xunsro(ncelet))
  do iel = 1, ncel
    xunsro(iel) = 1.d0/crom(iel)
  enddo
  call synsca(xunsro)

  if (iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then
    allocate(xdtsro(ncelet))
    do iel = 1, ncel
      xdtsro(iel) = dt(iel)/crom(iel)
    enddo
    call synsca(xdtsro)
  elseif (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then
    allocate(tpusro(6,ncelet))
    do iel = 1, ncel
      tpusro(1,iel) = vitenp(1,iel)*xunsro(iel)
      tpusro(2,iel) = vitenp(2,iel)*xunsro(iel)
      tpusro(3,iel) = vitenp(3,iel)*xunsro(iel)
      tpusro(4,iel) = vitenp(4,iel)*xunsro(iel)
      tpusro(5,iel) = vitenp(5,iel)*xunsro(iel)
      tpusro(6,iel) = vitenp(6,iel)*xunsro(iel)
    enddo
    call syntis(tpusro)
  endif

  ! Associate pointers to pressure diffusion coefficient
  viscap => xdtsro(:)
  if (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then
    vitenp => tpusro(:,:)
  endif

endif

!===============================================================================
! 2. Compute an approximated pressure increment if needed
!    that is when there are buoyancy terms (gravity and variable density)
!    with a free outlet.
!===============================================================================

! Standard initialization
do ifac = 1, nfac
  iflux(ifac) = 0.d0
enddo

do ifac = 1, nfabor
  coefa_dp(ifac) = 0.d0
  coefaf_dp(ifac) = 0.d0
  coefb_dp(ifac) = coefb_p(ifac)
  coefbf_dp(ifac) = coefbf_p(ifac)
  bflux(ifac) = 0.d0
enddo

! Compute a pseudo hydrostatic pressure increment stored
! in hydro_pres(.) with Homogeneous Neumann BCs everywhere
if (iphydr.eq.1.and.icalhy.eq.1) then

  ifcsor = isostd(nfabor+1)
  if (irangp.ge.0) then
    call parcmx (ifcsor)
  endif

  ! This computation is needed only if there are outlet faces
  if (ifcsor.le.0.and.iatmst.eq.0) then
    indhyd = 0
  else

    ! Work arrays for BCs
    allocate(coefap(ndimfb), cofafp(ndimfb))
    allocate(coefbp(ndimfb), cofbfp(ndimfb))

    do ifac = 1, nfabor

      iel = ifabor(ifac)

      if (iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then
        hint = dt(iel)/distb(ifac)
      ! Symmetric tensor diffusivity
      elseif (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then

        visci(1,1) = vitenp(1,iel)
        visci(2,2) = vitenp(2,iel)
        visci(3,3) = vitenp(3,iel)
        visci(1,2) = vitenp(4,iel)
        visci(2,1) = vitenp(4,iel)
        visci(2,3) = vitenp(5,iel)
        visci(3,2) = vitenp(5,iel)
        visci(1,3) = vitenp(6,iel)
        visci(3,1) = vitenp(6,iel)

        ! ||Ki.S||^2
        viscis = ( visci(1,1)*surfbo(1,ifac)       &
                 + visci(1,2)*surfbo(2,ifac)       &
                 + visci(1,3)*surfbo(3,ifac))**2   &
               + ( visci(2,1)*surfbo(1,ifac)       &
                 + visci(2,2)*surfbo(2,ifac)       &
                 + visci(2,3)*surfbo(3,ifac))**2   &
               + ( visci(3,1)*surfbo(1,ifac)       &
                 + visci(3,2)*surfbo(2,ifac)       &
                 + visci(3,3)*surfbo(3,ifac))**2

        ! IF.Ki.S
        fikis = ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,1)   &
                + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,1)   &
                + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,1)   &
                )*surfbo(1,ifac)                              &
              + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,2)   &
                + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,2)   &
                + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,2)   &
                )*surfbo(2,ifac)                              &
              + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,3)   &
                + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,3)   &
                + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,3)   &
                )*surfbo(3,ifac)

        distfi = distb(ifac)

        ! Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
        ! NB: eps =1.d-1 must be consistent with vitens.f90
        fikis = max(fikis, 1.d-1*sqrt(viscis)*distfi)

        hint = viscis/surfbn(ifac)/fikis

      endif

      ! LOCAL Neumann Boundary Conditions on the hydrostatic pressure
      !--------------------------------------------------------------

      qimp = 0.d0

      call set_neumann_scalar &
           !==================
         ( coefap(ifac), cofafp(ifac),             &
           coefbp(ifac), cofbfp(ifac),             &
           qimp        , hint )

    enddo

    ! External forces containing bouyancy force ONLY
    do iel = 1, ncel
      dronm1 = (croma(iel)-ro0)
      drom   = (crom(iel)-ro0)
      frchy(1,iel)  = dronm1*gx
      frchy(2,iel)  = dronm1*gy
      frchy(3,iel)  = dronm1*gz
      dfrchy(1,iel) = drom  *gx - frchy(1,iel)
      dfrchy(2,iel) = drom  *gy - frchy(2,iel)
      dfrchy(3,iel) = drom  *gz - frchy(3,iel)
    enddo

    ! Parallelism and periodicity treatment
    if (irangp.ge.0.or.iperio.eq.1) then
      call synvin(frchy)
      call synvin(dfrchy)
    endif

    call calhyd &
    !==========
    ( indhyd ,                                &
      !TODO
      !frchy, dfrchy,                          &
      frcxt  , dfrcxt ,                       &
      hydro_pres      , iflux  , bflux ,      &
      coefap , coefbp ,                       &
      cofafp , cofbfp ,                       &
      viscf  , viscb  ,                       &
      dam    , xam    ,                       &
      dphi   , rhs    ) !FIXME remove work arrays.

    ! Free memory
    deallocate(coefap, cofafp)
    deallocate(coefbp, cofbfp)

  endif

else

  indhyd = 0

endif

! Compute the BCs for the pressure increment
! (first we set the BCs of a standard pressure increment,
!  that are (A = 0, B_dp = B_p) for the gradient BCs
!  Then the A_dp is set thank to the pre-computed hydrostatic pressure
!  so that the pressure increment will be 0 on the reference outlet face.

if (iphydr.eq.1.or.iifren.eq.1) then

  if (indhyd.eq.1) then
    ifac0 = isostd(nfabor+1)
    if (ifac0.le.0) then
      phydr0 = 0.d0
    else
      iel0 = ifabor(ifac0)
      phydr0 = hydro_pres(iel0)                                     &
           +(cdgfbo(1,ifac0)-xyzcen(1,iel0))*dfrcxt(1 ,iel0) &
           +(cdgfbo(2,ifac0)-xyzcen(2,iel0))*dfrcxt(2 ,iel0) &
           +(cdgfbo(3,ifac0)-xyzcen(3,iel0))*dfrcxt(3 ,iel0)
    endif

    if (irangp.ge.0) then
      call parsom (phydr0)
    endif
  endif

  ! If hydrostatic pressure increment or free entrance Inlet
  if (indhyd.eq.1.or.iifren.eq.1) then

    do ifac = 1, nfabor
      iautof = 0
      ! automatic inlet/outlet face for atmospheric flow
      if (imeteo.gt.0) then
        iautof = iautom(ifac)
      endif

      if (isostd(ifac).eq.1.or.iatmst.eq.1.and.iautof.eq.1) then
        iel=ifabor(ifac)

        if (indhyd.eq.1) then
          coefa_dp(ifac) =  hydro_pres(iel)                               &
                         + (cdgfbo(1,ifac)-xyzcen(1,iel))*dfrcxt(1 ,iel)  &
                         + (cdgfbo(2,ifac)-xyzcen(2,iel))*dfrcxt(2 ,iel)  &
                         + (cdgfbo(3,ifac)-xyzcen(3,iel))*dfrcxt(3 ,iel)  &
                         -  phydr0
        endif

        ! Diffusive flux BCs
        if (iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then
          hint = dt(iel)/distb(ifac)

        ! Symmetric tensor diffusivity
        elseif (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then

          visci(1,1) = vitenp(1,iel)
          visci(2,2) = vitenp(2,iel)
          visci(3,3) = vitenp(3,iel)
          visci(1,2) = vitenp(4,iel)
          visci(2,1) = vitenp(4,iel)
          visci(2,3) = vitenp(5,iel)
          visci(3,2) = vitenp(5,iel)
          visci(1,3) = vitenp(6,iel)
          visci(3,1) = vitenp(6,iel)

          ! ||Ki.S||^2
          viscis = ( visci(1,1)*surfbo(1,ifac)       &
                   + visci(1,2)*surfbo(2,ifac)       &
                   + visci(1,3)*surfbo(3,ifac))**2   &
                 + ( visci(2,1)*surfbo(1,ifac)       &
                   + visci(2,2)*surfbo(2,ifac)       &
                   + visci(2,3)*surfbo(3,ifac))**2   &
                 + ( visci(3,1)*surfbo(1,ifac)       &
                   + visci(3,2)*surfbo(2,ifac)       &
                   + visci(3,3)*surfbo(3,ifac))**2

          ! IF.Ki.S
          fikis = ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,1)   &
                  + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,1)   &
                  + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,1)   &
                  )*surfbo(1,ifac)                              &
                + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,2)   &
                  + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,2)   &
                  + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,2)   &
                  )*surfbo(2,ifac)                              &
                + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,3)   &
                  + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,3)   &
                  + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,3)   &
                  )*surfbo(3,ifac)

          distfi = distb(ifac)

          ! Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
          ! NB: eps =1.d-1 must be consistent with vitens.f90
          fikis = max(fikis, 1.d-1*sqrt(viscis)*distfi)

          hint = viscis/surfbn(ifac)/fikis

        endif

        ! Free entrance boundary face (Bernoulli condition to link the pressure
        ! increment and the predicted velocity)
        if (itypfb(ifac).eq.ifrent) then

          ! Boundary mass flux of the predicted velocity
          bpmasf = vel(1, iel)*surfbo(1, ifac)    &
                 + vel(2, iel)*surfbo(2, ifac)    &
                 + vel(3, iel)*surfbo(3, ifac)

          ! Ingoing mass Flux, Bernoulli relation ship is used
          if (bpmasf.le.0.d0) then

            ! Head loss of the fluid outside the domain, between infinity and
            ! the entrance
            kpdc = b_head_loss(ifac)
            rho = brom(ifac)
            cfl = -(bmasfl(ifac)/surfbn(ifac)*dt(iel))    &
                / (2.d0*rho*distb(ifac))*(1.d0 + kpdc)

            pimp = - cvar_pr(iel)                                              &
                 - 0.5d0*(1.d0 + kpdc)*bmasfl(ifac)*bpmasf/surfbn(ifac)**2

            call set_convective_outlet_scalar &
                 !==================
               ( coefa_dp(ifac), coefaf_dp(ifac),             &
                 coefb_dp(ifac), coefbf_dp(ifac),             &
                 pimp        , cfl         , hint )

          else
            coefaf_dp(ifac) = - hint*coefa_dp(ifac)
          endif

        else
          coefaf_dp(ifac) = - hint*coefa_dp(ifac)
        endif

      endif
    enddo
  endif

endif

!===============================================================================
! 3. Building of the linear system to solve
!===============================================================================

! ---> Implicit term

do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo
! Implicit part of the cavitation source
if (icavit.gt.0.and.itscvi.eq.1) then
  do iel = 1, ncel
    rovsdt(iel) = rovsdt(iel) - cell_f_vol(iel)*dgdpca(iel)*(1.d0/rho2 - 1.d0/rho1)
  enddo
endif
! Strengthen the diagonal for Low Mach Algorithm
if (idilat.eq.3) then
  do iel = 1, ncel
    rovsdt(iel) = rovsdt(iel) + epsdp*cell_f_vol(iel)/dt(iel)
  enddo
endif

! ---> Face diffusivity
if (vcopt_p%idiff.ge.1) then

  ! Scalar diffusivity
  if (iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then

    if (ivofmt.ge.0) then
      imvisp = 1  ! VOF algorithm: continuity of the flux across internal faces
    else
      imvisp = imvisf
    endif

    call viscfa &
    !==========
   ( imvisp ,            &
     viscap ,            &
     viscf  , viscb  )

    if (vcopt_p%iwgrec.eq.1) then
      ! Weighting for gradient
      do iel = 1, ncel
        cpro_wgrec_s(iel) = viscap(iel)
      enddo
      call synsca(cpro_wgrec_s)
    endif

  ! Tensor diffusivity
  else if (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then

    ! Allocate temporary arrays
    allocate(weighf(2,nfac))
    allocate(weighb(ndimfb))

    iwarnp = vcopt_p%iwarni

    call vitens &
    !==========
   ( vitenp , iwarnp ,             &
     weighf , weighb ,             &
     viscf  , viscb  )

    if (vcopt_p%iwgrec.eq.1) then
      ! Weighting for gradient
      do iel = 1, ncel
        do isou = 1, 6
          cpro_wgrec_v(isou,iel) = vitenp(isou,iel)
        enddo
      enddo
      call syntis(cpro_wgrec_v)
    endif

  endif

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

iconvp = vcopt_p%iconv
idiffp = vcopt_p%idiff
ndircp = vcopt_p%ndircl

thetap = 1.d0
imucpp = 0

call matrix &
!==========
 ( iconvp , idiffp , ndircp , isym   ,                            &
   thetap , imucpp ,                                              &
   coefb_dp , coefbf_dp     , rovsdt ,                            &
   imasfl , bmasfl , viscf  , viscb  ,                            &
   rvoid  , dam    , xam    )

!===============================================================================
! 4. Mass flux initialization
!===============================================================================

! --- Flux de masse predit et premiere composante Rhie et Chow

! Allocate a work array for the gradient calculation
allocate(gradp(3,ncelet))

iccocg = 1
iprev  = 0
inc    = 1
call grdpor(inc)

! Pressure gradient
! NB: for the VOF algo. the weighting is automatically done
! thank to the iwgrec variable calculation option.
call field_gradient_potential(ivarfl(ipr), iprev, imrgra, inc,    &
                              iccocg, iphydr,                     &
                              frcxt, gradp)

if (iphydr.eq.1.and.iporos.eq.3) then

  call field_get_val_v_by_name("poro_div_duq", cpro_poro_div_duq)
  do iel = 1, ncel
    do isou = 1, 3
      trav(isou, iel) = gradp(isou, iel) - frcxt(isou, iel) &
        - cpro_poro_div_duq(isou, iel)
    enddo
  enddo

else if (iphydr.eq.1) then
  do iel = 1, ncel
    do isou = 1, 3
      trav(isou, iel) = gradp(isou, iel) - frcxt(isou, iel)
    enddo
  enddo
else
  do iel = 1, ncel
    do isou = 1, 3
      trav(isou,iel) = gradp(isou,iel)
    enddo
  enddo
endif

! --- Weakly compressible algorithm: semi analytic scheme
!     The RHS contains rho div(u*) and not div(rho u*)
!     so this term will be add afterwards
if (idilat.ge.4) then
  if (arak.gt.0.d0 .and. iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then
    do iel = 1, ncel
      ardtsr  = arak*(dt(iel)/crom(iel))
      do isou = 1, 3
        trav(isou,iel) = ardtsr*trav(isou,iel)
      enddo
    enddo
  else if (arak.gt.0.d0 .and. iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then
    do iel = 1, ncel
      arsr  = arak/crom(iel)

      trav(1,iel) = arsr*(                                 &
                           vitenp(1,iel)*trav(1,iel)      &
                         + vitenp(4,iel)*trav(2,iel)      &
                         + vitenp(6,iel)*trav(3,iel)      &
                         )
      trav(2,iel) = arsr*(                                 &
                           vitenp(4,iel)*trav(1,iel)      &
                         + vitenp(2,iel)*trav(2,iel)      &
                         + vitenp(5,iel)*trav(3,iel)      &
                         )
      trav(3,iel) = arsr*(                                 &
                           vitenp(6,iel)*trav(1,iel)      &
                         + vitenp(5,iel)*trav(2,iel)      &
                         + vitenp(3,iel)*trav(3,iel)      &
                         )

    enddo
  else
    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou = 1, 3
        trav(isou,iel) = 0.d0
      enddo
    enddo
  endif

! Standard algorithm
else
  if (arak.gt.0.d0 .and. iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then
    do iel = 1, ncel
      ardtsr  = arak*(dt(iel)/crom(iel))
      do isou = 1, 3
        trav(isou,iel) = vel(isou,iel) + ardtsr*trav(isou,iel)
      enddo
    enddo
  else if (arak.gt.0.d0 .and. iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then
    do iel = 1, ncel
      arsr  = arak/crom(iel)

      trav(1,iel) = vel(1,iel) + arsr*(                   &
                           vitenp(1,iel)*trav(1,iel)      &
                         + vitenp(4,iel)*trav(2,iel)      &
                         + vitenp(6,iel)*trav(3,iel)      &
                         )
      trav(2,iel) = vel(2,iel) + arsr*(                   &
                           vitenp(4,iel)*trav(1,iel)      &
                         + vitenp(2,iel)*trav(2,iel)      &
                         + vitenp(5,iel)*trav(3,iel)      &
                         )
      trav(3,iel) = vel(3,iel) + arsr*(                   &
                           vitenp(6,iel)*trav(1,iel)      &
                         + vitenp(5,iel)*trav(2,iel)      &
                         + vitenp(3,iel)*trav(3,iel)      &
                         )

    enddo
  else
    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou = 1, 3
        trav(isou, iel) = vel(isou, iel)
      enddo
    enddo
  endif
endif

! ---> Traitement du parallelisme et de la periodicite

if (irangp.ge.0.or.iperio.eq.1) then
  call synvin(trav)
endif

init   = 1
inc    = 1
! BCs will be taken into account after in idilat>=4
if (idilat.ge.4) inc = 0
iflmb0 = 1
if (iale.ge.1) iflmb0 = 0
nswrgp = vcopt_u%nswrgr
imligp = vcopt_u%imligr
iwarnp = vcopt_p%iwarni
epsrgp = vcopt_u%epsrgr
climgp = vcopt_u%climgr
itypfl = 1
if (ivofmt.ge.0.or.idilat.eq.4) itypfl = 0

call inimav &
!==========
 ( f_id0  , itypfl ,                                              &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom   , brom   ,                                              &
   trav   ,                                                       &
   coefav , coefbv ,                                              &
   imasfl , bmasfl )


! --- Projection aux faces des forces exterieures

if (iphydr.eq.1) then
  init   = 0
  inc    = 0

  call grdpor(inc)

  iccocg = 1
  nswrgp = vcopt_p%nswrgr
  imligp = vcopt_p%imligr
  iwarnp = vcopt_p%iwarni
  epsrgp = vcopt_p%epsrgr
  climgp = vcopt_p%climgr
  ircflp = vcopt_p%ircflu

  ! Scalar diffusivity
  if (iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then
    call projts &
    !==========
 ( init   , nswrgp ,                                              &
   dfrcxt ,                                                       &
   coefbf_p ,                                                     &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  ,                                              &
   viscap , viscap , viscap     )

  ! Tensor diffusivity
  else if (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then

    call projtv &
    !==========
  ( init   , nswrgp , ircflp ,                                     &
    dfrcxt ,                                                       &
    coefbf_p ,                                                     &
    viscf  , viscb  ,                                              &
    vitenp ,                                                       &
    weighf ,                                                       &
    imasfl , bmasfl )

  endif
endif

init   = 0
inc    = 1
call grdpor(inc)
iccocg = 1

!----------------------
! Rhie and Chow filter
!----------------------
if (arak.gt.0.d0) then

  allocate(ipro_visc(nfac))
  allocate(bpro_visc(nfabor))

  ! --- Prise en compte de Arak: la viscosite face est multipliee
  !       Le pas de temps aussi.
  do ifac = 1, nfac
    ipro_visc(ifac) = arak*viscf(ifac)
  enddo
  do ifac = 1, nfabor
    bpro_visc(ifac) = arak*viscb(ifac)
  enddo

  ! On annule la viscosite facette pour les faces couplees pour ne pas modifier
  ! le flux de masse au bord dans le cas d'un dirichlet de pression: la correction
  ! de pression et le filtre sont annules.
  if (nbrcpl.gt.0) then
    do ifac = 1, nfabor
      if (itypfb(ifac).eq.icscpd) then
        bpro_visc(ifac) = 0.d0
      endif
    enddo
  endif

  ! Scalar diffusivity
  !-------------------
  if (iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then

    allocate(cpro_visc(ncelet))

    do iel = 1, ncel
      cpro_visc(iel) = arak*viscap(iel)
    enddo

    nswrgp = vcopt_p%nswrgr
    imligp = vcopt_p%imligr
    iwarnp = vcopt_p%iwarni
    epsrgp = vcopt_p%epsrgr
    climgp = vcopt_p%climgr
    extrap = vcopt_p%extrag

    call itrmas &
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr ,     &
   0      , iwarnp ,                                                           &
   epsrgp , climgp , extrap ,                                                  &
   frcxt  ,                                                                    &
   cvar_pr,                                                                    &
   coefa_p, coefb_p, coefaf_p        , coefbf_p        ,                       &
   ipro_visc      , bpro_visc      ,                                           &
   cpro_visc ,                                                                 &
   imasfl , bmasfl )

    ! Projection du terme source pour oter la partie hydrostat de la pression
    if (iphydr.eq.1) then
      init   = 0
      inc    = 0
      iccocg = 1
      nswrgp = vcopt_p%nswrgr
      imligp = vcopt_p%imligr
      iwarnp = vcopt_p%iwarni
      epsrgp = vcopt_p%epsrgr
      climgp = vcopt_p%climgr

      ! A 0 boundary coefficient coefbf_dp is passed to projts
      ! to cancel boundary terms
      allocate(cofbfp(ndimfb))
      do ifac = 1,nfabor
        cofbfp(ifac) = 0.d0
      enddo

      call projts &
      !==========
 ( init   , nswrgp ,                                              &
   frcxt  ,                                                       &
   cofbfp ,                                                       &
   imasfl , bmasfl ,                                              &
   ipro_visc       , bpro_visc  ,                                 &
   cpro_visc, cpro_visc, cpro_visc    )

      deallocate(cofbfp)

    endif

    deallocate(cpro_visc)

  ! Tensor diffusivity
  !-------------------
  else if (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then

    allocate(cpro_vitenp(6, ncelet))
    do iel = 1, ncel
      cpro_vitenp(1,iel) = arak*vitenp(1,iel)
      cpro_vitenp(2,iel) = arak*vitenp(2,iel)
      cpro_vitenp(3,iel) = arak*vitenp(3,iel)
      cpro_vitenp(4,iel) = arak*vitenp(4,iel)
      cpro_vitenp(5,iel) = arak*vitenp(5,iel)
      cpro_vitenp(6,iel) = arak*vitenp(6,iel)
    enddo

    nswrgp = vcopt_p%nswrgr
    imligp = vcopt_p%imligr
    iwarnp = vcopt_p%iwarni
    epsrgp = vcopt_p%epsrgr
    climgp = vcopt_p%climgr
    extrap = vcopt_p%extrag
    ircflp = vcopt_p%ircflu

    call itrmav &
    !==========
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
   iphydr , 0      , iwarnp ,                                              &
   epsrgp , climgp , extrap ,                                              &
   frcxt  ,                                                                &
   cvar_pr,                                                                &
   coefa_p , coefb_p , coefaf_p , coefbf_p ,                               &
   ipro_visc        , bpro_visc ,                                          &
   cpro_vitenp      ,                                                      &
   weighf , weighb ,                                                       &
   imasfl , bmasfl )

    ! Projection du terme source pour oter la partie hydrostat de la pression
    if (iphydr.eq.1) then
      init   = 0
      inc    = 0
      nswrgp = vcopt_p%nswrgr
      imligp = vcopt_p%imligr
      iwarnp = vcopt_p%iwarni
      epsrgp = vcopt_p%epsrgr
      climgp = vcopt_p%climgr

      ! A 0 boundary coefficient coefbf_dp is passed to projtv
      ! to cancel boundary terms
      allocate(cofbfp(ndimfb))
      do ifac = 1, nfabor
        cofbfp(ifac) = 0.d0
      enddo

      call projtv &
      !==========
   ( init   , nswrgp , ircflp ,                                     &
     frcxt  ,                                                       &
     cofbfp ,                                                       &
     ipro_visc  , bpro_visc ,                                       &
     cpro_vitenp ,                                                  &
     weighf ,                                                       &
     imasfl, bmasfl )

      deallocate(cofbfp)

    endif

    deallocate(cpro_vitenp)

  endif

  deallocate(ipro_visc)
  deallocate(bpro_visc)

endif

!===============================================================================
! 5. Solving (Loop over the non-orthogonalities)
!===============================================================================

! --- Number of sweeps
nswmpr = vcopt_p%nswrsm

! --- Variables are set to 0
!       phi        is the increment of the pressure
!       dphi       is the increment of the increment between sweeps
!       cpro_divu       is the initial divergence of the predicted mass flux

do iel = 1, ncel
  phi(iel)  = 0.d0
  dphi(iel) = 0.d0
  phia(iel) = 0.d0
enddo

relaxp = vcopt_p%relaxv

! --- Initial divergence
init = 1

call divmas(init, imasfl , bmasfl , cpro_divu)

! --- Weakly compressible algorithm: semi analytic scheme
!     1. The RHS contains rho div(u*) and not div(rho u*)
!     2. Add dilatation source term to rhs
!     3. The mass flux is completed by rho u* . S

if (idilat.ge.4) then

  call field_get_val_s(iustdy(itsrho), cpro_tsrho)

  allocate(velflx(nfac), velflb(ndimfb))

  ! 1. The RHS contains rho div(u*) and not div(rho u*)
  init   = 1
  inc    = 1
  iflmb0 = 1
  if (iale.ge.1) iflmb0 = 0
  nswrgp = vcopt_u%nswrgr
  imligp = vcopt_u%imligr
  iwarnp = vcopt_p%iwarni
  epsrgp = vcopt_u%epsrgr
  climgp = vcopt_u%climgr

  itypfl = 0

  call inimav &
  !==========
 ( ivarfl(iu)      , itypfl ,                                     &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom, brom   ,                                                 &
   vel    ,                                                       &
   coefav , coefbv ,                                              &
   velflx , velflb )

  call divmas(init, velflx , velflb , res)

  if (idilat.eq.4) then
    do iel = 1, ncel
      cpro_divu(iel) = cpro_divu(iel) + res(iel)
    enddo
  else
    do iel = 1, ncel
      cpro_divu(iel) = cpro_divu(iel) + res(iel)*crom(iel)
    enddo
  endif

  ! 2. Add the dilatation source term D(rho)/Dt
  if (idilat.eq.4) then
    do iel = 1, ncel
      cpro_divu(iel) = cpro_divu(iel) &
                + cpro_tsrho(iel)/crom(iel)
    enddo
  else
    do iel = 1, ncel
      cpro_divu(iel) = cpro_divu(iel) + cpro_tsrho(iel)
    enddo
  endif

  ! 3. The mass flux is completed by u*.S (idilat=4)
  !                                  rho u* . S (idilat=5)
  init   = 0
  inc    = 1
  iflmb0 = 1
  if (iale.ge.1) iflmb0 = 0
  nswrgp = vcopt_u%nswrgr
  imligp = vcopt_u%imligr
  iwarnp = vcopt_p%iwarni
  epsrgp = vcopt_u%epsrgr
  climgp = vcopt_u%climgr

  itypfl = 1
  if (idilat.eq.4) itypfl = 0

  call inimav &
  !==========
 ( ivarfl(iu)      , itypfl ,                                     &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom, brom   ,                                                 &
   vel    ,                                                       &
   coefav , coefbv ,                                              &
   imasfl , bmasfl )

endif

! --- Masse source terms adding for volumic flow rate
if (ncesmp.gt.0) then
  do ii = 1, ncesmp
    iel = icetsm(ii)
    cpro_divu(iel) = cpro_divu(iel) -cell_f_vol(iel)*smacel(ii,ipr)
  enddo
endif

! --- Source term adding for condensation modelling
if (nfbpcd.gt.0) then
  do ii = 1, nfbpcd
    ifac= ifbpcd(ii)
    iel = ifabor(ifac)
    cpro_divu(iel) = cpro_divu(iel) - surfbn(ifac)*spcond(ii,ipr)
  enddo
endif

! --- volume Gamma source for metal mass structures
!     condensation modelling
if (icondv.eq.0) then
  allocate(surfbm(ncelet))
  surfbm(:) = 0.d0

  do ii = 1, ncmast
    iel= ltmast(ii)
    surfbm(iel) = s_metal*volume(iel)/voltot
    cpro_divu(iel) = cpro_divu(iel) - surfbm(iel)*svcond(iel,ipr)
  enddo

  deallocate(surfbm)
endif


! --- Source term associated to the mass aggregation
if (idilat.eq.2.or.idilat.eq.3) then

  ! Add source term
  do iel = 1, ncel
    drom = crom_eos(iel) - croma(iel)
    cpro_divu(iel) = cpro_divu(iel) + drom*cell_f_vol(iel)/dt(iel)
  enddo

endif

! ---> Termes sources Lagrangien
if (iilagr.eq.2 .and. ltsmas.eq.1) then
  do iel = 1, ncel
    cpro_divu(iel) = cpro_divu(iel) -tslagr(iel,itsmas)
  enddo
endif

! --- Cavitation source term
if (icavit.gt.0) then
  do iel = 1, ncel
    cpro_divu(iel) = cpro_divu(iel) - cell_f_vol(iel)*gamcav(iel)*(1.d0/rho2 - 1.d0/rho1)
  enddo
endif

! --- Initial right hand side
do iel = 1, ncel
  rhs(iel) = - cpro_divu(iel) - rovsdt(iel)*phi(iel)
enddo

! --- Right hand side residual
residu = sqrt(cs_gdot(ncel,rhs,rhs))

sinfo%rnsmbr = residu

! --- Norm resiudal
! Historical norm for the pressure step:
!       div(rho u* + dt gradP^(n))-Gamma
!       i.e.  RHS of the pressure + div(dt gradP^n) (otherwise there is a risk
!       a 0 norm at steady states...). Represents terms that pressure has to
!       balance.

do iel = 1, ncel
  dtsrom = dt(iel) / crom(iel)
  do isou = 1, 3
    trav(isou, iel) = dtsrom * gradp(isou,iel)
  enddo
enddo

!---> Parallelism
if (irangp.ge.0.or.iperio.eq.1) then
  call synvin(trav)
endif

! To save time, no space reconstruction
init   = 1
inc    = 1
iflmb0 = 1
if (iale.ge.1) iflmb0 = 0
nswrgp = 1
imligp = vcopt_u%imligr
iwarnp = vcopt_p%iwarni
epsrgp = vcopt_u%epsrgr
climgp = vcopt_u%climgr
itypfl = 1
! VOF algorithm: the pressure step corresponds to the
! correction of the volumetric flux, not the mass flux
if (idilat.ge.4.or.ivofmt.ge.0) itypfl = 0

call inimav &
(f_id0  , itypfl ,                                              &
 iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
 iwarnp ,                                                       &
 epsrgp , climgp ,                                              &
 crom   , brom   ,                                              &
 trav   ,                                                       &
 coefav , coefbv ,                                              &
 iflux  , bflux  )

init = 1
call divmas(init, iflux, bflux, res)

! Free memory
deallocate(iflux, bflux)

! --- Weakly compressible algorithm: semi analytic scheme
if (idilat.ge.4) then
  do iel = 1, ncel
    res(iel) = res(iel)*crom(iel)
  enddo
endif

! It is: div(dt/rho*rho grad P) + div(rho u*) - Gamma
! NB: if iphydr=1, div(rho u*) contains div(d fext).
do iel = 1, ncel
  res(iel) = res(iel) + cpro_divu(iel)
enddo

! --- Pressure norm
rnormp = sqrt(cs_gdot(ncel,res,res))

if (vcopt_p%iwarni.ge.2) then
  write(nfecra,1300)chaine(1:16) ,rnormp
endif
if (iterns.le.1) then
  sinfo%nbivar = 0
endif

! Pressure derive for the log
if (rnormp.lt.epzero) then
  sinfo%dervar = - sinfo%rnsmbr
else
  sinfo%dervar = sinfo%rnsmbr/rnormp
endif

isweep = 1

! Writing
if (vcopt_p%iwarni.ge.2) then
  write(nfecra,1400)chaine(1:16),isweep,residu, relaxp
endif

! Dynamic relaxation initialization
!----------------------------------
if (iswdyp.ge.1) then

  do iel = 1, ncelet
    adxkm1(iel) = 0.d0
    adxk(iel) = 0.d0
  enddo

  ! ||A.dx^0||^2 = 0
  nadxk = 0.d0

  rnorm2 = rnormp**2

  iccocg = 1
  init = 1
  inc  = 0
  if (iphydr.eq.1.or.iifren.eq.1) inc = 1
  nswrgp = vcopt_p%nswrgr
  imligp = vcopt_p%imligr
  iwarnp = vcopt_p%iwarni
  epsrgp = vcopt_p%epsrgr
  climgp = vcopt_p%climgr
  extrap = vcopt_p%extrag
  ircflp = vcopt_p%ircflu

  if (iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then

    call itrgrp &
    !==========
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr ,     &
   iwarnp ,                                                                    &
   epsrgp , climgp , extrap ,                                                  &
   dfrcxt ,                                                                    &
   dphi   ,                                                                    &
   coefa_dp  , coefb_dp  ,                                                     &
   coefaf_dp , coefbf_dp ,                                                     &
   viscf  , viscb  ,                                                           &
   viscap ,                                                                    &
   rhs0   )

  else if (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then

    call itrgrv &
    !==========
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
   iphydr , iwarnp ,                                                       &
   epsrgp , climgp , extrap ,                                              &
   dfrcxt ,                                                                &
   dphi   ,                                                                &
   coefa_dp  , coefb_dp  ,                                                 &
   coefaf_dp , coefbf_dp ,                                                 &
   viscf  , viscb  ,                                                       &
   vitenp ,                                                                &
   weighf , weighb ,                                                       &
   rhs0   )

  endif

endif

! Reconstruction loop (beginning)
!--------------------------------

do while (isweep.le.nswmpr.and.residu.gt.vcopt_p%epsrsm*rnormp)

  ! Solving on the increment dphi
  !-------------------------------
  if (iswdyp.eq.0) then
    do iel = 1, ncel
      dphi(iel) = 0.d0
    enddo
  else
    do iel = 1, ncel
      dphim1(iel) = dphi(iel)
      dphi(iel) = 0.d0
    enddo
  endif

  iwarnp = vcopt_p%iwarni
  epsilp = vcopt_p%epsilo

  ! Solver resiudal
  ressol = residu

  call sles_solve_native(ivarfl(ipr), '',                            &
                         isym, ibsize, iesize, dam, xam,             &
                         epsilp, rnormp, niterf, ressol, rhs, dphi)

  ! Dynamic relaxation of the system
  !---------------------------------
  if (iswdyp.ge.1) then

    ! Computation of the variable ralaxation coefficient

    !$omp parallel do
    do iel = 1, ncelet
      adxkm1(iel) = adxk(iel)
      adxk(iel) = - rhs0(iel)
    enddo

    iccocg = 1
    init = 0
    inc  = 0
    if (iphydr.eq.1.or.iifren.eq.1) inc = 1
    nswrgp = vcopt_p%nswrgr
    imligp = vcopt_p%imligr
    iwarnp = vcopt_p%iwarni
    epsrgp = vcopt_p%epsrgr
    climgp = vcopt_p%climgr
    extrap = vcopt_p%extrag

    if (iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then

      call itrgrp &
      !==========
   ( f_id0           , init   , inc    , imrgra ,              &
     iccocg , nswrgp , imligp , iphydr ,                       &
     iwarnp ,                                                  &
     epsrgp , climgp , extrap ,                                &
     dfrcxt ,                                                  &
     dphi   ,                                                  &
     coefa_dp  , coefb_dp  ,                                   &
     coefaf_dp , coefbf_dp ,                                   &
     viscf  , viscb  ,                                         &
     viscap ,                                                  &
     adxk   )

    else if (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then

      ircflp = vcopt_p%ircflu

      call itrgrv &
      !==========
   ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
     iphydr , iwarnp ,                                                       &
     epsrgp , climgp , extrap ,                                              &
     dfrcxt ,                                                                &
     dphi   ,                                                                &
     coefa_dp  , coefb_dp  ,                                                 &
     coefaf_dp , coefbf_dp ,                                                 &
     viscf  , viscb  ,                                                       &
     vitenp ,                                                                &
     weighf , weighb ,                                                       &
     adxk   )

    endif

    do iel = 1, ncel
      adxk(iel) = - adxk(iel)
    enddo

    ! ||E.dx^(k-1)-E.0||^2
    nadxkm1 = nadxk

    ! ||E.dx^k-E.0||^2
    nadxk = cs_gdot(ncel, adxk, adxk)

    ! < E.dx^k-E.0; r^k >
    paxkrk = cs_gdot(ncel, rhs, adxk)

    ! Relaxation with respect to dx^k and dx^(k-1)
    if (iswdyp.ge.2) then

      ! < E.dx^(k-1)-E.0; r^k >
      paxm1rk = cs_gdot(ncel, rhs, adxkm1)

      ! < E.dx^(k-1)-E.0; E.dx^k -E.0 >
      paxm1ax = cs_gdot(ncel, adxk, adxkm1)

      if (nadxkm1.gt.1.d-30*rnorm2.and.                    &
         (nadxk*nadxkm1-paxm1ax**2).gt.1.d-30*rnorm2) then
        beta = (paxkrk*paxm1ax - nadxk*paxm1rk)/(nadxk*nadxkm1-paxm1ax**2)
      else
        beta = 0.d0
      endif

    else
      beta = 0.d0
      paxm1ax = 1.d0
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
      write(nfecra,1200) chaine(1:16), isweep, alph, beta, &
                         paxkrk, nadxk, paxm1rk, nadxkm1, paxm1ax
    endif

  endif

  ! Update the increment of pressure
  !---------------------------------

  if (iswdyp.eq.0) then
    if (idtvar.ge.0.and.isweep.le.nswmpr.and.residu.gt.vcopt_p%epsrsm*rnormp) then
      do iel = 1, ncel
        phia(iel) = phi(iel)
        phi(iel) = phi(iel) + vcopt_p%relaxv*dphi(iel)
      enddo
    ! If it is the last sweep, update with the total increment
    else
      do iel = 1, ncel
        phia(iel) = phi(iel)
        phi(iel) = phi(iel) + dphi(iel)
      enddo
    endif
  elseif (iswdyp.eq.1) then
     do iel = 1, ncel
      phia(iel) = phi(iel)
      phi(iel) = phi(iel) + alph*dphi(iel)
    enddo
  elseif (iswdyp.ge.2) then
    do iel = 1, ncel
      phia(iel) = phi(iel)
      phi(iel) = phi(iel) + alph*dphi(iel) + beta*dphim1(iel)
    enddo
  endif

  ! --- Update the right hand side and update the residual
  !      rhs^{k+1} = - div(rho u^n) - D(dt, delta delta p^{k+1})
  !-------------------------------------------------------------

  iccocg = 1
  init = 1
  inc  = 0
  if (iphydr.eq.1.or.iifren.eq.1) inc = 1
  nswrgp = vcopt_p%nswrgr
  imligp = vcopt_p%imligr
  iwarnp = vcopt_p%iwarni
  epsrgp = vcopt_p%epsrgr
  climgp = vcopt_p%climgr
  extrap = vcopt_p%extrag

  if (iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then

    call itrgrp &
    !==========
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr ,     &
   iwarnp ,                                                                    &
   epsrgp , climgp , extrap ,                                                  &
   dfrcxt ,                                                                    &
   phi    ,                                                                    &
   coefa_dp  , coefb_dp  ,                                                     &
   coefaf_dp , coefbf_dp ,                                                     &
   viscf  , viscb  ,                                                           &
   viscap ,                                                                    &
   rhs    )

  else if (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then

    ircflp = vcopt_p%ircflu

    call itrgrv &
    !==========
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
   iphydr , iwarnp ,                                                       &
   epsrgp , climgp , extrap ,                                              &
   dfrcxt ,                                                                &
   phi    ,                                                                &
   coefa_dp  , coefb_dp  ,                                                 &
   coefaf_dp , coefbf_dp ,                                                 &
   viscf  , viscb  ,                                                       &
   vitenp ,                                                                &
   weighf , weighb ,                                                       &
   rhs    )

  endif

  do iel = 1, ncel
    rhs(iel) = - cpro_divu(iel) - rhs(iel) - rovsdt(iel)*phi(iel)
  enddo

  ! --- Convergence test
  residu = sqrt(cs_gdot(ncel,rhs,rhs))

  ! Writing
  sinfo%nbivar = sinfo%nbivar + niterf

  ! Writing
  if (vcopt_p%iwarni.ge.2) then
    write(nfecra,1400) chaine(1:16), isweep, residu, relaxp
    write(nfecra,1500) chaine(1:16), isweep, residu, rnormp, niterf
  endif

  isweep = isweep + 1

enddo
! --- Reconstruction loop (end)

! For logging
if (abs(rnormp).gt.0.d0) then
  sinfo%resvar = residu/rnormp
else
  sinfo%resvar = 0.d0
endif

! Writing
if(vcopt_p%iwarni.ge.1) then
  if (residu.le.vcopt_p%epsrsm*rnormp) then
    write(nfecra,1500) chaine(1:16), isweep-1, residu, rnormp, niterf

  ! Writing: non-convergence
  else if(isweep.gt.nswmpr) then
    write(nfecra,1600) chaine(1:16),nswmpr
  endif
endif

! Save convergence info
call field_set_key_struct_solving_info(ivarfl(ipr), sinfo)

! --- Compute the indicator, taken the volume into account (L2 norm)
!     or not
if(iescal(iesder).gt.0) then
  call field_get_val_s(iestim(iesder), c_estim_der)
  do iel = 1, ncel
    c_estim_der(iel) = abs(rhs(iel))/volume(iel)
  enddo
  if(iescal(iesder).eq.2) then
    do iel = 1, ncel
      c_estim_der(iel) = c_estim_der(iel)*sqrt(volume(iel))
    enddo
  endif
endif

! Update the mass flux
!---------------------

! On annule la viscosite facette pour les faces couplees pour ne pas modifier
! le flux de masse au bord dans le cas d'un dirichlet de pression: la correction
! de pression et le filtre sont annules.
if (nbrcpl.ge.0) then
  do ifac = 1, nfabor
    if (itypfb(ifac).eq.icscpd) then
      viscb(ifac) = 0.d0
    endif
  enddo
endif

iccocg = 1
init = 0
inc  = 0
! In case of hydrostatic pressure, inc is set to 1 to take explicit
! boundary conditions on the pressure (coefa)
if (iphydr.eq.1.or.iifren.eq.1) inc = 1
nswrgp = vcopt_p%nswrgr
imligp = vcopt_p%imligr
iwarnp = vcopt_p%iwarni
epsrgp = vcopt_p%epsrgr
climgp = vcopt_p%climgr
extrap = vcopt_p%extrag

if (iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then

  call itrmas &
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr ,     &
   0      , iwarnp ,                                                           &
   epsrgp , climgp , extrap ,                                                  &
   dfrcxt ,                                                                    &
   phia   ,                                                                    &
   coefa_dp  , coefb_dp  ,                                                     &
   coefaf_dp , coefbf_dp ,                                                     &
   viscf  , viscb  ,                                                           &
   viscap ,                                                                    &
   imasfl , bmasfl )

  ! The last increment is not reconstructed to fullfill exactly the continuity
  ! equation (see theory guide). The value of dfrcxt has no importance.
  iccocg = 0
  nswrgp = 0
  inc = 0

  call itrmas &
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr ,     &
   0      , iwarnp ,                                                           &
   epsrgp , climgp , extrap ,                                                  &
   dfrcxt ,                                                                    &
   dphi   ,                                                                    &
   coefa_dp  , coefb_dp  ,                                                     &
   coefaf_dp , coefbf_dp ,                                                     &
   viscf  , viscb  ,                                                           &
   viscap ,                                                                    &
   imasfl , bmasfl )

else if (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then

  ircflp = vcopt_p%ircflu

  call itrmav &
  !==========
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
   iphydr , 0      , iwarnp ,                                              &
   epsrgp , climgp , extrap ,                                              &
   dfrcxt ,                                                                &
   phia   ,                                                                &
   coefa_dp  , coefb_dp  ,                                                 &
   coefaf_dp , coefbf_dp ,                                                 &
   viscf  , viscb  ,                                                       &
   vitenp ,                                                                &
   weighf , weighb ,                                                       &
   imasfl , bmasfl )

  ! The last increment is not reconstructed to fullfill exactly the continuity
  ! equation (see theory guide). The value of dfrcxt has no importance.
  iccocg = 0
  nswrgp = 0
  inc = 0
  ircflp = 0

  call itrmav &
  !==========
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
   iphydr , 0      , iwarnp ,                                              &
   epsrgp , climgp , extrap ,                                              &
   dfrcxt ,                                                                &
   dphi   ,                                                                &
   coefa_dp  , coefb_dp  ,                                                 &
   coefaf_dp , coefbf_dp ,                                                 &
   viscf  , viscb  ,                                                       &
   vitenp ,                                                                &
   weighf , weighb ,                                                       &
   imasfl , bmasfl )

endif

!===============================================================================
! 6. Suppression of the mesh hierarchy
!===============================================================================

call sles_free_native(ivarfl(ipr), '')

!===============================================================================
! 7. Weakly compressible algorithm: semi analytic scheme
!    2nd step solving a convection diffusion equation
!===============================================================================

if (idilat.eq.5) then

  ! Allocate temporary arrays
  allocate(ddphi(ncelet))
  allocate(coefar(3,ndimfb), cofafr(3,ndimfb))
  allocate(coefbr(3,3,ndimfb), cofbfr(3,3,ndimfb))

  ! --- Convective flux: dt/rho grad(rho)
  inc = 1
  iccocg = 1
  nswrgp = vcopt_u%nswrgr
  imligp = vcopt_u%imligr
  iwarnp = vcopt_p%iwarni
  epsrgp = vcopt_u%epsrgr
  climgp = vcopt_u%climgr
  extrap = vcopt_u%extrag

  ! Dirichlet Boundary Condition on rho
  !------------------------------------

  allocate(coefa_rho(ndimfb), coefb_rho(ndimfb))

  do ifac = 1, nfabor
    coefa_rho(ifac) = brom(ifac)
    coefb_rho(ifac) = 0.d0
  enddo

  call gradient_s                                                 &
  (f_id0  , imrgra , inc    , iccocg , nswrgp , imligp , iwarnp , &
   epsrgp , climgp , extrap ,                                     &
   crom   ,                                                       &
   coefa_rho       , coefb_rho       ,                            &
   gradp  )

  deallocate(coefa_rho, coefb_rho)

  ! --- dt/rho * grad rho
  do iel = 1, ncel
    do isou = 1, 3
      trav(isou,iel) = gradp(isou,iel) * dt(iel) / crom(iel)
    enddo
  enddo

  ! --- (dt/rho * grad rho) . S

  init   = 1
  inc    = 1
  iflmb0 = 1
  if (iale.ge.1) iflmb0 = 0
  nswrgp = vcopt_u%nswrgr
  imligp = vcopt_u%imligr
  iwarnp = vcopt_p%iwarni
  epsrgp = vcopt_u%epsrgr
  climgp = vcopt_u%climgr
  extrap = vcopt_u%extrag

  itypfl = 0

  ! --- Viscosity
  call viscfa (imvisf, dt, viscf, viscb)

  ! --- Boundary Conditions for the convective flux
  do ifac = 1, nfabor

    iel = ifabor(ifac)

    ! Neumann Boundary Conditions
    !----------------------------

    qimpv(1) = 0.d0
    qimpv(2) = 0.d0
    qimpv(3) = 0.d0

    if (iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then
      hint = dt(iel)/distb(ifac)

    ! Symmetric tensor diffusivity
    else if (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then

      visci(1,1) = vitenp(1,iel)
      visci(2,2) = vitenp(2,iel)
      visci(3,3) = vitenp(3,iel)
      visci(1,2) = vitenp(4,iel)
      visci(2,1) = vitenp(4,iel)
      visci(2,3) = vitenp(5,iel)
      visci(3,2) = vitenp(5,iel)
      visci(1,3) = vitenp(6,iel)
      visci(3,1) = vitenp(6,iel)

      ! ||Ki.S||^2
      viscis = ( visci(1,1)*surfbo(1,ifac)       &
               + visci(1,2)*surfbo(2,ifac)       &
               + visci(1,3)*surfbo(3,ifac))**2   &
             + ( visci(2,1)*surfbo(1,ifac)       &
               + visci(2,2)*surfbo(2,ifac)       &
               + visci(2,3)*surfbo(3,ifac))**2   &
             + ( visci(3,1)*surfbo(1,ifac)       &
               + visci(3,2)*surfbo(2,ifac)       &
               + visci(3,3)*surfbo(3,ifac))**2

      ! IF.Ki.S
      fikis = ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,1)   &
              + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,1)   &
              + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,1)   &
              )*surfbo(1,ifac)                              &
            + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,2)   &
              + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,2)   &
              + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,2)   &
              )*surfbo(2,ifac)                              &
            + ( (cdgfbo(1,ifac)-xyzcen(1,iel))*visci(1,3)   &
              + (cdgfbo(2,ifac)-xyzcen(2,iel))*visci(2,3)   &
              + (cdgfbo(3,ifac)-xyzcen(3,iel))*visci(3,3)   &
              )*surfbo(3,ifac)

      distfi = distb(ifac)

      ! Take I" so that I"F= eps*||FI||*Ki.n when J" is in cell rji
      ! NB: eps =1.d-1 must be consistent with vitens.f90
      fikis = max(fikis, 1.d-1*sqrt(viscis)*distfi)

      hint = viscis/surfbn(ifac)/fikis

    endif

    call set_neumann_vector &
         !=================
       ( coefar(1,ifac)  , cofafr(1,ifac)  ,             &
         coefbr(1,1,ifac), cofbfr(1,1,ifac),             &
         qimpv           , hint             )

  enddo

  call inimav &
  !==========
 ( f_id0  , itypfl ,                                              &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom, brom   ,                                                 &
   trav   ,                                                       &
   coefar , coefbr ,                                              &
   velflx , velflb )

  ! --- Boundary condition for the pressure increment
  ! coefb, coefbf are those of the pressure
  allocate(coefa_dp2(ndimfb), coefaf_dp2(ndimfb))

  do ifac = 1, nfabor
   coefa_dp2(ifac) = 0.d0
   coefaf_dp2(ifac) = 0.d0
  enddo

  ! --- Convective source term
  do iel = 1, ncel
    rhs(iel) = 0.d0
  enddo

  ivar   = ipr
  f_id0  = -1
  iconvp = 1
  imasac = 1
  idiffp = 0
  nswrsp = 1
  imligp = vcopt_p%imligr
  ircflp = vcopt_p%ircflu
  ischcp = vcopt_p%ischcv
  isstpp = vcopt_p%isstpc
  inc    = 1
  iccocg = 1
  iwarnp = vcopt_p%iwarni
  imucpp = 0
  idftnp = vcopt_p%idften
  blencp = vcopt_p%blencv
  epsrgp = vcopt_p%epsrgr
  climgp = vcopt_p%climgr
  extrap = vcopt_p%extrag
  relaxp = vcopt_p%relaxv
  thetap = 1.d0
  ! all boundary convective flux with upwind
  icvflb = 0

  call bilsca &
  !==========
 ( idtvar , f_id0  , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   iwarnp , imucpp , idftnp , imasac ,                            &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   phi    , phi    ,                                              &
   coefa_dp2       , coefb_p, coefaf_dp2      , coefbf_p,         &
   velflx , velflb , viscf  , viscb  , rvoid  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   icvflb , ivoid  ,                                              &
   rhs   )

  ! --- Initialization of the variable to solve
  do iel = 1, ncel
    rovsdt(iel) = 340.d0/dt(iel) * cell_f_vol(iel)
    dphi(iel)   = 0.d0
    ddphi(iel)  = 0.d0
    rhs(iel)    = - rhs(iel)
  enddo

  ! --- Solve the convection diffusion equation

  idiffp = 1
  ! To reinforce the diagonal
  ndircp = 0
  nswrsp = vcopt_p%nswrsm
  nswrgp = vcopt_p%nswrgr
  imligp = vcopt_p%imligr
  ircflp = vcopt_p%ircflu
  ischcp = vcopt_p%ischcv
  isstpp = vcopt_p%isstpc
  iescap = 0
  imucpp = 0
  idftnp = vcopt_p%idften
  iswdyp = vcopt_p%iswdyn
  iwarnp = vcopt_p%iwarni
  blencp = vcopt_p%blencv
  epsilp = vcopt_p%epsilo
  epsrsp = vcopt_p%epsrsm
  epsrgp = vcopt_p%epsrgr
  climgp = vcopt_p%climgr
  extrap = vcopt_p%extrag
  relaxp = vcopt_p%relaxv
  thetap = vcopt_p%thetav
  ! all boundary convective flux with upwind
  icvflb = 0
  normp = -1.d0
  ! ivar = 0
  nomva0 = "Pr compress"

  ! --- Solve the convection diffusion equation

  call sles_push(ivarfl(ipr), "Pr compress")

  call codits &
  !==========
   ( idtvar , iterns , ivarfl(ivar)    , iconvp , idiffp , ndircp , &
     imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
     ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
     iwarnp , normp  ,                                              &
     blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
     relaxp , thetap ,                                              &
     dphi   , dphi   ,                                              &
     coefa_dp2       , coefb_p, coefaf_dp2      ,coefbf_p,          &
     velflx , velflb ,                                              &
     viscf  , viscb  , viscf  , viscb  , rvoid  ,                   &
     weighf , weighb ,                                              &
     icvflb , ivoid  ,                                              &
     rovsdt , rhs    , dphi   , ddphi  ,                            &
     rvoid  , rvoid  )

  call sles_pop(ivarfl(ipr))

  ! --- Update the increment of Pressure

  do iel = 1, ncel
    phi(iel) = phi(iel) + dphi(iel)
    ! Remove the last increment
    dphi(iel) = dphi(iel) - ddphi(iel)
  enddo

  ! --- Update the Mass flux

  init   = 0
  inc    = 1
  iccocg = 1
  nswrgp = vcopt_p%nswrgr
  imligp = vcopt_p%imligr
  iwarnp = vcopt_p%iwarni
  epsrgp = vcopt_p%epsrgr
  climgp = vcopt_p%climgr
  extrap = vcopt_p%extrag

  if (iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then
    call itrmas &
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr ,     &
   0      , iwarnp ,                                                           &
   epsrgp , climgp , extrap ,                                                  &
   dfrcxt ,                                                                    &
   dphi   ,                                                                    &
   coefa_dp2       , coefb_p, coefaf_dp2      ,coefbf_p,                       &
   viscf  , viscb  ,                                                           &
   dt     ,                                                                    &
   imasfl , bmasfl )

    ! The last increment is not reconstructed to fullfill exactly the continuity
    ! equation (see theory guide). The value of dfrcxt has no importance.
    iccocg = 0
    nswrgp = 0
    inc = 0

    call itrmas &
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydr ,     &
   0      , iwarnp ,                                                           &
   epsrgp , climgp , extrap ,                                                  &
   dfrcxt ,                                                                    &
   ddphi  ,                                                                    &
   coefa_dp2       , coefb_p, coefaf_dp2      ,coefbf_p,                       &
   viscf  , viscb  ,                                                           &
   dt     ,                                                                    &
   imasfl , bmasfl )

  else if (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then

    call itrmav &
   ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
     iphydr , 0      , iwarnp ,                                              &
     epsrgp , climgp , extrap ,                                              &
     dfrcxt ,                                                                &
     dphi   ,                                                                &
     coefa_dp2       , coefb_p, coefaf_dp2      ,coefbf_p,                   &
     viscf  , viscb  ,                                                       &
     vitenp ,                                                                &
     weighf , weighb ,                                                       &
     imasfl , bmasfl )

    ! The last increment is not reconstructed to fullfill exactly the continuity
    ! equation (see theory guide). The value of dfrcxt has no importance.
    iccocg = 0
    nswrgp = 0
    inc = 0

    call itrmav &
   ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , ircflp , &
     iphydr , 0      , iwarnp ,                                              &
     epsrgp , climgp , extrap ,                                              &
     dfrcxt ,                                                                &
     ddphi  ,                                                                &
     coefa_dp2       , coefb_p, coefaf_dp2      ,coefbf_p,                   &
     viscf  , viscb  ,                                                       &
     vitenp ,                                                                &
     weighf , weighb ,                                                       &
     imasfl , bmasfl )

  endif

  ! Free memory
  deallocate(ddphi)
  deallocate(coefa_dp2, coefaf_dp2)
  deallocate(coefar, coefbr)
  deallocate(cofafr, cofbfr)
  deallocate(velflx, velflb)

endif

!===============================================================================
! 8. Update the pressure field
!===============================================================================

if (idtvar.lt.0) then
  do iel = 1, ncel
    cvar_pr(iel) = cvar_pr(iel) + vcopt_p%relaxv*phi(iel)
  enddo
else
  do iel = 1, ncel
    cvar_pr(iel) = cvar_pr(iel) + phi(iel)
  enddo
endif

! Transformation of volumic mass fluxes into massic mass fluxes
if (idilat.eq.4) then

  do ifac = 1, nfabor
    bmasfl(ifac) = bmasfl(ifac) * brom(ifac)
  enddo

  do ifac = 1, nfac
    ii = ifacel(1, ifac)
    jj = ifacel(2, ifac)
    ! FIXME: should be coherent with the convective scheme of the species...
    rho = pond(ifac)*crom(ii)+(1.d0-pond(ifac))*crom(jj)
    imasfl(ifac) = imasfl(ifac) * rho
  enddo

endif

! Free memory
deallocate(dam, xam)
deallocate(trav)
deallocate(res, phia, dphi)
if (allocated(divu)) deallocate(divu)
deallocate(gradp)
deallocate(coefaf_dp, coefbf_dp)
deallocate(rhs, rovsdt)
if (allocated(weighf)) deallocate(weighf, weighb)
if (iswdyp.ge.1) deallocate(adxk, adxkm1, dphim1, rhs0)
if (icalhy.eq.1) deallocate(frchy, dfrchy, hydro_pres)
if (ivofmt.ge.0.or.idilat.eq.4) then
  if (allocated(xdtsro)) deallocate(xdtsro)
  if (allocated(xunsro)) deallocate(xunsro)
  if (allocated(tpusro)) deallocate(tpusro)
endif
if (allocated(cpro_rho_tc)) deallocate(cpro_rho_tc)
if (allocated(bpro_rho_tc)) deallocate(bpro_rho_tc)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1200 format ( &
 1X,A16,' Sweep: ',I5,' Dynamic relaxation: alpha = ',E12.5,' beta = ',E12.5,/,&
'    < dI^k  ; R^k > = ',E12.5,' ||dI^k  ||^2 = ',E12.5                     ,/,&
'    < dI^k-1; R^k > = ',E12.5,' ||dI^k-1||^2 = ',E12.5                     ,/,&
'   < dI^k-1; dI^k > = ',E12.5)
 1300 format(1X,A16,' : RESIDU DE NORMALISATION =', E14.6)
 1400 format(1X,A16,' : SWEEP = ',I5,' NORME SECOND MEMBRE = ',E14.6,  &
             ', RELAXP = ',E14.6)
 1500 format ( &
 1X,A16,' : Current reconstruction sweep = ',I5                     ,/,&
'           sweep residual = ',E12.5,', norm = ',E12.5              ,/,&
'           number of sweeps for solver = ',I5)
 1600 format( &
'@'                                                                 ,/,&
'@ @@ ATTENTION : ', A16,' ETAPE DE PRESSION'                       ,/,&
'@    ========='                                                    ,/,&
'@  Nombre d''iterations maximal ',I10   ,' atteint'                ,/,&
'@' )

#else

 1200 format ( &
 1X,A16,' Sweep: ',I5,' Dynamic relaxation: alpha = ',E12.5,' beta = ',E12.5,/,&
'    < dI^k  ; R^k > = ',E12.5,' ||dI^k  ||^2 = ',E12.5                     ,/,&
'    < dI^k-1; R^k > = ',E12.5,' ||dI^k-1||^2 = ',E12.5                     ,/,&
'   < dI^k-1; dI^k > = ',E12.5)
 1300 format(1X,A16,' : NORMED RESIDUALS = ', E14.6)
 1400 format(1X,A16,' : SWEEP = ',I5,' RIGHT HAND SIDE NORM = ',E14.6, &
             ', RELAXP = ',E14.6)
 1500 format ( &
 1X,A16,' : Current reconstruction sweep = ',I5                     ,/,&
'           sweep residual = ',E12.5,', norm = ',E12.5              ,/,&
'           number of sweeps for solver = ',I5)
 1600 format( &
'@'                                                                 ,/,&
'@ @@ WARNING: ', A16,' PRESSURE STEP'                              ,/,&
'@    ========'                                                     ,/,&
'@  Maximum number of iterations ',I10   ,' reached'                ,/,&
'@'                                                              )

#endif

!----
! End
!----

return

end subroutine resopv
