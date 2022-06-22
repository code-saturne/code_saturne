!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file navstv.f90
!>
!> \brief Solving of NS equations for incompressible or slightly compressible
!> flows for one time step. Both convection-diffusion and continuity steps are
!> performed.  The velocity components are solved together in once.
!>
!> Please refer to the
!> <a href="../../theory.pdf#navstv"><b>navstv</b></a> section
!> of the theory guide for more informations.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     iterns        index of the iteration on Navier-Stokes
!> \param[in]     icvrge        indicator of convergence
!> \param[in]     itrale        number of the current ALE iteration
!> \param[in]     isostd        indicator of standar outlet
!>                               +index of the reference face
!> \param[in]     dt            time step (per cell)
!> \param[in]     frcxt         external force generating the hydrostatic
!>                              pressure
!> \param[in]     trava         work array for pressure velocity coupling
!_______________________________________________________________________________


subroutine navstv &
 ( nvar   , nscal  , iterns , icvrge , itrale ,                   &
   isostd ,                                                       &
   dt     ,                                                       &
   frcxt  ,                                                       &
   trava  )

!===============================================================================

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
use paramx, only: isymet
use period
use ppppar
use ppthch
use ppincl
use cplsat
use mesh
use lagran, only: iilagr, ntersl
use rotation
use turbomachinery
use ptrglo
use field
use field_operator
use cavitation
use vof
use cs_c_bindings
use atincl, only: iatmst, iautom, imeteo
use cs_nz_condensation, only: nfbpcd, spcond, ifbpcd

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal  , iterns , icvrge , itrale

integer          isostd(nfabor+1)

double precision, pointer, dimension(:)   :: dt
double precision, pointer, dimension(:,:) :: frcxt
double precision, pointer, dimension(:,:) :: trava

! Local variables

integer          iccocg, inc, iel, iel1, iel2, ifac, imax, imaxt, imin, imint
integer          ii    , inod, itypfl, f_id, f_iddp
integer          isou, ivar, iitsm
integer          init, iautof
integer          iflmas, iflmab
integer          ivolfl_id, bvolfl_id
integer          iflmb0
integer          imrgrp, nswrgp, imligp, iwarnp
integer          nbrval, iappel
integer          ndircp, icpt
integer          numcpl
integer          f_dim , iflwgr
double precision rnorm , rnormt, rnorma, rnormi, vitnor
double precision dtsrom, unsrom, rhom, rovolsdt
double precision epsrgp, climgp, xyzmax(3), xyzmin(3)
double precision thetap, xdu, xdv, xdw
double precision rhofac, dtfac
double precision xnrdis, xnrtmp
double precision t1, t2, t3, t4
double precision visclc, visctc
double precision distbf, srfbnf, hint
double precision rnx, rny, rnz
double precision vr(3), vr1(3), vr2(3), vrn
double precision disp_fac(3)
double precision vol_fl_drhovol1, vol_fl_drhovol2

double precision, allocatable, dimension(:,:,:), target :: viscf
double precision, allocatable, dimension(:), target :: viscb
double precision, allocatable, dimension(:,:,:), target :: wvisfi
double precision, allocatable, dimension(:,:), target :: uvwk
double precision, dimension(:,:), pointer :: velk
double precision, allocatable, dimension(:), target :: wvisbi
double precision, allocatable, dimension(:), target :: cpro_rho_tc, bpro_rho_tc
double precision, allocatable, dimension(:) :: esflum, esflub
double precision, allocatable, dimension(:) :: intflx, bouflx
double precision, allocatable, dimension(:) :: secvif, secvib

double precision, allocatable, dimension(:,:), target :: gradp
double precision, dimension(:,:), pointer :: cpro_gradp
double precision, dimension(:), pointer :: coefa_dp, coefb_dp
double precision, dimension(:,:), pointer :: da_uu
double precision, dimension(:,:), pointer :: vel, vela
double precision, dimension(:,:,:), pointer :: viscfi
double precision, dimension(:), pointer :: viscbi
double precision, dimension(:,:), pointer :: dttens
double precision, dimension(:,:), pointer :: dfrcxt
double precision, dimension(:,:), pointer :: coefau, cofafu, claale
double precision, dimension(:,:,:), pointer :: coefbu, cofbfu, clbale
double precision, dimension(:), pointer :: coefa_p
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: ivolfl, bvolfl
double precision, dimension(:), pointer :: brom, broma, crom, croma, viscl, visct
double precision, dimension(:,:), pointer :: trav
double precision, dimension(:,:), pointer :: mshvel
double precision, dimension(:,:), pointer :: disale
double precision, dimension(:,:), pointer :: xyzno0
double precision, dimension(:), pointer :: porosi
double precision, dimension(:), pointer :: cvar_pr
double precision, dimension(:), pointer :: cpro_prtot, c_estim
double precision, dimension(:), pointer :: cvar_voidf, cvara_voidf
double precision, dimension(:), pointer :: cpro_rho_mass
double precision, dimension(:), pointer :: bpro_rho_mass
double precision, dimension(:), pointer :: brom_eos, crom_eos
double precision, dimension(:), pointer :: cpro_wgrec_s
double precision, dimension(:,:), pointer :: cpro_wgrec_v

type(var_cal_opt) :: vcopt_p, vcopt_u, vcopt

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_pressure_correction &
   ( iterns , nfbpcd , ncmast ,                                   &
     ifbpcd , ltmast , isostd ,                                   &
     vel    , da_uu  ,                                            &
     coefav , coefbv , coefa_dp        , coefb_dp ,               &
     spcond , svcond ,                                            &
     frcxt  , dfrcxt ,                                            &
     viscf  , viscb  )                                            &
    bind(C, name='cs_pressure_correction')

    use dimens, only: ndimfb
    use mesh

    implicit none

    ! Arguments

    integer, value :: iterns
    integer, value :: nfbpcd, ncmast

    integer          ifbpcd(nfbpcd)
    integer          ltmast(ncelet)
    integer          isostd(nfabor+1)

    double precision spcond(nfbpcd,*)
    double precision svcond(ncelet,*)
    double precision frcxt(3,ncelet), dfrcxt(3,ncelet)
    double precision viscf(nfac), viscb(ndimfb)
    double precision coefav(3,ndimfb)
    double precision coefbv(3,3,ndimfb)
    double precision vel(3,ncelet)
    double precision da_uu(6,ncelet)
    double precision coefa_dp(ndimfb)
    double precision coefb_dp(ndimfb)

  end subroutine cs_pressure_correction

  !=============================================================================

end interface

!===============================================================================

!===============================================================================
! 0. Initialization
!===============================================================================

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_p)

! Allocate temporary arrays for the velocity-pressure resolution
if (iand(vcopt_u%idften, ISOTROPIC_DIFFUSION).ne.0) then
  allocate(viscf(1, 1, nfac), viscb(ndimfb))
else if (iand(vcopt_u%idften, ANISOTROPIC_LEFT_DIFFUSION).ne.0) then
  allocate(viscf(3, 3, nfac), viscb(ndimfb))
endif

allocate(trav(3,ncelet))

! Allocate other arrays, depending on user options

call field_get_id("pressure_increment",f_iddp)

call field_get_coefa_s(f_iddp, coefa_dp)
call field_get_coefb_s(f_iddp, coefb_dp)

allocate(dfrcxt(3,ncelet))
if (iand(vcopt_u%idften, ISOTROPIC_DIFFUSION).ne.0) then
  if (itytur.eq.3.and.irijnu.eq.1) then
    allocate(wvisfi(1,1,nfac), wvisbi(ndimfb))
    viscfi => wvisfi(:,:,1:nfac)
    viscbi => wvisbi(1:ndimfb)
  else
    viscfi => viscf(:,:,1:nfac)
    viscbi => viscb(1:ndimfb)
  endif
else if (iand(vcopt_u%idften, ANISOTROPIC_LEFT_DIFFUSION).ne.0) then
  if (itytur.eq.3.and.irijnu.eq.1) then
    allocate(wvisfi(3,3,nfac), wvisbi(ndimfb))
    viscfi => wvisfi(1:3,1:3,1:nfac)
    viscbi => wvisbi(1:ndimfb)
  else
    viscfi => viscf(1:3,1:3,1:nfac)
    viscbi => viscb(1:ndimfb)
  endif
endif

if (ivisse.eq.1) then
  allocate(secvif(nfac),secvib(ndimfb))
endif

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)
call field_get_val_prev_v(ivarfl(iu), vela)

! Map some specific field arrays
if (idtten.ge.0) then
  call field_get_val_v(idtten, dttens)
else
  dttens => rvoid2
endif

call field_get_val_s(ivarfl(ipr), cvar_pr)

if (ivofmt.gt.0) then
  call field_get_val_s(ivarfl(ivolf2), cvar_voidf)
  call field_get_val_prev_s(ivarfl(ivolf2), cvara_voidf)
endif

! Initialize variables to avoid compiler warnings

ivar = 0
iflmas = 0
imax = 0

! pointer to velosity at sub iteration k for velocity-pressure inner iterations
if (nterup.gt.1) then

  allocate(uvwk(3, ncelet))
  !$omp parallel do private(isou)
  do iel = 1, ncel
    do isou = 1, 3
      uvwk(isou,iel) = vel(isou,iel)
    enddo
  enddo

  ! Compute the L2 velocity norm (it is zero at the first time step, so
  ! we recompute it)
  if (iterns.eq.1.or.xnrmu0.eq.0d0) then
    xnrtmp = 0.d0
    !$omp parallel do reduction(+:xnrtmp)
    do iel = 1, ncel
      xnrtmp = xnrtmp +(vel(1,iel)**2        &
                      + vel(2,iel)**2        &
                      + vel(3,iel)**2)       &
                      * cell_f_vol(iel)
    enddo
    xnrmu0 = xnrtmp
    if (irangp.ge.0) then
      call parsom (xnrmu0)
    endif
    ! En cas de couplage entre deux instances de code_saturne, on calcule
    ! la norme totale de la vitesse
    ! Necessaire pour que l'une des instances ne stoppe pas plus tot que les autres
    ! (il faudrait quand meme verifier les options numeriques, ...)
    do numcpl = 1, nbrcpl
      call tbrcpl ( numcpl, 1, 1, xnrmu0, xnrdis )
      xnrmu0 = xnrmu0 + xnrdis
    enddo
    xnrmu0 = sqrt(xnrmu0)
  endif

  ! On assure la periodicite ou le parallelisme de uvwk et la pression
  if (irangp.ge.0.or.iperio.eq.1) then
    call synvin(uvwk)
    call synsca(cvar_pr)
  endif

  velk => uvwk

else
  velk => vela
endif

! --- Physical quantities
call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

! Initialize timers
t1 = 0.d0
t2 = 0.d0
t3 = 0.d0
t4 = 0.d0

! Id of the mass flux
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)

! Pointers to the mass fluxes
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

! Pointers to properties
call field_get_val_s(icrom, crom_eos)
call field_get_val_s(ibrom, brom_eos)

if (irovar.eq.1.and.(idilat.gt.1.or.ivofmt.gt.0.or.ippmod(icompf).eq.3)) then
  ! If iterns = 1: this is density at time n
  call field_get_id("density_mass", f_id)
  call field_get_val_s(f_id, cpro_rho_mass)
  call field_get_id("boundary_density_mass", f_id)
  call field_get_val_s(f_id, bpro_rho_mass)

  ! Time interpolated density
  if (vcopt_u%thetav .lt. 1.d0 .and. itpcol .eq. 0) then
    call field_get_val_prev_s(icrom, croma)
    call field_get_val_prev_s(ibrom, broma)
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

! Pointers to BC coefficients
call field_get_coefa_v(ivarfl(iu), coefau)
call field_get_coefb_v(ivarfl(iu), coefbu)
call field_get_coefaf_v(ivarfl(iu), cofafu)
call field_get_coefbf_v(ivarfl(iu), cofbfu)

!===============================================================================
! 1. Prediction of the mass flux in case of Low Mach compressible algorithm
!===============================================================================

if ((idilat.eq.2.or.idilat.eq.3).and. &
    (ntcabs.gt.1.or.isuite.gt.0).and.ipredfl.ne.0) then

  call predfl(nvar, ncetsm, icetsm, dt, smacel)

endif

!===============================================================================
! 3. Pressure resolution and computation of mass flux for compressible flow
!===============================================================================

! Note, for the compressible algorithm written in pressure increment,
! this step is merged with the pressure correction step of the incompressible
! algorithm
if (ippmod(icompf).ge.0.and.ippmod(icompf).ne.3) then

  if(vcopt_p%iwarni.ge.1) then
    write(nfecra,1080)
  endif

  call cfmspr &
  ( nvar   , nscal  , iterns ,                                     &
    ncepdc , ncetsm , icepdc , icetsm , itypsm ,                   &
    dt     , vela   ,                                              &
    ckupdc , smacel )

endif

!===============================================================================
! 4. VoF: compute liquid-vapor mass transfer term (cavitating flows)
!===============================================================================

if (iand(ivofmt,VOF_MERKLE_MASS_TRANSFER).ne.0) then

  call field_get_val_s(iprtot, cpro_prtot)
  call cavitation_compute_source_term(cpro_prtot, cvara_voidf)

endif

!===============================================================================
! 5. Velocity prediction step
!===============================================================================

if (vcopt_u%iwarni.ge.1) then
  write(nfecra,1000)
endif

iappel = 1

allocate(da_uu(6,ncelet))

call predvv &
( iappel ,                                                       &
  nvar   , nscal  , iterns ,                                     &
  ncepdc , ncetsm ,                                              &
  icepdc , icetsm , itypsm ,                                     &
  dt     , vel    , vela   , velk   , da_uu  ,                   &
  tslagr , coefau , coefbu , cofafu , cofbfu ,                   &
  ckupdc , smacel , frcxt  ,                                     &
  trava  ,                   dfrcxt , dttens ,  trav  ,          &
  viscf  , viscb  , viscfi , viscbi , secvif , secvib )


! Bad cells regularisation
call cs_bad_cells_regularisation_vector(vel, 1)

! --- Sortie si pas de pression continuite
!       on met a jour les flux de masse, et on sort

if (iprco.le.0) then

  itypfl = 1
  init   = 1
  inc    = 1
  iflmb0 = 1
  if (iale.ge.1) iflmb0 = 0
  imrgrp = vcopt_u%imrgra
  nswrgp = vcopt_u%nswrgr
  imligp = vcopt_u%imligr
  iwarnp = vcopt_u%iwarni
  epsrgp = vcopt_u%epsrgr
  climgp = vcopt_u%climgr

  call inimav                                                     &
 ( ivarfl(iu)      , itypfl ,                                     &
   iflmb0 , init   , inc    , imrgrp , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom   , brom   ,                                              &
   vel    ,                                                       &
   coefau , coefbu ,                                              &
   imasfl , bmasfl )

  ! In the ALE framework, we add the mesh velocity
  if (iale.ge.1) then

    call field_get_val_v(ivarfl(iuma), mshvel)

    call field_get_val_v(fdiale, disale)
    call field_get_val_v_by_name("vtx_coord0", xyzno0)

    if (iflxmw.gt.0) then
      ! One temporary array needed for internal faces, in case some internal vertices
      !  are moved directly by the user
      allocate(intflx(nfac), bouflx(ndimfb))

      itypfl = 1
      init   = 1
      inc    = 1
      iflmb0 = 1
      call field_get_key_struct_var_cal_opt(ivarfl(iuma), vcopt)
      imrgrp = vcopt%imrgra
      nswrgp = vcopt%nswrgr
      imligp = vcopt%imligr
      iwarnp = vcopt%iwarni
      epsrgp = vcopt%epsrgr
      climgp = vcopt%climgr

      call field_get_coefa_v(ivarfl(iuma), claale)
      call field_get_coefb_v(ivarfl(iuma), clbale)

      call inimav &
    ( ivarfl(iuma)    , itypfl ,                                     &
      iflmb0 , init   , inc    , imrgrp , nswrgp , imligp ,          &
      iwarnp ,                                                       &
      epsrgp , climgp ,                                              &
      crom, brom,                                                    &
      mshvel ,                                                       &
      claale , clbale ,                                              &
      intflx , bouflx )
    endif

    ! Here we need of the opposite of the mesh velocity.
    do ifac = 1, nfabor
      ! Compute the mass flux using the nodes displacement
      if (iflxmw.eq.0) then
        disp_fac(1) = 0.d0
        disp_fac(2) = 0.d0
        disp_fac(3) = 0.d0
        icpt  = 0
        do ii = ipnfbr(ifac), ipnfbr(ifac+1)-1
          inod = nodfbr(ii)
          icpt = icpt + 1
          disp_fac(1) = disp_fac(1) + disale(1,inod) - (xyznod(1,inod)-xyzno0(1,inod))
          disp_fac(2) = disp_fac(2) + disale(2,inod) - (xyznod(2,inod)-xyzno0(2,inod))
          disp_fac(3) = disp_fac(3) + disale(3,inod) - (xyznod(3,inod)-xyzno0(3,inod))
        enddo
        iel = ifabor(ifac)
        bmasfl(ifac) = bmasfl(ifac) - brom(ifac) * (              &
              disp_fac(1) * surfbo(1,ifac)                        &
             +disp_fac(2) * surfbo(2,ifac)                        &
             +disp_fac(3) * surfbo(3,ifac) )/dt(iel)/icpt
      else
        bmasfl(ifac) = bmasfl(ifac) - bouflx(ifac)
      endif
    enddo

    do ifac = 1, nfac
      ! Compute the mass flux using the nodes displacement
      if (iflxmw.eq.0) then
        disp_fac(1) = 0.d0
        disp_fac(2) = 0.d0
        disp_fac(3) = 0.d0
        icpt  = 0
        do ii = ipnfac(ifac),ipnfac(ifac+1)-1
          inod = nodfac(ii)
          icpt = icpt + 1
          disp_fac(1) = disp_fac(1) + disale(1,inod) - (xyznod(1,inod)-xyzno0(1,inod))
          disp_fac(2) = disp_fac(2) + disale(2,inod) - (xyznod(2,inod)-xyzno0(2,inod))
          disp_fac(3) = disp_fac(3) + disale(3,inod) - (xyznod(3,inod)-xyzno0(3,inod))
        enddo
        ! For inner vertices, the mass flux due to the mesh displacement is
        !  recomputed from the nodes displacement
        iel1 = ifacel(1,ifac)
        iel2 = ifacel(2,ifac)
        dtfac = 0.5d0*(dt(iel1) + dt(iel2))
        rhofac = 0.5d0*(crom(iel1) + crom(iel2))
        imasfl(ifac) = imasfl(ifac) - rhofac*(                    &
              disp_fac(1) * surfac(1,ifac)                        &
             +disp_fac(2) * surfac(2,ifac)                        &
             +disp_fac(3) * surfac(3,ifac) )/dtfac/icpt
      else
        imasfl(ifac) = imasfl(ifac) - intflx(ifac)
      endif
    enddo

    if (iflxmw.gt.0) then
      ! Free memory
      deallocate(intflx, bouflx)
    endif
  endif

  ! Ajout de la vitesse du solide dans le flux convectif,
  ! si le maillage est mobile (solide rigide)
  ! En turbomachine, on connait exactement la vitesse de maillage a ajouter
  if (iturbo.ne.0) then
    !$omp parallel do private(iel1, iel2, rhofac, vr1, vr2)
    do ifac = 1, nfac
      iel1 = ifacel(1,ifac)
      iel2 = ifacel(2,ifac)
      if (irotce(iel1).ne.0 .or. irotce(iel2).ne.0) then

        rhofac = 0.5d0*(crom(iel1) + crom(iel2))
        call rotation_velocity(irotce(iel1), cdgfac(:,ifac), vr1)
        call rotation_velocity(irotce(iel2), cdgfac(:,ifac), vr2)

        imasfl(ifac) = imasfl(ifac) - 0.5d0 * rhofac*(      &
                          surfac(1,ifac)*(vr1(1) + vr2(1))  &
                        + surfac(2,ifac)*(vr1(2) + vr2(2))  &
                        + surfac(3,ifac)*(vr1(3) + vr2(3)) )
      endif
    enddo
    !$omp parallel do private(iel, rhofac, vr) if(nfabor > thr_n_min)
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      if (irotce(iel).ne.0) then

        rhofac = brom(ifac)
        call rotation_velocity(irotce(iel), cdgfbo(:,ifac), vr)

        bmasfl(ifac) = bmasfl(ifac) - rhofac*( surfbo(1,ifac)*vr(1) &
                                             + surfbo(2,ifac)*vr(2) &
                                             + surfbo(3,ifac)*vr(3) )
      endif
    enddo
  endif

  ! Free memory
  deallocate(viscf, viscb)
  deallocate(trav)
  deallocate(da_uu)
  deallocate(dfrcxt)
  if (allocated(wvisfi)) deallocate(wvisfi, wvisbi)
  if (allocated(uvwk)) deallocate(uvwk)
  if (allocated(secvif)) deallocate(secvif, secvib)
  if (allocated(cpro_rho_tc)) deallocate(cpro_rho_tc)
  if (allocated(bpro_rho_tc)) deallocate(bpro_rho_tc)
  return

endif

!===============================================================================
! 6. Update mesh for unsteady turbomachinery computations
!===============================================================================

if (iturbo.eq.2 .and. iterns.eq.1) then

  ! Update mesh

  call turbomachinery_update_mesh (ttcmob, rs_ell(1))

  call dmtmps(t1)

  if (ityint.eq.0) then

    do ifac = 1, nfabor
      ! --- To cancel the mass flux for symmetry BC
      if (itypfb(ifac).eq.isymet) then
        isympa(ifac) = 0
      else
        isympa(ifac) = 1
      endif
    enddo

    ! Scratch and resize temporary internal faces arrays

    deallocate(viscf)
    if (iand(vcopt_u%idften, ISOTROPIC_DIFFUSION).ne.0) then
      allocate(viscf(1, 1, nfac))
    else if (iand(vcopt_u%idften, ANISOTROPIC_LEFT_DIFFUSION).ne.0) then
      allocate(viscf(3, 3, nfac))
    endif

    if (allocated(wvisfi)) then
      deallocate(viscfi)

      if (vcopt_u%idften.eq.1) then
        if (itytur.eq.3.and.irijnu.eq.1) then
          allocate(wvisfi(1,1,nfac))
          viscfi => wvisfi(:,:,1:nfac)
        else
          viscfi => viscf(:,:,1:nfac)
        endif
      else if(vcopt_u%idften.eq.6) then
        if (itytur.eq.3.and.irijnu.eq.1) then
          allocate(wvisfi(3,3,nfac))
          viscfi => wvisfi(1:3,1:3,1:nfac)
        else
          viscfi => viscf(1:3,1:3,1:nfac)
        endif
      endif

    endif

    if (allocated(secvif)) then
      deallocate(secvif)
      allocate(secvif(nfac))
    endif

    ! Scratch, resize and initialize main internal faces properties array

    call turbomachinery_reinit_i_face_fields

    ! Update local pointers on "internal faces" fields

    call field_get_val_s(iflmas, imasfl)

    if (irangp.ge.0 .or. iperio.eq.1) then

      ! Resize auxiliary arrays (pointe module)

      call resize_aux_arrays

      ! Update turbomachinery module

      call turbomachinery_update

      ! Update field mappings ("owner" fields handled by turbomachinery_update)

      call fldtri
      call field_get_val_s_by_name('dt', dt)

      ! Resize other arrays related to the velocity-pressure resolution

      call resize_sym_tens_real_array(da_uu)
      call resize_vec_real_array(trav)
      call resize_vec_real_array(dfrcxt)

      ! Resize other arrays, depending on user options

      if (iilagr.gt.0 .and. ntersl.gt.0) &
        call resize_n_sca_real_arrays(ntersl, tslagr)

      if (iphydr.eq.1) then
        call field_get_val_v_by_name('volume_forces', frcxt)
      endif

      ! Update local pointers on "cells" fields

      call field_get_val_s(icrom, crom)
      call field_get_val_s(icrom, crom_eos)

      if (irovar.eq.1.and.(idilat.gt.1.or.ivofmt.gt.0.or.ippmod(icompf).eq.3)) then
        ! If iterns = 1: this is density at time n
        call field_get_id("density_mass", f_id)
        call field_get_val_s(f_id, cpro_rho_mass)

        ! Time interpolated density
        if (vcopt_u%thetav .lt. 1.d0 .and. itpcol .eq. 0) then

          call field_get_val_prev_s(icrom, croma)

          if (allocated(cpro_rho_tc)) deallocate(cpro_rho_tc)
          allocate(cpro_rho_tc(ncelet))

          do iel = 1, ncelet
            cpro_rho_tc(iel) =  vcopt_u%thetav * cpro_rho_mass(iel) &
                               + (1.d0 - vcopt_u%thetav) * croma(iel)
          enddo

          crom => cpro_rho_tc

        else
          crom => cpro_rho_mass
        endif
      endif

      call field_get_val_s(iviscl, viscl)
      call field_get_val_s(ivisct, visct)

      call field_get_val_v(ivarfl(iu), vel)
      call field_get_val_prev_v(ivarfl(iu), vela)

      call field_get_val_s(ivarfl(ipr), cvar_pr)

      if (idtten.ge.0) call field_get_val_v(idtten, dttens)

      if (ivofmt.gt.0) then
        call field_get_val_s(ivarfl(ivolf2), cvar_voidf)
        call field_get_val_prev_s(ivarfl(ivolf2), cvara_voidf)
      endif

      if (nterup.gt.1) then
        call resize_vec_real_array(velk)
        call resize_vec_real_array(trava)
      else
        velk => vela
      endif

    endif

  endif

  ! Update the Dirichlet wall boundary conditions for velocity (based on the
  ! solid body rotation on the new mesh).
  ! Note that the velocity BC update is made only if the user has not specified
  ! any specific Dirichlet condition for velocity.

  do ifac = 1, nfabor

    iel = ifabor(ifac)

    if (coftur(ifac).lt.rinfin*0.5d0) then

      ! --- Physical Propreties
      visclc = viscl(iel)
      visctc = visct(iel)

      ! --- Geometrical quantities
      distbf = distb(ifac)
      srfbnf = surfbn(ifac)

      ! Unit normal
      rnx = surfbo(1,ifac)/srfbnf
      rny = surfbo(2,ifac)/srfbnf
      rnz = surfbo(3,ifac)/srfbnf

      if (itytur.eq.3) then
        hint =   visclc         /distbf
      else
        hint = ( visclc+visctc )/distbf
      endif

      call rotation_velocity(irotce(iel), cdgfbo(:,ifac), vr)

      ! Gradient boundary conditions (Dirichlet)
      !-----------------------------
      vrn = vr(1)*rnx + vr(2)*rny + vr(3)*rnz

      coefau(1,ifac) = (1.d0-coftur(ifac))*(vr(1) - vrn*rnx) + vrn*rnx
      coefau(2,ifac) = (1.d0-coftur(ifac))*(vr(2) - vrn*rny) + vrn*rny
      coefau(3,ifac) = (1.d0-coftur(ifac))*(vr(3) - vrn*rnz) + vrn*rnz

      ! Flux boundary conditions (Dirichlet)
      !-------------------------

      cofafu(1,ifac) = -hfltur(ifac)*(vr(1) - vrn*rnx) - hint*vrn*rnx
      cofafu(2,ifac) = -hfltur(ifac)*(vr(2) - vrn*rny) - hint*vrn*rny
      cofafu(3,ifac) = -hfltur(ifac)*(vr(3) - vrn*rnz) - hint*vrn*rnz

    endif

  enddo

  call dmtmps(t2)

  rs_ell(2) = t2-t1

endif

!===============================================================================
! 7. Pressure correction step
!===============================================================================

if (vcopt_u%iwarni.ge.1) then
  write(nfecra,1200)
endif

if (ippmod(icompf).lt.0.or.ippmod(icompf).eq.3) then

  call cs_pressure_correction                &
    (iterns, nfbpcd, ncmast,                 &
     ifbpcd, ltmast, isostd,                 &
     vel, da_uu,                             &
     coefau, coefbu, coefa_dp, coefb_dp,     &
     spcond, svcond,                         &
     frcxt, dfrcxt,                          &
     viscf, viscb)

endif

! Bad cells regularisation
call cs_bad_cells_regularisation_scalar(cvar_pr)

!===============================================================================
! 8. Mesh velocity solving (ALE)
!===============================================================================

if (iale.ge.1) then

  if (itrale.gt.nalinf) then
    call cs_ale_solve_mesh_velocity(iterns, impale, ialtyb)
  endif

endif

!===============================================================================
! 9. Update of the fluid velocity field
!===============================================================================

if (ippmod(icompf).lt.0.or.ippmod(icompf).eq.3) then

  ! irevmc = 0: Update the velocity with the pressure gradient.

  if (irevmc.eq.0) then

    ! The predicted velocity is corrected by the cell gradient of the
    ! pressure increment.

    iccocg = 1
    inc = 0

    call grdpor(inc)

    if (iphydr.eq.1.or.iifren.eq.1) inc = 1

    ! Pressure increment gradient
    call field_get_id_try("pressure_increment_gradient", f_id)
    if (f_id.ge.0) then
      call field_get_val_v(f_id, cpro_gradp)
    else
      !Allocation
      allocate(gradp(3,ncelet))
      cpro_gradp => gradp
    endif

    if (ivofmt.ne.0) then
      call field_get_key_int(ivarfl(ipr), kwgrec, iflwgr)
      call field_get_dim(iflwgr, f_dim)
      if (f_dim.eq.1) then
        call field_get_val_s(iflwgr, cpro_wgrec_s)
        do iel = 1, ncel
          cpro_wgrec_s(iel) = dt(iel) / crom(iel)
        enddo
        call synsca(cpro_wgrec_s)
      else if (f_dim.eq.6) then
        call field_get_val_v(iflwgr, cpro_wgrec_v)
        do iel = 1, ncel
          do ii = 1, 6
            cpro_wgrec_v(ii,iel) = dttens(ii,iel) / crom(iel)
          enddo
        enddo
        call syntis(cpro_wgrec_v)
      endif
    endif

    if (iprcdo.eq.0) then
      call field_gradient_potential(f_iddp, 0, 0, inc,                   &
                                    iccocg, iphydr,                      &
                                    dfrcxt, cpro_gradp)
    endif
    ! Update the velocity field
    !--------------------------
    thetap = vcopt_p%thetav

    ! Specific handling of hydrostatic pressure
    !------------------------------------------
    if (iphydr.eq.1) then

      ! Scalar diffusion for the pressure
      if (iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then
        !$omp parallel do private(dtsrom, isou)
        do iel = 1, ncel
          dtsrom = thetap*dt(iel)/crom(iel)
          do isou = 1, 3
            vel(isou,iel) = vel(isou,iel)                            &
                 + dtsrom*(dfrcxt(isou, iel)-cpro_gradp(isou,iel))
          enddo
        enddo

      ! Tensorial diffusion for the pressure
      elseif (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then
        !$omp parallel do private(unsrom)
        do iel = 1, ncel
          unsrom = thetap/crom(iel)

          vel(1, iel) = vel(1, iel)                              &
               + unsrom*(                                        &
                 dttens(1,iel)*(dfrcxt(1, iel)-cpro_gradp(1,iel))     &
               + dttens(4,iel)*(dfrcxt(2, iel)-cpro_gradp(2,iel))     &
               + dttens(6,iel)*(dfrcxt(3, iel)-cpro_gradp(3,iel))     &
               )
          vel(2, iel) = vel(2, iel)                              &
               + unsrom*(                                        &
                 dttens(4,iel)*(dfrcxt(1, iel)-cpro_gradp(1,iel))     &
               + dttens(2,iel)*(dfrcxt(2, iel)-cpro_gradp(2,iel))     &
               + dttens(5,iel)*(dfrcxt(3, iel)-cpro_gradp(3,iel))     &
               )
          vel(3, iel) = vel(3, iel)                              &
               + unsrom*(                                        &
                 dttens(6,iel)*(dfrcxt(1 ,iel)-cpro_gradp(1,iel))     &
               + dttens(5,iel)*(dfrcxt(2 ,iel)-cpro_gradp(2,iel))     &
               + dttens(3,iel)*(dfrcxt(3 ,iel)-cpro_gradp(3,iel))     &
               )
        enddo
      endif

      ! Update of the Dirichlet boundary conditions on the
      ! pressure for the outlet
      call field_get_coefa_s(ivarfl(ipr), coefa_p)
      !$omp parallel do if(nfabor > thr_n_min) private(iautof)
      do ifac = 1, nfabor
        iautof = 0
        ! automatic inlet/outlet face for atmospheric flow
        if (imeteo.gt.0) then
          iautof = iautom(ifac)
        endif

        if (isostd(ifac).eq.1.or.iatmst.ge.1.and.iautof.ge.1) then
          coefa_p(ifac) = coefa_p(ifac) + coefa_dp(ifac)
        endif
      enddo


      ! Standard handling of hydrostatic pressure
      !------------------------------------------
    else

      ! Scalar diffusion for the pressure
      if (iand(vcopt_p%idften, ISOTROPIC_DIFFUSION).ne.0) then

        !$omp parallel do private(dtsrom, isou)
        do iel = 1, ncel
          dtsrom = thetap*dt(iel)/crom(iel)
          do isou = 1, 3
            vel(isou,iel) = vel(isou,iel) - dtsrom*cpro_gradp(isou,iel)
          enddo
        enddo

      ! Tensorial diffusion for the pressure
      else if (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then

        !$omp parallel do private(unsrom)
        do iel = 1, ncel
          unsrom = thetap/crom(iel)

          vel(1, iel) = vel(1, iel)                              &
                      - unsrom*(                                 &
                                 dttens(1,iel)*(cpro_gradp(1,iel))    &
                               + dttens(4,iel)*(cpro_gradp(2,iel))    &
                               + dttens(6,iel)*(cpro_gradp(3,iel))    &
                               )
          vel(2, iel) = vel(2, iel)                              &
                      - unsrom*(                                 &
                                 dttens(4,iel)*(cpro_gradp(1,iel))    &
                               + dttens(2,iel)*(cpro_gradp(2,iel))    &
                               + dttens(5,iel)*(cpro_gradp(3,iel))    &
                               )
          vel(3, iel) = vel(3, iel)                              &
                      - unsrom*(                                 &
                                 dttens(6,iel)*(cpro_gradp(1,iel))    &
                               + dttens(5,iel)*(cpro_gradp(2,iel))    &
                               + dttens(3,iel)*(cpro_gradp(3,iel))    &
                               )
        enddo

      endif
    endif

    !Free memory
    if (allocated(gradp)) deallocate(gradp)

  ! RT0 update from the mass fluxes
  else

    ! Initialization to 0
    do iel = 1, ncelet
      vel(1, iel) = 0.d0
      vel(2, iel) = 0.d0
      vel(3, iel) = 0.d0
    enddo

    ! vel = 1 / (rho Vol) SUM mass_flux (X_f - X_i)
    if (ivofmt.eq.0) then

      do ifac = 1, nfac

        iel1 = ifacel(1,ifac)
        iel2 = ifacel(2,ifac)

        vol_fl_drhovol1 = 0.d0
        ! If it is not a solid cell
        if (cell_is_active(iel1).eq.1) &
          vol_fl_drhovol1 = imasfl(ifac) / (crom(iel1) * cell_f_vol(iel1))


        vol_fl_drhovol2 = 0.d0
        ! If it is not a solid cell
        if (cell_is_active(iel2).eq.1) &
          vol_fl_drhovol2 = imasfl(ifac) / (crom(iel2) * cell_f_vol(iel2))

        vel(1, iel1) = vel(1, iel1) &
          + vol_fl_drhovol1 * (cdgfac(1, ifac) - xyzcen(1, iel1))
        vel(2, iel1) = vel(2, iel1) &
          + vol_fl_drhovol1 * (cdgfac(2, ifac) - xyzcen(2, iel1))
        vel(3, iel1) = vel(3, iel1) &
          + vol_fl_drhovol1 * (cdgfac(3, ifac) - xyzcen(3, iel1))

        vel(1, iel2) = vel(1, iel2) &
          - vol_fl_drhovol2 * (cdgfac(1, ifac) - xyzcen(1, iel2))
        vel(2, iel2) = vel(2, iel2) &
          - vol_fl_drhovol2 * (cdgfac(2, ifac) - xyzcen(2, iel2))
        vel(3, iel2) = vel(3, iel2) &
          - vol_fl_drhovol2 * (cdgfac(3, ifac) - xyzcen(3, iel2))

      enddo

      do ifac = 1, nfabor
        iel1 = ifabor(ifac)

        vol_fl_drhovol1 = 0.d0
        ! If it is not a solid cell
        if (cell_is_active(iel1).eq.1) &
          vol_fl_drhovol1 = bmasfl(ifac) / (crom(iel1) * cell_f_vol(iel1))

        vel(1, iel1) = vel(1, iel1) &
          + vol_fl_drhovol1 * (cdgfbo(1, ifac) - xyzcen(1, iel1))
        vel(2, iel1) = vel(2, iel1) &
          + vol_fl_drhovol1 * (cdgfbo(2, ifac) - xyzcen(2, iel1))
        vel(3, iel1) = vel(3, iel1) &
          + vol_fl_drhovol1 * (cdgfbo(3, ifac) - xyzcen(3, iel1))

      enddo

      ! VOF module
      ! vel = 1 / (Vol) SUM vol_flux (X_f - X_i)
    else

      ! Id of the volume flux
      call field_get_key_int(ivarfl(ivolf2), kimasf, ivolfl_id)
      call field_get_key_int(ivarfl(ivolf2), kbmasf, bvolfl_id)

      ! Pointers to the mass fluxes
      call field_get_val_s(ivolfl_id, ivolfl)
      call field_get_val_s(bvolfl_id, bvolfl)

      do ifac = 1, nfac

        iel1 = ifacel(1,ifac)
        iel2 = ifacel(2,ifac)

        vol_fl_drhovol1 = 0.d0
        ! If it is not a solid cell
        if (cell_is_active(iel1).eq.1) &
          vol_fl_drhovol1 = ivolfl(ifac) / cell_f_vol(iel1)

        vol_fl_drhovol2 = 0.d0
        ! If it is not a solid cell
        if (cell_is_active(iel2).eq.1) &
          vol_fl_drhovol2 = ivolfl(ifac) / cell_f_vol(iel2)

        vel(1, iel1) = vel(1, iel1) &
          + vol_fl_drhovol1 * (cdgfac(1, ifac) - xyzcen(1, iel1))
        vel(2, iel1) = vel(2, iel1) &
          + vol_fl_drhovol1 * (cdgfac(2, ifac) - xyzcen(2, iel1))
        vel(3, iel1) = vel(3, iel1) &
          + vol_fl_drhovol1 * (cdgfac(3, ifac) - xyzcen(3, iel1))

        vel(1, iel2) = vel(1, iel2) &
          - vol_fl_drhovol2 * (cdgfac(1, ifac) - xyzcen(1, iel2))
        vel(2, iel2) = vel(2, iel2) &
          - vol_fl_drhovol2 * (cdgfac(2, ifac) - xyzcen(2, iel2))
        vel(3, iel2) = vel(3, iel2) &
          - vol_fl_drhovol2 * (cdgfac(3, ifac) - xyzcen(3, iel2))

      enddo

      do ifac = 1, nfabor
        iel1 = ifabor(ifac)

        vol_fl_drhovol1 = 0.d0
        ! If it is not a solid cell
        if (cell_is_active(iel1).eq.1) &
          vol_fl_drhovol1 = bvolfl(ifac) / cell_f_vol(iel1)

        vel(1, iel1) = vel(1, iel1) &
          + vol_fl_drhovol1 * (cdgfbo(1, ifac) - xyzcen(1, iel1))
        vel(2, iel1) = vel(2, iel1) &
          + vol_fl_drhovol1 * (cdgfbo(2, ifac) - xyzcen(2, iel1))
        vel(3, iel1) = vel(3, iel1) &
          + vol_fl_drhovol1 * (cdgfbo(3, ifac) - xyzcen(3, iel1))

      enddo

    endif

    call synvin(vel)

  endif

  if (iphydr.eq.1) then

    ! Update external forces for the computation of the gradients
    !$omp parallel do
    do iel=1,ncel
      frcxt(1 ,iel) = frcxt(1 ,iel) * cell_is_active(iel) + dfrcxt(1 ,iel)
      frcxt(2 ,iel) = frcxt(2 ,iel) * cell_is_active(iel) + dfrcxt(2 ,iel)
      frcxt(3 ,iel) = frcxt(3 ,iel) * cell_is_active(iel) + dfrcxt(3 ,iel)
    enddo

    if (irangp.ge.0.or.iperio.eq.1) then
      call synvin(frcxt)
    endif

  endif

endif

! Bad cells regularisation
call cs_bad_cells_regularisation_vector(vel, 1)

! Mass flux initialization for VOF algorithm
if (ivofmt.gt.0) then
  do ifac = 1, nfac
    imasfl(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    bmasfl(ifac) = 0.d0
  enddo
endif

! In the ALE framework, we add the mesh velocity
if (iale.ge.1) then

  call field_get_val_v(ivarfl(iuma), mshvel)

  call field_get_val_v(fdiale, disale)
  call field_get_val_v_by_name("vtx_coord0", xyzno0)

  if (iflxmw.gt.0) then
    ! One temporary array needed for internal faces, in case some internal vertices
    !  are moved directly by the user
    allocate(intflx(nfac), bouflx(ndimfb))

    itypfl = 1
    init   = 1
    inc    = 1
    iflmb0 = 1
    call field_get_key_struct_var_cal_opt(ivarfl(iuma), vcopt)
    imrgrp = vcopt%imrgra
    nswrgp = vcopt%nswrgr
    imligp = vcopt%imligr
    iwarnp = vcopt%iwarni
    epsrgp = vcopt%epsrgr
    climgp = vcopt%climgr

    call field_get_coefa_v(ivarfl(iuma), claale)
    call field_get_coefb_v(ivarfl(iuma), clbale)

    call inimav &
  ( ivarfl(iuma)    , itypfl ,                                     &
    iflmb0 , init   , inc    , imrgrp , nswrgp , imligp ,          &
    iwarnp ,                                                       &
    epsrgp , climgp ,                                              &
    crom, brom,                                                    &
    mshvel ,                                                       &
    claale , clbale ,                                              &
    intflx , bouflx )
  endif

  ! Here we need of the opposite of the mesh velocity.
  !$omp parallel do private(disp_fac, icpt, ii, inod, iel) if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    ! Compute the mass flux using the nodes displacement
    if (iflxmw.eq.0) then
      disp_fac(1) = 0.d0
      disp_fac(2) = 0.d0
      disp_fac(3) = 0.d0
      icpt  = 0
      do ii = ipnfbr(ifac),ipnfbr(ifac+1)-1
        inod = nodfbr(ii)
        icpt = icpt + 1
        disp_fac(1) = disp_fac(1) + disale(1,inod) - (xyznod(1,inod)-xyzno0(1,inod))
        disp_fac(2) = disp_fac(2) + disale(2,inod) - (xyznod(2,inod)-xyzno0(2,inod))
        disp_fac(3) = disp_fac(3) + disale(3,inod) - (xyznod(3,inod)-xyzno0(3,inod))
      enddo
      iel = ifabor(ifac)
      bmasfl(ifac) = bmasfl(ifac) - brom(ifac) * (           &
         disp_fac(1) * surfbo(1,ifac)                        &
        +disp_fac(2) * surfbo(2,ifac)                        &
        +disp_fac(3) * surfbo(3,ifac) )/dt(iel)/icpt
    else
      bmasfl(ifac) = bmasfl(ifac) - bouflx(ifac)
    endif
  enddo

  do ifac = 1, nfac
    ! Compute the mass flux using the nodes displacement
    if (iflxmw.eq.0) then
      disp_fac(1) = 0.d0
      disp_fac(2) = 0.d0
      disp_fac(3) = 0.d0
      icpt  = 0
      do ii = ipnfac(ifac),ipnfac(ifac+1)-1
        inod = nodfac(ii)
        icpt = icpt + 1
        disp_fac(1) = disp_fac(1) + disale(1,inod) - (xyznod(1,inod)-xyzno0(1,inod))
        disp_fac(2) = disp_fac(2) + disale(2,inod) - (xyznod(2,inod)-xyzno0(2,inod))
        disp_fac(3) = disp_fac(3) + disale(3,inod) - (xyznod(3,inod)-xyzno0(3,inod))
      enddo

      ! For inner vertices, the mass flux due to the mesh displacement is
      !  recomputed from the nodes displacement
      iel1 = ifacel(1,ifac)
      iel2 = ifacel(2,ifac)
      dtfac = 0.5d0*(dt(iel1) + dt(iel2))
      rhofac = 0.5d0*(crom(iel1) + crom(iel2))
      imasfl(ifac) = imasfl(ifac) - rhofac*(                    &
            disp_fac(1) * surfac(1,ifac)                        &
           +disp_fac(2) * surfac(2,ifac)                        &
           +disp_fac(3) * surfac(3,ifac) )/dtfac/icpt
    else
      imasfl(ifac) = imasfl(ifac) - intflx(ifac)
    endif
  enddo

  if (iflxmw.gt.0) then
    ! Free memory
    deallocate(intflx, bouflx)
  endif
endif

!FIXME for me we should do that before predvv
! Ajout de la vitesse du solide dans le flux convectif,
! si le maillage est mobile (solide rigide)
! En turbomachine, on connait exactement la vitesse de maillage a ajouter
if (iturbo.ne.0) then

  call dmtmps(t3)

  !$omp parallel do private(iel1, iel2, rhofac, vr1, vr2)
  do ifac = 1, nfac
    iel1 = ifacel(1,ifac)
    iel2 = ifacel(2,ifac)
    if (irotce(iel1).ne.0 .or. irotce(iel2).ne.0) then

      rhofac = 0.5d0*(crom(iel1) + crom(iel2))
      call rotation_velocity(irotce(iel1), cdgfac(:,ifac), vr1)
      call rotation_velocity(irotce(iel2), cdgfac(:,ifac), vr2)

      imasfl(ifac) = imasfl(ifac) - 0.5d0 *rhofac*(         &
                          surfac(1,ifac)*(vr1(1) + vr2(1))  &
                        + surfac(2,ifac)*(vr1(2) + vr2(2))  &
                        + surfac(3,ifac)*(vr1(3) + vr2(3)) )
    endif
  enddo
  !$omp parallel do private(iel, rhofac, vr) if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    if (irotce(iel).ne.0) then

      rhofac = brom(ifac)
      call rotation_velocity(irotce(iel), cdgfbo(:,ifac), vr)

      bmasfl(ifac) = bmasfl(ifac) - rhofac*( surfbo(1,ifac)*vr(1) &
                                           + surfbo(2,ifac)*vr(2) &
                                           + surfbo(3,ifac)*vr(3) )
    endif
  enddo

  call dmtmps(t4)

  rs_ell(2) = rs_ell(2) + t4-t3

endif

!===============================================================================
! 10. VoF: void fraction solving and update the mixture density/viscosity
!      and mass flux (cs_pressure_correction solved the convective flux of
!      void fraction, divU)
!===============================================================================

if (ivofmt.gt.0) then

  ! Void fraction solving

  call resvoi(dt, iterns)

  ! Halo synchronization

  call synsca(cvar_voidf)

  ! Update mixture density/viscosity and mass flux

  call vof_update_phys_prop

  ! Verbosity

  if (mod(ntcabs,ntlist).eq.0.and.iterns.eq.nterup) then

    call vof_log_mass_budget

  endif

endif

! Update density (which is coherent with the mass)
!-------------------------------------------------

if (irovar.eq.1.and.(idilat.gt.1.or.ivofmt.gt.0.or.ippmod(icompf).eq.3)) then
  do iel = 1, ncelet
    cpro_rho_mass(iel) = crom_eos(iel)
  enddo

  do ifac = 1, nfabor
    bpro_rho_mass(ifac) = brom_eos(ifac)
  enddo
endif

!===============================================================================
! 11. Compute error estimators for correction step and the global algo
!===============================================================================

if (iestim(iescor).ge.0.or.iestim(iestot).ge.0) then

  ! Allocate temporary arrays
  allocate(esflum(nfac), esflub(nfabor))

  ! ---> ECHANGE DES VITESSES ET PRESSION EN PERIODICITE ET PARALLELISME

  !    Pour les estimateurs IESCOR et IESTOT, la vitesse doit etre echangee.

  !    Pour l'estimateur IESTOT, la pression doit etre echangee aussi.

  !    Cela ne remplace pas l'echange du debut de pas de temps
  !     a cause de cs_user_extra_operations qui vient plus tard et des calculs suite)


  ! --- Vitesse

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvin(vel)
  endif

  !  -- Pression

  if (iescal(iestot).gt.0) then

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(cvar_pr)
    endif

  endif

  ! ---> CALCUL DU FLUX DE MASSE DEDUIT DE LA VITESSE REACTUALISEE

  itypfl = 1
  init   = 1
  inc    = 1
  iflmb0 = 1
  if (iale.ge.1) iflmb0 = 0
  imrgrp = vcopt_u%imrgra
  nswrgp = vcopt_u%nswrgr
  imligp = vcopt_u%imligr
  iwarnp = vcopt_u%iwarni
  epsrgp = vcopt_u%epsrgr
  climgp = vcopt_u%climgr

  call inimav                                                     &
 ( ivarfl(iu)      , itypfl ,                                     &
   iflmb0 , init   , inc    , imrgrp , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom, brom,                                                    &
   vel    ,                                                       &
   coefau , coefbu ,                                              &
   esflum , esflub )


  ! ---> CALCUL DE L'ESTIMATEUR CORRECTION : DIVERGENCE DE ROM * U (N + 1)
  !                                          - GAMMA

  if (iestim(iescor).ge.0) then
    init = 1

    ! Allocate work arrays
    call divmas(init, esflum, esflub, c_estim)

    call field_get_val_s(iestim(iescor), c_estim)

    if (ncetsm.gt.0) then
      !$omp parallel do private(iel) if(ncetsm > thr_n_min)
      do iitsm = 1, ncetsm
        iel = icetsm(iitsm)
        c_estim(iel) = c_estim(iel)-cell_f_vol(iel)*smacel(iitsm,ipr)
      enddo
    endif

    if (iescal(iescor).eq.2) then
      !$omp parallel do
      do iel = 1, ncel
        c_estim(iel) = abs(c_estim(iel))
      enddo
    elseif (iescal(iescor).eq.1) then
      !$omp parallel do
      do iel = 1, ncel
        c_estim(iel) = abs(c_estim(iel)) / volume(iel)
      enddo
    endif

  endif

  ! ---> CALCUL DE L'ESTIMATEUR TOTAL

  if (iestim(iestot).ge.0) then

    !   INITIALISATION DE TRAV AVEC LE TERME INSTATIONNAIRE

    !$omp parallel do private(rovolsdt, isou)
    do iel = 1, ncel
      rovolsdt = crom(iel)*cell_f_vol(iel)/dt(iel)
      do isou = 1, 3
        trav(isou,iel) = rovolsdt * (vela(isou,iel) - vel(isou,iel))
      enddo
    enddo

    !   APPEL A PREDVV AVEC VALEURS AU PAS DE TEMPS COURANT
    iappel = 2
    call predvv &
 ( iappel ,                                                       &
   nvar   , nscal  , iterns ,                                     &
   ncepdc , ncetsm ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   dt     , vel    , vel    , velk   , da_uu  ,                   &
   tslagr , coefau , coefbu , cofafu , cofbfu ,                   &
   ckupdc , smacel , frcxt  ,                                     &
   trava  ,                   dfrcxt , dttens , trav   ,          &
   viscf  , viscb  , viscfi , viscbi , secvif , secvib )

  endif

  deallocate(esflum, esflub)
endif

!===============================================================================
! 12. Velocity/pressure inner iterations
!===============================================================================

if (nterup.gt.1) then

  ! Convergence test on U/P inner iterations, icvrge is 1 if converged
  icvrge = 1

  xnrtmp = 0.d0
  !$omp parallel do reduction(+:xnrtmp) private(xdu, xdv, xdw)
  do iel = 1, ncel
    xdu = vel(1,iel) - velk(1,iel)
    xdv = vel(2,iel) - velk(2,iel)
    xdw = vel(3,iel) - velk(3,iel)
    xnrtmp = xnrtmp +(xdu**2 + xdv**2 + xdw**2) * cell_f_vol(iel)
  enddo
  xnrmu = xnrtmp

  ! parallelism
  if (irangp.ge.0) call parsom (xnrmu)

  ! code-code coupling
  do numcpl = 1, nbrcpl
    call tbrcpl(numcpl, 1, 1, xnrmu, xnrdis)
    xnrmu = xnrmu + xnrdis
  enddo
  xnrmu = sqrt(xnrmu)

  ! Indicateur de convergence du point fixe
  if (xnrmu.ge.epsup*xnrmu0) icvrge = 0

endif

! Shift pressure field to set its spatial mean value to zero
! if there is no boundary faces with a Dirichlet condition on the pressure.
! Number of faces with Dirichlet condition for the pressure is:
! - ndircl if idiricl = 1
! - ndircl-1 if idircl = 0

if (vcopt_p%idircl.eq.1) then
  ndircp = vcopt_p%ndircl
else
  ndircp = vcopt_p%ndircl-1
endif
if (ndircp.le.0) then
  call field_set_volume_average(ivarfl(ipr), pred0)
endif

! Compute the total pressure (defined as a post-processed property).
! For the compressible module, the solved pressure is already the total pressure.
! NB: for Eddy Viscosity Models, TKE might be included in the solved pressure.

if (ippmod(icompf).lt.0) then
  call navstv_total_pressure
endif

!===============================================================================
! 13. Printing
!===============================================================================

if (vcopt_u%iwarni.ge.1) then

  write(nfecra,2000)

  rnorm = -1.d0
  do iel = 1, ncel
    rnorm = max(rnorm,abs(cvar_pr(iel)))
  enddo
  if (irangp.ge.0) call parmax (rnorm)

  write(nfecra,2100)rnorm

  rnorm = -1.d0
  imax = 1
  rnormt = -1.d0
  do iel = 1, ncel
    vitnor = sqrt(vel(1,iel)**2+vel(2,iel)**2+vel(3,iel)**2)
    if (vitnor.ge.rnormt) then
      rnormt = vitnor
      imaxt  = iel
    endif
  enddo
  if (rnormt .gt. rnorm) then
    rnorm = rnormt
    imax = imaxt
  endif

  xyzmax(1) = xyzcen(1,imax)
  xyzmax(2) = xyzcen(2,imax)
  xyzmax(3) = xyzcen(3,imax)

  if (irangp.ge.0) then
    nbrval = 3
    call parmxl (nbrval, rnorm, xyzmax)
  endif

  write(nfecra,2200) rnorm,xyzmax(1),xyzmax(2),xyzmax(3)

  rnorm = sqrt(vel(1,1)**2+vel(2,1)**2+vel(3,1)**2)
  imin = 1
  rnormt = rnorm
  do iel = 1, ncel
    vitnor = sqrt(vel(1,iel)**2+vel(2,iel)**2+vel(3,iel)**2)
    if (vitnor.le.rnormt) then
      rnormt = vitnor
      imint  = iel
    endif
  enddo
  if (rnormt .lt. rnorm) then
    rnorm = rnormt
    imin = imint
  endif

  xyzmin(1) = xyzcen(1,imin)
  xyzmin(2) = xyzcen(2,imin)
  xyzmin(3) = xyzcen(3,imin)

  if (irangp.ge.0) then
    nbrval = 3
    call parmnl (nbrval, rnorm, xyzmin)
  endif

  write(nfecra,2201) rnorm,xyzmin(1),xyzmin(2),xyzmin(3)


  ! With porosity
  if (iporos.ge.1) then
    call field_get_val_s(ipori, porosi)
    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(porosi)
    endif
  endif

  if (ivofmt.gt.0) then
    ! Id of the volume flux
    call field_get_key_int(ivarfl(ivolf2), kimasf, ivolfl_id)
    call field_get_key_int(ivarfl(ivolf2), kbmasf, bvolfl_id)

    ! Pointers to the mass fluxes
    call field_get_val_s(ivolfl_id, ivolfl)
    call field_get_val_s(bvolfl_id, bvolfl)
  endif

  rnorma = -grand
  rnormi =  grand
  do ifac = 1, nfac
    iel1 = ifacel(1,ifac)
    iel2 = ifacel(2,ifac)
    if (iporos.eq.1.or.iporos.eq.2) then
      rhom = (porosi(iel1)*crom(iel1)+porosi(iel2)*crom(iel2))*0.5d0
    else
      rhom = (crom(iel1)+crom(iel2))*0.5d0
    endif
    ! Deal with null fluid section
    if (suffan(ifac)/surfan(ifac).gt.epzero) then
      if (ivofmt.gt.0) then
        rnorm = abs(ivolfl(ifac))/suffan(ifac)
      else
        rnorm = abs(imasfl(ifac))/(suffan(ifac)*rhom)
      endif
    else
      rnorm = 0.d0
    endif
    rnorma = max(rnorma,rnorm)
    rnormi = min(rnormi,rnorm)
  enddo
  if (irangp.ge.0) then
    call parmax (rnorma)
    call parmin (rnormi)
  endif
  write(nfecra,2300)rnorma, rnormi
  rnorma = -grand
  rnormi =  grand
  do ifac = 1, nfabor
    if (ivofmt.gt.0) then
      ! Deal with null fluid section
      if (suffbn(ifac)/surfbn(ifac).gt.epzero) then
        rnorm = bvolfl(ifac)/(suffbn(ifac))
      else
        rnorm = 0.d0
      endif
    else
      if (iporos.eq.1.or.iporos.eq.2) then
        rnorm = bmasfl(ifac)/(surfbn(ifac)*brom(ifac)*porosi(ifabor(ifac)))
      else
        ! Deal with null fluid section
        if (suffbn(ifac)/surfbn(ifac).gt.epzero) then
          rnorm = bmasfl(ifac)/(suffbn(ifac)*brom(ifac))
        else
          rnorm = 0.d0
        endif
      endif
    endif
    rnorma = max(rnorma,rnorm)
    rnormi = min(rnormi,rnorm)
  enddo
  if (irangp.ge.0) then
    call parmax (rnorma)
    call parmin (rnormi)
  endif
  write(nfecra,2400)rnorma, rnormi

  rnorm = 0.d0
  !$omp parallel do reduction(+: rnorm) if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    rnorm = rnorm + bmasfl(ifac)
  enddo

  if (irangp.ge.0) call parsom (rnorm)

  write(nfecra,2500)rnorm

  write(nfecra,2001)

  if (nterup.gt.1) then
    if (icvrge.eq.0) then
      write(nfecra,2600) iterns
      write(nfecra,2601) xnrmu, xnrmu0, epsup
      write(nfecra,2001)
      if (iterns.eq.nterup) then
        write(nfecra,2603)
        write(nfecra,2001)
      endif
    else
      write(nfecra,2602) iterns
      write(nfecra,2601) xnrmu, xnrmu0, epsup
      write(nfecra,2001)
    endif
  endif

endif

if (iturbo.eq.2) then
  if (mod(ntcabs,ntlist).eq.0.and.iterns.eq.nterup) then
     write(nfecra,3000) rs_ell(1), rs_ell(1) + rs_ell(2)
     rs_ell(1) = 0.d0
     rs_ell(2) = 0.d0
  endif
endif

! Free memory
deallocate(viscf, viscb)
deallocate(trav)
deallocate(da_uu)
deallocate(dfrcxt)
if (allocated(wvisfi)) deallocate(wvisfi, wvisbi)
if (allocated(uvwk)) deallocate(uvwk)
if (allocated(secvif)) deallocate(secvif, secvib)
if (allocated(cpro_rho_tc)) deallocate(cpro_rho_tc)
if (allocated(bpro_rho_tc)) deallocate(bpro_rho_tc)

!--------
! Formats
!--------

 1000 format(/,                                                   &
'   ** SOLVING VELOCITY'                                       ,/,&
'      ----------------'                                       ,/)
 1080 format(/,                                                   &
'   ** SOLVING MASS BALANCE EQUATION                          ',/,&
'      -----------------------------                          ',/)
 1200 format(/,                                                   &
'   ** SOLVING CONTINUITY PRESSURE'                            ,/,&
'      ---------------------------'                            ,/)
 2000 format(/,' AFTER CONTINUITY PRESSURE',/,                    &
'-------------------------------------------------------------'  )
 2100 format(                                                           &
' Max. pressure',E12.4   ,' (max. absolute value)'             ,/)
 2200 format(                                                           &
' Max. velocity',E12.4   ,' in',3E11.3                         ,/)
 2201 format(                                                           &
' Min. velocity',E12.4   ,' in',3E11.3                         ,/)
 2300 format(                                                           &
' Max. velocity at interior face',E12.4   ,' ; min.',E12.4       )
 2400 format(                                                           &
' Max. velocity at boundary face',E12.4   ,' ; min.',E12.4       )
 2500 format(                                                           &
' Mass balance  at boundary  ',E14.6                             )
 2600 format(                                                           &
' Fixed point informations at iteration:',I10                  ,/)
 2601 format('norm = ',E12.4,' norm 0 = ',E12.4,' toler  = ',E12.4   ,/)
 2602 format(                                                           &
' Fixed point convergence at iteration ',I10                   ,/)
 2603 format(                                                           &
' Non convergence of fixed point for velocity pressure coupling' )
 2001 format(                                                           &
'-------------------------------------------------------------',/)
 3000 format(/,                                             &
'   ** INFORMATION ON UNSTEADY ROTOR/STATOR TREATMENT',/,&
'      ----------------------------------------------',/,&
' Time dedicated to mesh update (s):',F12.4,           /,&
' Global time                   (s):',F12.4,           /)

!----
! End
!----

return

end subroutine

!===============================================================================
! Local functions
!===============================================================================

!===============================================================================
! Function:
! ---------

!> \brief Update total pressure (defined as a post-processed property).
!
!> For the compressible module, the solved pressure is already
!> the total pressure.
!
!> Note: for Eddy Viscosity Models, the TKE may be included in the
!> solved pressure.

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________


subroutine navstv_total_pressure

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use atincl, only: iatmst, imomst
use numvar
use cstphy
use cstnum
use optcal
use pointe
use parall
use paramx, only: isymet
use ppincl
use mesh
use ptrglo
use field

!===============================================================================

implicit none

! Local variables

integer          iel
double precision xxp0 , xyp0 , xzp0

double precision, dimension(:), pointer :: cpro_rho
double precision, dimension(:,:), pointer :: cpro_momst
double precision, dimension(:), pointer :: cvar_pr
double precision, dimension(:), pointer :: cpro_prtot
double precision, dimension(:), pointer :: cvara_k

!===============================================================================
! 0. Initialization
!===============================================================================

if (ipr.lt.1 .or. iprtot.lt.0) return

call field_get_val_s(iprtot, cpro_prtot)
call field_get_val_s(ivarfl(ipr), cvar_pr)

xxp0   = xyzp0(1)
xyp0   = xyzp0(2)
xzp0   = xyzp0(3)

if (iatmst.eq.0) then
  do iel=1,ncel
    cpro_prtot(iel) =  cvar_pr(iel)                        &
                     + ro0*(  gx*(xyzcen(1,iel)-xxp0)      &
                            + gy*(xyzcen(2,iel)-xyp0)      &
                            + gz*(xyzcen(3,iel)-xzp0) )    &
                     + p0 - pred0
  enddo
else
  call field_get_val_v(imomst, cpro_momst)

  do iel=1,ncel
    cpro_prtot(iel) =   cvar_pr(iel)                             &
                      + ro0*(  gx*(xyzcen(1,iel)-xxp0)           &
                             + gy*(xyzcen(2,iel)-xyp0)           &
                             + gz*(xyzcen(3,iel)-xzp0))          &
                      + p0 - pred0                               &
                      - cpro_momst(1,iel)*(xyzcen(1,iel)-xxp0)   &
                      - cpro_momst(2,iel)*(xyzcen(2,iel)-xyp0)   &
                      - cpro_momst(3,iel)*(xyzcen(3,iel)-xzp0)
  enddo
endif

! For Eddy Viscosity Models, "2/3 rho k" is included in the solved pressure
if ((itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60).and. igrhok.ne.1) then
  call field_get_val_s(ivarfl(ik), cvara_k)
  call field_get_val_s(icrom, cpro_rho)
  do iel = 1, ncel
    cpro_prtot(iel) =   cpro_prtot(iel)                              &
                      - 2.d0 / 3.d0 * cpro_rho(iel) * cvara_k(iel)
  enddo
endif

!----
! End
!----

return

end subroutine navstv_total_pressure
