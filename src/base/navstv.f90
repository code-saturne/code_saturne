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
use atincl, only: iatmst, imomst, iautom, imeteo

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
integer          ii    , inod, itypfl, f_id
integer          isou, ivar, iitsm
integer          init, iautof
integer          iflmas, iflmab
integer          iflmb0
integer          nswrgp, imligp, iwarnp
integer          nbrval, iappel
integer          ndircp, icpt
integer          numcpl
integer          iflvoi, iflvob
double precision rnorm , rnormt, rnorma, rnormi, vitnor
double precision dtsrom, unsrom, rhom, rovolsdt
double precision epsrgp, climgp, extrap, xyzmax(3), xyzmin(3)
double precision thetap, xdu, xdv, xdw
double precision xxp0 , xyp0 , xzp0
double precision rhofac, dtfac
double precision xnrdis, xnrtmp
double precision t1, t2, t3, t4
double precision visclc, visctc
double precision distbf, srfbnf, hint
double precision rnx, rny, rnz
double precision vr(3), vr1(3), vr2(3), vrn
double precision disp_fac(3)
double precision mass_fl_drhovol1, mass_fl_drhovol2

double precision, allocatable, dimension(:,:,:), target :: viscf
double precision, allocatable, dimension(:), target :: viscb
double precision, allocatable, dimension(:,:,:), target :: wvisfi
double precision, allocatable, dimension(:,:), target :: uvwk
double precision, dimension(:,:), pointer :: velk
double precision, allocatable, dimension(:), target :: wvisbi
double precision, allocatable, dimension(:), target :: cpro_rho_tc, bpro_rho_tc
double precision, allocatable, dimension(:) :: phi
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: esflum, esflub
double precision, allocatable, dimension(:) :: intflx, bouflx
double precision, allocatable, dimension(:) :: secvif, secvib

double precision, dimension(:,:), allocatable :: gradp
double precision, dimension(:), allocatable :: coefa_dp, coefb_dp
double precision, dimension(:), allocatable :: xinvro
double precision, dimension(:,:), pointer :: grdphd
double precision, dimension(:,:), pointer :: vel, vela
double precision, dimension(:,:,:), pointer :: viscfi
double precision, dimension(:), pointer :: viscbi
double precision, dimension(:,:), pointer :: dttens
double precision, dimension(:,:), pointer :: dfrcxt
double precision, dimension(:,:), pointer :: coefau, cofafu, claale
double precision, dimension(:,:,:), pointer :: coefbu, cofbfu, clbale
double precision, dimension(:), pointer :: coefa_p
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: brom, broma, crom, croma, viscl, visct
double precision, dimension(:), pointer :: ivoifl, bvoifl
double precision, dimension(:), pointer :: coavoi, cobvoi
double precision, dimension(:,:), pointer :: trav
double precision, dimension(:,:), pointer :: mshvel
double precision, dimension(:,:), pointer :: disale
double precision, dimension(:,:), pointer :: disala
double precision, dimension(:,:), pointer :: cpro_momst
double precision, dimension(:), pointer :: porosi
double precision, dimension(:), pointer :: cvar_pr
double precision, dimension(:), pointer :: cpro_prtot, c_estim
double precision, dimension(:), pointer :: cvar_voidf, cvara_voidf
double precision, dimension(:), pointer :: cvara_k
double precision, dimension(:), pointer :: cpro_rho_mass
double precision, dimension(:), pointer :: bpro_rho_mass
double precision, dimension(:), pointer :: brom_eos, crom_eos

type(var_cal_opt) :: vcopt_p, vcopt_u, vcopt

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine resopv &
   ( nvar   , iterns , ncesmp , nfbpcd , ncmast ,                 &
     icetsm , ifbpcd , ltmast , isostd ,                          &
     dt     , vel    ,                                            &
     coefav , coefbv , coefa_dp        , coefb_dp ,               &
     smacel , spcond , svcond ,                                   &
     frcxt  , dfrcxt , tpucou ,                                   &
     viscf  , viscb  ,                                            &
     phi    , tslagr )

    use dimens, only: ndimfb
    use mesh

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
    double precision viscf(nfac), viscb(ndimfb)
    double precision phi(ncelet)
    double precision tslagr(ncelet,*)
    double precision coefav(3  ,ndimfb)
    double precision coefbv(3,3,ndimfb)
    double precision vel   (3  ,ncelet)
    double precision coefa_dp(ndimfb)
    double precision coefb_dp(ndimfb)

  end subroutine resopv

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

! Array for delta p gradient boundary conditions
allocate(coefa_dp(ndimfb), coefb_dp(ndimfb))

allocate(dfrcxt(3,ncelet))
if (iphydr.eq.2) then
  allocate(grdphd(ndim, ncelet))
else
  grdphd => rvoid2
endif
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

if (ivofmt.ge.0) then
  call field_get_val_s(ivarfl(ivolf2), cvar_voidf)
  call field_get_val_prev_s(ivarfl(ivolf2), cvara_voidf)
endif

! Allocate work arrays
allocate(w1(ncelet))

! Initialize variables to avoid compiler warnings

ivar = 0
iflmas = 0
imax = 0

! pointer to velosity at sub iteration k for PISO like algorithm
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
    ! En cas de couplage entre deux instances de Code_Saturne, on calcule
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

if (irovar.eq.1) then
  ! If iterns = 1: this is density at time n
  call field_get_id("density_mass", f_id)
  call field_get_val_s(f_id, cpro_rho_mass)
  call field_get_id("boundary_density_mass", f_id)
  call field_get_val_s(f_id, bpro_rho_mass)

  ! Time interpolated density
  if (vcopt_u%thetav .lt. 1.d0) then
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
! 2. Hydrostatic pressure prediction in case of Low Mach compressible algorithm
!===============================================================================

if (iphydr.eq.2) then

  call prehyd(grdphd, iterns)

endif

!===============================================================================
! 3. Pressure resolution and computation of mass flux for compressible flow
!===============================================================================

! Note, for the compressible algorithm written in pressure increment,
! this step is performed in the same time as the incompressible algorithm
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
! 4. Compute liquid-vapour mass transfer term for cavitating flows
!===============================================================================

if (icavit.ge.0) then

  call field_get_val_s(iprtot, cpro_prtot)

  call cavitation_compute_source_term (cpro_prtot, cvara_voidf)

endif

!===============================================================================
! 5. Velocity prediction step
!===============================================================================

if (vcopt_u%iwarni.ge.1) then
  write(nfecra,1000)
endif

iappel = 1

call predvv &
( iappel ,                                                       &
  nvar   , nscal  , iterns ,                                     &
  ncepdc , ncetsm ,                                              &
  icepdc , icetsm , itypsm ,                                     &
  dt     , vel    , vela   , velk   ,                            &
  tslagr , coefau , coefbu , cofafu , cofbfu ,                   &
  ckupdc , smacel , frcxt  , grdphd ,                            &
  trava  ,                   dfrcxt , dttens ,  trav  ,          &
  viscf  , viscb  , viscfi , viscbi , secvif , secvib ,          &
  w1     )

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
  nswrgp = vcopt_u%nswrgr
  imligp = vcopt_u%imligr
  iwarnp = vcopt_u%iwarni
  epsrgp = vcopt_u%epsrgr
  climgp = vcopt_u%climgr

  call inimav                                                     &
 ( ivarfl(iu)      , itypfl ,                                     &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
   iwarnp ,                                                       &
   epsrgp , climgp ,                                              &
   crom   , brom   ,                                              &
   vel    ,                                                       &
   coefau , coefbu ,                                              &
   imasfl , bmasfl )

  ! In the ALE framework, we add the mesh velocity
  if (iale.ge.1) then

    call field_get_val_v(ivarfl(iuma), mshvel)

    call field_get_val_prev_v(fdiale, disala)

    if (iflxmw.gt.0) then
      ! One temporary array needed for internal faces, in case some internal vertices
      !  are moved directly by the user
      allocate(intflx(nfac), bouflx(ndimfb))

      itypfl = 1
      init   = 1
      inc    = 1
      iflmb0 = 1
      call field_get_key_struct_var_cal_opt(ivarfl(iuma), vcopt)
      nswrgp = vcopt%nswrgr
      imligp = vcopt%imligr
      iwarnp = vcopt%iwarni
      epsrgp = vcopt%epsrgr
      climgp = vcopt%climgr

      call field_get_coefa_v(ivarfl(iuma), claale)
      call field_get_coefb_v(ivarfl(iuma), clbale)

      call inimav &
    ( ivarfl(iuma)    , itypfl ,                                     &
      iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
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
        do ii = ipnfbr(ifac),ipnfbr(ifac+1)-1
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
    !$omp parallel do private(iel1, iel2, dtfac, rhofac )
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
    !$omp parallel do private(iel, dtfac, rhofac) &
    !$omp          if(nfabor > thr_n_min)
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
  !--------------
  deallocate(coefa_dp, coefb_dp)
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

      ! Scratch and resize work arrays

      deallocate(w1)
      allocate(w1(ncelet))

      ! Resize auxiliary arrays (pointe module)

      call resize_aux_arrays

      ! Update turbomachinery module

      call turbomachinery_update

      ! Update field mappings ("owner" fields handled by turbomachinery_update)

      call fldtri
      call field_get_val_s_by_name('dt', dt)

      ! Resize other arrays related to the velocity-pressure resolution

      call resize_vec_real_array(trav)
      call resize_vec_real_array(dfrcxt)

      ! Resize other arrays, depending on user options

      if (iilagr.gt.0 .and. ntersl.gt.0) &
        call resize_n_sca_real_arrays(ntersl, tslagr)

      if (iphydr.eq.1) then
        call field_get_val_v_by_name('volume_forces', frcxt)
      elseif (iphydr.eq.2) then
        call resize_vec_real_array(grdphd)
      endif

      ! Update local pointers on "cells" fields

      call field_get_val_s(icrom, crom)

      call field_get_val_s(iviscl, viscl)
      call field_get_val_s(ivisct, visct)

      call field_get_val_v(ivarfl(iu), vel)
      call field_get_val_prev_v(ivarfl(iu), vela)

      call field_get_val_s(ivarfl(ipr), cvar_pr)

      if (idtten.ge.0) call field_get_val_v(idtten, dttens)

      if (ivofmt.ge.0) then
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

! Allocate temporary arrays for the pressure resolution
allocate(phi(ncelet))

if (ippmod(icompf).lt.0.or.ippmod(icompf).eq.3) then

  call resopv &
( nvar   , iterns , ncetsm , nfbpcd , ncmast ,                   &
  icetsm , ifbpcd , ltmast , isostd ,                            &
  dt     , vel    ,                                              &
  coefau , coefbu , coefa_dp        , coefb_dp ,                 &
  smacel , spcond , svcond ,                                     &
  frcxt  , dfrcxt , dttens ,                                     &
  viscf  , viscb  ,                                              &
  phi    , tslagr )

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

    ! Phi is the pressure increment

    ! ---> Periodicity and parallelism
    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(phi)
    endif

    iccocg = 1
    inc = 0

    call grdpor(inc)

    if (iphydr.eq.1.or.iifren.eq.1) inc = 1
    nswrgp = vcopt_p%nswrgr
    imligp = vcopt_p%imligr
    iwarnp = vcopt_p%iwarni
    epsrgp = vcopt_p%epsrgr
    climgp = vcopt_p%climgr
    extrap = vcopt_p%extrag

    !Allocation
    allocate(gradp(3, ncelet))

    if (ivofmt.lt.0) then
      call gradient_potential_s &
       (ivarfl(ipr)     , imrgra , inc    , iccocg , nswrgp , imligp , &
        iphydr , iwarnp , epsrgp , climgp , extrap ,                   &
        dfrcxt ,                                                       &
        phi    , coefa_dp        , coefb_dp        ,                   &
        gradp  )
    else
      allocate(xinvro(ncelet))
      do iel = 1, ncel
        xinvro(iel) = 1.d0/crom(iel)
      enddo

      call gradient_weighted_s &
      ( ivarfl(ipr)     , imrgra , inc    , iccocg , nswrgp , imligp , &
        iphydr, iwarnp  , epsrgp , climgp , extrap , dfrcxt ,          &
        phi    , xinvro , coefa_dp , coefb_dp ,                        &
        gradp  )

      deallocate(xinvro)
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
                 + dtsrom*(dfrcxt(isou, iel)-gradp(isou,iel))
          enddo
        enddo

      ! Tensorial diffusion for the pressure
      elseif (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then
        !$omp parallel do private(unsrom)
        do iel = 1, ncel
          unsrom = thetap/crom(iel)

          vel(1, iel) = vel(1, iel)                              &
               + unsrom*(                                        &
                 dttens(1,iel)*(dfrcxt(1, iel)-gradp(1,iel))     &
               + dttens(4,iel)*(dfrcxt(2, iel)-gradp(2,iel))     &
               + dttens(6,iel)*(dfrcxt(3, iel)-gradp(3,iel))     &
               )
          vel(2, iel) = vel(2, iel)                              &
               + unsrom*(                                        &
                 dttens(4,iel)*(dfrcxt(1, iel)-gradp(1,iel))     &
               + dttens(2,iel)*(dfrcxt(2, iel)-gradp(2,iel))     &
               + dttens(5,iel)*(dfrcxt(3, iel)-gradp(3,iel))     &
               )
          vel(3, iel) = vel(3, iel)                              &
               + unsrom*(                                        &
                 dttens(6,iel)*(dfrcxt(1 ,iel)-gradp(1,iel))     &
               + dttens(5,iel)*(dfrcxt(2 ,iel)-gradp(2,iel))     &
               + dttens(3,iel)*(dfrcxt(3 ,iel)-gradp(3,iel))     &
               )
        enddo
      endif

      ! Update external forces for the computation of the gradients
      !$omp parallel do
      do iel=1,ncel
        frcxt(1 ,iel) = frcxt(1 ,iel) + dfrcxt(1 ,iel)
        frcxt(2 ,iel) = frcxt(2 ,iel) + dfrcxt(2 ,iel)
        frcxt(3 ,iel) = frcxt(3 ,iel) + dfrcxt(3 ,iel)
      enddo
      if (irangp.ge.0.or.iperio.eq.1) then
        call synvin(frcxt)
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

        if (isostd(ifac).eq.1.or.iatmst.eq.1.and.iautof.eq.1) then
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
            vel(isou,iel) = vel(isou,iel) - dtsrom*gradp(isou,iel)
          enddo
        enddo

      ! Tensorial diffusion for the pressure
      else if (iand(vcopt_p%idften, ANISOTROPIC_DIFFUSION).ne.0) then

        !$omp parallel do private(unsrom)
        do iel = 1, ncel
          unsrom = thetap/crom(iel)

          vel(1, iel) = vel(1, iel)                              &
                      - unsrom*(                                 &
                                 dttens(1,iel)*(gradp(1,iel))    &
                               + dttens(4,iel)*(gradp(2,iel))    &
                               + dttens(6,iel)*(gradp(3,iel))    &
                               )
          vel(2, iel) = vel(2, iel)                              &
                      - unsrom*(                                 &
                                 dttens(4,iel)*(gradp(1,iel))    &
                               + dttens(2,iel)*(gradp(2,iel))    &
                               + dttens(5,iel)*(gradp(3,iel))    &
                               )
          vel(3, iel) = vel(3, iel)                              &
                      - unsrom*(                                 &
                                 dttens(6,iel)*(gradp(1,iel))    &
                               + dttens(5,iel)*(gradp(2,iel))    &
                               + dttens(3,iel)*(gradp(3,iel))    &
                               )
        enddo

      endif
    endif

    !Free memory
    deallocate(gradp)

    ! RT0 update from the mass fluxes
  else

    ! Initialization to 0
    do iel = 1, ncelet
      vel(1, iel) = 0.d0
      vel(2, iel) = 0.d0
      vel(3, iel) = 0.d0
    enddo

    ! vel = 1 / (rho Vol) SUM mass_flux (X_f - X_i)
    if (ivofmt.lt.0) then
      do ifac = 1, nfac

        iel1 = ifacel(1,ifac)
        iel2 = ifacel(2,ifac)

        mass_fl_drhovol1 = 0.d0
        ! If it is not a solid cell
        if (isolid(iporos, iel1).eq.0) &
          mass_fl_drhovol1 = imasfl(ifac) / (crom(iel1) * cell_f_vol(iel1))


        mass_fl_drhovol2 = 0.d0
        ! If it is not a solid cell
        if (isolid(iporos, iel2).eq.0) &
          mass_fl_drhovol2 = imasfl(ifac) / (crom(iel2) * cell_f_vol(iel2))

        vel(1, iel1) = vel(1, iel1) &
          + mass_fl_drhovol1 * (cdgfac(1, ifac) - xyzcen(1, iel1))
        vel(2, iel1) = vel(2, iel1) &
          + mass_fl_drhovol1 * (cdgfac(2, ifac) - xyzcen(2, iel1))
        vel(3, iel1) = vel(3, iel1) &
          + mass_fl_drhovol1 * (cdgfac(3, ifac) - xyzcen(3, iel1))

        vel(1, iel2) = vel(1, iel2) &
          - mass_fl_drhovol2 * (cdgfac(1, ifac) - xyzcen(1, iel2))
        vel(2, iel2) = vel(2, iel2) &
          - mass_fl_drhovol2 * (cdgfac(2, ifac) - xyzcen(2, iel2))
        vel(3, iel2) = vel(3, iel2) &
          - mass_fl_drhovol2 * (cdgfac(3, ifac) - xyzcen(3, iel2))

      enddo

      do ifac = 1, nfabor
        iel1 = ifabor(ifac)

        mass_fl_drhovol1 = 0.d0
        ! If it is not a solid cell
        if (isolid(iporos, iel1).eq.0) &
          mass_fl_drhovol1 = bmasfl(ifac) / (crom(iel1) * cell_f_vol(iel1))

        vel(1, iel1) = vel(1, iel1) &
          + mass_fl_drhovol1 * (cdgfbo(1, ifac) - xyzcen(1, iel1))
        vel(2, iel1) = vel(2, iel1) &
          + mass_fl_drhovol1 * (cdgfbo(2, ifac) - xyzcen(2, iel1))
        vel(3, iel1) = vel(3, iel1) &
          + mass_fl_drhovol1 * (cdgfbo(3, ifac) - xyzcen(3, iel1))

      enddo

      ! vel = 1 / (Vol) SUM vol_flux (X_f - X_i)
    else

      do ifac = 1, nfac

        iel1 = ifacel(1,ifac)
        iel2 = ifacel(2,ifac)

        mass_fl_drhovol1 = 0.d0
        ! If it is not a solid cell
        if (isolid(iporos, iel1).eq.0) &
          mass_fl_drhovol1 = imasfl(ifac) / cell_f_vol(iel1)


        mass_fl_drhovol2 = 0.d0
        ! If it is not a solid cell
        if (isolid(iporos, iel2).eq.0) &
          mass_fl_drhovol2 = imasfl(ifac) / cell_f_vol(iel2)

        vel(1, iel1) = vel(1, iel1) &
          + mass_fl_drhovol1 * (cdgfac(1, ifac) - xyzcen(1, iel1))
        vel(2, iel1) = vel(2, iel1) &
          + mass_fl_drhovol1 * (cdgfac(2, ifac) - xyzcen(2, iel1))
        vel(3, iel1) = vel(3, iel1) &
          + mass_fl_drhovol1 * (cdgfac(3, ifac) - xyzcen(3, iel1))

        vel(1, iel2) = vel(1, iel2) &
          - mass_fl_drhovol2 * (cdgfac(1, ifac) - xyzcen(1, iel2))
        vel(2, iel2) = vel(2, iel2) &
          - mass_fl_drhovol2 * (cdgfac(2, ifac) - xyzcen(2, iel2))
        vel(3, iel2) = vel(3, iel2) &
          - mass_fl_drhovol2 * (cdgfac(3, ifac) - xyzcen(3, iel2))

      enddo

      do ifac = 1, nfabor
        iel1 = ifabor(ifac)

        mass_fl_drhovol1 = 0.d0
        ! If it is not a solid cell
        if (isolid(iporos, iel1).eq.0) &
          mass_fl_drhovol1 = bmasfl(ifac) / cell_f_vol(iel1)

        vel(1, iel1) = vel(1, iel1) &
          + mass_fl_drhovol1 * (cdgfbo(1, ifac) - xyzcen(1, iel1))
        vel(2, iel1) = vel(2, iel1) &
          + mass_fl_drhovol1 * (cdgfbo(2, ifac) - xyzcen(2, iel1))
        vel(3, iel1) = vel(3, iel1) &
          + mass_fl_drhovol1 * (cdgfbo(3, ifac) - xyzcen(3, iel1))

      enddo

    endif

    call synvin(vel)

  endif

endif

! Bad cells regularisation
call cs_bad_cells_regularisation_vector(vel, 1)

! Mass flux initialization for VOF algorithm
if (ivofmt.ge.0) then
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
  call field_get_val_prev_v(fdiale, disala)

  if (iflxmw.gt.0) then
    ! One temporary array needed for internal faces, in case some internal vertices
    !  are moved directly by the user
    allocate(intflx(nfac), bouflx(ndimfb))

    itypfl = 1
    init   = 1
    inc    = 1
    iflmb0 = 1
    call field_get_key_struct_var_cal_opt(ivarfl(iuma), vcopt)
    nswrgp = vcopt%nswrgr
    imligp = vcopt%imligr
    iwarnp = vcopt%iwarni
    epsrgp = vcopt%epsrgr
    climgp = vcopt%climgr

    call field_get_coefa_v(ivarfl(iuma), claale)
    call field_get_coefb_v(ivarfl(iuma), clbale)

    call inimav &
  ( ivarfl(iuma)    , itypfl ,                                     &
    iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
    iwarnp ,                                                       &
    epsrgp , climgp ,                                              &
    crom, brom,                                                    &
    mshvel ,                                                       &
    claale , clbale ,                                              &
    intflx , bouflx )
  endif

  ! Here we need of the opposite of the mesh velocity.
  !$omp parallel do if(nfabor > thr_n_min)
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
!      and mass flux (resopv solved the convective flux of void fraction, divU)
!===============================================================================

if (ivofmt.ge.0) then

  ! Void fraction solving

  call resvoi(dt, iterns)

  ! Halo synchronization

  call synsca(cvar_voidf)

  ! Get the void fraction boundary conditions

  call field_get_coefa_s(ivarfl(ivolf2), coavoi)
  call field_get_coefb_s(ivarfl(ivolf2), cobvoi)

  ! Get the convective flux of the void fraction

  call field_get_key_int(ivarfl(ivolf2), kimasf, iflvoi)
  call field_get_key_int(ivarfl(ivolf2), kbmasf, iflvob)
  call field_get_val_s(iflvoi, ivoifl)
  call field_get_val_s(iflvob, bvoifl)

  ! Update mixture density/viscosity and mass flux

  call vof_update_phys_prop &
 ( cvar_voidf, coavoi, cobvoi, ivoifl, bvoifl, &
   crom, brom, imasfl, bmasfl )

  ! Verbosity

  if (mod(ntcabs,ntlist).eq.0.and.iterns.eq.nterup) then

    call field_get_val_prev_s(icrom, croma)

    call vof_print_mass_budget &
   ( crom, croma, brom, dt, imasfl, bmasfl )

  endif

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
  nswrgp = vcopt_u%nswrgr
  imligp = vcopt_u%imligr
  iwarnp = vcopt_u%iwarni
  epsrgp = vcopt_u%epsrgr
  climgp = vcopt_u%climgr

  call inimav                                                     &
 ( ivarfl(iu)      , itypfl ,                                     &
   iflmb0 , init   , inc    , imrgra , nswrgp , imligp ,          &
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
    call divmas(init, esflum, esflub, w1)

    call field_get_val_s(iestim(iescor), c_estim)

    if (ncetsm.gt.0) then
      !$omp parallel do private(iel) if(ncetsm > thr_n_min)
      do iitsm = 1, ncetsm
        iel = icetsm(iitsm)
        w1(iel) = w1(iel)-cell_f_vol(iel)*smacel(iitsm,ipr)
      enddo
    endif

    if (iescal(iescor).eq.2) then
      !$omp parallel do
      do iel = 1, ncel
        c_estim(iel) =  abs(w1(iel))
      enddo
    elseif (iescal(iescor).eq.1) then
      !$omp parallel do
      do iel = 1, ncel
        c_estim(iel) =  abs(w1(iel)) / volume(iel)
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
   dt     , vel    , vel    , velk   ,                            &
   tslagr , coefau , coefbu , cofafu , cofbfu ,                   &
   ckupdc , smacel , frcxt  , grdphd ,                            &
   trava  ,                   dfrcxt , dttens , trav   ,          &
   viscf  , viscb  , viscfi , viscbi , secvif , secvib ,          &
   w1     )

  endif

  deallocate(esflum, esflub)
endif

!===============================================================================
! 12. Loop on the velocity/Pressure coupling (PISO)
!===============================================================================

if (nterup.gt.1) then

  ! Convergence test on PISO-like algorithm, icvrge is 1 if converged
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
  call field_get_val_s(iprtot, cpro_prtot)
  xxp0   = xyzp0(1)
  xyp0   = xyzp0(2)
  xzp0   = xyzp0(3)

  if (iatmst.eq.0) then
    do iel=1,ncel
      cpro_prtot(iel)= cvar_pr(iel)               &
           + ro0*( gx*(xyzcen(1,iel)-xxp0)               &
           + gy*(xyzcen(2,iel)-xyp0)                     &
           + gz*(xyzcen(3,iel)-xzp0) )                   &
           + p0 - pred0
    enddo
  else
    call field_get_val_v(imomst, cpro_momst)

    do iel=1,ncel
      cpro_prtot(iel)= cvar_pr(iel)                      &
           + ro0*(gx*(xyzcen(1,iel)-xxp0)                &
           + gy*(xyzcen(2,iel)-xyp0)                     &
           + gz*(xyzcen(3,iel)-xzp0))                    &
           + p0 - pred0                                  &
           - cpro_momst(1,iel)*(xyzcen(1,iel)-xxp0)      &
           - cpro_momst(2,iel)*(xyzcen(2,iel)-xyp0)      &
           - cpro_momst(3,iel)*(xyzcen(3,iel)-xzp0)
    enddo
  endif
  ! For Eddy Viscosity Models, "2/3 rho k" is included in the solved pressure
  if ((itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60).and. igrhok.ne.1) then
    call field_get_val_s(ivarfl(ik), cvara_k)
    do iel = 1, ncel
      cpro_prtot(iel) = cpro_prtot(iel) - 2.d0 / 3.d0 * crom(iel) * cvara_k(iel)
    enddo
  endif
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
      rnorm = abs(imasfl(ifac))/(suffan(ifac)*rhom)
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
deallocate(phi)
deallocate(trav)
deallocate(dfrcxt)
deallocate(w1)
if (allocated(wvisfi)) deallocate(wvisfi, wvisbi)
if (allocated(uvwk)) deallocate(uvwk)
if (allocated(secvif)) deallocate(secvif, secvib)
if (allocated(cpro_rho_tc)) deallocate(cpro_rho_tc)
if (allocated(bpro_rho_tc)) deallocate(bpro_rho_tc)
if (iphydr.eq.2) deallocate(grdphd)

! Free memory
!--------------
deallocate(coefa_dp, coefb_dp)

!--------
! Formats
!--------
#if defined(_CS_LANG_FR)

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VITESSE                             ',/,&
'      --------------------------                             ',/)
 1080 format(/,                                                   &
'   ** RESOLUTION DE L''EQUATION DE MASSE                     ',/,&
'      ----------------------------------                     ',/)
 1200 format(/,                                                   &
'   ** RESOLUTION POUR LA PRESSION CONTINUITE                 ',/,&
'      --------------------------------------                 ',/)
 2000 format(/,' APRES PRESSION CONTINUITE',/,                    &
'-------------------------------------------------------------'  )
 2100 format(                                                           &
' Pression max.',E12.4   ,' (max. de la valeur absolue)       ',/)
 2200 format(                                                           &
' Vitesse  max.',E12.4   ,' en',3E11.3                         ,/)
 2201 format(                                                           &
' Vitesse  min.',E12.4   ,' en',3E11.3                         ,/)
 2300 format(                                                           &
' Vitesse  en face interne max.',E12.4   ,' ; min.',E12.4        )
 2400 format(                                                           &
' Vitesse  en face de bord max.',E12.4   ,' ; min.',E12.4        )
 2500 format(                                                           &
' Bilan de masse   au bord   ',E14.6                             )
 2600 format(                                                           &
' Informations Point fixe a l''iteration :',I10                ,/)
 2601 format('norme = ',E12.4,' norme 0 = ',E12.4,' toler  = ',E12.4 ,/)
 2602 format(                                                           &
' Convergence du point fixe a l''iteration ',I10               ,/)
 2603 format(                                                           &
' Non convergence du couplage vitesse pression par point fixe  ' )
 2001 format(                                                           &
'-------------------------------------------------------------',/)
 3000 format(/,                                                     &
'   ** INFORMATION SUR LE TRAITEMENT ROTOR/STATOR INSTATIONNAIRE',/,&
'      ---------------------------------------------------------',/,&
' Temps dedie a la mise a jour du maillage (s) :',F12.4,          /,&
' Temps total                              (s) :',F12.4,          /)

#else

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

#endif

!----
! End
!----

return

end subroutine
