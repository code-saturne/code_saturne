!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

!> \file phyvar.f90
!>
!> \brief This subroutine fills physical properties which are variable in time
!> (mainly the eddy viscosity).
!>
!> Some user subroutines are called which allows the setting of \f$ \rho \f$,
!> \f$ \mu \f$, etc.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     iterns        Navier-Stokes sub-iterations indicator:
!>                              - if strictly negative, indicate that this
!>                                function is called outside Navier-Stokes loop
!>                              - if positive, Navier-Stokes iteration number.
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________


subroutine phyvar &
 ( nvar   , nscal  ,                                              &
   iterns , dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use pointe
use albase
use parall
use period
use ppppar
use ppthch
use ppincl
use lagran
use mesh
use field
use field_operator
use cavitation
use vof
use cs_c_bindings
use cs_nz_condensation, only: nfbpcd, ifbpcd, ztpar, izzftcd, iztag1d
use cs_nz_tagmr, only: ztpar0

!===============================================================================

implicit none

! Arguments

integer          nvar, nscal, iterns
double precision dt(ncelet)

! Local variables

character(len=80) :: chaine
integer          ivar, iel, ifac, iscal, f_id
integer          ii, jj, iok, iok1, iok2, iisct, idfm, iggafm, iebdfm
integer          nn, isou, iz
integer          mbrom, ifcvsl, iscacp
integer          idftnp
integer          kturt, turb_flux_model, turb_flux_model_type

double precision vismax(nscamx), vismin(nscamx)
double precision varmn(4), varmx(4), ttke, visls_0
double precision xttkmg, xttdrb
double precision trrij,rottke
double precision, dimension(:), pointer :: field_s_v, field_s_b
double precision, dimension(:), pointer :: brom, crom
double precision, dimension(:), pointer :: cvar_k, cvar_ep
double precision, dimension(:,:), pointer :: cvar_rij
double precision, dimension(:), pointer :: sval
double precision, dimension(:,:), pointer :: visten, vistes, cpro_visma_v
double precision, dimension(:), pointer :: viscl, visct, cpro_vis
double precision, dimension(:), pointer :: cvar_voidf
double precision, dimension(:), pointer :: cpro_var, cpro_beta, cpro_visma_s
integer, dimension(:), pointer :: ifpt1d
double precision, dimension(:), pointer :: tppt1d
double precision, allocatable, dimension(:) :: ttmp

integer          ipass
data             ipass /0/
save             ipass

type(var_cal_opt) :: vcopt

!===============================================================================
! Interfaces
!===============================================================================

procedure() :: cs_physical_properties1, cs_physical_properties2
procedure() :: uiphyv, usphyv

interface

  subroutine cs_ht_convert_h_to_t_cells_solid() &
    bind(C, name='cs_ht_convert_h_to_t_cells_solid')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_ht_convert_h_to_t_cells_solid

  subroutine cs_les_mu_t_smago_dyn() &
    bind(C, name='cs_les_mu_t_smago_dyn')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_les_mu_t_smago_dyn

  subroutine cs_les_mu_t_smago_const() &
    bind(C, name='cs_les_mu_t_smago_const')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_les_mu_t_smago_const

  subroutine cs_les_mu_t_wale() &
    bind(C, name='cs_les_mu_t_wale')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_les_mu_t_wale

  subroutine cs_turbulence_ke_mu_t(phase_id) &
    bind(C, name='cs_turbulence_ke_mu_t')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: phase_id
  end subroutine cs_turbulence_ke_mu_t

  subroutine cs_turbulence_kw_mu_t(phase_id) &
    bind(C, name='cs_turbulence_kw_mu_t')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: phase_id
  end subroutine cs_turbulence_kw_mu_t

  subroutine cs_turbulence_rij_mu_t(phase_id) &
    bind(C, name='cs_turbulence_rij_mu_t')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: phase_id
  end subroutine cs_turbulence_rij_mu_t

  subroutine cs_turbulence_sa_mu_t() &
    bind(C, name='cs_turbulence_sa_mu_t')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_turbulence_sa_mu_t

  subroutine cs_turbulence_ml_mu_t() &
    bind(C, name='cs_turbulence_ml_mu_t')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_turbulence_ml_mu_t

  subroutine cs_turbulence_v2f_phi_mu_t() &
    bind(C, name='cs_turbulence_v2f_phi_mu_t')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_turbulence_v2f_phi_mu_t

  subroutine cs_turbulence_v2f_bl_v2k_mu_t() &
    bind(C, name='cs_turbulence_v2f_bl_v2k_mu_t')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_turbulence_v2f_bl_v2k_mu_t

  subroutine cs_turbulence_rij_compute_rusanov() &
    bind(C, name='cs_turbulence_rij_compute_rusanov')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_turbulence_rij_compute_rusanov

end interface

!===============================================================================
! Initializations
!===============================================================================

ipass = ipass + 1

!===============================================================================
! User settings
!===============================================================================

mbrom = 0

! Densities at boundaries are computed in cs_vof_compute_linear_rho_mu for VoF
if (ivofmt.gt.0) then
  mbrom = 1
endif

! First computation of physical properties for specific physics
! BEFORE the user
if (ippmod(iphpar).ge.1) then

  call cs_physical_properties1(mbrom)

endif

! - Interface code_saturne
!   ======================

call uiphyv()

call usphyv(nvar, nscal, mbrom, dt)

! C version

if (mbrom.eq.0 .and. nfabor.gt.0) then
  call field_get_val_s(ibrom, brom)
  brom(1) = -grand
endif

if (itherm.eq.2) then
  call cs_ht_convert_h_to_t_cells_solid
endif

call user_physical_properties()

if (mbrom.eq.0 .and. nfabor.gt.0) then
  if (brom(1) .gt. -grand) mbrom = 1
endif

! Finalization of physical properties for specific physics
! AFTER the user
if (ippmod(iphpar).ge.1) then
  call cs_physical_properties2
endif

! Boundary density based on adjacent cell value if not explicitely set.

if (mbrom.eq.0) then
  call field_get_val_s(icrom, crom)
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    brom(ifac) = crom(iel)
  enddo
endif

! Parallelism and periodicity
!
! In cs_solve_navier_stokes and visecv, we need rho in the halo
if (irangp.ge.0.or.iperio.eq.1) then
  call field_get_val_s(icrom, crom)
  call synsca(crom)
endif

! Only density may be updated in Navier Stokes loop
if (iterns.ge.1) return

!  Au premier pas de temps du calcul
!     Si on a indique que rho (visc) etait constant
!       et qu'on l'a modifie dans cs_user_physical_properties, ca ne va pas
!     On se sert de irovar (ivivar) pour ecrire et lire
!       rho (visc) dans le fichier suite

if (ntcabs.eq.ntpabs+1) then

  ! Masse volumique aux cellules et aux faces de bord
  iok1 = 0
  if (irovar.eq.0) then
    call field_get_val_s(icrom, crom)
    call field_get_val_s(ibrom, brom)
    do iel = 1, ncel
      if ( abs(crom(iel)-ro0   ).gt.epzero) then
        iok1 = 1
      endif
    enddo
    do ifac = 1, nfabor
      if ( abs(brom(ifac)-ro0   ).gt.epzero) then
        iok1 = 1
      endif
    enddo
  endif
  if (iok1.ne.0) then
    write(nfecra,9001)
  endif

  ! Viscosite moleculaire aux cellules
  iok2 = 0
  if (ivivar.eq.0) then
    call field_get_val_s(iviscl, viscl)
    do iel = 1, ncel
      if ( abs(viscl(iel)-viscl0).gt.epzero) then
        iok2 = 1
      endif
    enddo
  endif
  if (iok2.ne.0) then
    if ( ippmod(icompf) .ge. 0 ) then
      write(nfecra,9003)
    else
      write(nfecra,9002)
    endif
  endif

  if (iok1.ne.0.or.iok2.ne.0) then
    call csexit(1)
  endif

endif

!===============================================================================
! Compute the eddy viscosity
!===============================================================================

if (iturb.eq. 0) then

! Laminar
! =======

  call field_get_val_s(ivisct, visct)

  do iel = 1, ncel
    visct(iel) = 0.d0
  enddo

elseif (iturb.eq.10) then

! Mixing length model
! ===================

  call cs_turbulence_ml_mu_t

elseif (itytur.eq.2) then

! k-epsilon
! =========

  call field_get_val_s(ivisct, visct)
  call field_get_val_s(iviscl, viscl)
  call field_get_val_s(icrom, crom)
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)

  if (itytur.eq.2) then

    call cs_turbulence_ke_mu_t(-1)

  endif

elseif (itytur.eq.3) then

  call cs_turbulence_rij_mu_t(-1)

elseif (itytur.eq.4) then

  ! LES (Smagorinsky, dynamic Smagorinsky, or Wale)

  if (iturb.eq.40) then
    call cs_les_mu_t_smago_const()
  elseif (iturb.eq.41) then
    call cs_les_mu_t_smago_dyn()
  elseif (iturb.eq.42) then
    call cs_les_mu_t_wale()
  endif

elseif (itytur.eq.5) then

! v2f (phi-model and BL-v2/k)
! ===========================

  if (iturb.eq.50) then
    call cs_turbulence_v2f_phi_mu_t
  else if (iturb.eq.51) then
    call cs_turbulence_v2f_bl_v2k_mu_t
  endif

elseif (iturb.eq.60) then

! k-omega SST
! ===========

  call cs_turbulence_kw_mu_t(-1)

elseif (iturb.eq.70) then

! Spalart-Allmaras
! ================

  call cs_turbulence_sa_mu_t

endif

!===============================================================================
! Anisotropic turbulent viscosity (symmetric)
!===============================================================================
idfm = 0
iggafm = 0
iebdfm = 0

call field_get_key_id('turbulent_flux_model', kturt)

do iscal = 1, nscal
  call field_get_key_int(ivarfl(isca(iscal)), kturt, turb_flux_model)
  turb_flux_model_type = turb_flux_model / 10

  if (turb_flux_model_type.eq.3) idfm = 1
  if (turb_flux_model.eq.31) iebdfm = 1
  ! GGDH or AFM on current scalar
  ! and if DFM, GGDH on the scalar variance
  if (turb_flux_model_type.gt.0) iggafm = 1
enddo

if (idfm.eq.1 .or. itytur.eq.3 .and. idirsm.eq.1) then

  call field_get_val_v(ivsten, visten)

  if (itytur.eq.3) then
    call field_get_val_s(icrom, crom)
    call field_get_val_s(iviscl, viscl)

    call field_get_val_s(ivarfl(iep), cvar_ep)

    call field_get_val_v(ivarfl(irij), cvar_rij)

    ! EBRSM
    if (iturb.eq.32) then

      do iel = 1, ncel
        trrij = 0.5d0*(cvar_rij(1,iel)+cvar_rij(2,iel)+cvar_rij(3,iel))
        ttke  = trrij/cvar_ep(iel)
        ! Durbin scale
        xttkmg = xct*sqrt(viscl(iel)/crom(iel)/cvar_ep(iel))
        xttdrb = max(ttke,xttkmg)
        rottke  = csrij * crom(iel) * xttdrb * cell_is_active(iel)

        do isou = 1, 6
          visten(isou, iel) = rottke*cvar_rij(isou, iel)
        enddo
      enddo

      ! Other damping for EBDFM model (see F. Dehoux thesis)
      if (iebdfm.eq.1) then
        call field_get_val_v(ivstes, vistes) !FIXME one by scalar

        if (irijco.eq.1) then
          do iel = 1, ncel
            trrij = 0.5d0*(cvar_rij(1,iel)+cvar_rij(2,iel)+cvar_rij(3,iel))
            rottke  = csrij * crom(iel) * trrij / cvar_ep(iel) * cell_is_active(iel)

            do isou = 1, 6
              vistes(isou, iel) = rottke*cvar_rij(isou, iel)
            enddo
          enddo
        else
          do iel = 1, ncel
            trrij = 0.5d0*(cvar_rij(1,iel)+cvar_rij(2,iel)+cvar_rij(3,iel))
            ttke  = trrij/cvar_ep(iel)
            ! Durbin scale
            xttkmg = xct*sqrt(viscl(iel)/crom(iel)/cvar_ep(iel))
            xttdrb = max(ttke,xttkmg)
            !FIXME xttdrbt = xttdrb*sqrt((1.d0-alpha3)*PR/XRH + alpha3)
            rottke  = csrij * crom(iel) * xttdrb * cell_is_active(iel)

            do isou = 1, 6
              vistes(isou, iel) = rottke*cvar_rij(isou, iel)
            enddo
          enddo
        endif

      ! No damping with Durbing scale for the scalar
      else if (iggafm.eq.1) then
        call field_get_val_v(ivstes, vistes)

        do iel = 1, ncel
          trrij = 0.5d0*(cvar_rij(1,iel)+cvar_rij(2,iel)+cvar_rij(3,iel))
          rottke  = csrij * crom(iel) * trrij / cvar_ep(iel) * cell_is_active(iel)

          do isou = 1, 6
            vistes(isou, iel) = rottke*cvar_rij(isou, iel)
          enddo
        enddo

      endif

    ! LRR or SSG
    else
      do iel = 1, ncel
        trrij = 0.5d0*(cvar_rij(1,iel)+cvar_rij(2,iel)+cvar_rij(3,iel))
        rottke  = csrij * crom(iel) * trrij / cvar_ep(iel) * cell_is_active(iel)

        do isou = 1, 6
          visten(isou, iel) = rottke*cvar_rij(isou, iel)
        enddo
      enddo
    endif

  else

    do iel = 1, ncel
      visten(1,iel) = 0.d0
      visten(2,iel) = 0.d0
      visten(3,iel) = 0.d0
      visten(4,iel) = 0.d0
      visten(5,iel) = 0.d0
      visten(6,iel) = 0.d0
    enddo

  endif ! itytur = 3
endif

!===============================================================================
! Eddy viscosity correction for cavitating flows
!===============================================================================

if (iand(ivofmt,VOF_MERKLE_MASS_TRANSFER).ne.0.and.icvevm.eq.1) then
  if (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60 .or. iturb.eq.70) then

    call field_get_val_s(icrom, crom)
    call field_get_val_s(ivarfl(ivolf2), cvar_voidf)
    call field_get_val_s(ivisct, visct)

    call cavitation_correct_visc_turb (crom, cvar_voidf, visct)

  endif
endif

!===============================================================================
! User modification of the turbulent viscosity and symmetric tensor diffusivity
!===============================================================================

call cs_user_physical_properties_turb_viscosity_wrapper()

!==============================================================================
! Rusanov flux
!==============================================================================

call cs_turbulence_rij_compute_rusanov()

!===============================================================================
! Checking of the user values and put turbulent viscosity to 0 in
! disabled cells
!===============================================================================

! ---> Calcul des bornes des variables et impressions

! Indicateur d'erreur
iok = 0

nn = 3

! Retrieve values of some physical properties fields
call field_get_val_s(icrom, crom)
call field_get_val_s(iviscl, viscl)
if (icp.ge.0) nn  = 4
call field_get_val_s(ibrom, brom)

! Min et max sur les cellules
do ii = 1, nn
  ivar = 0
  if (ii.eq.1) then
    call field_get_val_s(icrom, sval)
    varmx(ii) = sval(1)
    varmn(ii) = sval(1)
    do iel = 2, ncel
      varmx(ii) = max(varmx(ii),sval(iel))
      varmn(ii) = min(varmn(ii),sval(iel))
    enddo
  else
    ivar = 1
    if (ii.eq.2) call field_get_val_s(iviscl, cpro_var)
    if (ii.eq.3) then
      call field_get_val_s(ivisct, cpro_var)
      ! Set turbulent viscosity to 0 in disabled cells
      do iel = 1, ncel
        cpro_var(iel) = cpro_var(iel) * cell_is_active(iel)
      enddo
    endif
    if (ii.eq.4) call field_get_val_s(icp, cpro_var)
    varmx(ii) = cpro_var(1)
    varmn(ii) = cpro_var(1)
    do iel = 2, ncel
      varmx(ii) = max(varmx(ii),cpro_var(iel))
      varmn(ii) = min(varmn(ii),cpro_var(iel))
    enddo
  endif
  if (ivar.gt.0) then
    if (irangp.ge.0) then
      call parmax (varmx(ii))
      call parmin (varmn(ii))
    endif
  endif
enddo

! Min et max sur les faces de bord (masse volumique uniquement)
ii   = 1
call field_get_val_s(ibrom, brom)
do ifac = 1, nfabor
  varmx(ii) = max(varmx(ii),brom(ifac) )
  varmn(ii) = min(varmn(ii),brom(ifac) )
enddo

if (irangp.ge.0) then
  call parmax (varmx(ii))
  call parmin (varmn(ii))
endif

! Writings
iok1 = 0
do ii = 1, nn
  if (ii.eq.1) call field_get_name(icrom, chaine)
  if (ii.eq.2) call field_get_name(iviscl, chaine)
  if (ii.eq.3) call field_get_name(ivisct, chaine)
  if (ii.eq.4) call field_get_name(icp, chaine)
  if (ipass.eq.1.or.varmn(ii).lt.0.d0) then
    if (iok1.eq.0) then
      write(nfecra,3010)
      iok1 = 1
    endif
    if ((ii.ne.3).or.(iturb.ne.0))                          &
         write(nfecra,3011)chaine(1:16),varmn(ii),varmx(ii)
  endif
enddo
if (iok1.eq.1) write(nfecra,3012)

! Verifications de valeur physique

! Masse volumique definie
ii = 1
call field_get_name(icrom, chaine)
if (varmn(ii).lt.0.d0) then
  write(nfecra,9011)chaine(1:16),varmn(ii)
  iok = iok + 1
endif

! Viscosite moleculaire definie
ii = 2
call field_get_name(iviscl, chaine)
if (varmn(ii).lt.0.d0) then
  write(nfecra,9011)chaine(1:16),varmn(ii)
  iok = iok + 1
endif

! Viscosite turbulente definie
! on ne clippe pas mu_t en modele LES dynamique, car on a fait
! un clipping sur la viscosite totale
ii = 3
call field_get_name(ivisct, chaine)
if (varmn(ii).lt.0.d0.and.iturb.ne.41) then
  write(nfecra,9012)varmn(ii)
  iok = iok + 1
endif

! Chaleur specifique definie
if (icp.ge.0) then
  ii = 4
  call field_get_name(icp, chaine)
  if (varmn(ii).lt.0.d0) then
    iisct = 0
    if (itherm.ne.0) iisct = 1
    do iscal = 1, nscal
      call field_get_key_int(ivarfl(isca(iscal)), kscacp, iscacp)
      if (iscacp.ne.0) then
        iisct = 1
      endif
    enddo
    if (iisct.eq.1) then
      write(nfecra,9011)chaine(1:16),varmn(ii)
      iok = iok + 1
    endif
  endif
endif

! ---> Calcul des bornes des scalaires et impressions

if (nscal.ge.1) then

  iok1 = 0
  do iscal = 1, nscal

    call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
    if (ifcvsl.ge.0) then
      call field_get_val_s(ifcvsl, cpro_vis)
    endif

    vismax(iscal) = -grand
    vismin(iscal) =  grand
    if (ifcvsl.ge.0) then
      do iel = 1, ncel
        vismax(iscal) = max(vismax(iscal),cpro_vis(iel))
        vismin(iscal) = min(vismin(iscal),cpro_vis(iel))
      enddo
      if (irangp.ge.0) then
        call parmax(vismax(iscal))
        call parmin(vismin(iscal))
      endif
    else
      call field_get_key_double(ivarfl(isca(iscal)), kvisl0, visls_0)
      vismax(iscal) = visls_0
      vismin(iscal) = visls_0
    endif

    ivar = isca(iscal)
    call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
    if (vcopt%iwarni.ge.1.or.ipass.eq.1.or.vismin(iscal).le.0.d0) then
      call field_get_label(ivarfl(ivar), chaine)
      if (iok1.eq.0) then
        write(nfecra,3110)
        iok1 = 1
      endif
      write(nfecra,3111) chaine(1:16),iscal,vismin(iscal),vismax(iscal)
    endif

  enddo
  if (iok1.eq.1) write(nfecra,3112)

  ! Verifications de valeur physique

  ! IOK a deja ete initialise

  do iscal = 1, nscal

    ivar = isca(iscal)

    if (vismin(iscal).lt.0.d0) then
      call field_get_label(ivarfl(ivar), chaine)
      write(nfecra,9111)chaine(1:16),iscal,vismin(iscal)
      iok = iok + 1
    endif

    if (ibeta.ge.0) then
      iok1 = 0
      call field_get_val_s(ibeta, cpro_beta)
      do iel = 1, ncel
        if (cpro_beta(iel).lt.0.d0) iok1 = 1
      enddo
      if (iok1.eq.1) write(nfecra,9013)
    endif
  enddo

endif

! ---> Calcul des bornes de viscosite de maillage en ALE

if (iale.ge.1 .and. ntcabs.eq.ntpabs+1) then

  call field_get_key_struct_var_cal_opt(ivarfl(iuma), vcopt)
  idftnp = vcopt%idften

  iok1 = 0
  if (iand(idftnp, ANISOTROPIC_LEFT_DIFFUSION).ne.0) then
    call field_get_val_v(ivisma, cpro_visma_v)
    do ii = 1, 6
      ! Min et max sur les cellules
      varmx(1) = cpro_visma_v(ii,1)
      varmn(1) = cpro_visma_v(ii,1)
      do iel = 2, ncel
        varmx(1) = max(varmx(1),cpro_visma_v(ii,iel))
        varmn(1) = min(varmn(1),cpro_visma_v(ii,iel))
      enddo
      if (irangp.ge.0) then
        call parmax (varmx(1))
        call parmin (varmn(1))
      endif

      ! Writings
      call field_get_name(ivisma, chaine)
      if (vcopt%iwarni.ge.1.or.ipass.eq.1.or.varmn(1).lt.0.d0) then
        if (iok1.eq.0) then
          write(nfecra,3210)
          iok1 = 1
        endif
        write(nfecra,3211)chaine(1:16),varmn(1),varmx(1)
      endif

      ! Verifications de valeur physique

      ! Viscosite de maillage definie
      call field_get_name(ivisma, chaine)
      if (varmn(1).le.0.d0) then
        write(nfecra,9211) varmn(1)
        iok = iok + 1
      endif

    enddo
  else if (iand(idftnp, ISOTROPIC_DIFFUSION).ne.0) then
    call field_get_val_s(ivisma, cpro_visma_s)
    ! Min et max sur les cellules
    varmx(1) = cpro_visma_s(1)
    varmn(1) = cpro_visma_s(1)
    do iel = 2, ncel
      varmx(1) = max(varmx(1),cpro_visma_s(iel))
      varmn(1) = min(varmn(1),cpro_visma_s(iel))
    enddo
    if (irangp.ge.0) then
      call parmax (varmx(1))
      call parmin (varmn(1))
    endif

    ! Writings
    call field_get_name(ivisma, chaine)
    if (vcopt%iwarni.ge.1.or.ipass.eq.1.or.varmn(1).lt.0.d0) then
      if (iok1.eq.0) then
        write(nfecra,3210)
        iok1 = 1
      endif
      write(nfecra,3211)chaine(1:16),varmn(1),varmx(1)
    endif

    ! Verifications de valeur physique

    ! Viscosite de maillage definie
    call field_get_name(ivisma, chaine)
    if (varmn(1).le.0.d0) then
      write(nfecra,9211) varmn(1)
      iok = iok + 1
    endif

  endif

  if (iok1.eq.1) write(nfecra,3212)

endif

! --->  arret eventuel

if (iok.ne.0) then
  write(nfecra,9999)iok
  call csexit (1)
endif

! Initialize boundary temperature if present and not initialized yet
!===================================================================

call field_get_id_try('boundary_temperature', f_id)
if (f_id .ge. 0) then
  call field_get_val_s(f_id, field_s_b)
  call field_get_id_try('temperature', f_id)
  if (f_id .ge. 0) then
    call field_get_val_s(f_id, field_s_v)
    do ifac = 1, nfabor
      if (field_s_b(ifac) .le. -grand) then
        iel = ifabor(ifac)
        field_s_b(ifac) = field_s_v(iel)
      endif
    enddo
  else if (itherm.eq.2) then
    call field_get_id_try('enthalpy', f_id)
    if (f_id .ge. 0) then
      call field_get_val_s(f_id, field_s_v)
      allocate(ttmp(ncelet))
      call cs_ht_convert_h_to_t_cells(field_s_v, ttmp);
      do ifac = 1, nfabor
        if (field_s_b(ifac) .le. -grand) then
          iel = ifabor(ifac)
          field_s_b(ifac) = ttmp(iel)
        endif
      enddo
      deallocate(ttmp)
    endif
  endif
  ! Last resort
  if (f_id .ge. 0) then
    do ifac = 1, nfabor
      if (field_s_b(ifac) .le. -grand) then
        iel = ifabor(ifac)
        field_s_b(ifac) = t0
      endif
    enddo
  endif
  ! For wall condensation, initialize to user-prescribed value
  if (icondb.ge.0) then
    call cs_1d_wall_thermal_get_faces(ifpt1d)
    call cs_1d_wall_thermal_get_temp(tppt1d)
    do ii = 1, nfbpcd
      ifac = ifbpcd(ii) + 1 ! C numbering
      iz  = izzftcd(ii) + 1 ! C numbering
      if (iztag1d(iz).eq.0) then
        field_s_b(ifac) = ztpar(iz)
      else if (iztag1d(iz).eq.1) then
        field_s_b(ifac) = ztpar0(iz)
      else
        do jj = 1, nfpt1d
          if (ifpt1d(jj) == ifac) then
            field_s_b(ifac) = tppt1d(jj)
          endif
        enddo
      endif
    enddo
  endif
endif

!--------
! Formats
!--------

 3010 format(                                                     &
' -----------------------------------------',                   /,&
' Property           Min. value  Max. value',                   /,&
' -----------------------------------------'                     )
 3011 format(                                                     &
 2x,    a16,      e12.4,      e12.4                              )
 3012 format(                                                     &
' -----------------------------------------',                   /)
 3110 format(                                                     &
' --- Diffusivity:',                                            /,&
' -----------------------------------------------',             /,&
' Scalar           Number  Min. value  Max. value',             /,&
' -----------------------------------------------'               )
 3111 format(                                                     &
 1x,    a16,   i7,      e12.4,      e12.4                        )
 3112 format(                                                     &
' -----------------------------------------------',             /)
 3210 format(                                                     &
' --- Mesh viscosity (ALE method)',                             /,&
' -----------------------------------------',                   /,&
' Property           Min. value  Max. value',                   /,&
' -----------------------------------------'                     )
 3211 format(                                                     &
 2x,    a16,      e12.4,      e12.4                              )
 3212 format(                                                     &
' -----------------------------------------',                   /)

 9001  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION',   /,&
'@    ========',                                                /,&
'@    INCOHERENCY BETWEEN PARAMETERS FOR THE DENSITY.',         /,&
'@',                                                            /,&
'@  The density has been declared constant',                    /,&
'@     (IROVAR=0) but its value has been modified',             /,&
'@     in cells or at boundary faces.',                         /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  Check the interface, cs_user_parameters.f90,',              /,&
'@     and cs_user_physical_properties',                        /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9002  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION',   /,&
'@    ========',                                                /,&
'@    INCOHERENCY BETWEEN PARAMETERS FOR',                      /,&
'@                                     THE MOLECULAR VISCOSITY',/,&
'@',                                                            /,&
'@  The molecular viscosity has been declared constant',        /,&
'@     (IVIVAR=0) but its value has been  modified in cells',   /,&
'@     or at boundary faces.',                                  /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  Check the interface, cs_user_parameters.f90,',              /,&
'@     and cs_user_physical_properties',                        /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9003  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION',   /,&
'@    ========',                                                /,&
'@    INCOHERENCY BETWEEN USCFPV AND USCFX2 FOR',               /,&
'@                                     THE MOLECULAR VISCOSITY',/,&
'@',                                                            /,&
'@  In the compressible module, the molecular viscosity is',    /,&
'@     constant by default (IVIVAR=0) and the value',    /,&
'@     of IVIVAR  has not been modified in uscfx2. Yet, its',   /,&
'@     value has been modified in cs_user_physical_properties.',/,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  Verify uscfx2 and cs_user_physical_properties.',            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9011  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION',   /,&
'@    ========',                                                /,&
'@    THE PHYSICAL PROPERTY ', a16  ,' HAS NOT BEEN',           /,&
'@                                          CORRECTLY DEFINED.',/,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  The physical property identified is variable and the',      /,&
'@    minimum reached is', e12.4                               ,/,&
'@  Verify that this property has been defined and',            /,&
'@    that the chosen law leads to correct values.',            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9012  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION',   /,&
'@    ========',                                                /,&
'@    THE TURBULENT VISCOSITY HAS NOT BEEN CORRECTLY DEFINED.', /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  The  minimum reached is', e12.4                            ,/,&
'@  Verify the density definition  and the turbulent viscosity',/,&
'@  modification in cs_user_physical_properties_turb_viscosity',/,&
'@  (if any).',                                                 /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9013  format( &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING : ABORT IN THE PHYSICAL QUANTITIES COMPUTATION',  /,&
'@    =========',                                               /,&
'@    INCOHERENCY BETWEEN PARAMETERS and the volumic thermal'  ,/,&
'@    expansion coefficient Beta'                              ,/,&
'@',                                                            /,&
'@  The density has been declared variable (IROVAR=1) but',     /,&
'@     the value of Beta has not been modified     ',           /,&
'@     in GUI or cs_user_physical_properties',                  /,&
'@',                                                            /,&
'@  The calculation will not be run',                           /,&
'@',                                                            /,&
'@  Check the interface or cs_user_physical_properties.',       /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@', /)


 9111  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION',   /,&
'@    ========',                                                /,&
'@    THE DIFFUSIVITY OF THE SCALAR ', a16                     ,/,&
'@       (SCALAR NUMBER', i10   ,') HAS NOT BEEN',              /,&
'@                                          CORRECTLY DEFINED.',/,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  The physical property identified is variable and the',      /,&
'@    minimum reached is', e12.4                               ,/,&
'@  Verify that this property has been defined in',             /,&
'@    cs_user_physical_properties and',                         /,&
'@    that the chosen law leads to correct values.',            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9211  format(                                                    &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION',   /,&
'@    ========',                                                /,&
'@    THE MESH VISCOSITY HAS NOT BEEN CORRECTLY DEFINED.',      /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  The  minimum reached is', e12.4,                            /,&
'@  Verify the definition of this property.',                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9999 format(                                                     &
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE PHYSICAL QUANTITIES COMPUTATION',   /,&
'@    ========',                                                /,&
'@    SOME PHYSICAL QUANTITIES HAVE INCORRECT VALUES',          /,&
'@',                                                            /,&
'@  The calculation will not be run (',i10,' errors).',         /,&
'@',                                                            /,&
'@  Refer to previous warnings for further information.',       /,&
'@  Verify cs_user_physical_properties.',                       /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

!----
! End
!----

return
end subroutine
