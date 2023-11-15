!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
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

!> \file iniva0.f90
!> \brief Computed variable initialization.
!> The time step, the indicator of wall distance computation are also
!> initialized just before reading a restart file or use the user
!> initializations.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     nscal         total number of scalars
!______________________________________________________________________________

subroutine iniva0 &
 ( nscal  )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use pointe
use entsor
use albase
use parall
use period
use ppppar
use ppthch
use ppincl
use cplsat
use field
use mesh
use cavitation
use vof
use cs_cf_bindings
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nscal

! Local variables

integer          iis   , iscal
integer          iel   , ifac
integer          ii    , jj    , idim, f_dim
integer          ifcvsl, iclvfl, kclvfl
integer          iflid, nfld, ifmaip, bfmaip, iflmas, iflmab
integer          kscmin
integer          f_type, idftnp
integer          keyvar
integer          f_id, kdflim

logical          have_previous

double precision clvfmn, visls_0, gravn2

double precision, dimension(:), pointer :: dt
double precision, dimension(:), pointer :: brom, crom, cpro_beta
double precision, dimension(:), pointer :: cofbcp
double precision, dimension(:), pointer :: porosi
double precision, dimension(:,:), pointer :: porosf
double precision, dimension(:), pointer :: field_s_v
double precision, dimension(:), pointer :: cpro_diff_lim
double precision, dimension(:), pointer :: cvar_pr
double precision, dimension(:), pointer :: viscl, visct, cpro_cp, cpro_prtot
double precision, dimension(:), pointer :: cpro_viscls, cproa_viscls, cvar_tempk
double precision, dimension(:), pointer :: cpro_visma_s
double precision, dimension(:), pointer :: mix_mol_mas
double precision, dimension(:,:), pointer :: cpro_visma_v
double precision, dimension(:,:), pointer :: xyzno0

type(var_cal_opt) :: vcopt_uma

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_gui_mesh_viscosity()  &
    bind(C, name='cs_gui_mesh_viscosity')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_mesh_viscosity

end interface

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

call field_get_key_id("variable_id", keyvar)

call field_get_val_s_by_name('dt', dt)

! Initialize temperature to reference temperature if present
call field_get_id_try('temperature', f_id)
if (f_id .ge. 0) then
  call field_get_val_s(f_id, field_s_v)
  do iel = 1, ncelet
    field_s_v(iel) = t0
  enddo
endif

! Initialize boundary temperature to "marker" if present
call field_get_id_try('boundary_temperature', f_id)
if (f_id .ge. 0) then
  call field_get_val_s(f_id, field_s_v)
  do iel = 1, nfabor
    field_s_v(iel) = -grand
  enddo
endif

call field_get_key_id("variance_clipping", kclvfl)

! Initialize variables to avoid compiler warnings

jj = 0

!===============================================================================
! 2. PAS DE TEMPS
!===============================================================================

! dt might be used on the halo cells during the ALE initialization
! otherwise dt is synchronized in the pressure correction step.
do iel = 1, ncelet
  dt(iel) = dtref
enddo

!===============================================================================
! 3.  INITIALISATION DES PROPRIETES PHYSIQUES
!===============================================================================

call field_get_val_s(ivarfl(ipr), cvar_pr)

!     Masse volumique
call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)

!     Masse volumique aux cellules (et au pdt precedent si ordre2 ou icalhy
!     ou algo. VOF)
do iel = 1, ncelet
  crom(iel)  = ro0
enddo

!     Masse volumique aux faces de bord (et au pdt precedent si ordre2)
do ifac = 1, nfabor
  brom(ifac) = ro0
enddo

! Note: for VOF or dilatable algorithms, density at twice previous time step is
! also stored and written here with "current to previous" function
call field_current_to_previous(icrom)
call field_current_to_previous(icrom)
call field_current_to_previous(ibrom)
call field_current_to_previous(ibrom)

! Boussinesq
if (ibeta.ge.0) then
  call field_get_val_s(ibeta, cpro_beta)
  do iel = 1, ncelet
    cpro_beta(iel) = -1.d0
  enddo
endif

! Molecular viscosity
call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

! Molecular viscosity at cells (and eventual previous value)
do iel = 1, ncelet
  viscl(iel) = viscl0
enddo
call field_current_to_previous(iviscl)

! Specific heat at cells (and eventual previous value)
if(icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)
  do iel = 1, ncelet
    cpro_cp(iel) = cp0
  enddo
  call field_current_to_previous(icp)
endif

! La pression totale sera initialisee a P0 + rho.g.r dans INIVAR
!  si l'utilisateur n'a pas fait d'initialisation personnelle
! Non valable en compressible
! For groundwater flows, the field of index iprtot is the
! pressure head (h = H - z). h is only used when gravity is taken
! into account

gravn2 = gx**2+gy**2+gz**2

if (ippmod(icompf).lt.0) then
  call field_get_val_s(iprtot, cpro_prtot)
  do iel = 1, ncelet
    cpro_prtot(iel) = - rinfin
  enddo
endif

! Initialization of mix_mol_mas with default values (air)
! (used in cs_cf_thermo_default_init)
if(ippmod(igmix).ge.0) then
  call field_get_val_s(igmxml, mix_mol_mas)
  do iel =1, ncelet
    mix_mol_mas(iel) = xmasmr
  enddo
endif

! Default initialisations for the compressible model
if (ippmod(icompf).ge.0) then
  ! In compressible, for now, the temperature is not solved but is a field of
  ! type variable anyway. The reference value has to be taken into account.
  call field_get_val_s(ivarfl(isca(itempk)), cvar_tempk)
  do iel = 1, ncelet
    cvar_tempk(iel) = t0
  enddo

  ! Default isochoric specific heat (cv0),
  ! total energy and density
  call cs_cf_thermo_default_init

  ! Default diffusivity for total energy

  call field_get_key_double(ivarfl(isca(itempk)), kvisl0, visls_0)
  visls_0 = visls_0 / cv0
  call field_set_key_double(ivarfl(isca(ienerg)), kvisl0, visls_0)
endif

! Diffusivite des scalaires
do iscal = 1, nscal
  call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
  call field_get_key_double(ivarfl(isca(iscal)), kvisl0, visls_0)
  ! Diffusivite aux cellules (et au pdt precedent si ordre2)
  if (ifcvsl.ge.0) then
    call field_get_val_s(ifcvsl, cpro_viscls)
    do iel = 1, ncelet
      cpro_viscls(iel) = visls_0
    enddo
    call field_have_previous(ifcvsl, have_previous)
    if (have_previous) then
      call field_get_val_prev_s(ifcvsl, cproa_viscls)
      do iel = 1, ncelet
        cproa_viscls(iel) = visls_0
      enddo
    endif
  endif
enddo

! Mesh viscosity for ALE
if (iale.ge.1) then

  call field_get_key_struct_var_cal_opt(ivarfl(iuma), vcopt_uma)
  idftnp = vcopt_uma%idften

  if (iand(idftnp, ANISOTROPIC_LEFT_DIFFUSION).ne.0) then
    call field_get_val_v(ivisma, cpro_visma_v)
    do iel = 1, ncelet
      do ii = 1, 3
        cpro_visma_v(ii  ,iel) = 1.d0
        cpro_visma_v(ii+3,iel) = 0.d0
      enddo
    enddo
  else if (iand(idftnp, ISOTROPIC_DIFFUSION).ne.0) then
    call field_get_val_s(ivisma, cpro_visma_s)
    do iel = 1, ncelet
      cpro_visma_s(iel) = 1.d0
    enddo
  endif

  call cs_gui_mesh_viscosity

endif

! Porosity
if (iporos.ge.1) then
  call field_get_val_s(ipori, porosi)
  if (compute_porosity_from_scan) then
    do iel = 1, ncelet
      porosi(iel) = 0.d0
    enddo
  else if (ibm_porosity_mode.gt.0) then
    call field_get_id_try('i_face_porosity', f_id)
    if (f_id .ge. 0) then
      call field_get_val_s(f_id, field_s_v)
      do iel = 1, nfac
        field_s_v(iel) = 1.d0
      enddo
    endif
    call field_get_id_try('b_face_porosity', f_id)
    if (f_id .ge. 0) then
      call field_get_val_s(f_id, field_s_v)
      do iel = 1, nfabor
        field_s_v(iel) = 1.d0
      enddo
    endif
  else
    do iel = 1, ncelet
      porosi(iel) = 1.d0
    enddo
  endif
  ! Tensorial porosity
  if (iporos.eq.2) then
    call field_get_val_v(iporf, porosf)
    do iel = 1, ncelet
      porosf(1, iel) = 1.d0
      porosf(2, iel) = 1.d0
      porosf(3, iel) = 1.d0
      porosf(4, iel) = 0.d0
      porosf(5, iel) = 0.d0
      porosf(6, iel) = 0.d0
    enddo
  endif
endif

!===============================================================================
! 4. INITIALISATION STANDARD DES VARIABLES DE CALCUL
!     On complete ensuite pour les variables turbulentes et les scalaires
!===============================================================================

!     On met la pression P* a PRED0
do iel = 1, ncelet
  cvar_pr(iel) = pred0
enddo

! initialize void fraction to minimum clipping value
if (ivofmt.gt.0) then

  call field_get_key_id("min_scalar_clipping", kscmin)
  call field_get_key_double(ivarfl(ivolf2), kscmin, clvfmn)

  call field_get_val_s(ivarfl(ivolf2), field_s_v)
  do iel = 1, ncelet
    field_s_v(iel) = clvfmn
  enddo

endif

call cs_turbulence_init_by_ref_quantities

!===============================================================================
! 5.  CLIPPING DES GRANDEURS SCALAIRES (SF K-EPS VOIR CI DESSUS)
!===============================================================================

if (nscal.gt.0) then

!    Clipping des scalaires non variance
  do iis = 1, nscal
    call field_get_dim(ivarfl(isca(iis)), f_dim)
    if(iscavr(iis).eq.0.and.f_dim.eq.1) then
      iscal = iis
      call clpsca(iscal)
    endif
  enddo

!     Clipping des variances qui sont clippees sans recours au scalaire
!        associe
  do iis = 1, nscal
    if (iscavr(iis).ne.0) then
      call field_get_key_int(ivarfl(isca(iis)), kclvfl, iclvfl)
      if (iclvfl.ne.1) then
        iscal = iis
        call clpsca(iscal)
      endif
    endif
  enddo

!     Clipping des variances qui sont clippees avec recours au scalaire
!        associe s'il est connu
  do iis = 1, nscal
    if (iscavr(iis).le.nscal.and.iscavr(iis).ge.1) then
      call field_get_key_int(ivarfl(isca(iis)), kclvfl, iclvfl)
      if (iclvfl.ne.1) then
        iscal = iis
        call clpsca(iscal)
      endif
    endif
  enddo

endif

!===============================================================================
! 6.  INITIALISATION DE CONDITIONS AUX LIMITES ET FLUX DE MASSE
!      NOTER QUE LES CONDITIONS AUX LIMITES PEUVENT ETRE UTILISEES DANS
!      PHYVAR, PRECLI
!===============================================================================

! Conditions aux limites

if (ienerg.gt.0) then
  call field_get_coefbc_s(ivarfl(isca(ienerg)), cofbcp)
  do ifac = 1, nfabor
    cofbcp(ifac) = 0.d0
  enddo
endif

! Boundary conditions

do ifac = 1, nfabor
  itrifb(ifac) = 0
enddo

! Old mass flux. We try not to do the same operation multiple times
! (for shared mass fluxes), without doing too complex tests.

call field_get_n_fields(nfld)

ifmaip = -1
bfmaip = -1

do iflid = 0, nfld - 1
  call field_get_type(iflid, f_type)
  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    ! Is this field not managed by CDO ? Not useful with CDO
    if (iand(f_type, FIELD_CDO)/=FIELD_CDO) then

      call field_get_key_int(iflid, kimasf, iflmas) ! interior mass flux
      call field_get_key_int(iflid, kbmasf, iflmab) ! boundary mass flux

      if (iflmas.ge.0 .and. iflmas.ne.ifmaip) then
        call field_current_to_previous(iflid)
        ifmaip = iflmas
      endif

      if (iflmab.ge.0 .and. iflmab.ne.bfmaip) then
        call field_current_to_previous(iflid)
        bfmaip = iflmab
      endif

    endif ! CDO ?
  endif ! VARIABLE ?

enddo

!===============================================================================
! 7.  INITIALISATIONS EN ALE
!===============================================================================

if (iale.ge.1) then
  call field_get_val_v_by_name("vtx_coord0", xyzno0)
  do ii = 1, nnod
    do idim = 1, 3
      xyzno0(idim,ii) = xyznod(idim,ii)
    enddo
  enddo
endif

!===============================================================================
! 8 Current to previous for variables
!===============================================================================

do iflid = 0, nfld - 1
  call field_get_type(iflid, f_type)
  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    ! Is this field not managed by CDO ? Perfomed elsewhere with CDO
    if (iand(f_type, FIELD_CDO)/=FIELD_CDO) then

      call field_current_to_previous(iflid)

    endif ! CDO ?
  endif ! VARIABLE ?
enddo


! Diffusion limiter initialization
call field_get_key_id("diffusion_limiter_id", kdflim)

do f_id = 0, nfld - 1

  call field_get_type(f_id, f_type)

  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    ! Is this field not managed by CDO ?
    if (iand(f_type, FIELD_CDO)/=FIELD_CDO) then

      call field_get_key_int(f_id, kdflim, iflid)

      if (iflid.ne.-1) then

        call field_get_val_s(iflid, cpro_diff_lim)

        do iel = 1, ncelet
          cpro_diff_lim(iel) = 1.d0
        enddo

      endif

    endif ! CDO ?
  endif ! VARIABLE ?
enddo

!----
! End
!----

return
end subroutine
