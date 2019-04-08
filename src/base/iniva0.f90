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
use darcy_module

!===============================================================================

implicit none

! Arguments

integer          nscal

! Local variables

integer          iis   , iscal
integer          iel   , ifac
integer          iclip , ii    , jj    , idim, f_dim
integer          ifcvsl
integer          iflid, nfld, ifmaip, bfmaip, iflmas, iflmab
integer          kscmin, kscmax
integer          f_type, idftnp
integer          keyvar
integer          f_id, kdflim

logical          have_previous

double precision xxk, xcmu, trii, clvfmn

double precision, dimension(:), pointer :: dt
double precision, dimension(:), pointer :: brom, crom
double precision, dimension(:), pointer :: cofbcp
double precision, dimension(:), pointer :: porosi
double precision, dimension(:,:), pointer :: porosf
double precision, dimension(:), pointer :: field_s_v
double precision, dimension(:), pointer :: cpro_diff_lim
double precision, dimension(:), pointer :: cvar_pr
double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_al
double precision, dimension(:), pointer :: cvar_phi, cvar_fb, cvar_omg, cvar_nusa
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cvar_r12, cvar_r13, cvar_r23
double precision, dimension(:,:), pointer :: cvar_rij
double precision, dimension(:), pointer :: viscl, visct, cpro_cp, cpro_prtot
double precision, dimension(:), pointer :: cpro_viscls, cproa_viscls, cvar_tempk
double precision, dimension(:), pointer :: cpro_visma_s
double precision, dimension(:), pointer :: mix_mol_mas
double precision, dimension(:,:), pointer :: cpro_visma_v

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

! Initialize variables to avoid compiler warnings

jj = 0

! En compressible, ISYMPA initialise (= 1) car utile dans le calcul
!     du pas de temps variable avant passage dans les C.L.

if ( ippmod(icompf).ge.0 ) then
  do ifac = 1, nfabor
    isympa(ifac) = 1
  enddo
endif

!===============================================================================
! 2. PAS DE TEMPS
!===============================================================================

! dt might be used on the halo cells during the ALE initialization
! otherwise dt is synchronized in the pressure correction step.
do iel = 1, ncelet
  dt (iel) = dtref
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
do iel = 1, ncel
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

! Moleacular viscosity
call field_get_val_s(iviscl, viscl)
call field_get_val_s(ivisct, visct)

! Molecular viscosity at cells (and eventual previous value)
do iel = 1, ncel
  viscl(iel) = viscl0
enddo
call field_current_to_previous(iviscl)

! Specific heat at cells (and eventual previous value)
if(icp.ge.0) then
  call field_get_val_s(icp, cpro_cp)
  do iel = 1, ncel
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
if ((ippmod(icompf).lt.0.and.ippmod(idarcy).lt.0).or.                          &
    (ippmod(idarcy).ge.0.and.darcy_gravity.ge.1)) then
  call field_get_val_s(iprtot, cpro_prtot)
  do iel = 1, ncel
    cpro_prtot(iel) = - rinfin
  enddo
endif

! Initialization of mix_mol_mas with default values (air)
! (used in cs_cf_thermo_default_init)
if(ippmod(igmix).ge.0) then
  call field_get_val_s(igmxml, mix_mol_mas)
  do iel =1, ncel
    mix_mol_mas(iel) = xmasmr
  enddo
endif

! Default initialisations for the compressible model
if (ippmod(icompf).ge.0) then
  ! In compressible, for now, the temperature is not solved but is a field of
  ! type variable anyway. The reference value has to be taken into account.
  call field_get_val_s(ivarfl(isca(itempk)), cvar_tempk)
  do iel = 1, ncel
    cvar_tempk(iel) = t0
  enddo

  ! Default isochoric specific heat (cv0),
  ! total energy and density
  call cs_cf_thermo_default_init

  ! Default diffusivity for total energy
  visls0(ienerg) = visls0(itempk)/cv0
endif

! Diffusivite des scalaires
do iscal = 1, nscal
  call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
  ! Diffusivite aux cellules (et au pdt precedent si ordre2)
  if (ifcvsl.ge.0) then
    call field_get_val_s(ifcvsl, cpro_viscls)
    do iel = 1, ncel
      cpro_viscls(iel) = visls0(iscal)
    enddo
    call field_have_previous(ifcvsl, have_previous)
    if (have_previous) then
      call field_get_val_prev_s(ifcvsl, cproa_viscls)
      do iel = 1, ncel
        cproa_viscls(iel) = visls0(iscal)
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
    do iel = 1, ncel
      do ii = 1, 3
        cpro_visma_v(ii  ,iel) = 1.d0
        cpro_visma_v(ii+3,iel) = 0.d0
      enddo
    enddo
  else if (iand(idftnp, ISOTROPIC_DIFFUSION).ne.0) then
    call field_get_val_s(ivisma, cpro_visma_s)
    do iel = 1, ncel
      cpro_visma_s(iel) = 1.d0
    enddo
  endif

  call cs_gui_mesh_viscosity
  call usvima

endif

! Porosity
if (iporos.ge.1) then
  call field_get_val_s(ipori, porosi)
  do iel = 1, ncelet
    porosi(iel) = 1.d0
    ! No solid cells
    isolid_0(iel) = 0
  enddo

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

else
  isolid_0(1) = 0
endif

!===============================================================================
! 4. INITIALISATION STANDARD DES VARIABLES DE CALCUL
!     On complete ensuite pour les variables turbulentes et les scalaires
!===============================================================================

!     On met la pression P* a PRED0
!$omp parallel do
do iel = 1, ncel
  cvar_pr(iel) = pred0
enddo

! initialize void fraction to minimum clipping value
if (ivofmt.gt.0) then

  call field_get_key_id("min_scalar_clipping", kscmin)
  call field_get_key_double(ivarfl(ivolf2), kscmin, clvfmn)

  call field_get_val_s(ivarfl(ivolf2), field_s_v)
  do iel = 1, ncel
    field_s_v(iel) = clvfmn
  enddo

endif

!===============================================================================
! 5. INITIALISATION DE K, RIJ ET EPS
!===============================================================================

!  Si UREF n'a pas ete donnee par l'utilisateur ou a ete mal initialisee
!    (valeur negative), on met les valeurs de k, Rij, eps et omega a
!    -10*GRAND. On testera ensuite si l'utilisateur les a modifiees dans
!    usiniv ou en lisant un fichier suite.

if(itytur.eq.2 .or. itytur.eq.5) then

  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)

  xcmu = cmu
  if (iturb.eq.50) xcmu = cv2fmu
  if (iturb.eq.51) xcmu = cpalmu

  if (uref.ge.0.d0) then
    do iel = 1, ncel
      cvar_k(iel) = 1.5d0*(0.02d0*uref)**2
      cvar_ep(iel) = cvar_k(iel)**1.5d0*xcmu/almax
    enddo

    iclip = 1
    call clipke(ncelet, ncel, iclip)

  else
    do iel = 1, ncel
      cvar_k(iel) = -grand
      cvar_ep(iel) = -grand
    enddo
  endif

  if (iturb.eq.50) then
    call field_get_val_s(ivarfl(iphi), cvar_phi)
    call field_get_val_s(ivarfl(ifb), cvar_fb)
    do iel = 1, ncel
      cvar_phi(iel) = 2.d0/3.d0
      cvar_fb(iel) = 0.d0
    enddo
  endif
  if (iturb.eq.51) then
    call field_get_val_s(ivarfl(ial), cvar_al)
    call field_get_val_s(ivarfl(iphi), cvar_phi)
    do iel = 1, ncel
      cvar_phi(iel) = 2.d0/3.d0
      cvar_al(iel) = 1.d0
    enddo
  endif

elseif(itytur.eq.3) then

  call field_get_val_s(ivarfl(iep), cvar_ep)

  if (irijco.eq.1) then
    call field_get_val_v(ivarfl(irij), cvar_rij)

    if (uref.ge.0.d0) then

      trii   = (0.02d0*uref)**2

      do iel = 1, ncel
        cvar_rij(1,iel) = trii
        cvar_rij(2,iel) = trii
        cvar_rij(3,iel) = trii
        cvar_rij(4,iel) = 0.d0
        cvar_rij(5,iel) = 0.d0
        cvar_rij(6,iel) = 0.d0
        xxk = 0.5d0*(cvar_rij(1,iel)+                             &
             cvar_rij(2,iel)+cvar_rij(3,iel))
        cvar_ep(iel) = xxk**1.5d0*cmu/almax
      enddo
      iclip = 1
      call clprij2(ncelet, ncel, iclip)

    else

      do iel = 1, ncel
        cvar_rij(1,iel) = -grand
        cvar_rij(2,iel) = -grand
        cvar_rij(3,iel) = -grand
        cvar_rij(4,iel) = -grand
        cvar_rij(5,iel) = -grand
        cvar_rij(6,iel) = -grand
        cvar_ep(iel)  = -grand
      enddo

    endif
  else
    call field_get_val_s(ivarfl(ir11), cvar_r11)
    call field_get_val_s(ivarfl(ir22), cvar_r22)
    call field_get_val_s(ivarfl(ir33), cvar_r33)
    call field_get_val_s(ivarfl(ir12), cvar_r12)
    call field_get_val_s(ivarfl(ir23), cvar_r23)
    call field_get_val_s(ivarfl(ir13), cvar_r13)

    if (uref.ge.0.d0) then

      trii   = (0.02d0*uref)**2

      do iel = 1, ncel
        cvar_r11(iel) = trii
        cvar_r22(iel) = trii
        cvar_r33(iel) = trii
        cvar_r12(iel) = 0.d0
        cvar_r13(iel) = 0.d0
        cvar_r23(iel) = 0.d0
        xxk = 0.5d0*(cvar_r11(iel)+                             &
             cvar_r22(iel)+cvar_r33(iel))
        cvar_ep(iel) = xxk**1.5d0*cmu/almax
      enddo
      iclip = 1
      call clprij(ncelet, ncel, iclip)

    else

      do iel = 1, ncel
        cvar_r11(iel) = -grand
        cvar_r22(iel) = -grand
        cvar_r33(iel) = -grand
        cvar_r12(iel) = -grand
        cvar_r13(iel) = -grand
        cvar_r23(iel) = -grand
        cvar_ep(iel)  = -grand
      enddo
    endif
  endif
  if(iturb.eq.32)then
    call field_get_val_s(ivarfl(ial), cvar_al)
    do iel = 1, ncel
      cvar_al(iel) = 1.d0
    enddo
  endif

elseif(iturb.eq.60) then

  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iomg), cvar_omg)

  if (uref.ge.0.d0) then

    do iel = 1, ncel
      cvar_k(iel) = 1.5d0*(0.02d0*uref)**2
      !     on utilise la formule classique eps=k**1.5/Cmu/ALMAX et omega=eps/Cmu/k
      cvar_omg(iel) = cvar_k(iel)**0.5d0/almax
    enddo
    !     pas la peine de clipper, les valeurs sont forcement positives

  else

    do iel = 1, ncel
      cvar_k(iel) = -grand
      cvar_omg(iel) = -grand
    enddo

  endif

elseif(iturb.eq.70) then

  call field_get_val_s(ivarfl(inusa), cvar_nusa)

  if (uref.ge.0.d0) then

    do iel = 1, ncel
      cvar_nusa(iel) = sqrt(1.5d0)*(0.02d0*uref)*almax
      !     on utilise la formule classique eps=k**1.5/Cmu/ALMAX
      !     et nusa=Cmu*k**2/eps
    enddo
    !     pas la peine de clipper, les valeurs sont forcement positives

  else

    do iel = 1, ncel
      cvar_nusa(iel) = -grand
    enddo

  endif

endif

!===============================================================================
! 6.  CLIPPING DES GRANDEURS SCALAIRES (SF K-EPS VOIR CI DESSUS)
!===============================================================================

if (nscal.gt.0) then

!    Clipping des scalaires non variance
  do iis = 1, nscal
    call field_get_dim(ivarfl(isca(iis)), f_dim)
    if(iscavr(iis).eq.0.and.f_dim.eq.1) then
      iscal = iis
      call clpsca(iscal)
      !==========
    endif
  enddo

!     Clipping des variances qui sont clippees sans recours au scalaire
!        associe
  do iis = 1, nscal
    if(iscavr(iis).ne.0.and.iclvfl(iis).ne.1) then
      iscal = iis
      call clpsca(iscal)
      !==========
    endif
  enddo

!     Clipping des variances qui sont clippees avec recours au scalaire
!        associe s'il est connu
  do iis = 1, nscal
    if (iscavr(iis).le.nscal.and.iscavr(iis).ge.1.and.iclvfl(iis).eq.1) then
      iscal = iis
      call clpsca(iscal)
      !==========
    endif
  enddo

endif

!===============================================================================
! 7.  INITIALISATION DE CONDITIONS AUX LIMITES ET FLUX DE MASSE
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
  itypfb(ifac) = 0
  itrifb(ifac) = 0
enddo

! Type symetrie : on en a besoin dans le cas du calcul des gradients
!     par moindres carres etendu avec extrapolation du gradient au bord
!     La valeur 0 permet de ne pas extrapoler le gradient sur les faces.
!     Habituellement, on evite l'extrapolation sur les faces de symetries
!     pour ne pas tomber sur une indetermination et une matrice 3*3 non
!     inversible dans les configurations 2D).
do ifac = 1, nfabor
  isympa(ifac) = 0
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
  endif

enddo

!===============================================================================
! 8.  INITIALISATIONS EN ALE
!===============================================================================

if (iale.ge.1) then
  do ii = 1, nnod
    impale(ii) = 0
  enddo
endif

if (iale.ge.1) then
  do ii = 1, nnod
    do idim = 1, 3
      xyzno0(idim,ii) = xyznod(idim,ii)
    enddo
  enddo
endif

!===============================================================================
! 9 Current to previous for variables
!===============================================================================

do iflid = 0, nfld - 1
  call field_get_type(iflid, f_type)
  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    call field_current_to_previous(iflid)
  endif
enddo


! Diffusion limiter initialization
call field_get_key_id("diffusion_limiter_id", kdflim)

do f_id = 0, nfld - 1

  call field_get_type(f_id, f_type)

  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then

    call field_get_key_int(f_id, kdflim, iflid)

    if (iflid.ne.-1) then

      call field_get_val_s(iflid, cpro_diff_lim)

      do iel = 1, ncelet
        cpro_diff_lim(iel) = 1.d0
      enddo

    endif
  endif
enddo

!----
! End
!----

return
end subroutine
