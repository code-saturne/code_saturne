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

!> \file varpos.f90
!> \brief Variables location initialization, according to calculation type
!> selected by the user.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine varpos

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use cstnum
use entsor
use albase
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use radiat
use mesh
use post
use field
use pointe, only:compute_porosity_from_scan, ibm_porosity_mode
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

integer          iscal , id, ityloc, itycat, pflag
integer          ii
integer          iok
integer          f_id, idftnp
integer          key_buoyant_id, is_buoyant_fld, key_restart_file
integer          keydri
integer          kturt, turb_flux_model, kisso2t
integer          ivar, iscdri, isso2t

double precision gravn2

character(len=80) :: f_name, s_label, s_name

type(var_cal_opt) :: vcopt

procedure() :: add_property_field, add_property_field_1d, hide_property
procedure() :: add_source_term_prev_field

!===============================================================================
! Initialization
!===============================================================================

! Key id for buoyant field (inside the Navier Stokes loop)
call field_get_key_id("is_buoyant", key_buoyant_id)
call field_get_key_id('turbulent_flux_model', kturt)
call field_get_key_id("scalar_time_scheme", kisso2t)
call field_get_key_id("restart_file", key_restart_file)

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)

! Determine itycor now that irccor is known (iturb/itytur known much earlier)
! type of rotation/curvature correction for turbulent viscosity models
if (irccor.eq.1.and.(itytur.eq.2.or.itytur.eq.5)) then
  itycor = 1
else if (irccor.eq.1.and.(iturb.eq.60.or.iturb.eq.70)) then
  itycor = 2
endif

pflag = POST_ON_LOCATION + POST_MONITOR

!===============================================================================
! Additional physical properties
!===============================================================================

if (ivofmt.gt.0) then
  ! variable density
  irovar = 1
  ivivar = 1
endif

! CP when variable
if (icp.ge.0) then
  call add_property_field_1d('specific_heat', 'Specific Heat', icp)

  call field_set_key_int(icp, keyvis, 1)
  call field_set_key_int(icp, keylog, 1)
endif

! ALE mesh viscosity
if (iale.ge.1) then
  call field_get_key_struct_var_cal_opt(ivarfl(iuma), vcopt)
  idftnp = vcopt%idften

  if (iand(idftnp, ISOTROPIC_DIFFUSION).ne.0) then
    call add_property_field('mesh_viscosity', 'Mesh Visc', 1, .false., ivisma)
  else if (iand(idftnp, ANISOTROPIC_LEFT_DIFFUSION).ne.0) then
    call add_property_field('mesh_viscosity', 'Mesh Visc', 6, .false., ivisma)
  endif
endif

if (irccor.eq.1) then
  if (idtvar.ge.0) then
    !> Strain rate tensor at the previous time step
    call add_property_field('strain_rate_tensor', 'Strain Rate Tensor', 6, &
                            .false., f_id)
    call hide_property(f_id)
  endif
endif

!===============================================================================
! Time-scheme related properties
!===============================================================================

! Dans le cas ou on a un schema en temps d'ordre 2, il faut aussi
!   prevoir les proprietes au temps n-1. Ce sera fait au dernier appel

! Dans le cas ou on calcule des moments, il faut en prevoir le nombre
!   et prevoir le nombre de tableaux necessaires pour stocker le
!   temps cumule de moyenne. On suppose que l'on peut s'aider
!   d'infos situees en tete de fichier suite (si on fait une
!   suite avec des moments non reinitialises).

! 1.1 PROPRIETES ADDITIONNELLES POUR LES ET SCHEMA EN TEMPS
! ---------------------------------------------------------

! Initialisations par defaut eventuelles et verifications
!   des options utilisees ci-dessous pour decider si l'on
!   reserve des tableaux supplementaires pour des grandeurs
!   au pas de temps precedent

iok = 0

!  Pression hydrostatique
if (iphydr.eq.1.and.icalhy.eq.-1) then
  gravn2 = gx**2+gy**2+gz**2
  if (gravn2.lt.epzero**2) then
    icalhy = 0
  else
    icalhy = 1
  endif
else
  icalhy = 0
endif

!   Global time stepping
!       for LES: 2nd order; 1st order otherwise
!       (2nd order forbidden for "coupled" k-epsilon)
if (ischtp.eq.-1) then
  if ((itytur.eq.4).or.(hybrid_turb.eq.4)) then
    ischtp = 2
  else
    ischtp = 1
  endif
endif

!   Schemas en temps : variables deduites
!   Schema pour le Flux de masse
if (istmpf.eq.-999) then
  if (ischtp.eq.1) then
    istmpf = 1
  else if (ischtp.eq.2) then
    istmpf = 2
  endif
endif

!   Collocated time scheme for gaz combustion
if (itpcol.eq.-1) then
  if (ischtp.eq.2.and.ippmod(islfm).ge.0) then
    itpcol = 1
  else
    itpcol = 0
  endif
endif

!     Termes sources NS,
if (isno2t.eq.-999) then
  if (ischtp.eq.1) then
    isno2t = 0
  else if (ischtp.eq.2) then
    !       Pour le moment par defaut on prend l'ordre 2
    isno2t = 1
  endif
endif
!     Termes sources turbulence (k-eps, Rij, v2f ou k-omega)
!     On n'autorise de changer ISTO2T qu'en Rij (sinon avec
!       le couplage k-eps/omega il y a pb)
if (isto2t.eq.-999) then
  if (ischtp.eq.1) then
    isto2t = 0
  else if (ischtp.eq.2) then
    !       Pour le moment par defaut on ne prend pas l'ordre 2
    !              ISTO2T = 1
    isto2t = 0
  endif
else if (itytur.eq.2.or.iturb.eq.50.or.iturb.ne.60) then
  write(nfecra,8132) iturb,isto2t
  iok = iok + 1
endif

do iscal = 1, nscal
  ! Termes sources Scalaires,
  call field_get_key_int(ivarfl(isca(iscal)), kisso2t, isso2t)
  if (isso2t.eq.-1) then
    if (ischtp.eq.1) then
      isso2t = 0
      call field_set_key_int(ivarfl(isca(iscal)), kisso2t, isso2t)
    else if (ischtp.eq.2) then
      ! Pour coherence avec Navier Stokes on prend l'ordre 2
      ! mais de toute facon qui dit ordre 2 dit LES et donc
      ! generalement pas de TS scalaire a interpoler.
      isso2t = 1
      call field_set_key_int(ivarfl(isca(iscal)), kisso2t, isso2t)
      if (iscal.eq.iscalt .and. iirayo.gt.0) then
        isso2t = 0
        call field_set_key_int(ivarfl(isca(iscal)), kisso2t, isso2t)
      end if
    endif
  endif

  call field_get_key_int(ivarfl(isca(iscal)), kturt, turb_flux_model)

  if (iscal.eq.iscalt) then
    if (turb_flux_model.gt.0.and.irovar.eq.1) then
      call add_property_field_1d('thermal_expansion', 'Beta', ibeta)
    endif
  endif

enddo

! Schemas en temps

!     Schema en temps global.
if (ischtp.ne. 1.and.ischtp.ne.2) then
  write(nfecra,8101) 'ISCHTP',ischtp
  iok = iok + 1
endif
if (ischtp.eq. 2.and.idtvar.ne.0) then
  write(nfecra,8111) ischtp,idtvar
  iok = iok + 1
endif
if (ischtp.eq. 2.and.itytur.eq.2) then
  write(nfecra,8112) ischtp,iturb
  iok = iok + 1
endif
if (ischtp.eq.1.and.itytur.eq.4) then
  write(nfecra,8113) ischtp,iturb
endif
if (ischtp.eq. 2.and.iturb.eq.50) then
  write(nfecra,8114) ischtp,iturb
  iok = iok + 1
endif
if (ischtp.eq. 2.and.iturb.eq.51.and.hybrid_turb.ne.4) then
  write(nfecra,8117) ischtp,iturb
  iok = iok + 1
endif
if (ischtp.eq. 2.and.iturb.eq.60.and.hybrid_turb.ne.4) then
  write(nfecra,8115) ischtp,iturb
  iok = iok + 1
endif
if (ischtp.eq. 2.and.iturb.eq.70) then
  write(nfecra,8116) ischtp,iturb
  iok = iok + 1
endif

! Schema en temps pour le flux de masse
if (istmpf.ne. 2.and.istmpf.ne.0.and.istmpf.ne. 1) then
  write(nfecra,8121) 'ISTMPF',istmpf
  iok = iok + 1
endif

! Schema en temps pour les termes sources de NS
if (isno2t.ne.0.and.isno2t.ne. 1.and.isno2t.ne.2) then
  write(nfecra,8131) 'ISNO2T',isno2t
  iok = iok + 1
endif
! Schema en temps pour les termes sources des grandeurs turbulentes
if (isto2t.ne.0.and.isto2t.ne. 1.and.isto2t.ne.2) then
  write(nfecra,8131) 'ISTO2T',isto2t
  iok = iok + 1
endif

do iscal = 1, nscal
  ! Schema en temps pour les termes sources des scalaires
  call field_get_key_int(ivarfl(isca(iscal)), kisso2t, isso2t)
  if (isso2t.ne.0.and.isso2t.ne. 1.and.isso2t.ne.2) then
    write(nfecra,8141) iscal,'ISSO2T',isso2t
    iok = iok + 1
  endif
enddo

! Stop si probleme
if (iok.gt.0) then
  call csexit(1)
endif

! add thermal expansion field for Boussinesq approximation
! if not already added
if (idilat.eq.0) then
  call field_get_id_try('thermal_expansion', ibeta)
  if (ibeta.lt.0) then
    call add_property_field_1d('thermal_expansion', 'Beta', ibeta)
  endif
endif

! Source term for weakly compressible algorithm (semi analytic scheme)
if (idilat.ge.4) then
  do iscal = 1, nscal
    id = ivarfl(isca(iscal))
    call field_get_name(id, s_name)
    call field_get_label(id, s_label)
    f_name  = trim(s_name) // '_dila_st'
    call add_property_field_1d(f_name, '', iustdy(iscal))
    id = iustdy(iscal)
    call field_set_key_int(id, keyvis, 0)
    ! Set restart file option for source terms
    call field_set_key_int(id, key_restart_file, RESTART_AUXILIARY)
  enddo
  itsrho = nscal + 1
  call add_property_field_1d('dila_st', '', iustdy(itsrho))
  id = iustdy(iscal)
  call field_set_key_int(id, keyvis, 0)
endif

! On a besoin d'un tableau pour les termes sources de Navier Stokes
!  a extrapoler. Ce tableau est NDIM dans le cas general et NDIM+1
!  si on extrapole aussi les termes sources de l equation sur le taux
!  de vide pour l'algo. VOF.
if (isno2t.gt.0) then
  call add_source_term_prev_field(ivarfl(iu))
  if (ivofmt.gt.0) then
    call add_source_term_prev_field(ivarfl(ivolf2))
  endif
endif

if (isto2t.gt.0) then
  ! The dimension of this array depends on turbulence model:
  if (itytur.eq.2) then
    call add_source_term_prev_field(ivarfl(ik))
    call add_source_term_prev_field(ivarfl(iep))
  else if (itytur.eq.3) then
    call add_source_term_prev_field(ivarfl(irij))
    call add_source_term_prev_field(ivarfl(iep))
    if (iturb.eq.32) then
      call add_source_term_prev_field(ivarfl(ial))
    endif
  else if (itytur.eq.5) then
    call add_source_term_prev_field(ivarfl(ik))
    call add_source_term_prev_field(ivarfl(iep))
    call add_source_term_prev_field(ivarfl(iphi))
    if (iturb.eq.50) then
      call add_source_term_prev_field(ivarfl(ifb))
    else if (iturb.eq.51) then
      call add_source_term_prev_field(ivarfl(ial))
    endif
  else if (iturb.eq.60) then
    call add_source_term_prev_field(ivarfl(ik))
    call add_source_term_prev_field(ivarfl(iomg))
  else if (iturb.eq.70) then
    call add_source_term_prev_field(ivarfl(inusa))
  endif
endif

! Proprietes des scalaires : termes sources pour theta schema
if (nscal.ge.1) then
  do ii = 1, nscal
    call field_get_key_int(ivarfl(isca(ii)), kisso2t, isso2t)
    if (isso2t.gt.0) then
      ! For buoyant scalars, save the current user source term
      call field_get_key_int(ivarfl(isca(ii)), key_buoyant_id, is_buoyant_fld)
      if (is_buoyant_fld.eq.1) then
        call add_source_term_field(ivarfl(isca(ii)))
      endif
      call add_source_term_prev_field(ivarfl(isca(ii)))
    endif
    ! Only usefull for Min/Max limiter
    call field_get_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt)
    if (vcopt%isstpc.eq.2) then
      call add_source_term_field(ivarfl(isca(ii)))
    endif
  enddo
endif

! Porosity
ityloc = 1 ! cells
itycat = FIELD_INTENSIVE + FIELD_PROPERTY

if (iporos.ge.1) then
  f_name = 'porosity'
  if (compute_porosity_from_scan .or. ibm_porosity_mode.gt.0) then
    !TODO move it to fldvar?
    call add_variable_field(f_name, f_name, 1, ivar)
    ipori = ivarfl(ivar)

    ! Pure convection equation (no time term)
    call field_get_key_struct_var_cal_opt(ipori, vcopt)
    vcopt%iconv = 1
    vcopt%blencv= 0.d0 ! Pure upwind
    vcopt%istat = 0
    vcopt%nswrsm = 1
    vcopt%idiff  = 0
    vcopt%idifft = 0
    vcopt%relaxv = 1.d0 ! No relaxation, even for steady algorithm.
    call field_set_key_struct_var_cal_opt(ipori, vcopt)

    ! Activate the drift for all scalars with key "drift" > 0
    iscdri = 1

    ! GNU function to return the value of iscdri
    ! with the bit value of iscdri at position
    ! 'DRIFT_SCALAR_ADD_DRIFT_FLUX' set to one
    iscdri = ibset(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)

    iscdri = ibset(iscdri, DRIFT_SCALAR_IMPOSED_MASS_FLUX)

    call field_set_key_int(ipori, keydri, iscdri)

  else
    call field_create(f_name, itycat, ityloc, 1, .false., ipori)
    call field_set_key_int(ipori, keylog, 1)
    call field_set_key_int(ipori, keyvis, pflag)
  endif

  f_name = 'cell_f_vol'
  call field_create(f_name,&
                    itycat,&
                    1,& ! location: cell
                    1,& ! dimension
                    .false.,&
                    f_id)

  if (iporos.eq.2) then
    f_name = 'tensorial_porosity'
    call field_create(f_name, itycat, ityloc, 6, .false., iporf)
  endif
  if (iporos.eq.3) then
    f_name = 'poro_div_duq'
    call field_create(f_name,&
                      itycat,&
                      1,& ! location: cell
                      3,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    f_name = 'i_poro_duq_0'
    call field_create(f_name,&
                      itycat,&
                      2,& ! location: inner faces
                      1,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    f_name = 'i_poro_duq_1'
    call field_create(f_name,&
                      itycat,&
                      2,& ! location: inner faces
                      1,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    f_name = 'b_poro_duq'
    call field_create(f_name,&
                      itycat,&
                      3,& ! location: boundary faces
                      1,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    f_name = 'i_f_face_normal'
    call field_create(f_name,&
                      itycat,&
                      2,& ! location: inner faces
                      3,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    f_name = 'i_f_face_surf'
    call field_create(f_name,&
                      itycat,&
                      2,& ! location: inner faces
                      1,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    f_name = 'b_f_face_normal'
    call field_create(f_name,&
                      itycat,&
                      3,& ! location: boundary faces
                      3,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    f_name = 'b_f_face_surf'
    call field_create(f_name,&
                      itycat,&
                      3,& ! location: boundary faces
                      1,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    f_name = 'b_f_face_cog'
    call field_create(f_name,&
                      itycat,&
                      3,& ! location: boundary faces
                      3,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    f_name = 'i_f_face_factor'
    call field_create(f_name,&
                      itycat,&
                      2,& ! location: inner faces
                      2,& ! dimension: 2 per face
                      .false.,&
                      f_id)
    call hide_property(f_id)

    f_name = 'b_f_face_factor'
    call field_create(f_name,&
                      itycat,&
                      3,& ! location: boundary faces
                      1,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    ! Interior faces weighting factor with new cell cog
    f_name = 'i_f_weight'
    call field_create(f_name,&
                      itycat,&
                      2,& ! location: inner faces
                      2,& ! dimension: 2 per face
                      .false.,&
                      f_id)

    ! Solid surface normal immersed in the cells
    f_name = 'c_w_face_normal'
    call field_create(f_name,&
                      itycat,&
                      1,& ! location: cell
                      3,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    ! Center of gravity of solid face immersed in the cells
    f_name = 'c_w_face_cog'
    call field_create(f_name,&
                      itycat,&
                      1,& ! location: cell
                      3,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    ! Solid surface of cells
    f_name = 'c_w_face_surf'
    call field_create(f_name,&
                      itycat,&
                      1,& ! location: cell
                      1,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    ! Distance between the centers of the cell and the solid face
    f_name = 'c_w_dist_inv'
    call field_create(f_name,&
                      itycat,&
                      1,& ! location: cell
                      1,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    ! Cell fluid center coordinates
    f_name = 'cell_f_cen'
    call field_create(f_name,&
                      itycat,&
                      1,& ! location: cell
                      3,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    ! Cell solid center coordinates
    f_name = 'cell_s_cen'
    call field_create(f_name,&
                      itycat,&
                      1,& ! location: cell
                      3,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    ! Porosity at internal faces
    f_name = 'i_face_porosity'
    call field_create(f_name,&
                      itycat,&
                      2,& ! location: inner faces
                      1,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

    ! Porosity at boundary faces
    f_name = 'b_face_porosity'
    call field_create(f_name,&
                      itycat,&
                      3,& ! location: boundary faces
                      1,& ! dimension
                      .false.,&
                      f_id)
    call hide_property(f_id)

  endif

endif

!===============================================================================
! Local time step and postprocessing fields
!===============================================================================

! Local time step

ityloc = 1 ! cells
itycat = FIELD_INTENSIVE

call field_create('dt', itycat, ityloc, 1, .false., id)
call field_set_key_str(id, keylbl, 'Local Time Step')
if (idtvar.gt.0) then
  if (idtvar.eq.2) then
    call field_set_key_int(id, keylog, 1)
    call field_set_key_int(id, keyvis, pflag)
  endif
endif

itycat = FIELD_INTENSIVE

! Transient velocity/pressure coupling, postprocessing field
! (variant used for computation is a tensorial field, not this one)

ncpdct = volume_zone_n_type_zones(VOLUME_ZONE_HEAD_LOSS)

if (ipucou.ne.0 .or. ncpdct.gt.0 .or. iporos.eq.2) then
  call field_create('dttens', itycat, ityloc, 6, .false., idtten)
  if (ipucou.ne.0 .or. ncpdct.gt.0) then
    call field_set_key_int(idtten, keyvis, POST_ON_LOCATION)
  endif
  call field_set_key_int(idtten, keylog, 1)
  if (iporos.eq.2) then
    call field_set_key_int(ivarfl(ipr), kwgrec, idtten)
  endif
endif

! Tensorial diffusivity

if (iporos.eq.2) then
  call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)
  vcopt%idften = ANISOTROPIC_LEFT_DIFFUSION
  call field_set_key_struct_var_cal_opt(ivarfl(iu), vcopt)
endif

! Diagonal cell tensor for the pressure solving when needed

if (ncpdct.gt.0.or.ipucou.eq.1.or.iporos.eq.2) then
  call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
  vcopt%idften = ANISOTROPIC_LEFT_DIFFUSION
  call field_set_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
endif

!===============================================================================
! Map to field pointers
!===============================================================================

call cs_field_pointer_map_base
call cs_field_pointer_map_boundary

return

!===============================================================================
! Formats
!===============================================================================

 8101 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 1 OR 2               ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verifier parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8111 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    WITH A SECOND ORDER SCHEME IN TIME: ISCHTP = ', I10      ,/,&
'@    IT IS NECESSARY TO USE A CONSTANT AND UNIFORM TIME STEP ',/,&
'@    BUT IDTVAR = ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verifier parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8112 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2ND ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    WITH K-EPSILON (ITURB = ',I10,' )'                       ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-epsilon.               ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8113 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   :      AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 1st ORDER SCHEME HAS BEEN IMPOSSED   (ISCHTP = ',I10   ,/,&
'@    FOR LES (ITURB = ',I10,' )'                              ,/,&
'@                                                            ',/,&
'@  The calculation will   be executed                        ',/,&
'@                                                            ',/,&
'@  It is recommended to verify  parameters.                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8114 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR PHI_FBAR (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-epsilon.               ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8117 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA                    ',/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR BL-V2/K  (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-epsilon.               ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8115 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR K-OMEGA   (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-omega.                 ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8116 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR SPALART   (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of Spalart-Allmaras.        ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8121 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1 OR 2            ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8131 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1 OR 2            ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8132 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  With the chosen turbulence model   , ITURB = ',I10         ,/,&
'@    the value of ISTO2T (extrapolation of the source terms  ',/,&
'@    for the turbulent variables) cannot be modified         ',/,&
'@    yet ISTO2T has been forced to ',I10                      ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8141 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITAL DATA FOR SCALARS    ',I10 ,/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1 OR 2            ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 5. End
!===============================================================================

return
end subroutine varpos

!===============================================================================

!> \brief add field defining previous source term values for a given field
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     f_id          base field id
!_______________________________________________________________________________

subroutine add_source_term_prev_field &
 ( f_id )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

! Arguments

integer, intent(in) :: f_id

! Local variables

character(len=64) :: f_name

integer :: type_flag, location_id, st_id, f_dim
logical :: has_previous

!===============================================================================

type_flag = FIELD_EXTENSIVE + FIELD_PROPERTY
location_id = 1 ! variables defined on cells
has_previous = .false.

! Define associated field

call field_get_dim(f_id, f_dim)
call field_get_name (f_id, f_name)

call field_create(trim(f_name)//'_prev_st', type_flag,               &
                  location_id, f_dim, has_previous, st_id)

call field_set_key_int(f_id, kstprv, st_id)

call hide_property(st_id)

return

end subroutine add_source_term_prev_field


!===============================================================================

!> \brief add field defining current source term values for a given field
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     f_id          base field id
!_______________________________________________________________________________

subroutine add_source_term_field &
 ( f_id )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

! Arguments

integer, intent(in) :: f_id

! Local variables

character(len=64) :: f_name

integer :: type_flag, location_id, st_id, f_dim

!===============================================================================

type_flag = FIELD_EXTENSIVE + FIELD_PROPERTY
location_id = 1 ! variables defined on cells

! Define asscociated field

call field_get_dim(f_id, f_dim)
call field_get_name (f_id, f_name)

call field_find_or_create(trim(f_name)//'_st', type_flag,       &
                            location_id, f_dim, st_id)

call field_set_key_int(f_id, kst, st_id)

return

end subroutine add_source_term_field
