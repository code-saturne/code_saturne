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

!> \file fldprp.f90
!> \brief Properties definition initialization, according to calculation type
!> selected by the user.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine fldprp

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
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

character(len=80) :: f_label, f_name
integer           :: ischcp
integer           :: idim1, idim3, idim6, iflid
integer           :: k_restart_id, key_n_restart_id
integer           :: type_flag, post_flag, location_id
integer           :: keypid
logical           :: has_previous

type(var_cal_opt) :: vcopt_u

!===============================================================================
! Interfaces
!===============================================================================

procedure() :: add_property_field, add_property_field_1d, hide_property
procedure() :: add_boundary_property_field_owner, ppprop

interface

  ! Interface to C function returning number of user-defined properties

  function cs_parameters_n_added_properties() result(n) &
    bind(C, name='cs_parameters_n_added_properties')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int)                                           :: n
  end function cs_parameters_n_added_properties

  ! Interface to C function building user-defined properties

  subroutine cs_parameters_create_added_properties() &
    bind(C, name='cs_parameters_create_added_properties')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_parameters_create_added_properties

  ! Interface to C function building properties

  subroutine cs_parameters_define_auxiliary_fields() &
    bind(C, name='cs_parameters_define_auxiliary_fields')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_parameters_define_auxiliary_fields

  !=============================================================================

end interface

!===============================================================================
! 0. Initialisation
!===============================================================================

call field_get_key_id("parent_field_id", keypid)

idim1 = 1
idim3 = 3
idim6 = 6

! Get key id for restart file
call field_get_key_id("restart_file", k_restart_id)
call field_get_key_id("restart_n_values", key_n_restart_id)

!===============================================================================
! 1. PROPRIETES PRINCIPALES
!===============================================================================

! --- Numerotation des proprietes presentes ici
!       Ceci depend du type de solveur branche derriere
!        (CP, Poly, Mono...)
!       Dans l'ideal, il y aurait donc plusieurs fldprp.

!       Pour le moment, on fait les hypotheses suivantes :
!         Il y a toujours, pour toutes les phases,  rho, viscl, visct
!         Il y a toujours la pression totale (sauf en compressible)
!         Lorsqu'elles sont variables, on a les proprietes suivantes :
!           . cp
!           . visls (par scalaire)
!           . csmago en LES dynamique
!         En ALE on a la viscosite de maillage
!         On a aussi les flux de masse porteurs :
!           . les variables u,v,w,p,turbulence sont portees par leur
!               phase (1 flux)
!           . les scalaires sont portes par leur phase (1 flux)
!           On suppose donc qu'il n'y a pas de scalaire en
!             taux de vide a convecter avec un flux particulier (ce
!             serait possible : ca rajoute un test, par exemple
!             if alpro...

!     ATTENTION : on s'arrange pour numeroter a la queue-leu-leu sans
!       trous les proprietes qui sont definies au centre des cellules
!       ceci permet ensuite de ne pas se fatiguer lors de la
!       construction de IPPPRO plus bas.
!      Cependant, pour les physiques particulieres, ce n'est pas le cas.

! Base properties, always present

call add_property_field_1d('density', 'Density', irom)
icrom = irom
! Postprocessed and in the log file by default, hidden in modini if not variable
call field_set_key_int(icrom, keylog, 1)
call field_set_key_int(icrom, keyvis, 1)

call add_boundary_property_field_owner('boundary_density', 'Boundary Density', &
                                       ibrom)

call add_property_field_1d('molecular_viscosity', 'Laminar Viscosity', iviscl)

call add_property_field_1d('turbulent_viscosity', 'Turb Viscosity', ivisct)
if (iturb.eq.0) then
  call hide_property(ivisct)
endif

! Hybrid RANS/LES function f_d is stored for Post Processing in hybrid_blend.
! If hybrid spatial scheme is activated for the velocity (ischcv=3)
! creation of the field hybrid_blend wihich contains the
! local blending factor for each cell

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
ischcp = vcopt_u%ischcv
if (ischcp.eq.3.or.hybrid_turb.gt.0) then
  call add_property_field_1d('hybrid_blend', 'Hybrid blending function', iflid)
end if

if (hybrid_turb.eq.3) then
  call add_property_field_1d('hybrid_sas_source_term', &
                             'SAS hybrid source term', iflid)
endif

if (hybrid_turb.eq.4) then
  call add_property_field_1d('k_tot',     'Energy total',     iflid)
  call add_property_field_1d('k_mod',     'Modelised Energy', iflid)
  call add_property_field_1d('k_res',     'Resolved Energy',  iflid)
  call add_property_field_1d('eps_mod',   'Mean Dissipation', iflid)
  if (iturb.eq.60) then
    call add_property_field_1d('omg_mod',  'Mean Specific Dissipation', iflid)
    call add_property_field_1d('f1_kwsst', 'Function F1 of k-omg SST',  iflid)
  end if
  call add_property_field_1d('htles_psi', 'Psi HTLES',          iflid)
  call add_property_field_1d('htles_r',   'Energy ratio',       iflid)
  call add_property_field_1d('htles_t',   'Time scale HTLES',   iflid)
  call add_property_field_1d('htles_icc', 'ICC coefficient',    iflid)
  call add_property_field_1d('htles_fs',  'Shielding function', iflid)
  call add_property_field_1d('Delta_max', 'Delta max',          iflid)

  ! Time averaged with exponential filtering, TODO use standard time moment
  call add_property_field_1d('vel_mag_mean','Mean velocity mag.',iflid)
  has_previous = .false.
  call add_property_field('velocity_mean', 'Vel Tavg', idim3, &
                          has_previous, iflid)
  ! Diagonal part of time moment of uiuj
  call add_property_field('ui2_mean', 'Vel Tavg', idim3, &
                          has_previous, iflid)

endif

if  (iturb.eq.60) then
  ! Square of the norm of the deviatoric part of the deformation rate
  ! tensor (\f$S^2=2S_{ij}^D S_{ij}^D\f$).
  call add_property_field_1d('s2', 'S2', iflid)
  call hide_property(iflid)

  ! Divergence of the velocity. More precisely, trace of the velocity gradient
  ! (and not a finite volume divergence term). Defined only for k-omega SST
  ! (because in this case it may be calculated at the same time as \f$S^2\f$)
  call add_property_field_1d('vel_gradient_trace', 'Vel. Gradient Trace', iflid)
  call hide_property(iflid)
endif

call add_property_field_1d('courant_number', 'CFL', icour)
if (idtvar.lt.0) then
  call hide_property(icour)
endif

if (ivofmt.gt.0) then
  call add_property_field_1d('volume_courant_number', 'CourantNbVol', iflid)
  if (idtvar.lt.0) then
    call hide_property(iflid)
  endif
endif
call add_property_field_1d('fourier_number', 'Fourier Number', ifour)
if (idtvar.lt.0) then
  call hide_property(ifour)
endif


! Total pressure is stored in property field of index iprtot
! if the compressible module is not enabled (otherwise Ptot=P*).
! For groundwater flows, this field is the pressure head (h = H - z),
! only used if the gravity is set.

if (ippmod(icompf).lt.0) then
  call add_property_field_1d('total_pressure', 'Total Pressure', iprtot)
  ! Save total pressure in auxiliary restart file
  call field_set_key_int(iprtot, k_restart_id, RESTART_AUXILIARY)
endif

! Cs^2 si on est en LES dynamique
if (iturb.eq.41) then
  call add_property_field_1d('smagorinsky_constant^2', 'Csdyn2', ismago)
else
  ismago = 0
endif

!     Numero max des proprietes ; ressert plus bas pour
!       ajouter celles relatives a la physique particuliere

! --- Modifications pour la physique particuliere

call ppprop

if (iand(ivofmt,VOF_FREE_SURFACE).ne.0) then
  idrift = 2
  has_previous = .false.
  call add_property_field('drift_velocity', 'Drift Velocity', idim3, &
                          has_previous, iflid)
endif

! --- Mesh displacement for ALE

if (iale.ge.1) then

  has_previous = .true.
  idim3 = 3
  f_name = 'mesh_displacement'
  f_label = 'Mesh displacement'
  type_flag = FIELD_PROPERTY
  post_flag = POST_ON_LOCATION
  location_id = 4 ! variables defined on vertices

  call field_create(f_name, type_flag, location_id, idim3, &
                    has_previous, iflid)
  call field_set_key_int(iflid, keyvis, post_flag)
  call field_set_key_int(iflid, keylog, 1)

  call field_set_key_str(iflid, keylbl, trim(f_label))

  call field_set_key_int(iflid, k_restart_id, RESTART_AUXILIARY)
  call field_set_key_int(iflid, key_n_restart_id, 2)

  has_previous = .false.
  idim3 = 3
  f_name = 'vtx_coord0'
  type_flag = FIELD_PROPERTY
  location_id = 4 ! variables defined on vertices

  call field_create(f_name, type_flag, location_id, idim3, &
                    has_previous, iflid)

  call hide_property(iflid)

endif

! Properties and other fields defined in C

call cs_parameters_define_auxiliary_fields

! User-defined properties

call cs_parameters_create_added_properties

! Set itemp if temperature is present as a property

if (itherm.eq.2 .and. itemp.eq.0) then
  call field_get_id_try('temperature', iflid)
  if (iflid.ge.0) itemp = iflid
endif

if (itherm.eq.4 .and. itemp.eq.0) then
  call field_get_id_try('temperature', iflid)
  if (iflid.ge.0) itemp = iflid
endif

! Map pointers

call cs_field_pointer_map_base
call cs_field_pointer_map_boundary

!====
! End
!====

return
end subroutine fldprp

!===============================================================================

!> \brief add field defining a one-dimensional property field defined on cells,
!>        with no previous time values and with default options
!
!> It is recommended not to define property names of more than 16
!> characters, to get a clear execution log (some advanced writing
!> levels take into account only the first 16 characters).
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     name          field name
!> \param[in]     label         field default label, or empty
!> \param[out]    f_id          field id
!_______________________________________________________________________________

subroutine add_property_field_1d &
 ( name, label, f_id )

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

procedure() :: add_property_field

! Arguments

character(len=*), intent(in) :: name, label
integer, intent(out)         :: f_id

! Local variables

integer  dim
logical  has_previous

!===============================================================================

has_previous = .false.
dim = 1

call add_property_field(name, label, dim, has_previous, f_id)

call field_set_key_int(f_id, keylog, 1)
call field_set_key_int(f_id, keyvis, 1)

return

end subroutine add_property_field_1d

!===============================================================================

!> \brief disable logging and postprocessing for a property field
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     f_id          field id
!_______________________________________________________________________________

subroutine hide_property &
 ( f_id )

!===============================================================================
! Module files
!===============================================================================

use entsor
use field

!===============================================================================

implicit none

! Arguments

integer, intent(in) :: f_id

! Local variables

!===============================================================================

call field_set_key_int(f_id, keyvis, 0)
call field_set_key_int(f_id, keylog, 0)

return

end subroutine hide_property

!===============================================================================

!> \brief add field defining a property field defined on cells,
!>        with default options
!
!> It is recommended not to define property names of more than 16
!> characters, to get a clear execution log (some advanced writing
!> levels take into account only the first 16 characters).
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     name          field name
!> \param[in]     label         field default label, or empty
!> \param[in]     dim           field dimension
!> \param[in]     has_previous  indicates if the field also has previous
!>                              time step values
!> \param[out]    f_id          matching field id
!_______________________________________________________________________________

subroutine add_property_field &
 ( name, label, dim, has_previous, f_id )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use post
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

character(len=*), intent(in) :: name, label
integer, intent(in)          :: dim
logical, intent(in)          :: has_previous
integer, intent(out)         :: f_id

! Local variables

integer  type_flag, location_id

character(len=len_trim(name)+1, kind=c_char) :: c_name
integer(c_int) :: c_type_flag
integer(c_int) :: c_location_id
integer(c_int) :: c_dim
logical(c_bool) :: c_has_previous
!===============================================================================

type_flag = FIELD_INTENSIVE + FIELD_PROPERTY
location_id = 1 ! variables defined on cells

! Test if the field has already been defined
call field_get_id_try(trim(name), f_id)
if (f_id .ge. 0) then
  write(nfecra,1000) trim(name)
  call csexit (1)
endif

! Create field
c_name = trim(name)//c_null_char
c_type_flag = type_flag
c_location_id = location_id
c_dim = dim

if (has_previous) then
  c_has_previous = .true.
else
  c_has_previous = .false.
endif

call cs_physical_property_define_from_field(c_name, c_type_flag, &
  c_location_id, dim, c_has_previous)

f_id = cs_physical_property_field_id_by_name(c_name)

call field_set_key_int(f_id, keyvis, 0)
call field_set_key_int(f_id, keylog, 1)

if (len(trim(label)).gt.0) then
  call field_set_key_str(f_id, keylbl, trim(label))
endif

return

!---
! Formats
!---

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ERROR:      STOP AT THE INITIAL DATA SETUP              ',/,&
'@    ======                                                  ',/,&
'@     FIELD: ', a, 'HAS ALREADY BEEN DEFINED.                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

end subroutine add_property_field

!===============================================================================

!> \brief add owner field defining a property field defined on boundary faces,
!>        with default options
!
!> It is recommended not to define property names of more than 16
!> characters, to get a clear execution log (some advanced writing
!> levels take into account only the first 16 characters).
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     name          field name
!> \param[in]     label         field default label, or empty
!> \param[out]    f_id          matching field id
!_______________________________________________________________________________

subroutine add_boundary_property_field_owner &
 ( name, label, f_id )

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

character(len=*), intent(in) :: name, label
integer, intent(out)         :: f_id

! Local variables

integer  type_flag, location_id, dim1
logical  has_previous

procedure() :: csexit

!===============================================================================

type_flag = FIELD_INTENSIVE + FIELD_PROPERTY
location_id = 3 ! variables defined on boundary faces
dim1 = 1
has_previous = .false.

! Test if the field has already been defined
call field_get_id_try(trim(name), f_id)
if (f_id .ge. 0) then
  write(nfecra,1000) trim(name)
  call csexit (1)
endif

! Create field

call field_create(name, type_flag, location_id, dim1, has_previous, f_id)

call field_set_key_int(f_id, keyvis, 0)
call field_set_key_int(f_id, keylog, 1)

if (len(trim(label)).gt.0) then
  call field_set_key_str(f_id, keylbl, trim(label))
endif

return

!---
! Formats
!---

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ERROR:      STOP AT THE INITIAL DATA SETUP              ',/,&
'@    ======                                                  ',/,&
'@     FIELD: ', a, 'HAS ALREADY BEEN DEFINED.                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

end subroutine add_boundary_property_field_owner

!===============================================================================
