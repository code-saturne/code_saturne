!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
use lagran
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use radiat
use ihmpre
use mesh
use post
use field
use cs_c_bindings
use darcy_module

!===============================================================================

implicit none

! Arguments

! Local variables

character(len=80) :: f_label, f_name, s_name, s_label
integer           :: ii, ivar, isorb, keysrb, igwfpr, keypre, ischcp
integer           :: idim1, idim3, idim6, iflid
integer           :: type_flag, post_flag, location_id
logical           :: has_previous

type(gwf_soilwater_partition) :: sorption_scal
type(var_cal_opt) :: vcopt_u

!===============================================================================
! Interfaces
!===============================================================================

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

  !=============================================================================

end interface

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
! Postprocessed and in the log file by default, hidden in modini if not variable.
call field_set_key_int(icrom, keylog, 1)
call field_set_key_int(icrom, keyvis, 1)

call add_boundary_property_field_owner('boundary_density', 'Boundary Density', &
                                       ibrom)

call add_property_field_1d('molecular_viscosity', 'Laminar Viscosity', iviscl)

call add_property_field_1d('turbulent_viscosity', 'Turb Viscosity', ivisct)
if (iturb.eq.0) then
  call hide_property(ivisct)
endif

! If hybrid spatial scheme is activated for the velocity (ischcv=3)
! creation of the field hybrid_blend wihich contains the
! local blending factor for each cell
call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
ischcp = vcopt_u%ischcv
if (ischcp.eq.3) then
  call add_property_field_1d('hybrid_blend', 'Hybrid blending function', iflid)
end if

if  (iturb.eq.60) then
  call add_property_field_1d('s2', 'S2', is2kw)
  call hide_property(is2kw)
  call add_property_field_1d('vel_gradient_trace', 'Vel. Gradient Trace', idivukw)
  call hide_property(idivukw)
  ! Hybrid RANS/LES function f_d is stored for Post Processing in hybrid_blend.
  ! If  hybrid spatial scheme is activated for the velocity (ischcv=3) f_d
  ! is used as blending factor and this field already exists
  if (iddes.eq.1.and.ischcp.ne.3) then
    call add_property_field_1d('hybrid_blend', 'Hybrid blending function', iflid)
  end if
endif

idim3 = 3
call add_property_field('grad_p', 'Presssure gradient', idim3, .false., igradp)

call add_property_field_1d('courant_number', 'CFL', icour)
call add_property_field_1d('fourier_number', 'Fourier Number', ifour)

! Total pressure is stored in property field of index iprtot
! if the compressible module is not enabled (otherwise Ptot=P*).
! For groundwater flows, this field is the pressure head (h = H - z),
! only used if the gravity is set.

if (ippmod(icompf).lt.0.and.ippmod(idarcy).lt.0) then
  call add_property_field_1d('total_pressure', 'Total Pressure', iprtot)
else if (ippmod(idarcy).ge.0.and.darcy_gravity.ge.1) then
  call add_property_field_1d('total_pressure', 'Pressure head', iprtot)
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

! --- Properties for Darcy module

if (ippmod(idarcy).eq.1) then

  has_previous = .true.
  idim1 = 1
  idim6 = 6
  f_name = 'saturation'
  f_label = 'Saturation'
  call add_property_field(f_name, f_label, idim1, has_previous, iflid)
  f_name = 'capacity'
  f_label = 'Capacity'
  call add_property_field(f_name, f_label, idim1, has_previous, iflid)
  f_name = 'permeability'
  f_label = 'Permeability'
  if (darcy_anisotropic_permeability.eq.0) then
    call add_property_field(f_name, f_label, idim1, has_previous, iflid)
  else
    call add_property_field(f_name, f_label, idim6, has_previous, iflid)
  endif
  f_name = 'soil_density'
  f_label = 'Soil density'
  call add_property_field(f_name, f_label, idim1, has_previous, iflid)

  call field_get_key_id("gwf_sorbed_concentration_id", keysrb)
  call field_get_key_id("gwf_precip_concentration_id", keypre)

  do ii = 1, nscal
    ivar = isca(ii)
    call field_get_key_struct_gwf_soilwater_partition(ivarfl(ivar), &
                                                      sorption_scal)
    call field_get_name(ivarfl(ivar), s_name)
    call field_get_name(ivarfl(ivar), s_label)

    f_name = trim(s_name)//'_kd'
    f_label = trim(s_label)//' Kd'
    call add_property_field(f_name, f_label, idim1, has_previous, &
                            sorption_scal%ikd)
    call hide_property(sorption_scal%ikd)
    f_name = trim(s_name)//'_delay'
    f_label = trim(s_label)//' delay'
    call add_property_field(f_name, f_label, idim1, has_previous, &
                            sorption_scal%idel)

    if (sorption_scal%kinetic.eq.1) then
      f_name = trim(s_name)//'_sorb_conc'
      f_label = trim(s_label)//' sorb conc'
      call add_property_field(f_name, f_label, idim1, has_previous, isorb)
      call field_set_key_int(ivarfl(ivar), keysrb, isorb)

      f_name = trim(s_name)//'_kplus'
      f_label = trim(s_label)//' kplus'
      call add_property_field(f_name, f_label, idim1, has_previous, &
                              sorption_scal%ikp)
      call hide_property(sorption_scal%ikp)
      f_name = trim(s_name)//'_kminus'
      f_label = trim(s_label)//' kminus'
      call add_property_field(f_name, f_label, idim1, has_previous, &
                              sorption_scal%ikm)
      call hide_property(sorption_scal%ikm)
    endif

    if (sorption_scal%imxsol.ge.0) then
      f_name = trim(s_name)//'_precip_conc'
      f_label = trim(s_label)//' precip conc'
      call add_property_field(f_name, f_label, idim1, has_previous, igwfpr)
      call field_set_key_int(ivarfl(ivar), keypre, igwfpr)

      f_name = trim(s_name)//'_solubility_index'
      f_label = trim(s_label)//' solubility index'
      call add_property_field(f_name, f_label, idim1, has_previous, &
                              sorption_scal%imxsol)
      call hide_property(sorption_scal%imxsol)
    endif

    call field_set_key_struct_gwf_soilwater_partition(ivarfl(ivar), &
                                                      sorption_scal)
  enddo

endif

! --- Mesh displacement for ALE

if (iale.eq.1) then

  has_previous = .true.
  idim3 = 3
  f_name = 'disale'
  f_label = 'Mesh displacement'
  type_flag = FIELD_PROPERTY
  post_flag = POST_ON_LOCATION
  location_id = 4 ! variables defined on vertices

  call field_create(f_name, type_flag, location_id, idim3, &
                    has_previous, fdiale)
  call field_set_key_int(fdiale, keyvis, post_flag)
  call field_set_key_int(fdiale, keylog, 1)

  call field_set_key_str(fdiale, keylbl, trim(f_label))

endif


! User-defined properties

call cs_parameters_create_added_properties

! Set itemp if temperature is present as a property

if (itherm.eq.2 .and. itemp.eq.0) then
  call field_get_id_try('temperature', iflid)
  if (iflid.ge.0) itemp = iflid
endif

! Map pointers

call cs_field_pointer_map_base
call cs_field_pointer_map_boundary

return

!===============================================================================
! 2. Formats
!===============================================================================

!===============================================================================
! 5. End
!===============================================================================

return
end subroutine fldprp

!===============================================================================

!> \brief add field defining a one-dimensional property field defined on cells,
!>        with no previous time values and with default options
!
!> It is recommended not to define property names of more than 16
!> characters, to get a clear execution listing (some advanced writing
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
!> characters, to get a clear execution listing (some advanced writing
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

!===============================================================================

implicit none

! Arguments

character(len=*), intent(in) :: name, label
integer, intent(in)          :: dim
logical, intent(in)          :: has_previous
integer, intent(out)         :: f_id

! Local variables

integer  type_flag, post_flag, location_id

!===============================================================================

type_flag = FIELD_INTENSIVE + FIELD_PROPERTY
post_flag = POST_ON_LOCATION
location_id = 1 ! variables defined on cells

! Test if the field has already been defined
call field_get_id_try(trim(name), f_id)
if (f_id .ge. 0) then
  write(nfecra,1000) trim(name)
  call csexit (1)
endif

! Create field

call field_create(name, type_flag, location_id, dim, has_previous, f_id)

call field_set_key_int(f_id, keyvis, post_flag)
call field_set_key_int(f_id, keylog, 1)

if (len(trim(label)).gt.0) then
  call field_set_key_str(f_id, keylbl, trim(label))
endif

return

!---
! Formats
!---

#if defined(_CS_LANG_FR)
 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ERREUR :    ARRET A L''ENTREE DES DONNEES               ',/,&
'@    ========                                                ',/,&
'@     LE CHAMP : ', a, 'EST DEJA DEFINI.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
#else
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
#endif

end subroutine add_property_field

!===============================================================================

!> \brief add owner field defining a property field defined on boundary faces,
!>        with default options
!
!> It is recommended not to define property names of more than 16
!> characters, to get a clear execution listing (some advanced writing
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

#if defined(_CS_LANG_FR)
 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ERREUR :    ARRET A L''ENTREE DES DONNEES               ',/,&
'@    ========                                                ',/,&
'@     LE CHAMP : ', a, 'EST DEJA DEFINI.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
#else
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
#endif

end subroutine add_boundary_property_field_owner

!===============================================================================
