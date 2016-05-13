!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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
use field
use cs_c_bindings
use darcy_module

!===============================================================================

implicit none

! Arguments

! Local variables

character(len=80) :: f_label, f_name, s_name
integer           :: ii
integer           :: ippok
integer           :: ipropp, idim1, idim3, idim6, iflid
integer           :: type_flag, location_id, ipp
logical           :: has_previous

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
! 0. INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings
ippok = 0

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
!           . cp    (par phase)
!           . visls (par scalaire)
!           . csmago (par phase) en LES dynamique
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

call add_property_field('density', 'Density', irom)
icrom = iprpfl(irom)

call add_boundary_property_field_owner('boundary_density', 'Boundary Density', &
                                       ibrom)

call add_property_field('molecular_viscosity', 'Laminar Viscosity', iviscl)

if (iturb.eq.0) then
  call add_property_field_hidden('turbulent_viscosity', 1, ivisct)
else
  call add_property_field('turbulent_viscosity', 'Turb Viscosity', ivisct)
endif

call add_property_field('courant_number', 'CFL', icour)
call add_property_field('fourier_number', 'Fourier Number', ifour)

! Total pressure is stored in property field of index iprtot
! if the compressible module is not enabled (otherwise Ptot=P*).
! For groundwater flows, this field is the pressure head (h = H - z),
! only used if the gravity is set.

if (ippmod(icompf).lt.0.and.ippmod(idarcy).lt.0) then
  call add_property_field('total_pressure', 'Total Pressure', iprtot)
else if (ippmod(idarcy).ge.0.and.darcy_gravity.ge.1) then
  call add_property_field('total_pressure', 'Pressure head', iprtot)
endif

! Cs^2 si on est en LES dynamique
if (iturb.eq.41) then
  call add_property_field('smagorinsky_constant^2', 'Csdyn2', ismago)
else
  ismago = 0
endif

!     Numero max des proprietes ; ressert plus bas pour
!       ajouter celles relatives a la physique particuliere

! --- Modifications pour la physique particuliere
!      des entiers NPROCE

!      Sauvegarde pour la physique particuliere de IPROP
!      afin d'initialiser les positions des variables d'etat
!      Attention IPROPP est le dernier numero affecte pour les proprietes.
ipropp = nproce

call ppprop
!==========

! --- Verification de NPROCE

if (nproce.gt.npromx) then
  write(nfecra,7200)nproce, npromx, nproce
  call csexit (1)
  !==========
endif

! --- Properties for Darcy module

if (ippmod(idarcy).eq.1) then

  has_previous = .true.
  idim1 = 1
  idim6 = 6
  f_name = 'saturation'
  f_label = 'Saturation'
  call add_property_field_owner(f_name, f_label, idim1, has_previous, iflid)
  f_name = 'capacity'
  f_label = 'Capacity'
  call add_property_field_owner(f_name, f_label, idim1, has_previous, iflid)
  f_name = 'permeability'
  f_label = 'Permeability'
  if (darcy_anisotropic_permeability.eq.0) then
    call add_property_field_owner(f_name, f_label, idim1, has_previous, iflid)
  else
    call add_property_field_owner(f_name, f_label, idim6, has_previous, iflid)
  endif
  do ii = 1, nscal

    if (isca(ii) .gt. 0) then

      call field_get_name(ivarfl(isca(ii)), s_name)

      f_name = trim(s_name)//'_delay'
      f_label = trim(s_name)//'_delay'
      call add_property_field_owner(f_name, f_label, idim1, has_previous, iflid)

    endif

  enddo

endif

! --- Mesh displacement for ALE

if (iale.eq.1) then

  has_previous = .true.
  idim3 = 3
  f_name = 'disale'
  f_label = 'Mesh displacement'
  type_flag = FIELD_PROPERTY
  location_id = 4 ! variables defined on vertices

  call field_create(f_name, type_flag, location_id, idim3, &
                    has_previous, fdiale)
  call field_set_key_int(fdiale, keyvis, 1)
  call field_set_key_int(fdiale, keylog, 1)

  ipp = field_post_id(fdiale)
  call field_set_key_int(fdiale, keyipp, ipp)

  call field_set_key_str(fdiale, keylbl, trim(f_label))

endif


! User-defined properties

call cs_parameters_create_added_properties

! Set itemp if temperature is present as a property

if (itherm.eq.2 .and. itemp.eq.0) then
  call field_get_id_try('temperature', iflid)
  if (iflid.ge.0) then
    do ii = 1, nproce
      if (iflid .eq. iprpfl(ii)) then
        itemp = ii
        exit
      endif
    enddo
  endif
endif

! Map pointers

call cs_field_pointer_map_base
call cs_field_pointer_map_boundary

return

!===============================================================================
! 2. Formats
!===============================================================================

#if defined(_CS_LANG_FR)

 7200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE PROPRIETES TROP GRAND                        ',/,&
'@                                                            ',/,&
'@  Le type de calcul defini                                  ',/,&
'@    correspond aux nombres de proprietes suivants           ',/,&
'@      au centre des cellules       : NPROCE = ',i10          ,/,&
'@  Le nombre de proprietes maximal prevu                     ',/,&
'@                      dans paramx.h est NPROMX = ',i10       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@  NPROMX doit valoir au moins ',i10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 7200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NUMBER OF VARIABLES TOO LARGE                          ',/,&
'@                                                            ',/,&
'@  The type of calculation defined                           ',/,&
'@    corresponds  to the following number of properties      ',/,&
'@      at the cell centers          : NPROCE = ',i10          ,/,&
'@  The maximum number of properties allowed                  ',/,&
'@                      in   paramx   is  NPROMX = ',i10       ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@  NPROMX must be at least     ',i10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!===============================================================================
! 5. End
!===============================================================================

return
end subroutine fldprp

!===============================================================================

!> \brief add field defining a hidden property field defined on cells
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     name          field name
!> \param[in]     dim           field dimension
!> \param[out]    iprop         matching field property id
!_______________________________________________________________________________

subroutine add_property_field_hidden &
 ( name, dim, iprop )

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

character(len=*), intent(in) :: name
integer, intent(in)          :: dim
integer, intent(out)         :: iprop

! Local variables

integer  id, type_flag, location_id, ii
logical  has_previous

!===============================================================================

type_flag = FIELD_INTENSIVE + FIELD_PROPERTY
location_id = 1 ! variables defined on cells
has_previous = .false.

! Test if the field has already been defined
call field_get_id_try(trim(name), id)
if (id .ge. 0) then
  write(nfecra,1000) trim(name)
  call csexit (1)
endif

! Create field

call field_create(name, type_flag, location_id, dim, has_previous, id)

call field_set_key_int(id, keyvis, 0)
call field_set_key_int(id, keylog, 0)

! Property number and mapping to field

iprop = nproce + 1
nproce = nproce + dim

call fldprp_check_nproce

do ii = 1, dim
  iprpfl(iprop + ii -1) = id
  ipproc(iprop + ii - 1) = iprop + ii - 1
enddo

! Postprocessing slots

do ii = 1, dim
  ipppro(iprop+ii-1) = 1
enddo

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

end subroutine add_property_field_hidden

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
!> \param[out]    iprop         matching field property id
!_______________________________________________________________________________

subroutine add_property_field &
 ( name, label, iprop )

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
integer, intent(out)         :: iprop

! Local variables

integer  type_flag, location_id, f_id, dim
logical  has_previous

!===============================================================================

type_flag = FIELD_INTENSIVE + FIELD_PROPERTY
location_id = 1 ! variables defined on cells
has_previous = .false.
dim = 1

! Test if the field has already been defined
call field_get_id_try(trim(name), f_id)
if (f_id .ge. 0) then
  write(nfecra,1000) trim(name)
  call csexit (1)
endif

! Create field

call field_create(name, type_flag, location_id, dim, has_previous, f_id)

call field_set_key_int(f_id, keyvis, 1)
call field_set_key_int(f_id, keylog, 1)

if (len(trim(label)).gt.0) then
  call field_set_key_str(f_id, keylbl, trim(label))
endif

! Property number and mapping to field

iprop = nproce + 1
nproce = nproce + 1

call fldprp_check_nproce

iprpfl(iprop) = f_id
ipproc(iprop) = iprop

! Postprocessing slots

ipppro(iprop) = field_post_id(f_id)

call field_set_key_int(f_id, keyipp, ipppro(iprop))

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

return

end subroutine add_property_field

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

subroutine hide_property_field &
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

integer  ipp, j

integer :: f_dim

!===============================================================================

call field_set_key_int(f_id, keyvis, 0)
call field_set_key_int(f_id, keylog, 0)

ipp = field_post_id(f_id)

if (ipp .gt. 1) then
  call field_get_dim(f_id, f_dim)
  do j = 1, f_dim
    ihisvr(ipp+j-1,1) = 0
  enddo
endif

return

end subroutine hide_property_field

!===============================================================================

!> \brief disable logging and postprocessing for a property
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iprop         property id
!_______________________________________________________________________________

subroutine hide_property &
 ( iprop )

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

integer, intent(in) :: iprop

! Local variables

integer  f_id

!===============================================================================

f_id = iprpfl(iprop)
call hide_property_field(f_id)

end subroutine hide_property

!===============================================================================

!> \brief check npromx is sufficient for the required number of properties.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine fldprp_check_nproce

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar

!===============================================================================

implicit none

! Arguments

! Local variables

if (nproce .gt. npromx) then
  write(nfecra,1000) nproce, npromx
  call csexit (1)
endif

return

!---
! Formats
!---

#if defined(_CS_LANG_FR)

 1000 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ERREUR :    ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    ========'                                                ,/,&
'@     NOMBRE DE PROPRIETES TROP GRAND'                        ,/,&
'@'                                                            ,/,&
'@  Le type de calcul defini'                                  ,/,&
'@    correspond a un nombre de proprietes NPROCE >= ', i10    ,/,&
'@  Le nombre de proprietes maximal prevu'                     ,/,&
'@                      dans paramx    est NPROMX  = ', i10    ,/,&
'@'                                                            ,/,&
'@  Le calcul ne sera pas execute.'                            ,/,&
'@'                                                            ,/,&
'@  Verifier les parametres'                                   ,/,&
'@'                                                            ,/,&
'@  Si NPROMX est augmente, le code doit etre reinstalle.'     ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#else

 1000 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ERROR:      STOP AT THE INITIAL DATA SETUP'              ,/,&
'@    ======'                                                  ,/,&
'@     NUMBER OF PROPERTIES TOO LARGE'                         ,/,&
'@'                                                            ,/,&
'@  The type of calculation defined'                           ,/,&
'@    corresponds to a number of properties NPROCE >= ', i10   ,/,&
'@  The maximum number of properties allowed'                  ,/,&
'@                      in   paramx     is  NPROMX  = ', i10   ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be executed'                        ,/,&
'@'                                                            ,/,&
'@  Verify   parameters.'                                      ,/,&
'@'                                                            ,/,&
'@  If NVARMX is increased, the code must be reinstalled.'     ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#endif

end subroutine fldprp_check_nproce

!===============================================================================

!> \brief add owner field defining a property field defined on cells,
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

subroutine add_property_field_owner &
 ( name, label, dim, has_previous, f_id )

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
integer, intent(in)          :: dim
logical, intent(in)          :: has_previous
integer, intent(out)         :: f_id

! Local variables

integer  type_flag, location_id, ipp

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

call field_create(name, type_flag, location_id, dim, has_previous, f_id)

call field_set_key_int(f_id, keyvis, 1)
call field_set_key_int(f_id, keylog, 1)

ipp = field_post_id(f_id)
call field_set_key_int(f_id, keyipp, ipp)

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

end subroutine add_property_field_owner

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
