!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

!> \file fldvar.f90
!> \brief Variables definition initialization, according to calculation type
!> selected by the user.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[out]    nmodpp        number of activated paricle physic models
!______________________________________________________________________________

subroutine fldvar &
( nmodpp )

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
use lagpar
use lagdim
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

!===============================================================================

implicit none

! Arguments

integer       nmodpp

! Local variables

integer       ipp   , ii
integer       iok   , ippok,  keycpl

!===============================================================================
! Interfaces
!===============================================================================

interface

  ! Interface to C function returning number of user-defined variables

  function cs_parameters_n_added_variables() result(n) &
    bind(C, name='cs_parameters_n_added_variables')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int)                                           :: n
  end function cs_parameters_n_added_variables

end interface

!===============================================================================
! 0. INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings
ippok = 0

call field_get_key_id('coupled', keycpl)

!===============================================================================
! CALCUL DE NSCAPP
! VERIFICATION DU NOMBRE DE SCALAIRES
! CONSTRUCTION DE ISCAPP
! CALCUL DE NSCAL

!  A la sortie de cette section, NSCAL, NSCAUS et NSCAPP sont connus.
!  On renseignera egalement ici les valeurs de ivisls
!    pour les scalaires physiques particulieres en question.
!  On en profite aussi pour remplir ITYTUR et ITYTURT puisque ITURB et ITURT
!    viennent d'etre definis.
!  On remplit aussi itycor puisque irccor, iturb et itytur viennent d'etre
!    definis.
!===============================================================================

! ---> Remplissage de ITYTUR
itytur = iturb/10

! ---> Remplissage de itycor :
! type de correction rotation/courbure pour les modeles de viscosite turbulente
if (irccor.eq.1.and.(itytur.eq.2.or.itytur.eq.5)) then
  itycor = 1
else if (irccor.eq.1.and.(iturb.eq.60.or.iturb.eq.70)) then
  itycor = 2
endif

! ---> Coherence modele
!     Rq : ATTENTION il faudrait renforcer le blindage

iok   = 0
nmodpp = 0
do ipp = 2, nmodmx
  if (ippmod(ipp).ne.-1) then
    nmodpp = nmodpp+1
    ippok = ipp
  endif
enddo
if (nmodpp.gt.1) then
  write(nfecra,6000)
  iok = iok + 1
endif

if (nmodpp.eq.1) then
  if (ippmod(ippok).lt.0 .or. ippmod(ippok).gt.5) then
    write(nfecra,6001)
    iok = iok + 1
  endif
endif

if (iok.ne.0) then
  call csexit (1)
  !==========
endif

! ---> On positionne l'indicateur global IPPMOD(IPHPAR)
!         0 : pas de physique particuliere
!         1 : physique particuliere enclenchee
!         2 : physique particuliere avec definition du coefficient
!             d'absorption par fichier parametrique pour le rayonnement
ippmod(iphpar) = 0
if (nmodpp.gt.0) then
  ippmod(iphpar) = 1
  if (ippmod(icompf).eq.-1 .and. ippmod(iatmos).eq.-1           &
                           .and. ippmod(iaeros).eq.-1)          &
    ippmod(iphpar) = 2
endif

! Define main variables
!======================

nvar = 0

! Velocity

call add_variable_field('velocity', 'Velocity', 3, iu)
call field_set_key_int(ivarfl(iu), keycpl, 1)

! All components point to same field
iv = iu + 1
iw = iv + 1

! Pressure

call add_variable_field('pressure', 'Pressure', 1, ipr)

if (ippmod(icompf).ge.0) then
  istat(ipr) = 1
else
  istat (ipr) = 0
endif
iconv (ipr) = 0

! Void fraction (cavitating flows)

if (icavit.ge.0) then
  call add_variable_field('void_fraction', 'Void Fraction', 1, ivoidf)
  idiff(ivoidf) = 0
endif

! --- Turbulence

if (itytur.eq.2) then
  call add_variable_field('k', 'Turb Kinetic Energy', 1, ik)
  call add_variable_field('epsilon', 'Turb Dissipation', 1, iep)
else if (itytur.eq.3) then
  call add_variable_field('r11', 'R11', 1, ir11)
  call add_variable_field('r22', 'R22', 1, ir22)
  call add_variable_field('r33', 'R33', 1, ir33)
  call add_variable_field('r12', 'R12', 1, ir12)
  call add_variable_field('r13', 'R13', 1, ir13)
  call add_variable_field('r23', 'R23', 1, ir23)
  call add_variable_field('epsilon', 'Turb Dissipation', 1, iep)
  if (iturb.eq.32) then
    call add_variable_field('alpha', 'Alphap', 1, ial)
    istat(ial)  = 0
    iconv(ial)  = 0
    ! For alpha, we always have a diagonal term, so do not shift the diagonal
    idircl(ial) = 0
  endif
else if (itytur.eq.5) then
  call add_variable_field('k', 'Turb Kinetic Energy', 1, ik)
  call add_variable_field('epsilon', 'Turb Dissipation', 1, iep)
  call add_variable_field('phi', 'Phi', 1, iphi)
  if (iturb.eq.50) then
    call add_variable_field('f_bar', 'f_bar', 1, ifb)
    istat(ifb)  = 0
    iconv(ifb)  = 0
    ! For fb, we always have a diagonal term, so do not shift the diagonal
    idircl(ifb) = 0
  else if (iturb.eq.51) then
    call add_variable_field('alpha', 'Alpha', 1, ial)
    istat(ial)  = 0
    iconv(ial)  = 0
    ! For alpha, we always have a diagonal term, so do not shift the diagonal
    idircl(ial) = 0
  endif
else if (iturb.eq.60) then
  call add_variable_field('k', 'Turb Kinetic Energy', 1, ik)
  call add_variable_field('omega', 'Omega', 1, iomg)
else if (iturb.eq.70) then
  call add_variable_field('nu_tilda', 'NuTilda', 1, inusa)
endif

! Mesh velocity with ALE

if (iale.eq.1) then

  call add_variable_field('mesh_velocity', 'Mesh Velocity', 3, iuma)
  call field_set_key_int(ivarfl(iuma), keycpl, 1)

  istat(iuma) = 0
  iconv(iuma) = 0
  imgr (iuma) = 1

  ! All components point to same field
  ivma = iuma + 1
  iwma = ivma + 1
  istat(ivma) = istat(iuma)
  iconv(ivma) = iconv(iuma)
  imgr (ivma) = imgr (iuma)
  istat(iwma) = istat(iuma)
  iconv(iwma) = iconv(iuma)
  imgr (iwma) = imgr (iuma)

endif

! Number of user variables

nscaus = cs_parameters_n_added_variables()

! ---> Lecture donnees thermochimie

call pplecd
!==========

! ---> Definition des variables

call ppvarp
!==========

! Thermal model with no specific physics

if (nmodpp.eq.0) then

  if (itherm .eq. 1) then
    call add_model_scalar_field('temperature', 'Temperature', iscalt)
  else if (itherm .eq. 2) then
    call add_model_scalar_field('enthalpy', 'Enthalpy', ihm)
    iscalt = ihm
  endif

  if (itherm.ne.0 .and. iihmpr.eq.1) then
    call uithsc(iscalt)
  endif

endif

! Initialise ivisls for specific physics fields other than variances

if (nscapp.gt.0) then
  do ii = 1, nscapp
    if (iscavr(iscapp(ii)).le.0 .and. ivisls(iscapp(ii)).lt.0) then
      ivisls(iscapp(ii)) = 0
    endif
  enddo
endif

call add_user_scalar_fields
!==========================

! ---> Verifications

iok = 0

if (nscaus.lt.0) then
  write(nfecra,6010) nscaus
  iok = iok + 1
endif

if (nscaus.gt.0 .or. nscapp.gt.0) then
  if ((nscaus+nscapp).gt.nscamx) then
    if (nscapp.le.0) then
      write(nfecra,6011) nscaus, nscamx, nscamx, nscaus
    else
      write(nfecra,6012) nscaus,nscapp,nscamx,nscamx-nscapp,nscaus+nscapp
    endif
    iok = iok + 1
  endif
endif

if (iok.ne.0) then
  call csexit (1)
  !==========
endif

return

!===============================================================================
! 5. FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)
 6000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     PLUSIEURS MODELES PHYSIQUES PARTICULIERES ACTIVES      ',/,&
'@                                                            ',/,&
'@  Un seul modele physique particuliere peut etre active a la',/,&
'@    fois.                                                   ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier les indicateurs de IPPMOD dans usppmo.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     SELECTION INCORRECTE DU MODELE PHYSIQUE PARTICULIERE   ',/,&
'@                                                            ',/,&
'@  Les valeurs des indicateurs du tableau IPPMOD ne sont pas ',/,&
'@    admissibles                                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier les indicateurs de IPPMOD dans usppmo.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE SCALAIRES ERRONE                             ',/,&
'@                                                            ',/,&
'@  Le nombre de scalaires utilisateur doit etre un entier    ',/,&
'@    positif ou nul. Il vaut ici   NSCAUS  = ',I10            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6011 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE SCALAIRES TROP GRAND                         ',/,&
'@                                                            ',/,&
'@  Le nombre de scalaires utilisateurs                       ',/,&
'@    demande                          est  NSCAUS = ',I10     ,/,&
'@  Le nombre de scalaires total                              ',/,&
'@    autorise   dans paramx           est  NSCAMX = ',I10     ,/,&
'@                                                            ',/,&
'@  La valeur maximale autorisee de NSCAUS                    ',/,&
'@                          est donc  NSCAMX        = ',I10    ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier NSCAUS.                                          ',/,&
'@                                                            ',/,&
'@  NSCAMX doit valoir au moins ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6012 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE SCALAIRES TROP GRAND                         ',/,&
'@                                                            ',/,&
'@  Le nombre de scalaires utilisateurs                       ',/,&
'@    demande                          est  NSCAUS = ',I10     ,/,&
'@  Le nombre de scalaires pour les physiques particulieres   ',/,&
'@    necessaire avec le modele choisi est  NSCAPP = ',I10     ,/,&
'@  Le nombre de scalaires total                              ',/,&
'@    autorise   dans paramx           est  NSCAMX = ',I10     ,/,&
'@                                                            ',/,&
'@  La valeur maximale autorisee de NSCAUS                    ',/,&
'@    avec le modele choisi est donc NSCAMX-NSCAPP = ',I10     ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier NSCAUS.                                          ',/,&
'@                                                            ',/,&
'@  NSCAMX doit valoir au moins ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 6000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     TOO MANY SPECIFIC PHYSICS MODULES ACTIVATED            ',/,&
'@                                                            ',/,&
'@  Only one specific physics module can be active for one    ',/,&
'@    given calculation.                                      ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Modify the indices of       IPPMOD in   usppmo.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     WRONG SELLECTION OF THE MODEL FOR SPECIFIC PHYSICS     ',/,&
'@                                                            ',/,&
'@  The values of the indices of the array IPPMOD are not     ',/,&
'@    admissible                                              ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Modify the indices of       IPPMOD in   usppmo.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     ERRONEOUS NUMBER OF SCALARS                            ',/,&
'@                                                            ',/,&
'@  The number of users scalars must be an integer either     ',/,&
'@   positive or zero. Here is      NSCAUS  = ',I10            ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6011 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NUMBER OF SCALARS TOO LARGE                            ',/,&
'@                                                            ',/,&
'@  The number of users scalars                               ',/,&
'@  requested                          is   NSCAUS = ',I10     ,/,&
'@  The total number of scalars                               ',/,&
'@    allowed    in   paramx           is   NSCAMX = ',I10     ,/,&
'@                                                            ',/,&
'@  The maximmum value allowed of   NSCAUS                    ',/,&
'@                          is in   NSCAMX        = ',I10      ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   NSCAUS.                                          ',/,&
'@                                                            ',/,&
'@  NSCAMX must be at least     ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6012 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NUMBER OF SCALARS TOO LARGE                            ',/,&
'@                                                            ',/,&
'@  The number of users scalars                               ',/,&
'@     requested                       is   NSCAUS = ',I10     ,/,&
'@  The number of scalars necessary for the specific physics'  ,/,&
'@    with the chosen model is              NSCAPP = ',I10     ,/,&
'@  The total number of scalars                               ',/,&
'@    allowed    in   paramx.h         is   NSCAMX = ',I10     ,/,&
'@                                                            ',/,&
'@  The maximum value allowed for  NSCAUS                     ',/,&
'@    with the chosen model is       NSCAMX-NSCAPP = ',I10     ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   NSCAUS.                                          ',/,&
'@                                                            ',/,&
'@  NSCAMX must be at least     ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!===============================================================================
! 5. FIN
!===============================================================================

return
end subroutine

!===============================================================================
! Local functions
!===============================================================================

!===============================================================================

!> \fn add_variable_field
!
!> \brief add field defining a general solved variable, with default options
!
!> It is recommended not to define variable names of more than 16
!> characters, to get a clear execution listing (some advanced writing
!> levels take into account only the first 16 characters).

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  name           field name
!> \param[in]  label          field default label, or empty
!> \param[in]  dim            field dimension
!> \param[out] ivar           variable number for defined field
!_______________________________________________________________________________

subroutine add_variable_field &
 ( name, label, dim, ivar )

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
integer, intent(out)         :: ivar

! Local variables

integer  id, ii, ipp
integer  type_flag, location_id
logical  interleaved, has_previous

integer, save :: keycpl = -1
integer, save :: keyvar = -1

type_flag = FIELD_INTENSIVE + FIELD_VARIABLE
location_id = 1         ! variables defined on cells
interleaved = .true.
has_previous = .true.

! Test if the field has already been defined
call field_get_id_try(trim(name), id)
if (id .ge. 0) then
  write(nfecra,1000) trim(name)
  call csexit (1)
endif

! Create field

if (keyvar.lt.0) then
  call field_get_key_id('coupled', keycpl)
  call field_get_key_id("variable_id", keyvar)
endif

call field_create(name, type_flag, location_id, dim, interleaved, has_previous, &
                  id)

call field_set_key_int(id, keyvis, 1)
call field_set_key_int(id, keylog, 1)

if (len(trim(label)).gt.0) then
  call field_set_key_str(id, keylbl, trim(label))
endif

ivar = nvar + 1
nvar = nvar + dim

! Check we have enough slots
call fldvar_check_nvar

ivarfl(ivar) = id
ipp = field_post_id(id)

call field_set_key_int(id, keyvar, ivar)

if (dim .gt. 1) then
  call field_set_key_int(id, keycpl, 1)
  do ii = 2, dim
    ivarfl(ivar + ii - 1) = id
  enddo
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

end subroutine add_variable_field

!===============================================================================

!> \function add_user_scalar_fields
!
!> \brief add fields defining user solved scalar variables,
!>        with default options
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine add_user_scalar_fields

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

! Local variables

integer  iscal, nfld1, nfld2
integer  dim, id, ipp
logical  interleaved

integer :: keyvar, keysca

!===============================================================================
! Interfaces
!===============================================================================

interface

  ! Interface to C function building user-defined variables

  subroutine cs_parameters_create_added_variables() &
    bind(C, name='cs_parameters_create_added_variables')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_parameters_create_added_variables

end interface

!===============================================================================

! Create fields

call field_get_n_fields(nfld1)

call cs_parameters_create_added_variables

call field_get_n_fields(nfld2)

! Now map those fields

iscal = 0

call field_get_key_id("scalar_id", keysca)
call field_get_key_id("variable_id", keyvar)

do id = nfld1, nfld2 - 1

  call field_get_dim(id, dim, interleaved)

  if (dim.ne.1) cycle ! fields of dimension > 1 may not be handled as scalars

  iscal = iscal + 1

  nvar = nvar + 1
  nscal = nscal + 1

  ! Check we have enough slots
  call fldvar_check_nvar

  isca(iscal) = nvar
  ivarfl(nvar) = id
  ipp = field_post_id(id)

  call field_set_key_int(id, keyvar, nvar)
  call field_set_key_int(id, keysca, iscal)

enddo

return

!---
! Formats
!---

end subroutine add_user_scalar_fields

!===============================================================================

!> \function add_model_scalar_field
!
!> \brief add field defining a non-user solved scalar variable,
!>        with default options
!
!> It is recommended not to define variable names of more than 16
!> characters, to get a clear execution listing (some advanced writing
!> levels take into account only the first 16 characters).
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  name           field name
!> \param[in]  label          field default label, or empty
!> \param[out] iscal          variable number for defined field
!_______________________________________________________________________________

subroutine add_model_scalar_field &
 ( name, label, iscal )

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
integer, intent(out)         :: iscal

! Local variables

integer  dim, id, ipp
integer  type_flag, location_id
logical  interleaved, has_previous

integer, save :: keyvar = -1
integer, save :: keysca = -1

type_flag = FIELD_INTENSIVE + FIELD_VARIABLE
dim = 1
location_id = 1 ! variables defined on cells
interleaved = .true.
has_previous = .true.

! Test if the field has already been defined
call field_get_id_try(trim(name), id)
if (id .ge. 0) then
  write(nfecra,1000) trim(name)
  call csexit (1)
endif

! Create field

if (keysca.lt.0) then
  call field_get_key_id("scalar_id", keysca)
  call field_get_key_id("variable_id", keyvar)
endif

call field_create(name, type_flag, location_id, dim, interleaved, has_previous, &
                  id)

call field_set_key_int(id, keyvis, 1)
call field_set_key_int(id, keylog, 1)

if (len(trim(label)).gt.0) then
  call field_set_key_str(id, keylbl, trim(label))
endif

nvar = nvar + 1
nscal = nscal + 1
nscapp = nscapp + 1
iscal = nscaus + nscapp

! Check we have enough slots
call fldvar_check_nvar
call fldvar_check_nscapp

isca(iscal) = nvar
iscapp(nscapp) = iscal
ivarfl(isca(iscal)) = id
ipp = field_post_id(id)

call field_set_key_int(id, keyvar, nvar)
call field_set_key_int(id, keysca, iscal)

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

end subroutine add_model_scalar_field

!===============================================================================
!> \function add_model_scalar_field
!
!> \brief add field defining a non-user solved scalar variable,
!>        with default options
!
!> It is recommended not to define variable names of more than 16
!> characters, to get a clear execution listing (some advanced writing
!> levels take into account only the first 16 characters).
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]  name           field name
!> \param[in]  label          field default label, or empty
!> \param[in]  dim            field dimension
!> \param[out] iscal          variable number for defined field
!_______________________________________________________________________________

subroutine add_model_field &
 ( name, label, dim, iscal )

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
integer, intent(out)         :: iscal

! Local variables

integer  ivar, id, ii, ipp
integer  type_flag, location_id,  keycpl
logical  interleaved, has_previous

integer, save :: keyvar = -1
integer, save :: keysca = -1

type_flag = FIELD_INTENSIVE + FIELD_VARIABLE
location_id = 1 ! variables defined on cells
interleaved = .true.
has_previous = .true.

! Test if the field has already been defined
call field_get_id_try(trim(name), id)
if (id .ge. 0) then
  write(nfecra,1000) trim(name)
  call csexit (1)
endif

! Create field

if (keysca.lt.0) then
  call field_get_key_id('coupled', keycpl)
  call field_get_key_id("scalar_id", keysca)
  call field_get_key_id("variable_id", keyvar)
endif

call field_create(name, type_flag, location_id, dim, interleaved, has_previous, &
                  id)

call field_set_key_int(id, keyvis, 1)
call field_set_key_int(id, keylog, 1)

if (len(trim(label)).gt.0) then
  call field_set_key_str(id, keylbl, trim(label))
endif

ivar = nvar + 1
nvar = nvar + dim
nscal = nscal + dim
iscal = nscaus + nscapp + 1
nscapp = nscapp + dim

! Check we have enough slots
call fldvar_check_nvar
call fldvar_check_nscapp

do ii = 1, dim
  isca(iscal + ii - 1) =  nvar - dim + ii
  ivarfl(isca(iscal + ii - 1)) = id
  iscapp(nscapp - dim + ii) = iscal + ii - 1
enddo

ipp = field_post_id(id)

call field_set_key_int(id, keyvar, nvar)
call field_set_key_int(id, keysca, iscal)

if (dim .gt. 1) then
  call field_set_key_int(id, keycpl, 1)
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

end subroutine add_model_field

!===============================================================================

!> \function fldvar_check_nvar

!> \brief check nvarmx is sufficient for the required number of variables.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine fldvar_check_nvar

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

if (nvar .gt. nvarmx) then
  write(nfecra,1000) nvar, nvarmx
  call csexit (1)
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
'@     NOMBRE DE VARIABLES TROP GRAND                         ',/,&
'@                                                            ',/,&
'@  Le type de calcul defini                                  ',/,&
'@    correspond a un nombre de variables NVAR   >= ', i10     ,/,&
'@  Le nombre de variables maximal prevu                      ',/,&
'@                      dans paramx   est NVARMX  = ', i10     ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres                                   ',/,&
'@                                                            ',/,&
'@  Si NVARMX est augmente, le code doit etre reinstalle.     ',/,&
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
'@     NUMBER OF VARIABLES TOO LARGE                          ',/,&
'@                                                            ',/,&
'@  The type of calculation defined                           ',/,&
'@    corresponds to a number of variables NVAR  >= ', i10     ,/,&
'@  The maximum number of variables allowed                   ',/,&
'@                      in   paramx   is  NVARMX  = ', i10     ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@  If NVARMX is increased, the code must be reinstalled.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

end subroutine fldvar_check_nvar

!===============================================================================

!> \function fldvar_check_nscapp
!
!> \brief check nscamx is sufficient for the required number of model scalars.
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine fldvar_check_nscapp

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

if ((nscaus+nscapp).gt.nscamx) then
  write(nfecra,1000) nscaus,nscamx,nscamx-nscaus
  call csexit (1)
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
'@     NOMBRE DE SCALAIRES TROP GRAND                         ',/,&
'@                                                            ',/,&
'@  Le nombre de scalaires utilisateurs                       ',/,&
'@    demande                          est  NSCAUS = ', i10    ,/,&
'@  Le nombre de scalaires total                              ',/,&
'@    autorise   dans paramx           est  NSCAMX = ', i10    ,/,&
'@                                                            ',/,&
'@  La valeur maximale possible de NSCAUS                     ',/,&
'@    avec le modele choisi est donc NSCAMX-NSCAUS = ', i10    ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier NSCAUS.                                          ',/,&
'@                                                            ',/,&
'@  Si NSCAMX est augmente, le code doit etre reinstalle.     ',/,&
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
'@     NUMBER OF SCALARS TOO LARGE                            ',/,&
'@                                                            ',/,&
'@  The number of users scalars                               ',/,&
'@     requested                       is   NSCAUS = ', i10    ,/,&
'@  The total number of scalars                               ',/,&
'@    allowed    in   paramx.h         est  NSCAMX = ', i10    ,/,&
'@                                                            ',/,&
'@  The maximum value possible for NSCAPP                     ',/,&
'@    with the chosen model is       NSCAMX-NSCAUS = ', i10    ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   NSCAUS.                                          ',/,&
'@                                                            ',/,&
'@  If NSCAMX is increased, the code must be reinstalled.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

end subroutine fldvar_check_nscapp

!===============================================================================
