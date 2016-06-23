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

!> \file initi1.f90
!> \brief Commons initialization.
!>
!------------------------------------------------------------------------------

subroutine initi1

!===============================================================================
! Module files
!===============================================================================

use paramx
use optcal
use entsor
use ihmpre
use ppincl, only: ippmod, nmodmx
use post
use cs_c_bindings
use field
use lagran
use ppincl
use cpincl
use dimens
use numvar, only: itempb
use radiat
use cs_fuel_incl

!===============================================================================

implicit none

! Arguments

! Local variables

integer          iok, ipp, nmodpp, imom, n_moments, f_id, f_type, nfld
integer          keyvar, ivar

double precision ttsuit, wtsuit

type(var_cal_opt) vcopt

!===============================================================================

interface

  subroutine gui_postprocess_fields()  &
      bind(C, name='cs_gui_postprocess_fields')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine gui_postprocess_fields

  subroutine gui_linear_solvers()  &
      bind(C, name='cs_gui_linear_solvers')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine gui_linear_solvers

  subroutine user_linear_solvers()  &
      bind(C, name='cs_user_linear_solvers')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine user_linear_solvers

end interface

!===============================================================================
! 1. Initialize modules before user.
!      . entsor
!      . dimens
!      . numvar
!      . pointe
!      . optcal
!      . mltgrd
!      . cstphy
!===============================================================================

call iniini

nmodpp = 0
do ipp = 2, nmodmx
  if (ippmod(ipp).ne.-1) then
    nmodpp = nmodpp+1
  endif
enddo

!===============================================================================
! 2. ENTREE DES DONNEES PAR L'UTILISATEUR
!      ET POSITIONNEMENT DES VARIABLES
!===============================================================================

call iniusi

call ppini1

!===============================================================================
! Map Fortran pointers to C global data
!===============================================================================

call elec_option_init

call cs_rad_transfer_options

! Additional fields

call addfld

! Time moments

call cs_gui_time_moments
call cs_user_time_moments

n_moments = cs_time_moment_n_moments()
do imom = 1, n_moments
  f_id = time_moment_field_id(imom)
  if (f_id.lt.0) cycle
  ipp = field_post_id(f_id)
enddo

! Postprocessing and logging

if (iihmpr.eq.1) then
  call csenso                                                   &
     ( nvppmx, ncapt,  nthist, frhist, iecaux,          &
       ihisvr, tplfmt, xyzcap )
endif

! Restart

ttsuit = -1.d0
wtsuit = -1.d0

call dflsui(ntsuit, ttsuit, wtsuit);

! Lagrangian model options

call lagran_init_map

call lagopt(isuite, iccvfg, iscalt, dtref)

call lagstati

!===============================================================================
! 3. DEFINITION DES COUPLAGES AVEC SYRTHES
!===============================================================================

! Le nombre de couplage SYRTHES doit etre connu avant MODINI a des fins
! de verification de coherence avec la definition des scalaires

if (iihmpr.eq.1) then
  call uisyrc
endif

call ussyrc

call ussatc

!===============================================================================
! 4. MODIFS APRES USINI1
!===============================================================================

call modini

!===============================================================================
! 5. Some additional fields and mappings
!===============================================================================

call fldini

call gui_postprocess_fields

call usipes(nmodpp)

call gui_linear_solvers
call user_linear_solvers

! Number of fields
call field_get_n_fields(nfld)
call field_get_key_id("variable_id", keyvar)

! Copy field calculation options into the field structure
do f_id = 0, nfld - 1

  call field_get_type(f_id, f_type)

  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    call field_get_key_int(f_id, keyvar, ivar)
    if (ivar.gt.0) then
      call field_get_key_struct_var_cal_opt(f_id, vcopt)
      vcopt%iwarni= iwarni(ivar)
      call field_set_key_struct_var_cal_opt(f_id, vcopt)
    endif
  endif
enddo

!===============================================================================
! 6. Coherency checks
!===============================================================================

iok = 0

call verini (iok)

if(iok.gt.0) then
  write(nfecra,9999)iok
  call csexit (1)
else
  write(nfecra,9998)
endif

#if defined(_CS_LANG_FR)

 9998 format(                                                   /,&
' Pas d erreur detectee lors de la verification des donnees'   ,/,&
'               (interface, cs_user_parameters.f90 et autres).',/)
 9999 format(                                                     &
'@'                                                            ,/,&
'@'                                                            ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    ========='                                               ,/,&
'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS' ,/,&
'@'                                                            ,/,&
'@  Le calcul ne sera pas execute (',i10,' erreurs).'          ,/,&
'@'                                                            ,/,&
'@  Se reporter aux impressions precedentes pour plus de'      ,/,&
'@    renseignements.'                                         ,/,&
'@  Verifier les donnees entrees dans l''interface'            ,/,&
'@    et dans les sous-programmes utilisateur.'                ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#else

 9998 format(                                                   /,&
' No error detected during the data verification'              ,/,&
'                          cs_user_parameters.f90 and others).',/)
 9999 format(                                                     &
'@'                                                            ,/,&
'@'                                                            ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION'                ,/,&
'@    ========'                                                ,/,&
'@    THE CALCULATION PARAMETERS ARE INCOHERENT OR INCOMPLET'  ,/,&
'@'                                                            ,/,&
'@  The calculation will not be run (',i10,' errors).'         ,/,&
'@'                                                            ,/,&
'@  See previous impressions for more informations.'           ,/,&
'@  Verify the provided data in the interface'                 ,/,&
'@    and in user subroutines.'                                ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#endif

!===============================================================================
! 7. Output
!===============================================================================

call impini

return
end subroutine
