!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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
use ppincl, only: ippmod, nmodmx, iatmos
use post
use cs_c_bindings
use field
use lagran
use cpincl
use dimens
use radiat
use cs_fuel_incl
use cdomod

!===============================================================================

implicit none

! Arguments

! Local variables

integer          iok, ipp, nmodpp

double precision ttsuit, wtsuit

!===============================================================================

interface

  subroutine gui_output()  &
      bind(C, name='cs_gui_output')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine gui_output

  subroutine user_finalize_setup_wrapper()  &
      bind(C, name='cs_user_finalize_setup_wrapper')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine user_finalize_setup_wrapper

  subroutine user_syrthes_coupling()  &
      bind(C, name='cs_user_syrthes_coupling')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine user_syrthes_coupling

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

if (ippmod(iatmos).ge.0) call cs_at_data_assim_initialize

! Lagrangian model options

call lagran_init_map

call lagopt(isuite, iccvfg, iscalt, dtref)

! Additional fields if not in CDO mode only

if (icdo.lt.2) then
  call addfld
endif

! Time moments

call cs_gui_time_moments
call cs_user_time_moments

! Postprocessing and logging

call gui_output

! Restart

ttsuit = -1.d0
wtsuit = -1.d0

call dflsui(ntsuit, ttsuit, wtsuit);

!===============================================================================
! 3. MODIFS APRES USINI1
!===============================================================================

! Do not call this routine if CDO mode only (default variables and properties
! are not defined anymore)
if (icdo.lt.2) then
   call modini
endif

!===============================================================================
! 4. Some additional fields and mappings
!===============================================================================

! Do not call this routine if CDO mode only (default variables and properties
! are not defined anymore)
if (icdo.lt.2) then
   call fldini
   call usipes(nmodpp)

   ! Avoid a second spurious call to this function
   ! Call in the C part if CDO is activated, i.e. when
   ! additional geometric quantities and connectivities are built
   if (icdo.lt.0) then
      call user_finalize_setup_wrapper
   endif
endif

!===============================================================================
! 5. Coherency checks
!===============================================================================

iok = 0

! No verification in CDO mode only. This done elsewhere
if (icdo.lt.2) then
   call verini (iok)
   call parameters_check
endif

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
