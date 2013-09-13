!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

!===============================================================================

implicit none

! Arguments

! Local variables

integer          iok, ipp, nmodpp
double precision ttsuit, wtsuit

!===============================================================================

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
!==========

nmodpp = 0
do ipp = 2, nmodmx
  if (ippmod(ipp).ne.-1) then
    nmodpp = nmodpp+1
  endif
enddo

!===============================================================================
! 2. ENTREE DES DONNEES PAR L'UTILISATEUR
!      ET POSITIONNEMENT DES VARIABLES (VARPOS)
!===============================================================================

call iniusi
!==========

call ppini1
!==========

call usipes(nmodpp)
!==========

ttsuit = -1.d0
wtsuit = -1.d0

call dflsui(ntsuit, ttsuit, wtsuit);
!==========

call rayopt
!==========

call lagopt
!==========

!===============================================================================
! 3. DEFINITION DES COUPLAGES AVEC SYRTHES
!===============================================================================

! Le nombre de couplage SYRTHES doit etre connu avant MODINI a des fins
! de verification de coherence avec la definition des scalaires

if (iihmpr.eq.1) then
  call uisyrc
  !==========
endif

call ussyrc
!==========

call ussatc
!==========

!===============================================================================
! 4. MODIFS APRES USINI1
!===============================================================================

call modini
!==========

!===============================================================================
! 5. Initial definition of fields
!===============================================================================

call fldini
!==========

call user_field_parameters
!=========================

call addfld

!===============================================================================
! 6. Coherency checks
!===============================================================================

iok = 0

call verini (iok)
!==========

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
!==========

!===============================================================================
! 7. Other initializations
!===============================================================================

call cs_log_init_moments(dtcmom)   ! Map dtcmom to C logging API
call cs_post_init_moments(dtcmom)  ! Map dtcmom to C post-processing API

return
end subroutine
