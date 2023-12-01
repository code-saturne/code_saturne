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

!===============================================================================
! Function :
! --------

!> \file lecamp.f90
!>
!> \brief Reading of main restart file.

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[inout]  oflmap        pointer to old field map
!_______________________________________________________________________________

subroutine lecamp &
 ( oflmap )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

use paramx
use dimens, only: nvar
use cstphy
use cstnum
use entsor
use optcal
use pointe
use numvar
use albase
use parall
use cplsat
use field
use atincl, only: init_at_chem
use atchem, only: ichemistry
use sshaerosol, only: iaerosol, CS_ATMO_AEROSOL_OFF
use turbomachinery
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

type(c_ptr)      oflmap

! Local variables

character(len=64) :: rubriq
character(len=2)  :: cindfp*2

character        ficsui*32
integer          ivar  , ivers(1), f_id
integer          ierror, itysup, nbval
integer          nberro, t_id
integer          nfmtru
integer          jale, jvolfl
integer          ival(1)
double precision rval(1)

logical(kind=c_bool) :: ncelok, nfaiok, nfabok, nsomok

type(c_ptr) :: rp
type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! 0. INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings

ivar = 0

! Memoire

!  ---> Banniere
write(nfecra,1000)

!     Longueur pour format de print
nfmtru = 36

!  --->  On code en chaine le numero des phases et scalaires

!     Indefini a 2 caracteres
cindfp='YY'

!===============================================================================
! 1. OUVERTURE DU FICHIER OU STOP
!===============================================================================

ficsui = 'main.csc'
call restart_create(ficsui, '', 0, rp)

! ---> Debut de la lecture
write(nfecra,1100)

!===============================================================================
! 2. ENTETES DU FICHIER SUITE OU STOP
!===============================================================================

!  --->  Rubrique "fichier suite ppal"
!        Pourrait porter le numero de version si besoin.

itysup = 0
nbval  = 1

call restart_read_int_t_compat(rp,                                      &
                               'code_saturne:checkpoint:main:version',  &
                               'version_fichier_suite_principal',       &
                               itysup, nbval, ivers, ierror)

if (ierror.ne.0) then
  ierror = cs_restart_check_if_restart_from_ncfd(rp)
  if (ierror.eq.0) then
    write(nfecra,9200) ficsui
    call csexit (1)
  endif
endif

!  --->  Tests sur les supports

call restart_check_base_location(rp,ncelok,nfaiok,nfabok,nsomok)

if (ncelok .eqv. .false.) then
  write(nfecra,9201)
  call csexit (1)
endif

! Inutile de tester les supports "faces internes" et "faces de bord"
! ils ne sont pas utilises ici

! Read field info

call restart_read_field_info(rp, oflmap)

! ---> Fin de la lecture des dimensions
write(nfecra,1299)

!===============================================================================
! 4. OPTIONS OU STOP
!===============================================================================

!  --->  Lecture des options

nberro = 0

!     Nombre de pas de temps, instant precedent

rubriq = 'nbre_pas_de_temps'
itysup = 0
nbval  = 1
call restart_read_section_int_t(rp,rubriq,itysup,nbval,ival,ierror)
ntpabs = ival(1) ! no direct read to avoid pointer issue

! If section doesnt exist, check if it is a restart from neptune:
if (ierror.ne.0) then
  rubriq = 'ntcabs'
  call restart_read_section_int_t(rp,rubriq,itysup,nbval,ival,ierror)
  ntpabs = ival(1) ! no direct read to avoid pointer issue
endif
nberro=nberro+ierror

rubriq = 'instant_precedent'
itysup = 0
nbval  = 1
call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
ttpabs = rval(1) ! no direct read to avoid pointer issue

! If section doesnt exist, check if it is a restart from neptune:
if (ierror.ne.0) then
  rubriq = 'ttcabs'
  call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
  ttpabs = rval(1) ! no direct read to avoid pointer issue
endif

nberro=nberro+ierror

! --->  Stop si erreur
if (nberro.ne.0) then
  write(nfecra,9400)
  call csexit (1)
endif

!     Methode ALE

nberro = 0

rubriq = 'methode_ALE'
itysup = 0
nbval  = 1
call restart_read_section_int_t(rp,rubriq,itysup,nbval,ival,ierror)
jale = ival(1)
nberro=nberro+ierror

! --->  Message si erreur (pas de stop pour compatibilite avec les fichiers anterieurs)
!       -> on n'affiche le message que si IALE>=1 (sinon RAS)
if (nberro.ne.0) then
  if (iale.ge.1) write(nfecra,9401)
  jale = 0
endif

!     VOF

nberro = 0

rubriq = 'vof'
itysup = 0
nbval  = 1
call restart_read_section_int_t(rp,rubriq,itysup,nbval,ival,ierror)
jvolfl = ival(1)
nberro=nberro+ierror

! warning if error found (backward compatibility)
! and only if VoF model is enabled
if (nberro.ne.0) then
  if (ivofmt.gt.0) write(nfecra,9405)
  jvolfl = -1
endif

! --->  Stop si pas de temps incoherent
if (ttmabs.ge.0) then
  if (ttpabs.gt.ttmabs) then
    write(nfecra,9411) ttpabs,ttmabs
    call csexit (1)
  endif
else if (ntpabs.gt.ntmabs) then
  write(nfecra,9410) ntpabs,ntmabs
  call csexit (1)
endif

! --->  Informations
write(nfecra,2410) ntpabs
write(nfecra,2411) ttpabs

! --->  Si le calcul precedent etait en ALE, on DOIT relire les
!         coordonnees des noeuds dans le fichier auxiliaire
if (iale.ge.1 .and. jale.ge.1) then
  if (ileaux.ne.1) then
    write(nfecra,9402)jale,iale,ileaux
    call csexit(1)
  endif
endif

!     Instant de maillage mobile precedent (rotor/stator)

nberro = 0

if (iturbo.ne.0) then
  call turbomachinery_restart_read(rp)
endif

! Fin de la lecture des options
write(nfecra,1499)

!===============================================================================
! 5. Read variables
!===============================================================================

call restart_read_variables(rp, oflmap, 0)

f_id = -1
do ivar = 1, nvar
  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
  if (vcopt%ibdtso.gt.1) then
    if (f_id.ne.ivarfl(ivar)) then
      ierror = 0
      f_id = ivarfl(ivar)
      do t_id = 1, vcopt%ibdtso - 1
        call restart_read_field_vals(rp, f_id, t_id, ierror)
        ierror = ierror + 1
      enddo
      if (ierror.gt.1) then
        vcopt%ibdtso = -vcopt%ibdtso
        call field_set_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
      endif
    endif
  endif
enddo

call restart_read_fields(rp, RESTART_MAIN)

!===============================================================================
! 6. LECTURE D'INFORMATIONS COMPLEMENTAIRES LEGERES
!===============================================================================

if (ichemistry.gt.0.or.iaerosol.ne.CS_ATMO_AEROSOL_OFF) then
  rubriq = 'atmospheric_chem'
  itysup = 0
  nbval  = 1
  ival(1) = 1
  call restart_read_section_int_t(rp,rubriq,itysup,nbval,ival,ierror)
  init_at_chem = ival(1)
  if (ierror.eq.0.and.init_at_chem.gt.0) then
    init_at_chem = 0
  endif
endif

!===============================================================================
! 7. FERMETURE DU FICHIER SUITE PRINCIPAL
!===============================================================================

call restart_destroy(rp)

write(nfecra,1799)

!===============================================================================
! 8. SORTIE
!===============================================================================

return

!===============================================================================
! 9. FORMATS
!===============================================================================

! --- ETAPES

 1000 format(/, 3x,'   READING THE MAIN RESTART FILE',/)
 1100 format(' Start reading'  )
 1299 format(' Reading dimensions complete'  )
 1499 format(' Reading options complete')
 1799 format(' Reading complete')

! --- INFORMATIONS

 2410 format                                                            &
 ('  Reading the previous time step number ',                    &
                      '(restarting computation)  NTPABS =   ',I10)
 2411 format                                                            &
 ('  Reading the previous time step number ',                    &
                      '(restarting computation)  TTPABS = ',E12.4)

! --- ERREURS

 9200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@      WRONG FILE TYPE                                       ',/,&
'@                                                            ',/,&
'@    The file ',A13      ,' does not look like a proper      ',/,&
'@      main restart file.                                    ',/,&
'@                                                            ',/,&
'@    The calculation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@    Please make sure the file used as a restart file        ',/,&
'@        actually is a correct main restart file.            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9201 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@      INCONSISTANT RESTART AND CHECKPOINT DATA              ',/,&
'@                                                            ',/,&
'@    The number of cells has changed                         ',/,&
'@                                                            ',/,&
'@    The calculation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@    Please make sure the file used as restart file does     ',/,&
'@        correspond to your case                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9400 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      ERROR AT READING THE TEMPORAL INFORMATION             ',/,&
'@                                                            ',/,&
'@    The computation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@    Please check the integrity of the file used as          ',/,&
'@        restart file                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9401 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : ERROR AT THE MAIN RESTART FILE READING        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      ERROR AT READING THE INDICATOR OF ALE METHOD          ',/,&
'@                                                            ',/,&
'@    The read restart file might come from a previous        ',/,&
'@      version of Code Saturne, without ALE.                 ',/,&
'@    The calculation will be executed but                    ',/,&
'@      ALE data will be reset.                               ',/,&
'@    Please check the integrity of the file used as          ',/,&
'@        restart file, however.                              ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9402 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@      IALE INDICATOR OF THE PREVIOUS CALCULATION = ',I10     ,/,&
'@      IALE INDICATOR OF THE CURRECT CALCULATION  = ',I10     ,/,&
'@                                                            ',/,&
'@    The coordinates of the mesh nodes need to be read again.',/,&
'@      They are stored in the auxiliary restart file.        ',/,&
'@    Therefore the ILEAUX indicator needs to be equal to 1.  ',/,&
'@    Its current value is ILEAUX = ',I10                     ,/, &
'@                                                            ',/,&
'@    The calculation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@    Please check the value of ILEAUX.                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9405 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : ERROR AT THE MAIN RESTART FILE READING        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      ERROR AT READING THE INDICATOR OF THE VOLUME OF FLUID ',/,&
'@                                                      METHOD',/,&
'@    The read restart file might come from a previous        ',/,&
'@      version of Code Saturne, without VOF.                 ',/,&
'@    The calculation will be executed but                    ',/,&
'@      Volume of Fluid method data will be reset.            ',/,&
'@    Please check the integrity of the file used as          ',/,&
'@        restart file, however.                              ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
9410 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@      NUMBER OF THE PREVIOUS TIME STEP  NTPABS = ',I10       ,/,&
'@      NUMBER OF TIME STEPS WANTED       NTMABS = ',I10       ,/,&
'@                                                            ',/,&
'@    The number of time steps (absolute) wanted, NTMABS,     ',/,&
'@      has to be greater or equal to than the number of      ',/,&
'@      time steps (absolute) already run, NTPABS.            ',/,&
'@                                                            ',/,&
'@    The calculation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@    Please check (increase) NTMABS.                         ',/,&
'@    Please make sure the file used as restart file does     ',/,&
'@          correspond to your case                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
9411 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : STOP AT THE MAIN RESTART FILE READING         ',/,&
'@    =========                                               ',/,&
'@      PREVIOUS TIME TTPABS = ',E12.4                         ,/,&
'@      TIME WANTED   TTMABS = ',E12.4                         ,/,&
'@                                                            ',/,&
'@    The number of time steps (absolute) wanted, NTMABS,     ',/,&
'@      has to be greater or equal to than the number of      ',/,&
'@      time steps (absolute) already run, NTPABS.            ',/,&
'@                                                            ',/,&
'@    The calculation cannot be executed.                     ',/,&
'@                                                            ',/,&
'@    Please check (increase) NTMABS.                         ',/,&
'@    Please make sure the file used as restart file does     ',/,&
'@          correspond to your case                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

end subroutine
