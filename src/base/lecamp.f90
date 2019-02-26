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
use siream, only: iaerosol
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
integer          jale, jcavit, jvolfl
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

ficsui = 'main'
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

!     Cavitation

nberro = 0

rubriq = 'cavitation'
itysup = 0
nbval  = 1
call restart_read_section_int_t(rp,rubriq,itysup,nbval,ival,ierror)
jcavit = ival(1)
nberro=nberro+ierror

! Message si erreur (pas de stop pour compatibilite avec fichiers anterieurs)
! -> on n'affiche le message que si ICAVIT>=0 (sinon RAS)
if (nberro.ne.0) then
  if (icavit.ge.0) write(nfecra,9404)
  jcavit = -1
endif

!     VOF

nberro = 0

rubriq = 'vof'
itysup = 0
nbval  = 1
call restart_read_section_int_t(rp,rubriq,itysup,nbval,ival,ierror)
jvolfl = ival(1)
nberro=nberro+ierror

! Message si erreur (pas de stop pour compatibilite avec fichiers anterieurs)
! -> on n'affiche le message que si IVOFMT>=0 (sinon RAS)
if (nberro.ne.0) then
  if (ivofmt.ge.0) write(nfecra,9405)
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

rubriq = 'instant_mobile_precedent'
itysup = 0
nbval  = 1
call restart_read_section_real_t(rp,rubriq,itysup,nbval,rval,ierror)
ttpmob = rval(1) ! no direct read to avoid pointer issue
nberro=nberro+ierror

! Message si erreur (pas de stop pour compatibilite avec fichiers anterieurs)
! -> on n'affiche le message que si iturbo=2 (sinon RAS)
if (nberro.ne.0) then
  if (iturbo.eq.2) write(nfecra,9403) ttpabs
  ttpmob = ttpabs
endif

! Information (uniquement si iturbo=2 et pas d affichage precedent)
if (iturbo.eq.2) then
  if (nberro.eq.0)  write(nfecra,2412) ttpmob
endif

call turbomachinery_restart_read(rp)

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

if (ichemistry.gt.0.or.iaerosol.gt.0) then
  rubriq = 'atmospheric_chem'
  itysup = 0
  nbval  = 1
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

#if defined(_CS_LANG_FR)

 1000 format(/, 3x,'   LECTURE DU FICHIER SUITE PRINCIPAL',/)
 1100 format(' Debut de la lecture')
 1299 format(' Fin de la lecture des dimensions')
 1499 format(' Fin de la lecture des options')
 1799 format(' Fin de la lecture')

#else

 1000 format(/, 3x,'   READING THE MAIN RESTART FILE',/)
 1100 format(' Start reading'  )
 1299 format(' Reading dimensions complete'  )
 1499 format(' Reading options complete')
 1799 format(' Reading complete')

#endif

! --- INFORMATIONS

#if defined(_CS_LANG_FR)

 2410 format                                                            &
 ('  Lecture du pas de temps precedent (suite) ',                &
                                                  'NTPABS = ',I10)
 2411 format                                                            &
 ('  Lecture du pas de temps precedent (suite) ',                &
                                                'TTPABS = ',E12.4)
 2412 format                                                            &
 ('  Lecture du temps de maillage mobile precedent (suite) ',    &
                                                'TTPMOB = ',E12.4)

#else

 2410 format                                                            &
 ('  Reading the previous time step number ',                    &
                      '(restarting computation)  NTPABS =   ',I10)
 2411 format                                                            &
 ('  Reading the previous time step number ',                    &
                      '(restarting computation)  TTPABS = ',E12.4)
 2412 format                                                            &
 ('  Reading the previous moving mesh moment ',                  &
                      '(restarting computation)  TTPMOB = ',E12.4)

#endif

! --- ERREURS

#if defined(_CS_LANG_FR)

 9200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@      TYPE DE FICHIER INCORRECT                             ',/,&
'@                                                            ',/,&
'@    Le fichier ',A13      ,' ne semble pas etre un fichier  ',/,&
'@      suite principal.                                      ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise correspond bien   ',/,&
'@        a un fichier suite principal.                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9201 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    Le nombre de cellules a ete modifie                     ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise correspond bien   ',/,&
'@        au cas traite.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9400 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@                                                            ',/,&
'@      ERREUR A LA LECTURE DES INFORMATIONS TEMPORELLES      ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise n''a pas ete      ',/,&
'@        endommage.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9401 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A LA LECTURE DU FICHIER SUITE        ',/,&
'@    =========                                      PRINCIPAL',/,&
'@                                                            ',/,&
'@      ERREUR A LA LECTURE DE L''INDICATEUR DE METHODE ALE   ',/,&
'@                                                            ',/,&
'@    Il se peut que le fichier suite relu corresponde a une  ',/,&
'@      version anterieure de Code_Saturne, sans methode ALE. ',/,&
'@    Le calcul sera execute en reinitialisant toutes les     ',/,&
'@      donnees ALE.                                          ',/,&
'@    Verifier neanmoins que le fichier suite utilise n''a    ',/,&
'@        pas ete endommage.                                  ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9402 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@      INDICATEUR IALE DU CALCUL PRECEDENT = ',I10            ,/,&
'@      INDICATEUR IALE DU CALCUL ACTUEL    = ',I10            ,/,&
'@                                                            ',/,&
'@    Les coordonnees des noeuds du maillage doivent etre     ',/,&
'@      relues. Elles sont stockees dans le fichier suite     ',/,&
'@      auxiliaire.                                           ',/,&
'@    L''indicateur ILEAUX doit donc etre positionne a 1.     ',/,&
'@    Il vaut ici ILEAUX = ',I10                               ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@    Verifier ILEAUX.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9403 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A LA LECTURE DU FICHIER SUITE        ',/,&
'@    =========                                      PRINCIPAL',/,&
'@                                                            ',/,&
'@      ERREUR A LA LECTURE DE L INSTANT DE MAILLAGE MOBILE  ',/,&
'@                                                   PRECEDENT',/,&
'@    Il se peut que le fichier suite relu corresponde a une  ',/,&
'@      version anterieure de Code_Saturne, sans couplage     ',/,&
'@      rotor/stator instationnaire.                          ',/,&
'@    Le calcul sera execute en initialisant l instant de     ',/,&
'@      maillage mobile precedent a TTCMOB = ',E12.4           ,/,&
'@    Verifier neanmoins que le fichier suite utilise n''a    ',/,&
'@        pas ete endommage.                                  ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9404 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A LA LECTURE DU FICHIER SUITE        ',/,&
'@    =========                                      PRINCIPAL',/,&
'@                                                            ',/,&
'@      ERREUR A LA LECTURE DE L''INDICATEUR DU MODELE DE     ',/,&
'@                                                  CAVITATION',/,&
'@                                                            ',/,&
'@    Il se peut que le fichier suite relu corresponde a une  ',/,&
'@      version anterieure de Code_Saturne, sans cavitation.  ',/,&
'@    Le calcul sera execute en reinitialisant toutes les     ',/,&
'@      donnees du modele de cavitation.                      ',/,&
'@    Verifier neanmoins que le fichier suite utilise n''a    ',/,&
'@        pas ete endommage.                                  ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9405 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A LA LECTURE DU FICHIER SUITE        ',/,&
'@    =========                                      PRINCIPAL',/,&
'@                                                            ',/,&
'@      ERREUR A LA LECTURE DE L''INDICATEUR DE LA METHODE    ',/,&
'@                                             VOLUME of FLUID',/,&
'@                                                            ',/,&
'@    Il se peut que le fichier suite relu corresponde a une  ',/,&
'@      version anterieure de Code_Saturne, sans methode      ',/,&
'@      Volume of Fluid.                                      ',/,&
'@    Le calcul sera execute en reinitialisant toutes les     ',/,&
'@      donnees de la methode Volume of Fluid.                ',/,&
'@    Verifier neanmoins que le fichier suite utilise n''a    ',/,&
'@        pas ete endommage.                                  ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9410 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@      NUMERO DU PAS DE TEMPS PRECEDENT NTPABS = ',I10        ,/,&
'@      NUMERO DU PAS DE TEMPS VISE      NTMABS = ',I10        ,/,&
'@                                                            ',/,&
'@    Le nombre de pas de temps (absolu) vise, NTMABS,        ',/,&
'@      doit etre superieur ou egal au nombre de pas de temps ',/,&
'@      (absolu) deja effectues, NTPABS.                      ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier (augmenter) NTMABS.                            ',/,&
'@    Verifier que le fichier suite utilise correspond bien   ',/,&
'@        au cas traite.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9411 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                      PRINCIPAL',/,&
'@      TEMPS PRECEDENT TTPABS = ',E12.4                      ,/,&
'@      TEMPS VISE      TTMABS = ',E12.4                       ,/,&
'@                                                            ',/,&
'@    Le nombre de pas de temps (absolu) vise, NTMABS,        ',/,&
'@      doit etre superieur ou egal au nombre de pas de temps ',/,&
'@      (absolu) deja effectues, NTPABS.                      ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier (augmenter) NTMABS.                            ',/,&
'@    Verifier que le fichier suite utilise correspond bien   ',/,&
'@        au cas traite.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
! NFMTRU = 36 pour A36

#else

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
 9403 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : ERROR AT THE MAIN RESTART FILE READING        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      ERROR AT READING THE PREVIOUS MOVING MESH MOMENT      ',/,&
'@                                                            ',/,&
'@    The read restart file might come from a previous        ',/,&
'@      version of Code Saturne, without unsteady             ',/,&
'@      rotor/stator coupling method.                         ',/,&
'@    The calculation will be executed with the previous      ',/,&
'@      moving mesh moment initialized to TTCMOB = ',E12.4     ,/,&
'@    Please check the integrity of the file used as          ',/,&
'@        restart file, however.                              ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9404 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING : ERROR AT THE MAIN RESTART FILE READING        ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@      ERROR AT READING THE INDICATOR OF THE CAVITATION MODEL',/,&
'@                                                            ',/,&
'@    The read restart file might come from a previous        ',/,&
'@      version of Code Saturne, without cavitation.          ',/,&
'@    The calculation will be executed but                    ',/,&
'@      cavitation model data will be reset.                  ',/,&
'@    Please check the integrity of the file used as          ',/,&
'@        restart file, however.                              ',/,&
'@                                                            ',/,&
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

#endif

end subroutine
