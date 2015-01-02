!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

!> \file elincl.f90
!> Module for electric arcs

module elincl

  !=============================================================================

  use paramx
  use ppthch

  implicit none

  !=============================================================================


  !--> DEFINITION DES PARAMETERS
  !    =========================

  !     PERMVI : Mu zero, permeabilite magnetique du vide H/m
  !     EPSZER : Epsilon zero, permittivite du vide F/m

  double precision permvi            , epszer
  parameter      ( permvi = 1.2566d-6, epszer = 8.854d-12 )

  !--> DONNEES EN COMMON POUR LE CHAUFFAGE EFFET JOULE
  !    ===============================================

  !     TH, NPOT et NPO sont deja dans ppthch.h

  ! ----- Fournies par l'utilisateur
  !       IENTM1       --> indicateur entree matiere premiere
  !       IELPH1       --> indicateur electrode phase 1
  !       IELPH2       --> indicateur electrode phase 2
  !       IELPH3       --> indicateur electrode phase 3
  !       IELNEU       --> indicateur electrode neutre
  !       ENH          --> tabulation enthalpie(temperature)
  !       USRHO        -->  - - - - - inverse masse volumique - - -
  !       SIG          -->  - - - - - conductivite  - - - -
  !       KAB          -->  - - - - - coeff absorption  - -
  !       VIS          -->  - - - - - viscosite - - - - - -
  !       LCP          -->  - - - - - Lambda/Cp

  integer, save ::           ientm1(ntypmx),ielph1(ntypmx),ielph2(ntypmx)
  integer, save ::           ielph3(ntypmx),ielneu(ntypmx)

  !       ENHEL        --> tabulation enthalpie      (temperature)
  !       RHOEL        -->  - - - - - masse volumique - - -
  !       CPEL         -->  - - - - - CP             - - -
  !       SIGEL        -->  - - - - - conductivite elec  - - - -
  !       XLABEL        -->  - - - - -  conductivite thermique  - -
  !       XKABEL        -->  - - - - -  coeff absorption  (pour Srad)- -
  !       VISEL        -->  - - - - - viscosite dynamique - - - - - -

  double precision, save ::  rhoel (ngazgm,npot), cpel  (ngazgm,npot)
  double precision, save ::  sigel (ngazgm,npot), visel (ngazgm,npot)
  double precision, save ::  xlabel(ngazgm,npot), xkabel(ngazgm,npot)

  ! CL sur les electrodes

  integer nelemx,nbtrmx
  parameter (nelemx = 1000 , nbtrmx = 100)

  integer, save ::        nbelec , nbtrf , ntfref
  integer, save ::        ielecc(nelemx),ielect(nelemx),ielecb(nelemx)

  integer, save ::       ibrpr(nbtrmx),ibrsec(nbtrmx)

  double precision, save :: tenspr(nbtrmx),rnbs(nbtrmx)
  double precision, save :: zr(nbtrmx)    ,zi(nbtrmx)

  double precision, save :: uroff(nbtrmx)    ,uioff(nbtrmx)

  !--> PARAMETRES POUR LA VERSION ARC ELECTRIQUE
  !    ========================================

  !     IXKABE : valeur lue dans le fichier dp_elec
  !             = 0 la derniere colonne du fichier est lue mais pas utilisee
  !             = 1 la derniere colonne du fivhier represente le coefficient
  !                 d'absorption
  !             = 2 la derniere colonne du fivhier represente le TS radiatif

  integer, save ::           ixkabe

  !    Grandeurs necessaires au claquage

  !      NTDCLA : iterration de debut du claquage
  !      ICLAQ  : indicateur pour savoir si on fait actuellement un claquage
  !                = 0 Pas de claquage
  !                = 1 Claquage
  !       XCLAQ ,YCLAQ ZCLAQ : Position de point de claquage

  integer, save ::           ntdcla , iclaq

  double precision, save ::  xclaq , yclaq , zclaq

  !--> DONNEES SUR LA CORRECTION DES VARIABLES ELECTRIQUES
  !    EN FONCTION D'UNE INTENSITE DE COURANT DONNEES
  !    ========================================

  !     IELCOR : = 0 pas de correction
  !              = 1 correction

  !     COUIMP : intensite de courant impose par l'utilisateur
  !                pour Arc Electrique
  !     PUISIM : puissance imposee pour Joule
  !     DPOT   : Delta du potentiel electrique entre l'Anode et la cathode
  !              (arc et Joule)
  !     COEJOU : coefficient de correction pour version Joule
  !     MODREC : modele de recalage de l'arc
  !              1 : modele cas general
  !              2 : modele avec un plan
  !     IZRECA : definition du plan de recalage
  !     IDRECA : defnition de la compsante a recaler
  !     CRIT_RECA : define criteria for recal

  double precision, save :: crit_reca(5)
  integer, save ::           ielcor, modrec, idreca
  integer, allocatable, dimension(:) :: izreca

  double precision, save ::  couimp , dpot , puisim , coejou, elcou

  !--> DONNEES POUR LES ESPECES AYANT UN IMPACT
  !    SUR LE PROBLEME ELECTRIQUE
  !    ========================================

  !     QESPEL : Charge massique des especes  C/kg
  !     SUSCEP : Susceptibilite (relation champ - mobilite) m2/s/V

  double precision, save ::   qespel(ngazgm), suscep(ngazgm)

  !=============================================================================

contains

  !=============================================================================

subroutine init_elec

  use mesh
  implicit none
  integer iel
  allocate(izreca(nfac))

  do iel = 1, nfac
    izreca(iel) = 0
  enddo

end subroutine init_elec

  !=============================================================================

subroutine finalize_elec

  implicit none
  deallocate(izreca)

end subroutine finalize_elec

end module elincl
