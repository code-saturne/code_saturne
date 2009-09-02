!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

!                             elincl.h

!===============================================================================

!            INCLUDE POUR LES VERSION ELECTRIQUES

!-------------------------------------------------------------------------------

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

integer           ientm1(ntypmx),ielph1(ntypmx),ielph2(ntypmx)
integer           ielph3(ntypmx),ielneu(ntypmx)
common / ichjou / ientm1        ,ielph1        ,ielph2       ,    &
                  ielph3        ,ielneu


!       ENHEL        --> tabulation enthalpie      (temperature)
!       RHOEL        -->  - - - - - masse volumique - - -
!       CPEL         -->  - - - - - CP             - - -
!       SIGEL        -->  - - - - - conductivite elec  - - - -
!       XLABEL        -->  - - - - -  conductivite thermique  - -
!       XKABEL        -->  - - - - -  coeff absorption  (pour Srad)- -
!       VISEL        -->  - - - - - viscosite dynamique - - - - - -

double precision  rhoel (ngazgm,npot), cpel  (ngazgm,npot)
double precision  sigel (ngazgm,npot), visel (ngazgm,npot)
double precision  xlabel(ngazgm,npot), xkabel(ngazgm,npot)
common / rchjou / rhoel              , cpel               ,       &
                  sigel              , visel              ,       &
                  xlabel             , xkabel


! CL sur les electrodes

integer nelemx,nbtrmx
parameter (nelemx = 1000 , nbtrmx = 100)

integer          nbelec , nbtrf , ntfref
common /eletrf / nbelec , nbtrf , ntfref

integer        ielecc(nelemx),ielect(nelemx),ielecb(nelemx)
common/eletrf/ielecc         ,ielect        ,ielecb

integer       ibrpr(nbtrmx),ibrsec(nbtrmx)
common/brtrsf/ibrpr        ,ibrsec

double precision tenspr(nbtrmx),rnbs(nbtrmx)
double precision zr(nbtrmx)    ,zi(nbtrmx)
common/crtrsf/   tenspr , rnbs , zr , zi

double precision uroff(nbtrmx)    ,uioff(nbtrmx)
common/offser/   uroff            ,uioff

!--> PARAMETRES POUR LA VERSION ARC ELECTRIQUE
!    ========================================

!     IXKABE : valeur lue dans le fichier dp_elec
!             = 0 la derniere colonne du fichier est lue mais pas utilisee
!             = 1 la derniere colonne du fivhier represente le coefficient
!                 d'absorption
!             = 2 la derniere colonne du fivhier represente le TS radiatif

integer           ixkabe
common / ioptel / ixkabe



!    Grandeurs necessaires au claquage

!      NTDCLA : iterration de debut du claquage
!      ICLAQ  : indicateur pour savoir si on fait actuellement un claquage
!                = 0 Pas de claquage
!                = 1 Claquage
!       XCLAQ ,YCLAQ ZCLAQ : Position de point de claquage

integer           ntdcla , iclaq
common / iclaqu / ntdcla , iclaq

double precision  xclaq , yclaq , zclaq
common / rclaqu / xclaq , yclaq , zclaq


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

integer           ielcor
common / iecorr / ielcor

double precision  couimp , dpot , puisim , coejou, elcou
common / recorr / couimp , dpot , puisim , coejou, elcou

!--> DONNEES POUR LES ESPECES AYANT UN IMPACT
!    SUR LE PROBLEME ELECTRIQUE
!    ========================================

!     QESPEL : Charge massique des especes  C/kg
!     SUSCEP : Susceptibilite (relation champ - mobilite) m2/s/V

double precision   qespel(ngazgm), suscep(ngazgm)
common / rdpbel /  qespel        , suscep


