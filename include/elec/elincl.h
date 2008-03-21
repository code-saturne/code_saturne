c@a
c@versb
C-----------------------------------------------------------------------
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2008 EDF S.A., France
C
C     contact: saturne-support@edf.fr
C
C     The Code_Saturne Kernel is free software; you can redistribute it
C     and/or modify it under the terms of the GNU General Public License
C     as published by the Free Software Foundation; either version 2 of
C     the License, or (at your option) any later version.
C
C     The Code_Saturne Kernel is distributed in the hope that it will be
C     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
C     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with the Code_Saturne Kernel; if not, write to the
C     Free Software Foundation, Inc.,
C     51 Franklin St, Fifth Floor,
C     Boston, MA  02110-1301  USA
C
C-----------------------------------------------------------------------
c@verse
C                             elincl.h
C
C***********************************************************************
C
C            INCLUDE POUR LES VERSION ELECTRIQUES
C
C-----------------------------------------------------------------------
C
C--> DEFINITION DES PARAMETERS
C    =========================
C
C     PERMVI : Mu zero, permeabilite magnetique du vide H/m
C     EPSZER : Epsilon zero, permittivite du vide F/m
C
      DOUBLE PRECISION PERMVI            , EPSZER
      PARAMETER      ( PERMVI = 1.2566D-6, EPSZER = 8.854D-12 )
C
C--> DONNEES EN COMMON POUR LE CHAUFFAGE EFFET JOULE
C    ===============================================
C
C     TH, NPOT et NPO sont deja dans ppthch.h
C
C ----- Fournies par l'utilisateur
C       IENTM1       --> indicateur entree matiere premiere
C       IELPH1       --> indicateur electrode phase 1
C       IELPH2       --> indicateur electrode phase 2
C       IELPH3       --> indicateur electrode phase 3
C       IELNEU       --> indicateur electrode neutre
C       ENH          --> tabulation enthalpie(temperature)
C       USRHO        -->  - - - - - inverse masse volumique - - -
C       SIG          -->  - - - - - conductivite  - - - -
C       KAB          -->  - - - - - coeff absorption  - -
C       VIS          -->  - - - - - viscosite - - - - - -
C       LCP          -->  - - - - - Lambda/Cp
C
      INTEGER           IENTM1(NTYPMX),IELPH1(NTYPMX),IELPH2(NTYPMX)
      INTEGER           IELPH3(NTYPMX),IELNEU(NTYPMX)
      COMMON / ICHJOU / IENTM1        ,IELPH1        ,IELPH2       ,
     &                  IELPH3        ,IELNEU
C
C
C       ENHEL        --> tabulation enthalpie      (temperature)
C       RHOEL        -->  - - - - - masse volumique - - -
C       CPEL         -->  - - - - - CP             - - -
C       SIGEL        -->  - - - - - conductivite elec  - - - -
C       XLABEL        -->  - - - - -  conductivite thermique  - -
C       XKABEL        -->  - - - - -  coeff absorption  (pour Srad)- -
C       VISEL        -->  - - - - - viscosite dynamique - - - - - -
C
      DOUBLE PRECISION  RHOEL (NGAZGM,NPOT), CPEL  (NGAZGM,NPOT)
      DOUBLE PRECISION  SIGEL (NGAZGM,NPOT), VISEL (NGAZGM,NPOT)
      DOUBLE PRECISION  XLABEL(NGAZGM,NPOT), XKABEL(NGAZGM,NPOT)
      COMMON / RCHJOU / RHOEL              , CPEL               ,
     &                  SIGEL              , VISEL              ,
     &                  XLABEL             , XKABEL
C
C
C CL sur les electrodes
C
      INTEGER NELEMX,NBTRMX
      PARAMETER (NELEMX = 1000 , NBTRMX = 100)
C
      INTEGER          NBELEC , NBTRF , NTFREF
      COMMON /ELETRF / NBELEC , NBTRF , NTFREF
C
      INTEGER        IELECC(NELEMX),IELECT(NELEMX),IELECB(NELEMX)
      COMMON/ELETRF/IELECC         ,IELECT        ,IELECB
C
      INTEGER       IBRPR(NBTRMX),IBRSEC(NBTRMX)
      COMMON/BRTRSF/IBRPR        ,IBRSEC
C
      DOUBLE PRECISION TENSPR(NBTRMX),RNBS(NBTRMX)
      DOUBLE PRECISION ZR(NBTRMX)    ,ZI(NBTRMX)
      COMMON/CRTRSF/   TENSPR , RNBS , ZR , ZI
C
      DOUBLE PRECISION UROFF(NBTRMX)    ,UIOFF(NBTRMX)
      COMMON/OFFSER/   UROFF            ,UIOFF
C
C--> PARAMETRES POUR LA VERSION ARC ELECTRIQUE
C    ========================================
C
C     IXKABE : valeur lue dans le fichier dp_elec
C             = 0 la derniere colonne du fichier est lue mais pas utilisee
C             = 1 la derniere colonne du fivhier represente le coefficient
C                 d'absorption
C             = 2 la derniere colonne du fivhier represente le TS radiatif
C
      INTEGER           IXKABE
      COMMON / IOPTEL / IXKABE
C
C
C
C    Grandeurs necessaires au claquage
C
C      NTDCLA : iterration de debut du claquage
C      ICLAQ  : indicateur pour savoir si on fait actuellement un claquage
C                = 0 Pas de claquage
C                = 1 Claquage
C       XCLAQ ,YCLAQ ZCLAQ : Position de point de claquage
C
      INTEGER           NTDCLA , ICLAQ
      COMMON / ICLAQU / NTDCLA , ICLAQ
C
      DOUBLE PRECISION  XCLAQ , YCLAQ , ZCLAQ
      COMMON / RCLAQU / XCLAQ , YCLAQ , ZCLAQ

C
C--> DONNEES SUR LA CORRECTION DES VARIABLES ELECTRIQUES
C    EN FONCTION D'UNE INTENSITE DE COURANT DONNEES
C    ========================================
C
C     IELCOR : = 0 pas de correction
C              = 1 correction
C
C     COUIMP : intensite de courant impose par l'utilisateur
C                pour Arc Electrique
C     PUISIM : puissance imposee pour Joule
C     DPOT   : Delta du potentiel electrique entre l'Anode et la cathode
C              (arc et Joule)
C     COEJOU : coefficient de correction pour version Joule
C
      INTEGER           IELCOR
      COMMON / IECORR / IELCOR
C
      DOUBLE PRECISION  COUIMP , DPOT , PUISIM , COEJOU
      COMMON / RECORR / COUIMP , DPOT , PUISIM , COEJOU
C
C--> DONNEES POUR LES ESPECES AYANT UN IMPACT
C    SUR LE PROBLEME ELECTRIQUE
C    ========================================
C
C     QESPEL : Charge massique des especes  C/kg
C     SUSCEP : Susceptibilite (relation champ - mobilite) m2/s/V
C
      DOUBLE PRECISION   QESPEL(NGAZGM), SUSCEP(NGAZGM)
      COMMON / RDPBEL /  QESPEL        , SUSCEP
C
C
c@z
