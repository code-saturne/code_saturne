c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS                  Code_Saturne version 1.3
C                      ------------------------
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
                        SUBROUTINE STRHIS
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , NCELET , NCEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , MODHIS ,
     &   IDEVEL , ITUSER , IA     ,
     &   RDEVEL , RTUSER , RA )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  ROUTINE D'ECRITURE DES HISTORIQUES
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! MODHIS       ! E  !  -> ! INDICATEUR VALANT 0,1 ou 2           !
CARGU !              !    !     ! 1,2 = ECRITURE INTERMEDIAIRE, FINALE |
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
CARGU ! RA           ! TR !  -- ! TABLEAU DES REELS                    !
CARGU !______________!____!_____!______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "parall.h"
      INCLUDE "alstru.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0, IDBRA0
      INTEGER          NDIM, NCELET, NCEL
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          MODHIS
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          NBNAME
      PARAMETER        (NBNAME=12)
      CHARACTER        NOMFIC*300, NENVAR*300
      CHARACTER*80     NAMEVR(NBNAME)
      INTEGER          II, JJ, ISTR, II1, II2, LPOS, INAM1, INAM2
      INTEGER          NBPDTE, JTCABS
      INTEGER          IDEBIA, IDEBRA
      DOUBLE PRECISION XTCABS
      DOUBLE PRECISION VARTMP(NSTRMX)
C
C NOMBRE DE PASSAGES DANS LA ROUTINE
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C=======================================================================
C 0. INITIALISATIONS LOCALES
C=======================================================================
C
C     Seul le processeur 0 ecrit, on n'ecrit rien si pas de structure
      IF (IRANGP.GT.0 .OR. NBSTRU.LE.0) RETURN
C
      IPASS = IPASS + 1
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      NAMEVR(1 ) = "deplacement X"
      NAMEVR(2 ) = "deplacement Y"
      NAMEVR(3 ) = "deplacement Z"
      NAMEVR(4 ) = "vitesse X"
      NAMEVR(5 ) = "vitesse Y"
      NAMEVR(6 ) = "vitesse Z"
      NAMEVR(7 ) = "acceleration X"
      NAMEVR(8 ) = "acceleration Y"
      NAMEVR(9 ) = "acceleration Z"
      NAMEVR(10) = "force X"
      NAMEVR(11) = "force Y"
      NAMEVR(12) = "force Z"

C--> Il n'y a pas eu d'historiques
      IF(IPASS.EQ.1.AND.MODHIS.EQ.2) RETURN
C
C=======================================================================
C 2. OUVERTURE DU FICHIER DE STOCKAGE hist.tmp
C=======================================================================
C
      IF(IPASS.EQ.1) THEN
        NOMFIC = ' '
        NOMFIC = EMPHIS
        CALL VERLON ( NOMFIC,II1,II2,LPOS)
C       ===========
C
        NOMFIC(II2+1:II2+11) = 'histstr.tmp'
        II2 = II2+11
        OPEN ( UNIT=IMPSTH(1), FILE=NOMFIC (II1:II2),
     &         STATUS='UNKNOWN', FORM='UNFORMATTED',
     &         ACCESS='SEQUENTIAL')
      ENDIF
C
C=======================================================================
C 3. ECRITURE DES RESULTATS dans le FICHIER DE STOCKAGE
C=======================================================================
C
      IF(MODHIS.EQ.0.OR.MODHIS.EQ.1) THEN
C
        WRITE(IMPSTH(1)) NTCABS, TTCABS, (XSTR  (1,II),II=1,NBSTRU)
        WRITE(IMPSTH(1)) NTCABS, TTCABS, (XSTR  (2,II),II=1,NBSTRU)
        WRITE(IMPSTH(1)) NTCABS, TTCABS, (XSTR  (3,II),II=1,NBSTRU)
        WRITE(IMPSTH(1)) NTCABS, TTCABS, (XPSTR (1,II),II=1,NBSTRU)
        WRITE(IMPSTH(1)) NTCABS, TTCABS, (XPSTR (2,II),II=1,NBSTRU)
        WRITE(IMPSTH(1)) NTCABS, TTCABS, (XPSTR (3,II),II=1,NBSTRU)
        WRITE(IMPSTH(1)) NTCABS, TTCABS, (XPPSTR(1,II),II=1,NBSTRU)
        WRITE(IMPSTH(1)) NTCABS, TTCABS, (XPPSTR(2,II),II=1,NBSTRU)
        WRITE(IMPSTH(1)) NTCABS, TTCABS, (XPPSTR(3,II),II=1,NBSTRU)
        WRITE(IMPSTH(1)) NTCABS, TTCABS, (FORSTR(1,II),II=1,NBSTRU)
        WRITE(IMPSTH(1)) NTCABS, TTCABS, (FORSTR(2,II),II=1,NBSTRU)
        WRITE(IMPSTH(1)) NTCABS, TTCABS, (FORSTR(3,II),II=1,NBSTRU)
C
      ENDIF
C
C
C=======================================================================
C 4. EN CAS DE SAUVEGARDE INTERMEDIAIRE OU FINALE,
C    TRANSMISSION DES INFORMATIONS DANS LES DIFFERENTS FICHIERS
C=======================================================================
C
C On sauve aussi au premier passage pour permettre une
C     verification des le debut du calcul
C
      IF(MODHIS.EQ.1.OR.MODHIS.EQ.2.OR.IPASS.EQ.1) THEN
C
C       --> nombre de pas de temps enregistres
C
        IF(MODHIS.EQ.2) THEN
          NBPDTE = IPASS - 1
        ELSE
          NBPDTE = IPASS
        ENDIF
C
C       --> ecriture un fichier par variable
        DO II = 1, NBNAME
C
C         --> nom du fichier
          NOMFIC = ' '
          NOMFIC = EMPHIS
          CALL VERLON ( NOMFIC,II1,II2,LPOS)
C         ===========
          NENVAR = 'str_'//NAMEVR(II)
          CALL VERLON(NENVAR,INAM1,INAM2,LPOS)
C         ===========
          CALL UNDSCR(INAM1,INAM2,NENVAR)
C         ===========
          NOMFIC(II2+1:II2+LPOS) = NENVAR(INAM1:INAM2)
          II2 = II2+LPOS
          NOMFIC(II2+1:II2+1) = '.'
          II2 = II2+1
          NENVAR = EXTHIS
          CALL VERLON(NENVAR,INAM1,INAM2,LPOS)
C         ===========
          CALL UNDSCR(INAM1,INAM2,NENVAR)
C         ===========
          NOMFIC(II2+1:II2+LPOS) = NENVAR(INAM1:INAM2)
          II2 = II2+LPOS
C         --> ouverture
          OPEN ( UNIT=IMPSTH(2), FILE=NOMFIC (II1:II2),
     &         STATUS='UNKNOWN', FORM='FORMATTED',
     &         ACCESS='SEQUENTIAL')
C         --> entete
          WRITE(IMPSTH(2),1000)NAMEVR(II),NBSTRU,NBPDTE,NBSTRU+2
          WRITE(IMPSTH(2),2000) (ISTR,ISTR=1,NBSTRU)
          WRITE(IMPSTH(2),2005)
C         --> impression des matrices de masse
          WRITE(IMPSTH(2),2001) ((XMSTRU(1,JJ,ISTR),JJ=1,3),
     &         ISTR=1,NBSTRU)
          WRITE(IMPSTH(2),2004) ((XMSTRU(2,JJ,ISTR),JJ=1,3),
     &         ISTR=1,NBSTRU)
          WRITE(IMPSTH(2),2004) ((XMSTRU(3,JJ,ISTR),JJ=1,3),
     &         ISTR=1,NBSTRU)
          WRITE(IMPSTH(2),2005)
C         --> impression des matrices d'amortissement
          WRITE(IMPSTH(2),2002) ((XCSTRU(1,JJ,ISTR),JJ=1,3),
     &         ISTR=1,NBSTRU)
          WRITE(IMPSTH(2),2004) ((XCSTRU(2,JJ,ISTR),JJ=1,3),
     &         ISTR=1,NBSTRU)
          WRITE(IMPSTH(2),2004) ((XCSTRU(3,JJ,ISTR),JJ=1,3),
     &         ISTR=1,NBSTRU)
          WRITE(IMPSTH(2),2005)
C         --> impression des matrices de raideur
          WRITE(IMPSTH(2),2003) ((XKSTRU(1,JJ,ISTR),JJ=1,3),
     &         ISTR=1,NBSTRU)
          WRITE(IMPSTH(2),2004) ((XKSTRU(2,JJ,ISTR),JJ=1,3),
     &         ISTR=1,NBSTRU)
          WRITE(IMPSTH(2),2004) ((XKSTRU(3,JJ,ISTR),JJ=1,3),
     &         ISTR=1,NBSTRU)
          WRITE(IMPSTH(2),2005)
C         --> impression de fin de section
          WRITE(IMPSTH(2),2006)
C
C
C         --> boucle sur les differents enregistrements
C             et les variables
          REWIND(IMPSTH(1))
          DO JJ = 1, NBPDTE
            DO II1 = 1, NBNAME
              READ(IMPSTH(1))
     &             JTCABS, XTCABS, (VARTMP(ISTR),ISTR=1,NBSTRU)
              IF(II1.EQ.II)
     &             WRITE(IMPSTH(2),3000)
     &             JTCABS, XTCABS, (VARTMP(ISTR),ISTR=1,NBSTRU)
            ENDDO
          ENDDO
C
C         --> fermeture fichier
          CLOSE(IMPSTH(2))
C
        ENDDO
C
      ENDIF
C
C
C----
C FORMATS
C----
 1000 FORMAT(
     &'# ---------------------------------------------------',/,
     &'#      FICHIER HISTORIQUE EN TEMPS                   ',/,
     &'#      VARIABLE ',A30                                 ,/,
     &'# ---------------------------------------------------',/,
     &'# ',/,
     &'#      NOMBRE DE STRUTURES       : ',I8,/,
     &'# ',/,
     &'#      NOMBRE D''ENREGISTREMENTS  : ',I8,/,
     &'# ',/,
     &'# ',/,
     &'# COLONNE 1        : NUMERO DU PAS DE TEMPS',/,
     &'#         2        : TEMPS PHYSIQUE (ou No pas de temps*DTREF',/,
     &'#                               en pas de temps non uniforme)',/,
     &'#         3 A ',I4,' : VALEUR POUR CHAQUE STRUCTURE          ',/,
     &'# ',/,
     &'# ---------------------------------------------------',/,
     &'# ')
 2000 FORMAT('# STRUCTURE    |',20(21X,I3,22X,'|'))
 2001 FORMAT('# MASSE        |',20(1X,G14.7,1X,G14.7,1X,G14.7,' |'))
 2002 FORMAT('# AMORTISSEMENT|',20(1X,G14.7,1X,G14.7,1X,G14.7,' |'))
 2003 FORMAT('# RAIDEUR      |',20(1X,G14.7,1X,G14.7,1X,G14.7,' |'))
 2004 FORMAT('#              |',20(1X,G14.7,1X,G14.7,1X,G14.7,' |'))
 2005 FORMAT('#')
 2006 FORMAT(
     &'# (dans le cas ou les caracteristiques des structures sont   ',/,
     &'#  variables, les valeurs ci-dessus correspondent au dernier ',/,
     &'#  pas de temps d''ecriture dans le fichier',/,
     &'#',/,
     &'# ---------------------------------------------------',/,
     &'# ')
 3000 FORMAT ( 1X,I7,1X,21(1X,E14.7))
C
      RETURN
      END
c@z
C
C
