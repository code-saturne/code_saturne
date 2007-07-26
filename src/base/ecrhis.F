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
C     Copyright (C) 1998-2007 EDF S.A., France
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
                        SUBROUTINE ECRHIS
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , NDIM   , NCELET , NCEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , MODHIS ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , RDEVEL , RTUSER , RA )
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
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE                !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! MODHIS       ! E  !  -> ! INDICATEUR VALANT 0,1 ou 2           !
CARGU !              !    !     ! 1,2 = ECRITURE INTERMEDIAIRE, FINALE |
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! XYZCEN       ! TR !  -> ! POINT ASSOCIES AUX VOLUMES DE CONTROL!
CARGU ! (NDIM,NCELET)!    !     !                                      !
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
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "parall.h"
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
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER        NOMFIC*300, NENVAR*300
      INTEGER          II, II1, II2, LPOS, INAM1, INAM2, LNG
      INTEGER          ICAP,NCAP,IPP,IRA,IPP2,NBPDTE, JTCABS
      INTEGER          IDIVDT, IXMSDT, IDEBIA, IDEBRA, IFINRA, IEL
      INTEGER          NBCAP(NVPPMX)
      DOUBLE PRECISION XTCABS,XYZTMP(3)
      DOUBLE PRECISION VARCAP(NCAPTM)
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
      IPASS = IPASS + 1
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C--> Il n'y a pas eu d'historiques ou s'il n'y a pas de capteur
      IF((IPASS.EQ.1.AND.MODHIS.EQ.2) .OR. NCAPT.EQ.0) RETURN
C
C=======================================================================
C 1. RECHERCHE DES NOEUDS PROCHES -> NODCAP
C=======================================================================
C
      IF(IPASS.EQ.1) THEN
C
        DO II = 1, NCAPT
          CALL FINDPT
C         ===========
     &   ( NCELET, NCEL, XYZCEN ,
     &     XYZCAP(1,II), XYZCAP(2,II), XYZCAP(3,II),
     &     NODCAP(II), NDRCAP(II))
        ENDDO
C
      ENDIF
C
C=======================================================================
C 2. OUVERTURE DU FICHIER DE STOCKAGE hist.tmp
C=======================================================================
C
      IF(IPASS.EQ.1 .AND. IRANGP.LE.0) THEN
        NOMFIC = ' '
        NOMFIC = EMPHIS
        CALL VERLON ( NOMFIC,II1,II2,LPOS)
C       ===========
C
        NOMFIC(II2+1:II2+8) = 'hist.tmp'
        II2 = II2+8
        OPEN ( UNIT=IMPHIS(1), FILE=NOMFIC (II1:II2),
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
        DO IPP = 2, NVPPMX
          IF(IHISVR(IPP,1).NE.0) THEN
            IRA = ABS(IPP2RA(IPP))
C
C     Pour les moments, il faut eventuellement diviser par le temps cumule
            IDIVDT = IPPMOM(IPP)
            IF(IDIVDT.EQ.0) THEN
              IXMSDT = IRA
            ELSE
              IXMSDT = IDEBRA
              IFINRA = IXMSDT + NCEL
              CALL RASIZE ('ECRHIS', IFINRA)
C             ===========
            ENDIF
            IF(IDIVDT.GT.0) THEN
              DO IEL = 1, NCEL
                RA(IXMSDT+IEL-1) = RA(IRA+IEL-1)/
     &               MAX(RA(IDIVDT+IEL-1),EPZERO)
              ENDDO
            ELSEIF(IDIVDT.LT.0) THEN
              DO IEL = 1, NCEL
                RA(IXMSDT+IEL-1) = RA(IRA+IEL-1)/
     &               MAX(DTCMOM(-IDIVDT),EPZERO)
              ENDDO
C           ELSE
C             RA(IXMSDT+IEL-1) = RA(IRA+IEL-1)
C             inutile car on a pose IXMSDT = IRA
            ENDIF
C
            IF(IHISVR(IPP,1).LT.0) THEN
              DO ICAP = 1, NCAPT
                IF (IRANGP.LT.0) THEN
                  VARCAP(ICAP) = RA(IXMSDT+NODCAP(ICAP)-1)
                ELSE
                  CALL PARHIS(NODCAP(ICAP), NDRCAP(ICAP),
C                 ===========
     &                        RA(IXMSDT), VARCAP(ICAP))
                ENDIF
              ENDDO
              NCAP = NCAPT
            ELSE
              DO ICAP = 1, IHISVR(IPP,1)
                IF (IRANGP.LT.0) THEN
                  VARCAP(ICAP) = RA(IXMSDT+NODCAP(IHISVR(IPP,ICAP+1))-1)
                ELSE
                  CALL PARHIS(NODCAP(IHISVR(IPP,ICAP+1)),
C                 ===========
     &                        NDRCAP(IHISVR(IPP,ICAP+1)),
     &                        RA(IXMSDT), VARCAP(ICAP))
                ENDIF
              ENDDO
              NCAP = IHISVR(IPP,1)
            ENDIF
            IF (IRANGP.LE.0) THEN
              WRITE(IMPHIS(1)) NTCABS, TTCABS, (VARCAP(ICAP),
     &                                           ICAP=1,NCAP)
            ENDIF
          ENDIF
        ENDDO
C
      ENDIF
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
C       --> nombre de capteur par variable
        DO IPP = 2, NVPPMX
          NBCAP(IPP) = IHISVR(IPP,1)
          IF(NBCAP(IPP).LT.0) NBCAP(IPP) = NCAPT
        ENDDO
C
C       --> ecriture un fichier par variable
C
        DO IPP = 2, NVPPMX
          IF(IHISVR(IPP,1).NE.0) THEN
C
            IF(IRANGP.LE.0) THEN
C           --> nom du fichier
              NOMFIC = ' '
              NOMFIC = EMPHIS
              CALL VERLON ( NOMFIC,II1,II2,LPOS)
C             ===========
              NENVAR = NOMVAR(IPP)
              CALL VERLON(NENVAR,INAM1,INAM2,LPOS)
C             ===========
              CALL UNDSCR(INAM1,INAM2,NENVAR)
C             ===========
              NOMFIC(II2+1:II2+LPOS) = NENVAR(INAM1:INAM2)
              II2 = II2+LPOS
              NOMFIC(II2+1:II2+1) = '.'
              II2 = II2+1
              NENVAR = EXTHIS
              CALL VERLON(NENVAR,INAM1,INAM2,LPOS)
C             ===========
              CALL UNDSCR(INAM1,INAM2,NENVAR)
C             ===========
              NOMFIC(II2+1:II2+LPOS) = NENVAR(INAM1:INAM2)
              II2 = II2+LPOS
C           --> ouverture
              OPEN ( UNIT=IMPHIS(2), FILE=NOMFIC (II1:II2),
     &               STATUS='UNKNOWN', FORM='FORMATTED',
     &               ACCESS='SEQUENTIAL')
C           --> entete
              WRITE(IMPHIS(2),100)
              WRITE(IMPHIS(2),101)
              WRITE(IMPHIS(2),102) NOMVAR(IPP)
              WRITE(IMPHIS(2),100)
              WRITE(IMPHIS(2),103)
              WRITE(IMPHIS(2),104)
              WRITE(IMPHIS(2),103)
            ENDIF
C
            IF(IHISVR(IPP,1).GT.0) THEN
              DO II=1,IHISVR(IPP,1)
                IF (IRANGP.LT.0 .OR.
     &              IRANGP.EQ.NDRCAP(IHISVR(IPP,II+1))) THEN
                  XYZTMP(1) = XYZCEN(1,NODCAP(IHISVR(IPP,II+1)))
                  XYZTMP(2) = XYZCEN(2,NODCAP(IHISVR(IPP,II+1)))
                  XYZTMP(3) = XYZCEN(3,NODCAP(IHISVR(IPP,II+1)))
                ENDIF
                IF (IRANGP.GE.0) THEN
                  LNG = 3
                  CALL PARBCR(NDRCAP(IHISVR(IPP,II+1)), LNG , XYZTMP)
C                 ===========
                ENDIF
                IF(IRANGP.LE.0) THEN
                  WRITE(IMPHIS(2),105) IHISVR(IPP,II+1),
     &                                 XYZTMP(1), XYZTMP(2), XYZTMP(3)
                ENDIF
              ENDDO
            ELSEIF(IHISVR(IPP,1).LT.0) THEN
              DO II=1,NCAPT
                IF (IRANGP.LT.0 .OR.
     &              IRANGP.EQ.NDRCAP(II)) THEN
                  XYZTMP(1) = XYZCEN(1,NODCAP(II))
                  XYZTMP(2) = XYZCEN(2,NODCAP(II))
                  XYZTMP(3) = XYZCEN(3,NODCAP(II))
                ENDIF
                IF (IRANGP.GE.0) THEN
                  LNG = 3
                  CALL PARBCR(NDRCAP(II), LNG , XYZTMP)
C                 ===========
                ENDIF
                IF(IRANGP.LE.0) THEN
                  WRITE(IMPHIS(2),105) II,
     &                                 XYZTMP(1), XYZTMP(2), XYZTMP(3)
                ENDIF
              ENDDO
            ENDIF
C
            IF(IRANGP.LE.0) THEN
C
              WRITE(IMPHIS(2),103)
              WRITE(IMPHIS(2),106) NBPDTE
              WRITE(IMPHIS(2),103)
C
              WRITE(IMPHIS(2),103)
              WRITE(IMPHIS(2),107)
              WRITE(IMPHIS(2),103)
C
              WRITE(IMPHIS(2),100)
              WRITE(IMPHIS(2),103)
C
C           --> boucle sur les differents enregistrements
C               et les variables
              REWIND(IMPHIS(1))
              DO II = 1, NBPDTE
                DO IPP2 = 2, NVPPMX
                  IF(IHISVR(IPP2,1).NE.0) THEN
                    READ(IMPHIS(1))
     &                JTCABS, XTCABS, (VARCAP(ICAP),ICAP=1,NBCAP(IPP2))
                    IF(IPP2.EQ.IPP)
     &                WRITE(IMPHIS(2),1000)
     &                JTCABS, XTCABS, (VARCAP(ICAP),ICAP=1,NBCAP(IPP))
                  ENDIF
                ENDDO
              ENDDO
C
C           --> fermeture fichier
              CLOSE(IMPHIS(2))
C
            ENDIF
C
          ENDIF
        ENDDO
C
      ENDIF
C
C=======================================================================
C 5. AFFICHAGES
C=======================================================================
C
 100  FORMAT ('# ---------------------------------------------------')
 101  FORMAT ('#      FICHIER HISTORIQUE EN TEMPS')
 102  FORMAT ('#      VARIABLE    ',A16)
 103  FORMAT ('# ')
 104  FORMAT ('#      POSITION DES CAPTEURS (colonne)')
 105  FORMAT ('# ',I6,')',3(1X,E14.7))
 106  FORMAT ('#      NOMBRE D''ENREGISTREMENTS :',I7)
 107  FORMAT (
     &'# COLONNE 1       : NUMERO DU PAS DE TEMPS ',/,
     &'#         2       : TEMPS PHYSIQUE (ou No pas de temps*DTREF ',/,
     &'#                               en pas de temps non uniforme)',/,
     &'#         3 A 100 : VALEUR AUX CAPTEURS')
 1000 FORMAT ( 1(1X,I7,1X),101(1X,E14.7))
C
      RETURN
      END
c@z
