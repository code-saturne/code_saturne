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
                        SUBROUTINE LAGHIS
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , NDIM   , NCELET , NCEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , MODHIS , NVLSTA ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , VOLUME , STATIS , STATIV ,
     &   RDEVEL , RTUSER , RA )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  ROUTINE D'ECRITURE DES HISTORIQUES POUR LE LAGRANGIEN
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
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0, IDBRA0
      INTEGER          NDIM, NCELET, NCEL
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE  , NVLSTA
      INTEGER          MODHIS
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
      DOUBLE PRECISION XYZCEN(NDIM,NCELET) , VOLUME(NCELET)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA)
      DOUBLE PRECISION STATIV(NCELET,NVLSTA-1)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER        NOMFIC*300, NENVAR*300
      INTEGER          II, II1, II2, LPOS, INAM1, INAM2, LNG
      INTEGER          ICAP,NCAP,IPP,IPP2,NBPDTE, JTCABS
      INTEGER          IDMOYD, IDEBIA, IDEBRA, IFINIA ,IFINRA
      INTEGER          IEL   , IVARL
      INTEGER          NBCAP(NVPPMX)
      INTEGER          IPAS  , ILPD1 , IL , ILFV1 , ICLA , ILTS1
C
      DOUBLE PRECISION XTCABS,XYZTMP(3)
      DOUBLE PRECISION VARCAP(NCAPTM)
      DOUBLE PRECISION DMOY
C
C NOMBRE DE PASSAGES DANS LA ROUTINE
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C***********************************************************************
C=======================================================================
C 0. INITIALISATIONS LOCALES
C=======================================================================
C
      IPASS = IPASS + 1
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C Test : Si il n'y a pas de capteur ====> On ne fait rien
C
      IF(NCAPT.LE.0) RETURN
C
C=======================================================================
C 2. OUVERTURE DU FICHIER DE STOCKAGE histla.tmp
C=======================================================================
C
      IF(IPASS.EQ.1 .AND. IRANGP.LE.0) THEN
        NOMFIC = ' '
        NOMFIC = EMPHIS
        CALL VERLON ( NOMFIC,II1,II2,LPOS)
C       ===========
C
        NOMFIC(II2+1:II2+10) = 'histla.tmp'
        II2 = II2+10
        OPEN ( UNIT=IMPLI1, FILE=NOMFIC (II1:II2),
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
        IFINIA = IDEBIA
        IDMOYD = IDEBRA
        IFINRA = IDMOYD + NCELET
        CALL RASIZE('LAGHIS',IFINRA)
C       ===========
C
        DO IPAS  = 1,1+NBCLST
C
C   Moyenne
C
        DO IL = 1,NVLSTA
C
          IVARL = (IPAS-1)*NVLSTA+IL
          ILPD1 = (IPAS-1)*NVLSTA+ILPD
          ILFV1 = (IPAS-1)*NVLSTA+ILFV
          ICLA  = IPAS -1
C
C Pour l'instant on fait des chrono sur toutes les variables Stat. Lag.
C et sur tout les capteurs
C
         IF ( IHSLAG(IVARL).GE. 1 ) THEN
C
           IF ( IVARL .NE. ILPD1 .AND. IVARL .NE. ILFV1 ) THEN
             DO IEL=1,NCEL
               IF ( STATIS(IEL,ILPD1) .GT. SEUIL ) THEN
                 RA(IDMOYD+IEL-1) = STATIS(IEL,IVARL)
     &                             /STATIS(IEL,ILPD1)
               ELSE
                 RA(IDMOYD+IEL-1) = 0.D0
               ENDIF
             ENDDO
C
           ELSE IF (IVARL.EQ.ILPD1) THEN
             DO IEL=1,NCEL
               RA(IDMOYD+IEL-1) = STATIS(IEL,IVARL)
             ENDDO
           ELSE
             DO IEL=1,NCEL
               IF (NPST.GT.0) THEN
                 RA(IDMOYD+IEL-1) = STATIS(IEL,IVARL)
     &                             /(DBLE(NPST)*VOLUME(IEL))
               ELSE
                 RA(IDMOYD+IEL-1) = 0.D0
               ENDIF
             ENDDO
           ENDIF
C
           DO ICAP = 1, NCAPT
             IF (IRANGP.LT.0) THEN
               IF ( IVARL .NE. ILPD1 .AND. IVARL .NE. ILFV1 ) THEN
                 VARCAP(ICAP) = RA(IDMOYD+NODCAP(ICAP)-1)
               ENDIF
             ELSE
               CALL PARHIS(NODCAP(ICAP), NDRCAP(ICAP),
C              ===========
     &                     RA(IDMOYD), VARCAP(ICAP))
             ENDIF
           ENDDO
           NCAP = NCAPT
C
           IF (IRANGP.LE.0) THEN
             WRITE(IMPLI1) NTCABS, TTCABS, (VARCAP(ICAP),
     &                                           ICAP=1,NCAP)
C

           ENDIF
C
         ENDIF
        ENDDO
C
C   Variance
C
        DO IL = 1,NVLSTA-1
C
         IVARL = (IPAS-1)*NVLSTA+IL
         ILPD1 = (IPAS-1)*NVLSTA+ILPD
         ILFV1 = (IPAS-1)*NVLSTA+ILFV
         ILTS1 = (IPAS-1)*NVLSTA+ILTS
         ICLA  = IPAS -1
C
C Pour l'instant on fait des chrono sur toutes les variables Stat. Lag.
C et sur tout les capteurs
C
         IF ( IHSLAG(IVARL).EQ. 2 ) THEN
           DO IEL = 1, NCEL
C
            IF ( IVARL.NE.ILFV ) THEN
              IF ( STATIS(IEL,ILPD1).GT.SEUIL ) THEN
                RA(IDMOYD+IEL-1) = STATIV(IEL,IVARL)/STATIS(IEL,ILPD1)
     &                       -( STATIS(IEL,IVARL)/STATIS(IEL,ILPD1)
     &                         *STATIS(IEL,IVARL)/STATIS(IEL,ILPD1) )
              ELSE
                RA(IDMOYD+IEL-1) = ZERO
              ENDIF
            ELSE
              IF ( STATIS(IEL,ILPD1).GT.SEUIL .AND. NPST.GT.0 ) THEN
                DMOY = STATIS(IEL,IVARL)
     &                /(DBLE(NPST)*VOLUME(IEL))
                RA(IDMOYD+IEL-1) = STATIV(IEL,IVARL)
     &                             /( DBLE(NPST) * VOLUME(IEL) )
     &                             -DMOY*DMOY
C
              ELSE IF ( STATIS(IEL,ILPD1).GT.SEUIL .AND.
     &                  IPLAS.GE.IDSTNT                  ) THEN
                DMOY =  STATIS(IEL,IVARL) / VOLUME(IEL)
                RA(IDMOYD+IEL-1) = STATIV(IEL,ILFV) / VOLUME(IEL)
     &                            -DMOY*DMOY
              ELSE
                RA(IDMOYD+IEL-1) = ZERO
              ENDIF
            ENDIF
            RA(IDMOYD+IEL-1) = SQRT( MAX(ZERO,RA(IDMOYD+IEL-1)))
          ENDDO
C
          DO ICAP = 1, NCAPT
            IF (IRANGP.LT.0) THEN
              VARCAP(ICAP) = RA(IDMOYD+NODCAP(ICAP)-1)
            ELSE
              CALL PARHIS(NODCAP(ICAP), NDRCAP(ICAP),
C             ===========
     &                    RA(IDMOYD), VARCAP(ICAP))
            ENDIF
          ENDDO
          NCAP = NCAPT
C
          IF (IRANGP.LE.0) THEN
            WRITE(IMPLI1) NTCABS, TTCABS, (VARCAP(ICAP),
     &                                         ICAP=1,NCAP)
          ENDIF
C
         ENDIF
        ENDDO
C
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
        DO IPP = 1, 2*NVLSTA
          NBCAP(IPP) = NCAPT
        ENDDO
C
C       --> ecriture un fichier par variable
C
C
        DO IPAS  = 1,1+NBCLST
C
         DO IPP = 1, 2*NVLSTA
C
          IF ( IPP .LE. NVLSTA) THEN
            IVARL = (IPAS-1)*NVLSTA+IL
            ILPD1 = (IPAS-1)*NVLSTA+ILPD
            ILFV1 = (IPAS-1)*NVLSTA+ILFV
            ILTS1 = (IPAS-1)*NVLSTA+ILTS
            ICLA  = IPAS -1
          ELSE
            IVARL = (IPAS-1)*NVLSTA+(IL-NVLSTA)
            ILPD1 = (IPAS-1)*NVLSTA+ILPD
            ILFV1 = (IPAS-1)*NVLSTA+ILFV
            ILTS1 = (IPAS-1)*NVLSTA+ILTS
            ICLA  = IPAS -1
          ENDIF
C
          IF ( (IPP .LE. NVLSTA .AND. IHSLAG(IPP).GE. 1 )
     &        .OR.
     &         (IPP .GT. NVLSTA .AND. (IPP-NVLSTA).NE.ILPD
     &                     .AND. IHSLAG(IPP-NVLSTA).EQ. 2 ) )THEN
C
            IF(IRANGP.LE.0) THEN
C           --> nom du fichier
              NOMFIC = ' '
              NOMFIC = EMPHIS
              CALL VERLON (NOMFIC,II1,II2,LPOS)
C             ===========
C
              IF ( IPAS.EQ.1 ) THEN
                IF ( IPP.LE.NVLSTA ) THEN
                  NENVAR = NOMLAG(IPP)
                ELSE
                  NENVAR = NOMLAV(IPP-NVLSTA)
                ENDIF
              ELSE
                IF ( IPP.LE.NVLSTA ) THEN
                  WRITE(NENVAR,'(A8,A4,I3)')
     &                     NOMLAG(IPP),'_grp',ICLA
                ELSE
                  WRITE(NENVAR,'(A8,A4,I3)')
     &                     NOMLAV(IPP-NVLSTA),'_grp',ICLA
                ENDIF
              ENDIF
              CALL VERLON(NENVAR,INAM1,INAM2,LPOS)
C             ===========
              CALL UNDSCR(INAM1,INAM2,NENVAR)
C             ===========
              NOMFIC(II2+1:II2+4)='Lag_'
              NOMFIC(II2+4+1:II2+4+LPOS) = NENVAR(INAM1:INAM2)
              II2 = II2+4+LPOS
              NOMFIC(II2+1:II2+1) = '.'
              II2 = II2+1
              NENVAR = EXTHIS
              CALL VERLON(NENVAR,INAM1,INAM2,LPOS)
C             ===========
              CALL UNDSCR(INAM1,INAM2,NENVAR)
C             ===========
              NOMFIC(II2+1:II2+LPOS) = NENVAR(INAM1:INAM2)
              II2 = II2+LPOS
C             --> ouverture
              OPEN ( UNIT=IMPLI2, FILE=NOMFIC (II1:II2),
     &               STATUS='UNKNOWN', FORM='FORMATTED',
     &               ACCESS='SEQUENTIAL')
C             --> entete
              WRITE(IMPLI2,100)
              WRITE(IMPLI2,101)
              WRITE(IMPLI2,102) NOMLAG(IPP)
              WRITE(IMPLI2,100)
              WRITE(IMPLI2,103)
              WRITE(IMPLI2,104)
              WRITE(IMPLI2,103)
            ENDIF
C
            DO II=1,NCAPT
              IF (IRANGP.LT.0 .OR.
     &          IRANGP.EQ.NDRCAP(II)) THEN
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
                WRITE(IMPLI2,105) II,
     &                            XYZTMP(1), XYZTMP(2), XYZTMP(3)
              ENDIF
            ENDDO
C
            IF(IRANGP.LE.0) THEN
C
              WRITE(IMPLI2,103)
              WRITE(IMPLI2,106) NBPDTE
              WRITE(IMPLI2,103)
C
              WRITE(IMPLI2,103)
              WRITE(IMPLI2,107)
              WRITE(IMPLI2,103)
C
              WRITE(IMPLI2,100)
              WRITE(IMPLI2,103)
C
C             --> boucle sur les differents enregistrements
C                et les variables
              REWIND(IMPLI1)
              DO II = 1, NBPDTE
                DO IPP2 = 1, 2*NVLSTA
                  IF ( (IPP2 .LE.NVLSTA
     &                   .AND. IHSLAG(IPP2) .GE. 1 )
     &                .OR.
     &                 (IPP2 .GT.NVLSTA
     &                   .AND. IHSLAG(IPP2-NVLSTA) .EQ. 2 ) ) THEN
C
                  READ(IMPLI1)
     &            JTCABS, XTCABS, (VARCAP(ICAP),ICAP=1,NBCAP(IPP2))
C
                  IF(IPP2.EQ.IPP)
     &              WRITE(IMPLI2,1000)
     &               JTCABS, XTCABS, (VARCAP(ICAP),ICAP=1,NBCAP(IPP))
C
                  ENDIF
                ENDDO
              ENDDO
C
C         --> fermeture fichier
              CLOSE(IMPLI2)
C
            ENDIF
C
          ENDIF
C
        ENDDO
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C 5. EN CAS DE SAUVEGARDE FINALE, DESTRUCTION DU TMP
C=======================================================================
C
CMO      IF(MODHIS.EQ.2) THEN
CMO        NOMFIC = ' '
CMO        NOMFIC = EMPHIS
CMO        CALL VERLON ( NOMFIC,II1,II2,LPOS)
CMOC       ===========
CMO        NOMFIC(II2+1:II2+8) = 'hist.tmp'
CMO        II2 = II2+8
CMO        LPOS = LPOS+8
CMOC
CMO        NENVAR = ' '
CMO        NENVAR = 'rm '
CMO        NENVAR(4:4+LPOS) = NOMFIC ( II1:II2 )
CMO        CALL SYSTEM(NENVAR)
CMOC       ===========
CMO      ENDIF
C
C=======================================================================
C 6. AFFICHAGES
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
C
C----
C FIN
C----
C
      END
c@z
