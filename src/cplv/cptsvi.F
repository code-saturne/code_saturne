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
                        SUBROUTINE CPTSVI
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NUMTRA ,
     &   RTP    , PROPCE , VOLUME ,
     &   SMBRS , ROVSDT ,
     &   XF1M   , XF2M   ,
     &   W1     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC ROUTINE PHYSIQUE PARTICULIERE : FLAMME CP
CFONC   TERMES DE TRANSFERT DE MASSE ENTRE LA PHASE CONTINUE
CFONC   ET LA PHASE DISPERSEE RELATIF A LA VARIANCE DES DIFFERENTS
CFONC   TRACEURS (BILANS EXPLICITE ET IMPLICITE)
CFONC
CFONC   ATTENTION CETTE SUBROUTINE N'EST OPERATIONNELLE QUE
CFONC               POUR NUMTRA = 3 ou 4
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NUMTRA       ! E  !  -> ! NUMERO DU TRACEUR CONCERNE (1,2,3,4) !
CARGU ! RTP         ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT PRECEDENT)      !
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! SMBRS(NCELET)! TR ! <-  ! SECOND MEMBRE EXPLICITE              !
CARGU ! ROVSDT(NCELET! TR ! <-  ! PARTIE DIAGONALE IMPLICITE           !
CARGU ! XF1M, XF2M   ! TR ! <-> ! SOMME DE F1,F2,F3,F4 SUR L'ENSEMBLE  !
CARGU !              !    !     ! SUR L'ENSEMBLE DES CHARBONS ET       !
CARGU !              !    !     ! DES CLASSES                          !
CARGU ! W1(NCELET)   ! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
C-----------------------------------------------------------------------
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
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET , NCEL   , NUMTRA
C
      DOUBLE PRECISION RTP(NCELET,*), PROPCE(NCELET,*)
      DOUBLE PRECISION VOLUME(NCELET)
      DOUBLE PRECISION SMBRS(NCELET), ROVSDT(NCELET)
      DOUBLE PRECISION XF1M(NCELET)  , XF2M(NCELET)
      DOUBLE PRECISION W1(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IEL    , ICLA   , ICHA
      INTEGER          IPHAS  , IPCROM , IXCHCL , IXCKCL , IXNPCL
      INTEGER          IPCGD1 , IPCGD2 , IPCGHT , IPCDIA
C
      DOUBLE PRECISION GAMDEV1 , GAMDEV2 , GMDEVT , GAMHET  , D2S3
      DOUBLE PRECISION FDEV(4) , FHET(4) , FSD(4)  , FSH(4) , FTRAC(4)
      DOUBLE PRECISION DIAMDV , DIAMHT , X2 , XMP , DIAMAD
C
C***********************************************************************
C
C=======================================================================
C 1. CALCULS PRELIMINAIRES
C=======================================================================
C
C --> Calcul de X1 = 1. - SOMME(X2) dans W1
C
      DO IEL = 1, NCEL
        W1(IEL) = 1.D0
      ENDDO
C
      DO ICLA = 1, NCLACP
        IXCHCL = ISCA(IXCH(ICLA))
        IXCKCL = ISCA(IXCK(ICLA))
        IXNPCL = ISCA(INP(ICLA ))
        DO IEL = 1, NCEL
          X2 = RTP(IEL,IXCHCL)+RTP(IEL,IXCKCL)
     &       + RTP(IEL,IXNPCL)*XMASH(ICLA)
          W1(IEL) = W1(IEL) - X2
        ENDDO
      ENDDO
C
      IPHAS = 1
      IPCROM = IPPROC(IROM(IPHAS))
C
C Calcul de F1 et F2 sur l'ensemble des charbons
C
      DO IEL = 1, NCEL
        XF1M(IEL) = 0.D0
        XF2M(IEL) = 0.D0
      ENDDO
      DO ICHA = 1, NCHARB
        DO IEL = 1, NCEL
          XF1M(IEL) = XF1M(IEL) + RTP(IEL,ISCA(IF1M(ICHA)))
          XF2M(IEL) = XF2M(IEL) + RTP(IEL,ISCA(IF2M(ICHA)))
        ENDDO
      ENDDO
C
      D2S3 = 2.D0/3.D0
C
C=======================================================================
C 2. CALCULS DES BILANS EXPLICITE ET IMPLICITE
C=======================================================================
C
      DO ICLA = 1, NCLACP
C
        IXCHCL = ISCA(IXCH(ICLA))
        IXCKCL = ISCA(IXCK(ICLA))
        IXNPCL = ISCA(INP(ICLA ))
        IPCGD1 = IPPROC(IGMDV1(ICLA))
        IPCGD2 = IPPROC(IGMDV2(ICLA))
        IPCGHT = IPPROC(IGMHET(ICLA))
        IPCDIA = IPPROC(IDIAM2(ICLA))
C
        DO IEL = 1, NCEL
C
C --> Calculs preliminaires
C
          FTRAC(1) = XF1M(IEL)           / W1(IEL)
          FTRAC(2) = XF2M(IEL)           / W1(IEL)
          FTRAC(3) = RTP(IEL,ISCA(IF3M)) / W1(IEL)
          FTRAC(4) = 1.D0-FTRAC(1)-FTRAC(2)-FTRAC(3)
C
          XMP = XMP0(ICLA)*RTP(IEL,IXNPCL)
          X2  = RTP(IEL,IXCHCL) + RTP(IEL,IXCKCL) +
     &          RTP(IEL,IXNPCL)*XMASH(ICLA)
C
          IF ( XMP .GT. EPSICP .AND.
     &         X2  .GT. EPSICP        ) THEN
C
C --> Prise en compte des TS interfaciaux relatifs
C     au phenomene de devolatilisation
C
            IF (RTP(IEL,IXCHCL).GT.EPSICP) THEN
C
C     On prend comme diametre le diametre initial de la particule
C     la devolatilisation ayant lieu a diametre constant
C
              DIAMDV = DIAM20(ICLA)
              FSD(1) = 1.D0 - (1.D0-FTRAC(1))
     &             * EXP( ( RTP(IEL,IXCHCL)*PROPCE(IEL,IPCGD1)  )
     &                   /( 2.D0*PI*2.77D-4*DIAMDV
     &                      *RTP(IEL,IXNPCL)*PROPCE(IEL,IPCROM) ) )
C
              FSD(2) = 1.D0 - (1.D0-FTRAC(2))
     &             * EXP( ( RTP(IEL,IXCHCL)*PROPCE(IEL,IPCGD2)  )
     &                   /( 2.D0*PI*2.77D-4*DIAMDV
     &                      *RTP(IEL,IXNPCL)*PROPCE(IEL,IPCROM) ) )
C
              FSD(3) = FTRAC(3) * (1.D0 -   FSD(1) -   FSD(2))
     &                        / (1.D0 - FTRAC(1) - FTRAC(2))
              FSD(4) = FTRAC(4) * (1.D0 -   FSD(1) -   FSD(2))
     &                        / (1.D0 - FTRAC(1) - FTRAC(2))
C
              GAMDEV1 = - PROPCE(IEL,IPCROM) * RTP(IEL,IXCHCL)
     &                  * PROPCE(IEL,IPCGD1)
              GAMDEV2 = - PROPCE(IEL,IPCROM) * RTP(IEL,IXCHCL)
     &                  * PROPCE(IEL,IPCGD2)
              GMDEVT = GAMDEV1+GAMDEV2
C
              FDEV(1) = GAMDEV1 / GMDEVT
              FDEV(2) = GAMDEV2 / GMDEVT
              FDEV(3) = ZERO
              FDEV(4) = ZERO
C
              IF (  (FSD(NUMTRA)-FTRAC(NUMTRA))
     &           *(2.D0*FDEV(NUMTRA)-FSD(NUMTRA)-FTRAC(NUMTRA))
     &          .GT. EPSICP ) THEN
                SMBRS(IEL) = SMBRS(IEL) + GMDEVT
     &                    *VOLUME(IEL)*(FSD(NUMTRA)-FTRAC(NUMTRA))
     &                    *(2.D0*FDEV(NUMTRA)-FSD(NUMTRA)-FTRAC(NUMTRA))
c                ROVSDT(IEL) = ROVSDT(IEL) + ZERO
              ENDIF
C
            ENDIF
C
C --> Prise en compte des TS interfaciaux relatifs
C     au phenomene de combustion heterogene
C
C     On prend comme diametre le diametre du coke
C
            DIAMHT = PROPCE(IEL,IPCDIA)
            DIAMAD = PROPCE(IEL,IPCDIA) / DIAM20(ICLA)
            FSH(3) = 1.D0
            IF ( DIAMAD.GT.EPSICP ) THEN
              FSH(3) = 1.D0 - (1.D0-FTRAC(3))
     &               * EXP( ( (RTP(IEL,IXCKCL)**D2S3)
     &                        *PROPCE(IEL,IPCGHT)       )
     &                     /( 2.D0*PI*2.77D-4*DIAMHT
     &                        *RTP(IEL,IXNPCL)*PROPCE(IEL,IPCROM) ) )
            ENDIF
C
            FSH(1) = FTRAC(1) * (1.D0-FSH(3))/(1.D0-FTRAC(3))
            FSH(2) = FTRAC(2) * (1.D0-FSH(3))/(1.D0-FTRAC(3))
            FSH(4) = FTRAC(4) * (1.D0-FSH(3))/(1.D0-FTRAC(3))
C
            GAMHET = - PROPCE(IEL,IPCROM) * (RTP(IEL,IXCKCL)**D2S3)
     &                 * PROPCE(IEL,IPCGHT)
C
            FHET(1) = ZERO
            FHET(2) = ZERO
            FHET(3) = 1.D0
            FHET(4) = ZERO
C
            IF (  (FSH(NUMTRA)-FTRAC(NUMTRA))
     &           *(2.D0*FHET(NUMTRA)-FSH(NUMTRA)-FTRAC(NUMTRA))
     &         .GT. EPSICP ) THEN
              SMBRS(IEL) = SMBRS(IEL) + GAMHET
     &                    *VOLUME(IEL)*(FSH(NUMTRA)-FTRAC(NUMTRA))
     &                    *(2.D0*FHET(NUMTRA)-FSH(NUMTRA)-FTRAC(NUMTRA))
c              ROVSDT(IEL) = ROVSDT(IEL) + ZERO
            ENDIF
C
          ENDIF
C
        ENDDO
C
      ENDDO
C
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
