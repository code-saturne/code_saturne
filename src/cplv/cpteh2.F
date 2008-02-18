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
                        SUBROUTINE CPTEH2
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   ,
     &   RTP    , PROPCE ,
     &   EH0    , EH1     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC CALCUL DE LA TEMPERATURE DES PARTICULES
CFONC   EN FONCTION DE L'ENTHALPIE DU SOLIDE ET DES CONCENTRATIONS
CFONC
c@fonce
C                             ARGUMENTS
c@argub
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT)        !
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! EH0          ! TR ! <-> ! TABLEAU REEL DE TRAVAIL              !
CARGU ! EH1          ! TR ! <-> ! TABLEAU REEL DE TRAVAIL              !
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
C     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C**********************************************************************
C     DONNEES EN COMMON
C**********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
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
      INTEGER          NCELET, NCEL
      DOUBLE PRECISION RTP(NCELET,*), PROPCE(NCELET,*)
      DOUBLE PRECISION EH0(NCELET), EH1(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          I      , ICLA   , ICHA   , ICEL
      INTEGER          IPCTE1 , IPCTE2
      INTEGER          IHFLT2
      DOUBLE PRECISION H2     , X2     , XCH    , XCK
      DOUBLE PRECISION XASH   , XNP    , XTES   , XWAT
C
C***********************************************************************
C RQ IMPORTANTE : On suppose pour l'instant que H2 = H02 + CP2(T2-TREF)
C
      IHFLT2 = 1
C
C=======================================================================
C 1. CALCULS PRELIMINAIRES
C=======================================================================
C
C --- Initialisation des tableaux
C
      DO ICEL = 1, NCEL
        EH0(ICEL) = ZERO
        EH1(ICEL) = ZERO
      ENDDO
C
C --- Initialisation de T2 a T1
C
      IPCTE1 = IPPROC(ITEMP1)
      DO ICLA = 1, NCLACP
        IPCTE2 = IPPROC(ITEMP2(ICLA))
        DO ICEL = 1, NCEL
          PROPCE(ICEL,IPCTE2) = PROPCE(ICEL,IPCTE1)
        ENDDO
      ENDDO
C
C=======================================================================
C 2. CALCUL DE LA TEMPERATURE DES PARTICULES
C=======================================================================
C
      IF ( IHFLT2.EQ.0 ) THEN
C
C --> H2 fonction lineaire de T2
C
        DO ICLA = 1, NCLACP
          IPCTE2 = IPPROC(ITEMP2(ICLA))
          ICHA = ICHCOR(ICLA)
          DO ICEL = 1, NCEL
            PROPCE(ICEL,IPCTE2) =
     &            (RTP(ICEL,ISCA(IH2(ICLA)))-H02CH(ICHA))
     &            / CP2CH(ICHA) + TREFTH
          ENDDO
        ENDDO
C
      ELSE
C
C --> H2 tabule
C
        DO ICLA = 1, NCLACP
C
          IPCTE2 = IPPROC(ITEMP2(ICLA))
C
          I = NPOC-1
          DO ICEL = 1, NCEL
            XCH  = RTP(ICEL,ISCA(IXCH(ICLA)))
            XNP  = RTP(ICEL,ISCA(INP(ICLA)))
            XCK  = RTP(ICEL,ISCA(IXCK(ICLA)))
            XASH = XMASH(ICLA)*XNP
            IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
              XWAT = RTP(ICEL,ISCA(IXWT(ICLA)))
            ELSE
              XWAT = 0.D0
            ENDIF
C
            X2   = XCH + XCK + XASH + XWAT
C
            IF ( X2 .GT. EPSICP*100.D0 ) THEN
C
              H2   = RTP(ICEL,ISCA(IH2(ICLA)))/X2
C
              XTES = XMP0(ICLA)*XNP
C
              IF ( XTES.GT.EPSICP .AND. X2.GT.EPSICP*100.D0 ) THEN
                EH1(ICEL) = XCH /X2 * EHSOLI(ICH(ICHCOR(ICLA) ),I+1)
     &                    + XCK /X2 * EHSOLI(ICK(ICHCOR(ICLA) ),I+1)
     &                    + XASH/X2 * EHSOLI(IASH(ICHCOR(ICLA)),I+1)
     &                    + XWAT/X2 * EHSOLI(IWAT(ICHCOR(ICLA)),I+1)
                IF ( H2.GE.EH1(ICEL) ) THEN
                  PROPCE(ICEL,IPCTE2) = THC(I+1)
C
                ENDIF
              ENDIF
C
            ENDIF
C
          ENDDO
C
          I = 1
          DO ICEL = 1, NCEL
            XCH  = RTP(ICEL,ISCA(IXCH(ICLA)))
            XNP  = RTP(ICEL,ISCA(INP(ICLA)))
            XCK  = RTP(ICEL,ISCA(IXCK(ICLA)))
            XASH = XMASH(ICLA)*XNP
            IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
              XWAT = RTP(ICEL,ISCA(IXWT(ICLA)))
            ELSE
              XWAT = 0.D0
            ENDIF
C
            X2   = XCH + XCK + XASH + XWAT
C
            IF ( X2 .GT. EPSICP*100.D0 ) THEN
C
              H2   = RTP(ICEL,ISCA(IH2(ICLA)))/X2
C
              XTES = XMP0(ICLA)*XNP
C
              IF ( XTES.GT.EPSICP .AND. X2.GT.EPSICP*100.D0 ) THEN
                EH0(ICEL) = XCH /X2 * EHSOLI(ICH(ICHCOR(ICLA) ),I)
     &                    + XCK /X2 * EHSOLI(ICK(ICHCOR(ICLA) ),I)
     &                    + XASH/X2 * EHSOLI(IASH(ICHCOR(ICLA)),I)
     &                    + XWAT/X2 * EHSOLI(IWAT(ICHCOR(ICLA)),I)
                IF ( H2.LE.EH0(ICEL) ) THEN
                  PROPCE(ICEL,IPCTE2) = THC(I)
                ENDIF
              ENDIF
C
            ENDIF
C
          ENDDO
C
          DO I = 1, NPOC-1
            DO ICEL = 1, NCEL
              XCH  = RTP(ICEL,ISCA(IXCH(ICLA)))
              XNP  = RTP(ICEL,ISCA(INP(ICLA)))
              XCK  = RTP(ICEL,ISCA(IXCK(ICLA)))
              XASH = XMASH(ICLA)*XNP
              IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
                XWAT = RTP(ICEL,ISCA(IXWT(ICLA)))
              ELSE
                XWAT = 0.D0
              ENDIF
C
              X2   = XCH + XCK + XASH + XWAT
C
              IF ( X2 .GT. EPSICP*100.D0 ) THEN
C
                H2   = RTP(ICEL,ISCA(IH2(ICLA)))/X2
C
                XTES = XMP0(ICLA)*XNP
C
                IF ( XTES.GT.EPSICP .AND. X2.GT.EPSICP*100.D0 ) THEN
                  EH0(ICEL) = XCH /X2 * EHSOLI(ICH(ICHCOR(ICLA) ),I  )
     &                      + XCK /X2 * EHSOLI(ICK(ICHCOR(ICLA) ),I  )
     &                      + XASH/X2 * EHSOLI(IASH(ICHCOR(ICLA)),I  )
     &                      + XWAT/X2 * EHSOLI(IWAT(ICHCOR(ICLA)),I  )
C
                  EH1(ICEL) = XCH /X2 * EHSOLI(ICH(ICHCOR(ICLA) ),I+1)
     &                      + XCK /X2 * EHSOLI(ICK(ICHCOR(ICLA) ),I+1)
     &                      + XASH/X2 * EHSOLI(IASH(ICHCOR(ICLA)),I+1)
     &                      + XWAT/X2 * EHSOLI(IWAT(ICHCOR(ICLA)),I+1)
C
                  IF ( H2.GE.EH0(ICEL) .AND. H2.LE.EH1(ICEL) ) THEN
                    PROPCE(ICEL,IPCTE2) = THC(I) + (H2-EH0(ICEL)) *
     &                    (THC(I+1)-THC(I))/(EH1(ICEL)-EH0(ICEL))
                  ENDIF
                ENDIF
C
              ENDIF
C
            ENDDO
          ENDDO
C
        ENDDO
C
      ENDIF
C
C----
C FIN
C----
C
      RETURN
      END
c@z
