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
                        SUBROUTINE FUTSVI
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NUMTRA ,
     &   RTP    , PROPCE , VOLUME ,
     &   SMBRS  , ROVSDT ,
     &   W1     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC ROUTINE PHYSIQUE PARTICULIERE : FLAMME FUEL
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
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT PRECEDENT)      !
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! SMBRS(NCELET)! TR ! <-  ! SECOND MEMBRE EXPLICITE              !
CARGU ! ROVSDT(NCELET! TR ! <-  ! PARTIE DIAGONALE IMPLICITE           !
CARGU ! XHLF         ! TR ! <-> ! SOMME DE F1,F3,F4 SUR LE FUEL        !
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
      INCLUDE "fuincl.h"
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
      DOUBLE PRECISION XHLF(NCELET)
      DOUBLE PRECISION W1(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IEL    , IPHAS
      INTEGER          IPCROM
      INTEGER          IXCKCL , IXNPCL , IPCGHT , IPCDIA
C
      DOUBLE PRECISION GAMVAP , GAMHET  , D2S3
      DOUBLE PRECISION FVAP(4), FSD(4)  , FSH(4) , FTRAC(4)
      DOUBLE PRECISION FHET(4)
      DOUBLE PRECISION DIAMHT , DIAMAD
C
      INTEGER          IPCTE1 , IPCTE2
      DOUBLE PRECISION T2MT1
C
C***********************************************************************
C
C=======================================================================
C 1. CALCULS PRELIMINAIRES
C=======================================================================
C
C ---- W1 = X1
C
      DO IEL = 1, NCEL
        W1(IEL) = 1.D0 - RTP(IEL,ISCA(IYFOL))
      ENDDO
C
      D2S3 = 2.D0/3.D0
C
C=======================================================================
C 2. CALCULS DES BILANS EXPLICITE ET IMPLICITE
C=======================================================================
C
C
      IPHAS  = 1
      IPCROM = IPPROC(IROM(IPHAS))
      IXCKCL = ISCA(IFHTF)
      IXNPCL = ISCA(ING)
      IPCGHT = IPPROC(IGMHTF)
      IPCDIA = IPPROC(IDIAM3)
      IPCTE2 = IPPROC(ITEMP3)
      IPCTE1 = IPPROC(ITEMP1)
C
      DO IEL = 1, NCEL
C
        IF ( RTP(IEL,ISCA(ING)) .GT. EPSIFL ) THEN
C
C --> Calculs preliminaires
C
          FTRAC(1) = RTP(IEL,ISCA(IFVAP)) / W1(IEL)
          FTRAC(3) = RTP(IEL,ISCA(IFHTF)) / W1(IEL)
          FTRAC(4) = 1.D0-FTRAC(1)-FTRAC(3)
C
C --> Prise en compte des TS interfaciaux relatifs
C     au phenomene d'evaporation
C
c        GAMVAP = PROPCE(IEL,IPPROC(IGMEVA))*PROPCE(IEL,IPCROM)
C
          T2MT1  = PROPCE(IEL,IPCTE2)-PROPCE(IEL,IPCTE1)
          GAMVAP = -PROPCE(IEL,IPPROC(IGMEVA))*T2MT1
     &                                        *PROPCE(IEL,IPCROM)
C
          FSD(1) = 0.D0
          FSD(3) = 0.D0
          FSD(4) = 0.D0
C
          FVAP(1) = 1.D0
          FVAP(3) = 0.D0
          FVAP(4) = 0.D0
C
          IF ( (FSD(NUMTRA)-FTRAC(NUMTRA))
     &        *(2.D0*FVAP(NUMTRA)-FSD(NUMTRA)-FTRAC(NUMTRA))
     &            .GT. EPSIFL ) THEN
            SMBRS(IEL)  = SMBRS(IEL) + GAMVAP
     &                   *VOLUME(IEL)*(FSD(NUMTRA)-FTRAC(NUMTRA))
     &                   *(2.D0*FVAP(NUMTRA)-FSD(NUMTRA)-FTRAC(NUMTRA))
          ENDIF
C
C --> Prise en compte des TS interfaciaux relatifs
C     au phenomene de combustion heterogene
C
          GAMHET = -PROPCE(IEL,IPCROM)*PROPCE(IEL,IPCGHT)
C
          DIAMHT = ( ( RTP(IEL,ISCA(IYFOL))
     &               /(RTP(IEL,ISCA(ING))*RHO0FL)
     &                -PI*(DINIKF**3)*XINKF/6.D0  )
     &                 *6.D0/(PI*(1.D0-XINKF)) )
          DIAMHT = ABS(DIAMHT)**(1.D0/3.D0)
C
          DIAMAD = PROPCE(IEL,IPCDIA) / DINIKF
C
          FSH(3) = 1.D0
          IF ( DIAMAD.GT.EPSIFL ) THEN
            FSH(3) = 1.D0 - (1.D0-FTRAC(3))
     &              * EXP( PROPCE(IEL,IPCGHT)
     &                    /( 2.D0*PI*2.77D-4*DIAMHT
     &                      *RTP(IEL,IXNPCL)*PROPCE(IEL,IPCROM) ) )
          ENDIF
C
          FSH(1) = FTRAC(1) * (1.D0-FSH(3))/(1.D0-FTRAC(3))
          FSH(4) = FTRAC(4) * (1.D0-FSH(3))/(1.D0-FTRAC(3))
C
          FHET(1) = ZERO
          FHET(3) = 1.D0
          FHET(4) = ZERO
C
          IF (  (FSH(NUMTRA)-FTRAC(NUMTRA))
     &         *(2.D0*FHET(NUMTRA)-FSH(NUMTRA)-FTRAC(NUMTRA))
     &           .GT. EPSIFL ) THEN
            SMBRS(IEL)  = SMBRS(IEL) + GAMHET
     &                   *VOLUME(IEL)*(FSH(NUMTRA)-FTRAC(NUMTRA))
     &                   *(2.D0*FHET(NUMTRA)-FSH(NUMTRA)-FTRAC(NUMTRA))
          ENDIF
C
        ENDIF
C
      ENDDO
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
