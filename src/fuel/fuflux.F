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
                        SUBROUTINE FUFLUX
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NCELET , NCEL   ,
     &   RTPA   , PROPCE , VOLUME ,
     &   W1     , W2     , W3     ,
     &   RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CALCUL DES TERMES DE TRANSFERT DE MASSE ENTRE LA PHASE CONTINUE
CFONC ET LA PHASE DISPERSEE
CFONC
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! RTPA         ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT PRECEDENT)      !
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! W1, W2, W3   ! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU ! RA(*)        ! TR !  -  ! MACRO TABLEAU REEL                   !
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
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
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
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NCELET , NCEL
C
      DOUBLE PRECISION RTPA(NCELET,*), PROPCE(NCELET,*)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION VOLUME(NCELET)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IEL    , IPHAS  , ICLA
      INTEGER          IPCROM , IPCTE1 , IPCTE2 , IPCRO2 , IPCDIA
      INTEGER          IPCGEV , IPCGHT , IPCYOX
      INTEGER          IPCVST,IPCVSL,IPCCP,IPCHGL
C
      DOUBLE PRECISION XNG,XNUSS
      DOUBLE PRECISION PPARO2 , XDFFLI , XDFEXT , XDFTOT0 , XDFTOT1
      DOUBLE PRECISION DIACKA, XUASH
      DOUBLE PRECISION DCOKE , SURF
C
C=======================================================================
C 1. INITIALISATIONS ET CALCULS PRELIMINAIRES
C=======================================================================
C
C --- Initialisation memoire
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C --- Initialisation des termes de transfert de masse
C
      DO ICLA = 1, NCLAFU
        IPCGEV = IPPROC(IGMEVA(ICLA))
        IPCGHT = IPPROC(IGMHTF(ICLA))
        DO IEL = 1, NCEL
          PROPCE(IEL,IPCGEV) = ZERO
          PROPCE(IEL,IPCGHT) = ZERO
        ENDDO
      ENDDO
C
C --- Initialisation des tableaux de travail
C
      DO IEL = 1, NCEL
        W1(IEL) = ZERO
        W2(IEL) = ZERO
        W3(IEL) = ZERO
      ENDDO
C
C --- Calcul de la masse volumique du melange gazeux
C
      IPHAS  = 1
      IPCROM = IPPROC(IROM(IPHAS))
      IPCTE1 = IPPROC(ITEMP1)
      IPCYOX = IPPROC(IYM1(IO2))
C
C --- Numero des grandeurs physiques (voir usclim)
      IPCROM = IPPROC(IROM(IPHAS))
      IPCVST = IPPROC(IVISCT(IPHAS))
C
C --> Terme source pour l'enthalpie du liquide
C
      DO ICLA = 1, NCLAFU
C
        IPCRO2 = IPPROC(IROM3 (ICLA))
        IPCDIA = IPPROC(IDIAM3(ICLA))
        IPCTE2 = IPPROC(ITEMP3(ICLA))
        IPCGHT = IPPROC(IGMHTF(ICLA))
        IPCHGL = IPPROC(IH1HLF(ICLA))
C
C
C ---- Contribution aux bilans explicite et implicite
C        des echanges par diffusion moleculaire
C        6 Lambda Nu / diam**2 / Rho2 * Rho * (T1-T2)
C ------ Calcul de lambda dans W1
C
        XNUSS = 2.D0
        DO IEL = 1, NCEL
          IF ( IVISLS(IHM).GT.0 ) THEN
            IPCVSL = IPPROC(IVISLS(IHM))
            IF ( ICP(IPHAS).GT.0 ) THEN
              IPCCP   = IPPROC(ICP(IPHAS))
              W1(IEL) = PROPCE(IEL,IPCVSL) * PROPCE(IEL,IPCCP)
            ELSE
              W1(IEL) = PROPCE(IEL,IPCVSL) * CP0(IPHAS)
            ENDIF
          ELSE
            IF ( ICP(IPHAS).GT.0 ) THEN
              IPCCP   = IPPROC(ICP(IPHAS))
              W1(IEL) = VISLS0(IHM) * PROPCE(IEL,IPCCP)
            ELSE
              W1(IEL) = VISLS0(IHM) * CP0(IPHAS)
            ENDIF
          ENDIF
        ENDDO
C
C----Contribution aux bilans explicite et implicite des
C    echanges par diffusion moleculaire
C        6 Lambda Nu / diam**2 / Rho2 * Rho
C    le diametre est en mm  donc on multiplie par 1.D-3
C    pour l'avoir en m
C
        DO IEL = 1, NCEL
C
C
          IF ( RTPA(IEL,ISCA(ING(ICLA))) .GT. EPSIFL .AND.
     &           PROPCE(IEL,IPCTE1).GT. PROPCE(IEL,IPCTE2) ) THEN
C
C
c          PROPCE(IEL,IPCHGL) = 6.D0 * W1(IEL) * XNUSS
c     &                         /((PROPCE(IEL,IPCDIA)*1.D-3)**2)
c     &                         /PROPCE(IEL,IPCRO2)
C
            PROPCE(IEL,IPCHGL) = W1(IEL)*XNUSS*RTPA(IEL,ISCA(ING(ICLA)))
     &                          *PI*PROPCE(IEL,IPCDIA)*1.D6
C
          ELSE
            PROPCE(IEL,IPCHGL) =0.D0
          ENDIF
C
        ENDDO
C
      ENDDO
C
C=======================================================================
C 2. TRANSFERTS DE MASSE PAR EVAPORATION
C=======================================================================
C
      DO ICLA = 1, NCLAFU
C
        IPCRO2 = IPPROC(IROM3 (ICLA))
        IPCDIA = IPPROC(IDIAM3(ICLA))
        IPCTE2 = IPPROC(ITEMP3(ICLA))
        IPCGEV = IPPROC(IGMEVA(ICLA))
        IPCHGL = IPPROC(IH1HLF(ICLA))
C
C FO & PPl 16/09/05  Version avec parametres en dur, en attendant les
C                    lectures et inclusion en include
        TEVAP1 = 150.D0 + TKELVI
        TEVAP2 = 450.D0 + TKELVI
C
C
C PPl 161205 On répartit le flux interfacial d'enthalpie entre
C            l'échauffement de la goutte et l'évaporation
C
        DO IEL = 1, NCEL
C
C --- Transfert de masse du a l'evaporation
C
          PROPCE(IEL,IPCGEV) = ZERO
C
C         IF (PROPCE(IEL,IPCTE2)    .GE. TEVAP1 .AND.
C    &       PROPCE(IEL,IPCTE2)    .LE. TEVAP2 .AND.
C    &       RTPA(IEL,ISCA(IYFOL(ICLA))) .GE. EPSIFL      ) THEN
C          PROPCE(IEL,IPCGEV) = PROPCE(IEL,IPPROC(IH1HLF)) /
C     &    ( RTPA(IEL,ISCA(IYFOL)) * CP2FOL * (TEVAP2-TEVAP1)  /
C     &      ( ( RTPA(IEL,ISCA(IYFOL)) + RTPA(IEL,ISCA(IFVAP)) )
C     &    * (1.d0-FKC))                                     + HRFVAP )
C PPl 161205
C
C PPl 090106 On teste plutôt les diamètres car désormais les
C        gouttes finissent trop fines
C
          IF ( PROPCE(IEL,IPCTE2)    .GT. TEVAP1 .AND.
     &         PROPCE(IEL,IPCDIA)    .GT. DINIKF(ICLA) .AND.
     &         RTPA(IEL,ISCA(IYFOL(ICLA))) .GT. EPSIFL       ) THEN
            PROPCE(IEL,IPCGEV) = PROPCE(IEL,IPCHGL)
     &                        /( HRFVAP + CP2FOL*(TEVAP2-TEVAP1) )
          ENDIF
C
        ENDDO
C
      ENDDO
C
C=======================================================================
C 3. CALCUL DE RHO_COKE MOYEN
C    On suppose pour le calcul de la masse volumique du coke que
C    l'evaporation a lieu a volume de coke constant
C    Par la suite, on suppose (pour commencer ?) que la
C    combustion hétérogène a lieu à volume constant =>
C    à masse volumique décroissante
C=======================================================================
C
C --- Initialisation
C
      RHOKF  = RHO0FL
C
C --- Calcul de la masse volumique moyenne du coke
C
C=======================================================================
C 4. TRANSFERTS DE MASSE PAR COMBUSTION HETEROGENE
C=======================================================================
C
      DO ICLA = 1, NCLAFU
C
        IPCRO2 = IPPROC(IROM3 (ICLA))
        IPCDIA = IPPROC(IDIAM3(ICLA))
        IPCTE2 = IPPROC(ITEMP3(ICLA))
        IPCGEV = IPPROC(IGMEVA(ICLA))
        IPCGHT = IPPROC(IGMHTF(ICLA))
C
        DO IEL = 1, NCEL
C
          IF ( PROPCE(IEL,IPCDIA) .LE. DINIKF(ICLA) .AND.
     &         PROPCE(IEL,IPCDIA) .GT. DINIIN(ICLA) .AND.
     &         RTPA(IEL,ISCA(IYFOL(ICLA))) .GT. EPSIFL ) THEN
C
            XNG   = RTPA(IEL,ISCA(ING(ICLA)))*1.D9
C
            XUASH = XNG*(1.D0-XINFOL)
     &             *(RHO0FL*PI*((DINIFL(ICLA)*1.D-3)**2))/6.D0
C
C --- Calcul de la pression partielle en oxygene (atm)
C                                                 ---
C       PO2 = RHO1*RR*T*YO2/MO2
C
            PPARO2 = PROPCE(IEL,IPPROC(IROM1))*RR*PROPCE(IEL,IPCTE1)
     &              *PROPCE(IEL,IPCYOX)/WMOLE(IO2)
            PPARO2 = PPARO2 / PREFTH
C
C --- Calcul de Dcoke en metres
C
            DCOKE = ( ( RTPA(IEL,ISCA(IYFOL(ICLA)))
     &              /(RTPA(IEL,ISCA(ING(ICLA)))*RHO0FL)
     &              -PI*(DINIKF(ICLA)**3)*XINKF/6.D0  )
     &               *6.D0/(PI*(1.D0-XINKF)) )**(1.D0/3.D0)
     &            *1.D-3
           IF ( DCOKE .LT. 0.D0 ) THEN
             WRITE(NFECRA,*) 'erreur Dcoke = ',Dcoke,IEL
             CALL CSEXIT(1)
           ENDIF
C
C --- Coefficient de cinetique chimique de formation de CO
C       en (kg.m-2.s-1.atm(-n))
C
            XDFFLI = AHETFL*EXP(-EHETFL*4185.D0
     &              /(RR*PROPCE(IEL,IPCTE1)))
C
C --- Coefficient de diffusion en  (Kg/m2/s/atm) : XDFEXT
C     Coefficient global pour n=0.5 en (kg/m2/s) : XDFTOT0
C     Coefficient global pour n=1   en (Kg/m2/s) : XDFTOT1
C
            DIACKA = DCOKE/(DINIKF(ICLA)*1.D-3)
            IF ( DIACKA .GT. EPSIFL ) THEN
              XDFEXT = 2.53D-7*((PROPCE(IEL,IPCTE1))**0.75D0)
     &                / DCOKE*2.D0
              XDFTOT1 = PPARO2 / ( 1.D0/XDFFLI + 1.D0/XDFEXT )
              XDFTOT0 = -(XDFFLI**2)/(2.D0*XDFEXT**2)+(PPARO2*XDFFLI**2
     &                  +(XDFFLI**4)/(2.D0*XDFEXT**2))**0.5D0
            ELSE
              XDFTOT1 = XDFFLI*PPARO2
              XDFTOT0 = XDFFLI*PPARO2**0.5D0
            ENDIF
C
C     Surface
C
            SURF = PI*(DCOKE**2)
C
C --- Calcul de PROPCE(IEL,IPCGHT) = - COXCK*XDFTOT0*PPARO2*XNP < 0
C --- ou        PROPCE(IEL,IPCGHT) = - COXCK*XDFTOT1*PPARO2*XNP < 0
C
            IF (IOFHET.EQ.1) THEN
c             PROPCE(IEL,IPCGHT) = - XDFTOT1*COXCK*XNG
              PROPCE(IEL,IPCGHT) = - XDFTOT1*SURF*XNG
            ELSE
c             PROPCE(IEL,IPCGHT) = - XDFTOT0*COXCK*XNG
              PROPCE(IEL,IPCGHT) = - XDFTOT0*SURF*XNG
            ENDIF
C
          ELSE
            PROPCE(IEL,IPCGHT) = 0.D0
          ENDIF
C
        ENDDO
C
      ENDDO
C
C=======================================================================
C FORMATS
C--------
C----
C FIN
C----
C
      RETURN
      END
c@z
