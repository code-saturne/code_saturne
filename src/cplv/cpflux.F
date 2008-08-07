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
                        SUBROUTINE CPFLUX
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
      INCLUDE "ppincl.h"
      INCLUDE "ppcpfu.h"
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
      INTEGER          IEL    , IPHAS  , ICHA   , ICLA
      INTEGER          IPCROM , IPCTE1 , IPCTEM , IPCRO2 , IPCDIA
      INTEGER          IPCGD1 , IPCGD2 , IPCGCH , IPCGHT , IPCYOX
      INTEGER          IPCSEC
      INTEGER          IPCVSL , IPCCP, IROMF , IPCTE2
      INTEGER          NPOIN1,NPOIN2,NPOIN3,NPOIN4,NPOIN5,NPOIN6
      INTEGER          NPOIN63, NPOIN65
      INTEGER          NPYV, MODNTL
C
      DOUBLE PRECISION X2     , XCH    , XCK    , XASH   , XNP , XUASH
      DOUBLE PRECISION PPARO2 , XDFCHI , XDFEXT , XDFTOT0 , XDFTOT1
      DOUBLE PRECISION DEVTO1(NCHARM)  , DEVTO2(NCHARM) , COXCK ,  DEN
      DOUBLE PRECISION DIACKA
      DOUBLE PRECISION DP , LV, YVS, YV , TEBL , SHRD , XMEAU, XMGAZ
      DOUBLE PRECISION XNUSS, TLIMIT , TMINI
      DOUBLE PRECISION TMIN, TMAX, YVMIN, YVMAX, YYMAX
C
      INTEGER          IPYCO2
      DOUBLE PRECISION PPRCO2
C
C***********************************************************************
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
      DO ICLA = 1, NCLACP
        IPCGD1 = IPPROC(IGMDV1(ICLA))
        IPCGD2 = IPPROC(IGMDV2(ICLA))
        IPCGCH = IPPROC(IGMDCH(ICLA))
        IPCGHT = IPPROC(IGMHET(ICLA))
        DO IEL = 1, NCEL
          PROPCE(IEL,IPCGD1) = ZERO
          PROPCE(IEL,IPCGD2) = ZERO
          PROPCE(IEL,IPCGCH) = ZERO
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
C
C ---- W1 = Somme (X2i)
C      W2 = Somme (X2i/Rho2i)
C
      DO ICLA = 1, NCLACP
        IPCRO2 = IPPROC(IROM2(ICLA))
        DO IEL = 1, NCEL
          XCK  = RTPA(IEL,ISCA(IXCK(ICLA)))
          XCH  = RTPA(IEL,ISCA(IXCH(ICLA)))
          XASH = RTPA(IEL,ISCA(INP (ICLA)))*XMASH(ICLA)
          X2   = XCH + XCK + XASH
C
C         Prise en compte de l'humidite
C
          IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
            X2 = X2+RTPA(IEL,ISCA(IXWT(ICLA)))
          ENDIF
C
          W1(IEL) = W1(IEL) + X2
          W2(IEL) = W2(IEL) + X2 / PROPCE(IEL,IPCRO2)
        ENDDO
      ENDDO
C
C ---- W3 : Rho 1
C
      DO IEL = 1, NCEL
        W3(IEL) = (1.D0-W1(IEL)) / (1.D0/PROPCE(IEL,IPCROM)-W2(IEL))
      ENDDO
C
C
C
C=======================================================================
C 2. TRANSFERTS DE MASSE PAR DEVOLATILISATION
C=======================================================================
C
      IPCTE1 = IPPROC(ITEMP1)
C
      DO ICLA = 1, NCLACP
C
        IPCGD1 = IPPROC(IGMDV1(ICLA))
        IPCGD2 = IPPROC(IGMDV2(ICLA))
        IPCGCH = IPPROC(IGMDCH(ICLA))
        IPCTEM = IPPROC(ITEMP2(ICLA))
C
        DO IEL = 1, NCEL
C
C --- Transfert de masse du au degagemnt des MV legeres (s-1) < 0
C
          PROPCE(IEL,IPCGD1) = -Y1CH(ICHCOR(ICLA))*A1CH(ICHCOR(ICLA))
     &      * EXP(-E1CH(ICHCOR(ICLA))/(RR*PROPCE(IEL,IPCTEM)))
C
C --- Transfert de masse du au degagemnt des MV lourdes (s-1) < 0
C
          PROPCE(IEL,IPCGD2) = -Y2CH(ICHCOR(ICLA))*A2CH(ICHCOR(ICLA))
     &      * EXP(-E2CH(ICHCOR(ICLA))/(RR*PROPCE(IEL,IPCTEM)))
C
C --- Taux de disparition du charbon reactif (s-1) < 0
C
          PROPCE(IEL,IPCGCH) = - A1CH(ICHCOR(ICLA))
     &      * EXP(-E1CH(ICHCOR(ICLA))/(RR*PROPCE(IEL,IPCTEM)))
     &                         - A2CH(ICHCOR(ICLA))
     &      * EXP(-E2CH(ICHCOR(ICLA))/(RR*PROPCE(IEL,IPCTEM)))
C
        ENDDO
C
      ENDDO
C
C
C=======================================================================
C 3. CALCUL DE RHO_COKE MOYEN POUR CHAQUE CHARBON
C    On suppose pour le calcul de la masse volumique du coke que
C    la devolatilisation a lieu a volume constant
C=======================================================================
C
C --- Initialisation
C
      DO ICHA = 1, NCHARB
        DEVTO1(ICHA) = ZERO
        DEVTO2(ICHA) = ZERO
        RHOCK(ICHA)  = RHO0CH(ICHA)
      ENDDO
C
C --- Calcul de l'integrale de GMDEV1 et GMDEV2 pour chaque charbon
C
      IPHAS = 1
      IPCROM = IPPROC(IROM(IPHAS))
      DO ICLA = 1, NCLACP
        IPCGD1 = IPPROC(IGMDV1(ICLA))
        IPCGD2 = IPPROC(IGMDV2(ICLA))
        DO IEL = 1, NCEL
          XCH = RTPA(IEL,ISCA(IXCH(ICLA)))
          DEVTO1(ICHCOR(ICLA)) = DEVTO1(ICHCOR(ICLA)) -
     &      ( PROPCE(IEL,IPCGD1)*XCH*PROPCE(IEL,IPCROM)
     &        *VOLUME(IEL) )
          DEVTO2(ICHCOR(ICLA)) = DEVTO2(ICHCOR(ICLA)) -
     &      ( PROPCE(IEL,IPCGD2)*XCH*PROPCE(IEL,IPCROM)
     &        *VOLUME(IEL) )
        ENDDO
        IF(IRANGP.GE.0) THEN
          CALL PARSOM(DEVTO1(ICHCOR(ICLA)))
          CALL PARSOM(DEVTO2(ICHCOR(ICLA)))
        ENDIF
      ENDDO
C
C --- Calcul de la masse volumique moyenne du coke
C
      DO ICHA = 1, NCHARB
        DEN = Y2CH(ICHA)*DEVTO1(ICHA)+Y1CH(ICHA)*DEVTO2(ICHA)
        IF ( DEN.GT.EPSICP ) THEN
          RHOCK(ICHA) = RHO0CH(ICHA) * ( 1.D0 -
     &     ( Y1CH(ICHA)*Y2CH(ICHA)*(DEVTO1(ICHA)+DEVTO2(ICHA))/DEN ) )
        ENDIF
      ENDDO
C
C=======================================================================
C 4. TRANSFERTS DE MASSE PAR COMBUSTION HETEROGENE AVEC O2
C=======================================================================
C
      IPCYOX = IPPROC(IYM1(IO2))
      IPCTE1 = IPPROC(ITEMP1)
C
      DO ICLA = 1, NCLACP
C
        IPCGHT = IPPROC(IGMHET(ICLA))
        IPCDIA = IPPROC(IDIAM2(ICLA))
        ICHA   = ICHCOR(ICLA)
        IPCTEM = IPPROC(ITEMP2(ICLA))
C
        DO IEL = 1, NCEL
C
          XNP   = RTPA(IEL,ISCA(INP(ICLA)))
          XUASH = XNP*(1.D0-XASHCH(ICHA))*XMP0(ICLA)
C
C --- Calcul de la pression partielle en oxygene (atm)
C                                                 ---
C       PO2 = RHO1*RR*T*YO2/MO2
C
          PPARO2 = W3(IEL)*RR*PROPCE(IEL,IPCTE1)
     &           * PROPCE(IEL,IPCYOX)/WMOLE(IO2)
          PPARO2 = PPARO2 / PREFTH
C
C --- Coefficient de cinetique chimique de formation de CO
C       en (kg.m-2.s-1.atm(-n))
C
          XDFCHI = AHETCH(ICHCOR(ICLA))
     &      * EXP(-EHETCH(ICHCOR(ICLA))*4185.D0
     &             / (RR * PROPCE(IEL,IPCTEM)) )
C
C --- Coefficient de diffusion en  (Kg/m2/s/atm) : XDFEXT
C     Coefficient global pour n=0.5 en (kg/m2/s) : XDFTOT0
C     Coefficient global pour n=1   en (Kg/m2/s) : XDFTOT1
C
          DIACKA = PROPCE(IEL,IPCDIA)/DIAM20(ICLA)
          IF ( DIACKA .GT. EPSICP ) THEN
            XDFEXT = 2.53D-7*((PROPCE(IEL,IPCTEM))**0.75D0)
     &              / PROPCE(IEL,IPCDIA)*2.D0
            XDFTOT1 = PPARO2 / ( 1.D0/XDFCHI + 1.D0/XDFEXT )
            XDFTOT0 = -(XDFCHI**2)/(2.D0*XDFEXT)+(PPARO2*XDFCHI**2
     &                + (XDFCHI**4)/(4.D0*XDFEXT**2))**0.5D0
          ELSE
            XDFTOT1 = XDFCHI*PPARO2
            XDFTOT0 = XDFCHI*PPARO2**0.5D0
          ENDIF
C
C Rq AE : Pour l'instant, on se limite a ce test
C         L'introduction d'une correlation de soufflage viendra
C         ulterieurement.
C
C --- Calcul de COXCK tq : Se = COXCK * XCK**(2/3)
C
          COXCK = ZERO
          IF ( XUASH.GT.EPSICP ) THEN
            COXCK = PI*(DIAM20(ICLA)**2)*
     &            (RHO20(ICLA)/(RHOCK(ICHA)*XUASH))**(2.D0/3.D0)
          ENDIF
C
C --- Calcul de PROPCE(IEL,IPCGHT) = - COXCK*XDFTOT0*PPARO2*XNP < 0
C --- ou  PROPCE(IEL,IPCGHT) = - COXCK*XDFTOT1*PPARO2*XNP < 0
C
          IF (IOCHET(ICHA).EQ.1) THEN
            PROPCE(IEL,IPCGHT) = - XDFTOT1*COXCK*XNP
          ELSE
            PROPCE(IEL,IPCGHT) = - XDFTOT0*COXCK*XNP
          ENDIF
C
        ENDDO
C
      ENDDO
C
C=======================================================================
C 5. TRANSFERTS DE MASSE PAR COMBUSTION HETEROGENE AVEC CO2
C=======================================================================
C
      IF ( IHTCO2 .EQ. 1) THEN
C
        IPYCO2 = IPPROC(IYM1(ICO2))
        IPCTE1 = IPPROC(ITEMP1)
C
        DO ICLA = 1, NCLACP
C
          IPCGHT = IPPROC(IGHCO2(ICLA))
          IPCDIA = IPPROC(IDIAM2(ICLA))
          ICHA   = ICHCOR(ICLA)
          IPCTEM = IPPROC(ITEMP2(ICLA))
C
          DO IEL = 1, NCEL
C
            XNP   = RTPA(IEL,ISCA(INP(ICLA)))
            XUASH = XNP*(1.D0-XASHCH(ICHA))*XMP0(ICLA)
C
C --- Calcul de la pression partielle en oxygene (atm)
C                                                 ---
C       PCO2 = RHO1*RR*T*YO2/MO2
C
            PPRCO2 = W3(IEL)*RR*PROPCE(IEL,IPCTE1)
     &              *PROPCE(IEL,IPYCO2)/WMOLE(ICO2)
            PPRCO2 = PPRCO2 / PREFTH
C
C --- Coefficient de cinetique chimique de formation de CO
C       en (kg.m-2.s-1.atm(-n))
C
            XDFCHI = AHETC2(ICHCOR(ICLA))
     &        * EXP(-EHETC2(ICHCOR(ICLA))*4185.D0
     &               / (RR * PROPCE(IEL,IPCTEM)) )
C
C --- Coefficient de diffusion en  (Kg/m2/s/atm) : XDFEXT
C     Coefficient global pour n=0.5 en (kg/m2/s) : XDFTOT0
C     Coefficient global pour n=1   en (Kg/m2/s) : XDFTOT1
C
            DIACKA = PROPCE(IEL,IPCDIA)/DIAM20(ICLA)
            IF ( DIACKA .GT. EPSICP ) THEN
              XDFEXT = 2.53D-7*((PROPCE(IEL,IPCTEM))**0.75D0)
     &                / PROPCE(IEL,IPCDIA)*2.D0
              XDFTOT1 = PPRCO2 / ( 1.D0/XDFCHI + 1.D0/XDFEXT )
              XDFTOT0 = -(XDFCHI**2)/(2.D0*XDFEXT)+(PPRCO2*XDFCHI**2
     &                  +(XDFCHI**4)/(4.D0*XDFEXT**2))**0.5D0
            ELSE
              XDFTOT1 = XDFCHI*PPRCO2
              XDFTOT0 = XDFCHI*PPRCO2**0.5D0
            ENDIF
C
C Rq AE : Pour l'instant, on se limite a ce test
C         L'introduction d'une correlation de soufflage viendra
C         ulterieurement.
C
C --- Calcul de COXCK tq : Se = COXCK * XCK**(2/3)
C
            COXCK = ZERO
            IF ( XUASH.GT.EPSICP ) THEN
              COXCK = PI*(DIAM20(ICLA)**2)*
     &              (RHO20(ICLA)/(RHOCK(ICHA)*XUASH))**(2.D0/3.D0)
            ENDIF
C
C --- Calcul de PROPCE(IEL,IPCGHT) = - COXCK*XDFTOT0*PPRCO2*XNP < 0
C --- ou  PROPCE(IEL,IPCGHT) = - COXCK*XDFTOT1*PPRCO2*XNP < 0
C
            IF (IOETC2(ICHA).EQ.1) THEN
              PROPCE(IEL,IPCGHT) = - XDFTOT1*COXCK*XNP
            ELSE
              PROPCE(IEL,IPCGHT) = - XDFTOT0*COXCK*XNP
            ENDIF
C
          ENDDO
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C 6. TRANSFERTS DE MASSE LORS DE LA PHASE DE SECHAGE
C=======================================================================
C
      IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
C
C      Chaleur Latente en J/kg
        LV     = 2.263D+6
        TEBL   = 100.D0+TKELVI
C
        XNUSS = 2.D0
        SHRD  = 2.D0
        XMEAU = 0.018D0
C
        TLIMIT = 302.24D0
        TMINI   = TLIMIT*(1.D0-TLIMIT/(LV*XMEAU))
        write(NFECRA,*) ' TMIN = ',TMINI
C
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
        IF(NTLIST.GT.0) THEN
          MODNTL = MOD(NTCABS,NTLIST)
        ELSEIF(NTLIST.EQ.-1.AND.NTCABS.EQ.NTMABS) THEN
          MODNTL = 0
        ELSE
          MODNTL = 1
        ENDIF
C
        DO ICLA = 1, NCLACP
C
          NPOIN1 = 0
          NPOIN2 = 0
          NPOIN3 = 0
          NPOIN4 = 0
          NPOIN5 = 0
          NPOIN6 = 0
          NPOIN63 = 0
          NPOIN65 = 0
C
          ICHA   = ICHCOR(ICLA)
C
          IROMF  = IPPROC(IROM1)
          IPCTE1 = IPPROC(ITEMP1)
          IPCTE2 = IPPROC(ITEMP2(ICLA))
          IPCSEC = IPPROC(IGMSEC(ICLA))
          IPCRO2 = IPPROC(IROM2(ICLA))
C
C -------- Calcul du diametre des particules dans W2
C          On calcule le d20 = (A0.D0**2+(1-A0)*DCK**2)**0.5
C
C
          TMAX = -1.D+20
          TMIN =  1.D+20
          YVMIN = 1.D+20
          YVMAX =-1.D+20
C
          NPYV  = 0
          YYMAX = 0.D0
C
          DO IEL = 1, NCEL
C
            PROPCE(IEL,IPCSEC) = 0.D0
C
            IF ( RTPA(IEL,ISCA(IXWT(ICLA))) .GT. EPSICP  ) THEN
C
              NPOIN1 = NPOIN1 + 1
C
              XNP = RTPA(IEL,ISCA(INP(ICLA)))
C
C             Calcul du diametre des particules dans W2
C             On calcule le d20 = (A0.D0**2+(1-A0)*DCK**2)**0.5
C
              DP  = ( XASHCH(ICHA)*DIAM20(ICLA)**2 +
     &                 (1.D0-XASHCH(ICHA))*PROPCE(IEL,IPCDIA)**2
     &                  )**0.5
C
c              IF ( PROPCE(IEL,IPCTE2).LT. TEBL ) THEN
C
                NPOIN2 = NPOIN2 + 1
C
                IF ( PROPCE(IEL,IPCTE2) .GT. TLIMIT )THEN
                  XMGAZ = PROPCE(IEL,IPPROC(IMMEL))
                  YVS   = XMEAU/XMGAZ
     &                   *EXP( LV*XMEAU
     &                        *(1.D0/TEBL-1.D0/PROPCE(IEL,IPCTE2))
     &                        /RR )
C
                ELSE
                  XMGAZ = PROPCE(IEL,IPPROC(IMMEL))
                  YVS   = XMEAU/XMGAZ
     &                   *EXP( LV*XMEAU
     &                        *(1.D0/TEBL-1.D0/TLIMIT)
     &                        /RR )
     &                   *(LV*XMEAU*(PROPCE(IEL,IPCTE2)-TMINI))
     &                   /(TLIMIT*TLIMIT)
                ENDIF
C
C                YV    = 2.D0*YVS/3.D0
C
                YV = YVS+ (1.D0/3.D0)
     &                   *(PROPCE(IEL,IPPROC(IYM1(IH2O)))-YVS)
C
               IF ( YV .LT. 0.D0 ) THEN
                 write(NFECRA,*) YV,YVS,PROPCE(IEL,IPPROC(IYM1(IH2O)))
                 write(NFECRA,*) PROPCE(IEL,IPCTE2),TMINI
                 call csexit(1)
               ENDIF
               IF ( YV .GE. 1.D0) THEN
                 YYMAX = MAX(YYMAX,YV)
                 NPYV  = NPYV +1
                 YV    = 0.99D0
               ENDIF
C
                YVMIN=MIN(YVMIN,YV)
                YVMAX=MAX(YVMAX,YV)
C
                IF ( YV .GT. EPSICP .AND. YV .LT. 1.D0) THEN
C
                  NPOIN3 = NPOIN3 + 1
C
                  PROPCE(IEL,IPCSEC) = PI*DP*PROPCE(IEL,IROMF)
     &                                *DIFTL0*SHRD*XNP
     &                                *LOG(1.D0/(1.D0-YV))
                  if ( PROPCE(IEL,IPCSEC) .LT. 0.D0 ) THEN
                    PROPCE(IEL,IPCSEC) = 0.D0
                    NPOIN63 = NPOIN63 + 1
                  endif
                ELSE
C
                  NPOIN4 = NPOIN4 + 1
C
                  PROPCE(IEL,IPCSEC) = 0.D0
                ENDIF
C
c              ELSE
C
c                NPOIN5 = NPOIN5 + 1
C
c                PROPCE(IEL,IPCSEC) = W1(IEL)*XNUSS*XNP
c     &                              *PI*DP
c     &                              *( PROPCE(IEL,IPCTE1)
c     &                                -PROPCE(IEL,IPCTE2))/LV
C
C
c                PROPCE(IEL,IPCSEC) = 6.D0*W1(IEL)*XNUSS
c     &                              /(DP**2)
c     &                              / PROPCE(IEL,IPCRO2)
c     &                              *( PROPCE(IEL,IPCTE1)
c     &                                -PROPCE(IEL,IPCTE2))/LV
C
c                IF ( PROPCE(IEL,IPCSEC) .LT. 0.D0 ) THEN
c                  PROPCE(IEL,IPCSEC) = 0.D0
c                  NPOIN65 = NPOIN65 + 1
c                endif
C
c               ENDIF
C
            ELSE
              PROPCE(IEL,IPCSEC) = 0.D0
            ENDIF
C
            TMAX = MAX(PROPCE(IEL,IPCSEC),TMAX)
            TMIN = MIN(PROPCE(IEL,IPCSEC),TMIN)
C
C
C
          ENDDO
C
        IF ( IRANGP .GE. 0 ) THEN
          CALL PARMIN(TMIN)
          CALL PARMAX(TMAX)
          CALL PARMIN(YVMIN)
          CALL PARMAX(YVMAX)
          CALL PARCPT(NPYV)
          CALL PARMAX(YYMAX)
        ENDIF

C
          IF (MODNTL.EQ.0) THEN
            WRITE(NFECRA,*) ' POUR LA CLASSE = ',ICLA,NCEL
            WRITE(NFECRA,*) ' NBRE DE PTS XWAT >0   =',NPOIN1
            WRITE(NFECRA,*) '             T < TEBL  =',NPOIN2
            WRITE(NFECRA,*) '             YV     >0 =',NPOIN3
            WRITE(NFECRA,*) '             YV     <0 =',NPOIN4
            WRITE(NFECRA,*) '             T >= TEBL =',NPOIN5
            WRITE(NFECRA,*) '             ANNUL G   =',
     &           NPOIN63+NPOIN65
            WRITE(NFECRA,*) '  ANNUL G   T < TEBL=',
     &           NPOIN63
            WRITE(NFECRA,*) '  ANNUL G   T >= TEBL=',
     &           NPOIN65
            WRITE(NFECRA,*) ' MIN MAX TS = ',TMIN,TMAX
            WRITE(NFECRA,*) ' MIN MAX YV = ',YVMIN,YVMAX
C
            WRITE(NFECRA,*) ' CLIPPING DE YV EN MAX ',NPYV,YYMAX
          ENDIF
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C FORMATS
C----
C
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
