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
                        SUBROUTINE FUCYNO
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   ,
     &   INDPDF ,
     &   RTP    , PROPCE ,
     &   F1M    , F3M    , F4M   , F1CL  , F3CL , F4CL ,
     &   F4M1   , F4M2   , D4CL  , D4F4  , HREC , F4S3 ,
     &   FUEL1  , FUEL3  , OXYD  , PROD1 , PROD2  ,
     &   XINER  , XH2S   , XSO2   )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CALCUL DE :
CFONC
CFONC           K1 exp(-E1/RT)   (conversion HCN en NO)
CFONC           K2 exp(-E2/RT)   (conversion HCN en N2)
CFONC           K3 exp(-E3/RT)   (NO thermique)
CFONC
CFONC   VALEURS CELLULES
CFONC   ----------------
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT)        !
CARGU ! FVAP         ! TR !  -> ! MOYENNE DU TRACEUR 1 MVl [FOVm+CO]   !
CARGU ! FHET         ! TR !  -> ! MOYENNE DU TRACEUR 3 C hÈterogËne    !
CARGU ! F4P2M        ! TR !  -> ! VARIANCE DU TRACEUR 4 (air)          !
CARGU ! INDPDF       ! TE !  -> ! PASSAGE PAR LES PDF                  !
CARGU ! FUEL1        ! TR !  <- ! FRACTION MASSIQUE FOVm               !
CARGU ! FUEL3        ! TR !  <- ! FRACTION MASSIQUE CO                 !
CARGU ! OXYD         ! TR !  <- ! FRACTION MASSIQUE O2                 !
CARGU ! PROD1        ! TR !  <- ! FRACTION MASSIQUE CO2                !
CARGU ! PROD2        ! TR !  <- ! FRACTION MASSIQUE H2O                !
CARGU ! XINER        ! TR !  <- ! FRACTION MASSIQUE N2                 !
CARGU ! XH2S         ! TR !  <- ! FRACTION MASSIQUE H2S                !
CARGU ! XSO2         ! TR !  <- ! FRACTION MASSIQUE SO2                !
CARGU ! F4M1         ! TR ! <-> ! BORNE MINIMUM                        !
CARGU ! F4M2         ! TR ! <-> ! BORNE MAX                            !
CARGU ! D4CL         ! TR ! <-> ! AMPLITUDE DU PIC DE DIRAC EN F4CL    !
CARGU ! D4F4         ! TR ! <-> ! AMPLITUDE DU PIC DE DIRAC EN 1       !
CARGU !              !    !     ! (AIR PUR)                            !
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
      INCLUDE "parall.h"
      INCLUDE "pointe.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "fuincl.h"
      INCLUDE "ppincl.h"
      INCLUDE "ppcpfu.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET , NCEL
      INTEGER          INDPDF(NCELET)
C
      DOUBLE PRECISION RTP(NCELET,*), PROPCE(NCELET,*)
      DOUBLE PRECISION F1M(NCELET) , F3M(NCELET)  , F4M(NCELET)
      DOUBLE PRECISION F1CL(NCELET), F3CL(NCELET) , F4CL(NCELET)
C
      DOUBLE PRECISION F4M1(NCELET) , F4M2(NCELET) , D4CL(NCELET)
      DOUBLE PRECISION D4F4(NCELET) , HREC(NCELET) , F4S3(NCEL)
C
      DOUBLE PRECISION FUEL1(NCELET), FUEL3(NCELET)
      DOUBLE PRECISION OXYD(NCELET) , PROD1(NCELET), PROD2(NCELET)
      DOUBLE PRECISION XINER(NCELET)
      DOUBLE PRECISION XH2S(NCELET) , XSO2(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IEL , ICLA , I
C
      INTEGER          N1 , N2 , N3 , N4 , N5 , N6 , N7 , N8 , N9
      INTEGER          N10, N11
      INTEGER          NBRINT
      INTEGER          IPDF1 , IPDF2 , IPDF3
      INTEGER          IEXP1 , IEXP2 , IEXP3
      PARAMETER        (NBRINT = 11)
      INTEGER          INTTMP(NBRINT)
      INTEGER          NPART
      PARAMETER        (NPART = 200 )
C
      DOUBLE PRECISION WMO2
C
      DOUBLE PRECISION EE1,EE2,EE3,KK1,KK2,KK3
      DOUBLE PRECISION TS3 , TS3DEN , TS3NUM , U , T , TS3S , DIRAC
      DOUBLE PRECISION Q , R , S , LPA , LRI , TMPGAZ , TFUEL , XMX2
      DOUBLE PRECISION BB1 , BB2 , BB3 , BB4 , WMCO , XO2 , BB
      DOUBLE PRECISION GH1, GH2 ,GH3 ,GH4 ,GH5,GH6, GH7 , TG
      DOUBLE PRECISION GH8, GH9 ,GH10,GH11
      DOUBLE PRECISION DGS , RSF4 , RSCL , BMOY , P , SPDF
      DOUBLE PRECISION VAL(NPART+1),TT(NPART+1) , GS(NPART+1)
C
C***********************************************************************
C
C
C=======================================================================
C 0. INITIALISATIONS
C=======================================================================
C
C
C=======================================================================
C 1. CALCULS PRELIMINAIRES
C=======================================================================
C
C --- Masses molaires
C
      WMO2   = WMOLE(IO2)
C
C --- pointeurs
C
      IEXP1 = IPPROC(IGHCN1)
      IEXP2 = IPPROC(IGHCN2)
      IEXP3 = IPPROC(IGNOTH)
C
C Parametres des lois d'Arrhenius
C
      KK1 = 3.0E12
      EE1 = 3.0E4
      KK2 = 1.2E10
      EE2 = 3.35E4
      KK3 = 3.4E12
      EE3 = 6.69E4
C
C Pour les termes, indicateur de calcul par PDF ou non
C       = 1 --> passage par pdf
C       = 0 --> on ne passe pas par les pdf
C
      IPDF1 = 0
      IPDF2 = 0
      IPDF3 = 0
C
C ---- Initialisation
C
      DO IEL = 1, NCEL
c        PROPCE(IEL,IEXP1)  = ZERO
c        PROPCE(IEL,IEXP2)  = ZERO
c        PROPCE(IEL,IEXP3)  = ZERO
      ENDDO
C
C=======================================================================
C 3. CALCUL SANS LES PDF
C=======================================================================
C
      DO IEL = 1, NCEL
C
        TG  = PROPCE(IEL,IPPROC(ITEMP1))
        PROPCE(IEL,IEXP1)  = KK1*EXP(-EE1/TG)
C
      ENDDO
C
      DO IEL = 1, NCEL
C
        TG  = PROPCE(IEL,IPPROC(ITEMP1))
        XO2 = OXYD(IEL)*PROPCE(IEL,IPPROC(IMMEL))
     &       /WMO2
        IF ( XO2 .GT. 0.018D0 ) THEN
          BB = 0.D0
        ELSE IF ( XO2 .LT. 0.0025D0 ) THEN
          BB= 1.D0
        ELSE
          BB=(0.018D0-XO2)/(0.018D0-0.0025D0)
        ENDIF
C
        PROPCE(IEL,IEXP2)  = KK2*EXP(-EE2/TG)*(XO2**BB)
C
      ENDDO
C
      DO IEL = 1, NCEL
C
        TG  = PROPCE(IEL,IPPROC(ITEMP1))
        XO2 = OXYD(IEL)*PROPCE(IEL,IPPROC(IMMEL))
     &       /WMO2
C
        PROPCE(IEL,IEXP3)  = KK3*EXP(-EE3/TG)*(XO2**0.5D0)
C
      ENDDO
C
C=======================================================================
C 3. CALCUL AVEC LES PDF
C=======================================================================
C
      IF ( IPDF1 .EQ. 1 .OR. IPDF2 .EQ. 1 .OR. IPDF3 .EQ. 1 ) THEN
C
        DO IEL=1,NCEL
C
          IF ( INDPDF(IEL).EQ.1 ) THEN
C
C         Calcul de Tfioul moyen
C
            XMX2  = 0.D0
            TFUEL = 0.D0
C
            DO ICLA=1,NCLAFU
              XMX2 = XMX2 + RTP(IEL,ISCA(IYFOL(ICLA)))
            ENDDO
C
            IF ( XMX2 .GT. 0.D0 ) THEN
              DO ICLA=1,NCLAFU
                TFUEL = TFUEL + RTP(IEL,ISCA(IYFOL(ICLA)))
     &                         *PROPCE(IEL,IPPROC(ITEMP3(ICLA)))
              ENDDO
              TFUEL = TFUEL/XMX2
            ELSE
              TFUEL = PROPCE(IEL,IPPROC(ITEMP1))
            ENDIF
C
C  On recupere la valeur de TAIR
C
            TAIRE = RTP(IEL,ISCA(ITAIRE))
C
c On initialise par les temperatures Tair (en f4) et Tcl (TFUEL) aux extrÅÈmitÅÈs
C
            DIRAC  = D4CL(IEL)*TFUEL + D4F4(IEL)*TAIRE
C
c1 On recupere la valeur de la temperature moyenne
C
            TMPGAZ = PROPCE(IEL,IPPROC(ITEMP1))
C
c1 Integration premier intervalle, entre CL et S3 (stoechio derniÅËre rÅÈaction)
C
c1 Parametres de la droite entre CL et S3
C
c1 Detection des bornes d'integration (celles de la pdf) du domaine pauvre
C
            IF(TFUEL.GT.EPSICP .AND. TMPGAZ.GT.EPSICP) THEN
c
              IF(F4S3(IEL).LE.F4M1(IEL)) THEN
                P=((F4M1(IEL)+F4M2(IEL))/2.-F4S3(IEL))/(1.-F4S3(IEL))
                TS3S=((TMPGAZ-DIRAC)/HREC(IEL)
     &                              /(F4M2(IEL)-F4M1(IEL))-TAIRE*P)
     &                              /(1.D0-P)
              ELSE IF(F4S3(IEL).GE.F4M2(IEL)) THEN
                TS3S=((TMPGAZ-DIRAC)/HREC(IEL)
     &               /(F4M2(IEL)-F4M1(IEL))-TFUEL)
     &               *(F4S3(IEL)-F4CL(IEL))
     &               /((F4M1(IEL)+F4M2(IEL))/2.D0-F4CL(IEL))+TFUEL
              ELSE
                TS3S=(2.*(TMPGAZ-DIRAC)/HREC(IEL)
     &               - TFUEL*(F4S3(IEL)-F4M1(IEL))**2
     &                /(F4S3(IEL)-F4CL(IEL))
     &               - TAIRE*(F4M2(IEL)-F4S3(IEL))**2
     &                /(1.D0-F4S3(IEL)))
     &      /((2.*F4M2(IEL)-F4M1(IEL)-F4S3(IEL))
     &               +(F4M1(IEL)-F4CL(IEL))*(F4S3(IEL)-F4M1(IEL))
     &                                     /(F4S3(IEL)-F4CL(IEL))
     &               -(F4M2(IEL)-F4S3(IEL))**2/(1.D0-F4S3(IEL)))
              ENDIF
C
              BB1 = MAX( F4M1(IEL) , F4CL(IEL) )
              BB2 = MIN( F4M2(IEL) , F4S3(IEL) )
              BB3 = MAX( F4M1(IEL) , F4S3(IEL) )
              BB4 = MIN( F4M2(IEL) , 1.D0 )
C
C On definit quelques parametres intermediaires pour simplifier le code
C
              IF (BB2 .GT. BB1 ) THEN
                LRI = HREC(IEL)
              ELSE
                LRI = 0.D0
              ENDIF
              IF (BB4 .GT. BB3 ) THEN
                LPA = HREC(IEL)
              ELSE
                LPA = 0.D0
              ENDIF
C
              P = F4S3(IEL) - F4CL(IEL)
              Q = BB2**2 - BB1**2
              R = BB2 - BB1
              S = 1.D0 - F4S3(IEL)
              T = BB4**2 - BB3**2
              U = BB4 - BB3
C
c On calcul de TS3
C
c           TS3=(TMPGAZ-DIRAC+TFUEL*(((LRI*Q)/(2.D0*P))-((LRI*F4S3*R)/P))
c     &         - TAIRE*( ((LPA*T)/(2.D0*S)) +(LPA*U) -((LPA*U)/S)))
c     &         /
c     &         (((LRI*Q)/(2.D0*P))+(LRI*R)-((LRI*F4S3*R)/P)-((LPA*T)
c     &         /(2.D0*S)) +((LPA*U)/S) )
C
              GH1 = -LRI*Q/(P*2.D0)
              GH2 =  LRI*R*F4S3(IEL)/P
              GH3 =  LPA*T/(2.D0*S)
              GH4 =  LPA*U
              GH5 =  LPA*U/S
              GH6 =  LRI*Q/(P*2.D0)
              GH7 =  LRI*R
              GH8 =  F4S3(IEL)*LRI*R/P
              GH9 = LPA*T/(S*2.D0)
              GH10= LPA*U/S
              GH11 = TMPGAZ - DIRAC
C
              TS3NUM = GH11 - TFUEL*(GH1+GH2) - TAIRE*(GH3+GH4-GH5)
              TS3DEN = GH6 + GH7 -GH8 - GH9 + GH10
C

              TS3 = TS3NUM / TS3DEN
              IF(ABS(TS3-TS3S).GT.1.D0) THEN
                WRITE(NFECRA,*) 'TS3 paul-TS3 sandro ', IEL,TS3,TS3S
              ENDIF
            ENDIF

            SPDF = HREC(IEL) * (F4M2(IEL) - F4M1(IEL))
C
C Integration
C
            XO2 = OXYD(IEL)*PROPCE(IEL,IPPROC(IMMEL))
     &           /WMO2
C
            DO I = 1, NPART+1
              GS(I) = F4M1(IEL) + DBLE(I-1)/DBLE(NPART)
     &                           *(F4M2(IEL)-F4M1(IEL))
C
              IF( GS(I) .LE. F4S3(IEL) ) THEN
                TT(I) = (TS3-TFUEL)/(F4S3(IEL)-F4CL(IEL))* GS(I)
     &               + TS3 - F4S3(IEL)*(TS3 - TFUEL)
     &                                / ( F4S3(IEL) - F4CL(IEL) )
C
              ELSE
                TT(I) = (TAIRE -TS3)/(1.D0-F4S3(IEL))*GS(I)
     &               + TAIRE - (TAIRE - TS3)/(1.D0-F4S3(IEL))
C
              ENDIF
            ENDDO
C
            IF(XO2.GT.0.018D0) THEN
              BMOY=0.D0
            ELSE IF(XO2 .LT. 0.0025D0) THEN
              BMOY=1.D0
            ELSE
              BMOY=(0.018D0-XO2)/(0.018D0-0.0025D0)
            ENDIF
C
C      DGS est le pas d'integration
C
            DGS = ( F4M2(IEL) - F4M1(IEL) ) / DBLE(NPART)
C
C Calcul de K1*EXP(-E1/T)
C
            IF ( IPDF1 .EQ. 1 ) THEN
C
              RSF4 = KK1*EXP(-EE1/TAIRE)*D4F4(IEL)
              RSCL = KK1*EXP(-EE1/TFUEL)*D4CL(IEL)
C
              DO I = 1, NPART+1
                VAL(I) = KK1 * EXP(-EE1/TT(I)) * HREC(IEL)
              ENDDO
C
              PROPCE(IEL,IEXP1)= RSF4 + RSCL
C
              DO I = 1, NPART
                PROPCE(IEL,IEXP1) = PROPCE(IEL,IEXP1)
     &                           +0.5D0*DGS*(VAL(I)+VAL(I+1))
              ENDDO
C
            ENDIF
C
C  Calcul de K2*EXP(-E2/T)
C
            IF ( IPDF2 .EQ. 1 ) THEN
C
              IF ( XO2 .GT. 0.D0 ) THEN
                PROPCE(IEL,IEXP2) = KK2*EXP(-EE2/TAIRE)
     &                                 *D4F4(IEL)*XO2**BMOY
     &                             +KK2*EXP(-EE2/TFUEL)
     &                                 *D4CL(IEL)*XO2**BMOY
C
                DO I = 1, NPART+1
                  VAL(I) = KK2*EXP(-EE2/TT(I))*HREC(IEL)
                ENDDO
C
                DO I = 1, NPART
                  PROPCE(IEL,IEXP2) = PROPCE(IEL,IEXP2)
     &                               +0.5D0*DGS*(VAL(I)+VAL(I+1))
     &                                     *XO2**BMOY
                ENDDO
              ELSE
                PROPCE(IEL,IEXP2) = 0.D0
              ENDIF
C
            ENDIF
C
C  Calcul de K3*EXP(-E3/T)
C
            IF ( IPDF3 .EQ. 1 ) THEN
C
              IF ( XO2 .GT. 0.D0 ) THEN
                PROPCE(IEL,IEXP3) = KK3*EXP(-EE3/TAIRE)
     &                                 *D4F4(IEL)*XO2**(0.5D0)
     &                             +KK3*EXP(-EE3/TFUEL)
     &                                 *D4CL(IEL)*XO2**(0.5D0)
C
                DO I = 1, NPART+1
                  VAL(I) = KK3*EXP(-EE3/TT(I))*HREC(IEL)
                ENDDO
C
                DO I = 1, NPART
                  PROPCE(IEL,IEXP3)= PROPCE(IEL,IEXP3)
     &                     + 0.5D0*DGS*(VAL(I)+VAL(I+1))*XO2**(0.5D0)
                ENDDO
              ELSE
                PROPCE(IEL,IEXP3)= 0.D0
              ENDIF
            ENDIF
C
          ENDIF
C
        ENDDO
C
      ENDIF
C
      RETURN
      END
c@z
