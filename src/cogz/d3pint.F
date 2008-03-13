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
                        SUBROUTINE D3PINT
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , INDPDF ,
     &   DIRMIN , DIRMAX , FDEB   , FFIN   , HREC   ,
     &   FM     , HM     , P      ,
     &   PROPCE ,
     &   W1      )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC ROUTINE PHYSIQUE PARTICULIERE : FLAMME DE DIFFUSION
CFONC  Integration des variables thermodynamiques en fonction de
CFONC   la fraction de melange
CFONC   Rq : Il serait judicieux de ponderer l'integration de la
CFONC            temperature par les CP
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! INDPDF       ! TI !  -> ! INDICTEUR PASSAGE OU NON PAR LES PDF !
CARGU ! DIRMIN       ! TR !  -> ! PDF : DIRAC EN FMIN                  !
CARGU ! DIRMAX       ! TR !  -> ! PDF : DIRAC EN FMAX                  !
CARGU ! FDEB         ! TR !  -> ! PDF : ABSCISSE DEBUT RECTANGLE       !
CARGU ! FFIN         ! TR !  -> ! PDF : ABSCISSE FIN RECTANGLE         !
CARGU ! HREC         ! TR !  -> ! PDF : HAUTEUR RECTANGLE              !
CARGU ! FM           ! TR !  -> ! FRACTION DE MELANGE MOYENNE          !
CARGU ! HM           ! TR !  -> ! ENTHALPIE MASSIQUE MOYENNE           !
CARGU !              !    !     !  SI ECOULEMENT PERMEATIQUE           !
CARGU ! P            ! TR !  -> ! PRESSION                             !
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     ! CELLULES ( Concentrations, Temp. )   !                        !
CARGU ! W1           ! TR !  -  ! TABLEAU DE TAVAIL                    !
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
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C
C***********************************************************************
C ARGUMENTS
C
      INTEGER          NCELET, NCEL
      INTEGER          INDPDF(NCELET)
      DOUBLE PRECISION DIRMIN(NCELET), DIRMAX(NCELET)
      DOUBLE PRECISION FDEB(NCELET), FFIN(NCELET), HREC(NCELET)
      DOUBLE PRECISION FM(NCELET), HM(NCELET), P(NCELET)
      DOUBLE PRECISION PROPCE(NCELET,*), W1(NCELET)
C
C
C VARIABLES LOCALES
C
      INTEGER          ICEL, ICG, IPHAS
      INTEGER          IH, IF, JH, JF, IPCROM
      INTEGER          IPCSCA, IPCTEM, IPCKAB, IPCT4, IPCT3
      DOUBLE PRECISION AA1, BB1, AA2, BB2, F1, F2, A, B, FMINI, FMAXI
      DOUBLE PRECISION U, V, C, D, TEMSMM, FSIR
C
C
C***********************************************************************
C
      INTEGER       IPASS
      DATA          IPASS /0/
      SAVE          IPASS

C
C***********************************************************************
C
C=======================================================================
C 0. ON COMPTE LES PASSAGES
C=======================================================================
C
      IPASS = IPASS + 1
C
C
C=======================================================================
C 1. INTEGRATION DES NGAZG FRACTIONS MASSIQUES D'ESPECES GLOBALES
C=======================================================================
C
C ---> En flamme de diffusion chimie 3 points :
C      - il n'y a qu'une seule reaction globale (IR= )
C      - le taux de melange f varie entre 0 et 1
      FSIR = FS(1)
      FMINI = ZERO
      FMAXI = 1.D0
C
      DO ICEL = 1, NCEL
C
        DO ICG = 1, NGAZG
C
C --->  Determination des parametres des droites Variables(f)
C       Y = A + B F
C       En flamme de diffusion, il n'y a qu'une seule reaction globale
C         (IR=1)
C       Par definition les fractions massiques des especes globales
C         sont alors
C
C         F       0  FS  1
C         YFUEL   0  0   1
C         YOXYD   1  0   0
C         YPROD   0  1   0
C
          IF ( ICG.EQ.1 ) THEN
C         Fuel
            AA1 =  ZERO
            BB1 =  ZERO
            AA2 = -FSIR/(1.D0-FSIR)
            BB2 =  1.D0/(1.D0-FSIR)

          ELSEIF ( ICG.EQ.2 ) THEN
C         Oxydant
            AA1 =  1.D0
            BB1 = -1.D0/FSIR
            AA2 =  ZERO
            BB2 =  ZERO
          ELSEIF ( ICG.EQ.3 ) THEN
C         Produits
            AA1 =  ZERO
            BB1 =  1.D0/FSIR
            AA2 =  1.D0/(1.D0-FSIR)
            BB2 = -1.D0/(1.D0-FSIR)
          ENDIF
C
          IPCSCA = IPPROC(IYM(ICG))
C
          IF ( INDPDF(ICEL) .EQ. 1 ) THEN
C
C ---> Integration de la PDF
C
            PROPCE(ICEL,IPCSCA) = DIRMIN(ICEL) * ( AA1 + BB1 * FMINI )
     &                          + DIRMAX(ICEL) * ( AA2 + BB2 * FMAXI )
            IF ( FDEB(ICEL).LT.FSIR ) THEN
              F1 = FDEB(ICEL)
              F2 = MIN( FSIR,FFIN(ICEL) )
              PROPCE(ICEL,IPCSCA) = PROPCE(ICEL,IPCSCA)
     &             + HREC(ICEL)*(F2-F1)*(AA1+BB1*5.D-1*(F2+F1))
            ENDIF
            IF ( FFIN(ICEL).GT.FSIR ) THEN
              F1 = MAX(FSIR,FDEB(ICEL))
              F2 = FFIN(ICEL)
              PROPCE(ICEL,IPCSCA) = PROPCE(ICEL,IPCSCA)
     &             + HREC(ICEL)*(F2-F1)*(AA2+BB2*5.D-1*(F2+F1))
            ENDIF
          ELSE
C
C ---> Degenerescence sur la valeur moyenne
C
            IF ( FM(ICEL).LE.FSIR ) THEN
              PROPCE(ICEL,IPCSCA) = AA1+BB1*FM(ICEL)
            ELSE
              PROPCE(ICEL,IPCSCA) = AA2+BB2*FM(ICEL)
            ENDIF
C
          ENDIF
C
        ENDDO
C
      ENDDO
C
C
C=======================================================================
C 2. DETERMINATION LOCALE DE L'ENTHALPIE DES GAZ STOECHIOMETRIQUES
C    BRULES EN PERMEATIQUE (non adiab)
C    (stocke dans le tableau W1)
C=======================================================================
C
C ---> Calcul de HSTOE dans W1
C
C ---- Initialisation
C
      DO ICEL = 1, NCEL
        W1(ICEL) = HSTOEA
      ENDDO
C
      IF ( IPPMOD(ICOD3P).EQ.1 ) THEN
C
        CALL D3PHST
C       ===========
     &  ( NCELET , NCEL    , INDPDF ,
     &    DIRMIN , DIRMAX  , FDEB   , FFIN   , HREC   ,
     &    FM     , HM      ,
     &    W1      )
C
      ENDIF
C
C
C=======================================================================
C 3. INTEGRATION a) DE LA TEMPERATURE
C                b) DU COEFFICIENT D'ABSORPTION si rayonnement
C                c) DES TERME T^4 et T^3 si rayonnement
C                d) DE LA MASSE VOLUMIQUE
C=======================================================================
C
C
C ---> Positions des variables, coefficients
C
      IPCTEM = IPPROC(ITEMP)
      IF ( IRAYPP.GT.0 ) THEN
        IPCKAB = IPPROC(ICKABS)
        IPCT4  = IPPROC(IT4M)
        IPCT3  = IPPROC(IT3M)
      ENDIF
C
      IPHAS = 1
      IPCROM = IPPROC(IROM(IPHAS))
C
      DO ICEL = 1, NCEL
C
        IF ( INDPDF(ICEL) .EQ. 1 ) THEN
C
C ---> Integration de la PDF
C
          IH = 1
          DO JH = 1,(NMAXH-1)
            IF ( W1(ICEL).GT.HH(JH+1) .AND. W1(ICEL).LE.HH(JH) )
     &           IH = JH
          ENDDO
          IF ( W1(ICEL) .GE. HH(1)     ) IH = 1
          IF ( W1(ICEL) .LE. HH(NMAXH) ) IH = NMAXH-1
          PROPCE(ICEL,IPCTEM) = DIRMIN(ICEL)*TINOXY +
     &                          DIRMAX(ICEL)*TINFUE
          TEMSMM = DIRMIN(ICEL)/WMOLG(2)*TINOXY
     &           + DIRMAX(ICEL)/WMOLG(1)*TINFUE
          IF ( IRAYPP.GT.0 ) THEN
            PROPCE(ICEL,IPCKAB) =
     &        DIRMIN(ICEL)*CKABSG(2)  + DIRMAX(ICEL)*CKABSG(1)
            PROPCE(ICEL,IPCT4) =
     &        DIRMIN(ICEL)*TINOXY**4 + DIRMAX(ICEL)*TINFUE**4
            PROPCE(ICEL,IPCT3) =
     &        DIRMIN(ICEL)*TINOXY**3 + DIRMAX(ICEL)*TINFUE**3
          ENDIF
          IF = 1
          DO JF = 1, (NMAXF-1)
            IF ( FDEB(ICEL).GE.FF(JF) .AND.
     &           FDEB(ICEL).LT.FF(JF+1) ) IF = JF
          ENDDO
          IF ( FDEB(ICEL) .LE. FF(1)     ) IF = 1
          IF ( FDEB(ICEL) .GE. FF(NMAXF) ) IF = NMAXF-1
          F2 = ZERO
          F1 = FDEB(ICEL)
          DO WHILE ( (FFIN(ICEL)-F2).GT.EPZERO )
            F2 = MIN(FF(IF+1),FFIN(ICEL))
C          Dans le tableau TFH,
C           on extrait sur chaque ligne i : T = Ai+Bi*F
C           et on construit pour la valeur courante de HSTOE (W1)
C                               T = A+B*F
            AA1 = TFH(IF,IH)
            BB1 = (TFH(IF+1,IH)-TFH(IF,IH))/(FF(IF+1)-FF(IF))
            AA2 = TFH(IF,IH+1)
            BB2 = (TFH(IF+1,IH+1)-TFH(IF,IH+1))/(FF(IF+1)-FF(IF))
            A = AA1 + (W1(ICEL)-HH(IH))*(AA2-AA1)/(HH(IH+1)-HH(IH))
            B = BB1 + (W1(ICEL)-HH(IH))*(BB2-BB1)/(HH(IH+1)-HH(IH))
            A = A - B*FF(IF)
C
C ----- Calcul de la temperature par integration
C
            PROPCE(ICEL,IPCTEM) = PROPCE(ICEL,IPCTEM)
     &         + HREC(ICEL)*(F2-F1)*(A+B*(F1+F2)/2.D0)
C
C ----- Preparation aux calculs du coefficient d'absorption
C                               de T^4 et de T^3
C         Cote pauvre
C           UNSMM = (FS-F)/FS / WMOLG(2)+ F/FS / WMOLG(3)
C           CKABS = (FS-F)/FS * CKABSG(2) + F/FS * CKABSG(3)
C         Cote riche
C           UNSMM = (F-FS)/(1-FS)/WMOLG(1) + (1-F)/(1-FS)/WMOLG(3)
C           CKABS = (F-FS)/(1-FS)*CKABSG(1)  + (1-F)/(1-FS)*CKABSG(3)
C         Partout
C           UNSMM = c + df
C           CKABS = u + vF
C           TEMSMM = T*UNSMM = (c+df)*(a+bf) = ca +(cb+ad)f + bd f^2
C           T^4 = (a+bf)^4
C               = a4 + 4a3b f + 6a2b2 f^2 + 4ab3 f^3 + b4 f^4
C           T^3 = (a+bf)^3
C               = a3 + 3a2b f + 3ab2 f^2 + b3 f^3

C
            IF ( F1.LT.FSIR ) THEN
C         On a demarre cote pauvre
              C =   1.D0/WMOLG(2)
              D = (-1.D0/WMOLG(2)+1.D0/WMOLG(3))/FSIR
            ELSE
C         On termine cote riche (en commencant avec f1=fs)
              C = (  -FSIR/WMOLG(1)+1.D0/WMOLG(3))/(1.D0-FSIR)
              D = (   1.D0/WMOLG(1)-1.D0/WMOLG(3))/(1.D0-FSIR)
            ENDIF
C
            IF ( IRAYPP.GT.0 ) THEN
              IF ( F1.LT.FSIR ) THEN
C         On a demarre cote pauvre
                U =   CKABSG(2)
                V = (-CKABSG(2)+ CKABSG(3))/FSIR
              ELSE
C         On termine cote riche (en commencant avec f1=fs)
                U = (-FSIR*CKABSG(1)+ CKABSG(3))/(1.D0-FSIR)
                V = (      CKABSG(1)- CKABSG(3))/(1.D0-FSIR)
              ENDIF
C
C ----- Calcul du coefficient d'absorption
C           et des termes T^4 et de T^3 (si rayonnement)
C
              PROPCE(ICEL,IPCKAB) = PROPCE(ICEL,IPCKAB) +
     &          HREC(ICEL)*( U*(F2-F1) + V*(F2**2-F1**2)*0.5D0 )
C
              PROPCE(ICEL,IPCT4) = PROPCE(ICEL,IPCT4) +
     &          HREC(ICEL)*
     &      (      A**4            * (F2-F1)
     &   +   (4.D0*A**3  *B      ) * (F2**2-F1**2)/2.D0
     &   +   (6.D0*(A**2)*(B**2) ) * (F2**3-F1**3)/3.D0
     &   +   (4.D0*A     *(B**3) ) * (F2**4-F1**4)/4.D0
     &   +   (            (B**4) ) * (F2**5-F1**5)/5.D0  )
C
              PROPCE(ICEL,IPCT3) = PROPCE(ICEL,IPCT3) +
     &          HREC(ICEL)*
     &      (      (A**3)          * (F2-F1)
     &   +   (3.D0*(A**2)*B      ) * (F2**2-F1**2)/2.D0
     &   +   (3.D0*A     *(B**2) ) * (F2**3-F1**3)/3.D0
     &   +   (            (B**3) ) * (F2**4-F1**4)/4.D0  )
C
            ENDIF
C
C ----- Calcul du terme Temperature/masse molaire
C
            TEMSMM = TEMSMM + HREC(ICEL)*
     &        ( A*C       * (F2-F1)
     &        + (C*B+A*D) * (F2**2-F1**2)/2.D0
     &        +  B*D      * (F2**3-F1**3)/3.D0 )
C
            IF = IF+1
            F1 = F2
          ENDDO
C
        ELSE
C
C ---> Degenerescence sur la valeur moyenne
C
          IH = 1
          DO JH = 1, (NMAXH-1)
            IF ( W1(ICEL).GT.HH(JH+1) .AND. W1(ICEL).LE.HH(JH) )
     &            IH = JH
          ENDDO
          IF ( W1(ICEL) .GE. HH(1)     ) IH =1
          IF ( W1(ICEL) .LE. HH(NMAXH) ) IH =NMAXH-1
          IF = 1
          DO JF = 1, (NMAXF-1)
            IF ( FM(ICEL).GE.FF(JF) .AND. FM(ICEL).LT.FF(JF+1) )
     &           IF = JF
          ENDDO
          IF ( FM(ICEL) .LE. FF(1)     ) IF = 1
          IF ( FM(ICEL) .GE. FF(NMAXF) ) IF = NMAXF-1
          AA1 = TFH(IF,IH)
          BB1 = (TFH(IF+1,IH)-TFH(IF,IH))/(FF(IF+1)-FF(IF))
          AA2 = TFH(IF,IH+1)
          BB2 = (TFH(IF+1,IH+1)-TFH(IF,IH+1))/(FF(IF+1)-FF(IF))
          A  = AA1 + (W1(ICEL)-HH(IH))*(AA2-AA1)/(HH(IH+1)-HH(IH))
          B  = BB1 + (W1(ICEL)-HH(IH))*(BB2-BB1)/(HH(IH+1)-HH(IH))
          A  = A - B*FF(IF)
C
C ----- Calcul de la temperature a partir de la valeur moyenne
C
          PROPCE(ICEL,IPCTEM) = A+B*FM(ICEL)
C
          IF ( FM(ICEL).LT.FSIR ) THEN
C         On a demarre cote pauvre
            C =   1.D0/WMOLG(2)
            D = (-1.D0/WMOLG(2)+1.D0/WMOLG(3))/FSIR
          ELSE
C         On termine cote riche (en commencant avec f1=fs)
            C = (  -FSIR/WMOLG(1)+1.D0/WMOLG(3))/(1.D0-FSIR)
            D = (   1.D0/WMOLG(1)-1.D0/WMOLG(3))/(1.D0-FSIR)
           ENDIF
C
          IF ( IRAYPP.GT.0 ) THEN
            IF ( FM(ICEL).LT.FSIR ) THEN
C         On a demarre cote pauvre
              U =   CKABSG(2)
              V = (-CKABSG(2)+ CKABSG(3))/FSIR
            ELSE
C         On termine cote riche (en commencant avec f1=fs)
              U = (-FSIR*CKABSG(1)+ CKABSG(3))/(1.D0-FSIR)
              V = (      CKABSG(1)- CKABSG(3))/(1.D0-FSIR)
            ENDIF
C
C ----- Calcul du coefficient d'absorption
C         et des termes T^4 et de T^3
C         a partir de la valeur moyenne (si rayonnement)
C
            PROPCE(ICEL,IPCKAB) = U + V*FM(ICEL)
            PROPCE(ICEL,IPCT4) = A**4
     &       + (4.D0*(A**3)*B      ) * FM(ICEL)
     &       + (6.D0*(A**2)*(B**2) ) * FM(ICEL)**2
     &       + (4.D0*A     *(B**3) ) * FM(ICEL)**3
     &       + (            (B**4) ) * FM(ICEL)**4
C
            PROPCE(ICEL,IPCT3) = A**3
     &       + ( 3.D0*(A**2)*B      ) * FM(ICEL)
     &       + ( 3.D0*A     *(B**2) ) * FM(ICEL)**2
     &       + (             (B**3) ) * FM(ICEL)**3
C
          ENDIF
C
C ----- Calcul du terme Temperature/masse molaire
C
          TEMSMM = A*C +(C*B+A*D)*FM(ICEL) + B*D*FM(ICEL)**2
C
        ENDIF
C
C ---> Calcul de la masse volumique
C
        IF (IPASS.GT.1.OR.(ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1)) THEN
          PROPCE(ICEL,IPCROM) = SRROM*PROPCE(ICEL,IPCROM)
     &                        + (1.D0-SRROM)*
     &                        ( P0(IPHAS)/(RR*TEMSMM) )
        ENDIF
C
      ENDDO
C
C
      RETURN
      END
C
c@z
