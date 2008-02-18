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
                        SUBROUTINE PPINII
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  INITIALISATION PAR DEFAUT DES COMMONS PHYSIQUE PARTICULIERE
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
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
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "fuincl.h"
      INCLUDE "ppincl.h"
      INCLUDE "elincl.h"
C
C***********************************************************************
C
C VARIABLES LOCALES
C
      INTEGER         IPP, ICLA, ICHA, IESP, IGG, IT, IR, IH, IF, IZONE
      INTEGER         ISOL, IGE, IAT , III
      INTEGER         IDIRAC
C
C
C***********************************************************************
C
C=======================================================================
C 1. REMPLISSAGE INCLUDE ppincl.h
C                INCLUDE GENERAL PROPRE A LA PHYSIQUE PARTICULIERE
C=======================================================================
C
C ---> Initialisation tableau IPPMOD
C
      DO IPP = 1, NMODMX
        IPPMOD(IPP) = -1
      ENDDO
C
C ---> Initialisation pour le rayonnement
C
      IRAYPP = 0
C
C ---> Initialisation pour la combustion gaz
C       Variables transportees
      IFM    = 0
      IFP2M  = 0
      IYGFM  = 0
      ICM    = 0
      ICP2M  = 0
      IFPCPM = 0
      IYFM   = 0
      IYFP2M = 0
      ICOYFP = 0
C       Variables algebriques ou d'etat
      NSALPP = 0
      NSALTO = 0
      ITEMP  = 0
      DO IGG = 1, NGAZGM
        IYM(IGG) = 0
      ENDDO
      ICKABS = 0
      IT4M   = 0
      IT3M   = 0
      DO IDIRAC = 1, NDRACM
        IRHOL (IDIRAC) = 0
        ITEML (IDIRAC) = 0
        IFMEL (IDIRAC) = 0
        IFMAL (IDIRAC) = 0
        IAMPL (IDIRAC) = 0
        ITSCL (IDIRAC) = 0
        IMAML (IDIRAC) = 0
      ENDDO
C
C ---> Initialisation pour la combustion du charbon
C       Variables transportees
      DO ICHA = 1, NCHARB
        IF1M(ICHA) = 0
        IF2M(ICHA) = 0
      ENDDO
      IF3M   = 0
      IF4P2M = 0
      IF5M   = 0
      IYCO2  = 0
      DO ICLA = 1, NCLCPM
        IXCK(ICLA)    = 0
        IXCH(ICLA)   = 0
        INP(ICLA)    = 0
        IH2(ICLA)    = 0
        IXWT(ICLA)   = 0
      ENDDO
C       Variables algebriques ou d'etat
      ITEMP1 = 0
      DO IGE = 1, NGAZEM
        IYM1(IGE) = 0
      ENDDO
      IMMEL = 0
      DO ICLA = 1, NCLCPM
        IX2(ICLA)    = 0
        ITEMP2(ICLA) = 0
        IROM2(ICLA)  = 0
        IDIAM2(ICLA) = 0
        IGMDCH(ICLA) = 0
        IGMDV1(ICLA) = 0
        IGMDV2(ICLA) = 0
        IGMHET(ICLA) = 0
        IGMSEC(ICLA) = 0
      ENDDO
C
C ---> Initialisation pour la combustion fuel
C       Variables transportees
C
      IHM     = 0
      ING     = 0
      IYFOL   = 0
      IHLF    = 0
      IFVAP   = 0
      IFHTF   = 0
      IF4P2M  = 0
C
C       Variables algebriques ou d'etat
C
      ITEMP1 = 0
      DO IGE = 1, NGAZEM
        IYM1(IGE) = 0
      ENDDO
C
      IX3     = 0
      ITEMP3  = 0
      IROM3   = 0
      IDIAM3  = 0
      IGMEVA  = 0
      IGMHTF  = 0
C
C ---> Initialisation pour l'effet Joule , Arc electrique et conduction ionique
C
      IPOTR     = 0
      IPOTI     = 0
      IPOTVA(1) = 0
      IPOTVA(2) = 0
      IPOTVA(3) = 0
      DO IESP=1,NGAZGM
        IYCOEL(IESP) = 0
      ENDDO
C
C ---> Initialisation de l'enthalpie commune a tous
C
      IHM       = 0
C
C ---> Coefficient de relation de la masse volumique
C      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
      SRROM = 1.D0
C
C
C=======================================================================
C 2. REMPLISSAGE INCLUDE ppthch.h
C                INCLUDE THERMOCHIMIE POUR LA PHYSIQUE PARTICULIERE
C=======================================================================
C
C
C ---> Initialisation Common / TCHPPI /
C
      NPO   = 0
      NGAZE = 0
      NGAZG = 0
      NATO  = 0
      NRGAZ = 0
C
C ---> Initialisation Common / TCHPPR /
C
      DO IT = 1, NPOT
        TH(IT) = ZERO
      ENDDO
C
      DO IR = 1, NRGAZM
        FS(IR) = ZERO
      ENDDO
C
      DO IGG = 1, NGAZGM
        DO IT = 1, NPOT
          EHGAZG(IGG,IT) = ZERO
        ENDDO
        DO IR = 1, NRGAZM
          STOEG(IGG,IR) = ZERO
        ENDDO
        WMOLG(IGG) = ZERO
        CKABSG(IGG)= ZERO
      ENDDO
C
      CKABS1 = ZERO
C
      DIFTL0 = ZERO
C
      DO IGE = 1, NGAZEM
        DO IT = 1, NPOT
          EHGAZE(IGE,IT) = ZERO
        ENDDO
        WMOLE(IGE) = ZERO
      ENDDO
C
      DO IAT = 1, NATOM
        WMOLAT(IAT) = ZERO
      ENDDO
C
      XCO2 = ZERO
      XH2O = ZERO
C
C
C=======================================================================
C 3. REMPLISSAGE INCLUDE coincl.h
C                INCLUDE POUR LA PHYSIQUE PARTICULIERE RELATIF A
C                LA COMBUSTION GAZ
C=======================================================================
C
C ---> Modele de flamme de diffusion (chimie 3 points)
C
      NMAXH = 0
      NMAXF = 0
      TINOXY = ZERO
      TINFUE = ZERO
      DO IZONE = 1, NOZPPM
        IENTOX(IZONE) = 0
        IENTFU(IZONE) = 0
      ENDDO
      HINFUE = ZERO
      HINOXY = ZERO
      HSTOEA = ZERO
      DO IH = 1, NMAXHM
        HH(IH) = ZERO
      ENDDO
      DO IF = 1, NMAXFM
        FF(IF)= ZERO
        DO IH = 1, NMAXHM
          TFH(IF,IH) = ZERO
        ENDDO
      ENDDO
C
C ---> Modele de flamme de premelange (modele EBU et LWC)
C      On prend 300K pour temperature des gaz frais.
C        En suite de calcul c'est ecrase par le fichier suite
C        Cette valeur ne sert que lors du premier appel a ebuphy
C          (ensuite, la valeur imposee dans les CL prend le pas)
C
      CEBU  = 2.5D0
      VREF  = ZERO
      LREF  = ZERO
      TA    = ZERO
      TSTAR = ZERO
      FRMEL = ZERO
      TGF   = 300.D0
      DO IZONE = 1, NOZPPM
        IENTGF(IZONE) = 0
        IENTGB(IZONE) = 0
        QIMP(IZONE)   = ZERO
        FMENT(IZONE)  = ZERO
        TKENT(IZONE)  = ZERO
      ENDDO
      HGF   = ZERO
      TGBAD = ZERO
C
      FMIN = ZERO
      FMAX = 1.D0
      HMIN = ZERO
      HMAX = ZERO
      COEFF1 = ZERO
      COEFF2 = ZERO
      COEFF3 = ZERO
C
C=======================================================================
C 4. REMPLISSAGE INCLUDE cpincl.h
C                INCLUDE POUR LA PHYSIQUE PARTICULIERE RELATIF A
C                LA COMBUSTION CP
C=======================================================================
C
C ---> Donnees relatives au charbon
C
C ----   Par charbon
      NCHARB = 0
      NSOLID = 0
C
      DO ICHA = 1, NCHARM
C
        NCLPCH(ICHA) = 0
C
        CCH(ICHA)    = ZERO
        HCH(ICHA)    = ZERO
        OCH(ICHA)    = ZERO
        ALPHA(ICHA)  = ZERO
        BETA(ICHA)   = ZERO
        PCICH(ICHA)  = ZERO
        RHO0CH(ICHA) = ZERO
C
        CCK(ICHA)    = ZERO
        HCK(ICHA)    = ZERO
        OCK(ICHA)    = ZERO
        RHOCK(ICHA)  = ZERO
        GAMMA(ICHA)  = ZERO
        DELTA(ICHA)  = ZERO
        PCICK(ICHA)  = ZERO
C
        XASHCH(ICHA) = ZERO
        CPASHC(ICHA) = ZERO
        H0ASHC(ICHA) = ZERO
C
        XWATCH(ICHA) = ZERO
C
        IY1CH(ICHA)  = 0
        Y1CH(ICHA)   = ZERO
        A1CH(ICHA)   = ZERO
        E1CH(ICHA)   = ZERO
        IY2CH(ICHA)  = 0
        Y2CH(ICHA)   = ZERO
        A2CH(ICHA)   = ZERO
        E2CH(ICHA)   = ZERO
C
        AHETCH(ICHA) = ZERO
        EHETCH(ICHA) = ZERO
C
        IOCHET(ICHA) = 0
C
        ICH(ICHA)    = 0
        ICK(ICHA)    = 0
        IASH(ICHA)   = 0
        IWAT(ICHA)   = 0
C
      ENDDO
C
      ALPHAM = ZERO
      BETAM  = ZERO
C
      DO ISOL = 1, NSOLIM
        DO IT = 1, NPOT
          EHSOLI(ISOL,IT) = ZERO
        ENDDO
        WMOLS(ISOL) = ZERO
      ENDDO
C
C ----   Par classe
C
      NCLACP = 0
      DO ICLA = 1, NCLCPM
        ICHCOR(ICLA) = 0
        DIAM20(ICLA) = ZERO
        DIA2MN(ICLA) = ZERO
        RHO20(ICLA)  = ZERO
        RHO2MN(ICLA) = ZERO
        XMP0(ICLA)   = ZERO
        XMASH(ICLA)  = ZERO
      ENDDO
C
C ---> Definition des Pointeurs du tableau TBMCR utilise dans cpphy1.F
C      et les sous-programmes appeles
C
        DO ICHA = 1, NCHARM
          IF1MC(ICHA) = 0
          IF2MC(ICHA) = 0
        ENDDO
        IX1MC   = 0
        IX2MC   = 0
        ICHX1F1 = 0
        ICHX2F2 = 0
        ICOF1   = 0
        ICOF2   = 0
C
C ---> Donnees relatives a la combustion des especes gazeuses
C
      DO ICHA = 1, NCHARM
        ICHX1C(ICHA) = 0
        ICHX2C(ICHA) = 0
      ENDDO
      ICHX1 = 0
      ICHX2 = 0
      ICO   = 0
      IO2   = 0
      ICO2  = 0
      IH2O  = 0
      IN2   = 0
C
      XSI   = 3.76D0
      F3MAX =  2.D0*0.012D0 / (2.D0*0.028D0 + XSI*0.028D0)
C
      DO ICHA = 1, NCHARM
        CHX1(ICHA) = ZERO
        CHX2(ICHA) = ZERO
        A1(ICHA)   = ZERO
        A2(ICHA)   = ZERO
        B1(ICHA)   = ZERO
        B2(ICHA)   = ZERO
      ENDDO
C
C
C ---> Donnees complementaires relatives au calcul de rho sur les facettes
C      d'entree
C
      DO IZONE = 1, NOZPPM
        IENTAT(IZONE) = 0
        IENTCP(IZONE) = 0
        TIMPAT(IZONE) = ZERO
        DO ICLA = 1, NCLCPM
          X20(IZONE,ICLA) = ZERO
        ENDDO
      ENDDO
C
C=======================================================================
C 5. REMPLISSAGE INCLUDE fuincl.h
C                INCLUDE POUR LA PHYSIQUE PARTICULIERE RELATIF A
C                LA COMBUSTION FUEL
C=======================================================================
C
C ---> Donnees relatives au fuel
C
      NSOLID = 0
C
      CFOL    = ZERO
      HFOL    = ZERO
      OFOL     = ZERO
      PCIFOL  = ZERO
      RHO0FL = ZERO
C
      CKF    = ZERO
      HKF    = ZERO
      OKF    = ZERO
      RHOKF  = ZERO
      PCIKF = ZERO
C
      IYFOL  = 0
      YFOL   = ZERO
      AFOL   = ZERO
      EFOL   = ZERO
C
      AHETFL = ZERO
      EHETFL = ZERO
C
      IOFHET = 0
C
      IFOL   = 0
      IKF   = 0
C
      DO ISOL = 1, NSOLIM
        DO IT = 1, NPOT
          EHSOLI(ISOL,IT) = ZERO
        ENDDO
        WMOLS(ISOL) = ZERO
      ENDDO
C
C ---> Donnees relatives a la combustion des especes gazeuses
C
      IFOV  = 0
      ICO   = 0
      IO2   = 0
      ICO2  = 0
      IH2O  = 0
      IN2   = 0
C
      XSI   = 3.76D0
      FVAPMX =  2.D0*0.012D0 / (2.D0*0.028D0 + XSI*0.028D0)
C
      FOV  = ZERO
      A  = ZERO
      B  = ZERO
C
C
C ---> Donnees complementaires relatives au calcul de rho sur les facettes
C      d'entree
C
      DO IZONE = 1, NOZPPM
        IENTAT(IZONE)  = 0
        IENTFL(IZONE) = 0
        TIMPAT(IZONE)  = ZERO
        TIMPFL(IZONE) = ZERO
      ENDDO
C
C=======================================================================
C 6. REMPLISSAGE INCLUDE elincl.h
C                INCLUDE POUR LA PHYSIQUE PARTICULIERE RELATIF A
C                EFFET JOULE, ARC ELECTRIQUE, CONDUCTION IONIQUE
C=======================================================================
C
      DO III=1,NTYPMX
        IENTM1(III)=0
        IELPH1(III)=0
        IELPH2(III)=0
        IELPH3(III)=0
        IELNEU(III)=0
      ENDDO
C
      DO IT=1,NPOT
        DO IESP=1,NGAZGM
          RHOEL (IESP,IT) = 0.D0
          CPEL  (IESP,IT) = 0.D0
          SIGEL (IESP,IT) = 0.D0
          XLABEL(IESP,IT) = 0.D0
          VISEL (IESP,IT) = 0.D0
          XKABEL(IESP,IT) = 0.D0
        ENDDO
      ENDDO
C
      DO IESP=1,NGAZGM
        QESPEL(IESP) = 0.D0
        SUSCEP(IESP) = 0.D0
      ENDDO
C
      RETURN
      END
C
c@z
