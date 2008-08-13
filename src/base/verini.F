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
                        SUBROUTINE VERINI
C                       *****************
C     -------------------------------------------------------------
     & ( IOK    )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  VERIFICATION DES PARAMETRES DE CALCUL
CFONC    APRES INTERVENTION UTILISATEUR
CFONC    (COMMONS)
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
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "cstnum.h"
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "albase.h"
      INCLUDE "alstru.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppincl.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "radiat.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IOK
C
C VARIABLES LOCALES
C
      CHARACTER        CHAINE*80, CHAIN2*80
      INTEGER          II    , IIS   , JJ    , IISCT , IPHAS
      INTEGER          ISCAL , IEST  , IIESCA, IVAR
      INTEGER          IUIPH , IVIPH , IWIPH, IKIPH, IEIPH
      INTEGER          ITRBPH, NBSCCP
      INTEGER          IPP   , IPPM  , IMGROK, NBCCOU
      INTEGER          IOKPRE, INDEST, ITESTS, IIIDEF, ISTOP
      INTEGER          IRESOP, IPOLOP
      DOUBLE PRECISION TESTTH
      DOUBLE PRECISION ARAKFR
C
C***********************************************************************
C
C=======================================================================
C 1. ENTREES SORTIES entsor.h : formats 1000
C=======================================================================
C
C --- Dimension
C
      IF (NDIM.NE.3) THEN
        WRITE(NFECRA,1100)NDIM
        IOK = IOK + 1
      ENDIF
C
C --- Suite, Chrono, Historiques, Listing
C
      IF(NTSUIT.LT.-1) THEN
        WRITE(NFECRA,1200) 'NTSUIT (Periode   Sortie Suite  )',NTSUIT
        IOK = IOK + 1
      ENDIF
C
      IF(NTCHR.NE.-1.AND.NTCHR.LE.0) THEN
        WRITE(NFECRA,1210) 'NTCHR  (Periode   Sortie Chrono.)',NTCHR
        IOK = IOK + 1
      ENDIF
C
      DO IPP = 2, NVPPMX
        IF(ICHRVR(IPP).NE.1.AND.ICHRVR(IPP).NE.0) THEN
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,1220)CHAINE(1:8),IPP,ICHRVR(IPP)
          IOK = IOK + 1
        ENDIF
      ENDDO
C
      IF(NCAPT.LT.0.OR.NCAPT.GT.NCAPTM) THEN
        WRITE(NFECRA,1230)NCAPTM,NCAPT
        IOK = IOK + 1
      ENDIF
C
      IF(NTHIST.LE.0.AND.NTHIST.NE.-1) THEN
        WRITE(NFECRA,1210) 'NTHIST (Periode   Sortie Histo. )',NTHIST
        IOK = IOK + 1
      ENDIF
C
      IF(NTHSAV.LT.-1) THEN
        WRITE(NFECRA,1200) 'NTHSAV (Periode   Save   Histo. )',NTHSAV
        IOK = IOK + 1
      ENDIF
C
      DO IPP = 2, NVPPMX
        IF( IHISVR(IPP,1).GT.NCAPT.OR.
     &     (IHISVR(IPP,1).LT.0.AND.IHISVR(IPP,1).NE.-1) ) THEN
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,1240)CHAINE(1:8),IPP,NCAPT,IHISVR(IPP,1)
          IOK = IOK + 1
        ENDIF
      ENDDO
C
      DO IPP = 2, NVPPMX
        IF( (IHISVR(IPP,1).GT.0.AND.IHISVR(IPP,1).LT.NCAPT)) THEN
          DO JJ = 1, IHISVR(IPP,1)
            IF(IHISVR(IPP,JJ+1).LE.0.OR.IHISVR(IPP,JJ+1).GT.NCAPT) THEN
              CHAINE=NOMVAR(IPP)
              WRITE(NFECRA,1250)
     &          CHAINE(1:8),IPP,JJ+1,NCAPT,IHISVR(IPP,JJ+1)
              IOK = IOK + 1
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
      DO IPP = 2, NVPPMX
        IF( ILISVR(IPP).NE.0.AND.ILISVR(IPP).NE.1) THEN
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,1260)CHAINE(1:8),IPP,ILISVR(IPP)
          IOK = IOK + 1
        ENDIF
      ENDDO
C
      DO IPP = 2, NVPPMX
        IF(ITRSVR(IPP).NE.0) THEN
          IF(ITRSVR(IPP).LE.0.OR.ITRSVR(IPP).GT.NVAR) THEN
            CHAINE=NOMVAR(IPP)
            WRITE(NFECRA,1270)CHAINE(1:8),IPP,NVAR,ITRSVR(IPP)
            IOK = IOK + 1
          ENDIF
        ENDIF
      ENDDO
C
      IF(NTLIST.NE.-1.AND.NTLIST.LE.0) THEN
        WRITE(NFECRA,1210) 'NTLIST (Periode   Sortie Listing)',NTLIST
        IOK = IOK + 1
      ENDIF
C
C --- Post traitement automatique (bord)
C
      IF(IPSTDV.NE.1.AND.
     &   MOD(IPSTDV,IPSTYP).NE.0.AND.
     &   MOD(IPSTDV,IPSTCL).NE.0.AND.
     &   MOD(IPSTDV,IPSTFT).NE.0) THEN
        WRITE(NFECRA,1300) 'IPSTDV',
     &                     'IPSTYP',IPSTYP,
     &                     'IPSTCL',IPSTCL,
     &                     'IPSTFT',IPSTFT,
     &                     'IPSTFO',IPSTFO,
     &                     'IPSTDV',IPSTDV
        IOK = IOK + 1
      ENDIF
C
C
C=======================================================================
C 2. OPTIONS DU CALCUL : TABLEAUX DE optcal.h : formats 2000
C=======================================================================
C
C --- Dimensions
C
      IF(NPHAS.LT.0.OR.NPHAS.GT.NPHSMX) THEN
        WRITE(NFECRA,2000)'NPHAS ',NPHSMX,NPHAS
        IOK = IOK + 1
      ENDIF
      IF(NSCAL.LT.0.OR.NSCAL.GT.NSCAMX) THEN
        WRITE(NFECRA,2000)'NSCAL ',NSCAMX,NSCAL
        IOK = IOK + 1
      ENDIF
      IF(NSCAUS.LT.0.OR.NSCAUS.GT.NSCAMX) THEN
        WRITE(NFECRA,2000)'NSCAUS',NSCAMX,NSCAUS
        IOK = IOK + 1
      ENDIF
      IF(NSCAPP.LT.0.OR.NSCAPP.GT.NSCAMX) THEN
        WRITE(NFECRA,2000)'NSCAPP',NSCAMX,NSCAPP
        IOK = IOK + 1
      ENDIF
      IF(NVAR.LT.0.OR.NVAR.GT.NVARMX) THEN
        WRITE(NFECRA,2000)'NVAR  ',NVARMX,NVAR
        IOK = IOK + 1
      ENDIF
C
C --- Rho et visc constants ou variables
C
      DO IPHAS = 1, NPHAS
        IF(IROVAR(IPHAS).NE.0.AND.IROVAR(IPHAS).NE.1) THEN
          WRITE(NFECRA,2201)IPHAS, 'IROVAR',IROVAR(IPHAS)
          IOK = IOK + 1
        ENDIF
        IF(IVIVAR(IPHAS).NE.0.AND.IVIVAR(IPHAS).NE.1) THEN
          WRITE(NFECRA,2201)IPHAS, 'IVIVAR',IVIVAR(IPHAS)
          IOK = IOK + 1
        ENDIF
      ENDDO
C
C --- Definition des equations, schema en temps, schema convectif
C
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          IF( (ICONV (II).NE.0.AND.ICONV (II).NE.1).OR.
     &        (ISTAT (II).NE.0.AND.ISTAT (II).NE.1).OR.
     &        (IDIFF (II).NE.0.AND.IDIFF (II).NE.1).OR.
     &        (IDIFFT(II).NE.0.AND.IDIFFT(II).NE.1).OR.
     &        (THETAV(II).GT.1.D0.OR.THETAV(II).LT.0.D0).OR.
     &        (BLENCV(II).GT.1.D0.OR.BLENCV(II).LT.0.D0).OR.
     &        (ISCHCV(II).NE.0.AND.ISCHCV(II).NE.1).OR.
     &        (ISSTPC(II).NE.0.AND.ISSTPC(II).NE.1)   ) THEN
            CHAINE=NOMVAR(IPP)
            WRITE(NFECRA,2100) CHAINE(1:8),
     &                         ISTAT(II),ICONV(II),
     &                         IDIFF(II),IDIFFT(II),
     &                         THETAV(II),BLENCV(II),ISCHCV(II),
     &                         ISSTPC(II)
            IOK = IOK + 1
          ENDIF
        ENDIF
      ENDDO
C
C
C     Extrap de rho : necessairement rho variable
      DO IPHAS = 1, NPHAS
        IF(IROVAR(IPHAS).EQ.0.AND.IROEXT(IPHAS).GT.0) THEN
          WRITE(NFECRA,2005)IPHAS,IROEXT(IPHAS),IROVAR(IPHAS)
          IOK = IOK + 1
        ENDIF
      ENDDO
C
C     Coherence des vitesses
C       Pour le moment, theta est fixe automatiquement dans modini
C       On conserve quand meme le test pour plus tard puisqu'il est ecrit.
      DO IPHAS = 1, NPHAS
        IF(ABS(THETAV(IV(IPHAS))-THETAV(IU(IPHAS))).GT.EPZERO.OR.
     &     ABS(THETAV(IW(IPHAS))-THETAV(IU(IPHAS))).GT.EPZERO) THEN
          WRITE(NFECRA,2111) IPHAS,THETAV(IU(IPHAS)),THETAV(IV(IPHAS)),
     &                             THETAV(IW(IPHAS))
          IOK = IOK + 1
        ENDIF
      ENDDO
C
C
C     Theta pression : vaut 1
C       Pour le moment, theta est fixe automatiquement dans modini
C         (on ne devrait donc jamais voir cet affichage)
C       On conserve quand meme le test pour plus tard puisqu'il est ecrit.
      DO IPHAS = 1, NPHAS
        JJ = IPR(IPHAS)
        IF(ABS(THETAV(JJ)-1.0D0).GT.EPZERO) THEN
          IPP    = IPPRTP(JJ)
          CHAINE=NOMVAR(IPP)
          WRITE(NFECRA,2112) IPHAS,THETAV(JJ)
          IOK = IOK + 1
        ENDIF
      ENDDO
C
C     Flux de masse et proprietes
C       Pour le moment, theta est fixe automatiquement dans modini
C       Donc pas de test sur sa valeur ici
C
C     On verifie que toutes les phases sont traitees pareil
      TESTTH = THETFL(1)
      ITESTS = ISTMPF(1)
      DO IPHAS = 1, NPHAS
C
        IF(ABS(TESTTH-THETFL(IPHAS)).GT.EPZERO) THEN
          WRITE(NFECRA,2113) IPHAS,TESTTH,THETFL(IPHAS)
          IOK = IOK + 1
        ENDIF
        IF(ITESTS.NE.ISTMPF(IPHAS)) THEN
          WRITE(NFECRA,2114) IPHAS,ITESTS,ISTMPF(IPHAS)
          IOK = IOK + 1
        ENDIF
      ENDDO
C
C     En LES il y a des verification de coherence supplementaires.
C       (simple avertissement si on s'ecarte des choix std)
C        mais stop si on fait plus de 5% d'upwind
C     Schema centre sans/avec test de pente, nwsrsm
      DO IPHAS = 1, NPHAS
        IF(ITYTUR(IPHAS).EQ.4) THEN
          DO II = 1,3
            IF(II.EQ.1) JJ = IU(IPHAS)
            IF(II.EQ.2) JJ = IV(IPHAS)
            IF(II.EQ.3) JJ = IW(IPHAS)
            IPP    = IPPRTP(JJ)
            CHAINE=NOMVAR(IPP)
            IF(ABS(THETAV(JJ)-0.5D0).GT.EPZERO) THEN
              WRITE(NFECRA,2121) IPHAS,CHAINE(1:8),THETAV(JJ)
            ENDIF
            IF (BLENCV(JJ).LT.0.95D0) THEN
              WRITE(NFECRA,2127) IPHAS,CHAINE(1:8),BLENCV(JJ)
              IOK = IOK + 1
            ELSEIF(ABS(BLENCV(JJ)-1.D0).GT.EPZERO) THEN
              WRITE(NFECRA,2122) IPHAS,CHAINE(1:8),BLENCV(JJ)
            ENDIF
            IF(ISSTPC(JJ).EQ.0) THEN
              WRITE(NFECRA,2123) IPHAS,CHAINE(1:8),ISSTPC(JJ)
            ENDIF
          ENDDO
        ENDIF
        IF(ITYTUR(IPHAS).EQ.4.OR.ISCHTP(IPHAS).EQ.2) THEN
          DO II = 1,3
            IF(II.EQ.1) JJ = IU(IPHAS)
            IF(II.EQ.2) JJ = IV(IPHAS)
            IF(II.EQ.3) JJ = IW(IPHAS)
            IPP    = IPPRTP(JJ)
            CHAINE=NOMVAR(IPP)
            IIIDEF = 10
            IF(NSWRSM(JJ).NE.IIIDEF) THEN
              WRITE(NFECRA,2125) IPHAS,CHAINE(1:8),IIIDEF,NSWRSM(JJ)
            ENDIF
          ENDDO
          JJ = IPR(IPHAS)
          IPP    = IPPRTP(JJ)
          CHAINE=NOMVAR(IPP)
          IIIDEF = 5
          IF(NSWRSM(JJ).NE.IIIDEF) THEN
            WRITE(NFECRA,2125) IPHAS,CHAINE(1:8),IIIDEF,NSWRSM(JJ)
          ENDIF
        ENDIF
      ENDDO
      DO II = 1, NSCAL
        IPHAS = IPHSCA(II)
        IF(ITYTUR(IPHAS).EQ.4) THEN
          JJ    = ISCA(II)
          IPP   = IPPRTP(JJ)
          CHAINE=NOMVAR(IPP)
          IF(ABS(THETAV(JJ)-0.5D0).GT.EPZERO) THEN
            WRITE(NFECRA,2121) IPHAS,CHAINE(1:8),THETAV(JJ)
          ENDIF
          IF (BLENCV(JJ).LT.0.95D0) THEN
            WRITE(NFECRA,2127) IPHAS,CHAINE(1:8),BLENCV(JJ)
            IOK = IOK + 1
          ELSEIF(ABS(BLENCV(JJ)-1.D0).GT.EPZERO) THEN
            WRITE(NFECRA,2122) IPHAS,CHAINE(1:8),BLENCV(JJ)
          ENDIF
          IF(ISSTPC(JJ).EQ.1) THEN
            WRITE(NFECRA,2124) IPHAS,CHAINE(1:8),ISSTPC(JJ)
          ENDIF
        ENDIF
        IF(ITYTUR(IPHAS).EQ.4.OR.ISCHTP(IPHAS).EQ.2) THEN
          IIIDEF = 10
          IF(NSWRSM(JJ).NE.IIIDEF) THEN
            WRITE(NFECRA,2125) IPHAS,CHAINE(1:8),IIIDEF,NSWRSM(JJ)
          ENDIF
        ENDIF
      ENDDO
C
C     Test du theta de la viscosite secondaire, du flux de masse et
C     de la viscosite par rapport a celui de la vitesse
      DO IPHAS = 1, NPHAS
        JJ = IU(IPHAS)
        IF( ABS(THETAV(JJ)-1.D0).LT.EPZERO.AND.
     &       (ISTMPF(IPHAS).EQ.2.OR.
     &        ISNO2T(IPHAS).NE.0.OR.
     &        ISTO2T(IPHAS).NE.0.OR.
     &        IROEXT(IPHAS).NE.0.OR.
     &        IVIEXT(IPHAS).NE.0.OR.
     &        ICPEXT(IPHAS).NE.0   ) ) THEN
          WRITE(NFECRA,2131) IPHAS,THETAV(JJ),
     &         ISTMPF(IPHAS),ISNO2T(IPHAS),ISTO2T(IPHAS),
     &         IROEXT(IPHAS),IVIEXT(IPHAS),ICPEXT(IPHAS)
        ENDIF
        IF( ABS(THETAV(JJ)-0.5D0).LT.EPZERO.AND.
     &       (ISTMPF(IPHAS).NE.2.OR.
     &        ISNO2T(IPHAS).NE.1.OR.
     &        ISTO2T(IPHAS).NE.1.OR.
     &        IROEXT(IPHAS).NE.1.OR.
     &        IVIEXT(IPHAS).NE.1.OR.
     &        ICPEXT(IPHAS).NE.1   ) ) THEN
          WRITE(NFECRA,2132) IPHAS,THETAV(JJ),
     &         ISTMPF(IPHAS),ISNO2T(IPHAS),ISTO2T(IPHAS),
     &         IROEXT(IPHAS),IVIEXT(IPHAS),ICPEXT(IPHAS)
        ENDIF
      ENDDO
      DO ISCAL = 1, NSCAL
        IPHAS = IPHSCA(ISCAL)
        IF(ISSO2T(ISCAL).NE.ISNO2T(IPHAS))THEN
          WRITE(NFECRA,2133) ISCAL,ISSO2T(ISCAL),IPHAS,ISNO2T(IPHAS)
        ENDIF
        IF(IVSEXT(ISCAL).NE.IVIEXT(IPHAS))THEN
          WRITE(NFECRA,2134) ISCAL,IVSEXT(ISCAL),IPHAS,IVIEXT(IPHAS)
        ENDIF
      ENDDO
C
C     Test du theta de la diffusivite des scalaires et de Cp : ils doivent etre
C       variables en (en espace) si on les extrapole (en temps) (...)
      DO IPHAS = 1, NPHAS
        IF( ICPEXT(IPHAS).GT.0 .AND. ICP(IPHAS).LE.0 ) THEN
          WRITE(NFECRA,2135) IPHAS, ICPEXT(IPHAS), ICP(IPHAS)
          IOK = IOK + 1
        ENDIF
      ENDDO
      DO ISCAL = 1, NSCAL
        IF( IVSEXT(ISCAL).GT.0 .AND. IVISLS(ISCAL).LE.0 ) THEN
          WRITE(NFECRA,2136) ISCAL, IVSEXT(ISCAL), IVISLS(ISCAL)
          IOK = IOK + 1
        ENDIF
      ENDDO

C
C
C     Pour les tests suivants : Utilise-t-on un estimateur d'erreur ?
      INDEST = 0
      DO IPHAS = 1, NPHAS
        DO IEST = 1, NESTMX
          IIESCA = IESCAL(IEST,IPHAS)
          IF(IIESCA.GT.0) THEN
            INDEST = 1
          ENDIF
        ENDDO
      ENDDO
C
C     Estimateurs incompatibles avec calcul a champ de vitesse
C       fige (on ne fait rien, sauf ecrire des betises dans le listing)
       IF(INDEST.EQ.1.AND.ICCVFG.EQ.1) THEN
         WRITE(NFECRA,2137)
         IOK = IOK + 1
       ENDIF
C
C     A priori, pour le moment, l'ordre 2 en temps
C       (rho, visc, termes sources N.S, theta vitesse)
C     est incompatible avec
C       - estimateurs
C       - ipucou
C       - iphydr et icalhy
C       - dt variable en espace ou en temps et stationnaire
C     Ici on s'arrete si on n'est pas dans le cas du schema std
      DO IPHAS = 1, NPHAS
        IUIPH = IU(IPHAS)
        IVIPH = IV(IPHAS)
        IWIPH = IW(IPHAS)
        IF( (ABS(THETAV(IUIPH)-1.0D0).GT.1.D-3).OR.
     &      (ABS(THETAV(IVIPH)-1.0D0).GT.1.D-3).OR.
     &      (ABS(THETAV(IWIPH)-1.0D0).GT.1.D-3).OR.
     &      (    THETSN(IPHAS)       .GT.0.D0 ).OR.
     &      (    ISNO2T(IPHAS)       .GT.0    ).OR.
     &      (    THETRO(IPHAS)       .GT.0.D0 ).OR.
     &      (    IROEXT(IPHAS)       .GT.0    ).OR.
     &      (    THETVI(IPHAS)       .GT.0.D0 ).OR.
     &      (    IVIEXT(IPHAS)       .GT.0    )    ) THEN
           IF(INDEST.EQ.1.OR.IPUCOU.EQ.1.OR.
     &        IPHYDR.EQ.1.OR.ICALHY.EQ.1.OR.
     &        IDTVAR.EQ.1.OR.IDTVAR.EQ.2.OR.IDTVAR.LT.0) THEN
             WRITE(NFECRA,2140) IPHAS,
     &            THETAV(IUIPH),THETAV(IVIPH),THETAV(IWIPH),
     &            ISNO2T(IPHAS),THETSN(IPHAS),
     &            IROEXT(IPHAS),THETRO(IPHAS),
     &            IVIEXT(IPHAS),THETVI(IPHAS)
             IOK = IOK + 1
           ENDIF
         ENDIF
       ENDDO
C
C
C     Iterations sur navsto
C     Doit etre un entier superieur ou egal a 1
C     Pour le moment, on interdit NTERUP > 1 et
C       estimateurs, matrices poids, pression hydrostatique
C       et algo stationnaire
C
      IF(NTERUP.LE.0) THEN
        WRITE(NFECRA,3100) 'NTERUP', NTERUP
        IOK = IOK + 1
      ENDIF
C
      IF(NTERUP.GT.1) THEN
C
        IF(IPUCOU.EQ.1.OR.INDEST.EQ.1.OR.
     &       IPHYDR.EQ.1.OR.ICALHY.EQ.1.OR.IDTVAR.EQ.-1) THEN
          WRITE(NFECRA,2141) NTERUP
          IOK = IOK + 1
        ENDIF
C
      ENDIF
C
C     A priori, pour le moment, l'ordre 2 en temps
C       n'est pas pris en compte en k-eps, v2f ou k-omega couple : on s'arrete
      DO IPHAS = 1, NPHAS
        IF (ITYTUR(IPHAS).EQ.2 .AND.IKECOU(IPHAS).EQ.1) THEN
          IF((    THETST(IPHAS)       .GT.0.D0 ).OR.
     &       (    ISTO2T(IPHAS)       .GT.0    ).OR.
     &       (ABS(THETAV(IK (IPHAS))-1.0D0).GT.EPZERO).OR.
     &       (ABS(THETAV(IEP(IPHAS))-1.0D0).GT.EPZERO) ) THEN
            WRITE(NFECRA,2142)IPHAS,ITURB(IPHAS),IKECOU(IPHAS),
     &           THETST(IPHAS),ISTO2T(IPHAS),
     &           THETAV(IK (IPHAS)),THETAV(IEP(IPHAS))
            IOK = IOK + 1
          ENDIF
        ENDIF
        IF (ITURB(IPHAS).EQ.50.AND.IKECOU(IPHAS).EQ.1) THEN
          IF((    THETST(IPHAS)       .GT.0.D0 ).OR.
     &       (    ISTO2T(IPHAS)       .GT.0    ).OR.
     &       (ABS(THETAV(IK  (IPHAS))-1.0D0).GT.EPZERO).OR.
     &       (ABS(THETAV(IEP (IPHAS))-1.0D0).GT.EPZERO).OR.
     &       (ABS(THETAV(IPHI(IPHAS))-1.0D0).GT.EPZERO).OR.
     &       (ABS(THETAV(IFB (IPHAS))-1.0D0).GT.EPZERO) ) THEN
            WRITE(NFECRA,2143)IPHAS,ITURB(IPHAS),IKECOU(IPHAS),
     &           THETST(IPHAS),ISTO2T(IPHAS),
     &           THETAV(IK  (IPHAS)),THETAV(IEP (IPHAS)),
     &           THETAV(IPHI(IPHAS)),THETAV(IFB (IPHAS))
            IOK = IOK + 1
          ENDIF
        ENDIF
        IF (ITURB(IPHAS).EQ.60.AND.IKECOU(IPHAS).EQ.1) THEN
          IF((    THETST(IPHAS)       .GT.0.D0 ).OR.
     &       (    ISTO2T(IPHAS)       .GT.0    ).OR.
     &       (ABS(THETAV(IK  (IPHAS))-1.0D0).GT.EPZERO).OR.
     &       (ABS(THETAV(IOMG(IPHAS))-1.0D0).GT.EPZERO) ) THEN
            WRITE(NFECRA,2144)IPHAS,ITURB(IPHAS),IKECOU(IPHAS),
     &           THETST(IPHAS),ISTO2T(IPHAS),
     &           THETAV(IK  (IPHAS)),THETAV(IOMG(IPHAS))
            IOK = IOK + 1
          ENDIF
        ENDIF
      ENDDO
C
C     A priori, pour le moment, l'ordre 2 en temps
C       (rho, visc, cp, termes sources N.S., Turb., Scal., theta)
C     est incompatible avec les physiques particulieres
C     Ici on s'arrete si on n'est pas dans le cas du schema std
C
      IPPM = 0
      DO IPP = 1, NMODMX
        IF(IPPMOD(IPP).NE.-1) THEN
          IPPM = 1
        ENDIF
      ENDDO
C
      IF(IPPM.GT.0) THEN
        ISTOP = 0
        DO IVAR = 1, NVAR
          IF( (ABS(THETAV(IVAR)-1.0D0).GT.1.D-3) ) ISTOP = 1
        ENDDO
        DO IPHAS = 1, NPHAS
          IF((    THETSN(IPHAS)       .GT.0.D0 ).OR.
     &       (    ISNO2T(IPHAS)       .GT.0    ).OR.
     &       (    THETRO(IPHAS)       .GT.0.D0 ).OR.
     &       (    IROEXT(IPHAS)       .GT.0    ).OR.
     &       (    THETVI(IPHAS)       .GT.0.D0 ).OR.
     &       (    IVIEXT(IPHAS)       .GT.0    ).OR.
     &       (    THETCP(IPHAS)       .GT.0.D0 ).OR.
     &       (    ICPEXT(IPHAS)       .GT.0    )    ) ISTOP = 1
        ENDDO
        DO ISCAL = 1, NSCAL
          IF((    THETSS(ISCAL)       .GT.0.D0 ).OR.
     &       (    ISSO2T(ISCAL)       .GT.0    ).OR.
     &       (    THETVS(IPHAS)       .GT.0.D0 ).OR.
     &       (    IVSEXT(IPHAS)       .GT.0    )    ) ISTOP = 1
        ENDDO
C
        IF(ISTOP.NE.0) THEN
          WRITE(NFECRA,2145)
          IOK = IOK + 1
        ENDIF
      ENDIF
C
C     A priori, pour le moment, l'ordre 2 en temps
C       n'est pas pris en compte pour les termes issus du Lagrangien.
C       On pourrait le signaler et continuer : on s'arrete.
      IF(IILAGR .EQ. 2) THEN
        DO IPHAS = 1, NPHAS
          IF((    THETSN(IPHAS)       .GT.0.D0 ).OR.
     &      (    ISNO2T(IPHAS)       .GT.0    ).OR.
     &       (    THETST(IPHAS)       .GT.0.D0 ).OR.
     &       (    ISTO2T(IPHAS)       .GT.0    ) ) THEN
            WRITE(NFECRA,2146)IPHAS,THETSN(IPHAS),ISNO2T(IPHAS),
     &                              THETST(IPHAS),ISTO2T(IPHAS)
            IOK = IOK + 1
          ENDIF
        ENDDO
        DO ISCAL = 1, NSCAL
          IF (ISCSTH(ISCAL).EQ.1.OR.ISCSTH(ISCAL).EQ.2) THEN
            IF((    THETSS(ISCAL)       .GT.0.D0 ).OR.
     &         (    ISSO2T(ISCAL)       .GT.0    )) THEN
              WRITE(NFECRA,2147)
     &    'lagrangien ',ISCAL,THETSS(ISCAL),ISSO2T(ISCAL),'uslag1'
              IOK = IOK + 1
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C     A priori, pour le moment, l'ordre 2 en temps
C       n'est pas pris en compte pour les termes issus du rayonnement.
C       On pourrait le signaler et continuer : on s'arrete.
      IF (IIRAYO.GT.0) THEN
        DO IPHAS = 1, NPHAS
          DO ISCAL = 1, NSCAL
            IF (ISCAL.EQ.ISCALT(IPHAS)) THEN
              IF((    THETSS(ISCAL)       .GT.0.D0 ).OR.
     &           (    ISSO2T(ISCAL)       .GT.0    )) THEN
                WRITE(NFECRA,2147)
     &         'rayonnement',ISCAL,THETSS(ISCAL),ISSO2T(ISCAL),'usray1'
                IOK = IOK + 1
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
C --- Algorithme stationnaire
      IF (IDTVAR.LT.0) THEN
        DO IPP = 2, NVPPMX
          II = ITRSVR(IPP)
          IF (II.GE.1) THEN
            IF (RELAXV(II).GT.1D0.OR.RELAXV(II).LT.0D0) THEN
              CHAINE=NOMVAR(IPP)
              WRITE(NFECRA,2148) CHAINE(1:8),RELAXV(II)
              IOK = IOK + 1
            ENDIF
          ENDIF
        ENDDO
        DO IPHAS = 1, NPHAS
          IF((RELAXV(IV(IPHAS)).NE.RELAXV(IU(IPHAS)))
     &   .OR.(RELAXV(IW(IPHAS)).NE.RELAXV(IU(IPHAS))) ) THEN
            WRITE(NFECRA,2149) RELAXV(IU(IPHAS)),RELAXV(IV(IPHAS)),
     &           RELAXV(IW(IPHAS))
            IOK = IOK + 1
          ENDIF
        ENDDO
C       L'algorithme stationnaire n'est pas compatible avec le module Lagrangien
        IF (IILAGR.NE.0) THEN
          WRITE(NFECRA,2150) IILAGR
          IOK = IOK + 1
        ENDIF
C       L'algorithme stationnaire n'est pas compatible avec la LES
        DO IPHAS = 1, NPHAS
          IF (ITYTUR(IPHAS).EQ.4) THEN
            WRITE(NFECRA,2151) ITURB(IPHAS)
            IOK = IOK + 1
          ENDIF
        ENDDO
      ENDIF
C
C --- Reconstruction des gradients
      IF(IMRGRA.NE.0.AND.IMRGRA.NE.1.AND.
     &   IMRGRA.NE.2.AND.IMRGRA.NE.3.AND.
     &   IMRGRA.NE.4 ) THEN
        WRITE(NFECRA,2205) 'IMRGRA',IMRGRA
        IOK = IOK + 1
      ENDIF
C
C On verifie l'angle de non orthogonalite de selection du
C   voisinage etendu dans le cas du moindre carre qui l'utilise
C
      IF(IMRGRA.EQ.3) THEN
        IF(ANOMAX.GT.PI*0.5D0.OR.ANOMAX.LT.0.D0) THEN
          WRITE(NFECRA,2206) ANOMAX, IMRGRA
        ENDIF
      ENDIF
C
C Extrapolation : indetermination possible par mc,
C     necessitant un traitement particulier dans gradmc,
C     pour lequel on fait certaines hypotheses
      IF(IMRGRA.EQ.1.OR.IMRGRA.EQ.2.OR.IMRGRA.EQ.3) THEN
        IF((NPHAS.GT.1).OR.
     &     ((ABS(EXTRAG(IPR(1))-1.D0).GT.EPZERO).AND.
     &      (ABS(EXTRAG(IPR(1))     ).GT.EPZERO)  )) THEN
          WRITE(NFECRA,2207) IMRGRA, NPHAS, EXTRAG(IPR(1))
          IOK = IOK + 1
        ENDIF
      ENDIF
C
C Les nombres de sweeps n'ont pas a etre verifies :
C  ce sont simplement des entiers (negatifs si on veut etre sur de ne
C  *jamais* entrer dans les boucles)
C
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          IF(IMLIGR(II).GT.1) THEN
            CHAINE=NOMVAR(IPP)
            WRITE(NFECRA,2300) CHAINE(1:8),II,IMLIGR(II)
            IOK = IOK + 1
          ENDIF
        ENDIF
      ENDDO
C
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          IF(IRCFLU(II).NE.1.AND.IRCFLU(II).NE.0) THEN
            CHAINE=NOMVAR(IPP)
            WRITE(NFECRA,2310) CHAINE(1:8),II,IRCFLU(II)
            IOK = IOK + 1
          ENDIF
        ENDIF
      ENDDO
C Non reconstruction des flux en SOLU n'a pas de sens pour la convection
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          IF(IRCFLU(II).EQ.0.AND.ISCHCV(II).EQ.0.AND.
     &                           BLENCV(II).NE.0.D0) THEN
            CHAINE=NOMVAR(IPP)
            WRITE(NFECRA,2311) CHAINE(1:8),II,IRCFLU(II),II,ISCHCV(II)
            IOK = IOK + 1
          ENDIF
        ENDIF
      ENDDO
C
C Il n'y a pas besoin de test sur les epsilons
C   Ce sont simplement des reels
C   Une valeur negative indique qu'on veut atteindre
C   le nombre d'iterations maximal
C
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          IF(CLIMGR(II).LT.1.D0) THEN
            CHAINE=NOMVAR(IPP)
            WRITE(NFECRA,2320) CHAINE(1:8),II,CLIMGR(II)
            IOK = IOK + 1
          ENDIF
        ENDIF
      ENDDO
C
C
C EXTRAG non nul permis uniquement pour la pression.
C        et dans ce cas egal a 1
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          IF(ABS(EXTRAG(II)     ).GE.EPZERO) THEN
            IOKPRE = 0
            DO IPHAS = 1, NPHAS
              IF (II.EQ.IPR(IPHAS)) THEN
                IOKPRE = 1
                IF(ABS(EXTRAG(II)-1.D0).GE.EPZERO) THEN
                  CHAINE=NOMVAR(IPP)
                  WRITE(NFECRA,2330) CHAINE(1:8),II,EXTRAG(II)
                  IOK = IOK + 1
                ENDIF
              ENDIF
            ENDDO
            IF(IOKPRE.EQ.0) THEN
              CHAINE=NOMVAR(IPP)
              WRITE(NFECRA,2331) CHAINE(1:8),II,EXTRAG(II)
              IOK = IOK + 1
            ENDIF
          ENDIF
        ENDIF
      ENDDO
C
C
C --- Solveurs iteratifs
C
C Il n'y a pas besoin de test sur les epsilons
C   Ce sont simplement des reels
C   Une valeur negative indique qu'on veut atteindre
C   le nombre d'iterations maximal
C Il n'y a pas besoin de test sur le nombre d'iterations
C   Ce sont simplement des entiers
C   Une valeur negative indique au'on veut sortir de suite
C
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          IF(IRESOL(II).NE.-1) THEN
            IRESOP = MOD(IRESOL(II),1000)
            IPOLOP = (IRESOL(II)-IRESOP)/1000
            IF ((IRESOP.NE.0.AND.
     &           IRESOP.NE.1.AND.IRESOP.NE.2).OR.
     &          (IRESOP.EQ.1.AND.IPOLOP.NE.0)) THEN
              CHAINE=NOMVAR(IPP)
              WRITE(NFECRA,2400) CHAINE(1:8),II,IRESOL(II)
              IOK = IOK + 1
            ENDIF
          ENDIF
          IF (IDIRCL(II).NE.0.AND.IDIRCL(II).NE.1) THEN
            CHAINE=NOMVAR(IPP)
            WRITE(NFECRA,2401) CHAINE(1:8),II,IDIRCL(II)
            IOK = IOK + 1
          ENDIF
        ENDIF
      ENDDO
C
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          IF ((IRESOL(II).EQ.0.AND.ICONV(II).EQ.1).OR.
     &        (IRESOL(II).EQ.1.AND.ICONV(II).EQ.0)) THEN
            CHAINE=NOMVAR(IPP)
            WRITE(NFECRA,2410) CHAINE(1:8),II,IRESOL(II),ICONV(II)
          ENDIF
        ENDIF
      ENDDO
C
C --- Le multigrille sera verifie d'un seul bloc plus bas.
C
C
C --- Suite de calcul
C
      IF(ISUITE.NE.0.AND.ISUITE.NE.1) THEN
        WRITE(NFECRA,2200) 'ISUITE',ISUITE
        IOK = IOK + 1
      ENDIF
      IF(ILEAUX.NE.0.AND.ILEAUX.NE.1) THEN
        WRITE(NFECRA,2200) 'ILEAUX',ILEAUX
        IOK = IOK + 1
      ENDIF
      IF(IECAUX.NE.0.AND.IECAUX.NE.1) THEN
        WRITE(NFECRA,2200) 'IECAUX',IECAUX
        IOK = IOK + 1
      ENDIF
C En LES, on previent que ce n'est pas malin de ne pas relire le fichier
C   auxiliaire
      IF(IECAUX.EQ.0.OR.ILEAUX.EQ.0) THEN
        DO IPHAS = 1, NPHAS
          IF(ITYTUR(IPHAS).EQ.4) THEN
            WRITE(NFECRA,2420) IPHAS,ITURB(IPHAS),ILEAUX,IECAUX
          ENDIF
        ENDDO
      ENDIF
C
C --- Reperage du temps et marche en temps
C
C Le nombre de pas de temps pourrait bien etre negatif : pas de test.
C
      IF(INPDT0.NE.0.AND.INPDT0.NE.1) THEN
        WRITE(NFECRA,2200) 'INPDT0',INPDT0
        IOK = IOK + 1
      ENDIF
C
      IF(IDTVAR.LT.-1.OR.IDTVAR.GT.2) THEN
        WRITE(NFECRA,2500) 'IDTVAR',IDTVAR
        IOK = IOK + 1
      ENDIF
C
      IF(IDTVAR.GT.0.AND.VARRDT.LT.0.D0) THEN
        WRITE(NFECRA,2510)'VARRDT', VARRDT
        IOK = IOK + 1
      ENDIF
C
      IF(DTREF .LT.0.D0) THEN
        WRITE(NFECRA,2510)'DTREF ', DTREF
        IOK = IOK + 1
      ENDIF
C
      IF(DTMIN.LE.0.D0   .OR. DTMAX.LE.0.D0 .OR.
     &   DTMIN.GT.DTMAX                       ) THEN
        WRITE(NFECRA,2520) DTMIN, DTMAX
        IF(IDTVAR.GT.0) THEN
          IOK = IOK + 1
        ENDIF
      ENDIF
C
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          IF(CDTVAR(II).LE.0.D0) THEN
            CHAINE=NOMVAR(IPP)
            WRITE(NFECRA,2530) CHAINE(1:8),II,CDTVAR(II)
            IOK = IOK + 1
          ENDIF
        ENDIF
      ENDDO
C
      IF(IPTLRO.NE.0.AND.IPTLRO.NE.1) THEN
        WRITE(NFECRA,2200) 'IPTLRO',IPTLRO
        IOK = IOK + 1
      ENDIF
C
C     Si on est en pas de temps constant, on ne touche pas le pas de temps,
C     mais on indique juste les depassements de critere local.
      IF(IDTVAR.EQ.0.AND.IPTLRO.EQ.1) THEN
        WRITE(NFECRA,2540) IPTLRO,IDTVAR
      ENDIF
      IF(IDTVAR.EQ.-1.AND.IPTLRO.EQ.1) THEN
        WRITE(NFECRA,2541) IPTLRO,IDTVAR
      ENDIF
C
C --- Turbulence
C
C    Modele
C
      DO IPHAS = 1, NPHAS
        ITRBPH = ITURB(IPHAS)
        IF ( ITRBPH.NE. 0.AND.ITRBPH.NE.10.AND.ITRBPH.NE.20.AND.
     &       ITRBPH.NE.21.AND.ITRBPH.NE.30.AND.ITRBPH.NE.31.AND.
     &       ITRBPH.NE.40.AND.ITRBPH.NE.41.AND.ITRBPH.NE.50.AND.
     &       ITRBPH.NE.60) THEN
          WRITE(NFECRA,2600) IPHAS,'ITURB  ',ITRBPH
          IOK = IOK + 1
        ENDIF
      ENDDO
C
C     Methode des vortex pour la LES
C
      IF (IVRTEX.NE.0 .AND.IVRTEX.NE.1) THEN
        WRITE(NFECRA,2200) 'IVRTEX ',IVRTEX
        IOK = IOK + 1
      ENDIF
C On impose qu'il n'y ait qu'une seule phase
      IF (IVRTEX.EQ.1 .AND. NPHAS.GT.1) THEN
        WRITE(NFECRA,2605)NPHAS,IVRTEX
        IOK = IOK + 1
      ENDIF
      IPHAS = 1
      IF(IVRTEX.EQ.1.AND.ITYTUR(IPHAS).NE.4) THEN
        WRITE(NFECRA,2606)ITYTUR,IVRTEX
        IVRTEX = 0
      ENDIF
C
C    Nb de variables
C
      IF(NSCAL.GE.1) THEN
        DO IPHAS = 1, NPHAS
          IF(ISCALT(IPHAS).GT.NSCAL) THEN
            WRITE(NFECRA,2610)IPHAS,
     &                 'NUMERO DU SCALAIRE TEMPERATURE ',ISCALT(IPHAS),
     &                 'NOMBRE DE SCALAIRES            ',NSCAL
            IOK = IOK + 1
          ENDIF
          IF(  (NVAR.LT. 4+NSCAL               ) .OR.
     &         (NVAR.LT. 6+NSCAL.AND.ITYTUR(IPHAS).EQ.2).OR.
     &         (NVAR.LT.11+NSCAL.AND.ITYTUR(IPHAS).EQ.3).OR.
     &         (NVAR.LT. 8+NSCAL.AND.ITURB(IPHAS).EQ.50).OR.
     &         (NVAR.LT. 6+NSCAL.AND.ITURB(IPHAS).EQ.60)    ) THEN
            WRITE(NFECRA,2610)IPHAS,
     &                 'NOMBRE DE VARIABLES            ',NVAR,
     &                 'NOMBRE DE SCALAIRES            ',NSCAL
            IOK = IOK + 1
          ENDIF
        ENDDO
      ENDIF
C
      DO IPHAS = 1, NPHAS
C
        IF(IDEUCH(IPHAS).LT.0.OR.IDEUCH(IPHAS).GT.2) THEN
          WRITE(NFECRA,2211)IPHAS, 'IDEUCH',IDEUCH(IPHAS)
          IOK = IOK + 1
        ENDIF
        IF (IDEUCH(IPHAS).NE.0 .AND.
     &       (ITURB(IPHAS).EQ.0 .OR. ITURB(IPHAS).EQ.10 .OR.
     &       ITYTUR(IPHAS).EQ.4) ) THEN
           WRITE(NFECRA,2209)IPHAS,ITURB(IPHAS),IDEUCH(IPHAS)
           IOK = IOK + 1
        ENDIF
        IF(ILOGPO(IPHAS).NE.0.AND.ILOGPO(IPHAS).NE.1) THEN
          WRITE(NFECRA,2201)IPHAS, 'ILOGPO',ILOGPO(IPHAS)
          IOK = IOK + 1
        ENDIF
C
C      Specifique k-epsilon, v2f et k-omega
C
       IF(ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50
     &       .OR. ITURB(IPHAS).EQ.60) THEN
          IF( (NVAR.LE.5.AND.ITYTUR(IPHAS).EQ.2) .OR.
     &        (NVAR.LE.7.AND.ITURB(IPHAS).EQ.50) .OR.
     &        (NVAR.LE.5.AND.ITURB(IPHAS).EQ.60) ) THEN
            WRITE(NFECRA,2610)IPHAS,
     &                 'NOMBRE DE VARIABLES            ',NVAR,
     &                 'OPTION POUR LA TURBULENCE      ',ITURB(IPHAS)
            IOK = IOK + 1
          ENDIF
C     Le choix de ICLKEP n'est possible qu'en k-eps ou v2f
          IF (ITURB(IPHAS).NE.60) THEN
            IF(ICLKEP(IPHAS).NE.0.AND.ICLKEP(IPHAS).NE.1) THEN
              WRITE(NFECRA,2201) IPHAS,'ICLKEP',ICLKEP(IPHAS)
              IOK = IOK + 1
            ENDIF
          ENDIF
          IF(IKECOU(IPHAS).NE.0.AND.IKECOU(IPHAS).NE.1) THEN
            WRITE(NFECRA,2201)IPHAS, 'IKECOU',IKECOU(IPHAS)
            IOK = IOK + 1
          ENDIF
C     En k-eps a prod lin et en v2f on force IKECOU a 0
          IF (IKECOU(IPHAS).EQ.1 .AND.
     &         (ITURB(IPHAS).EQ.21 .OR. ITURB(IPHAS).EQ.50)) THEN
            WRITE(NFECRA,2208)IPHAS,ITURB(IPHAS),IKECOU(IPHAS)
            IOK = IOK + 1
          ENDIF
C     En stationnaire on force IKECOU a 0
          IF (IKECOU(IPHAS).NE.0.AND.IDTVAR.LT.0) THEN
            WRITE(NFECRA,2210)IPHAS,IKECOU(IPHAS)
            IOK = IOK + 1
          ENDIF

          IF(IGRHOK(IPHAS).NE.0.AND.IGRHOK(IPHAS).NE.1) THEN
            WRITE(NFECRA,2201)IPHAS, 'IGRHOK',IGRHOK(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF(IGRAKE(IPHAS).NE.0.AND.IGRAKE(IPHAS).NE.1) THEN
            WRITE(NFECRA,2201)IPHAS, 'IGRAKE',IGRAKE(IPHAS)
            IOK = IOK + 1
          ENDIF
c        IF( IGRAKE.EQ.1.AND.(GX**2+GY**2+GZ**2).LE.EPZERO**2 ) THEN
c          WRITE(NFECRA,2620)'IGRAKE',IGRAKE,GX,GY,GZ
c          IOK = IOK + 1
c        ENDIF
          IF(NSCAL.GT.0) THEN
            IF(ISCALT(IPHAS).LE.0.AND.
     &           (GX**2+GY**2+GZ**2).GE.EPZERO**2) THEN
              WRITE(NFECRA,2621) IPHAS,GX,GY,GZ,ISCALT(IPHAS)
              IF(IGRAKE(IPHAS).EQ.1) THEN
                WRITE(NFECRA,2622)IPHAS, 'IGRAKE',IGRAKE(IPHAS)
              ENDIF
CMO            IOK = IOK + 1
            ENDIF
          ENDIF
C
C     Si RELAXV(IK) a ete modifie par l'utilisateur mais que IKECOU n'est
C     pas egal a 0, on previent que ce sera sans effet
C     Sinon on verifie que RELAXV(IK) est bien compris entre 0 et 1
C     (en stationnaire cela a deja ete fait plus haut)
          IKIPH = IK(IPHAS)
          IF (ITYTUR(IPHAS).EQ.6) THEN
            IEIPH = IOMG(IPHAS)
          ELSE
            IEIPH = IEP(IPHAS)
          ENDIF
          IF ( (ABS(RELAXV(IKIPH)+999.D0).GT.EPZERO .OR.
     &          ABS(RELAXV(IEIPH)+999.D0).GT.EPZERO ) .AND.
     &          IKECOU(IPHAS).NE.0) WRITE(NFECRA,2623) IPHAS,
     &            RELAXV(IKIPH),RELAXV(IEIPH)
          IF (IKECOU(IPHAS).EQ.0 .AND. IDTVAR.GE.0) THEN
            IF(RELAXV(IKIPH).GT.1.D0.OR.RELAXV(IKIPH).LT.0.D0 .OR.
     &         RELAXV(IEIPH).GT.1.D0.OR.RELAXV(IEIPH).LT.0.D0) THEN
              WRITE(NFECRA,2624) IPHAS,RELAXV(IKIPH),RELAXV(IEIPH)
              IOK = IOK + 1
            ENDIF
          ENDIF
C
        ENDIF
C
C     Specifique Rij-epsilon
C
        IF(ITYTUR(IPHAS).EQ.3) THEN
          IF(NVAR.LE.10) THEN
            WRITE(NFECRA,2610)IPHAS,
     &                 'NOMBRE DE VARIABLES            ',NVAR,
     &                 'OPTION POUR LA TURBULENCE      ',ITURB(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF(IRIJNU(IPHAS).NE.0.AND.IRIJNU(IPHAS).NE.1) THEN
            WRITE(NFECRA,2201)IPHAS,'IRIJNU',IRIJNU(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF(IRIJRB(IPHAS).NE.0.AND.IRIJRB(IPHAS).NE.1) THEN
            WRITE(NFECRA,2201)IPHAS, 'IRIJRB',IRIJRB(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF (ITURB(IPHAS).EQ.30) THEN
C     echo de paroi et implicitation speciale de la diffusion de epsilon
C     seulement en Rij standard
            IF(IRIJEC(IPHAS).NE.0.AND.IRIJEC(IPHAS).NE.1) THEN
              WRITE(NFECRA,2201)IPHAS, 'IRIJEC',IRIJEC(IPHAS)
              IOK = IOK + 1
            ENDIF
            IF(IDIFRE(IPHAS).NE.0.AND.IDIFRE(IPHAS).NE.1) THEN
              WRITE(NFECRA,2201)IPHAS, 'IDIFRE',IDIFRE(IPHAS)
              IOK = IOK + 1
            ENDIF
          ENDIF
          IF(IGRARI(IPHAS).NE.0.AND.IGRARI(IPHAS).NE.1) THEN
            WRITE(NFECRA,2201)IPHAS, 'IGRARI',IGRARI(IPHAS)
            IOK = IOK + 1
          ENDIF
c        IF( IGRARI.EQ.1.AND.(GX**2+GY**2+GZ**2).LE.EPZERO**2 ) THEN
c          WRITE(NFECRA,2620)'IGRARI',IGRARI,GX,GY,GZ
c          IOK = IOK + 1
c        ENDIF
          IF(ICLSYR(IPHAS).NE.0.AND.ICLSYR(IPHAS).NE.1) THEN
            WRITE(NFECRA,2201) IPHAS,'ICLSYR',ICLSYR(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF(ICLPTR(IPHAS).NE.0.AND.ICLPTR(IPHAS).NE.1) THEN
            WRITE(NFECRA,2201)IPHAS, 'ICLPTR',ICLPTR(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF(NSCAL.GT.0) THEN
            IF(ISCALT(IPHAS).LE.0.AND.
     &           (GX**2+GY**2+GZ**2).GE.EPZERO**2) THEN
              WRITE(NFECRA,2621)IPHAS,GX,GY,GZ,ISCALT(IPHAS)
              IF(IGRARI(IPHAS).EQ.1) THEN
                WRITE(NFECRA,2622)IPHAS, 'IGRARI',IGRARI(IPHAS)
              ENDIF
CMO            IOK = IOK + 1
            ENDIF
          ENDIF
        ENDIF
C
C     Specifique LES
C
        IF(ITYTUR(IPHAS).EQ.4) THEN
          IF(IDRIES(IPHAS).NE.1.AND.IDRIES(IPHAS).NE.0) THEN
            WRITE(NFECRA,2201) IPHAS,'IDRIES',IDRIES(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF(IDRIES(IPHAS).NE.0.AND.ITURB(IPHAS).EQ.41) THEN
            WRITE(NFECRA,2630) IPHAS,IDRIES(IPHAS),ITURB(IPHAS)
            IOK = IOK + 1
          ENDIF
C         La reduction du voisinage etendu peut degrader
C         les resultats du modele dynamique en LES
          IF(ITURB(IPHAS).EQ.41.AND.IMRGRA.EQ.3) THEN
            WRITE(NFECRA,2607) IPHAS, ITURB(IPHAS), IMRGRA
          ENDIF
        ENDIF
C
      ENDDO
C
C --- Stokes
C
C
      IF(IPRCO .NE.0.AND.IPRCO .NE.1) THEN
        WRITE(NFECRA,2200) 'IPRCO ',IPRCO
        IOK = IOK + 1
      ENDIF
      IF(IPRCO.EQ.1) THEN
        DO IPHAS = 1, NPHAS
          IF(IREVMC(IPHAS).NE.0.AND.IREVMC(IPHAS).NE.1.AND.
     &                              IREVMC(IPHAS).NE.2) THEN
            WRITE(NFECRA,2211) IPHAS,'IREVMC',IREVMC(IPHAS)
            IOK = IOK + 1
          ENDIF
          ARAKFR = ARAK(IPHAS)
          IF (IDTVAR.LT.0) ARAKFR=ARAKFR*RELAXV(IU(IPHAS))
          IF(ARAKFR.GT.1.D0 .OR. ARAKFR.LT.0.D0) THEN
            WRITE(NFECRA,2640) IPHAS,'ARAK  ',ARAKFR
            IOK = IOK + 1
          ENDIF
          IF( RELAXV(IPR(IPHAS)).GT.1D0 .OR.
     &        RELAXV(IPR(IPHAS)).LT.0D0     ) THEN
            WRITE(NFECRA,2625) IPHAS,RELAXV(IPR(IPHAS))
            IOK = IOK + 1
          ENDIF
        ENDDO
      ENDIF
C
C --- Couplage U-P
C
      IF (IPUCOU.NE.0.AND.IPUCOU.NE.1) THEN
        WRITE(NFECRA,2200) 'IPUCOU ',IPUCOU
        IOK = IOK + 1
      ENDIF
C
C Incompatibilite pour le moment des matrices poids avec un theta schema
C Si theta n'est pas egal a 1 pour la vitesse (incompatibilite du
C pas de temps variable aussi)
C
      DO IPHAS = 1, NPHAS
        JJ = IU(IPHAS)
        IF((ABS(THETAV(JJ)-1.0D0).GT.EPZERO).AND.
     &    ((IDTVAR.NE.0).OR.(IPUCOU.EQ.1))) THEN
          WRITE(NFECRA,2204) THETAV(JJ),IDTVAR,IPUCOU
        ENDIF
      ENDDO
C
C --- Prise en compte de la pression hydrostatique
C
C     IPHYDR et ICALHY sont testes dans varpos
C       (valeur 0 ou 1)
C
      IF (IPHYDR.EQ.1.AND.NPHAS.GT.1) THEN
        WRITE(NFECRA,2202) IPHYDR,NPHAS
        IOK = IOK +1
      ENDIF
C
C --- Champ de vitesse fige
C
      IF (ICCVFG.NE.0.AND.ICCVFG.NE.1) THEN
        WRITE(NFECRA,2200) 'ICCVFG ',ICCVFG
        IOK = IOK + 1
      ENDIF
C
C --- Interpolation face des viscosites
C
      IF(IMVISF.NE.0.AND.IMVISF.NE.1) THEN
        WRITE(NFECRA,2200) 'IMVISF',IMVISF
        IOK = IOK + 1
      ENDIF
C
C --- Traitement de la temperature pour couplage SYRTHES
C     Verification du nombre de couplages
C
      IF(ITBRRB.NE.0 .AND. ITBRRB.NE.1       ) THEN
        WRITE(NFECRA,2200) 'ITBRRB',ITBRRB
        IOK = IOK + 1
      ENDIF
C
C     On regarde si ICPSYR a des valeurs realistes
      IF(NSCAL.GT.0) THEN
        DO ISCAL = 1, NSCAL
          IF(ICPSYR(ISCAL).NE.0.AND.ICPSYR(ISCAL).NE.1) THEN
            CHAINE=NOMVAR(IPPRTP(ISCA(ISCAL)))
            WRITE(NFECRA,2650)CHAINE(1:8),'ICPSYR',ISCAL,ICPSYR(ISCAL)
            IOK = IOK + 1
          ENDIF
        ENDDO
      ENDIF
C
C     On compte le nombre de scalaires couples
      NBSCCP = 0
      IF(NSCAL.GT.0) THEN
        DO ISCAL = 1, NSCAL
          NBSCCP = NBSCCP+ICPSYR(ISCAL)
        ENDDO
      ENDIF
C
C     On regarde s'il y a du couplage
C
      CALL NBCSYR (NBCCOU)
C     ===========
C
C     S'il n'y a pas de couplage
      IF(NBCCOU.EQ.0) THEN
C
C       et qu'il n'y a pas zero scalaire couple, on s'arrete
        IF(NBSCCP.NE.0) THEN
          WRITE(NFECRA,2660)NBSCCP,NSCAL
          IOK = IOK + 1
        ENDIF
C
C     Sinon, s'il y a du couplage
      ELSE
C
C       et qu'il n'y a pas un et un seul scalaire couple, on s'arrete
        IF(NBSCCP.NE.1) THEN
          WRITE(NFECRA,2661)NSCAL,NBSCCP
          IOK = IOK + 1
        ENDIF
C
C          que le scalaire couple n'est pas la temperature, on s'arrete
C          attention : tout est pret, mais il n'y a pas eu de valid.
C          en outre, ca permet de bien verifier que l'utilisateur
C          ne s'est pas trompe dans 99% des cas d'utilisation
C         (en compressible, on couple l'energie)
        DO ISCAL = 1, NSCAL
          IF(ICPSYR(ISCAL).EQ.1) THEN
            IF(IPPMOD(ICOMPF).LT.0) THEN
              IF(ABS(ISCSTH(ISCAL)).NE.1) THEN
                WRITE(NFECRA,2662)ISCAL,ISCSTH(ISCAL)
                IOK = IOK + 1
              ENDIF
            ELSE
              IF(ISCSTH(ISCAL).NE.3) THEN
                WRITE(NFECRA,2663)ISCAL,ISCAL,ISCSTH(ISCAL)
                IOK = IOK + 1
              ENDIF
            ENDIF
          ENDIF
        ENDDO
C
      ENDIF
C
C
C --- Estimateurs  d'erreur pour Navier-Stokes
C
      DO IPHAS = 1, NPHAS
        DO IEST = 1, NESTMX
          IIESCA = IESCAL(IEST,IPHAS)
          IF (IIESCA.NE.0.AND.IIESCA.NE.1.AND.IIESCA.NE.2) THEN
            WRITE(NFECRA,2664) IPHAS,IEST,IEST,IPHAS,IIESCA,
     &                         IESPRE,IESDER,IESCOR,IESTOT
            IOK = IOK + 1
          ENDIF
        ENDDO
      ENDDO
C
C
C --- Distance a la paroi
C
      IF(INEEDY.EQ.1) THEN
C
        IF(ABS(ICDPAR).NE.1.AND.ABS(ICDPAR).NE.2) THEN
          WRITE(NFECRA,2700) ICDPAR
          IOK = IOK + 1
        ENDIF
        IF(NITMAY.LT.1) THEN
          WRITE(NFECRA,3100) 'NITMAY',NITMAY
          IOK = IOK + 1
        ENDIF
        IF(IMLIGY.GT.1) THEN
          WRITE(NFECRA,2750) 'IMLIGY',IMLIGY
          IOK = IOK + 1
        ENDIF
        IF(IRCFLY.NE.1.AND.IRCFLY.NE.0) THEN
          WRITE(NFECRA,2200) 'IRCFLY',IRCFLY
          IOK = IOK + 1
        ENDIF
        IF(ISCHCY.NE.1.AND.ISCHCY.NE.0) THEN
          WRITE(NFECRA,2200) 'ISCHCY',ISCHCY
          IOK = IOK + 1
        ENDIF
        IF(ISSTPY.NE.1.AND.ISSTPY.NE.0) THEN
          WRITE(NFECRA,2200) 'ISSTPY',ISSTPY
          IOK = IOK + 1
        ENDIF
        IF(IMGRPY.NE.1.AND.IMGRPY.NE.0) THEN
          WRITE(NFECRA,2200) 'IMGRPY',IMGRPY
          IOK = IOK + 1
        ENDIF
        IF(NTCMXY.LT.1) THEN
          WRITE(NFECRA,3100) 'NTCMXY',NTCMXY
          IOK = IOK + 1
        ENDIF
C
        IF(BLENCY.GT.1.D0.OR.BLENCY.LT.0.D0) THEN
          WRITE(NFECRA,2710) 'BLENCY',BLENCY
          IOK = IOK + 1
        ENDIF
        IF(CLIMGY.LT.1.D0) THEN
          WRITE(NFECRA,2720) 'CLIMGY',CLIMGY
          IOK = IOK + 1
        ENDIF
        IF(ABS(EXTRAY-1.D0).GT.EPZERO.AND.ABS(EXTRAY).GT.EPZERO) THEN
          WRITE(NFECRA,2730) 'EXTRAY',EXTRAY
          IOK = IOK + 1
        ENDIF
        IF(COUMXY.LE.0.D0) THEN
          WRITE(NFECRA,2740) 'COUMXY',COUMXY
          IOK = IOK + 1
        ENDIF
        IF(YPLMXY.LE.0.D0) THEN
          WRITE(NFECRA,2740) 'YPLMXY',YPLMXY
          IOK = IOK + 1
        ENDIF
C
      ENDIF
C
C
C=======================================================================
C 2. MULTIGRILLE : TABLEAUX DU MULTIGRILLE : formats 3000
C=======================================================================
C
C --- Options generales
C
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          IF(IMGR(II).NE.0.AND.IMGR(II).NE.1) THEN
            CHAINE=NOMVAR(IPP)
            WRITE(NFECRA,3000) CHAINE(1:8),'IMGR  ',II,IMGR(II)
            IOK = IOK + 1
          ENDIF
          DO IPHAS = 1, NPHAS
            IF(IMGR(II).EQ.1.AND.II.NE.IPR(IPHAS)) THEN
              WRITE(NFECRA,3001)
              IOK = IOK + 1
            ENDIF
          ENDDO
          IF(IMGR(II).EQ.1) THEN
            IF(NCYMAX(II).LE.0) THEN
              WRITE(NFECRA,3010) CHAINE(1:8),'NCYMAX',II,NCYMAX(II)
              IOK = IOK + 1
            ENDIF
            IF(NITMGF(II).LE.0) THEN
              WRITE(NFECRA,3010) CHAINE(1:8),'NITMGF',II,NITMGF(II)
              IOK = IOK + 1
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      IMGROK = 0
      DO IPP = 2, NVPPMX
        II = ITRSVR(IPP)
        IF(II.GE.1) THEN
          IF(IMGR(II).EQ.1) THEN
            IMGROK = 1
          ENDIF
        ENDIF
      ENDDO
C
C --- Options specifiques de niveau plus eleve
C
      IF(IMGROK.EQ.1.AND.NCEGRM.LE.0)THEN
        WRITE(NFECRA,3100)'NCEGRM',NCEGRM
        IOK = IOK + 1
      ENDIF
      IF(IMGROK.EQ.1.AND.NGRMAX.LE.0)THEN
        WRITE(NFECRA,3100)'NGRMAX',NGRMAX
        IOK = IOK + 1
      ENDIF
C
C=======================================================================
C 3. TABLEAUX DE cstphy.h : formats 4000
C=======================================================================
C
      DO IPHAS = 1, NPHAS
C
C --- Constantes physiques de chaque phase
C
        IF(RO0(IPHAS)   .LT.0D0) THEN
          WRITE(NFECRA,2511)IPHAS,'RO0   ', RO0(IPHAS)
          IOK = IOK + 1
        ENDIF
        IF(VISCL0(IPHAS).LT.0D0) THEN
          WRITE(NFECRA,2511)IPHAS,'VISCL0', VISCL0(IPHAS)
          IOK = IOK + 1
        ENDIF
C
C --- Turbulence
C
C    On a besoin de UREF si on initialise la turbulence (debut de calcul
C      en turbulence ou suite laminaire->turbulent)
C    Ici on met juste un avertissement, car sans UREF l'utilisateur peut
C      initialiser la turbulence a la main. Un test complementaire sera fait
C      dans inivar.
        IF(ITYTUR(IPHAS).EQ.2.OR.ITYTUR(IPHAS).EQ.3
     &       .OR.ITURB(IPHAS).EQ.50.OR.ITURB(IPHAS).EQ.60) THEN
          IF(UREF(IPHAS)  .LT.0.D0) THEN
            WRITE(NFECRA,4100)IPHAS, UREF(IPHAS)
          ENDIF
        ENDIF
C
        IF(ITURB(IPHAS).EQ.10) THEN
          IF(XLOMLG(IPHAS).LE.0.D0) THEN
            WRITE(NFECRA,2511)IPHAS,'XLOMLG', XLOMLG(IPHAS)
            IOK = IOK + 1
          ENDIF
        ENDIF
C
C     LES
        IF(ITYTUR(IPHAS).EQ.4) THEN
          IF(XLESFL(IPHAS).LT.0.D0) THEN
            WRITE(NFECRA,2511) IPHAS,'XLESFL', XLESFL(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF(ALES  (IPHAS).LT.0.D0) THEN
            WRITE(NFECRA,2511) IPHAS,'ALES  ', ALES(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF(BLES  (IPHAS).LT.0.D0) THEN
            WRITE(NFECRA,2511) IPHAS,'BLES  ', BLES(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF(CSMAGO(IPHAS).LT.0.D0) THEN
            WRITE(NFECRA,2511)IPHAS, 'CSMAGO', CSMAGO(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF(IDRIES(IPHAS).EQ.1.AND.CDRIES(IPHAS).LT.0) THEN
            WRITE(NFECRA,2511) IPHAS,'CDRIES', CDRIES(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF(ITURB(IPHAS).EQ.41) THEN
            IF(XLESFD(IPHAS).LT.0.D0) THEN
              WRITE(NFECRA,2511) IPHAS,'XLESFD', XLESFD(IPHAS)
              IOK = IOK + 1
            ENDIF
            IF(SMAGMX(IPHAS).LT.0.D0) THEN
              WRITE(NFECRA,2511) IPHAS,'SMAGMX', SMAGMX(IPHAS)
              IOK = IOK + 1
            ENDIF
          ENDIF
        ENDIF
C
      ENDDO
C
C --- Scalaires
C
      IF(NSCAL.GT.0) THEN
C
C     Scalaire passif, temperature, enthalpie, energie
        DO II = 1, NSCAL
          IF (ISCSTH(II).LT.-1.OR.ISCSTH(II).GT.3) THEN
            CHAINE=NOMVAR(IPPRTP(ISCA(II)))
            WRITE(NFECRA,4300)CHAINE(1:8),II,ISCSTH(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C     Phase porteuse
        DO II = 1, NSCAL
          IF (IPHSCA(II).GT.NPHAS.OR.IPHSCA(II).LT.1) THEN
            CHAINE=NOMVAR(IPPRTP(ISCA(II)))
            WRITE(NFECRA,4310)CHAINE(1:8),II,NPHAS,IPHSCA(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C     Scalaire associe dans le cas des variances
        DO II = 1, NSCAL
          IF (ISCAVR(II).GT.NSCAL.OR.ISCAVR(II).LT.0) THEN
            CHAINE=NOMVAR(IPPRTP(ISCA(II)))
            WRITE(NFECRA,4320)CHAINE(1:8),II,NSCAL,ISCAVR(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C     On verifie que le scalaire associe a une fluctuation n'est
C     pas une fluctuation
        DO II = 1, NSCAL
          IF (ISCAVR(II).GT.0) THEN
            IF (ISCAVR(ISCAVR(II)).GT.0) THEN
              CHAINE=NOMVAR(IPPRTP(ISCA(II)))
              CHAIN2=NOMVAR(IPPRTP(ISCA(ISCAVR(II))))
              WRITE(NFECRA,4321)CHAINE(1:8),CHAIN2(1:8),II,ISCAVR(II),
     &             ISCAVR(II),ISCAVR(ISCAVR(II))
              IOK = IOK + 1
            ENDIF
          ENDIF
        ENDDO
C
C     Mode de clipping dans le cas des variances
C       pour les scalaires non variance, on n'utilise pas ICLVFL.
C       On demande donc a l'utilisateur de ne pas y toucher
C       (ca permet d'etre sur qu'il sait ce qu'il fait)
        DO II = 1, NSCAL
          IF (ISCAVR(II).LE.NSCAL.AND.ISCAVR(II).GT.0) THEN
            IF(ICLVFL(II).NE.0.AND.
     &         ICLVFL(II).NE.1.AND.ICLVFL(II).NE.2) THEN
              CHAINE=NOMVAR(IPPRTP(ISCA(II)))
              WRITE(NFECRA,4330)CHAINE(1:8),II,ICLVFL(II)
              IOK = IOK + 1
            ENDIF
          ELSEIF (ISCAVR(II).EQ.0) THEN
            IF(ICLVFL(II).NE.-1) THEN
              CHAINE=NOMVAR(IPPRTP(ISCA(II)))
              WRITE(NFECRA,4331)CHAINE(1:8),II,ICLVFL(II)
              IOK = IOK + 1
            ENDIF
          ENDIF
        ENDDO
C
C     Valeur de la diffusivite positive (si cste)
C        si pas cste, on verifiera apres usphyv
        DO II = 1, NSCAL
          IF (IVISLS(II).LE.0.AND.VISLS0(II).LT.0D0) THEN
            CHAINE=NOMVAR(IPPRTP(ISCA(II)))
            WRITE(NFECRA,4340)CHAINE(1:8),II,II,IVISLS(II),VISLS0(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C     Valeur du sigma positif
        DO II = 1, NSCAL
          IF (SIGMAS(II).LE.0D0) THEN
            CHAINE=NOMVAR(IPPRTP(ISCA(II)))
            WRITE(NFECRA,4350)CHAINE(1:8),II,SIGMAS(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C     Si on n'utilise pas la borne inf de clipping, on demande
C      a l'utilisateur de ne pas y toucher (ca permet d'etre sur
C      qu'il sait ce qu'il fait)
        DO II = 1, NSCAL
          IF(ISCAVR(II).GT.0.AND.ISCAVR(II).LE.NSCAL.AND.
     &       ICLVFL(II).NE.2.AND.SCAMIN(II).NE.-GRAND) THEN
            CHAINE=NOMVAR(IPPRTP(ISCA(II)))
            WRITE(NFECRA,4360)CHAINE(1:8),II,SCAMIN(II),II,ICLVFL(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C     Si on n'utilise pas la borne sup de clipping, on demande
C      a l'utilisateur de ne pas y toucher (ca permet d'etre sur
C      qu'il sait ce qu'il fait)
        DO II = 1, NSCAL
          IF(ISCAVR(II).GT.0.AND.ISCAVR(II).LE.NSCAL.AND.
     &       ICLVFL(II).NE.2.AND.SCAMAX(II).NE.GRAND) THEN
            CHAINE=NOMVAR(IPPRTP(ISCA(II)))
            WRITE(NFECRA,4361)CHAINE(1:8),II,SCAMAX(II),II,ICLVFL(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C     Valeur de la borne sup de clipping si on l'utilise
        DO II = 1, NSCAL
          IF(ISCAVR(II).GT.0.AND.ISCAVR(II).LE.NSCAL.AND.
     &       ICLVFL(II).EQ.2.AND.SCAMAX(II).LE.0.D0) THEN
            CHAINE=NOMVAR(IPPRTP(ISCA(II)))
            WRITE(NFECRA,4370)CHAINE(1:8),II,SCAMAX(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C     Warning sur Rvarfl < 0
        DO II = 1, NSCAL
          IF(ISCAVR(II).GT.0.AND.ISCAVR(II).LE.NSCAL.AND.
     &                           RVARFL(II).LE.0.D0) THEN
            CHAINE=NOMVAR(IPPRTP(ISCA(II)))
            WRITE(NFECRA,4380)CHAINE(1:8),II,RVARFL(II)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C
C  Rien a verifier sur SCAMIN SCAMAX
C
C     Si CP0 est utilise (resolution d'un scalaire en temperature
C       et CP constant), il doit etre positif
        DO IPHAS = 1, NPHAS
          IF(ICP(IPHAS).EQ.0) THEN
            IF (CP0(IPHAS).LT.0.D0) THEN
              IISCT = 0
              DO IIS = 1, NSCAL
                IF (IPHSCA(IIS).EQ.IPHAS.AND.
     &              ABS(ISCSTH(IIS)).EQ.1) THEN
                  IISCT = 1
                ENDIF
              ENDDO
              IF (IISCT.EQ.1) THEN
                WRITE(NFECRA,2511)IPHAS,'CP0   ',CP0(IPHAS)
                IOK = IOK + 1
              ENDIF
            ENDIF
          ENDIF
        ENDDO
C
      ENDIF
C
C=======================================================================
C 4. TABLEAUX DE period.h : formats 5000
C=======================================================================
C
C --- periodicite de rotation incompatible avec couplage
C       renforce vitesse pression et l'ALE
      IF(IPEROT.GT.0.AND. (IPUCOU.NE.0.OR.IALE.NE.0)) THEN
        WRITE(NFECRA,5003)IPERIO,IPUCOU,IALE
        IOK = IOK + 1
      ENDIF
C
C --- periodicite incompatible avec le mode de calcul
C       direct de la distance a la paroi
      IF(IPERIO.EQ.1.AND.INEEDY.EQ.1.AND.ABS(ICDPAR).EQ.2) THEN
        WRITE(NFECRA,5005)IPERIO, ICDPAR
        IOK = IOK + 1
      ENDIF
C
C
cC --- periodicite incompatible avec le rayonnement DOM
c      IF(IPERIO.EQ.1.AND.IIRAYO.GT.0) THEN
c        DO IPHAS = 1, NPHAS
c          IF(IRAYON(IPHAS).EQ.1) THEN
c            WRITE(NFECRA,5007) IPERIO , IPHAS, IRAYON(IPHAS)
c            IOK = IOK + 1
c          ENDIF
c        ENDDO
c      ENDIF
C
C --- periodicite de rotation incompatible avec le rayonnement DOM
      IF(IPEROT.GT.0.AND.IIRAYO.GT.0) THEN
        DO IPHAS = 1, NPHAS
          IF(IRAYON(IPHAS).EQ.1) THEN
            WRITE(NFECRA,5008) IPERIO , IPHAS, IRAYON(IPHAS)
            IOK = IOK + 1
          ENDIF
        ENDDO
      ENDIF
C
C --- periodicite de rotation douteuse avec rij
C      (et donc a fortiori avec ordre 2 sur Rij)
      IF(IPEROT.GT.0) THEN
        DO IPHAS = 1, NPHAS
          IF(ITYTUR(IPHAS).EQ.3) THEN
            WRITE(NFECRA,5009)IPERIO,ITURB(IPHAS),IPHAS
C            IOK = IOK + 1
          ENDIF
        ENDDO
      ENDIF
C
C --- periodicite de rotation douteuse avec ordre 2 vitesse
      IF(IPEROT.GT.0) THEN
        DO IPHAS = 1, NPHAS
          IUIPH = IU(IPHAS)
          IVIPH = IV(IPHAS)
          IWIPH = IW(IPHAS)
          IF( (ABS(THETAV(IUIPH)-0.5D0).LT.1.D-3).OR.
     &        (ABS(THETAV(IVIPH)-0.5D0).LT.1.D-3).OR.
     &        (ABS(THETAV(IWIPH)-0.5D0).LT.1.D-3)) THEN
            WRITE(NFECRA,5010)IPERIO,IPHAS,
     &        THETAV(IUIPH),THETAV(IVIPH),THETAV(IWIPH)
C            IOK = IOK + 1
          ENDIF
        ENDDO
      ENDIF
C
C
C=======================================================================
C 5. TABLEAUX DE parall.h : formats 6000 (limitations)
C=======================================================================
C
C --- parallelisme incompatible avec lagrangien
      IF(IRANGP.GE.0.AND.IILAGR.NE.0) THEN
        WRITE(NFECRA,6002)IRANGP,IILAGR
        IOK = IOK + 1
      ENDIF
C
C --- parallelisme incompatible avec le mode de calcul
C       direct de la distance a la paroi
      IF(IRANGP.GE.0.AND.INEEDY.EQ.1.AND.ABS(ICDPAR).EQ.2) THEN
        WRITE(NFECRA,6005)IPERIO, ICDPAR
        IOK = IOK + 1
      ENDIF
C
C=======================================================================
C 6. METHODE ALE (albase.h, alstru.h) : formats 7000
C=======================================================================
C
      IF (IALE.NE.0 .AND. IALE.NE.1) THEN
        WRITE(NFECRA,7000)IALE
        IOK = IOK + 1
      ENDIF
C
      IF (IALE.EQ.1) THEN
C
        IF (NALINF.LT.0) THEN
          WRITE(NFECRA,7010)NALINF
          IOK = IOK + 1
        ENDIF
C
        IF (ALPNMK.LT.0.D0 .OR. ALPNMK.GT.1.D0  .OR.
     &       GAMNMK.LT.0.D0 .OR. GAMNMK.GT.1.D0  .OR.
     &       BETNMK.LT.0.D0 .OR. BETNMK.GT.0.5D0 ) THEN
          WRITE(NFECRA,7020)ALPNMK,BETNMK,GAMNMK
          IOK = IOK + 1
        ENDIF
C
        IF (NALIMX.LE.0) THEN
          WRITE(NFECRA,7030)NALIMX
          IOK = IOK + 1
        ENDIF
        IF (EPALIM.LE.0.D0) THEN
          WRITE(NFECRA,7040)EPALIM
          IOK = IOK + 1
        ENDIF
C
        IF (ITALIN.NE.-999 .AND. ITALIN.NE.0 .AND. ITALIN.NE.1) THEN
          WRITE(NFECRA,7050)ITALIN
          IOK = IOK + 1
        ENDIF

      ENDIF
C
C=======================================================================
C 7. COMPRESSIBLE : formats 8000
C=======================================================================
C
      IF(IPPMOD(ICOMPF).GE.0) THEN
        DO IPHAS = 1, NPHAS
          IF(T0(IPHAS).LE.0.D0.OR.P0(IPHAS).LE.0.D0) THEN
            WRITE(NFECRA,8000)IPHAS,T0(IPHAS),P0(IPHAS)
            IOK = IOK + 1
          ENDIF
        ENDDO
      ENDIF
C
C=======================================================================
C 8. FORMATS VERIFICATION
C=======================================================================
C
 1100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    NDIM   DOIT ETRE UN ENTIER EGAL A 3                     ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute                            ',/,
     &'@                                                            ',/,
     &'@  La dimension de l espace est 3, meme pour les calculs     ',/,
     &'@    physiquement bidimensionnels.                           ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@  Verifier le fichier maillage.                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ',A33,                          ' DOIT ETRE UN ENTIER   ',/,
     &'@    SUPERIEUR OU EGAL A -1                                  ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1210 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ',A33,                          ' DOIT ETRE UN ENTIER   ',/,
     &'@    STRICTEMENT POSITIF OU EGAL A -1                        ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1220 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    ICHRVR(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0 OU 1    ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  ICHRVR indique si la variable doit etre incluse dans les  ',/,
     &'@    fichiers de post-traitement                             ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1230 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    NCAPT  DOIT ETRE UN ENTIER INFERIEUR OU EGAL A ',I10     ,/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  NCAPT  est le nombre de sondes utilisees pour les         ',/,
     &'@    historiques.                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1240 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    IHISVR(',I10   ,',1) DOIT ETRE UN ENTIER EGAL A -1      ',/,
     &'@      POSITIF OU NUL ET INFERIEUR OU EGAL A NCAPT =',I10     ,/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  IHISVR(I,1) indique les sondes a utiliser pour la variable',/,
     &'@    I (-1 signifiant que toutes sont utilisees)             ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1250 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    IHISVR(',I10   ,',',I10   ,') DOIT ETRE UN ENTIER       ',/,
     &'@      STRICTEMENT POSITIF ET                                ',/,
     &'@      INFERIEUR OU EGAL A NCAPT = ',I10                      ,/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute                            ',/,
     &'@                                                            ',/,
     &'@  IHISVR(I,j+1) indique le numero de la jieme sonde a       ',/,
     &'@    utiliser pour la variable a post-traiter numero I       ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1260 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    ILISVR(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0 OU 1    ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute                            ',/,
     &'@                                                            ',/,
     &'@  ILISVR(I) indique si la variable I sera suivie lors des   ',/,
     &'@    impressions dans le listing                             ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1270 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    ITRSVR(',I10   ,') DOIT ETRE UN ENTIER COMPRIS ENTRE    ',/,
     &'@      0 ET NVAR=',I10                                        ,/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1300 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LE PARAMETRE ',A6,' DOIT ETRE UN ENTIER MULTIPLE DE     ',/,
     &'@      DE CERTAINS DES ENTIERS SUIVANTS :                    ',/,
     &'@      ',A6,' = ',I10                                         ,/,
     &'@      ',A6,' = ',I10                                         ,/,
     &'@      ',A6,' = ',I10                                         ,/,
     &'@      ',A6,' = ',I10                                         ,/,
     &'@    IL VAUT ICI ',A6,' = ', I10                              ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute                            ',/,
     &'@                                                            ',/,
     &'@  Ce parametre precise les variables supplementaires        ',/,
     &'@    a post-traiter.                                         ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN ENTIER                              ',/,
     &'@      STRICTEMENT POSITIF ET INFERIEUR OU EGAL A ',I10       ,/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2005 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@    ON DEMANDE UNE EXTRAPOLATION TEMPORELLE DE RHO AVEC     ',/,
     &'@      IROEXT(IPHAS) = ',I10                                  ,/,
     &'@    CECI EST INCOMPATIBLE AVEC RHO CONSTANT                 ',/,
     &'@      IROVAR(IPHAS) = ',I10                                  ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     PARAMETRES DU SCHEMA NUMERIQUE POUR LA VARIABLE ',A8    ,/,
     &'@                                                            ',/,
     &'@ Parametre               ISTAT     ICONV                    ',/,
     &'@ Valeurs acceptees     0 ou  1   0 ou  1                    ',/,
     &'@ Valeurs entrees ici',I10     ,I10                           ,/,
     &'@                                                            ',/,
     &'@ Parametre                         IDIFF     IDIFFT         ',/,
     &'@ Valeurs acceptees               0 ou  1   0 ou  1          ',/,
     &'@ Valeurs entrees ici',10X     ,I10      ,I10                 ,/,
     &'@                                                            ',/,
     &'@ Parametre               THETAV   BLENCV    ISCHCV    ISSTPC',/,
     &'@ Valeurs acceptees     [0.; 1.] [0.; 1.]   0 ou  1   0 ou  1',/,
     &'@ Valeurs entrees ici      ',E14.5     ,E14.5,I10  ,I10       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2111 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@   INCOMPATIBILITE POUR LE SCHEMA EN TEMPS                  ',/,
     &'@                                                            ',/,
     &'@   Schema en temps pour la vitesse phase ',I10               ,/,
     &'@      THETA n''a pas la meme valeur pour les 3 composantes  ',/,
     &'@                                                            ',/,
     &'@ Parametre THETAV              U          V          W      ',/,
     &'@ Valeurs entrees ici ',E10.2,E10.2,E10.2                     ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2112 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@   DONNEES NON ADMISSIBLES POUR LE SCHEMA EN TEMPS          ',/,
     &'@                                                            ',/,
     &'@  LE PARAMETRE THETAV POUR LA PRESSION DOIT VALOIR 1        ',/,
     &'@                                                            ',/,
     &'@  Pour la phase ', I10,' il vaut ici ',E14.5                 ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2113 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@   DONNEES NON ADMISSIBLES POUR LE SCHEMA EN TEMPS          ',/,
     &'@                                                            ',/,
     &'@   PHASE ',I10   ,' EN L.E.S.                               ',/,
     &'@   LE PARAMETRE THETFL DU SCHEMA EN TEMPS POUR LE FLUX DE   ',/,
     &'@     MASSE EST DIFFERENT DE CELUI DE LA PHASE 1 ',E14.5      ,/,
     &'@     THETFL A ETE IMPOSE ICI A ',E14.5                       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute                             ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2114 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@   DONNEES NON ADMISSIBLES POUR LE SCHEMA EN TEMPS          ',/,
     &'@                                                            ',/,
     &'@   PHASE ',I10   ,'                                         ',/,
     &'@   LE PARAMETRE ISTMPF DU SCHEMA EN TEMPS POUR LE FLUX DE   ',/,
     &'@     MASSE EST DIFFERENT DE CELUI DE LA PHASE 1 ',I10        ,/,
     &'@     ISTMPF A ETE IMPOSE ICI A ',I10                         ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute                             ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2121 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,
     &'@    *********                                               ',/,
     &'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,
     &'@                                                            ',/,
     &'@   PHASE ',I10   ,' EN L.E.S.                               ',/,
     &'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE THETAV DU SCHEMA ',/,
     &'@     EN TEMPS DE LA VARIABLE ',A8  ,' EST 0.5               ',/,
     &'@     THETAV A ETE IMPOSE ICI A ',E14.5                       ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute                                    ',/,
     &'@                                                            ',/,
     &'@  Il est conseille de verifier les parametres donnes via    ',/,
     &'@  l''interface ou usini1.                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2122 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,
     &'@    *********                                               ',/,
     &'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,
     &'@                                                            ',/,
     &'@   PHASE ',I10   ,' EN L.E.S.                               ',/,
     &'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE BLENCV DU SCHEMA ',/,
     &'@     CONVECTIF DE LA VARIABLE ',A8  ,' EST 1.0              ',/,
     &'@     BLENCV A ETE IMPOSE ICI A ',E14.5                       ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute                                    ',/,
     &'@                                                            ',/,
     &'@  Il est conseille de verifier les parametres donnes via    ',/,
     &'@  l''interface ou usini1.                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2123 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,
     &'@    *********                                               ',/,
     &'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,
     &'@                                                            ',/,
     &'@   PHASE ',I10   ,' EN L.E.S.                               ',/,
     &'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE ISSTPC DU SCHEMA ',/,
     &'@     CONVECTIF DE LA VARIABLE ',A8  ,' EST 1                ',/,
     &'@     ISSTPC A ETE IMPOSE ICI A ',I10                         ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute                                    ',/,
     &'@                                                            ',/,
     &'@  Il est conseille de verifier les parametres donnes via    ',/,
     &'@  l''interface ou usini1.                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2124 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,
     &'@    *********                                               ',/,
     &'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,
     &'@                                                            ',/,
     &'@   PHASE ',I10   ,' EN L.E.S.                               ',/,
     &'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE ISSTPC DU SCHEMA ',/,
     &'@     CONVECTIF DE LA VARIABLE ',A8  ,' EST 0                ',/,
     &'@     ISSTPC A ETE IMPOSE ICI A ',I10                         ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute                                    ',/,
     &'@                                                            ',/,
     &'@  Il est conseille de verifier les parametres donnes via    ',/,
     &'@  l''interface ou usini1.                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2125 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,
     &'@    *********                                               ',/,
     &'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,
     &'@                                                            ',/,
     &'@   PHASE ',I10   ,' ORDRE 2 EN TEMPS OU LES                 ',/,
     &'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE NSWRSM POUR      ',/,
     &'@     LA VARIABLE ',A8  ,' EST  ',I10                        ,/,
     &'@     NSWRSM A ETE IMPOSE ICI A ',I10                         ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute                                    ',/,
     &'@                                                            ',/,
     &'@  Il est conseille de verifier les parametres donnes via    ',/,
     &'@  l''interface ou usini1.                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2127 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,
     &'@                                                            ',/,
     &'@   PHASE ',I10   ,' EN L.E.S.                               ',/,
     &'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE BLENCV DU SCHEMA ',/,
     &'@     CONVECTIF DE LA VARIABLE ',A8  ,' EST 1.0              ',/,
     &'@     BLENCV A ETE IMPOSE ICI A ',E14.5                       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute                             ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2131 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,
     &'@    *********                                               ',/,
     &'@   CHOIX DU SCHEMA EN TEMPS                                 ',/,
     &'@                                                            ',/,
     &'@     PHASE ',I10                                             ,/,
     &'@     LE SCHEMA EN TEMPS POUR LA VITESSE EST D ORDRE 1       ',/,
     &'@       (THETAV = ',E10.2 ,')                                ',/,
     &'@     CERTAINS TERMES SONT CEPENDANT PRIS A L''ORDRE 2 AVEC  ',/,
     &'@       LES CHOIX SUIVANTS :                                 ',/,
     &'@                                                            ',/,
     &'@ Parametres       ISTMPF ISNO2T ISTO2T IROEXT IVIEXT ICPEXT ',/,
     &'@ Valeurs entrees ',6I7                                       ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute.                                   ',/,
     &'@                                                            ',/,
     &'@  Il est conseille de verifier les parametres donnes via    ',/,
     &'@  l''interface ou usini1.                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2132 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,
     &'@    *********                                               ',/,
     &'@   CHOIX DU SCHEMA EN TEMPS                                 ',/,
     &'@                                                            ',/,
     &'@     PHASE ',I10                                             ,/,
     &'@     LE SCHEMA EN TEMPS POUR LA VITESSE EST D ORDRE 2       ',/,
     &'@       (THETAV = ',E10.2 ,')                                ',/,
     &'@     CERTAINS TERMES SONT CEPENDANT PRIS A L''ORDRE 1 AVEC  ',/,
     &'@       LES CHOIX SUIVANTS :                                 ',/,
     &'@                                                            ',/,
     &'@ Parametres       ISTMPF ISNO2T ISTO2T IROEXT IVIEXT ICPEXT ',/,
     &'@ Valeurs entrees ',6I7                                       ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute.                                   ',/,
     &'@                                                            ',/,
     &'@  Il est conseille de verifier les parametres donnes via    ',/,
     &'@  l''interface ou usini1.                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2133 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,
     &'@    *********                                               ',/,
     &'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,
     &'@                                                            ',/,
     &'@   SCALAIRE ',I10,' ISSO2T = ',I10                           ,/,
     &'@     EST DIFFERENT DE ISNO2T DE LA PHASE IPHAS = ',I10       ,/,
     &'@     ISNO2T(IPHAS) = ',I10                                   ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute                                    ',/,
     &'@                                                            ',/,
     &'@  Il est conseille de verifier les parametres donnes via    ',/,
     &'@  l''interface ou usini1.                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2134 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,
     &'@    *********                                               ',/,
     &'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,
     &'@                                                            ',/,
     &'@   SCALAIRE ',I10,' IVSEXT = ',I10                           ,/,
     &'@     EST DIFFERENT DE IVIEXT DE LA PHASE IPHAS = ',I10       ,/,
     &'@     IVIEXT(IPHAS) = ',I10                                   ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute                                    ',/,
     &'@                                                            ',/,
     &'@  Il est conseille de verifier les parametres donnes via    ',/,
     &'@  l''interface ou usini1.                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2135 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,
     &'@    *********                                               ',/,
     &'@   CHOIX INCOMPATIBLE POUR LE SCHEMA EN TEMPS               ',/,
     &'@                                                            ',/,
     &'@   Phase IPHAS = ',I10                                       ,/,
     &'@     La  chaleur massique est extrapolee en temps avec      ',/,
     &'@       ICPEXT(IPHAS) = ',I10                                 ,/,
     &'@     Pour cela, elle doit etre variable, or                 ',/,
     &'@       ICP(IPHAS)    = ',I10                                 ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute                             ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1 ',/,
     &'@    - desactiver le choix d''extrapolation de Cp en temps   ',/,
     &'@      ou                                                    ',/,
     &'@    - imposer Cp variable                                   ',/,
     &'@         (et le renseigner alors dans usphyv)               ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2136 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,
     &'@    *********                                               ',/,
     &'@   CHOIX INCOMPATIBLE POUR LE SCHEMA EN TEMPS               ',/,
     &'@                                                            ',/,
     &'@   Scalaire ISCAL = ',I10                                    ,/,
     &'@     La  diffusivite      est extrapolee en temps avec      ',/,
     &'@       IVSEXT(ISCAL) = ',I10                                 ,/,
     &'@     Pour cela, elle doit etre variable, or                 ',/,
     &'@       IVISLS(ISCAL) = ',I10                                 ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute                             ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1 ',/,
     &'@    - desactiver le choix d''extrapolation en temps         ',/,
     &'@                                     de la diffusivite      ',/,
     &'@      ou                                                    ',/,
     &'@    - imposer la diffusivite variable                       ',/,
     &'@         (et la renseigner alors dans usphyv)               ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2137 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,
     &'@    *********                                               ',/,
     &'@   CHOIX INCOMPATIBLE POUR LES ESTIMATEURS D''ERREUR        ',/,
     &'@                                                            ',/,
     &'@  On a active un ou plusieurs estimateurs d''erreur pour    ',/,
     &'@    Navier-Stokes dans un calcul a champ de vitesse fige.   ',/,
     &'@    Le ou les estimateurs ne seront pas calcules.           ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute                             ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou        ',/,
     &'@    usini1 :                                                ',/,
     &'@      desactiver les estimateurs d erreur ou                ',/,
     &'@               le calcul a champ de vitesse fige (ICCVFG)   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2140 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Pour la phase ',I10   ,' on souhaite utiliser un schema en',/,
     &'@    temps d''ordre 2 :                                      ',/,
     &'@      U,V,W : THETA = ',3E12.4                               ,/,
     &'@      Termes sources Navier-Stokes: ISNO2T = ',I10           ,/,
     &'@                                    THETSN = ',E12.4         ,/,
     &'@      Masse volumique             : IROEXT = ',I10           ,/,
     &'@                                    THETRO = ',E12.4         ,/,
     &'@      Viscosite                   : IVIEXT = ',I10           ,/,
     &'@                                    THETVI = ',E12.4         ,/,
     &'@  La version actuelle ne le permet pas lorsque l''une des   ',/,
     &'@    options suivantes a ete activee (c''est le cas ici) :   ',/,
     &'@    - utilisation d''un estimateur d''erreur (IESCAL)       ',/,
     &'@    - couplage instationnaire (IPUCOU)                      ',/,
     &'@    - prise en compte specifique de la pression             ',/,
     &'@      hydrostatique (IPHYDR et ICALHY)                      ',/,
     &'@    - pas de temps variable en temps ou en espace  ou       ',/,
     &'@      algorithme stationnaire (IDTVAR)                      ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2141 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,
     &'@                                                            ',/,
     &'@  On souhaite utiliser un couplage                          ',/,
     &'@    vitesse-pression par point fixe (NTERUP = ',I10          ,/,
     &'@  La version actuelle ne le permet pas lorsque l''une des   ',/,
     &'@    options suivantes a ete activee (c''est le cas ici) :   ',/,
     &'@    - utilisation d''un estimateur d''erreur (IESCAL)       ',/,
     &'@    - couplage instationnaire (IPUCOU)                      ',/,
     &'@    - prise en compte specifique de la pression             ',/,
     &'@      hydrostatique (IPHYDR et ICALHY)                      ',/,
     &'@    - algorithme stationnaire (IDTVAR=-1)                   ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2142 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Pour la phase ',I10   ,' avec le modele de turbulence     ',/,
     &'@    k-epsilon (ITURB = ',I10 ,') couple (IKECOU = ',I10,') :',/,
     &'@    la version courante ne permet pas de traiter les        ',/,
     &'@    equations du modele k-epsilon a l''ordre 2 en temps avec',/,
     &'@    couplage.                                               ',/,
     &'@    Une ou plusieurs valeurs parmi les suivantes ne sont    ',/,
     &'@    donc pas permises :                                     ',/,
     &'@                                                            ',/,
     &'@       THETST    ISTO2T     THETA K   THETA EPS             ',/,
     &'@ ',     E12.4,      I10,      E12.4,      E12.4              ,/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2143 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Pour la phase ',I10   ,' avec le modele de turbulence     ',/,
     &'@    v2f (ITURB = ',I10   ,') couple (IKECOU = ',I10,    ') :',/,
     &'@    la version courante ne permet pas de traiter les        ',/,
     &'@    equations du modele v2f a l''ordre 2 en temps avec      ',/,
     &'@    couplage.                                               ',/,
     &'@    Une ou plusieurs valeurs parmi les suivantes ne sont    ',/,
     &'@    donc pas permises :                                     ',/,
     &'@                                                            ',/,
     &'@       THETST    ISTO2T     THETA K   THETA EPS             ',/,
     &'@ ',     E12.4,      I10,      E12.4,      E12.4              ,/,
     &'@     THETA PHI    THETA FB                                  ',/,
     &'@ ',      E12.4,      E12.4                                   ,/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2144 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Pour la phase ',I10   ,' avec le modele de turbulence     ',/,
     &'@    k-omega (ITURB = ',I10   ,') couple (IKECOU = ',I10,') :',/,
     &'@    la version courante ne permet pas de traiter les        ',/,
     &'@    equations du modele k-omega a l''ordre 2 en temps avec  ',/,
     &'@    couplage.                                               ',/,
     &'@    Une ou plusieurs valeurs parmi les suivantes ne sont    ',/,
     &'@    donc pas permises :                                     ',/,
     &'@                                                            ',/,
     &'@       THETST    ISTO2T     THETA K   THETA OMEGA           ',/,
     &'@ ',     E12.4,      I10,      E12.4,      E12.4              ,/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2145 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  La version actuelle ne permet pas de modifier le schema en',/,
     &'@    temps lorsqu''une physique particuliere est activee     ',/,
     &'@    (combustion, charbon, electrique).                      ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface, usini1   ',/,
     &'@    et usppmo.                                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2146 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Phase ',I10   ,' les termes sources provenant du module   ',/,
     &'@    Lagrangien ne sont pas traites a l''ordre 2 en temps    ',/,
     &'@    dans la version courante malgre le choix utilisateur    ',/,
     &'@    suivant :                                               ',/,
     &'@                                                            ',/,
     &'@     Pour Navier-Stokes    Pour la turbulence               ',/,
     &'@       THETSN    ISNO2T      THETST    ISTO2T               ',/,
     &'@ ',     E12.4,      I10,      E12.4,      I10                ,/,
     &'@                                                            ',/,
     &'@  (Les autres termes sources pourraient etre traites a      ',/,
     &'@   l''ordre 2)                                              ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface, usini1   ',/,
     &'@  et uslag1.                                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2147 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Les termes sources provenant du module ',A11               ,/,
     &'@    ne sont pas traites a l''ordre 2 en temps               ',/,
     &'@    dans la version courante malgre le choix utilisateur    ',/,
     &'@    suivant :                                               ',/,
     &'@                                                            ',/,
     &'@       Pour le scalaire ',I10                                ,/,
     &'@       THETSS    ISSO2T                                     ',/,
     &'@ ',     E12.4,      I10                                      ,/,
     &'@                                                            ',/,
     &'@  (Les autres termes sources pourraient etre traites a      ',/,
     &'@   l''ordre 2)                                              ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface, usini1   ',/,
     &'@  et ',A6                                                    ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2148 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     ALGORITHME STATIONNAIRE                                ',/,
     &'@     COEFFICIENT DE RELAXATION POUR LA VARIABLE ',A8         ,/,
     &'@                                                            ',/,
     &'@ RELAXV doit etre un reel compris entre 0 et 1              ',/,
     &'@ Il vaut ici ',E14.5                                         ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2149  FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@   INCOMPATIBILITE ALGORITHME STATIONNAIRE                  ',/,
     &'@                                                            ',/,
     &'@   Coefficient de relaxation vitesse phase ',I10             ,/,
     &'@      RELAXV n''a pas la meme valeur pour les 3 composantes ',/,
     &'@                                                            ',/,
     &'@ Parametre RELAXV              U          V          W      ',/,
     &'@ Valeurs entrees ici ',E10.2,E10.2,E10.2                     ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2150  FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@   L''ALGORITHME STATIONNAIRE N''EST PAS COMPATIBLE AVEC    ',/,
     &'@    LE MODULE LAGRANGIEN QUI EST UNE APPROCHE INSTATIONNAIRE',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  L''indicateur IILAGR a ete positionne a ',I10              ,/,
     &'@    dans uslag1 (module lagrangien active pour IILAGR>0).   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2151  FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@   L''ALGORITHME STATIONNAIRE N''EST PAS COMPATIBLE AVEC    ',/,
     &'@    LA L.E.S. QUI EST UNE MODELISATION INSTATIONNAIRE DE    ',/,
     &'@    LA TURBULENCE                                           ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  L''indicateur ITURB a ete positionne a ',I10               ,/,
     &'@    dans usini1.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0 OU 1                ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2201 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0 OU 1                ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2202 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LA PRISE EN COMPTE EXPLICITE DE LA PRESSION             ',/,
     &'@    HYDROSTATIQUE N''EST IMPLANTEE QUE POUR LE MONOPHASIQUE ',/,
     &'@    IPHYDR VAUT ICI ',I10                                    ,/,
     &'@    ET NPHAS        ',I10                                    ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C 2203 FORMAT(
C     &'@                                                            ',/,
C     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
C     &'@                                                            ',/,
C     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
C     &'@    *********                                               ',/,
C     &'@  ON DEMANDE LA PRISE EN COMPTE EXPLICITE DE LA PRESSION    ',/,
C     &'@    HYDROSTATIQUE POUR LES CONDITIONS DE SORTIE EN PRESSION ',/,
C     &'@    AVEC ACCELERATION DE LA PESANTEUR NULLE                 ',/,
C     &'@                                                            ',/,
C     &'@  ICALHY VAUT ICI ',I10                                      ,/,
C     &'@  G      VAUT ICI ',3E14.5                                   ,/,
C     &'@                                                            ',/,
C     &'@  Le calcul ne sera pas execute.                            ',/,
C     &'@                                                            ',/,
C     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
C     &'@                                                            ',/,
C     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
C     &'@                                                            ',/)
 2204 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@   DONNEES NON ADMISSIBLES POUR LE SCHEMA EN TEMPS          ',/,
     &'@                                                            ',/,
     &'@  ON DEMANDE LA PRISE UN SCHEMA EN TEMPS POUR LA VITESSE    ',/,
     &'@    D''ORDRE 2 AVEC UN PAS DE TEMPS NON CONSTANT, UN        ',/,
     &'@    ALGORITHME STATIONNAIRE OU LES MATRICES POIDS           ',/,
     &'@                                                            ',/,
     &'@  THETAV VAUT ICI ',E14.5,' POUR LA VITESSE                 ',/,
     &'@  ALORS QUE IDTVAR VAUT ',I10                                ,/,
     &'@  ET IPUCOU VAUT        ',I10                                ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2205 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0, 1, 2, 3 OU 4       ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2206 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ANOMAX DOIT ETRE UN REEL POSITIF OU NUL ET              ',/,
     &'@                             INFERIEUR OU EGAL A PI/2       ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  On demande la reconstruction des gradients par moindres   ',/,
     &'@    carres sur voisinage etendu reduit (IMRGRA = ',I10  ,').',/,
     &'@    Le critere est base sur l''angle de non orthogonalite   ',/,
     &'@    des faces ANOMAX qui doit etre fourni en radians et     ',/,
     &'@    compris dans les bornes indiquees ci-dessus.            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2207 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    L''UTILISATION DE LA METHODE DE CALCUL DE GRADIENT PAR  ',/,
     &'@      MOINDRES CARRES EST IMPOSSIBLE AVEC                   ',/,
     &'@        NPHAS SUPERIEUR A 1   OU                            ',/,
     &'@        EXTRAG(IPR(1)) DIFFERENT DE 0 ET 1                  ',/,
     &'@                                                            ',/,
     &'@    ON A ICI                                                ',/,
     &'@        IMRGRA         = ',I10                               ,/,
     &'@        NPHAS          = ',I10                               ,/,
     &'@        EXTRAG(IPR(1)) = ',E14.5                             ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2208 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@    EN K-EPS PROD LIN (ITURB=21) ET EN V2F (ITURB=50)       ',/,
     &'@    IKECOU DOIT ETRE EGAL A 0                               ',/,
     &'@    ITURB  VAUT ICI ',I10                                    ,/,
     &'@    IKECOU VAUT ICI ',I10                                    ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2209 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LE MODELE DE PAROI A DEUX ECHELLES (IDEUCH=1 OU 2)      ',/,
     &'@    EST INCOMPATIBLE AVEC UN CALCUL EN LAMINAIRE, EN        ',/,
     &'@    LONGUEUR DE MELANGE OU EN L.E.S.                        ',/,
     &'@    POUR LA PHASE ',I10                                      ,/,
     &'@    ON A ICI ITURB=',I10                                     ,/,
     &'@         ET IDEUCH=',I10                                     ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2210 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    L''ALGORITHME STATIONNAIRE EST INCOMPATIBLE AVEC LE     ',/,
     &'@    COUPLAGE DES TERMES SOURCES EN K-EPS, V2F OU K-OMEGA    ',/,
     &'@                                                            ',/,
     &'@    POUR LA PHASE   ',I10                                    ,/,
     &'@    ON A ICI IKECOU=',I10                                    ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2211 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0, 1 OU 2             ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2300 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    IMLIGR(',I10   ,') DOIT ETRE UN ENTIER                  ',/,
     &'@      INFERIEUR OU EGAL A 1                                 ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  IMLIGR(I) indique le mode de limitation des gradients     ',/,
     &'@    pour la variable I                                      ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2310 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    IRCFLU(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0 OU 1    ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  IRCFLU(I) indique si les flux sont reconstruits           ',/,
     &'@    pour la variable I                                      ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2311 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    IRCFLU(',I10   ,') = ',I10   ,' EST INCOMPATIBLE AVEC   ',/,
     &'@    ISCHCV(',I10   ,') = ',I10                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  IRCFLU(I) = 0 indique que les flux ne sont pas            ',/,
     &'@    reconstruits pour la variable I.                        ',/,
     &'@  ISCHCV(I) = 0 (schema SOLU) demande une reconstruction    ',/,
     &'@    pour les flux convectifs.                               ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2320 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    CLIMGR(',I10   ,') DOIT ETRE UN REEL                    ',/,
     &'@      SUPERIEUR OU EGAL A 1                                 ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  CLIMGR(I) est le coefficient de limitation des gradients  ',/,
     &'@    pour la variable I                                      ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2330 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    EXTRAG(',I10   ,') DOIT ETRE UN REEL EGAL A 0 OU 1      ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2331 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    EXTRAG(',I10   ,') DOIT ETRE NUL                        ',/,
     &'@      (Des valeurs non nulles sont autorisees pour la       ',/,
     &'@       pression uniquement)                                 ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2400 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    IRESOL(',I10   ,') DOIT ETRE UN ENTIER EGAL             ',/,
     &'@                                     A -1 OU A IPOL*1000+ J ',/,
     &'@      AVEC IPOL LE DEGRE DU POLYNOME DE PRECONDITIONNEMENT  ',/,
     &'@        ET J    = 0 POUR GRADIENT CONJUGUE                  ',/,
     &'@                = 1 POUR JACOBI   (IPOL = 0 DANS CE CAS)    ',/,
     &'@                = 2 POUR BI-CGSTAB                          ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  IRESOL(I) indique le solveur lineaire a utiliser          ',/,
     &'@    pour la variable I                                      ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2401 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    IDIRCL(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0 OU 1    ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  IDIRCL(I) indique si le code doit decaler la diagonale de ',/,
     &'@    la matrice de la variable I en l''absence de Dirichlet  ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2410 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    RISQUE D ECHEC A LA RESOLUTION DU SYSTEME LINEAIRE      ',/,
     &'@    IRESOL(',I10   ,') = ',I10                               ,/,
     &'@      ET LA VARIABLE EST CONVECTEE (ICONV = ',I10,')        ',/,
     &'@                                                            ',/,
     &'@  Le calcul sera engage.                                    ',/,
     &'@                                                            ',/,
     &'@  Le solveur iteratif choisi peut ne pas converger sur le   ',/,
     &'@    systeme lineaire resultant du type de probleme considere',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2420 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE    ',I10                                           ,/,
     &'@    RISQUE DE PERTE D''INFORMATION EN CALCUL SUITE          ',/,
     &'@                                                            ',/,
     &'@  Le calcul sera engage.                                    ',/,
     &'@                                                            ',/,
     &'@  Un modele de LES a ete active par ITURB(IPHAS) = ',I10     ,/,
     &'@    mais on a desactive l''ecriture ou la lecture du fichier',/,
     &'@    suite auxiliaire :                                      ',/,
     &'@    ILEAUX = ',I10   ,'    IECAUX = ',I10                    ,/,
     &'@  Bien que ce fichier ne soit pas necessaire a la poursuite ',/,
     &'@    d''un calcul, il contient neanmoins des informations    ',/,
     &'@    qui permettent d''eviter les perturbations numeriques   ',/,
     &'@    au moment des suites.                                   ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2500 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN ENTIER EGAL A -1, 0, 1 OU 2         ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2510 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN REEL POSITIF                        ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2511 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@    ',A6,' DOIT ETRE UN REEL POSITIF                        ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2520 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     DTMIN DOIT ETRE SUPERIEUR OU EGAL A 0. ET              ',/,
     &'@                     INFERIEUR OU EGAL A DTMAX              ',/,
     &'@      ICI DTMIN = ',E14.5     ,' ET DTMAX = ',E14.5          ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Quand le pas de temps n est pas uniforme et constant,     ',/,
     &'@    les reels DTMIN et DTMAX bornent ses variations.        ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2530 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    CDTVAR(',I10   ,') DOIT ETRE UN REEL STRICTEMENT POSITIF',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  CDTVAR(I) est le coefficient multiplicatif applique au pas',/,
     &'@    de temps pour la resolution de la variable I.           ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2540 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@  ON DEMANDE UNE LIMITATION DU PAS DE TEMPS LIEE AUX EFFETS ',/,
     &'@    DE DENSITE (IPTLRO = ',I10   ,') AVEC UNE OPTION DE     ',/,
     &'@    PAS DE TEMPS FIXE (IDTVAR = ',I10   ,')                 ',/,
     &'@                                                            ',/,
     &'@  Le calcul sera engage, mais le pas de temps ne sera pas   ',/,
     &'@    clippe. Le code indiquera juste les eventuels           ',/,
     &'@    depassements de critere local.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2541 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@  ON DEMANDE UNE LIMITATION DU PAS DE TEMPS LIEE AUX EFFETS ',/,
     &'@    DE DENSITE (IPTLRO = ',I10   ,') AVEC UN ALGORITHME     ',/,
     &'@    STATIONNAIRE (IDTVAR = ',I10   ,')                      ',/,
     &'@                                                            ',/,
     &'@  Le calcul sera engage, l''option IPTLRO ignoree.          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2600 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0, 10, 20, 21, 30, 31,',/,
     &'@    40, 41, 50, OU 60'                                       ,/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2605 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LA METHODE DES VORTEX NE PEUT ETRE ACTIVEE QUE POUR UN  ',/,
     &'@    CALCUL MONOPHASIQUE.                                    ',/,
     &'@    ON A ICI                                                ',/,
     &'@    NPHAS =',I10                                             ,/,
     &'@    IVRETX=',I10                                             ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2606 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LA METHODE DES VORTEX NE CONCERNE QUE LES CALCULS EN LES',/,
     &'@    ON A ICI                                                ',/,
     &'@    ITURB(1)=',I10                                           ,/,
     &'@    IVRETX  =',I10                                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute en ignorant le mot clef IVRTEX.    ',/,
     &'@    (il est repositione a 0)                                ',/,
     &'@                                                            ',/,
     &'@  Il est conseille de verifier les parametres donnes via    ',/,
     &'@  l''interface ou usini1.                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2607 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ON DEMANDE LA REDUCTION DU VOISINAGE ETENDU POUR LE     ',/,
     &'@    CALCUL DES GRADIENTS PAR MOINDRE CARRES, POUR UN        ',/,
     &'@    CALCUL EN L.E.S AVEC LE MODELE DYNAMIQUE.               ',/,
     &'@    MODELE DYNAMIQUE.                                       ',/,
     &'@      PHASE  ',I10                                           ,/,
     &'@      ITURB =',I10                                           ,/,
     &'@      IMRGRA=',I10                                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera engage.                                    ',/,
     &'@                                                            ',/,
     &'@  Le calcul de la moyenne locale de la constante de         ',/,
     &'@    Smagorinsky dynamique peut etre degrade.                ',/,
     &'@                                                            ',/,
     &'@  Il est conseille :                                        ',/,
     &'@    - soit d''utiliser une methode de calculs des gradients ',/,
     &'@      par moindre carres sur voisinage etendu complet       ',/,
     &'@      (IMRGRA = 2)                                          ',/,
     &'@    - soit de calculer sa propre moyenne de la constante    ',/,
     &'@      dynamique dans la subroutine USSMAG.                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2610 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10   ,' : DONNEES INCOHERENTES                 ',/,
     &'@    ',A31,I10                                                ,/,
     &'@    ',A31,I10                                                ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
c 2620 FORMAT
c     & (' ATTENTION : ON DEMANDE ',A6,' = ',I10,/,
c     &  '                  AVEC GRAVITE = ',3E14.5)
 2621 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@  ON DEMANDE LA PRISE EN COMPTE DE L''ACCELERATION DE LA    ',/,
     &'@    PESANTEUR ',3E14.5                                       ,/,
     &'@    SANS RESOLUTION D''UNE VARIABLE TEMPERATURE OU ENERGIE  ',/,
     &'@    (ISCALT = ',I10   ,')                                   ',/,
     &'@                                                            ',/,
     &'@  Le calcul sera engage.                                    ',/,
     &'@                                                            ',/,
     &'@  Il n''y a pas d''incompatibilite a prendre en compte      ',/,
     &'@    l''acceleration de la pesanteur sans effets thermiques, ',/,
     &'@    mais, souvent, la masse volumique depend essentiellement',/,
     &'@    de la temperature et la combinaison des options du      ',/,
     &'@    present calcul est caracteristique d''un oubli.         ',/,
     &'@                                                            ',/,
     &'@  Il est conseille de verifier les parametres donnes via    ',/,
     &'@  l''interface ou usini1.                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2622 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@  ON DEMANDE LA PRISE EN COMPTE DES TERMES DE GRAVITE DANS  ',/,
     &'@    LES EQUATIONS DE LA TURBULENCE (',A6,' = ',I10   ,')    ',/,
     &'@    SANS RESOLUTION D''UNE VARIABLE TEMPERATURE OU ENERGIE  ',/,
     &'@                                                            ',/,
     &'@  Le calcul sera engage.                                    ',/,
     &'@                                                            ',/,
     &'@  Si des effets de gravite sont recherches, il convient de  ',/,
     &'@    s assurer que la masse volumique est variable.          ',/,
     &'@  Le nombre de Prandtl turbulent sera pris egal a 1.        ',/,
     &'@  Elle peut varier en fonction d autres grandeurs que       ',/,
     &'@    la temperature ou l enthalpie ; si c est le cas, ce     ',/,
     &'@    message pourra etre ignore ; sinon, verifier usini1     ',/,
     &'@    ou imposer une variation de masse volumique dans usphyv.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2623 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@                                                            ',/,
     &'@ Le coefficient RELAXV des variables de la turbulence a ete ',/,
     &'@ modifie alors que IKECOU ne vaut pas 0. Il vaut ici        ',/,
     &'@ - pour k                  : ',E12.4                         ,/,
     &'@ - pour epsilon (ou omega) : ',E12.4                         ,/,
     &'@                                                            ',/,
     &'@ La modification sera sans effet (RELAXV n''est utile que   ',/,
     &'@ si IKECOU=0)                                               ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2624 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@                                                            ',/,
     &'@ Le coefficient RELAXV des variables de la turbulence doit  ',/,
     &'@ etre un reel compris entre 0 et 1. Il vaut ici :           ',/,
     &'@ - pour k                  : ',E12.4                         ,/,
     &'@ - pour epsilon (ou omega) : ',E12.4                         ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2625 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@                                                            ',/,
     &'@ Le coefficient RELAXV de la pression doit etre un reel     ',/,
     &'@ compris entre 0 et 1. Il vaut ici : ',E12.4                 ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2630 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@                                                            ',/,
     &'@  On demande la prise en compte de l''amortissement de      ',/,
     &'@    Van Driest (IDRIES(IPHAS) = ', I10                       ,/,
     &'@    avec le modele dynamique (ITURB(IPHAS) = ',I10           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas realise.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2640 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@    ',A6,' DOIT ETRE UN REEL INCLUS DANS L''INTERVALLE [0;1]',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2650 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',A8                                            ,/,
     &'@    ',A6,'(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0 OU 1    ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  ICPSYR(I) est l indicateur de couplage du scalaire I avec ',/,
     &'@    SYRTHES.                                                ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2660 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    INCOHERENCE POUR LE COUPLAGE SYRTHES                    ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Aucun couplage avec SYRTHES n''a ete defini.              ',/,
     &'@  Le nombre de scalaires couples est cependant ',I10         ,/,
     &'@    (Le nombre de scalaires total est ici ',I10   ,')       ',/,
     &'@  Verifier le couplage SYRTHES-Noyau.                       ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2661 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    INCOHERENCE POUR LE COUPLAGE SYRTHES                    ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Un couplage avec SYRTHES a ete defini.                    ',/,
     &'@  Le couplage avec SYRTHES necessite de disposer d''un      ',/,
     &'@    scalaire couple (et d''un seul).                        ',/,
     &'@  Le nombre de scalaires total   est ici ',I10               ,/,
     &'@  Le nombre de scalaires couples est ici ',I10               ,/,
     &'@  Verifier le couplage SYRTHES-Noyau.                       ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2662 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    INCOHERENCE POUR LE COUPLAGE SYRTHES                    ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Un couplage avec SYRTHES a ete defini.                    ',/,
     &'@  Si le code est couple a SYRTHES, le scalaire couple doit  ',/,
     &'@    etre la temperature.                                    ',/,
     &'@  Le scalaire couple est ici le scalaire ',I10               ,/,
     &'@    Ce n''est pas une temperature car                       ',/,
     &'@                    ISCSTH(',I10   ,') = ',I10               ,/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@  Pour coupler un scalaire qui n''est pas la temperature,   ',/,
     &'@    contacter l''equipe de developpement.                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2663 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    INCOHERENCE POUR LE COUPLAGE SYRTHES                    ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Un couplage avec SYRTHES a ete defini.                    ',/,
     &'@  En compressible, si le code est couple a SYRTHES, le      ',/,
     &'@    scalaire couple doit etre l''energie.                   ',/,
     &'@  Le scalaire couple est ici le scalaire ',I10               ,/,
     &'@    Ce n''est pas l''energie car                            ',/,
     &'@                    ISCSTH(',I10   ,') = ',I10               ,/,
     &'@                                                            ',/,
     &'@  Contacter l''equipe de developpement.                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2664 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    CHOIX DU CALCUL DES ESTIMATEURS D''ERREUR               ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Pour la phase ',I10   ,' l''indicateur IESCAL relatif a   ',/,
     &'@    l''estimateur d''erreur numero IEST = ',I10   ,' pour   ',/,
     &'@    Navier-Stokes doit etre un entier egal a 0, 1 ou 2.     ',/,
     &'@  Il vaut ici : IESCAL(',I10  ,',',I10   ,') = ',I10         ,/,
     &'@                                                            ',/,
     &'@  Rq. : les valeurs possibles de IEST sont :                ',/,
     &'@        IESPRE = ',I10                                       ,/,
     &'@        IESDER = ',I10                                       ,/,
     &'@        IESCOR = ',I10                                       ,/,
     &'@        IESTOT = ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2700 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    CHOIX DU MODE DE CALCUL DE LA DISTANCE A LA PAROI       ',/,
     &'@                                                            ',/,
     &'@  ICDPAR DOIT ETRE UN ENTIER EGAL A -2, -1, 1 ou 2          ',/,
     &'@  IL VAUT ICI ', I10                                         ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2710 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN REEL INCLUS DANS LINTERVALLE [0.;1.]',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2720 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN REEL SUPERIEUR OU EGAL A 1.         ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2730 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN REEL EGAL A 0. OU A 1.              ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2740 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN REEL STRICTEMENT POSITIF.           ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2750 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN ENTIER INFERIEUR OU EGAL A 1        ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    ',A6,'(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0 OU 1    ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@  L algorithme de resolution multigrille algebrique         ',/,
     &'@  n est disponible que pour la variable pression.           ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VARIABLE ',A8                                            ,/,
     &'@    ',A6,'(',I10   ,') DOIT ETRE UN ENTIER                  ',/,
     &'@      STRICTEMENT POSITIF                                   ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN ENTIER STRICTEMENT POSITIF          ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 4100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ENTREE DES DONNEES                          ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@    LA VITESSE DE REFERENCE UREF N''A PAS ETE INITIALISEE   ',/,
     &'@    OU A ETE MAL INITIALISEE (VALEUR NEGATIVE).             ',/,
     &'@    ELLE VAUT ICI ',E14.5                                    ,/,
     &'@                                                            ',/,
     &'@  Le calcul pourra etre execute si la turbulence est        ',/,
     &'@  initialisee a partir d''un fichier suite de calcul ou par ',/,
     &'@  la routine usiniv.                                        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 4300 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',A8                                            ,/,
     &'@    ISCSTH(',I10   ,') DOIT ETRE UN ENTIER EGAL A -1,0,1,2,3',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4310 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',A8                                            ,/,
     &'@    IPHSCA(',I10   ,') DOIT ETRE UN ENTIER                  ',/,
     &'@      STRICTEMENT POSITIF ET                                ',/,
     &'@      INFERIEUR OU EGAL A NPHAS = ',I10                      ,/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  IPHSCA est le numero de la phase porteuse du scalaire.    ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4320 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',A8                                            ,/,
     &'@    ISCAVR(',I10   ,') DOIT ETRE UN ENTIER                  ',/,
     &'@      POSITIF OU NUL ET                                     ',/,
     &'@      INFERIEUR OU EGAL A NSCAL = ',I10                      ,/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Si ISCAVR(I) est nul, le scalaire I n est pas une variance',/,
     &'@  Si ISCAVR(I) est positif, le scalaire I est une variance :',/,
     &'@    il s agit de la variance des fluctuations du scalaire J ',/,
     &'@    dont le numero est ISCAVR(I)                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4321 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LE SCALAIRE ',A8  ,'EST DEFINI COMME LA FLUCTUATION     ',/,
     &'@    DU SCALAIRE ',A8                                         ,/,
     &'@    (ISCAVR(',I10   ,') = ',I10   ,'),                      ',/,
     &'@    QUI EST LUI-MEME DEFINI COMME UNE FLUCTUATION           ',/,
     &'@    (ISCAVR(',I10   ,') = ',I10   ,').                      ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Si ISCAVR(I) est positif, le scalaire I est une variance :',/,
     &'@    il s agit de la variance des fluctuations du scalaire J ',/,
     &'@    dont le numero est ISCAVR(I) et on a donc forcement     ',/,
     &'@    ISCAVR(J) = 0                                           ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4330 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',A8                                            ,/,
     &'@    ICLVFL(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0, 1 OU 2 ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  ICLVFL(I) indique le mode de clipping du scalaire I       ',/,
     &'@    lorsqu il s agit d une variance de fluctuations.        ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4331 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',A8                                            ,/,
     &'@    ICLVFL(',I10   ,') N EST UTILISE QUE POUR LES VARIANCES ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@    ALORS QUE LE SCALAIRE N EST PAS UNE VARIANCE.           ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  ICLVFL(I) indique le mode de clipping du scalaire I       ',/,
     &'@    lorsqu il s agit d une variance de fluctuations.        ',/,
     &'@    Il n est pas utilise pour les autres scalaires.         ',/,
     &'@  L utilisateur est invite a ne pas modifier ICLVFL pour    ',/,
     &'@    les scalaires qui ne sont pas des variances.            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4340 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',A8                                            ,/,
     &'@    VISLS0(',I10   ,') DOIT ETRE UN REEL POSITIF            ',/,
     &'@      QUAND IVISLS(',I10   ,') = ',I10                       ,/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  VISLS0(I) est le coefficient de diffusion moleculaire du  ',/,
     &'@    scalaire et doit etre positif quand ivisls est different',/,
     &'@    de 1 (dans ce cas, un coefficient variable est donne    ',/,
     &'@    dans usphyv).                                           ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4350 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',A8                                            ,/,
     &'@    SIGMAS(',I10   ,') DOIT ETRE UN REEL POSITIF            ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  SIFMAS(I) est le nombre de Prandtl turbulent associe      ',/,
     &'@    au scalaire I.                                          ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4360 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',A8                                            ,/,
     &'@    SCAMIN(',I10   ,') VAUT ICI ',E14.5                      ,/,
     &'@      AVEC ICLVFL(',I10   ,') = ',I10                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  SCAMIN(I) est la valeur minimale acceptee pour le         ',/,
     &'@    scalaire I. Lorsque le scalaire est une variance        ',/,
     &'@    (ISCAVR(I) > 0) la valeur de SCAMIN n est prise en      ',/,
     &'@    compte que si ICLVFL(I) = 2                             ',/,
     &'@  Si l utilisateur souhaite effectivement que le            ',/,
     &'@    scalaire I (en fait, une variance) soit limite a SCAMIN ',/,
     &'@    (positif) il faut imposer ICLVFL = 2 dans usini1.       ',/,
     &'@  Si l utilisateur souhaite utiliser l option ICLVFL = 1    ',/,
     &'@    il est invite a ne pas modifier SCAMIN dans usini1.     ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4361 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',A8                                            ,/,
     &'@    SCAMAX(',I10   ,') VAUT ICI ',E14.5                      ,/,
     &'@      AVEC ICLVFL(',I10   ,') = ',I10                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  SCAMAX(I) est la valeur maximale acceptee pour le         ',/,
     &'@    scalaire I. Lorsque le scalaire est une variance        ',/,
     &'@    (ISCAVR(I) > 0) la valeur de SCAMAX n est prise en      ',/,
     &'@    compte que si ICLVFL(I) = 2                             ',/,
     &'@  Si l utilisateur souhaite effectivement que le            ',/,
     &'@    scalaire I (en fait, une variance) soit limite a SCAMAX ',/,
     &'@    (positif) il faut imposer ICLVFL = 2 dans usini1.       ',/,
     &'@  Si l utilisateur souhaite utiliser l option ICLVFL = 1    ',/,
     &'@    il est invite a ne pas modifier SCAMAX dans usini1.     ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4370 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',A8                                            ,/,
     &'@    SCAMAX(',I10   ,') DOIT ETRE UN REEL STRICTEMENT POSITIF',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  SCAMAX(I) est la valeur maximale acceptee pour le         ',/,
     &'@    scalaire I, ici une variance                            ',/,
     &'@  Avec ICLVFL(I) = 2, la valeur de SCAMAX doit donc etre    ',/,
     &'@   strictment positive.                                     ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4380 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',A8                                            ,/,
     &'@    RVARFL(',I10   ,') DOIT ETRE UN REEL STRICTEMENT POSITIF',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  RVARFL(I) est le coefficient R pour le scalaire I (qui est',/,
     &'@    une variance) intervenant dans le terme de dissipation :',/,
     &'@    - (1/R) rho scalaire epsilon/k                          ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 5003 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LA PERIODICITE DE ROTATION N''EST PAS COMPATIBLE AVEC LE',/,
     &'@      COUPLAGE VITESSE PRESSION RENFORCE  OU LA METHODE     ',/,
     &'@      ALE DANS LA VERSION COURANTE                          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  La variable COMMANDE_PERIO a ete renseignee dans le       ',/,
     &'@    lanceur (la periodicite a ete activee, ce qui se traduit',/,
     &'@    par IPERIO = ',I10,   ')                                ',/,
     &'@    et certaines periodicites sont de rotation.             ',/,
     &'@  L''indicateur IPUCOU a ete positionne a ',I10              ,/,
     &'@    dans l''interface ou usini1 (couplage renforce pour     ',/,
     &'@    IPUCOU=1).                                              ',/,
     &'@  L''indicateur IALE a ete positionne a ',I10                ,/,
     &'@    dans l''interface ou usini1 (methode activee si IALE=1) ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 5005 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    MODE DE CALCUL DE LA DISTANCE A LA PAROI INCOMPATIBLE   ',/,
     &'@      AVEC LA PERIODICITE                                   ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  La variable COMMANDE_PERIO a ete renseignee dans le       ',/,
     &'@    lanceur (la periodicite a ete activee, ce qui se traduit',/,
     &'@    par IPERIO = ',I10,   ').                               ',/,
     &'@  Les parametres de calcul specifies necessitent le calcul  ',/,
     &'@    la distance a la paroi (Rij-eps LRR avec echo de paroi, ',/,
     &'@    LES avec amortissement de van Driest ou k-omega SST).   ',/,
     &'@  Le mode de calcul de la distance a la paroi defini par    ',/,
     &'@    ICDPAR = ',I10,   ' ne permet pas de prendre en compte  ',/,
     &'@    la periodicite.                                         ',/,
     &'@                                                            ',/,
     &'@  Utiliser ICDPAR = 1 ou -1.                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
c 5007 FORMAT(
c     &'@                                                            ',/,
c     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
c     &'@                                                            ',/,
c     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
c     &'@    *********                                               ',/,
c     &'@    LA PERIODICITE             N''EST PAS COMPATIBLE AVEC LE',/,
c     &'@      RAYONNEMENT SEMI-TRANSPARENT  (ORDONNEES DISCRETES)   ',/,
c     &'@      DANS LA VERSION COURANTE                              ',/,
c     &'@                                                            ',/,
c     &'@  Le calcul ne peut etre execute.                           ',/,
c     &'@                                                            ',/,
c     &'@  L''indicateur IPERIO a ete positionne a ',I10              ,/,
c     &'@    dans usini1 (periodicite activee pour IPERIO=1).        ',/,
c     &'@  L''indicateur IRAYON(',I10,') a ete positionne a ',I10     ,/,
c     &'@    dans usray1.                                            ',/,
c     &'@                                                            ',/,
c     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
c     &'@                                                            ',/)
 5008 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LA PERIODICITE DE ROTATION N''EST PAS COMPATIBLE AVEC LE',/,
     &'@      RAYONNEMENT SEMI-TRANSPARENT  (ORDONNEES DISCRETES)   ',/,
     &'@      DANS LA VERSION COURANTE                              ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  La variable COMMANDE_PERIO a ete renseignee dans le       ',/,
     &'@    lanceur (la periodicite a ete activee, ce qui se traduit',/,
     &'@    par IPERIO = ',I10,   ')                                ',/,
     &'@    et certaines periodicites sont de rotation.             ',/,
     &'@  L''indicateur IRAYON(',I10,') a ete positionne a ',I10     ,/,
     &'@    dans l''interface ou usray1.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 5009 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    DES DEFAUTS PEUVENT SE RENCONTRER LORS DE L''UTILISATION',/,
     &'@      DE LA PERIODICITE DE ROTATION EN RIJ-EPSILON.         ',/,
     &'@                                                            ',/,
     &'@  La variable COMMANDE_PERIO a ete renseignee dans le       ',/,
     &'@    lanceur (la periodicite a ete activee, ce qui se traduit',/,
     &'@    par IPERIO = ',I10,   ')                                ',/,
     &'@    et certaines periodicites sont de rotation.             ',/,
     &'@  L''indicateur ITURB(IPHAS) a ete positionne a ',I10        ,/,
     &'@    pour la phase IPHAS = ',I10   ,' dans usini1.           ',/,
     &'@                                                            ',/,
     &'@  Le calcul peut etre execute.                              ',/,
     &'@    Les defauts eventuels evoques proviennent de la prise en',/,
     &'@    compte de la rotation du tenseur de viscosite orthotrope',/,
     &'@    Il a cependant en general une influence faible de sorte ',/,
     &'@    que les tests menes jusqu''a present n''ont pas fait    ',/,
     &'@    apparaitre de probleme.                                 ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 5010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    DES DEFAUTS PEUVENT SE RENCONTRER LORS DE L''UTILISATION',/,
     &'@      DE LA PERIODICITE DE ROTATION AVEC UN SCHEMA D''ORDRE ',/,
     &'@      DEUX EN TEMPS POUR LA VITESSE.                        ',/,
     &'@                                                            ',/,
     &'@  La variable COMMANDE_PERIO a ete renseignee dans le       ',/,
     &'@    lanceur (la periodicite a ete activee, ce qui se traduit',/,
     &'@    par IPERIO = ',I10,   ')                                ',/,
     &'@    et certaines periodicites sont de rotation.             ',/,
     &'@  Les indicateurs THETAV des trois composantes Ux, Uy, Uz   ',/,
     &'@    de la vitesse pour la phase IPHAS = ',I10                ,/,
     &'@    ont ete positionnes (dans usini1 ou par defaut suite aux',/,
     &'@    options de calcul selectionnees) aux valeurs suivantes :',/,
     &'@    THETAV(IU(IPHAS))  THETAV(IV(IPHAS))  THETAV(IW(IPHAS)) ',/,
     &'@        ',E14.5     ,'     ',E14.5     ,'     ',E14.5        ,/,
     &'@                                                            ',/,
     &'@  Le calcul peut etre execute.                              ',/,
     &'@    Les defauts eventuels evoques proviennent de la prise en',/,
     &'@    compte de la rotation du tenseur gradient de vitesse.   ',/,
     &'@    En effet, on met en oeuvre une methode explicite qui est',/,
     &'@    susceptible de faire chuter l''ordre du schema en temps.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 6002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LE PARALLELISME N''EST PAS COMPATIBLE AVEC LE MODULE    ',/,
     &'@      LAGRANGIEN DANS LA VERSION COURANTE                   ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Le processeur courant est de rang ',I10                    ,/,
     &'@  L''indicateur IILAGR a ete positionne a ',I10              ,/,
     &'@    dans uslag1 (module lagrangien active pour IILAGR>0).   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 6005 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    MODE DE CALCUL DE LA DISTANCE A LA PAROI INCOMPATIBLE   ',/,
     &'@      AVEC LE PARALLELISME                                  ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Le processeur courant est de rang ',I10                    ,/,
     &'@  Les parametres de calcul specifies necessitent le calcul  ',/,
     &'@    la distance a la paroi (Rij-eps LRR avec echo de paroi, ',/,
     &'@    LES avec amortissement de van Driest ou k-omega SST).   ',/,
     &'@  Le mode de calcul de la distance a la paroi defini par    ',/,
     &'@    ICDPAR = ',I10,   ' ne permet pas de prendre en compte  ',/,
     &'@    le parallelisme.                                        ',/,
     &'@                                                            ',/,
     &'@  Utiliser ICDPAR = 1 ou -1.                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    INDICATEUR DE METHODE ALE                               ',/,
     &'@                                                            ',/,
     &'@  IALE DOIT VALOIR 0 OU 1                                   ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usalin.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    NOMBRE D''ITERATIONS D''INITIALISATION DU FLUIDE EN ALE ',/,
     &'@                                                            ',/,
     &'@  NALINF DOIT ETRE UN ENTIER POSITIF                        ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usalin.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    COEFFICIENTS DE LA METHODE DE NEWMARK NON VALIDES       ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  ALPNMK doit etre compris entre 0 et 1                     ',/,
     &'@  BETNMK doit etre compris entre 0 et 1/2                   ',/,
     &'@  GAMNMK doit etre compris entre 0 et 1                     ',/,
     &'@  On a ici :                                                ',/,
     &'@                                                            ',/,
     &'@       ALPNMK      BETNMK      GAMNMK                       ',/,
     &'@ ',     E12.4,      E12.4,      E12.4                        ,/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface, usini1   ',/,
     &'@    ou usalin.                                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7030 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    NOMBRE D''ITERATIONS MAX DE COUPLAGE IMPLICITE EN ALE   ',/,
     &'@                                                            ',/,
     &'@  NALIMX DOIT ETRE UN ENTIER POSITIF                        ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usalin.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7040 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PRECISION DU COUPLAGE IMPLICITE EN ALE                  ',/,
     &'@                                                            ',/,
     &'@  EPALIM DOIT ETRE UN REEL STRICTEMENT POSITIF              ',/,
     &'@    IL VAUT ICI ',E14.5                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usalin.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7050 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ITERATION D''INITIALISATION DE L''ALE                   ',/,
     &'@                                                            ',/,
     &'@  ITALIN DOIT VALOIR 0 OU 1                                 ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes dans usalin.               ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********   MODULE COMPRESSIBLE                         ',/,
     &'@                                                            ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@    T0 ET P0 DOIVENT ETRE DES REELS STRICTEMENT POSITIFS    ',/,
     &'@    ILS VALENT ICI :                                        ',/,
     &'@                   T0 = ',E14.5                              ,/,
     &'@                   P0 = ',E14.5                              ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C=======================================================================
C 8. SORTIE
C=======================================================================
C
      RETURN
      END
c@z
