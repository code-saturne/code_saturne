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
                        SUBROUTINE VARPOS
C                       *****************
C     -------------------------------------------------------------
     &(NMODPP )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  INIT DES POSITIONS DES VARIABLES SELON
CFONC    LE TYPE DE CALCUL INDIQUE PAR L'UTILISATEUR
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NMODPP       ! E  ! <-  ! NOMBRE DE MODELES PHYS.PART. ACTIVES !
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
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "albase.h"
      INCLUDE "parall.h"
      INCLUDE "radiat.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagdim.h"
      INCLUDE "lagran.h"
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
      INTEGER       NMODPP
C
C VARIABLES LOCALES
C
C
      CHARACTER     RUBRIQ*64, CMOY4*4, CINDFM*4
      INTEGER       IVAR  , IPP   , ISCAL , IPHAS , IPROP , IPROFL
      INTEGER       II    , JJ    , KK    , LL
      INTEGER       IOK   , IPPOK , IPPPST
      INTEGER       IFLUM , ICONDL
      INTEGER       IEST  , IIFLAA, IVISPH, IPROFA, IPROPP
      INTEGER       IMOM  , JMOM  , IMOLD , JMOLD , JBMOMT, JTCABS
      INTEGER       IIPLUS, IIMOIN, JMOMOK, IDTO  , JDTO
      INTEGER       IERROR, IRTYP , ITYSUP, NBVAL
      INTEGER       IDFFIN, IDFMJI, IMOMR , ILECEC
      INTEGER       ILSMOM, IVERS , INMFIN
      INTEGER       IMPAMX
      INTEGER       NCELOK, NFAIOK, NFABOK, NSOMOK
      INTEGER       NFMTMO, NBERRO
      INTEGER       IDTOLD(NBMOMX)
C
      DOUBLE PRECISION GRAVN2
C
      INTEGER       IPASS
      DATA          IPASS /0/
      SAVE          IPASS
C
      INTEGER       NPRMAX
      DATA          NPRMAX /0/
      SAVE          NPRMAX
C
      INTEGER       NPPMAX
      DATA          NPPMAX /0/
      SAVE          NPPMAX
C
C=======================================================================
C 0. INITIALISATIONS
C=======================================================================
C
      IPASS = IPASS + 1
C
C     Nombre max pour les formats choisis
      NFMTMO = 9999
C     Indefini a 4 caracteres
      CINDFM = 'YYYY'
C
C=======================================================================
C 1. PREMIER APPEL : VERIFICATION DU NOMBRE DE PHASES
C                    RETURN
C=======================================================================
C
      IF(IPASS.EQ.1) THEN
C
        IOK = 0
C
        IF(NPHAS.LE.0) THEN
          WRITE(NFECRA,5000) NPHAS
          IOK = IOK + 1
        ENDIF
        IF(NPHAS.GT.NPHSMX) THEN
          WRITE(NFECRA,5001) NPHAS, NPHSMX, NPHAS
          IOK = IOK + 1
        ENDIF
C
        IF(IOK.NE.0) THEN
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
        RETURN
C
      ENDIF
C
C=======================================================================
C 2. SECOND APPEL : CALCUL DE NSCAPP
C                   VERIFICATION DU NOMBRE DE SCALAIRES
C                   CONSTRUCTION DE ISCAPP
C                   CALCUL DE NSCAL
C                   RETURN
C
C  C'est juste avant ce second appel que les modeles de combustion
C    auront ete renseignes. C'est dans la section ci-dessous qu'on en
C    en deduira NSCAPP (avant meme les verifications).
C  A la sortie de cette section, NSCAL, NSCAUS et NSCAPP sont connus.
C  On renseignera egalement ici les valeurs de IPHSCA, ISCAVR, IVISLS
C    pour les scalaires physiques particulieres en question.
C  On en profite aussi pour remplir ITYTUR puisque ITURB vient d'etre
C    defini.
C=======================================================================
C
      IF(IPASS.EQ.2) THEN
C
C ---> Remplissage de ITYTUR
        DO IPHAS = 1, NPHAS
          ITYTUR(IPHAS) = ITURB(IPHAS)/10
        ENDDO
C
C ---> Coherence modele
C     Rq : ATTENTION il faudrait renforcer le blindage
C
        IOK   = 0
        NMODPP = 0
        DO IPP = 1, NMODMX
          IF ( IPPMOD(IPP).NE.-1 ) THEN
            NMODPP = NMODPP+1
            IPPOK = IPP
          ENDIF
        ENDDO
        IF ( NMODPP.GT.1 ) THEN
          WRITE(NFECRA,6000)
          IOK = IOK + 1
        ENDIF
C
        IF ( NMODPP.EQ.1 ) THEN
          IF ( IPPMOD(IPPOK).LT.0 .OR. IPPMOD(IPPOK).GT.5 ) THEN
            WRITE(NFECRA,6001)
            IOK = IOK + 1
          ENDIF
        ENDIF
C
        IF(IOK.NE.0) THEN
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C
C ---> Lecture donnees thermochimie
C
        CALL PPLECD
C       ===========
C
C ---> Calcul de NSCAPP
C
        CALL PPCSCA
C       ===========
C
C
C ---> Verifications
C
        IOK = 0
C
        IF(NSCAUS.LT.0) THEN
          WRITE(NFECRA,6010) NSCAUS
          IOK = IOK + 1
        ENDIF
C
        IF(NSCAUS.GT.0 .OR. NSCAPP.GT.0) THEN
          IF((NSCAUS+NSCAPP).GT.NSCAMX) THEN
C
            IF(NSCAPP.LE.0) THEN
              WRITE(NFECRA,6011)
     &             NSCAUS,       NSCAMX,NSCAMX       ,NSCAUS
            ELSE
              WRITE(NFECRA,6012)
     &             NSCAUS,NSCAPP,NSCAMX,NSCAMX-NSCAPP,NSCAUS+NSCAPP
            ENDIF
            IOK = IOK + 1
          ENDIF
        ENDIF
C
        IF(IOK.NE.0) THEN
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C ---> Calcul de ISCAPP et NSCAL
C      On prefere que l'identite porte sur les scalaires utilisateurs,
C        ca minimisera peut etre des erreurs utilisateur
C
        ISCAL = MAX(0,NSCAUS)
        IF (NSCAPP.GT.0) THEN
          DO II = 1, NSCAPP
            ISCAL = ISCAL + 1
            ISCAPP(II) = ISCAL
          ENDDO
C
          CALL PPVARP
C         ===========
C
        ENDIF
C
        NSCAL = ISCAL
C
        RETURN
C
      ENDIF
C
C
C=======================================================================
C 3. TROISIEME APPEL : VERIFICATIONS ET
C        POSITIONNEMENT DES VARIABLES  : IPR, IU ... ISCA, NVAR
C                       ET DES PROPRIETES PPALES
C=======================================================================
C
      IF(IPASS.EQ.3) THEN
C
C
C ---> 3.1 VERIFICATIONS
C      -----------------
C
        IOK = 0
C
C ---  NPHAS et NSCAL ont deja ete verifies, mais on ne sait jamais.
C
        IF(NPHAS.LE.0) THEN
          WRITE(NFECRA,5000) NPHAS
          IOK = IOK + 1
        ENDIF
        IF(NPHAS.GT.NPHSMX) THEN
          WRITE(NFECRA,5001) NPHAS, NPHSMX, NPHAS
          IOK = IOK + 1
        ENDIF
C
        IF(NSCAL.LT.0) THEN
          WRITE(NFECRA,7010) NSCAL, NSCAUS, NSCAPP
          IOK = IOK + 1
        ENDIF
        IF(NSCAL.GT.NSCAMX) THEN
          WRITE(NFECRA,7011) NSCAL, NSCAMX, NSCAUS, NSCAPP, NSCAL
          IOK = IOK + 1
        ENDIF
C
C --- IPHSCA(ISCAL) doit etre compris entre 0 et NPHAS.
C
        IF(NSCAUS.GT.0) THEN
          DO II = 1, NSCAUS
            ISCAL = II
            IF(IPHSCA(ISCAL).GT.NPHAS.OR.IPHSCA(ISCAL).LT.0) THEN
              WRITE(NFECRA,7020) ISCAL, ISCAL, II,II,
     &             IPHSCA(ISCAL), NPHAS
              IOK = IOK + 1
            ENDIF
          ENDDO
        ENDIF
        IF(NSCAPP.GT.0) THEN
          DO II = 1, NSCAPP
            ISCAL = ISCAPP(II)
            IF(IPHSCA(ISCAL).GT.NPHAS.OR.IPHSCA(ISCAL).LT.0) THEN
              WRITE(NFECRA,7021) ISCAL, ISCAL, II,II,
     &             IPHSCA(ISCAL), NPHAS
              IOK = IOK + 1
            ENDIF
          ENDDO
        ENDIF
C
C --- ISCAVR(ISCAL) doit etre compris entre 0 et NSCAL.
C
        IF(NSCAUS.GT.0) THEN
          DO II = 1, NSCAUS
            ISCAL = II
            IF(ISCAVR(ISCAL).GT.NSCAL.OR.ISCAVR(ISCAL).LT.0) THEN
              WRITE(NFECRA,7030) ISCAL,II,II,ISCAVR(ISCAL),NSCAL
              IOK = IOK + 1
            ENDIF
          ENDDO
        ENDIF
        IF(NSCAPP.GT.0) THEN
          DO II = 1, NSCAPP
            ISCAL = ISCAPP(II)
            IF(ISCAVR(ISCAL).GT.NSCAL.OR.ISCAVR(ISCAL).LT.0) THEN
              WRITE(NFECRA,7031) ISCAL,II,II,ISCAVR(ISCAL),NSCAL
              IOK = IOK + 1
            ENDIF
          ENDDO
        ENDIF
C
C
C --- IVISLS(ISCAL) doit etre non initialise pour les variances
C     Il prend la valeur du scalaire associe
C     Tous les tests qui suivent sont utilises pour le message
C       d'erreur eventuel.
C
        IF(NSCAUS.GT.0) THEN
          DO JJ = 1, NSCAUS
            II    = JJ
            ISCAL = ISCAVR(II)
C       Si on a une variance avec ivisls initialise : erreur
            IF    ( (ISCAL.GT.0.AND.ISCAL.LE.NSCAL).AND.
     &               IVISLS(II).NE.-1                        ) THEN
              LL = 0
              DO KK = 1, NSCAUS
                IF(       KK .EQ.ISCAL) LL = KK
              ENDDO
              DO KK = 1, NSCAPP
                IF(ISCAPP(KK).EQ.ISCAL) LL = -KK
              ENDDO
              IF(LL.GT.0) THEN
                WRITE(NFECRA,7040)
     &               II,II,JJ,ISCAL,LL,JJ,ISCAL,JJ,IVISLS(ISCAL)
              ELSE
                WRITE(NFECRA,7041)
     &               II,II,JJ,ISCAL,-LL,JJ,ISCAL,JJ,IVISLS(ISCAL)
              ENDIF
              IOK = IOK + 1
C     Si on n'a pas une variance std mais que ivisls est incorrect : erreur
            ELSEIF( (ISCAL.LE.0 .OR.ISCAL.GT.NSCAL).AND.
     &             (IVISLS(II).NE.0.AND.IVISLS(II).NE.1) ) THEN
              WRITE(NFECRA,7050) II,JJ,JJ,IVISLS(II)
              IOK = IOK + 1
            ENDIF
          ENDDO
        ENDIF
C
        IF(NSCAPP.GT.0) THEN
          DO JJ = 1, NSCAPP
            II    = ISCAPP(JJ)
            ISCAL = ISCAVR(II)
C       Si on a une variance avec ivisls initialise : erreur
            IF    ( (ISCAL.GT.0.AND.ISCAL.LE.NSCAL).AND.
     &               IVISLS(II).NE.-1                        ) THEN
              LL = 0
              DO KK = 1, NSCAUS
                IF(       KK .EQ.ISCAL) LL = KK
              ENDDO
              DO KK = 1, NSCAPP
                IF(ISCAPP(KK).EQ.ISCAL) LL = -KK
              ENDDO
              IF(LL.GT.0) THEN
                WRITE(NFECRA,7042)
     &               II,II,JJ,ISCAL,LL,JJ,ISCAL,JJ,IVISLS(ISCAL)
              ELSE
                WRITE(NFECRA,7043)
     &               II,II,JJ,ISCAL,-LL,JJ,ISCAL,JJ,IVISLS(ISCAL)
              ENDIF
              IOK = IOK + 1
C       Si on n'a pas une variance std mais que ivisls est incorrect : erreur
            ELSEIF( (ISCAL.LE.0 .OR.ISCAL.GT.NSCAL).AND.
     &              (IVISLS(II).NE.0.AND.IVISLS(II).NE.1) ) THEN
              WRITE(NFECRA,7051)II,JJ,JJ,IVISLS(II)
              IOK = IOK + 1
            ENDIF
          ENDDO
        ENDIF
C
C       On initialise les ivisls des variances
        IF(NSCAL.GT.0) THEN
          DO II = 1, NSCAL
            ISCAL = ISCAVR(II)
            IF(ISCAL.GT.0.AND.ISCAL.LE.NSCAL) THEN
              IVISLS(II) = IVISLS(ISCAL)
            ENDIF
          ENDDO
        ENDIF
C
C ---> IPHSCA
C      Pour les variances de fluctuations, les valeurs de IPHSCA
C        ne doivent pas avoir ete modifiees par l'utilisateur
C        Elles sont prises egales aux valeurs correspondantes
C        pour le scalaire associe.
C
        IF(NSCAUS.GT.0) THEN
          DO JJ = 1, NSCAUS
            II    = JJ
            ISCAL = ISCAVR(II)
            IF(ISCAL.GT.0.AND.ISCAL.LE.NSCAL)THEN
              IF(IPHSCA(II).EQ.0) THEN
                IPHSCA(II) = IPHSCA(ISCAL)
              ELSE
                LL = 0
                DO KK = 1, NSCAUS
                  IF(       KK .EQ.ISCAL) LL = KK
                ENDDO
                DO KK = 1, NSCAPP
                  IF(ISCAPP(KK).EQ.ISCAL) LL = -KK
                ENDDO
                IF(LL.GT.0) THEN
                  WRITE(NFECRA,7060)II,
     &                 II,JJ,ISCAL,LL,JJ,ISCAL,
     &                 JJ,IPHSCA(ISCAL)
                ELSE
                  WRITE(NFECRA,7061)II,
     &                 II,JJ,ISCAL,-LL,JJ,ISCAL,
     &                 JJ,IPHSCA(ISCAL)
                ENDIF
                IOK = IOK + 1
              ENDIF
            ENDIF
          ENDDO
        ENDIF
C
        IF(NSCAPP.GT.0) THEN
          DO JJ = 1, NSCAPP
            II    = ISCAPP(JJ)
            ISCAL = ISCAVR(II)
            IF(ISCAL.GT.0.AND.ISCAL.LE.NSCAL)THEN
              IF(IPHSCA(II).EQ.0) THEN
                IPHSCA(II) = IPHSCA(ISCAL)
              ELSE
                LL = 0
                DO KK = 1, NSCAUS
                  IF(       KK .EQ.ISCAL) LL = KK
                ENDDO
                DO KK = 1, NSCAPP
                  IF(ISCAPP(KK).EQ.ISCAL) LL = -KK
                ENDDO
                IF(LL.GT.0) THEN
                  WRITE(NFECRA,7062)II,
     &                 II,JJ,ISCAL,LL,JJ,ISCAL,
     &                 JJ,IPHSCA(ISCAL)
                ELSE
                  WRITE(NFECRA,7063)II,
     &                 II,JJ,ISCAL,-LL,JJ,ISCAL,
     &                 JJ,IPHSCA(ISCAL)
                ENDIF
                IOK = IOK + 1
              ENDIF
            ENDIF
          ENDDO
        ENDIF
C
C ---> VISCOSITE ALE
        IF (IALE.EQ.1) THEN
          IF (IORTVM.NE.0 .AND. IORTVM.NE.1) THEN
            WRITE(NFECRA,7070) IORTVM
            IOK = IOK + 1
          ENDIF
        ENDIF
C --- On s'arrete si quelque chose s'est mal passe
C
        IF(IOK.NE.0) THEN
          CALL CSEXIT (1)
        ENDIF
C
C
C ---> 3.2 POSITIONNEMENT DES VARIABLES  : IPR, IU ... ISCA, NVAR
C      --------------------------------
C
        IVAR = 0
C
        DO IPHAS = 1, NPHAS
C
C --- Pression : supposons ici qu'il n'y a qu'une pression
C       quelque soit le nombre de phases.
C     Le reste du code devrait se preter a l'extension, du fait
C       qu'il n'y a qu'ici qu'on fait cette hypothese (evidemment, il
C       faudrait ecrire un schema a deux pressions et remplacer navsto,
C       mais c'est un "detail"...relativement a ce qui nous occupe dans
C       le present sous-programme)
C
          IF(IPHAS.EQ.1) THEN
            IVAR          = IVAR + 1
            IPR   (IPHAS) = IVAR
          ELSE
            IPR   (IPHAS) = IPR(1)
          ENDIF
C
C --- Vitesse
          IVAR          = IVAR + 1
          IU    (IPHAS) = IVAR
          IVAR          = IVAR + 1
          IV    (IPHAS) = IVAR
          IVAR          = IVAR + 1
          IW    (IPHAS) = IVAR
C
C --- Turbulence
          IF (ITYTUR(IPHAS).EQ.2) THEN
            IVAR          = IVAR + 1
            IK    (IPHAS) = IVAR
            IVAR          = IVAR + 1
            IEP   (IPHAS) = IVAR
          ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
            IVAR          = IVAR + 1
            IR11  (IPHAS) = IVAR
            IVAR          = IVAR + 1
            IR22  (IPHAS) = IVAR
            IVAR          = IVAR + 1
            IR33  (IPHAS) = IVAR
            IVAR          = IVAR + 1
            IR12  (IPHAS) = IVAR
            IVAR          = IVAR + 1
            IR13  (IPHAS) = IVAR
            IVAR          = IVAR + 1
            IR23  (IPHAS) = IVAR
            IVAR          = IVAR + 1
            IEP   (IPHAS) = IVAR
          ELSEIF(ITURB(IPHAS).EQ.50) THEN
            IVAR          = IVAR + 1
            IK    (IPHAS) = IVAR
            IVAR          = IVAR + 1
            IEP   (IPHAS) = IVAR
            IVAR          = IVAR + 1
            IPHI  (IPHAS) = IVAR
            IVAR          = IVAR + 1
            IFB   (IPHAS) = IVAR
          ELSEIF(ITURB(IPHAS).EQ.60) THEN
            IVAR          = IVAR + 1
            IK    (IPHAS) = IVAR
            IVAR          = IVAR + 1
            IOMG  (IPHAS) = IVAR
          ENDIF
C
        ENDDO
C
C --- Scalaires
        IF(NSCAPP.GE.1) THEN
          DO JJ = 1, NSCAPP
            II       = ISCAPP(JJ)
            IVAR     = IVAR + 1
            ISCA(II) = IVAR
          ENDDO
        ENDIF
        IF(NSCAUS.GE.1) THEN
          DO JJ = 1, NSCAUS
            II       = JJ
            IVAR     = IVAR + 1
            ISCA(II) = IVAR
          ENDDO
        ENDIF
C
C --- Vitesse de maillage en ALE
        IF (IALE.EQ.1) THEN
          IVAR = IVAR + 1
          IUMA = IVAR
          IVAR = IVAR + 1
          IVMA = IVAR
          IVAR = IVAR + 1
          IWMA = IVAR
        ENDIF
C
C --- Nombre total de variables
        NVAR = IVAR
C
C --- Verification de NVAR
C
        IF(NVAR.GT.NVARMX) THEN
          WRITE(NFECRA,7100)NVAR,NVARMX,NVAR
          CALL CSEXIT (1)
        ENDIF
C
C
C --- Maintenant on peut faire ceci :
C
        DO IPHAS = 1, NPHAS
          ISTAT (IPR(IPHAS)) = 0
          ICONV (IPR(IPHAS)) = 0
          IF (ITURB(IPHAS).EQ.50) THEN
            ISTAT(IFB(IPHAS))  = 0
            ICONV(IFB(IPHAS))  = 0
C     Pour fb, on sait qu'on a un terme diagonal, meme si ISTAT=0,
C       donc on ne decalera pas la diagonale
            IDIRCL(IFB(IPHAS)) = 0
          ENDIF
        ENDDO
        IF (IALE.EQ.1) THEN
          ISTAT(IUMA) = 0
          ICONV(IUMA) = 0
          ISTAT(IVMA) = 0
          ICONV(IVMA) = 0
          ISTAT(IWMA) = 0
          ICONV(IWMA) = 0
        ENDIF
C
C
C ---> 3.3 POSITIONNEMENT DES PROPRIETES PRINCIPALES
C      --------------------------------
C
C --- Numerotation des proprietes presentes ici
C       Ceci depend du type de solveur branche derriere
C        (CP, Poly, Mono...)
C       Dans l'ideal, il y aurait donc plusieurs varpos.
C
C       Pour le moment, on fait les hypotheses suivantes :
C         Il y a toujours, pour toutes les phases,  rho, viscl, visct
C         Il y a toujours la pression totale (sauf en compressible)
C         Lorsqu'elles sont variables, on a les proprietes suivantes :
C           . cp    (par phase)
C           . visls (par scalaire)
C           . csmago (par phase) en LES dynamique
C         En ALE on a la viscosite de maillage
C         On a aussi les flux de masse porteurs :
C           . les variables u,v,w,p,turbulence sont portees par leur
C               phase (1 flux)
C           . les scalaires sont portes par leur phase (1 flux)
C           On suppose donc qu'il n'y a pas de scalaire en
C             taux de vide a convecter avec un flux particulier (ce
C             serait possible : ca rajoute un test, par exemple
C             if alpro...
C
C     ATTENTION : on s'arrange pour numeroter a la queue-leu-leu sans
C       trous les proprietes qui sont definies au centre des cellules
C       ceci permet ensuite de ne pas se fatiguer lors de la
C       construction de IPPPRO plus bas.
C      Cependant, pour les physiques particulieres, ce n'est pas le cas.
C
C
        IPROP = 0
C
C   Proprietes des phases : proprietes toujours presentes
        DO IPHAS = 1, NPHAS
          IPROP         = IPROP + 1
          IROM  (IPHAS) = IPROP
          IPROP         = IPROP + 1
          IVISCL(IPHAS) = IPROP
          IPROP         = IPROP + 1
          IVISCT(IPHAS) = IPROP
          IPROP         = IPROP + 1
          ICOUR (IPHAS) = IPROP
          IPROP         = IPROP + 1
          IFOUR (IPHAS) = IPROP
        ENDDO
C
C  Pression totale stockee dans IPRTOT, si on n'est pas en compressible
C  (sinon Ptot=P* !)
        IF (IPPMOD(ICOMPF).LT.0) THEN
          DO IPHAS = 1, NPHAS
            IF (IPHAS.EQ.1) THEN
              IPROP         = IPROP + 1
              IPRTOT(1)     = IPROP
            ELSE
              IPRTOT(IPHAS) = IPRTOT(1)
            ENDIF
          ENDDO
        ENDIF
C
C  Proprietes des phases : CP s'il est variable
        DO IPHAS = 1, NPHAS
          IF(ICP(IPHAS).NE.0) THEN
            IPROP         = IPROP + 1
            ICP   (IPHAS) = IPROP
          ENDIF
        ENDDO
C
C  Proprietes des phases : Cs^2 si on est en LES dynamique
        DO IPHAS = 1, NPHAS
          IF(ITURB(IPHAS).EQ.41) THEN
            IPROP         = IPROP + 1
            ISMAGO(IPHAS) = IPROP
          ELSE
            ISMAGO(IPHAS) = -1
          ENDIF
        ENDDO
C
C  Viscosite de maillage en ALE
        IF (IALE.EQ.1) THEN
          IPROP     = IPROP + 1
          IVISMA(1) = IPROP
C     si la viscosite est isotrope, les trois composantes pointent
C       au meme endroit
          IF (IORTVM.EQ.0) THEN
            IVISMA(2) = IPROP
            IVISMA(3) = IPROP
          ELSE
            IPROP     = IPROP + 1
            IVISMA(2) = IPROP
            IPROP     = IPROP + 1
            IVISMA(3) = IPROP
          ENDIF
        ENDIF
C
C   Proprietes des phases : estimateurs d'erreur
        DO IPHAS = 1, NPHAS
          DO IEST = 1, NESTMX
            IPROP              = IPROP + 1
            IESTIM(IEST,IPHAS) = IPROP
          ENDDO
        ENDDO
C
C
C   Proprietes des scalaires : VISCLS si elle est variable
C     On utilisera IVISLS comme suit :
C       Pour le scalaire II
C         si IVISLS(II) = 0    : visc = VISLS0(II)
C         si IVISLS(II) .GT. 0 : visc = PROPCE(IEL ,IPPROC(IVISLS(II)))
C     Ceci permet de ne pas reserver des tableaux vides pour les
C       scalaires a viscosite constante
        IF(NSCAL.GE.1) THEN
          DO II = 1, NSCAL
            IF(IVISLS(II).NE.0) THEN
              IPROP      = IPROP + 1
              IVISLS(II) = IPROP
            ENDIF
          ENDDO
        ENDIF
C
C
C   Proprietes des variables : flux de masse porteur
C
        IPROFL = IPROP
        DO IPHAS = 1, NPHAS
          IPROP               = IPROP + 1
          IF(IPHAS.EQ.1) THEN
            IFLUMA(IPR (IPHAS)) = IPROP
          ENDIF
          IFLUMA(IU  (IPHAS)) = IPROP
          IFLUMA(IV  (IPHAS)) = IPROP
          IFLUMA(IW  (IPHAS)) = IPROP
          IF(ITYTUR(IPHAS).EQ.2) THEN
            IFLUMA(IK  (IPHAS)) = IPROP
            IFLUMA(IEP (IPHAS)) = IPROP
          ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
            IFLUMA(IR11(IPHAS)) = IPROP
            IFLUMA(IR22(IPHAS)) = IPROP
            IFLUMA(IR33(IPHAS)) = IPROP
            IFLUMA(IR12(IPHAS)) = IPROP
            IFLUMA(IR13(IPHAS)) = IPROP
            IFLUMA(IR23(IPHAS)) = IPROP
            IFLUMA(IEP (IPHAS)) = IPROP
          ELSEIF(ITURB(IPHAS).EQ.50) THEN
            IFLUMA(IK  (IPHAS)) = IPROP
            IFLUMA(IEP (IPHAS)) = IPROP
            IFLUMA(IPHI(IPHAS)) = IPROP
            IFLUMA(IFB (IPHAS)) = IPROP
          ELSEIF(ITURB(IPHAS).EQ.60) THEN
            IFLUMA(IK  (IPHAS)) = IPROP
            IFLUMA(IOMG(IPHAS)) = IPROP
          ENDIF
        ENDDO
        DO ISCAL = 1, NSCAL
          IFLUMA(ISCA(ISCAL)) = IFLUMA(IU(IPHSCA(ISCAL)))
        ENDDO
        IF (IALE.EQ.1) THEN
          IFLUMA(IUMA) = IFLUMA(IPR(1))
          IFLUMA(IVMA) = IFLUMA(IPR(1))
          IFLUMA(IWMA) = IFLUMA(IPR(1))
        ENDIF
C     Nombre total de flux de masse
C       IPROFL ressert plus bas.
        NFLUMA = IPROP - IPROFL
C
C     Numero max des proprietes ; ressert plus bas pour
C       ajouter celles relatives a la physique particuliere
        NPRMAX = IPROP
C
C
C --- Positionnement dans les tableaux PROPCE, PROFA, PROFB
C
C   Au centre des cellules (tout sauf les flux de masse)
C     Pour les fluctuations, le pointeur de la diffusivite
C     envoie directement sur la diffusivite du scalaire associe.
C
C  On positionne en meme temps les pointeurs IPPPRO pour
C     le post traitement des proprietes physiques definies
C     aux cellules afin de ne pas en oublier.
C     Les pointeurs ont ete initialises a 1 (poubelle).
C       dans iniini (avant usini1).
C     IPPPST commence a 2 car 1 est une poubelle.
C  Attention, IPPPST ressert plus bas.
C
C
C   - Repere pour la position des grandeurs posttraitees.
C
        IPPPST = 1
C
C   - Pointeurs post pour les variables.
C
        DO IVAR = 1, NVAR
          IPPPST        = IPPPST + 1
          IPPRTP(IVAR)  = IPPPST
        ENDDO
C
C
C   - Numero des proprietes physiques definies au centre de cellules.
C     Et pointeur pour le post.
C      (attention, la viscosite turbulente n'est pas post traitable
C       pour les calculs laminaires, car elle n'existe pas)
C      (attention, en module electrique, on utilise la meme
C       conductivite electrique pour le potentiel reel et le potentiel
C       imaginaire
C
        IPROP = 0
        DO IPHAS = 1, NPHAS
C
          IPROP                 = IPROP  + 1
          IPPROC(IROM  (IPHAS)) = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
          IPROP                 = IPROP  + 1
          IPPROC(IVISCL(IPHAS)) = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
          IPROP                 = IPROP  + 1
          IPPROC(IVISCT(IPHAS)) = IPROP
          IF (ITURB(IPHAS).EQ.0) THEN
            IPPPRO(IPROP)         = 1
          ELSE
            IPPPST                = IPPPST + 1
            IPPPRO(IPROP)         = IPPPST
          ENDIF
          IPROP                 = IPROP  + 1
          IPPROC(ICOUR (IPHAS)) = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
          IPROP                 = IPROP  + 1
          IPPROC(IFOUR (IPHAS)) = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
          IPROP                 = IPROP  + 1
          IPPROC(IPRTOT(IPHAS)) = IPROP
          IPPPST                = IPPPST + 1
          IPPPRO(IPROP)         = IPPPST
C
          IF(ICP(IPHAS).GT.0) THEN
            IPROP                 = IPROP + 1
            IPPROC(ICP   (IPHAS)) = IPROP
            IPPPST                = IPPPST + 1
            IPPPRO(IPROP)         = IPPPST
          ENDIF
C
          IF(ISMAGO(IPHAS).NE.-1) THEN
            IPROP                 = IPROP  + 1
            IPPROC(ISMAGO(IPHAS)) = IPROP
            IPPPST                = IPPPST + 1
            IPPPRO(IPROP)         = IPPPST
          ENDIF
C
          IF (IALE.EQ.1) THEN
            IPROP             = IPROP + 1
            IPPROC(IVISMA(1)) = IPROP
            IPPPST            = IPPPST + 1
            IPPPRO(IPROP)     = IPPPST
            IF (IORTVM.EQ.1) THEN
              IPROP             = IPROP + 1
              IPPROC(IVISMA(2)) = IPROP
              IPPPST            = IPPPST + 1
              IPPPRO(IPROP)     = IPPPST
              IPROP             = IPROP + 1
              IPPROC(IVISMA(3)) = IPROP
              IPPPST            = IPPPST + 1
              IPPPRO(IPROP)     = IPPPST
            ENDIF
          ENDIF
C
          DO IEST = 1, NESTMX
            IF(IESCAL(IEST,IPHAS).GT.0) THEN
              IPROP                      = IPROP + 1
              IPPROC(IESTIM(IEST,IPHAS)) = IPROP
              IPPPST                     = IPPPST + 1
              IPPPRO(IPROP)              = IPPPST
            ENDIF
          ENDDO

        ENDDO
C
C     Conductivite electrique imaginaire :
C     La conductivite reelle et imaginaire sont dans le meme tableau.
C       En Joule, on ne reserve donc pas de propriete "viscosite"
C       pour le potentiel imaginaire.
C       Intervention 1/2
C
        IF(IPPMOD(IELJOU).EQ.2 .OR. IPPMOD(IELJOU).EQ.4) THEN
          DO II = 1, NSCAL
            IF(IVISLS(II).GT.0.AND.II.NE.IPOTI) THEN
              IF(ISCAVR(II).LE.0) THEN
                IPROP                 = IPROP + 1
                IPPROC(IVISLS(II)  )  = IPROP
                IPPPST                = IPPPST + 1
                IPPPRO(IPROP)         = IPPPST
              ENDIF
            ENDIF
          ENDDO
        ELSE
C
          DO II = 1, NSCAL
            IF(IVISLS(II).GT.0) THEN
              IF(ISCAVR(II).LE.0) THEN
                IPROP                 = IPROP + 1
                IPPROC(IVISLS(II)  )  = IPROP
                IPPPST                = IPPPST + 1
                IPPPRO(IPROP)         = IPPPST
              ENDIF
            ENDIF
          ENDDO
C
       ENDIF
C
        DO II = 1, NSCAL
          IF(IVISLS(II).GT.0) THEN
            IF(ISCAVR(II).GT.0.AND.ISCAVR(II).LE.NSCAL) THEN
              IPPROC(IVISLS(II)  )  = IPPROC(IVISLS(ISCAVR(II))  )
            ENDIF
          ENDIF
        ENDDO
C
C     Conductivite electrique imaginaire :
C     La conductivite reelle et imaginaire sont dans le meme tableau.
C       En Joule, le pointeur sur la  "viscosite" du potentiel imaginaire
C         renvoie sur celle du potentiel reel.
C       Intervention 2/2
C
        IF(IPPMOD(IELJOU).EQ.2 .OR. IPPMOD(IELJOU).EQ.4) THEN
          IPPROC(IVISLS(IPOTI)  )  = IPPROC(IVISLS(IPOTR)  )
        ENDIF
C
        NPROCE = IPROP
C
C   Au centre des faces de bord (rho et flux de masse)
C
        IPROP = 0
        DO IPHAS = 1, NPHAS
          IPROP                 = IPROP + 1
          IPPROB(IROM  (IPHAS)) = IPROP
        ENDDO
        DO IFLUM = 1, NFLUMA
          IPROP                 = IPROP + 1
          IPPROB(IPROFL+IFLUM)  = IPROP
        ENDDO
        NPROFB = IPROP
C
C   Au centre des faces internes (flux de masse)
C
        IPROP = 0
        DO IFLUM = 1, NFLUMA
          IPROP                 = IPROP + 1
          IPPROF(IPROFL+IFLUM)  = IPROP
        ENDDO
        NPROFA = IPROP
C
C
C --- Modifications pour la physique particuliere
C      des entiers NPROCE, NPROFA, NPROFB
C
C      Sauvegarde pour la physique particuliere de IPROP
C      afin d'initialiser les positions des variables d'etat
C      Attention IPROPP est le dernier numero affecte pour les proprietes.
C                IPPPST est le rang de la derniere grandeur definie aux
C                       cellules pour le post traitement
        IPROPP = NPRMAX
        CALL PPPROP (IPROPP,IPPPST)
C       ===========
C
C     On sauve IPROPP modifie pour les passages ulterieurs dans varpos
        NPRMAX = IPROPP
C     On sauve IPPPST modifie pour les passages ulterieurs dans varpos
        NPPMAX = IPPPST
C
C
C --- Verification de NPROCE, NPROFA, NPROFB
C
        IF(NPROCE.GT.NPROMX.OR.
     &     NPROFA.GT.NPROMX.OR.NPROFB.GT.NPROMX) THEN
          WRITE(NFECRA,7200)NPROCE, NPROFA, NPROFB, NPROMX,
     &         MAX(MAX(NPROCE,NPROFA),NPROFB)
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
        RETURN
C
      ENDIF
C
C=======================================================================
C 4. QUATRIEME APPEL :
C        POSITIONNEMENT DES PROPRIETES POUR LE SCHEMA EN TEMPS,
C                                           LES MOMENTS ET FIN
C=======================================================================
C
C     Dans le cas ou on a un schema en temps d'ordre 2, il faut aussi
C       prevoir les proprietes au temps n-1. Ce sera fait au dernier appel
C
C     Dans le cas ou on calcule des moments, il faut en prevoir le nombre
C       et prevoir le nombre de tableaux necessaires pour stocker le
C       temps cumule de moyenne. On suppose que l'on peut s'aider
C       d'infos situees en tete de fichier suite (si on fait une
C       suite avec des moments non reinitialises).
C
      IF(IPASS.EQ.4) THEN
C
C
C ---> 4.1 PROPRIETES ADDITIONNELLES POUR LES ET SCHEMA EN TEMPS
C      ---------------------------------------------------------
C
C --- Initialisations par defaut eventuelles et verifications
C       des options utilisees ci-dessous pour decider si l'on
C       reserve des tableaux supplementaires pour des grandeurs
C       au pas de temps precedent
C
        IOK = 0
C
C     Pression hydrostatique
        IF(IPHYDR.EQ.0) THEN
          ICALHY = 0
        ELSEIF(IPHYDR.EQ.1) THEN
          GRAVN2 = GX**2+GY**2+GZ**2
          IF (GRAVN2.LT.EPZERO**2) THEN
            ICALHY = 0
          ELSE
            ICALHY = 1
          ENDIF
        ENDIF
C
C     Schemas en temps
C         en LES : Ordre 2 ; sinon Ordre 1
C         (en particulier, ordre 2 impossible en k-eps couple)
        DO IPHAS = 1, NPHAS
          IF(ISCHTP(IPHAS).EQ.-999) THEN
            IF(ITYTUR(IPHAS).EQ.4) THEN
              ISCHTP(IPHAS) = 2
            ELSE
              ISCHTP(IPHAS) = 1
            ENDIF
          ENDIF
        ENDDO
C
C     Schemas en temps : variables deduites
        DO IPHAS = 1, NPHAS
C     Schema pour le Flux de masse
          IF(ISTMPF(IPHAS).EQ.-999) THEN
            IF(ISCHTP(IPHAS).EQ.1) THEN
              ISTMPF(IPHAS) = 1
            ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
              ISTMPF(IPHAS) = 2
            ENDIF
          ENDIF
C     Masse volumique
          IF(IROEXT(IPHAS).EQ.-999) THEN
            IF(ISCHTP(IPHAS).EQ.1) THEN
              IROEXT(IPHAS) = 0
            ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
C       Pour le moment par defaut on ne prend pas l'ordre 2
C              IROEXT(IPHAS) = 1
              IROEXT(IPHAS) = 0
            ENDIF
          ENDIF
C     Viscosite
          IF(IVIEXT(IPHAS).EQ.-999) THEN
            IF(ISCHTP(IPHAS).EQ.1) THEN
              IVIEXT(IPHAS) = 0
            ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
C       Pour le moment par defaut on ne prend pas l'ordre 2
C              IVIEXT(IPHAS) = 1
              IVIEXT(IPHAS) = 0
            ENDIF
          ENDIF
C     Chaleur massique
          IF(ICPEXT(IPHAS).EQ.-999) THEN
            IF(ISCHTP(IPHAS).EQ.1) THEN
              ICPEXT(IPHAS) = 0
            ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
C       Pour le moment par defaut on ne prend pas l'ordre 2
C              ICPEXT(IPHAS) = 1
              ICPEXT(IPHAS) = 0
            ENDIF
          ENDIF
C     Termes sources NS,
          IF(ISNO2T(IPHAS).EQ.-999) THEN
            IF(ISCHTP(IPHAS).EQ.1) THEN
              ISNO2T(IPHAS) = 0
C            ELSEIF(ISCHTP(IPHAS).EQ.2.AND.IVISSE(IPHAS).EQ.1) THEN
            ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
C       Pour le moment par defaut on prend l'ordre 2
              ISNO2T(IPHAS) = 1
C              ISNO2T(IPHAS) = 0
            ENDIF
          ENDIF
C     Termes sources turbulence (k-eps, Rij, v2f ou k-omega)
C     On n'autorise de changer ISTO2T qu'en Rij (sinon avec
C       le couplage k-eps/omega il y a pb)
          IF(ISTO2T(IPHAS).EQ.-999) THEN
            IF(ISCHTP(IPHAS).EQ.1) THEN
              ISTO2T(IPHAS) = 0
            ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
C       Pour le moment par defaut on ne prend pas l'ordre 2
C              ISTO2T(IPHAS) = 1
              ISTO2T(IPHAS) = 0
            ENDIF
          ELSE IF( ITYTUR(IPHAS).EQ.2.OR.ITURB(IPHAS).EQ.50
     &         .OR.ITURB(IPHAS).NE.60) THEN
            WRITE(NFECRA,8132) IPHAS,ITURB(IPHAS),ISTO2T(IPHAS)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
        DO ISCAL = 1, NSCAL
C     Termes sources Scalaires,
          IPHAS = IPHSCA(ISCAL)
          IF(ISSO2T(ISCAL).EQ.-999) THEN
            IF(ISCHTP(IPHAS).EQ.1) THEN
              ISSO2T(ISCAL) = 0
            ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
C       Pour coherence avec Navier Stokes on prend l'ordre 2
C       mais de toute facon qui dit ordre 2 dit LES et donc
C       generalement pas de TS scalaire a interpoler.
              ISSO2T(ISCAL) = 1
C              ISSO2T(ISCAL) = 0
            ENDIF
          ENDIF
C     Diffusivite scalaires
          IF(IVSEXT(ISCAL).EQ.-999) THEN
            IPHAS = IPHSCA(ISCAL)
            IF(ISCHTP(IPHAS).EQ.1) THEN
              IVSEXT(ISCAL) = 0
            ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
C       Pour le moment par defaut on ne prend pas l'ordre 2
C              IVSEXT(ISCAL) = 1
              IVSEXT(ISCAL) = 0
            ENDIF
          ENDIF
        ENDDO
C
C
C     Pression hydrostatique
        IF (IPHYDR.NE.0.AND.IPHYDR.NE.1) THEN
          WRITE(NFECRA,8021) 'IPHYDR ',IPHYDR
          IOK = IOK + 1
        ENDIF
C
C     Viscosite secondaire
        DO IPHAS = 1, NPHAS
          IVISPH = IVISSE(IPHAS)
          IF (IVISPH.NE.0.AND.IVISPH.NE.1) THEN
            WRITE(NFECRA,8022) IPHAS,'IVISSE ',IVISPH
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C     Schemas en temps
        DO IPHAS = 1, NPHAS
C
C     Schema en temps global.
          IF(ISCHTP(IPHAS).NE. 1.AND.ISCHTP(IPHAS).NE.2) THEN
            WRITE(NFECRA,8101) IPHAS,'ISCHTP',ISCHTP(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF(ISCHTP(IPHAS).EQ. 2.AND.IDTVAR.NE.0) THEN
            WRITE(NFECRA,8111) IPHAS,ISCHTP(IPHAS),IDTVAR
            IOK = IOK + 1
          ENDIF
          IF(ISCHTP(IPHAS).EQ. 2.AND.ITYTUR(IPHAS).EQ.2) THEN
            WRITE(NFECRA,8112) IPHAS,ISCHTP(IPHAS),ITURB(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF(ISCHTP(IPHAS).EQ.1.AND.ITYTUR(IPHAS).EQ.4) THEN
            WRITE(NFECRA,8113) IPHAS,ISCHTP(IPHAS),ITURB(IPHAS)
          ENDIF
          IF(ISCHTP(IPHAS).EQ. 2.AND.ITURB(IPHAS).EQ.50) THEN
            WRITE(NFECRA,8114) IPHAS,ISCHTP(IPHAS),ITURB(IPHAS)
            IOK = IOK + 1
          ENDIF
          IF(ISCHTP(IPHAS).EQ. 2.AND.ITURB(IPHAS).EQ.60) THEN
            WRITE(NFECRA,8115) IPHAS,ISCHTP(IPHAS),ITURB(IPHAS)
            IOK = IOK + 1
          ENDIF
C
C     Schema en temps pour le flux de masse
          IF(ISTMPF(IPHAS).NE. 2.AND.ISTMPF(IPHAS).NE.0.AND.
     &       ISTMPF(IPHAS).NE. 1) THEN
            WRITE(NFECRA,8121) IPHAS,'ISTMPF',ISTMPF(IPHAS)
            IOK = IOK + 1
          ENDIF
C
C     Schema en temps pour les termes sources de NS
          IF(ISNO2T(IPHAS).NE.0.AND.
     &       ISNO2T(IPHAS).NE. 1.AND.ISNO2T(IPHAS).NE.2) THEN
            WRITE(NFECRA,8131) IPHAS,'ISNO2T',ISNO2T(IPHAS)
            IOK = IOK + 1
          ENDIF
C     Schema en temps pour les termes sources des grandeurs
C     turbulentes
          IF(ISTO2T(IPHAS).NE.0.AND.
     &       ISTO2T(IPHAS).NE. 1.AND.ISTO2T(IPHAS).NE.2) THEN
            WRITE(NFECRA,8131) IPHAS,'ISTO2T',ISTO2T(IPHAS)
            IOK = IOK + 1
          ENDIF
C
C     Schema en temps pour la masse volumique
          IF(IROEXT(IPHAS).NE.0.AND.
     &       IROEXT(IPHAS).NE. 1.AND.IROEXT(IPHAS).NE.2) THEN
            WRITE(NFECRA,8131) IPHAS,'IROEXT',IROEXT(IPHAS)
            IOK = IOK + 1
          ENDIF
C
C     Schema en temps pour la viscosite
          IF(IVIEXT(IPHAS).NE.0.AND.
     &       IVIEXT(IPHAS).NE. 1.AND.IVIEXT(IPHAS).NE.2) THEN
            WRITE(NFECRA,8131) IPHAS,'IVIEXT',IVIEXT(IPHAS)
            IOK = IOK + 1
          ENDIF
C
C     Schema en temps pour la chaleur specifique
          IF(ICPEXT(IPHAS).NE.0.AND.
     &       ICPEXT(IPHAS).NE. 1.AND.ICPEXT(IPHAS).NE.2) THEN
            WRITE(NFECRA,8131) IPHAS,'ICPEXT',ICPEXT(IPHAS)
            IOK = IOK + 1
          ENDIF
C
        ENDDO
C
        DO ISCAL = 1, NSCAL
C     Schema en temps pour les termes sources des scalaires
          IF(ISSO2T(ISCAL).NE.0.AND.
     &       ISSO2T(ISCAL).NE. 1.AND.ISSO2T(ISCAL).NE.2) THEN
            WRITE(NFECRA,8141) ISCAL,'ISSO2T',ISSO2T(ISCAL)
            IOK = IOK + 1
          ENDIF
C     Schema en temps pour la viscosite
          IF(IVSEXT(ISCAL).NE.0.AND.
     &       IVSEXT(ISCAL).NE. 1.AND.IVSEXT(ISCAL).NE.2) THEN
            WRITE(NFECRA,8141) ISCAL,'IVSEXT',IVSEXT(ISCAL)
            IOK = IOK + 1
          ENDIF
        ENDDO
C
C     Stop si probleme
        IF(IOK.GT.0) THEN
          CALL CSEXIT(1)
        ENDIF
C
C
C --- Reprise du dernier numero de propriete
        IPROP  = NPRMAX
C
C --- Numeros de propriete
        DO IPHAS = 1, NPHAS
C     On a besoin de la masse volumique si on l'extrapole ou si ICALHY
          IF(IROEXT(IPHAS).GT.0.OR.ICALHY.EQ.1) THEN
            IPROP         = IPROP + 1
            IROMA (IPHAS) = IPROP
          ENDIF
C     Dans le cas d'une extrapolation de la viscosite totale
          IF(IVIEXT(IPHAS).GT.0) THEN
            IPROP         = IPROP + 1
            IVISLA(IPHAS) = IPROP
            IPROP         = IPROP + 1
            IVISTA(IPHAS) = IPROP
          ENDIF
        ENDDO
C     Proprietes des phases : CP s'il est variable
        DO IPHAS = 1, NPHAS
          IF(ICP(IPHAS).NE.0) THEN
            IF(ICPEXT(IPHAS).GT.0) THEN
              IPROP         = IPROP + 1
              ICPA  (IPHAS) = IPROP
            ENDIF
          ENDIF
        ENDDO
C     On a besoin d'un tableau pour les termes sources de Navier Stokes
C       a extrapoler. Ce tableau est NDIM
        DO IPHAS = 1, NPHAS
          IF(ISNO2T(IPHAS).GT.0) THEN
            IPROP         = IPROP + 1
            ITSNSA(IPHAS) = IPROP
          ENDIF
          IF(ISTO2T(IPHAS).GT.0) THEN
            IPROP         = IPROP + 1
            ITSTUA(IPHAS) = IPROP
          ENDIF
        ENDDO
C     Proprietes des scalaires : termes sources pour theta schema
C       et VISCLS si elle est variable
        IF(NSCAL.GE.1) THEN
          DO ISCAL = 1, NSCAL
            IF(ISSO2T(ISCAL).GT.0) THEN
              IPROP         = IPROP + 1
              ITSSCA(ISCAL) = IPROP
            ENDIF
            IF(IVISLS(ISCAL).NE.0) THEN
              IF(IVSEXT(ISCAL).GT.0) THEN
                IPROP         = IPROP + 1
                IVISSA(ISCAL) = IPROP
              ENDIF
            ENDIF
          ENDDO
        ENDIF
C     Proprietes des variables : flux de masse porteur
C       On en ajoute un (et un seul) par phase s'il existe une phase
C         qui en a besoin.
C       On est donc dans l'hypothese implicite qu'il n'y a qu'un seul
C         flux de masse pour toutes les variables
C         et ceci n'est pas optimal si il y a plusieurs phases traitees
C         differemment.
C       On suppose que si une phase en a besoin, on en ajoute donc
C         autant qu'il y a de flux de masse, cad NFLUMA
C
C     On les initialise a -1 pour iniva0
        DO IVAR = 1, NVARMX
          IFLUAA(IVAR) = -1
        ENDDO
C     On regarde s'il y en a besoin
        IIFLAA = 0
        DO IPHAS = 1, NPHAS
          IF(ISTMPF(IPHAS).NE.1) IIFLAA = 1
        ENDDO
C     On les affecte
        IPROFA = IPROP
        IF(IIFLAA.EQ.1) THEN
          DO IPHAS = 1, NPHAS
            IF(IPHAS.EQ.1) THEN
              IPROP               = IPROP + 1
              IFLUAA(IPR (IPHAS)) = IPROP
            ENDIF
            IFLUAA(IU  (IPHAS)) = IPROP
            IFLUAA(IV  (IPHAS)) = IPROP
            IFLUAA(IW  (IPHAS)) = IPROP
            IF(ITYTUR(IPHAS).EQ.2) THEN
              IFLUAA(IK  (IPHAS)) = IPROP
              IFLUAA(IEP (IPHAS)) = IPROP
            ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
              IFLUAA(IR11(IPHAS)) = IPROP
              IFLUAA(IR22(IPHAS)) = IPROP
              IFLUAA(IR33(IPHAS)) = IPROP
              IFLUAA(IR12(IPHAS)) = IPROP
              IFLUAA(IR13(IPHAS)) = IPROP
              IFLUAA(IR23(IPHAS)) = IPROP
              IFLUAA(IEP (IPHAS)) = IPROP
            ELSEIF(ITURB(IPHAS).EQ.50) THEN
              IFLUAA(IK  (IPHAS)) = IPROP
              IFLUAA(IEP (IPHAS)) = IPROP
              IFLUAA(IPHI(IPHAS)) = IPROP
              IFLUAA(IFB (IPHAS)) = IPROP
            ELSEIF(ITURB(IPHAS).EQ.60) THEN
              IFLUAA(IK  (IPHAS)) = IPROP
              IFLUAA(IOMG(IPHAS)) = IPROP
            ENDIF
          ENDDO
          DO ISCAL = 1, NSCAL
            IFLUAA(ISCA(ISCAL)) = IFLUAA(IU(IPHSCA(ISCAL)))
          ENDDO
        ENDIF
C
C --- Sauvegarde du dernier numero de propriete
        NPRMAX = IPROP
C
C
C --- Reprise du dernier NPROCE et du dernier NPPMAX
        IPROP                 = NPROCE
        IPPPST                = NPPMAX
C
C --- Positionnement des PROPCE
        DO IPHAS = 1, NPHAS
C
C     Variables schema en temps
          IF(IROEXT(IPHAS).GT.0.OR.ICALHY.EQ.1) THEN
            IPROP                 = IPROP  + 1
            IPPROC(IROMA (IPHAS)) = IPROP
          ENDIF
          IF(IVIEXT(IPHAS).GT.0) THEN
            IPROP                 = IPROP  + 1
            IPPROC(IVISLA(IPHAS)) = IPROP
          ENDIF
          IF(IVIEXT(IPHAS).GT.0) THEN
            IPROP                 = IPROP  + 1
            IPPROC(IVISTA(IPHAS)) = IPROP
          ENDIF
          IF(ICPEXT(IPHAS).GT.0) THEN
            IPROP                 = IPROP + 1
            IPPROC(ICPA  (IPHAS)) = IPROP
          ENDIF
          IF(ISNO2T(IPHAS).GT.0) THEN
            IPROP                 = IPROP + 1
            IPPROC(ITSNSA(IPHAS)) = IPROP
C     Ce tableau est NDIM :
            IPROP                 = IPROP + NDIM-1
          ENDIF
          IF(ISTO2T(IPHAS).GT.0) THEN
            IPROP                 = IPROP + 1
            IPPROC(ITSTUA(IPHAS)) = IPROP
C     Ce tableau est 2, 7 ou 4 selon le modele de turbulence :
            IF    (ITYTUR(IPHAS).EQ.2) THEN
              IPROP                 = IPROP + 2-1
            ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
              IPROP                 = IPROP + 7-1
            ELSEIF(ITURB(IPHAS).EQ.50) THEN
              IPROP                 = IPROP + 4-1
            ENDIF
          ENDIF
        ENDDO
C
        DO II = 1, NSCAL
C Termes source des scalaires pour theta schema
          IF(ISSO2T(II).GT.0) THEN
            IPROP                 = IPROP + 1
            IPPROC(ITSSCA(II))    = IPROP
          ENDIF
C
          IF(IVISLS(II).GT.0) THEN
            IF(ISCAVR(II).LE.0) THEN
              IF(IVSEXT(II).GT.0) THEN
                IPROP                 = IPROP + 1
                IPPROC(IVISSA(II)  )  = IPROP
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        DO II = 1, NSCAL
          IF(IVISLS(II).GT.0) THEN
            IF(ISCAVR(II).GT.0.AND.ISCAVR(II).LE.NSCAL) THEN
              IF(IVSEXT(II).GT.0) THEN
                IPPROC(IVISSA(II)  )  = IPPROC(IVISSA(ISCAVR(II))  )
              ENDIF
            ENDIF
          ENDIF
        ENDDO
C
C --- Sauvegarde du dernier NPROCE et du dernier NPPMAX
        NPROCE = IPROP
        NPPMAX = IPPPST
C
C
C --- Reprise du dernier NPROFB
        IPROP                 = NPROFB
C
C --- Positionnement des PROPFB
        DO IPHAS = 1, NPHAS
C     Variables schema en temps : rhoa (pas pour icalhy)
          IF(IROEXT(IPHAS).GT.0) THEN
            IPROP                 = IPROP  + 1
            IPPROB(IROMA (IPHAS)) = IPROP
          ENDIF
        ENDDO
C     Variables schema en temps : flux de masse A
        IF(IIFLAA.EQ.1) THEN
          DO IFLUM = 1, NFLUMA
            IPROP                 = IPROP + 1
            IPPROB(IPROFA+IFLUM)  = IPROP
          ENDDO
        ENDIF
C
C --- Sauvegarde du dernier NPROFB
        NPROFB = IPROP
C
C
C --- Reprise du dernier NPROFA
        IPROP                 = NPROFA
C
C --- Positionnement des PROPFA
        IF(IIFLAA.EQ.1) THEN
          DO IFLUM = 1, NFLUMA
            IPROP                 = IPROP + 1
            IPPROF(IPROFA+IFLUM)  = IPROP
          ENDDO
        ENDIF
C
C --- Sauvegarde du dernier NPROFA
        NPROFA = IPROP
C
C
C ---> 4.2 CALCUL DE LA TAILLE DU TABLEAU DES TEMPS CUMULES POUR LES MOMENTS
C      ---------------------------------------------------------------------
C
C     Pour verification des definitions de moments
        IOK = 0
C
C --- Calcul du nombre de moments definis (et verif qu'il n y a pas de trous)
        NBMOMT = 0
        INMFIN = 0
        DO IMOM = 1, NBMOMX
C     Si on n'est pas a la fin de la liste
          IF(INMFIN.EQ.0) THEN
C       Si il y en a, ca en fait en plus
            IF(IDFMOM(1,IMOM).NE.0) THEN
              NBMOMT = NBMOMT + 1
C       Si il n'y en a pas, c'est la fin de la liste
            ELSE
              INMFIN = 1
            ENDIF
C     Si on est a la fin de la liste, il n'en faut plus
          ELSE
            IF(IDFMOM(1,IMOM).NE.0) THEN
              IOK = IOK + 1
            ENDIF
          ENDIF
        ENDDO
C
        IF(IOK.NE.0) THEN
          WRITE(NFECRA,8200)NBMOMT+1,IDFMOM(1,NBMOMT+1),NBMOMT,NBMOMT
          DO IMOM = 1, NBMOMX
            WRITE(NFECRA,8201)IMOM,IDFMOM(1,IMOM)
          ENDDO
          WRITE(NFECRA,8202)
        ENDIF
C
C --- Verification de IDFMOM
        IOK = 0
        DO IMOM = 1, NBMOMX
          IDFFIN = 0
          DO JJ = 1, NDGMOX
            IDFMJI = IDFMOM(JJ,IMOM)
            IF(IDFFIN.EQ.0) THEN
              IF(IDFMJI.LT.-NPRMAX) THEN
                IOK = IOK + 1
                WRITE(NFECRA,8210)JJ,IMOM,IDFMJI,NPRMAX
              ELSEIF(IDFMJI.GT.NVAR) THEN
                IOK = IOK + 1
                WRITE(NFECRA,8211)JJ,IMOM,IDFMJI,NVAR
              ELSEIF(IDFMJI.LT.0) THEN
                IF (IPPROC(-IDFMJI).LE.0) THEN
                  IOK = IOK + 1
                  WRITE(NFECRA,8212)JJ,IMOM,IDFMJI,-IDFMJI,
     &                 IPPROC(-IDFMJI)
                ENDIF
              ELSEIF(IDFMJI.EQ.0) THEN
                IDFFIN = 1
              ENDIF
            ELSE
              IF(IDFMJI.NE.0) THEN
                IOK = IOK + 1
                WRITE(NFECRA,8213)IMOM,JJ,IDFMJI
              ENDIF
            ENDIF
          ENDDO
        ENDDO
C
C --- Verification de NTDMOM (>0)
        DO IMOM = 1, NBMOMT
          IF(NTDMOM(IMOM).LT.0) THEN
            IOK = IOK + 1
            WRITE(NFECRA,8214)IMOM,NTDMOM(IMOM)
          ENDIF
        ENDDO
C
        IF(IOK.NE.0) THEN
          CALL CSEXIT(1)
        ENDIF
C
C
C --- Indicateur pour juger de l'utilite de consulter le fichier suite
C     Il n'est pas indispensable de consulter avec succes le fichier suite
C       lorsqu'il n'y a aucun moment a relire ou qu'on ne fait pas de suite
C       (si ILSMOM=0)
C       ie ILSMOM = 0 == (ISUITE=0 OR IMOOLD(IMOM)=-1 pour tout IMOM)
C       ie ILSMOM = 1 == (ISUITE=1 AND IMOOLD(IMOM) different de -1 pour un IMOM)
        ILSMOM = 0
        IF(ISUITE.EQ.1.AND.NBMOMT.GT.0) THEN
          DO IMOM = 1, NBMOMT
            IF(IMOOLD(IMOM).NE.-1) THEN
              ILSMOM = 1
            ENDIF
          ENDDO
        ENDIF
C
C --- Lecture du fichier suite (debut : info sur les moments et sur le ntpabs)
        IF(ILSMOM.EQ.1) THEN
C
C     Ouverture
C        (ILECEC=1:lecture)
          ILECEC = 1
          CALL OPNSUI(FICAMX,LEN(FICAMX),ILECEC,IMPAMX,IERROR)
C         ===========
          IF (IERROR.NE.0) THEN
            WRITE(NFECRA,8300) FICAMX
            CALL CSEXIT (1)
          ENDIF
C
C     Tests : type de fichier et dimension des supports (puisqu'on y
C       relira des infos, autant s'arreter au plus vite
C       s'il n'est pas correct)
C
C     Type (fichier auxiliaire)
C
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'version_fichier_suite_auxiliaire'
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                IVERS,IERROR)
C
          IF (IERROR.NE.0) THEN
            WRITE(NFECRA,8301)FICAMX
            CALL CSEXIT (1)
          ENDIF
C
C     Supports (nombre de cellules)
C       Ce n'est pas possible de tester ceci en parallele ici
C       car on n'a pas encore appele parge1
C
          IF(IRANGP.LT.0) THEN
            CALL TSTSUI(IMPAMX,NCELOK,NFAIOK,NFABOK,NSOMOK)
C           ===========
            IF (NCELOK.EQ.0) THEN
              WRITE(NFECRA,8302)
              CALL CSEXIT (1)
            ENDIF
          ENDIF
C
C
          NBERRO=0
C
C     Nb de moments
C       Si le nombre de moments n'a pas pu etre relu, on devra s'arreter
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 1
          RUBRIQ = 'nombre_moyennes_temps'
          JBMOMT = 0
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                JBMOMT,IERROR)
          NBERRO=NBERRO+IERROR
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,8311)
          ENDIF
C
C     Nombre de pas de temps final du calcul precedent
          RUBRIQ = 'nbre_pas_de_temps'
          ITYSUP = 0
          NBVAL  = 1
          IRTYP  = 1
          JTCABS = 0
          CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,IRTYP,
     &                JTCABS,IERROR)
          NBERRO=NBERRO+IERROR
          IF(IERROR.NE.0) THEN
            WRITE(NFECRA,8312)
          ENDIF
C
          IF(NBERRO.NE.0) THEN
            CALL CSEXIT (1)
          ENDIF
C
        ENDIF
C
C
C --- Avec JTCABS et JBMOMT
C     on complete et on verifie la correspondance nouveaux -> anciens moments
C
C       Si suite et il existe un IMOOLD different de -1 :
C                      si IMOOLD mal renseigne erreur
C                         (valeur non admissible ou demande de relecture
C                             avec NTDBMO pas coherent)
C                      si IMOOLD pas renseigne pour NTDMOM > JTCABS : initialise
C                                              pour IMOM < JBMOMT+1 : on relit
C                                              pour IMOM > JBMOMT : initialise
C                      si on pointe deux fois sur le meme : erreur
C                                  (simplifie la lecture lecamo)
        IOK = 0
C
        IF(ILSMOM.EQ.1) THEN
          DO IMOM = 1, NBMOMT
            IF(IMOOLD(IMOM).GT.JBMOMT.OR.IMOOLD(IMOM).EQ.0
     &                             .OR.IMOOLD(IMOM).LT.-2) THEN
              WRITE(NFECRA,8400) IMOM,IMOOLD(IMOM),JBMOMT
              IOK = IOK + 1
            ELSEIF(IMOOLD(IMOM).LE.JBMOMT.AND.IMOOLD(IMOM).GT.0
     &                             .AND.NTDMOM(IMOM).GT.JTCABS) THEN
              WRITE(NFECRA,8401) IMOM,IMOOLD(IMOM),NTDMOM(IMOM),JTCABS
              IOK = IOK + 1
            ELSEIF(IMOOLD(IMOM).EQ.-2.AND.NTDMOM(IMOM).GT.JTCABS) THEN
              IMOOLD(IMOM)=-1
            ELSEIF(IMOOLD(IMOM).EQ.-2.AND.IMOM.LE.JBMOMT) THEN
              IMOOLD(IMOM)=IMOM
            ELSEIF(IMOOLD(IMOM).EQ.-2.AND.IMOM.GT.JBMOMT) THEN
              IMOOLD(IMOM)=-1
            ENDIF
          ENDDO
          DO IMOM = 1, NBMOMT
            IF(IMOOLD(IMOM).GT.0) THEN
              DO JMOM = 1, IMOM-1
                IF(IMOOLD(IMOM).EQ.IMOOLD(JMOM)) THEN
                  WRITE(NFECRA,8402) IMOM,JMOM,IMOOLD(IMOM)
                  IOK = IOK + 1
                ENDIF
              ENDDO
            ENDIF
          ENDDO
C
C       Si pas suite : si IMOOLD pas renseigne, on initialise
C                      sinon erreur
        ELSEIF(ISUITE.EQ.0) THEN
          DO IMOM = 1, NBMOMT
            IF(IMOOLD(IMOM).NE.-2) THEN
              WRITE(NFECRA,8403) ISUITE,IMOM,IMOOLD(IMOM)
              IOK = IOK + 1
            ELSE
              IMOOLD(IMOM)=-1
            ENDIF
          ENDDO
        ENDIF
C
        IF(IOK.NE.0) THEN
          CALL CSEXIT (1)
        ENDIF
C
C --- D'autres informations sont lues dans le fichier suite pour
C     l'initialisation du numero dans le fichier suite du cumul temporel
C       associe a chaque moyenne (IDTOLD)
C     IDTOLD(IMOM) donne pour la moyenne (du calcul courant) IMOM
C       le numero dans le fichier suite du cumul temporel associe
C       a la moyenne du calcul precedent a laquelle correspond IMOM (ouf!)
C
C     Initialisation
        DO IMOM = 1, NBMOMX
          IDTOLD(IMOM) = 0
        ENDDO
C
C     Si on ne lit rien dans le fichier suite, IDTOLD restera nul
        IF(ILSMOM.EQ.1) THEN
C
          NBERRO = 0
C
C     On ne lit des choses que s'il y a des moments stockes
          IF ( JBMOMT.GT.0 ) THEN
C
C        On lit les infos pour tous les moments courants
            DO IMOM = 1, NBMOMT
C           On s'interesse uniquement aux moments qui suivnt d'anciens moments
              IF(IMOOLD(IMOM).GT.0) THEN
C            Si le numero de l'ancien moment est incompatible avec le format
C              on cherche numero_cumul_temps_momentYYYY qui n existe pas
C              ca genere une erreur
                IF(IMOM.LE.NFMTMO) THEN
                  WRITE(CMOY4,'(I4.4)')IMOOLD(IMOM)
                ELSE
                  CMOY4 = CINDFM
                ENDIF
                ITYSUP = 0
                NBVAL  = 1
                IRTYP  = 1
                RUBRIQ = 'numero_cumul_temps_moment'//CMOY4
                CALL LECSUI(IMPAMX,RUBRIQ,LEN(RUBRIQ),ITYSUP,NBVAL,
     &                      IRTYP,IDTOLD(IMOM),IERROR)
                NBERRO=NBERRO+IERROR
                IF(IDTOLD(IMOM).EQ.0.OR.IERROR.NE.0) THEN
                  WRITE(NFECRA,8313)IMOM,IMOOLD(IMOM)
                  IOK = IOK + 1
                ENDIF
              ENDIF
            ENDDO
          ENDIF
C
          IF(NBERRO.NE.0) THEN
            CALL CSEXIT (1)
          ENDIF
C
C     Fermeture du fichier suite auxiliaire
          CALL CLSSUI(IMPAMX,IERROR)
C
          IF (IERROR.NE.0) THEN
            WRITE(NFECRA,8390) FICAMX
          ENDIF
C
        ENDIF
C

C --- Remplissage de IDTMOM (pointeur dans DTCMOM pour chaque moment)
C       > 0 : pointe sur un tableau NCEL dans PROPCE (DT cumule non uniforme)
C             remplissage sans trous de 1 a n
C       < 0 : pointe sur une case d'un tableau NPROMX  (DT cumule uniforme)
C             remplissage sans trous de -1 a -p
C
        IIPLUS = 0
        IIMOIN = 0
C
        DO IMOM = 1, NBMOMX
          IDTMOM(IMOM) = 0
        ENDDO
C
C       Pour les moments reinitialises (suite) ou initialises (non suite)
C         un nouveau tableau, sauf si leur calcul commence au meme instant
        DO IMOM = 1, NBMOMT
          IMOLD = IMOOLD(IMOM)
C        Si on (re)initialise IMOM
          IF(IMOLD.EQ.-1) THEN
C          On cherche si on en a deja vu (re)initialisees au meme moment
            IMOMR = 0
            DO JMOM = 1, IMOM-1
              JMOLD = IMOOLD(JMOM)
              IF(JMOLD.EQ.-1.AND.NTDMOM(JMOM).EQ.NTDMOM(IMOM)) THEN
                IMOMR = JMOM
              ENDIF
            ENDDO
C          Si oui : on utilise le meme tableau
            IF(IMOMR.GT.0) THEN
              IDTMOM(IMOM) = IDTMOM(IMOMR)
C          Si non : on a besoin d'un nouveau tableau ou reel
            ELSE
              IF(IDTVAR.EQ.2) THEN
                IIPLUS = IIPLUS + 1
                IDTMOM(IMOM) = IIPLUS
              ELSE
                IIMOIN = IIMOIN - 1
                IDTMOM(IMOM) = IIMOIN
              ENDIF
            ENDIF
          ENDIF
        ENDDO
C
C       Pour les moments IMOM relus dans IMOLD
C         (ie suite + IMOOLD non egal a -1)
        IF(ISUITE.EQ.1) THEN
          DO IMOM = 1, NBMOMT
            IMOLD = IMOOLD(IMOM)
            IF(IMOLD.GT.0) THEN
C         On regarde si le DTcumule du IMOLD a deja ete vu
C            (indicateur JMOMOK)
              IDTO = IDTOLD(IMOM)
              JMOMOK = 0
              DO JMOM = 1, IMOM-1
                JMOLD = IMOOLD(JMOM)
                IF(JMOLD.GT.0) THEN
                  JDTO = IDTOLD(JMOM)
                  IF(JDTO.EQ.IDTO) THEN
                    JMOMOK = JMOM
                  ENDIF
                ENDIF
              ENDDO
C         Si on le voit pour la premiere fois, ca en fait un de plus
              IF(JMOMOK.EQ.0) THEN
C           Si DT non uniforme dans le present calcul ou que
C              DT cumule etait non uniforme dans l'ancien, on a un
C              DT cumule non uniforme (sinon, il est uniforme)
                IF(IDTVAR.EQ.2.OR.IDTO.GT.0) THEN
                  IIPLUS = IIPLUS + 1
                  IDTMOM(IMOM) = IIPLUS
                ELSE
                  IIMOIN = IIMOIN - 1
                  IDTMOM(IMOM) = IIMOIN
                ENDIF
C         Si on l'a deja rencontre, on pointe au meme endroit
              ELSE
                IDTMOM(IMOM) = IDTMOM(JMOMOK)
              ENDIF
            ENDIF
          ENDDO
        ENDIF
C
C       Verification de IDTMOM : normalement, jamais ca plante ici.
        IOK = 0
        DO IMOM = 1, NBMOMT
          IF(IDTMOM(IMOM).EQ.0) THEN
            IOK = IOK + 1
            WRITE(NFECRA,8410)IMOM, IDTMOM(IMOM)
          ENDIF
        ENDDO
        IF(IOK.NE.0) THEN
          CALL CSEXIT (1)
        ENDIF
C
C --- Calcul du nombre de tableaux NCEL "temps cumule"
        NBDTCM = 0
        DO IMOM = 1, NBMOMT
          NBDTCM=MAX(IDTMOM(IMOM),NBDTCM)
        ENDDO
C
C ---> 4.3 POSITIONNEMENT DANS PROPCE DES MOMENTS ET DU TEMPS CUMULE
C      -------------------------------------------------------------
C
C --- Reprise du dernier numero de propriete
        IPROP  = NPRMAX
C
C --- Numeros de propriete
        DO IMOM = 1, NBMOMT
          IPROP  = IPROP + 1
          ICMOME(IMOM) = IPROP
        ENDDO
        DO II = 1, NBDTCM
          IPROP  = IPROP + 1
          ICDTMO(II) = IPROP
        ENDDO
C
C --- Sauvegarde du dernier numero de propriete
        NPRMAX = IPROP
C
C --- Reprise des derniers NPROCE et NPPMAX (PROPCE et POST-TRAITEMENT)
        IPROP                 = NPROCE
        IPPPST                = NPPMAX
C
C --- Positionnement
        DO IMOM = 1, NBMOMT
          IPROP                = IPROP + 1
          IPPROC(ICMOME(IMOM)) = IPROP
          IPPPST               = IPPPST + 1
          IPPPRO(IPROP)        = IPPPST
        ENDDO
        DO II = 1, NBDTCM
          IPROP                = IPROP + 1
          IPPROC(ICDTMO(II))   = IPROP
          IPPPST               = IPPPST + 1
          IPPPRO(IPROP)        = IPPPST
        ENDDO
C
C --- Sauvegarde du dernier NPROCE et NPPMAX
        NPROCE = IPROP
        NPPMAX = IPPPST
C
C
C
C ---> 4.4  POSITIONNEMENT DES CONDITIONS AUX LIMITES
C      ---------------------------------------------------------------------
C
C --- Numerotation des tableaux NFABOR de type COEFA/COEFB presents ici
C       On suppose pour le moment que seules les variables de calcul
C         disposent de conditions aux limites. Ceci pourra ensuite etre
C         etendu aux variables pysiques en declarant un pointeur du type
C         de ICLRTP (par exemple ICLPRO)
C       On suppose que les seules variables qui ont 2 types de cl sont
C         les variables vitesse en k-epsilon, k-omega et en LES, et la
C         pression si IPHYDR=1
C
        ICONDL = 0
        DO IVAR = 1, NVAR
          ICONDL = ICONDL + 1
          ICLRTP(IVAR,ICOEF ) = ICONDL
          ICLRTP(IVAR,ICOEFF) = ICONDL
        ENDDO
        DO IPHAS = 1, NPHAS
          IF( ITYTUR(IPHAS).EQ.2 .OR. ITYTUR(IPHAS).EQ.4
     &         .OR. ITURB(IPHAS).EQ.60 ) THEN
            IVAR = IU(IPHAS)
            ICONDL = ICONDL + 1
            ICLRTP(IVAR,ICOEFF) = ICONDL
            IVAR = IV(IPHAS)
            ICONDL = ICONDL + 1
            ICLRTP(IVAR,ICOEFF) = ICONDL
            IVAR = IW(IPHAS)
            ICONDL = ICONDL + 1
            ICLRTP(IVAR,ICOEFF) = ICONDL
          ENDIF
        ENDDO
        IF (IPHYDR.EQ.1) THEN
          DO IPHAS = 1, NPHAS
            ICONDL = ICONDL + 1
            ICLRTP(IPR(IPHAS),ICOEFF) = ICONDL
          ENDDO
        ENDIF
C
C Compressible
        IF (IPPMOD(ICOMPF).GE.0) THEN
          DO IPHAS = 1, NPHAS
            ICONDL = ICONDL + 1
            ICLRTP(ISCA(IENERG(IPHAS)),ICOEFF) = ICONDL
          ENDDO
        ENDIF
C
        NCOFAB = ICONDL
C
C
C ---> 4.5 POINTEURS POST-PROCESSING / LISTING / HISTORIQUES / CHRONOS
C      ---------------------------------------------------------------------
C
C --- Les pointeurs ont ete initialises a 1 (poubelle).
C       dans iniini (avant usini1).
C
C     On posttraitera les variables localisees au centre des cellules.
C
C      IPPRTP(IVAR)  pour RTP a ete complete plus haut.
C
C      IPPPRO(IPPROC(II)) pour PROPCE a ete complete plus haut
C        au fur et a mesure (voir ppprop en particulier)
C
C      Le rang de la derniere propriete pour le post est IPPPST.
C
C
        IF(IDTVAR.NE.0) THEN
          IPPPST = IPPPST + 1
          IPPDT  = IPPPST
        ENDIF
C
        IF(IPUCOU.EQ.1) THEN
          IPPPST = IPPPST + 1
          IPPTX  = IPPPST
          IPPPST = IPPPST + 1
          IPPTY  = IPPPST
          IPPPST = IPPPST + 1
          IPPTZ  = IPPPST
        ENDIF
C
C
C Verification de la limite sur IPPPST
C
        IF(IPPPST.GT.NVPPMX) THEN
          WRITE(NFECRA,8900)IPPPST,NVPPMX
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
        RETURN
C
      ENDIF
C
C=======================================================================
C 5. CINQUIEME APPEL :
C        RESERVATION D'UNE PLACE DANS PROPCE SI RAYONNEMENT
C        ET LAGRANGIEN AVEC THERMIQUE DES PARTICULES
C=======================================================================
C
      IF (IPASS.EQ.5) THEN
C
        IF ( IIRAYO.EQ.1 .AND.
     &       IILAGR.GT.0 .AND.
     &       IPHYLA.EQ.1 .AND.
     &       ITPVAR.EQ.1       ) THEN
C
          IPROP  = NPRMAX
          ILUMN  = IPROP + 1
          NPRMAX = IPROP
C
          IPROP = NPROCE
          IPROP = IPROP + 1
          IPPROC(ILUMN) = IPROP
          NPROCE = IPROP
C
C
C --- Verification de NPROCE, NPROFA, NPROFB
C
          IF(NPROCE.GT.NPROMX.OR.
     &       NPROFA.GT.NPROMX.OR.NPROFB.GT.NPROMX) THEN
            WRITE(NFECRA,7200)NPROCE, NPROFA, NPROFB, NPROMX,
     &           MAX(MAX(NPROCE,NPROFA),NPROFB)
            CALL CSEXIT (1)
C           ===========
          ENDIF
C
        ENDIF
C
        RETURN
C
      ENDIF
C
C
C=======================================================================
C 6. FORMATS
C=======================================================================
C
 5000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     NOMBRE DE PHASES ERRONE                                ',/,
     &'@                                                            ',/,
     &'@  Le nombre de phases doit etre un entier strictement       ',/,
     &'@    positif. Il vaut ici                NPHAS  = ',I10       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 5001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     NOMBRE DE PHASES TROP GRAND                            ',/,
     &'@                                                            ',/,
     &'@  Le nombre de phases                                       ',/,
     &'@    - demande          dans usini1   est NPHAS  = ',I10      ,/,
     &'@    - maximal autorise dans paramx.h est NPHSMX = ',I10      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@  NPHSMX doit valoir au moins ',I10                          ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 6000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     PLUSIEURS MODELES PHYSIQUES PARTICULIERES ACTIVES      ',/,
     &'@                                                            ',/,
     &'@  Un seul modele physique particuliere peut etre active a la',/,
     &'@    fois.                                                   ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Modifier les indicateurs de IPPMOD dans usppmo.           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 6001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     SELECTION INCORRECTE DU MODELE PHYSIQUE PARTICULIERE   ',/,
     &'@                                                            ',/,
     &'@  Les valeurs des indicateurs du tableau IPPMOD ne sont pas ',/,
     &'@    admissibles                                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Modifier les indicateurs de IPPMOD dans usppmo.           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 6010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     NOMBRE DE SCALAIRES ERRONE                             ',/,
     &'@                                                            ',/,
     &'@  Le nombre de scalaires utilisateur doit etre un entier    ',/,
     &'@    positif ou nul. Il vaut ici   NSCAUS  = ',I10            ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 6011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     NOMBRE DE SCALAIRES TROP GRAND                         ',/,
     &'@                                                            ',/,
     &'@  Le nombre de scalaires utilisateurs                       ',/,
     &'@    demande    dans usini1           est  NSCAUS = ',I10     ,/,
     &'@  Le nombre de scalaires total                              ',/,
     &'@    autorise   dans paramx.h         est  NSCAMX = ',I10     ,/,
     &'@                                                            ',/,
     &'@  La valeur maximale autorisee de NSCAUS                    ',/,
     &'@                          est donc  NSCAMX        = ',I10    ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier NSCAUS dans usini1.                              ',/,
     &'@                                                            ',/,
     &'@  NSCAMX doit valoir au moins ',I10                          ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 6012 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     NOMBRE DE SCALAIRES TROP GRAND                         ',/,
     &'@                                                            ',/,
     &'@  Le nombre de scalaires utilisateurs                       ',/,
     &'@    demande    dans usini1           est  NSCAUS = ',I10     ,/,
     &'@  Le nombre de scalaires pour les physiques particulieres   ',/,
     &'@    necessaire avec le modele choisi est  NSCAPP = ',I10     ,/,
     &'@  Le nombre de scalaires total                              ',/,
     &'@    autorise   dans paramx.h         est  NSCAMX = ',I10     ,/,
     &'@                                                            ',/,
     &'@  La valeur maximale autorisee de NSCAUS                    ',/,
     &'@    avec le modele choisi est donc NSCAMX-NSCAPP = ',I10     ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier NSCAUS dans usini1.                              ',/,
     &'@                                                            ',/,
     &'@  NSCAMX doit valoir au moins ',I10                          ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     NOMBRE DE SCALAIRES ERRONE                             ',/,
     &'@                                                            ',/,
     &'@  Le nombre de scalaires doit etre un entier                ',/,
     &'@    positif ou nul. Il vaut ici   NSCAL   = ',I10            ,/,
     &'@    Remarque : NSCAUS = ',I10                                ,/,
     &'@               NSCAPP = ',I10                                ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     NOMBRE DE SCALAIRES TROP GRAND                         ',/,
     &'@                                                            ',/,
     &'@  Le nombre de scalaires necessaire  est  NSCAL  = ',I10     ,/,
     &'@  Le nombre de scalaires total                              ',/,
     &'@    autorise   dans paramx.h         est  NSCAMX = ',I10     ,/,
     &'@                                                            ',/,
     &'@    Remarque : NSCAUS = ',I10                                ,/,
     &'@               NSCAPP = ',I10                                ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@  NSCAMX doit valoir au moins ',I10                          ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE PORTEUSE INCORRECTE POUR LE SCALAIRE ',I10         ,/,
     &'@                                                            ',/,
     &'@  La phase porteuse du scalaire ',I10                        ,/,
     &'@          (scalaire utilisateur ',I10   ,') indiquee dans   ',/,
     &'@    usini1 est IPHSCA(',I10   ,') = ',I10                    ,/,
     &'@  Elle devrait etre comprise entre 0 et NPHAS = ',I10        ,/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier IPHSCA dans usini1.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7021 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE PORTEUSE INCORRECTE POUR LE SCALAIRE ',I10         ,/,
     &'@                                                            ',/,
     &'@  La phase porteuse du scalaire ',I10                        ,/,
     &'@          (scalaire physique particuliere ',I10   ,')       ',/,
     &'@           est IPHSCA(ISCAPP(',I10   ,')) = ',I10            ,/,
     &'@  Elle devrait etre comprise entre 0 et NPHAS = ',I10        ,/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier IPHSCA.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7030 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ASSOCIE A UNE VARIANCE INCORRECT               ',/,
     &'@                                                            ',/,
     &'@  Le scalaire ',I10                                          ,/,
     &'@    (scalaire utilisateur           ',I10   ,') represente  ',/,
     &'@     une variance puisque                                   ',/,
     &'@    ISCAVR(',I10   ,') vaut ',I10   ,' (non nul)            ',/,
     &'@  Les valeurs de ISCAVR doivent cependant etre              ',/,
     &'@    superieures ou egales a                  0              ',/,
     &'@    inferieures ou egales a NSCAL = ',I10                    ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier ISCAVR dans usini1.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7031 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ASSOCIE A UNE VARIANCE INCORRECT               ',/,
     &'@                                                            ',/,
     &'@  Le scalaire ',I10                                          ,/,
     &'@    (scalaire physique particuliere ',I10   ,') represente  ',/,
     &'@     une variance puisque    ',/,
     &'@    ISCAVR(ISCAPP(',I10   ,')) vaut ',I10   ,' (non nul)    ',/,
     &'@  Les valeurs de ISCAVR doivent cependant etre              ',/,
     &'@    superieures ou egales a                  0              ',/,
     &'@    inferieures ou egales a NSCAL = ',I10                    ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier ISCAVR dans usini1.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7040 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',I10                                           ,/,
     &'@                                                            ',/,
     &'@  Le scalaire ',I10                                          ,/,
     &'@    (scalaire utilisateur           ',I10   ,') represente  ',/,
     &'@    la variance des fluctuations du scalaire ',I10           ,/,
     &'@    (scalaire utilisateur           ',I10   ,') puisque     ',/,
     &'@    ISCAVR(',I10   ,') vaut ',I10   ,' (non nul)            ',/,
     &'@                                                            ',/,
     &'@  L''indicateur de diffusivite IVISLS(',I10  ,')            ',/,
     &'@    ne doit pas etre renseigne.                             ',/,
     &'@  Il sera pris automatiquement egal a celui du scalaire     ',/,
     &'@    associe, soit ',I10                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier IVISLS dans usini1.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7041 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',I10                                           ,/,
     &'@                                                            ',/,
     &'@  Le scalaire ',I10                                          ,/,
     &'@    (scalaire utilisateur           ',I10   ,') represente  ',/,
     &'@    la variance des fluctuations du scalaire ',I10           ,/,
     &'@    (scalaire physique particuliere ',I10   ,') puisque     ',/,
     &'@    ISCAVR(',I10   ,') vaut ',I10   ,' (non nul)            ',/,
     &'@                                                            ',/,
     &'@  L''indicateur de diffusivite IVISLS(',I10  ,')            ',/,
     &'@    ne doit pas etre renseigne.                             ',/,
     &'@  Il sera pris automatiquement egal a celui du scalaire     ',/,
     &'@    associe, soit ',I10                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier IVISLS dans usini1.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7042 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',I10                                           ,/,
     &'@                                                            ',/,
     &'@  Le scalaire ',I10                                          ,/,
     &'@    (scalaire physique particuliere ',I10   ,') represente  ',/,
     &'@    la variance des fluctuations du scalaire ',I10           ,/,
     &'@    (scalaire utilisateur           ',I10   ,') puisque     ',/,
     &'@    ISCAVR(ISCAPP(',I10   ,')) vaut ',I10   ,' (non nul)    ',/,
     &'@                                                            ',/,
     &'@  L''indicateur de diffusivite IVISLS(ISCAPP(',I10  ,'))    ',/,
     &'@    ne doit pas etre renseigne.                             ',/,
     &'@  Il sera pris automatiquement egal a celui du scalaire     ',/,
     &'@    associe, soit ',I10                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier IVISLS.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7043 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',I10                                           ,/,
     &'@                                                            ',/,
     &'@  Le scalaire ',I10                                          ,/,
     &'@    (scalaire physique particuliere ',I10   ,') represente  ',/,
     &'@    la variance des fluctuations du scalaire ',I10           ,/,
     &'@    (scalaire physique particuliere ',I10   ,') puisque     ',/,
     &'@    ISCAVR(ISCAPP(',I10   ,')) vaut ',I10   ,' (non nul)    ',/,
     &'@                                                            ',/,
     &'@  L''indicateur de diffusivite IVISLS(ISCAPP(',I10  ,'))    ',/,
     &'@    ne doit pas etre renseigne.                             ',/,
     &'@  Il sera pris automatiquement egal a celui du scalaire     ',/,
     &'@    associe, soit ',I10                                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier IVISLS.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7050 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',I10                                           ,/,
     &'@                                                            ',/,
     &'@  L''indicateur de diffusivite variable                     ',/,
     &'@    du scalaire utilisateur           ',I10                  ,/,
     &'@    IVISLS(',I10   ,')                                      ',/,
     &'@    doit etre un entier egal a 0 ou 1.                      ',/,
     &'@    Il vaut ici ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier IVISLS dans usini1.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7051 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',I10                                           ,/,
     &'@                                                            ',/,
     &'@  L''indicateur de diffusivite variable                     ',/,
     &'@    du scalaire physique particuliere ',I10                  ,/,
     &'@    IVISLS(ISCAPP(',I10   ,'))                              ',/,
     &'@    doit etre un entier egal a 0 ou 1.                      ',/,
     &'@    Il vaut ici ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier IVISLS.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7060 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',I10   ,' NE PAS MODIFIER LA PHASE PORTEUSE   ',/,
     &'@                                                            ',/,
     &'@  Le scalaire ',I10                                          ,/,
     &'@    (scalaire utilisateur           ',I10   ,') represente  ',/,
     &'@    la variance des fluctuations du scalaire ',I10           ,/,
     &'@    (scalaire utilisateur           ',I10   ,') puisque     ',/,
     &'@    ISCAVR(',I10   ,') vaut ',I10   ,' (non nul)            ',/,
     &'@                                                            ',/,
     &'@  La phase porteuse IPHSCA(',I10   ,') du scalaire          ',/,
     &'@    ne doit pas etre renseignee.                            ',/,
     &'@  Elle sera automatiquement prise identique a la phase      ',/,
     &'@    porteuse du scalaire associe, soit ',I10                 ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier IPHSCA dans usini1.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7061 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',I10   ,' NE PAS MODIFIER LA PHASE PORTEUSE   ',/,
     &'@                                                            ',/,
     &'@  Le scalaire ',I10                                          ,/,
     &'@    (scalaire utilisateur           ',I10   ,') represente  ',/,
     &'@    la variance des fluctuations du scalaire ',I10           ,/,
     &'@    (scalaire physique particuliere ',I10   ,') puisque     ',/,
     &'@    ISCAVR(',I10   ,') vaut ',I10   ,' (non nul)            ',/,
     &'@                                                            ',/,
     &'@  La phase porteuse IPHSCA(',I10   ,') du scalaire          ',/,
     &'@    ne doit pas etre renseignee.                            ',/,
     &'@  Elle sera automatiquement prise identique a la phase      ',/,
     &'@    porteuse du scalaire associe, soit ',I10                 ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier IPHSCA dans usini1.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7062 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',I10   ,' NE PAS MODIFIER LA PHASE PORTEUSE   ',/,
     &'@                                                            ',/,
     &'@  Le scalaire ',I10                                          ,/,
     &'@    (scalaire physique particuliere ',I10   ,') represente  ',/,
     &'@    la variance des fluctuations du scalaire ',I10           ,/,
     &'@    (scalaire utilisateur           ',I10   ,') puisque     ',/,
     &'@    ISCAVR(ISCAPP(',I10   ,')) vaut ',I10   ,' (non nul)    ',/,
     &'@                                                            ',/,
     &'@  La phase porteuse IPHSCA(ISCAPP(',I10   ,')) du scalaire  ',/,
     &'@    ne doit pas etre renseignee.                            ',/,
     &'@  Elle sera automatiquement prise identique a la phase      ',/,
     &'@    porteuse du scalaire associe, soit ',I10                 ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier IPHSCA.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7063 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',I10   ,' NE PAS MODIFIER LA PHASE PORTEUSE   ',/,
     &'@                                                            ',/,
     &'@  Le scalaire ',I10                                          ,/,
     &'@    (scalaire physique particuliere ',I10   ,') represente  ',/,
     &'@    la variance des fluctuations du scalaire ',I10           ,/,
     &'@    (scalaire physique particuliere ',I10   ,') puisque     ',/,
     &'@    ISCAVR(ISCAPP(',I10   ,')) vaut ',I10   ,' (non nul)    ',/,
     &'@                                                            ',/,
     &'@  La phase porteuse IPHSCA(ISCAPP(',I10   ,')) du scalaire  ',/,
     &'@    ne doit pas etre renseignee.                            ',/,
     &'@  Elle sera automatiquement prise identique a la phase      ',/,
     &'@    porteuse du scalaire associe, soit ',I10                 ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier IPHSCA.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7070 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR IORTVM NE PEUT PRENDRE QUE LES VALEURS    ',/,
     &'@      0 OU 1.                                               ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usalin.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     NOMBRE DE VARIABLES TROP GRAND                         ',/,
     &'@                                                            ',/,
     &'@  Le type de calcul defini dans usini1                      ',/,
     &'@    correspond a un nombre de variables NVAR   = ',I10       ,/,
     &'@  Le nombre de variables maximal prevu                      ',/,
     &'@                      dans paramx.h est NVARMX = ',I10       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@  NVARMX doit valoir au moins ',I10                          ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     NOMBRE DE PROPRIETES TROP GRAND                        ',/,
     &'@                                                            ',/,
     &'@  Le type de calcul defini dans usini1                      ',/,
     &'@    correspond aux nombres de proprietes suivants           ',/,
     &'@      au centre des cellules       : NPROCE = ',I10          ,/,
     &'@      au centre des faces internes : NPROFA = ',I10          ,/,
     &'@      au centre des faces de bord  : NPROFB = ',I10          ,/,
     &'@  Le nombre de proprietes maximal prevu                     ',/,
     &'@                      dans paramx.h est NPROMX = ',I10       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@  NPROMX doit valoir au moins ',I10                          ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8021 FORMAT(
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
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8022 FORMAT(
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
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8101 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 1 ou 2                ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8111 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,
     &'@    *********                                               ',/,
     &'@    AVEC UN SCHEMA EN TEMPS D ORDRE 2 : ISCHTP = ', I10      ,/,
     &'@    IL FAUT UTILISER UN PAS DE TEMPS CONSTANT ET UNIFORME   ',/,
     &'@    OR IDTVAR = ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8112 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,
     &'@    *********                                               ',/,
     &'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,
     &'@    EN K-EPSILON (ITURB = ',I10,' )'                         ,/,
     &'@                                                            ',/,
     &'@   La version courante ne supporte pas l ordre 2 avec le    ',/,
     &'@   couplage des termes sources du k-epsilon.                ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Modifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8113 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES PHASE ',I10    ,/,
     &'@    *********                                               ',/,
     &'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 1 (ISCHTP = ',I10   ,/,
     &'@    EN LES (ITURB = ',I10,' )'                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute.                                   ',/,
     &'@                                                            ',/,
     &'@  Il est conseille de verifier usini1.                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8114 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,
     &'@    *********                                               ',/,
     &'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,
     &'@    EN V2F       (ITURB = ',I10,' )'                         ,/,
     &'@                                                            ',/,
     &'@   La version courante ne supporte pas l''ordre 2 avec le   ',/,
     &'@   couplage des termes sources du k-epsilon.                ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Modifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8115 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,
     &'@    *********                                               ',/,
     &'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,
     &'@    EN K-OMEGA   (ITURB = ',I10,' )'                         ,/,
     &'@                                                            ',/,
     &'@   La version courante ne supporte pas l''ordre 2 avec le   ',/,
     &'@   couplage des termes sources du k-omega.                  ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Modifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8121 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN ENTIER EGAL A  0, 1 OU 2            ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8131 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0, 1 OU 2             ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8132 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@  Avec le modele de turbulence choisi, ITURB = ',I10         ,/,
     &'@    la valeur de ISTO2T (extrapolation des termes sources   ',/,
     &'@    pour les variables turbulentes) ne doit pas etre modifie',/,
     &'@    or ISTO2T a ete force a ',I10                            ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8141 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES SCALAIRE ',I10 ,/,
     &'@    *********                                               ',/,
     &'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0, 1 OU 2             ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,
     &'@    *********                                               ',/,
     &'@      SUR LA LISTE DES MOYENNES TEMPORELLES                 ',/,
     &'@                                                            ',/,
     &'@    La valeur de IDFMOM(1,',I10   ,' est ',I10               ,/,
     &'@      ceci indique que l''on a defini ',I10   ,' moyennes   ',/,
     &'@      temporelles en renseignant les IMOM = ',I10            ,/,
     &'@      premieres cases du tableau IDFMOM(.,IMOM).            ',/,
     &'@    Les cases suivantes devraient etre nulles.              ',/,
     &'@                                                            ',/,
     &'@    Ce n''est cependant pas le cas :                        ',/,
     &'@                                                            ',/,
     &'@        IMOM    IDFMOM(1,IMOM)                              ',/,
     &'@  ----------------------------                              '  )
 8201 FORMAT(
     &'@  ',I10   ,'        ',     I10                                 )
 8202 FORMAT(
     &'@  ----------------------------                              ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas  execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier usini1.                                        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8210 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,
     &'@    *********                                               ',/,
     &'@      SUR LES VARIABLES COMPOSANT LES MOYENNES TEMPORELLES  ',/,
     &'@                                                            ',/,
     &'@    IDFMOM(',I10   ,',',I10   ,') = ',I10                    ,/,
     &'@      Les valeurs negatives renvoient a des proprietes      ',/,
     &'@      physiques, or il n y en a que NPRMAX = ', I10          ,/,
     &'@    La valeur de IDFMOM est donc erronee.                   ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier usini1.                                        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8211 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,
     &'@    *********                                               ',/,
     &'@      SUR LES VARIABLES COMPOSANT LES MOYENNES TEMPORELLES  ',/,
     &'@                                                            ',/,
     &'@    IDFMOM(',I10   ,',',I10   ,') = ',I10                    ,/,
     &'@      Les valeurs positives renvoient a des variables de    ',/,
     &'@      calcul, or il n y en a que NVAR   = ', I10             ,/,
     &'@    La valeur de IDFMOM est donc erronee.                   ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier usini1.                                        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8212 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,
     &'@    *********                                               ',/,
     &'@      SUR LES VARIABLES COMPOSANT LES MOYENNES TEMPORELLES  ',/,
     &'@                                                            ',/,
     &'@    La valeur                                               ',/,
     &'@      IDFMOM(',I10   ,',',I10   ,') = ',I10                  ,/,
     &'@      n''est pas une propriete associee aux cellules        ',/,
     &'@      (IPPROC(',I10   ,') = ',I10   ,')                     ',/,
     &'@    La valeur de IDFMOM est donc erronee.                   ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier usini1.                                        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8213 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,
     &'@    *********                                               ',/,
     &'@      SUR LES VARIABLES COMPOSANT LES MOYENNES TEMPORELLES  ',/,
     &'@                                                            ',/,
     &'@    Le tableau IDFMOM(JJ,IMOM) pour IMOM = ',I10             ,/,
     &'@      doit etre renseigne continuement. Or ici,             ',/,
     &'@      IDFMOM(',I10,',IMOM) est non nul (=',I10   ,')        ',/,
     &'@      alors qu il existe II < JJ pour lequel                ',/,
     &'@      IDFMOM(II,IMOM) est nul.                              ',/,
     &'@    La valeur de IDFMOM est donc erronee.                   ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier usini1.                                        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8214 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,
     &'@    *********                                               ',/,
     &'@  SUR L''INSTANT DE DEBUT DE CALCUL DES MOYENNES TEMPORELLES',/,
     &'@                                                            ',/,
     &'@    La variable NTDMOM(IMOM)   pour IMOM = ',I10             ,/,
     &'@      doit etre renseignee dans usini1 pour indiquer        ',/,
     &'@      a partir de quel pas de temps (absolu) doit etre      ',/,
     &'@      calculee la moyenne temporelle IMOM correspondante.   ',/,
     &'@      NTDMOM(IMOM) vaur ici ',I10                            ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier usini1.                                        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8300 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@      ERREUR A L''OUVERTURE DU FICHIER SUITE AUXILIAIRE     ',/,
     &'@                                                            ',/,
     &'@    Pour permettre de realiser une suite de calcul en       ',/,
     &'@      prenant en compte les moyennes temporelles,           ',/,
     &'@      on cherche a relire le fichier suite auxiliaire.      ',/,
     &'@    Une erreur se produit a son ouverture.                  ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier l''existence et le nom (',A13,') du            ',/,
     &'@        fichier suite dans le repertoire de travail.        ',/,
     &'@    Il est possible de s''affranchir du fichier suite       ',/,
     &'@      auxiliaire en reinitialisant les moyennes (IMOOLD)    ',/,
     &'@      ou en ne calculant pas de moyennes.                   ',/,
     &'@      Voir alors usini1.                                    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8301 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@      TYPE DE FICHIER INCORRECT                             ',/,
     &'@                                                            ',/,
     &'@    Le fichier ',A13      ,' ne semble pas etre un fichier  ',/,
     &'@      suite auxiliaire.                                     ',/,
     &'@                                                            ',/,
     &'@    Pour permettre de realiser une suite de calcul en       ',/,
     &'@      prenant en compte les moyennes temporelles,           ',/,
     &'@      on cherche a relire le fichier suite auxiliaire.      ',/,
     &'@    Une erreur se produit a sa lecture.                     ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier que le fichier suite utilise correspond bien   ',/,
     &'@        a un fichier suite auxiliaire.                      ',/,
     &'@    Il est possible de s''affranchir du fichier suite       ',/,
     &'@      auxiliaire en reinitialisant les moyennes (IMOOLD)    ',/,
     &'@      ou en ne calculant pas de moyennes.                   ',/,
     &'@      Voir alors usini1.                                    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8302 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,
     &'@                                                            ',/,
     &'@    Le nombre de cellules a ete modifie                     ',/,
     &'@                                                            ',/,
     &'@    Pour permettre de realiser une suite de calcul en       ',/,
     &'@      prenant en compte les moyennes temporelles,           ',/,
     &'@      on cherche a relire le fichier suite auxiliaire.      ',/,
     &'@    Une erreur se produit a sa lecture.                     ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Il est possible de s''affranchir du fichier suite       ',/,
     &'@      auxiliaire en reinitialisant les moyennes (IMOOLD)    ',/,
     &'@      ou en ne calculant pas de moyennes.                   ',/,
     &'@      Voir alors usini1.                                    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8311 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********  AUXILIAIRE (varpos)                          ',/,
     &'@                                                            ',/,
     &'@    Pour permettre de realiser une suite de calcul en       ',/,
     &'@      prenant en compte les moyennes temporelles,           ',/,
     &'@      on cherche a relire le fichier suite auxiliaire.      ',/,
     &'@    Erreur a la lecture du nombre de moyennes               ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier le fichier suite auxiliaire.                   ',/,
     &'@    Il est possible de s''affranchir du fichier suite       ',/,
     &'@      auxiliaire en reinitialisant les moyennes (IMOOLD)    ',/,
     &'@      ou en ne calculant pas de moyennes.                   ',/,
     &'@      Voir alors usini1.                                    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8312 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********  AUXILIAIRE (varpos)                          ',/,
     &'@                                                            ',/,
     &'@    Pour permettre de realiser une suite de calcul en       ',/,
     &'@      prenant en compte les moyennes temporelles,           ',/,
     &'@      on cherche a relire le fichier suite auxiliaire.      ',/,
     &'@    Erreur a la lecture du numero du pas de temps precedent ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier le fichier suite auxiliaire.                   ',/,
     &'@    Il est possible de s''affranchir du fichier suite       ',/,
     &'@      auxiliaire en reinitialisant les moyennes (IMOOLD)    ',/,
     &'@      ou en ne calculant pas de moyennes.                   ',/,
     &'@      Voir alors usini1.                                    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8313 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,
     &'@    *********  AUXILIAIRE (varpos)                          ',/,
     &'@                                                            ',/,
     &'@    Pour permettre de realiser une suite de calcul en       ',/,
     &'@      prenant en compte les moyennes temporelles,           ',/,
     &'@      on cherche a relire le fichier suite auxiliaire.      ',/,
     &'@    Erreur a la lecture du numero du cumul temporel         ',/,
     &'@      de la moyenne                  ', I10                  ,/,
     &'@      associee a l''ancienne moyenne ', I10                  ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Verifier le fichier suite auxiliaire.                   ',/,
     &'@    Il est possible de s''affranchir du fichier suite       ',/,
     &'@      auxiliaire en reinitialisant les moyennes (IMOOLD)    ',/,
     &'@      ou en ne calculant pas de moyennes.                   ',/,
     &'@      Voir alors usini1.                                    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8390 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR A LA FERMETURE DU FICHIER SUITE      ',/,
     &'@    *********   AUXILIAIRE (varpos)                         ',/,
     &'@                                                            ',/,
     &'@    Probleme sur le fichier de nom (',A13,')                ',/,
     &'@                                                            ',/,
     &'@    Le calcul se poursuit...                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8400 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE (VARPOS)',/,
     &'@    *********                                     AUXILIAIRE',/,
     &'@      ERREUR A LA LECTURE DU FICHIER SUITE                  ',/,
     &'@                                                            ',/,
     &'@    On souhaite faire correspondre la moyenne   ',I10        ,/,
     &'@            du present calcul avec la moyenne   ',I10        ,/,
     &'@            du calcul precedent, or, le numero des anciennes',/,
     &'@            moyennes  doit etre strictement positif et      ',/,
     &'@            inferieur ou egal a                 ',I10        ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut etre execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier les valeurs de IMOOLD dans usini1.             ',/,
     &'@    Verifier que le fichier suite est le bon.               ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8401 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,
     &'@    *********                                               ',/,
     &'@      SUR LA CORRESPONDANCE AVEC LES ANCIENNES MOYENNES     ',/,
     &'@                                                            ',/,
     &'@    On souhaite initialiser la nouvelle moyenne IMOM = ',I10 ,/,
     &'@      en relisant la moyenne IMOOLD(IMOM) = ',I10            ,/,
     &'@      dans le fichier suite.                                ',/,
     &'@    Or on a specifie que le pas de temps initial pour le    ',/,
     &'@      calcul de la moyenne IMOM etait NTDMOM(IMOM) = ',I10   ,/,
     &'@      et le calcul stocke dans le fichier suite correspond  ',/,
     &'@      au pas de temps = ',I10                                ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Modifier la valeur de  NTDMOM dans usini1               ',/,
     &'@      verifier que le fichier suite est le bon.             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8402 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,
     &'@    *********                                               ',/,
     &'@      SUR LA CORRESPONDANCE AVEC LES ANCIENNES MOYENNES     ',/,
     &'@                                                            ',/,
     &'@    Deux moyennes distinctes ',I10,' et ', I10               ,/,
     &'@      sont initialises avec la meme ',I10                    ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Verifier IMOOLD dans usini1      .                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8403 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,
     &'@    *********                                               ',/,
     &'@      SUR LA CORRESPONDANCE AVEC LES ANCIENNES MOYENNES     ',/,
     &'@                                                            ',/,
     &'@    On ne souhaite pas faire un calcul suite puisque        ',/,
     &'@      ISUITE = ',I10                                         ,/,
     &'@      mais on a renseigne le tableau IMOOLD de              ',/,
     &'@      correspondance nouvelles -> anciennes moyennes :      ',/,
     &'@      IMOOLD(',I10,') = ',I10                                ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Ne pas modifier IMOOLD dans usini1.                     ',/,
     &'@      ou realiser un calcul suite.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8410 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA CONSTRUCTION DE IDTMOM  (varpos) ',/,
     &'@    *********                                               ',/,
     &'@      CALCUL DES MOYENNES TEMPORELLES                       ',/,
     &'@                                                            ',/,
     &'@    IDTMOM(',I10,') = ', I10                                 ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8900 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     NOMBRE DE VARIABLES A SUIVRE TROP GRAND                ',/,
     &'@                                                            ',/,
     &'@  Le type de calcul defini dans usini1                      ',/,
     &'@    correspond a un nombre de variables a suivre dans       ',/,
     &'@    le listing et le post-processing egal a      ',I10       ,/,
     &'@  Le nombre de variables a suivre maximal prevu             ',/,
     &'@                      dans paramx.h est NVPPMX = ',I10       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@  Contacter l assistance.                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
C
C=======================================================================
C 5. FIN
C=======================================================================
C
      RETURN
      END
c@z
