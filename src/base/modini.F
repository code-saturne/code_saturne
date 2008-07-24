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
                        SUBROUTINE MODINI
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  MODIFICATION DES PARAMETRES DE CALCUL
CFONC     APRES INTERVENTION UTILISATEUR
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
C
C***********************************************************************
C
C ARGUMENTS
C
C
C VARIABLES LOCALES
C
      INTEGER          II, JJ, IVAR, IPHAS, IOK, IEST, IMOM, IKW
      INTEGER          ICOMPT, IPP, NBCCOU, NN
      INTEGER          NSCACP, ISCAL
      DOUBLE PRECISION RELXSP
C
C***********************************************************************
C
C Indicateur erreur (0 : pas d'erreur)
      IOK = 0
C
C=======================================================================
C 1. ENTREES SORTIES entsor.h
C=======================================================================
C
C ---> Niveau d'impression listing
C       Non initialise -> standard
      DO II = 1, NVARMX
        IF(IWARNI(II).EQ.-10000) THEN
          IWARNI(II) = 0
        ENDIF
      ENDDO
C
C---> Rayonnement
C       Format des fichiers suite
C
      IF(IFOAVR.EQ.-1) THEN
        IFOAVR = IFOAVA
      ENDIF
C
C---> Lagrangien
C       Format des fichiers suite
C
      IF(IFOAVL.EQ.-1) THEN
        IFOAVL = IFOAVA
      ENDIF
      IF(IFOVLS.EQ.-1) THEN
        IFOVLS = IFOAVA
      ENDIF
C
C---> Variables de calcul ITRSVR = ivar. Sinon 0 (valeur d'initialisation).
C
      DO IVAR = 1, NVAR
        ITRSVR(IPPRTP(IVAR)) = IVAR
      ENDDO
C
C---> sorties chrono?
C     Sauf mention contraire de l'utilisateur, on sort a la fin les
C        variables de calcul, la viscosite, rho, le pas de temps s'il
C        est variable, les estimateurs s'ils sont actives, les moments
C        s'il y en a et la viscosite de maillage en ALE.
C
      DO II = 2, NVPPMX
        IF(ITRSVR(II).GE.1.AND.ICHRVR(II).EQ.-999) THEN
          ICHRVR(II) = 1
        ENDIF
      ENDDO
      DO IPHAS = 1, NPHAS
        IPP = IPPPRO(IPPROC(IROM(IPHAS)))
        IF(                       ICHRVR(IPP).EQ.-999) ICHRVR(IPP) = 1
        IPP = IPPPRO(IPPROC(IVISCT(IPHAS)))
        IF( (ITURB(IPHAS).EQ.10 .OR. ITYTUR(IPHAS).EQ.2
     &       .OR. ITURB(IPHAS).EQ.50 .OR. ITURB(IPHAS).EQ.60)
     &       .AND.ICHRVR(IPP).EQ.-999) ICHRVR(IPP) = 1
        IF (IDTVAR.LT.0) THEN
          ICHRVR(IPPPRO(IPPROC(ICOUR(IPHAS)))) = 0
          ICHRVR(IPPPRO(IPPROC(IFOUR(IPHAS)))) = 0
        ENDIF
        DO IEST = 1, NESTMX
          IF(IESCAL(IEST,IPHAS).GT.0) THEN
            IPP = IPPPRO(IPPROC(IESTIM(IEST,IPHAS)))
            IF(                     ICHRVR(IPP).EQ.-999) ICHRVR(IPP) = 1
          ENDIF
        ENDDO
      ENDDO
      IPP = IPPDT
      IF (IDTVAR.LT.0) ICHRVR(IPP) = 0
      IF(IDTVAR.EQ.2.AND.ICHRVR(IPP).EQ.-999) ICHRVR(IPP) = 1
      IF(IPUCOU.NE.1) THEN
        ICHRVR(IPPTX) = 0
        ICHRVR(IPPTY) = 0
        ICHRVR(IPPTZ) = 0
      ENDIF
C
      IF(NBMOMT.GT.0) THEN
        DO IMOM = 1, NBMOMT
          IPP = IPPPRO(IPPROC(ICMOME(IMOM)))
          IF(ICHRVR(IPP).EQ.-999) ICHRVR(IPP) = 1
        ENDDO
      ENDIF
      IF (IALE.EQ.1) THEN
        NN = 1
        IF (IORTVM.EQ.1) NN = 3
        DO II = 1, NN
          IPP = IPPPRO(IPPROC(IVISMA(II)))
          IF(ICHRVR(IPP).EQ.-999) ICHRVR(IPP) = 1
        ENDDO
      ENDIF
C
      DO II = 1, NVPPMX
        IF(ICHRVR(II).EQ.-999) THEN
          ICHRVR(II) = 0
        ENDIF
      ENDDO
C
      ICOMPT = 0
      DO II = 2, NVPPMX
        IF(ICHRVR(II).EQ.1) ICOMPT = ICOMPT+1
      ENDDO
      IF(ICOMPT.EQ.0) NTCHR = -1
C
C
C---> sorties historiques ?
C      Si une valeur non modifiee par l'utilisateur (=-999)
C        on la met a sa valeur par defaut
C      On sort toutes les variables a tous les pas de temps par defaut
C      IHISVR nb de sonde et numero par variable (-999 non initialise)
C             -1 : toutes les sondes
C      NTHIST = -1 : on ne sort pas d'historiques
C      NTHIST =  n : on sort des historiques tous les n pas de temps
C      NTHSAV = -1 : on sauvegarde a la fin uniquement
C      NTHSAV =  0 : periode par defaut (voir caltri)
C             > 0  : periode
C
      DO II = 2, NVPPMX
        IF(ITRSVR(II).GE.1.AND.IHISVR(II,1).EQ.-999) THEN
          IHISVR(II,1) = -1
        ENDIF
      ENDDO
      IF (IDTVAR.LT.0) IHISVR(IPPDT ,1) = 0
      IF(IDTVAR.EQ.2.AND.IHISVR(IPPDT ,1).EQ.-999) IHISVR(IPPDT ,1) = -1
      IF(IPUCOU.NE.1) THEN
        IHISVR(IPPTX,1) = 0
        IHISVR(IPPTY,1) = 0
        IHISVR(IPPTZ,1) = 0
      ENDIF
      DO IPHAS = 1, NPHAS
        IPP = IPPPRO(IPPROC(IVISCT(IPHAS)))
        IF( (ITURB(IPHAS).EQ.10 .OR. ITYTUR(IPHAS).EQ.2
     &       .OR. ITURB(IPHAS).EQ.50 .OR. ITURB(IPHAS).EQ.60)
     &       .AND.IHISVR(IPP,1).EQ.-999) IHISVR(IPP,1) = -1
        IF (IDTVAR.LT.0) THEN
          IHISVR(IPPPRO(IPPROC(ICOUR(IPHAS))),1) = 0
          IHISVR(IPPPRO(IPPROC(IFOUR(IPHAS))),1) = 0
        ENDIF
      ENDDO
      IF(NBMOMT.GT.0) THEN
        DO IMOM = 1, NBMOMT
          IPP = IPPPRO(IPPROC(ICMOME(IMOM)))
          IF(IHISVR(IPP,1).EQ.-999) IHISVR(IPP,1) = -1
        ENDDO
      ENDIF
C
      DO II = 1, NVPPMX
        IF(IHISVR(II,1).EQ.-999) THEN
          IHISVR(II,1) = 0
        ENDIF
      ENDDO
C
C     Si on est en ALE, on a un test equivalent dans strini.F
      IF (IALE.EQ.0) THEN
        ICOMPT = 0
        DO II = 2, NVPPMX
          IF(IHISVR(II,1).NE.0) ICOMPT = ICOMPT+1
        ENDDO
C
        IF(ICOMPT.EQ.0.OR.NCAPT.EQ.0) THEN
          NTHIST = -1
        ENDIF
      ENDIF
C
C ---> Nom des variables
C
      DO IPHAS = 1, NPHAS
C
        IF(NOMVAR(IPPRTP(IPR   (IPHAS))) .EQ.' ') THEN
          WRITE(NOMVAR(IPPRTP(IPR   (IPHAS))),'(A6,I2.2)')'PresPh',IPHAS
        ENDIF
        IF(NOMVAR(IPPRTP(IU    (IPHAS))) .EQ.' ') THEN
          WRITE(NOMVAR(IPPRTP(IU    (IPHAS))),'(A6,I2.2)')'VitesX',IPHAS
        ENDIF
        IF(NOMVAR(IPPRTP(IV    (IPHAS))) .EQ.' ') THEN
          WRITE(NOMVAR(IPPRTP(IV    (IPHAS))),'(A6,I2.2)')'VitesY',IPHAS
        ENDIF
        IF(NOMVAR(IPPRTP(IW    (IPHAS))) .EQ.' ') THEN
          WRITE(NOMVAR(IPPRTP(IW    (IPHAS))),'(A6,I2.2)')'VitesZ',IPHAS
        ENDIF
        IF(ITYTUR(IPHAS).EQ.2) THEN
          IF(NOMVAR(IPPRTP(IK    (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IK    (IPHAS))),'(A6,I2.2)')
     &                                                    'EnTurb',IPHAS
          ENDIF
          IF(NOMVAR(IPPRTP(IEP   (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IEP   (IPHAS))),'(A6,I2.2)')
     &                                                    'Dissip',IPHAS
          ENDIF
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
          IF(NOMVAR(IPPRTP(IR11  (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IR11  (IPHAS))),'(A6,I2.2)')
     &                                                    'R11pha',IPHAS
          ENDIF
          IF(NOMVAR(IPPRTP(IR22  (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IR22  (IPHAS))),'(A6,I2.2)')
     &                                                    'R22pha',IPHAS
          ENDIF
          IF(NOMVAR(IPPRTP(IR33  (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IR33  (IPHAS))),'(A6,I2.2)')
     &                                                    'R33pha',IPHAS
          ENDIF
          IF(NOMVAR(IPPRTP(IR12  (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IR12  (IPHAS))),'(A6,I2.2)')
     &                                                    'R12pha',IPHAS
          ENDIF
          IF(NOMVAR(IPPRTP(IR13  (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IR13  (IPHAS))),'(A6,I2.2)')
     &                                                    'R13pha',IPHAS
          ENDIF
          IF(NOMVAR(IPPRTP(IR23  (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IR23  (IPHAS))),'(A6,I2.2)')
     &                                                    'R23pha',IPHAS
          ENDIF
          IF(NOMVAR(IPPRTP(IEP   (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IEP   (IPHAS))),'(A6,I2.2)')
     &                                                    'Dissip',IPHAS
          ENDIF
        ELSEIF(ITURB(IPHAS).EQ.50) THEN
          IF(NOMVAR(IPPRTP(IK    (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IK    (IPHAS))),'(A6,I2.2)')
     &                                                    'EnTurb',IPHAS
          ENDIF
          IF(NOMVAR(IPPRTP(IEP   (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IEP   (IPHAS))),'(A6,I2.2)')
     &                                                    'Dissip',IPHAS
          ENDIF
          IF(NOMVAR(IPPRTP(IPHI  (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IPHI  (IPHAS))),'(A6,I2.2)')
     &                                                    'phipha',IPHAS
          ENDIF
          IF(NOMVAR(IPPRTP(IFB   (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IFB   (IPHAS))),'(A6,I2.2)')
     &                                                    'fbarre',IPHAS
          ENDIF
        ELSEIF(ITURB(IPHAS).EQ.60) THEN
          IF(NOMVAR(IPPRTP(IK    (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IK    (IPHAS))),'(A6,I2.2)')
     &                                                    'EnTurb',IPHAS
          ENDIF
          IF(NOMVAR(IPPRTP(IOMG  (IPHAS))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPRTP(IOMG  (IPHAS))),'(A5,I2.2)')
     &                                                    'Omega',IPHAS
          ENDIF
        ENDIF
C
        IF(NOMVAR(IPPPRO(IPPROC(IROM  (IPHAS)))) .EQ.' ') THEN
          WRITE(NOMVAR(IPPPRO(IPPROC(IROM  (IPHAS)))),'(A6,I2.2)')
     &                                                    'MasVol',IPHAS
        ENDIF
        IF(NOMVAR(IPPPRO(IPPROC(IVISCT(IPHAS)))) .EQ.' ') THEN
          WRITE(NOMVAR(IPPPRO(IPPROC(IVISCT(IPHAS)))),'(A6,I2.2)')
     &                                                    'VisTur',IPHAS
        ENDIF
        IF(NOMVAR(IPPPRO(IPPROC(IVISCL(IPHAS)))) .EQ.' ') THEN
          WRITE(NOMVAR(IPPPRO(IPPROC(IVISCL(IPHAS)))),'(A6,I2.2)')
     &                                                    'VisMol',IPHAS
        ENDIF
        IF (ISMAGO(IPHAS).GT.0) THEN
          IF(NOMVAR(IPPPRO(IPPROC(ISMAGO(IPHAS)))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPPRO(IPPROC(ISMAGO(IPHAS)))),'(A6,I2.2)')
     &                                                    'Csdyn2',IPHAS
          ENDIF
        ENDIF
        IF(ICP   (IPHAS).GT.0) THEN
          IF(NOMVAR(IPPPRO(IPPROC(ICP   (IPHAS)))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPPRO(IPPROC(ICP   (IPHAS)))),'(A6,I2.2)')
     &                                                    'ChalSp',IPHAS
          ENDIF
        ENDIF
        IF(IESCAL(IESPRE,IPHAS).GT.0) THEN
          IPP = IPPPRO(IPPROC(IESTIM(IESPRE,IPHAS)))
          IF(NOMVAR(IPP) .EQ.' ') THEN
            WRITE(NOMVAR(IPP),'(A5,I1,I2.2)')
     &                             'EsPre',IESCAL(IESPRE,IPHAS),IPHAS
          ENDIF
        ENDIF
        IF(IESCAL(IESDER,IPHAS).GT.0) THEN
          IPP = IPPPRO(IPPROC(IESTIM(IESDER,IPHAS)))
          IF(NOMVAR(IPP) .EQ.' ') THEN
            WRITE(NOMVAR(IPP),'(A5,I1,I2.2)')
     &                             'EsDer',IESCAL(IESDER,IPHAS),IPHAS
          ENDIF
        ENDIF
        IF(IESCAL(IESCOR,IPHAS).GT.0) THEN
          IPP = IPPPRO(IPPROC(IESTIM(IESCOR,IPHAS)))
          IF(NOMVAR(IPP) .EQ.' ') THEN
            WRITE(NOMVAR(IPP),'(A5,I1,I2.2)')
     &                             'EsCor',IESCAL(IESCOR,IPHAS),IPHAS
          ENDIF
        ENDIF
        IF(IESCAL(IESTOT,IPHAS).GT.0) THEN
          IPP = IPPPRO(IPPROC(IESTIM(IESTOT,IPHAS)))
          IF(NOMVAR(IPP) .EQ.' ') THEN
            WRITE(NOMVAR(IPP),'(A5,I1,I2.2)')
     &                             'EsTot',IESCAL(IESTOT,IPHAS),IPHAS
          ENDIF
        ENDIF
C
      ENDDO
C
      DO JJ = 1, NSCAUS
        II = JJ
        IF(NOMVAR(IPPRTP(ISCA  (II   ))) .EQ.' ') THEN
          WRITE(NOMVAR(IPPRTP(ISCA  (II   ))),'(A5,I3.3)')'Scaus' ,II
        ENDIF
      ENDDO
      DO JJ = 1, NSCAPP
        II = ISCAPP(JJ)
        IF(NOMVAR(IPPRTP(ISCA  (II   ))) .EQ.' ') THEN
          WRITE(NOMVAR(IPPRTP(ISCA  (II   ))),'(A5,I3.3)')'Scapp' ,II
        ENDIF
      ENDDO
C
      IF(NBMOMT.GT.0) THEN
        DO IMOM = 1, NBMOMT
          IPP = IPPPRO(IPPROC(ICMOME(IMOM)))
          IF(NOMVAR(IPP) .EQ.' ') THEN
            WRITE(NOMVAR(IPP),'(A6,I2.2)')'MoyTps',IMOM
          ENDIF
        ENDDO
      ENDIF
C
      IF(NOMVAR(IPPDT                ) .EQ.' ') THEN
        WRITE(NOMVAR(IPPDT                ),'(A8     )')'PasDeTmp'
      ENDIF
C
      IF(NOMVAR(IPPTX                ) .EQ.' ') THEN
        WRITE(NOMVAR(IPPTX                ),'(A8     )')'Tx      '
      ENDIF
      IF(NOMVAR(IPPTY                ) .EQ.' ') THEN
        WRITE(NOMVAR(IPPTY                ),'(A8     )')'Ty      '
      ENDIF
      IF(NOMVAR(IPPTZ                ) .EQ.' ') THEN
        WRITE(NOMVAR(IPPTZ                ),'(A8     )')'Tz      '
      ENDIF
C
      IF (IALE.EQ.1) THEN
        IF(NOMVAR(IPPRTP(IUMA)) .EQ.' ') THEN
          WRITE(NOMVAR(IPPRTP(IUMA)),'(A8)')'VitmailX'
        ENDIF
        IF(NOMVAR(IPPRTP(IVMA)) .EQ.' ') THEN
          WRITE(NOMVAR(IPPRTP(IVMA)),'(A8)')'VitmailY'
        ENDIF
        IF(NOMVAR(IPPRTP(IWMA)) .EQ.' ') THEN
          WRITE(NOMVAR(IPPRTP(IWMA)),'(A8)')'VitmailZ'
        ENDIF
        IF (IORTVM.EQ.0) THEN
          IF(NOMVAR(IPPPRO(IPPROC(IVISMA(1)))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPPRO(IPPROC(IVISMA(1)))),'(A8)')'ViscMail'
          ENDIF
        ELSE
          IF(NOMVAR(IPPPRO(IPPROC(IVISMA(1)))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPPRO(IPPROC(IVISMA(1)))),'(A8)')'ViscMaiX'
          ENDIF
          IF(NOMVAR(IPPPRO(IPPROC(IVISMA(2)))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPPRO(IPPROC(IVISMA(2)))),'(A8)')'ViscMaiY'
          ENDIF
          IF(NOMVAR(IPPPRO(IPPROC(IVISMA(3)))) .EQ.' ') THEN
            WRITE(NOMVAR(IPPPRO(IPPROC(IVISMA(3)))),'(A8)')'ViscMaiZ'
          ENDIF
        ENDIF
      ENDIF
C
C ---> Sorties listing
C
      DO IPHAS = 1, NPHAS
        IPP = IPPPRO(IPPROC(IROM  (IPHAS)))
        IF(IROVAR(IPHAS).EQ.1.AND.ILISVR(IPP).EQ.-999) ILISVR(IPP) = 1
        IPP = IPPPRO(IPPROC(IVISCT(IPHAS)))
        IF( (ITURB(IPHAS).EQ.10 .OR. ITYTUR(IPHAS).EQ.2
     &       .OR. ITURB(IPHAS).EQ.50 .OR. ITURB(IPHAS).EQ.60)
     &       .AND.ILISVR(IPP).EQ.-999) ILISVR(IPP) = 1
        IPP = IPPPRO(IPPROC(ICOUR(IPHAS)))
        IF (ILISVR(IPP).EQ.-999 .OR. IDTVAR.LT.0) ILISVR(IPP) = 0
        IPP = IPPPRO(IPPROC(IFOUR(IPHAS)))
        IF (ILISVR(IPP).EQ.-999 .OR. IDTVAR.LT.0) ILISVR(IPP) = 0
        IF(IESCAL(IESPRE,IPHAS).GT.0) THEN
          IPP = IPPPRO(IPPROC(IESTIM(IESPRE,IPHAS)))
          IF(                     ILISVR(IPP).EQ.-999) ILISVR(IPP) = 1
        ENDIF
        IF(IESCAL(IESDER,IPHAS).GT.0) THEN
          IPP = IPPPRO(IPPROC(IESTIM(IESDER,IPHAS)))
          IF(                     ILISVR(IPP).EQ.-999) ILISVR(IPP) = 1
        ENDIF
        IF(IESCAL(IESCOR,IPHAS).GT.0) THEN
          IPP = IPPPRO(IPPROC(IESTIM(IESCOR,IPHAS)))
          IF(                     ILISVR(IPP).EQ.-999) ILISVR(IPP) = 1
        ENDIF
        IF(IESCAL(IESTOT,IPHAS).GT.0) THEN
          IPP = IPPPRO(IPPROC(IESTIM(IESTOT,IPHAS)))
          IF(                     ILISVR(IPP).EQ.-999) ILISVR(IPP) = 1
        ENDIF
      ENDDO
C
      IF(NBMOMT.GT.0) THEN
        DO IMOM = 1, NBMOMT
          IPP = IPPPRO(IPPROC(ICMOME(IMOM)))
          IF(ILISVR(IPP).EQ.-999) ILISVR(IPP) = 1
        ENDDO
      ENDIF
C
      IPP = IPPDT
      IF (IDTVAR.LT.0) ILISVR(IPP) = 0
      IF(ILISVR(IPP).EQ.-999) ILISVR(IPP)  = 1
      IF(IPUCOU.NE.1 .OR. IDTVAR.LT.0) THEN
        ILISVR(IPPTX) = 0
        ILISVR(IPPTY) = 0
        ILISVR(IPPTZ) = 0
      ENDIF
C
      DO II = 2, NVPPMX
        IF(ITRSVR(II).GE.1.AND.ILISVR(II).EQ.-999) THEN
          ILISVR(II) = 1
        ENDIF
      ENDDO
      DO II = 1, NVPPMX
        IF(ILISVR(II).EQ.-999) THEN
          ILISVR(II) = 0
        ENDIF
      ENDDO
C
C
C
C=======================================================================
C 2. POSITION DES VARIABLES DE numvar.h
C=======================================================================
C
C ---> Reperage des variables qui disposeront de deux types de CL
C
C     Fait dans varpos.
C     Si l'utilisateur y a touche ensuite, on risque l'incident.
C
C=======================================================================
C 3. OPTIONS DU CALCUL : TABLEAUX DE optcal.h
C=======================================================================
C
C ---> Schema en temps
C
C
      DO IPHAS = 1, NPHAS
C
C   -- Flux de masse
        IF(ABS(THETFL(IPHAS)+999.D0).GT.EPZERO) THEN
          WRITE(NFECRA,1001) IPHAS,ISTMPF(IPHAS)
          IOK = IOK + 1
        ELSEIF(ISTMPF(IPHAS).EQ.0) THEN
          THETFL(IPHAS) = 0.D0
        ELSEIF(ISTMPF(IPHAS).EQ.2) THEN
          THETFL(IPHAS) = 0.5D0
        ENDIF
C
C    -- Proprietes physiques
        IF(ABS(THETRO(IPHAS)+999.D0).GT.EPZERO) THEN
          WRITE(NFECRA,1011) IPHAS,'IROEXT',IROEXT(IPHAS),'THETRO'
          IOK = IOK + 1
        ELSEIF(IROEXT(IPHAS).EQ.0) THEN
          THETRO(IPHAS) = 0.0D0
        ELSEIF(IROEXT(IPHAS).EQ.1) THEN
          THETRO(IPHAS) = 0.5D0
        ELSEIF(IROEXT(IPHAS).EQ.2) THEN
          THETRO(IPHAS) = 1.D0
        ENDIF
        IF(ABS(THETVI(IPHAS)+999.D0).GT.EPZERO) THEN
          WRITE(NFECRA,1011) IPHAS,'IVIEXT',IVIEXT(IPHAS),'THETVI'
          IOK = IOK + 1
        ELSEIF(IVIEXT(IPHAS).EQ.0) THEN
          THETVI(IPHAS) = 0.0D0
        ELSEIF(IVIEXT(IPHAS).EQ.1) THEN
          THETVI(IPHAS) = 0.5D0
        ELSEIF(IVIEXT(IPHAS).EQ.2) THEN
          THETVI(IPHAS) = 1.D0
        ENDIF
        IF(ABS(THETCP(IPHAS)+999.D0).GT.EPZERO) THEN
          WRITE(NFECRA,1011) IPHAS,'ICPEXT',ICPEXT(IPHAS),'THETCP'
          IOK = IOK + 1
        ELSEIF(ICPEXT(IPHAS).EQ.0) THEN
          THETCP(IPHAS) = 0.0D0
        ELSEIF(ICPEXT(IPHAS).EQ.1) THEN
          THETCP(IPHAS) = 0.5D0
        ELSEIF(ICPEXT(IPHAS).EQ.2) THEN
          THETCP(IPHAS) = 1.D0
        ENDIF
C
C    -- Termes sources NS
        IF(ABS(THETSN(IPHAS)+999.D0).GT.EPZERO) THEN
          WRITE(NFECRA,1011) IPHAS,'ISNO2T',ISNO2T(IPHAS),'THETSN'
          IOK = IOK + 1
        ELSEIF(ISNO2T(IPHAS).EQ.1) THEN
          THETSN(IPHAS) = 0.5D0
        ELSEIF(ISNO2T(IPHAS).EQ.2) THEN
          THETSN(IPHAS) = 1.D0
        ELSEIF(ISNO2T(IPHAS).EQ.0) THEN
          THETSN(IPHAS) = 0.D0
        ENDIF
C
C    -- Termes sources grandeurs turbulentes
        IF(ABS(THETST(IPHAS)+999.D0).GT.EPZERO) THEN
          WRITE(NFECRA,1011) IPHAS,'ISTO2T',ISTO2T(IPHAS),'THETST'
          IOK = IOK + 1
        ELSEIF(ISTO2T(IPHAS).EQ.1) THEN
          THETST(IPHAS) = 0.5D0
        ELSEIF(ISTO2T(IPHAS).EQ.2) THEN
          THETST(IPHAS) = 1.D0
        ELSEIF(ISTO2T(IPHAS).EQ.0) THEN
          THETST(IPHAS) = 0.D0
        ENDIF
C
      ENDDO
C
      DO ISCAL = 1, NSCAL
C    -- Termes sources des scalaires
        IF(ABS(THETSS(ISCAL)+999.D0).GT.EPZERO) THEN
          WRITE(NFECRA,1021) ISCAL,'ISSO2T',ISSO2T(ISCAL),'THETSS'
          IOK = IOK + 1
        ELSEIF(ISSO2T(ISCAL).EQ.1) THEN
          THETSS(ISCAL) = 0.5D0
        ELSEIF(ISSO2T(ISCAL).EQ.2) THEN
          THETSS(ISCAL) = 1.D0
        ELSEIF(ISSO2T(ISCAL).EQ.0) THEN
          THETSS(ISCAL) = 0.D0
        ENDIF
C    -- Diffusivite des scalaires
        IF(ABS(THETVS(ISCAL)+999.D0).GT.EPZERO) THEN
          WRITE(NFECRA,1021) ISCAL,'IVSEXT',IVSEXT(ISCAL),'THETVS'
          IOK = IOK + 1
        ELSEIF(IVSEXT(ISCAL).EQ.0) THEN
          THETVS(ISCAL) = 0.D0
        ELSEIF(IVSEXT(ISCAL).EQ.1) THEN
          THETVS(ISCAL) = 0.5D0
        ELSEIF(IVSEXT(ISCAL).EQ.2) THEN
          THETVS(ISCAL) = 1.D0
        ENDIF
      ENDDO
C
C     Ici on interdit que l'utilisateur fixe lui meme THETAV, par securite
C       mais on pourrait le laisser faire
C       (enlever le IOK, modifier le message et les tests dans verini)
      DO IPHAS = 1, NPHAS
C
C     Vitesse pression (la pression est prise sans interp)
        IF(ABS(THETAV(IU (IPHAS))+999.D0).GT.EPZERO.OR.
     &     ABS(THETAV(IV (IPHAS))+999.D0).GT.EPZERO.OR.
     &     ABS(THETAV(IW (IPHAS))+999.D0).GT.EPZERO.OR.
     &     ABS(THETAV(IPR(IPHAS))+999.D0).GT.EPZERO) THEN
          WRITE(NFECRA,1031) IPHAS,'VITESSE-PRESSION ','THETAV'
          IOK = IOK + 1
        ELSEIF(ISCHTP(IPHAS).EQ.1) THEN
          THETAV(IU (IPHAS)) = 1.D0
          THETAV(IV (IPHAS)) = 1.D0
          THETAV(IW (IPHAS)) = 1.D0
          THETAV(IPR(IPHAS)) = 1.D0
        ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
          THETAV(IU (IPHAS)) = 0.5D0
          THETAV(IV (IPHAS)) = 0.5D0
          THETAV(IW (IPHAS)) = 0.5D0
          THETAV(IPR(IPHAS)) = 1.D0
        ENDIF
C
C     Turbulence (en k-eps : ordre 1)
        IF(ITYTUR(IPHAS).EQ.2) THEN
          IF(ABS(THETAV(IK (IPHAS))+999.D0).GT.EPZERO.OR.
     &       ABS(THETAV(IEP(IPHAS))+999.D0).GT.EPZERO) THEN
            WRITE(NFECRA,1031) IPHAS,'VARIABLES   K-EPS','THETAV'
            IOK = IOK + 1
          ELSEIF(ISCHTP(IPHAS).EQ.1) THEN
            THETAV(IK (IPHAS)) = 1.D0
            THETAV(IEP(IPHAS)) = 1.D0
          ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
C     pour le moment, on ne peut pas passer par ici (cf varpos)
            THETAV(IK (IPHAS)) = 0.5D0
            THETAV(IEP(IPHAS)) = 0.5D0
          ENDIF
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
          IF(ABS(THETAV(IR11(IPHAS))+999.D0).GT.EPZERO.OR.
     &       ABS(THETAV(IR22(IPHAS))+999.D0).GT.EPZERO.OR.
     &       ABS(THETAV(IR33(IPHAS))+999.D0).GT.EPZERO.OR.
     &       ABS(THETAV(IR12(IPHAS))+999.D0).GT.EPZERO.OR.
     &       ABS(THETAV(IR13(IPHAS))+999.D0).GT.EPZERO.OR.
     &       ABS(THETAV(IR23(IPHAS))+999.D0).GT.EPZERO.OR.
     &       ABS(THETAV(IEP (IPHAS))+999.D0).GT.EPZERO) THEN
            WRITE(NFECRA,1031) IPHAS,'VARIABLES  RIJ-EP','THETAV'
            IOK = IOK + 1
          ELSEIF(ISCHTP(IPHAS).EQ.1) THEN
            THETAV(IR11(IPHAS)) = 1.D0
            THETAV(IR22(IPHAS)) = 1.D0
            THETAV(IR33(IPHAS)) = 1.D0
            THETAV(IR12(IPHAS)) = 1.D0
            THETAV(IR13(IPHAS)) = 1.D0
            THETAV(IR23(IPHAS)) = 1.D0
            THETAV(IEP (IPHAS)) = 1.D0
          ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
            THETAV(IR11(IPHAS)) = 0.5D0
            THETAV(IR22(IPHAS)) = 0.5D0
            THETAV(IR33(IPHAS)) = 0.5D0
            THETAV(IR12(IPHAS)) = 0.5D0
            THETAV(IR13(IPHAS)) = 0.5D0
            THETAV(IR23(IPHAS)) = 0.5D0
            THETAV(IEP (IPHAS)) = 0.5D0
          ENDIF
        ELSEIF(ITURB(IPHAS).EQ.50) THEN
          IF(ABS(THETAV(IK  (IPHAS))+999.D0).GT.EPZERO.OR.
     &       ABS(THETAV(IEP (IPHAS))+999.D0).GT.EPZERO.OR.
     &       ABS(THETAV(IPHI(IPHAS))+999.D0).GT.EPZERO.OR.
     &       ABS(THETAV(IFB (IPHAS))+999.D0).GT.EPZERO) THEN
            WRITE(NFECRA,1031) IPHAS,'VARIABLES     V2F','THETAV'
            IOK = IOK + 1
          ELSEIF(ISCHTP(IPHAS).EQ.1) THEN
            THETAV(IK  (IPHAS)) = 1.D0
            THETAV(IEP (IPHAS)) = 1.D0
            THETAV(IPHI(IPHAS)) = 1.D0
            THETAV(IFB (IPHAS)) = 1.D0
          ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
C     pour le moment, on ne peut pas passer par ici (cf varpos)
            THETAV(IK  (IPHAS)) = 0.5D0
            THETAV(IEP (IPHAS)) = 0.5D0
            THETAV(IPHI(IPHAS)) = 0.5D0
            THETAV(IFB (IPHAS)) = 0.5D0
          ENDIF
        ELSEIF(ITURB(IPHAS).EQ.60) THEN
          IF(ABS(THETAV(IK  (IPHAS))+999.D0).GT.EPZERO.OR.
     &       ABS(THETAV(IOMG(IPHAS))+999.D0).GT.EPZERO ) THEN
            WRITE(NFECRA,1031) IPHAS,'VARIABLES K-OMEGA','THETAV'
            IOK = IOK + 1
          ELSEIF(ISCHTP(IPHAS).EQ.1) THEN
            THETAV(IK  (IPHAS)) = 1.D0
            THETAV(IOMG(IPHAS)) = 1.D0
          ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
C     pour le moment, on ne peut pas passer par ici (cf varpos)
            THETAV(IK  (IPHAS)) = 0.5D0
            THETAV(IOMG(IPHAS)) = 0.5D0
          ENDIF
        ENDIF
C
      ENDDO
C
C     Scalaires
      DO ISCAL = 1, NSCAL
        IPHAS = IPHSCA(ISCAL)
        IVAR  = ISCA(ISCAL)
        IF(ABS(THETAV(IVAR)+999.D0).GT.EPZERO) THEN
          WRITE(NFECRA,1041) IPHAS,'SCALAIRE',ISCAL,'THETAV'
          IOK = IOK + 1
        ELSEIF(ISCHTP(IPHAS).EQ.1) THEN
          THETAV(IVAR) = 1.D0
        ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
          THETAV(IVAR) = 0.5D0
        ENDIF
      ENDDO
C
C     Vitesse de maillage en ALE
      IF (IALE.EQ.1) THEN
        IPHAS = 1
        IF(ABS(THETAV(IUMA)+999.D0).GT.EPZERO.OR.
     &     ABS(THETAV(IVMA)+999.D0).GT.EPZERO.OR.
     &     ABS(THETAV(IWMA)+999.D0).GT.EPZERO) THEN
          WRITE(NFECRA,1032) 'THETAV'
          IOK = IOK + 1
        ELSEIF(ISCHTP(IPHAS).EQ.1) THEN
          THETAV(IUMA) = 1.D0
          THETAV(IVMA) = 1.D0
          THETAV(IWMA) = 1.D0
        ELSEIF(ISCHTP(IPHAS).EQ.2) THEN
C     pour le moment, on ne peut pas passer par ici (cf varpos)
          THETAV(IUMA) = 0.5D0
          THETAV(IVMA) = 0.5D0
          THETAV(IWMA) = 0.5D0
        ENDIF
      ENDIF
C
C ---> ISSTPC
C        Si l'utilisateur n'a rien specifie pour le test de pente (=-999),
C        On impose 1 (ie sans) pour la vitesse en LES
C                  0 (ie avec) sinon
C
      DO IPHAS = 1, NPHAS
        IF(ITYTUR(IPHAS).EQ.4) THEN
          II = IU(IPHAS)
          IF(ISSTPC(II).EQ.-999) ISSTPC(II) = 1
          II = IV(IPHAS)
          IF(ISSTPC(II).EQ.-999) ISSTPC(II) = 1
          II = IW(IPHAS)
          IF(ISSTPC(II).EQ.-999) ISSTPC(II) = 1
          DO JJ = 1, NSCAL
            II = ISCA(JJ)
            IF(ISSTPC(II).EQ.-999) ISSTPC(II) = 0
          ENDDO
        ENDIF
      ENDDO
C
      DO II = 1, NVARMX
        IF (ISSTPC(II).EQ.-999) THEN
          ISSTPC(II) = 0
        ENDIF
      ENDDO
C
C ---> BLENCV
C        Si l'utilisateur n'a rien specifie pour le schema convectif
C                  1 (ie centre) pour les vitesses et
C                                      les scalaires utilisateurs
C                  0 (ie upwind pur) pour le reste
C   (en particulier, en L.E.S. toutes les variables sont donc en centre)
C
      DO IPHAS = 1, NPHAS
        II = IU(IPHAS)
        IF(ABS(BLENCV(II)+999.D0).LT.EPZERO) BLENCV(II) = 1.D0
        II = IV(IPHAS)
        IF(ABS(BLENCV(II)+999.D0).LT.EPZERO) BLENCV(II) = 1.D0
        II = IW(IPHAS)
        IF(ABS(BLENCV(II)+999.D0).LT.EPZERO) BLENCV(II) = 1.D0
      ENDDO
      DO JJ = 1, NSCAUS
        II = ISCA(JJ)
        IF(ABS(BLENCV(II)+999.D0).LT.EPZERO) BLENCV(II) = 1.D0
      ENDDO
C
      DO II = 1, NVARMX
        IF (ABS(BLENCV(II)+999.D0).LT.EPZERO) THEN
          BLENCV(II) = 0.D0
        ENDIF
      ENDDO
C
C
C ---> NSWRSM, EPSRSM ET EPSILO
C        Si l'utilisateur n'a rien specifie  (NSWRSM=-999),
C        On impose
C           a l'ordre 1 :
C                  2  pour la pression
C                  1  pour les autres variables
C                  on initialise EPSRSM a 1.D-8
C                  on initialise EPSILO a 1.D-8
C           a l'ordre 2 :
C                  5  pour la pression
C                  10 pour les autres variables
C                  on initialise EPSRSM a 1.D-5
C                  on initialise EPSILO a 1.D-5
C     Attention aux tests dans verini
C
      DO IPHAS = 1, NPHAS
        IF(ISCHTP(IPHAS).EQ.2) THEN
          II = IPR(IPHAS)
          IF(NSWRSM(II).EQ.-999) NSWRSM(II) = 5
          IF(ABS(EPSRSM(II)+999.D0).LT.EPZERO) EPSRSM(II) = 1.D-5
          IF(ABS(EPSILO(II)+999.D0).LT.EPZERO) EPSILO(II) = 1.D-5
          II = IU(IPHAS)
          IF(NSWRSM(II).EQ.-999) NSWRSM(II) = 10
          IF(ABS(EPSRSM(II)+999.D0).LT.EPZERO) EPSRSM(II) = 1.D-5
          IF(ABS(EPSILO(II)+999.D0).LT.EPZERO) EPSILO(II) = 1.D-5
          II = IV(IPHAS)
          IF(NSWRSM(II).EQ.-999) NSWRSM(II) = 10
          IF(ABS(EPSRSM(II)+999.D0).LT.EPZERO) EPSRSM(II) = 1.D-5
          IF(ABS(EPSILO(II)+999.D0).LT.EPZERO) EPSILO(II) = 1.D-5
          II = IW(IPHAS)
          IF(NSWRSM(II).EQ.-999) NSWRSM(II) = 10
          IF(ABS(EPSRSM(II)+999.D0).LT.EPZERO) EPSRSM(II) = 1.D-5
          IF(ABS(EPSILO(II)+999.D0).LT.EPZERO) EPSILO(II) = 1.D-5
          DO JJ = 1, NSCAL
            II = ISCA(JJ)
            IF(NSWRSM(II).EQ.-999) NSWRSM(II) = 10
            IF(ABS(EPSRSM(II)+999.D0).LT.EPZERO) EPSRSM(II) = 1.D-5
            IF(ABS(EPSILO(II)+999.D0).LT.EPZERO) EPSILO(II) = 1.D-5
          ENDDO
        ENDIF
        II = IPR(IPHAS)
        IF(NSWRSM(II).EQ.-999) NSWRSM(II) = 2
      ENDDO
C
      DO II = 1, NVARMX
        IF (NSWRSM(II).EQ.-999) NSWRSM(II) = 1
        IF(ABS(EPSRSM(II)+999.D0).LT.EPZERO) EPSRSM(II) = 1.D-8
        IF(ABS(EPSILO(II)+999.D0).LT.EPZERO) EPSILO(II) = 1.D-8
      ENDDO
C
C ---> ANOMAX
C        Si l'utilisateur n'a rien specifie pour l'angle de non
C          orthogonalite pour la selection du voisinage etendu,
C        On impose pi/4 avec gradmc + support etendu
C               et rien sinon
C
      IF (IMRGRA.EQ.3) THEN
        IF (ANOMAX.LE.-GRAND) THEN
          ANOMAX = PI*0.25D0
        ENDIF
      ENDIF
C
C ---> IMLIGR
C        Si l'utilisateur n'a rien specifie pour la limitation des
C          gradients (=-999),
C        On impose -1 avec gradrc (pas de limitation)
C               et  1 avec gradmc (limitation)
C
      IF (IMRGRA.EQ.0.OR.IMRGRA.EQ.4) THEN
        DO II = 1, NVARMX
          IF (IMLIGR(II).EQ.-999) THEN
            IMLIGR(II) = -1
          ENDIF
        ENDDO
      ELSEIF (IMRGRA.EQ.1.OR.IMRGRA.EQ.2.OR.IMRGRA.EQ.3) THEN
        DO II = 1, NVARMX
          IF (IMLIGR(II).EQ.-999) THEN
            IMLIGR(II) = 1
          ENDIF
        ENDDO
      ENDIF
C
C ---> DTMIN DTMAX CDTVAR
C
C
      IF(DTMIN.LE.-GRAND) THEN
        DTMIN = 0.1D0*DTREF
      ENDIF
      IF(DTMAX.LE.-GRAND) THEN
        DTMAX = 1.0D3*DTREF
      ENDIF
C
      DO IPHAS = 1, NPHAS
C
C     Ici, ce n'est pas grave pour le moment,
C      etant entendu que ces coefs ne servent pas
C      s'ils servaient, attention dans le cas a plusieurs phases avec
C      une seule pression : celle ci prend le coef de la derniere phase
        CDTVAR(IV (IPHAS)) = CDTVAR(IU(IPHAS))
        CDTVAR(IW (IPHAS)) = CDTVAR(IU(IPHAS))
        CDTVAR(IPR(IPHAS)) = CDTVAR(IU(IPHAS))
C
        IF(ITYTUR(IPHAS).EQ.2) THEN
          CDTVAR(IEP (IPHAS)) = CDTVAR(IK  (IPHAS))
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
          CDTVAR(IR22(IPHAS)) = CDTVAR(IR11(IPHAS))
          CDTVAR(IR33(IPHAS)) = CDTVAR(IR11(IPHAS))
          CDTVAR(IR12(IPHAS)) = CDTVAR(IR11(IPHAS))
          CDTVAR(IR13(IPHAS)) = CDTVAR(IR11(IPHAS))
          CDTVAR(IR23(IPHAS)) = CDTVAR(IR11(IPHAS))
          CDTVAR(IEP (IPHAS)) = CDTVAR(IR11(IPHAS))
        ELSEIF(ITURB(IPHAS).EQ.50) THEN
          CDTVAR(IEP (IPHAS)) = CDTVAR(IK  (IPHAS))

          CDTVAR(IPHI(IPHAS)) = CDTVAR(IK  (IPHAS))
c     CDTVAR(IFB) est en fait inutile car pas de temps dans l'eq de f_barre
          CDTVAR(IFB (IPHAS)) = CDTVAR(IK  (IPHAS))
        ELSEIF(ITURB(IPHAS).EQ.60) THEN
          CDTVAR(IOMG(IPHAS)) = CDTVAR(IK  (IPHAS))
        ENDIF
C
      ENDDO
C
C ---> IDEUCH, YPLULI
C      En laminaire, longueur de melange et LES, une echelle de vitesse
C      Sinon, 2 echelles, sauf si l'utilisateur choisit 1 echelle.
C      On a initialise IDEUCH a -999 pour voir si l'utilisateur essaye
C        de choisir deux echelles quand ce n'est pas possible et le
C        prevenir dans la section verification.
C
      DO IPHAS = 1, NPHAS
        IF(IDEUCH(IPHAS).EQ.-999) THEN
          IF(ITURB(IPHAS).EQ. 0.OR.
     &       ITURB(IPHAS).EQ.10.OR.
     &       ITYTUR(IPHAS).EQ.4 ) THEN
            IDEUCH(IPHAS) = 0
          ELSE
            IDEUCH(IPHAS) = 1
          ENDIF
        ENDIF
C     Pour YPLULI, 1/XKAPPA est la valeur qui assure la continuite de la derivee entre
C     la zone lineaire et la zone logarithmique. Dans le cas des lois de paroi invariantes,
C     on utilise la valeur de continuite du profil de vitesse, 10,88.
C     Pour la L.E.S., on remet 10,88, afin d'eviter des clic/clac quand on est a la limite
C     (en modele a une echelle en effet, YPLULI=1/XKAPPA ne permet pas forcement de calculer
C      u* de maniere totalement satisfaisante)
        IF (YPLULI(IPHAS).LT.-GRAND) THEN
          IF (IDEUCH(IPHAS).EQ.2 .OR. ITYTUR(IPHAS).EQ.4 ) THEN
            YPLULI(IPHAS) = 10.88D0
          ELSE
            YPLULI(IPHAS) = 1.D0/XKAPPA
          ENDIF
        ENDIF
      ENDDO
C
C
C ---> Van Driest
      DO IPHAS = 1, NPHAS
        IF(IDRIES(IPHAS).EQ.-1) THEN
C   On met 1 en supposant qu'en periodicite ou parallele on utilise le
C     mode de calcul de la distance a la paroi qui les prend en charge
C     (ICDPAR=+/-1, valeur par defaut)
          IF(ITURB(IPHAS).EQ.40) THEN
            IDRIES(IPHAS) = 1
          ELSEIF(ITURB(IPHAS).EQ.41) THEN
            IDRIES(IPHAS) = 0
          ENDIF
        ENDIF
      ENDDO
C
C
C ---> ICPSYR
C      Si l'utilisateur n'a pas modifie ICPSYR, on prend par defaut :
C        s'il n y a pas de couplage
C          0 pour tous les scalaires
C        sinon
C          1 pour le scalaire ISCALT s'il existe
C          0 pour les autres
C      Les modifs adequates devront etre ajoutees pour les physiques
C        particulieres
C      Les tests de coherence seront faits dans verini.
C
      IF(NSCAL.GT.0) THEN
C
C     On regarde s'il y a du couplage
C
        CALL NBCSYR (NBCCOU)
C       ===========
C
C     S'il y a du couplage
        IF (NBCCOU .NE. 0) THEN
C
C       On compte le nombre de scalaires couples
          NSCACP = 0
          DO ISCAL = 1, NSCAL
            IF(ICPSYR(ISCAL).EQ.1) THEN
              NSCACP = NSCACP + 1
            ENDIF
          ENDDO
C
C       Si l'utilisateur n'a pas couple de scalaire,
          IF(NSCACP.EQ.0) THEN
C
C         On couple le scalaire temperature de la premiere phase
            DO IPHAS = 1, NPHAS
              IF(ISCALT(IPHAS).GT.0.AND.ISCALT(IPHAS).LE.NSCAL) THEN
                ICPSYR(ISCALT(IPHAS)) = 1
                GOTO 100
              ENDIF
            ENDDO
 100        CONTINUE
C
          ENDIF
C
        ENDIF
C
C     Pour tous les autres scalaires, non renseignes pas l'utilisateur
C       on ne couple pas
        DO ISCAL = 1, NSCAMX
          IF(ICPSYR(ISCAL).EQ.-999) THEN
            ICPSYR(ISCAL) = 0
          ENDIF
        ENDDO
C
      ENDIF
C
C
C ---> ISCSTH
C      Si l'utilisateur n'a pas modifie ISCSTH, on prend par defaut :
C        scalaire passif  pour les scalaires autres que ISCALT
C      Les modifs adequates devront etre ajoutees pour les physiques
C        particulieres
C      Noter en outre que, par defaut, si on choisit temperature
C        elle est en K (ceci n'est utile que pour le rayonnement et les pp)
C
C         =-10: non renseigne
C         =-1 : temperature en C
C         = 0 : passif
C         = 1 : temperature en K
C         = 2 : enthalpie
C
C
C
      IF(NSCAL.GT.0) THEN
        DO II = 1, NSCAL
          IF(ISCSTH(II).EQ.-10)THEN
            IPHAS = IPHSCA(II)
            IF(II.NE.ISCALT(IPHAS)) THEN
              ISCSTH(II) = 0
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C ---> ICALHY
C      Calcul de la pression hydrostatique en sortie pour les conditions de
C        Dirichlet sur la pression. Se deduit de IPHYDR et de la valeur de
C        la gravite (test assez arbitraire sur la norme).
C      ICALHY est initialise a -1 (l'utilisateur peut avoir force
C        0 ou 1 et dans ce cas, on ne retouche pas)
C
      IF(ICALHY.NE.-1.AND.ICALHY.NE.0.AND.ICALHY.NE.1) THEN
        WRITE(NFECRA,1061)ICALHY
        IOK = IOK + 1
      ENDIF
C
C
C ---> IDGMOM
C      Calcul du degre des moments
C
      DO IMOM = 1, NBMOMX
        IDGMOM(IMOM) = 0
      ENDDO
      DO IMOM = 1, NBMOMT
        DO II = 1, NDGMOX
          IF(IDFMOM(II,IMOM).NE.0) IDGMOM(IMOM) = IDGMOM(IMOM) + 1
        ENDDO
      ENDDO
C
C ---> ICDPAR
C      Calcul de la distance a la paroi. En standard, on met ICDPAR a -1, au cas
C      ou les faces de bord auraient change de type d'un calcul a l'autre. En k-omega,
C      il faut la distance a la paroi pour une suite propre, donc on initialise a 1 et
C      on avertit (dans verini).
      IKW = 0
      DO IPHAS = 1, NPHAS
        IF (ITURB(IPHAS).EQ.60) IKW = 1
      ENDDO
      IF (ICDPAR.EQ.-999) THEN
        ICDPAR = -1
        IF (IKW.EQ.1) ICDPAR = 1
        IF (ISUITE.EQ.1 .AND. IKW.EQ.1) WRITE(NFECRA,2000)
      ENDIF
      IF (ICDPAR.EQ.-1 .AND. IKW.EQ.1 .AND. ISUITE.EQ.1)
     &     WRITE(NFECRA,2001)
C
C ---> INEEDY, IMLIGY
C      Calcul de la distance a la paroi
C       (une seule phase ...)
C
      INEEDY = 0
      DO IPHAS = 1, NPHAS
        IF((ITURB(IPHAS).EQ.30.AND.IRIJEC(IPHAS).EQ.1).OR.
     &     (ITYTUR(IPHAS).EQ.4.AND.IDRIES(IPHAS).EQ.1).OR.
     &      ITURB(IPHAS).EQ.60 ) THEN
          INEEDY = 1
        ENDIF
      ENDDO
C
      IF (IMRGRA.EQ.0 .OR. IMRGRA.EQ.4) THEN
        IF (IMLIGY.EQ.-999) THEN
          IMLIGY = -1
        ENDIF
      ELSEIF (IMRGRA.EQ.1.OR.IMRGRA.EQ.2.OR.IMRGRA.EQ.3) THEN
        IF (IMLIGY.EQ.-999) THEN
          IMLIGY = 1
        ENDIF
      ENDIF
C
C     Warning : non initialise => comme la vitesse
      IF(IWARNY.EQ.-999) THEN
        IWARNY = IWARNI(IU(1))
      ENDIF
C
C
C ---> IKECOU
C     En k-eps prod lin, v2f ou k-omega, on met IKECOU a 0 par defaut,
C     sinon on le laisse a 1
C     Dans verini on bloquera le v2f et le k-eps prod lin si IKECOU.NE.0
C     On bloquera aussi le stationnaire si IKECOU.NE.0
      DO IPHAS = 1, NPHAS
        IF (IKECOU(IPHAS).EQ.-999) THEN
          IF (IDTVAR.LT.0) THEN
            IKECOU(IPHAS) = 0
          ELSE IF (ITURB(IPHAS).EQ.21 .OR. ITURB(IPHAS).EQ.50
     &        .OR. ITURB(IPHAS).EQ.60 ) THEN
            IKECOU(IPHAS) = 0
          ELSE
            IKECOU(IPHAS) = 1
          ENDIF
        ENDIF
      ENDDO
C
C ---> RELAXV
      IF (IDTVAR.LT.0) THEN
        RELXSP = 1.D0-RELXST
        IF (RELXSP.LE.EPZERO) RELXSP = RELXST
        DO IPHAS = 1, NPHAS
          IF (ABS(RELAXV(IPR(IPHAS))+999.D0).LE.EPZERO)
     &         RELAXV(IPR(IPHAS)) = RELXSP
        ENDDO
        DO II = 1, NVARMX
          IF (ABS(RELAXV(II)+999.D0).LE.EPZERO) RELAXV(II) = RELXST
        ENDDO
      ELSE
        DO IPHAS = 1, NPHAS
          IF ( IKECOU(IPHAS).EQ.0) THEN
            IF (ITYTUR(IPHAS).EQ.2 .OR. ITYTUR(IPHAS).EQ.5) THEN
              IF (ABS(RELAXV(IK(IPHAS))+999.D0).LT.EPZERO)
     &             RELAXV(IK(IPHAS)) = 0.7D0
              IF (ABS(RELAXV(IEP(IPHAS))+999.D0).LT.EPZERO)
     &             RELAXV(IEP(IPHAS)) = 0.7D0
            ELSE IF (ITYTUR(IPHAS).EQ.6) THEN
              IF (ABS(RELAXV(IK(IPHAS))+999.D0).LT.EPZERO)
     &             RELAXV(IK(IPHAS)) = 0.7D0
              IF (ABS(RELAXV(IOMG(IPHAS))+999.D0).LT.EPZERO)
     &             RELAXV(IOMG(IPHAS)) = 0.7D0
            ENDIF
          ENDIF
          IF (ABS(RELAXV(IPR(IPHAS))+999.D0).LT.EPZERO)
     &             RELAXV(IPR(IPHAS)) = 1.D0
        ENDDO
      ENDIF
C
C ---> SPECIFIQUE STATIONNAIRE
      IF (IDTVAR.LT.0) THEN
        DTREF = 1.D0
        DTMIN = 1.D0
        DTMAX = 1.D0
        DO II = 1, NVARMX
          ISTAT(II) = 0
        ENDDO
        DO IPHAS = 1, NPHAS
          ARAK(IPHAS) = ARAK(IPHAS)/MAX(RELAXV(IU(IPHAS)),EPZERO)
        ENDDO
      ENDIF
C
C ---> INEEDF
C     Si on a demande un posttraitement des efforts aux bords, on
C     les calcule !
      IF (MOD(IPSTDV,IPSTFO).EQ.0) THEN
        INEEDF = 1
      ENDIF
C     Si on est en ALE, par defaut on calcule les efforts aux bords
C     (un test eventuel sur la presence de structures viendrait trop
C     tard)
      IF (IALE.EQ.1) INEEDF = 1
C
C=======================================================================
C 4. TABLEAUX DE cstphy.h
C=======================================================================
C
C ---> Constantes
C    Ca fait un calcul en double, mais si qqn a bouge cmu, apow, bpow,
C     ca servira.
C
      CPOW    = APOW**(2.D0/(1.D0-BPOW))
      DPOW    = 1.D0/(1.D0+BPOW)
      CMU025 = CMU**0.25D0
C
C ---> ICLVFL
C      Si l'utilisateur n'a pas modifie ICLVFL, on prend par defaut :
C        0 pour les variances
C      Les modifs adequates devront etre ajoutees pour les physiques
C        particulieres
C
      DO ISCAL = 1, NSCAL
        IF(ISCAVR(ISCAL).GT.0) THEN
          IF(ICLVFL(ISCAL).EQ.-1) THEN
            ICLVFL(ISCAL) = 0
          ENDIF
        ENDIF
      ENDDO
C
C
C ---> VISLS0 (IPHSCA, IVISLS ont ete verifies dans varpos)
C      Pour les variances de fluctuations, les valeurs du tableau
C        precedent ne doivent pas avoir ete modifiees par l'utilisateur
C        Elles sont prises egales aux valeurs correspondantes pour le
C        scalaire associe.
C
      IF(NSCAL.GT.0) THEN
        DO II = 1, NSCAL
          ISCAL = ISCAVR(II)
          IF(ISCAL.GT.0.AND.ISCAL.LE.NSCAL)THEN
            IF(VISLS0(II).LT.-GRAND) THEN
              VISLS0(II) = VISLS0(ISCAL)
            ELSE
              WRITE(NFECRA,1071)II,
     &          II,ISCAL,II,ISCAL,
     &          II,II,ISCAL,VISLS0(ISCAL)
              IOK = IOK + 1
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C ---> XYZP0 : reference pour la pression hydrostatique
C      On considere que l'utilisateur a specifie la reference
C      a partir du moment ou il a specifie une coordonnee.
C      Pour les coordonnees non specifiees, on met 0.
C
      DO IPHAS = 1, NPHAS
        DO II = 1, 3
          IF (XYZP0(II,IPHAS).GT.-0.5D0*RINFIN) THEN
            IXYZP0(IPHAS) = 1
          ELSE
            XYZP0(II,IPHAS) = 0.D0
          ENDIF
        ENDDO
      ENDDO
C
C=======================================================================
C 5. ELEMENTS DE albase.h
C=======================================================================
C
      IF (IALE.EQ.1) THEN
        IF (ISUITE.EQ.0 .AND. ITALIN.EQ.-999 ) ITALIN = 1
      ELSE
        ITALIN = 0
      ENDIF
C
C=======================================================================
C 6. COEFFICIENTS DE alstru.h
C=======================================================================
C
      IF (BETNMK.LT.-0.5D0*GRAND)
     &     BETNMK = (1.D0-ALPNMK)**2/4.D0
      IF (GAMNMK.LT.-0.5D0*GRAND)
     &     GAMNMK = (1.D0-2.D0*ALPNMK)/2.D0
      IF (AEXXST.LT.-0.5D0*GRAND) AEXXST = 0.5D0
      IF (BEXXST.LT.-0.5D0*GRAND) BEXXST = 0.0D0
      IF (CFOPRE.LT.-0.5D0*GRAND) CFOPRE = 2.0D0
C
C=======================================================================
C 6. STOP SI PB
C=======================================================================
C
      IF(IOK.NE.0) THEN
        CALL CSEXIT (1)
      ENDIF
C
C
 1001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',   I10,' ISTMPF = ',   I10                       ,/,
     &'@    THETFL SERA INITIALISE AUTOMATIQUEMENT.                 ',/,
     &'@    NE PAS LE MODIFIER DANS usini1.                         ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',   I10,' ',A6,' = ',   I10                       ,/,
     &'@    ',A6,' SERA INITIALISE AUTOMATIQUEMENT.                 ',/,
     &'@    NE PAS LE MODIFIER DANS usini1.                         ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1021 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',   I10,' ',A6,' = ',   I10                    ,/,
     &'@    ',A6,' SERA INITIALISE AUTOMATIQUEMENT.                 ',/,
     &'@    NE PAS LE MODIFIER DANS usini1.                         ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1031 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',   I10,' ',A17                                   ,/,
     &'@    ',A6,' SERA INITIALISE AUTOMATIQUEMENT.                 ',/,
     &'@    NE PAS LE MODIFIER DANS usini1.                         ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1032 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    VITESSE DE MAILLAGE EN ALE                              ',/,
     &'@    ',A6,' SERA INITIALISE AUTOMATIQUEMENT.                 ',/,
     &'@    NE PAS LE MODIFIER DANS usini1.                         ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1041 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',   I10,' ',A8,' ',I10                            ,/,
     &'@    ',A6,' SERA INITIALISE AUTOMATIQUEMENT.                 ',/,
     &'@    NE PAS LE MODIFIER DANS usini1.                         ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1061 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    ICALHY doit etre un entier egal a 0 ou 1                ',/,
     &'@                                                            ',/,
     &'@  Il vaut ici ',I10                                          ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1071 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE ',I10   ,' NE PAS MODIFIER LA DIFFUSIVITE      ',/,
     &'@                                                            ',/,
     &'@  Le scalaire ',I10   ,' represente la variance des         ',/,
     &'@    fluctuations du scalaire ',I10                           ,/,
     &'@                             (ISCAVR(',I10   ,') = ',I10     ,/,
     &'@  La diffusivite VISLS0(',I10   ,') du scalaire ',I10        ,/,
     &'@    ne doit pas etre renseignee :                           ',/,
     &'@    elle sera automatiquement prise egale a la diffusivite  ',/,
     &'@    du scalaire ',I10   ,' soit ',E14.5                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@  Le modele de turbulence k-omega a ete choisi. Pour gerer  ',/,
     &'@    correctement la suite de calcul, l''indicateur ICDPAR a ',/,
     &'@    ete mis a 1 (relecture de la distance a la paroi dans le',/,
     &'@    fichier suite).                                         ',/,
     &'@  Si cette initialisation pose probleme (modification du    ',/,
     &'@    nombre et de la position des faces de paroi depuis le   ',/,
     &'@    calcul precedent), forcer ICDPAR=-1 dans usini1 (il peut',/,
     &'@    alors y avoir un leger decalage dans la viscosite       ',/,
     &'@    turbulente au premier pas de temps).                    ',/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute.                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@  Le modele de turbulence k-omega a ete choisi, avec        ',/,
     &'@    l''option de recalcul de la distance a la paroi         ',/,
     &'@    (ICDPAR=-1). Il se peut qu''il y ait un leger decalage  ',/,
     &'@    dans la viscosite turbulente au premier pas de temps.   ',/,
     &'@                                                            ',/,
     &'@  Le calcul sera execute.                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      RETURN
      END
c@z
