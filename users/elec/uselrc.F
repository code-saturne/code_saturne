c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS
C
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
                        SUBROUTINE USELRC
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , VISCF  , VISCB  ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  ROUTINE PHYSIQUE PARTICULIERE POUR LE MODULE ELECTRIQUE
CFONC
CFONC              CALCULS DU COEFFICIENT DE RECALAGE
CFONC                POUR LES VARIABLES ELECTIQUES
CFONC              RECALAGE DES VARIABLES ELECTRIQUES
CFONC                EN FONCTION DE CE COEFFICIENT
CFONC
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE                !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NFML         ! E  !  -> ! NOMBRE DE FAMILLES D ENTITES         !
CARGU ! NPRFML       ! E  !  -> ! NOMBRE DE PROPRIETESE DES FAMILLES   !
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE SOMMETS                    !
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IFACEL       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMFBR       ! TE !  -> ! NUMERO DE FAMILLE D'UNE FACE DE BORD !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMCEL       ! TE !  -> ! NUMERO DE FAMILLE D'UNE CELLULE      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! IPRFML       ! TE !  -> ! PROPRIETES D'UNE FAMILLE             !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! ITYPSM       ! TE !  -> ! TYPE DE SOURCE DE MASSE POUR LES     !
CARGU ! (NCESMP,NVAR)!    !     !  VARIABLES (cf. USTSMA)              !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! XYZCEN       ! TR !  -> ! POINT ASSOCIES AUX VOLUMES DE CONTROL!
CARGU ! (NDIM,NCELET !    !     !                                      !
CARGU ! SURFAC       ! TR !  -> ! VECTEUR SURFACE DES FACES INTERNES   !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! SURFBO       ! TR !  -> ! VECTEUR SURFACE DES FACES DE BORD    !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! CDGFAC       ! TR !  -> ! CENTRE DE GRAVITE DES FACES INTERNES !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! CDGFBO       ! TR !  -> ! CENTRE DE GRAVITE DES FACES DE BORD  !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS (OPTIONNEL)    !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! SMACEL       ! TR !  -> ! VALEUR DES VARIABLES ASSOCIEE A LA   !
CARGU ! (NCESMP,*   )!    !     !  SOURCE DE MASSE                     !
CARGU !              !    !     !  POUR IVAR=IPR, SMACEL=FLUX DE MASSE !
CARGU ! VISCF(NFAC)  ! TR !  -  ! TABLEAU DE TRAVAIL    FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! TABLEAU DE TRAVAIL    FACES DE BORD  !
CARGU ! W1..9(NCELET ! TR !  -  ! TABLEAU DE TRAVAIL    CELLULES       !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
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
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "elincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE)
      INTEGER          ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IEL    , IFAC   , IUTILE
      INTEGER          IPCEFJ , IPCDC1 , IPCDC2 , IPCDC3 , IPCSIG
      INTEGER          IPDCRP , IDIMVE
C
      DOUBLE PRECISION SOMJE , COEPOA , COEFAV , COEPOT
      DOUBLE PRECISION EMAX  , AIEX   , AMEX
      DOUBLE PRECISION RAYO  , ECONS  , Z1     , Z2   , POSI
      DOUBLE PRECISION ELCOU , DTJ    , DTJM   , DELHSH , CDTJ , CPMX
      DOUBLE PRECISION XELEC , YELEC  , ZELEC
C
C***********************************************************************
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
C=======================================================================
C 2.  ARC ELECTRIQUE
C=======================================================================
C
C
C
      IF ( IPPMOD(IELARC).GE.1 ) THEN
C
C       3 Exemples : pour activer l'un des 3 mettre IUTILE= 1, 2 ou 3
C       Par defaut IUTILE = 1
C
        IUTILE = 1
C
C 2.1 : 1er exemple : cas general
C ===============================
C
        IF ( IUTILE .EQ. 1) THEN
C
C       CALCUL DU COEFFICIENT DE RECALAGE
C       -------------------------------
C
C  Calcul de l'integrale sur le Volume de J.E
C     (c'est forcement positif ou nul)
C
          IPCEFJ = IPPROC(IEFJOU)
          SOMJE = 0.D0
          DO IEL = 1, NCEL
            SOMJE = SOMJE+PROPCE(IEL,IPCEFJ)*VOLUME(IEL)
          ENDDO
C
          IF(IRANGP.GE.0) THEN
            CALL PARSOM (SOMJE)
          ENDIF
C
          COEPOT = COUIMP*DPOT/MAX(SOMJE,EPZERO)
C
          COEPOA = COEPOT
C
C  On impose COEPOT >= 0.75 et COEPOT <= 1.5
C
          IF ( COEPOT .GT. 1.50D0 ) COEPOT = 1.50D0
          IF ( COEPOT .LT. 0.75D0 ) COEPOT = 0.75D0
C
          WRITE(NFECRA,1000)COEPOA,COEPOT
 1000     FORMAT(/,
     & ' Courant impose/Courant= ',E14.5,', Coeff. recalage= ',E14.5)
C
C       RECALAGE DES VARIABLES ELECTRIQUES
C       ---------------------------------------
C
C        Valeur de DPOT
C        --------------
C
          DPOT = DPOT*COEPOT
C
C        Potentiel Electrique (on pourrait eviter ; c'est pour le post)
C        --------------------
C
          DO IEL = 1, NCEL
            RTP(IEL,ISCA(IPOTR)) = RTP(IEL,ISCA(IPOTR))*COEPOT
          ENDDO
C
C
C        Densite de courant (sert pour A et pour jXB)
C        ------------------
C
          IF(IPPMOD(IELARC).GE.1 ) THEN
            DO IDIMVE = 1, NDIMVE
              IPDCRP = IPPROC(IDJR(IDIMVE))
              DO IEL = 1, NCEL
                PROPCE(IEL,IPDCRP) = PROPCE(IEL,IPDCRP) * COEPOT
              ENDDO
            ENDDO
          ENDIF
C
C        Effet Joule (sert pour H au pas de temps suivant)
C        -----------
C
          IPCEFJ = IPPROC(IEFJOU)
          DO IEL = 1, NCEL
            PROPCE(IEL,IPCEFJ) = PROPCE(IEL,IPCEFJ)*COEPOT**2
          ENDDO
C
C Fin 1er exemple
C
        ELSE IF ( IUTILE .EQ. 2) THEN
C
C 2.2 : 2eme exemple : Autre methode de recalage
C ==============================================
C    Ceci est un cas particulier et doit etre adapte en fonction
C    du cas et du maillage (intervenir aussi dans uselcl)
C
C        Calcul de l'integrale sur le Volume de J.E
C        -----------------------------------
C        (c'est forcement positif ou nul)
C
          IPCEFJ = IPPROC(IEFJOU)
          SOMJE = 0.D0
          DO IEL = 1, NCEL
            SOMJE = SOMJE+PROPCE(IEL,IPCEFJ)*VOLUME(IEL)
          ENDDO
C
          IF(IRANGP.GE.0) THEN
            CALL PARSOM (SOMJE)
          ENDIF
C
          IF (SOMJE .NE. 0) THEN
            COEPOT = COUIMP*DPOT/MAX(SOMJE,EPZERO)
          ENDIF
          WRITE(NFECRA,1001) COUIMP,DPOT,SOMJE
C
C        Calcul de l'intensite du courant d'arc
C        --------------------------------------
C          Calcul de l'integrale de J sur une surface plane
C          perpendiculaire a l'axe de l'arc
C
C       ATTENTION : changer la valeur des tests sur CDGFAC(3,IFAC)
C                   en fonction du maillage
C
          IPCDC3 = IPPROC(IDJR(3))
          ELCOU  = 0.D0
          DO IFAC = 1, NFAC
            IF( SURFAC(1,IFAC).EQ.0.D0 .AND. SURFAC(2,IFAC).EQ.0.D0
     &           .AND. CDGFAC(3,IFAC) .LT. 0.7D-2
     &           .AND. CDGFAC(3,IFAC) .GT. 0.65D-2 ) THEN
              IEL = IFACEL(1,IFAC)
              ELCOU = ELCOU + PROPCE(IEL,IPCDC3) * SURFAC(3,IFAC)
            ENDIF
          ENDDO
C
          IF(IRANGP.GE.0) THEN
            CALL PARSOM (ELCOU)
          ENDIF
C
          IF ( ABS(ELCOU).GE.1.D-06 ) THEN
            ELCOU=ABS(ELCOU)
          ELSE
            ELCOU=0.D0
          ENDIF
          IF(ELCOU.NE.0.D0) COEPOA = COUIMP/ELCOU
          COEPOT = COEPOA
C
          WRITE(NFECRA,*) ' ELCOU = ',ELCOU
C
          DTJ = 1.D15
          DTJM =DTJ
          DELHSH = 0.D0
          CDTJ= 2.0D2
C
          DO IEL = 1, NCEL
            IF(PROPCE(IEL,IPPROC(IROM(1))).NE.0.D0)
     &           DELHSH =  PROPCE(IEL,IPCEFJ) * DT(IEL)
     &           /PROPCE(IEL,IPPROC(IROM(1)))
C
            IF(DELHSH.NE.0.D0) THEN
              DTJM= RTP(IEL,ISCA(ISCALT(1)))/DELHSH
            ELSE
              DTJM= DTJ
            ENDIF
            DTJM=ABS(DTJM)
            DTJ =MIN(DTJ,DTJM)
          ENDDO
          IF(IRANGP.GE.0) THEN
            CALL PARMIN (DTJ)
          ENDIF
          WRITE(NFECRA,*) ' DTJ = ',DTJ
C
          CPMX= SQRT(CDTJ*DTJ)
          COEPOT=CPMX
          IF(NTCABS.GT.5) THEN
            IF(COEPOA.GE.1.05D0 .AND. COEPOT.LE.CPMX) THEN
              COEPOT=CPMX
            ELSE
              COEPOT=COEPOA
            ENDIF
          ENDIF
C
          WRITE(NFECRA,1008)CPMX,COEPOA,COEPOT
          WRITE(NFECRA,1009)ELCOU,DPOT*COEPOT
C
C        RECALAGE DES VARIABLES ELECTRIQUES
C        ----------------------------------
C
C         Valeur de DPOT
C         --------------
C
          DPOT = DPOT*COEPOT
C
C         Potentiel Electrique (on pourrait eviter ; c'est pour le post)
C         --------------------
C
          DO IEL = 1, NCEL
            RTP(IEL,ISCA(IPOTR)) = RTP(IEL,ISCA(IPOTR))*COEPOT
          ENDDO
C
C
C      Densite de courant (sert pour A et pour jXB)
C      ------------------
C
          IF(IPPMOD(IELARC).GE.1 ) THEN
            DO IDIMVE = 1, NDIMVE
              DO IEL = 1, NCEL
                IPDCRP = IPPROC(IDJR(IDIMVE))
                PROPCE(IEL,IPDCRP) = PROPCE(IEL,IPDCRP) * COEPOT
              ENDDO
            ENDDO
          ENDIF
C
C      Effet Joule (sert pour H au pas de temps suivant)
C      -----------
C
          IPCEFJ = IPPROC(IEFJOU)
          DO IEL = 1, NCEL
            PROPCE(IEL,IPCEFJ) = PROPCE(IEL,IPCEFJ)*COEPOT**2
          ENDDO
C
        ELSE IF ( IUTILE .EQ. 3) THEN
C
C 2.3 :  3eme exemple : cas avec claquage
C =======================================
C    Ceci est un cas particulier et doit etre adapte en fonction
C    du cas et du maillage (intervenir aussi dans uselcl)
C
C
C        Utilisation d'une rampe d'intensite
C        -----------------------------------
C
          IF ( NTCABS.LE.200 ) THEN
            COUIMP = 200.D0
          ENDIF
C
          IF ( NTCABS.GT.200.AND.NTCABS.LE.400 ) THEN
            COUIMP = 200.D0 + 2 * (NTCABS-200)
          ENDIF
C
          IF ( NTCABS.GT.400 ) THEN
            COUIMP = 600.D0
          ENDIF
C
C        UTILISANT D'UN CLAQUAGE AUTO
C        ----------------------------
C
          IF(NTCABS.LE.400.or.NTCABS.EQ.NTPABS+1) ICLAQ = 0
C
          ECONS = 1.5D5
C
C        ON REPERE SI IL Y A CLAQUAGE ET SI OUI OU
C        -----------------------------------------
C
          IF(NTCABS.GE.400 .AND. ICLAQ .EQ. 0 ) THEN
            AMEX = 1.D30
            AIEX = -1.D30
            EMAX = 0.D0
C
C     les composantes du champ electrique : J/SIGMA
C
            IPCDC1 = IPPROC(IDJR(1))
            IPCDC2 = IPPROC(IDJR(2))
            IPCDC3 = IPPROC(IDJR(3))
            IPCSIG = IPPROC(IVISLS(IPOTR))
C
            DO IEL = 1, NCEL
C
              XELEC = PROPCE(IEL,IPCDC1)/PROPCE(IEL,IPCSIG)
              YELEC = PROPCE(IEL,IPCDC2)/PROPCE(IEL,IPCSIG)
              ZELEC = PROPCE(IEL,IPCDC3)/PROPCE(IEL,IPCSIG)
C
              W1(IEL) = SQRT ( XELEC**2 + YELEC**2 + ZELEC**2 )
C
              AMEX =  MIN(AMEX,W1(IEL))
              AIEX =  MAX(AIEX,W1(IEL))
              IF( W1(IEL) .GE. ECONS) THEN
                WRITE(NFECRA,*) 'claquage ', NTCABS, W1(IEL)
                ICLAQ = 1
                NTDCLA = NTCABS
                IF(W1(IEL).GT.EMAX) THEN
                  XCLAQ =  XYZCEN(1,IEL)
                  YCLAQ =  XYZCEN(2,IEL)
                  ZCLAQ =  XYZCEN(3,IEL)
                  EMAX = W1(IEL)
                ENDIF
              ENDIF
            ENDDO
C
            WRITE(NFECRA,*)
            WRITE(NFECRA,*) ' NT min et max de E = ',NTCABS,AMEX,AiEX
            WRITE(NFECRA,*)
            WRITE(NFECRA,*) XCLAQ,YCLAQ,ZCLAQ,NTDCLA
C
          ENDIF
C
C        SI IL Y A CLAQUAGE : ON IMPOSE COLONNE CHAUDE DU CENTRE VERS
C        LE POINT DE CLAQUAGE
C        =============================================================
C
          IF(ICLAQ .EQ. 1) THEN
            IF(NTCABS.LE.NTDCLA+30) THEN
              Z1 = ZCLAQ - 3.D-4
              IF(Z1.LE.0.D0) Z1 = 0.D0
              Z2 = ZCLAQ + 3.D-4
              IF(Z2.GE.2.D-2) Z2 = 2.D-2
C
              DO IEL = 1, NCEL
C
                IF( XYZCEN(3,IEL).GE.Z1 .AND. XYZCEN(3,IEL).LE.Z2) THEN
                  RAYO = SQRT((XCLAQ*XYZCEN(1,IEL)-YCLAQ*XYZCEN(2,IEL)
     &                 /SQRT(XCLAQ**2+YCLAQ**2))**2+(XYZCEN(3,IEL)
     &                 -ZCLAQ)**2)
                  POSI=XCLAQ*XYZCEN(1,IEL)
                  IF( RAYO.LE.5D-4 .AND. POSI.GE.0D0 ) THEN
C                    RTP(IEL,ISCA(IHM)) = 16.D6
                    RTP(IEL,ISCA(IHM)) = 8.D7
                  ENDIF
                ENDIF
              ENDDO
            ELSE
              ICLAQ = 0
            ENDIF
          ENDIF
C
C        Calcul de l'integrale sur le Volume de J.E
C        -----------------------------------
C        (c'est forcement positif ou nul)
C
          IPCEFJ = IPPROC(IEFJOU)
          SOMJE = 0.D0
          DO IEL = 1, NCEL
            SOMJE = SOMJE+PROPCE(IEL,IPCEFJ)*VOLUME(IEL)
          ENDDO
C
          IF(IRANGP.GE.0) THEN
            CALL PARSOM (SOMJE)
          ENDIF
C
          IF (SOMJE .NE. 0) THEN
            COEPOT = COUIMP*DPOT/MAX(SOMJE,EPZERO)
          ENDIF
          WRITE(NFECRA,1001) COUIMP,DPOT,SOMJE
C
C        Calcul de l'intensite du courant d'arc
C        --------------------------------------
C          Calcul de l'integrale de J sur une surface plane
C          perpendiculaire a l'axe de l'arc
C
C       ATTENTION : changer la valeur des tests sur CDGFAC(3,IFAC)
C                   en fonction du maillage
C
          IPCDC3 = IPPROC(IDJR(3))
          ELCOU = 0.D0
          DO IFAC = 1, NFAC
            IF( SURFAC(1,IFAC).EQ.0.D0 .AND. SURFAC(2,IFAC).EQ.0.D0
     &           .AND. CDGFAC(3,IFAC) .GT. 0.05D-2
     &           .AND. CDGFAC(3,IFAC) .LT. 0.08D-2 ) THEN
              IEL = IFACEL(1,IFAC)
              ELCOU = ELCOU + PROPCE(IEL,IPCDC3) * SURFAC(3,IFAC)
            ENDIF
          ENDDO
C
          IF(IRANGP.GE.0) THEN
            CALL PARSOM (ELCOU)
          ENDIF
C
          IF ( ABS(ELCOU).GE.1.D-06 ) THEN
            ELCOU=ABS(ELCOU)
          ELSE
            ELCOU=0.D0
          ENDIF
          IF(ELCOU.NE.0.D0) COEPOA = COUIMP/ELCOU
          COEPOT = COEPOA
C
          WRITE(NFECRA,*) ' ELCOU = ',ELCOU
C
          DTJ = 1.D15
          DTJM =DTJ
          DELHSH = 0.D0
          CDTJ= 2.0D2
C
          DO IEL = 1, NCEL
            IF(PROPCE(IEL,IPPROC(IROM(1))).NE.0.D0)
     &           DELHSH =  PROPCE(IEL,IPCEFJ) * DT(IEL)
     &           /PROPCE(IEL,IPPROC(IROM(1)))
C
            IF(DELHSH.NE.0.D0) THEN
              DTJM= RTP(IEL,ISCA(ISCALT(1)))/DELHSH
            ELSE
              DTJM= DTJ
            ENDIF
            DTJM=ABS(DTJM)
            DTJ =MIN(DTJ,DTJM)
          ENDDO
C
          IF(IRANGP.GE.0) THEN
            CALL PARMIN (DTJ)
          ENDIF
C
          CPMX= SQRT(CDTJ*DTJ)
          COEPOT=CPMX
          IF(NTCABS.GT.5) THEN
            IF(COEPOA.GE.1.05D0 .AND. COEPOT.LE.CPMX) THEN
              COEPOT=CPMX
            ELSE
              COEPOT=COEPOA
            ENDIF
          ENDIF
C
          WRITE(NFECRA,1008)CPMX,COEPOA,COEPOT
          WRITE(NFECRA,1009)ELCOU,DPOT*COEPOT
C
C        RECALAGE DES VARIABLES ELECTRIQUES
C        ----------------------------------
C
C         Valeur de DPOT
C         --------------
C
          DPOT = DPOT*COEPOT
C
C         Potentiel Electrique (on pourrait eviter ; c'est pour le post)
C         --------------------
C
          DO IEL = 1, NCEL
            RTP(IEL,ISCA(IPOTR)) = RTP(IEL,ISCA(IPOTR))*COEPOT
          ENDDO
C
C
C      Densite de courant (sert pour A et pour jXB)
C      ------------------
C
          IF(IPPMOD(IELARC).GE.1 ) THEN
            DO IDIMVE = 1, NDIMVE
              DO IEL = 1, NCEL
                IPDCRP = IPPROC(IDJR(IDIMVE))
                PROPCE(IEL,IPDCRP) = PROPCE(IEL,IPDCRP) * COEPOT
              ENDDO
            ENDDO
          ENDIF
C
C      Effet Joule (sert pour H au pas de temps suivant)
C      -----------
C
          IPCEFJ = IPPROC(IEFJOU)
          DO IEL = 1, NCEL
            PROPCE(IEL,IPCEFJ) = PROPCE(IEL,IPCEFJ)*COEPOT**2
          ENDDO
C
        ELSE
          WRITE(NFECRA,5000) IUTILE
          CALL CSEXIT(1)
        ENDIF
      ENDIF
C
C=======================================================================
C 3.  EFFET JOULE
C=======================================================================
C
      IF ( IPPMOD(IELJOU).GE.1 ) THEN
C
C 3.1  CALCUL DU COEFFICIENT DE RECALAGE
C --------------------------------------
C
C  Calcul de l'integrale sur le Volume de J.E
C     (c'est forcement positif ou nul)
C
        IPCEFJ = IPPROC(IEFJOU)
        SOMJE = 0.D0
        DO IEL = 1, NCEL
          SOMJE = SOMJE+PROPCE(IEL,IPCEFJ)*VOLUME(IEL)
        ENDDO
C
        IF(IRANGP.GE.0) THEN
          CALL PARSOM (SOMJE)
        ENDIF
C
        COEPOT = SQRT(PUISIM/MAX(SOMJE,EPZERO))
C
        COEFAV = COEPOT
C
C  On impose COEF >= 0.75 et COEF <= 1.5
C
        IF ( COEPOT .GT. 1.50D0 ) COEPOT = 1.50D0
        IF ( COEPOT .LT. 0.75D0 ) COEPOT = 0.75D0
C
        WRITE(NFECRA,2000)COEFAV,COEJOU
 2000   FORMAT(/,
     & ' Puissance impose/Somme jE= ',E14.5,', Coeff. recalage= ',E14.5)
C
C
C 3.2  RECALAGE DES VARIABLES JOULE
C ---------------------------------
C
C       Valeur de DPOT (au cas ou utile)
C       --------------
C
        DPOT = DPOT*COEPOT
C
C       Coefficient correcteur COEJOU cumule
C       ------------------------------------
C
        COEJOU = COEJOU*COEPOT
C
C       Potentiel Electrique (on pourrait eviter ; c'est pour le post)
C       --------------------
C
        IF ( IPPMOD(IELJOU).NE.3 .AND. IPPMOD(IELJOU).NE.4 ) THEN
          DO IEL = 1, NCEL
            RTP(IEL,ISCA(IPOTR)) = RTP(IEL,ISCA(IPOTR))*COEPOT
          ENDDO
        ENDIF
C
C      Potentiel complexe (on pourrait eviter ; c'est pour le post)
C      -----------------
C
        IF ( IPPMOD(IELJOU).EQ.2 ) THEN
          DO IEL = 1, NCEL
            RTP(IEL,ISCA(IPOTI)) = RTP(IEL,ISCA(IPOTI))*COEPOT
          ENDDO
        ENDIF
C
C
C      Effet Joule (sert pour H au pas de temps suivant)
C      -----------
C
        IPCEFJ = IPPROC(IEFJOU)
        DO IEL = 1, NCEL
          PROPCE(IEL,IPCEFJ) = PROPCE(IEL,IPCEFJ)*COEPOT**2
        ENDDO
C
      ENDIF
C
C--------
C FORMATS
C--------
C
 1001  FORMAT(/, ' Courant impose= ',E14.5, /,
     &              ' Dpot= ',E14.5,/,
     &              ' Somje= ',E14.5)
C
 1008  FORMAT(/,' Cpmx   = ',E14.5,/,
     &          ' COEPOA = ',E14.5,/,
     &          ' COEPOT = ',E14.5)
C
 1009  FORMAT(/,' Courant calcule = ',E14.5,/,
     &          ' Dpot recale     = ',E14.5)
C
 5000  FORMAT(/,' ERREUR DANS USELRC :',/,
     &          ' VALEUR NON PERMISE DE IUTILE ',/,
     &          ' VERIFIER VOS DONNEES ')
C
C----
C FIN
C----
C
      RETURN
      END
c@z
