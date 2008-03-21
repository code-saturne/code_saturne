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
                        SUBROUTINE RIJECH
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPHAS  , IVAR   , ISOU   , IPP    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , PRODUC , SMBR   ,
     &   COEFAX , COEFBX ,
     &   PRODUK , W2     , W3     , W4     , EPSK   , W6     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC TERMES D'ECHO DE PAROI
CFONC    POUR Rij
CFONC VAR  = R11 R22 R33 R12 R13 R23
CFONC ISOU =  1   2   3   4   5   6
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
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE PHASE                      !
CARGU ! IVAR         ! E  !  -> ! NUMERO DE VARIABLE                   !
CARGU ! ISOU         ! E  !  -> ! NUMERO DE PASSAGE                    !
CARGU ! IPP          ! E  !  -> ! NUMERO DE VARIABLE POUR SORTIES POST !
CARGU ! IFACEL       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMFBR       ! TE !  -> ! NUMERO DE FAMILLE D'UNE FACE DE BORD !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMCEL       ! TE !  -> ! NUMERO DE FAMILLE D'UNE CELLULE      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! IPRFML       ! TE !  -> ! PROPRIETES D'UNE FAMILLE             !
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
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
CARGU ! PRODUC       ! TR !  -> ! PRODUCTION                           !
CARGU !  (6,NCELET)  !    !     !                                      !
CARGU ! SMBR(NCELET  ! TR ! <-> ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! COEFAX,COEFBX! TR !  -  ! TABLEAU DE TRAVAIL POUR LES COND.    !
CARGU !  (NFABOR)    !    !     !    AUX LIMITES DE LA DIST. PAROI     !
CARGU ! PRODUK(NCELET! TR !  -  ! TABLEAU DE TRAVAIL PRODUCTION        !
CARGU ! EPSK  (NCELET! TR !  -  ! TABLEAU DE TRAVAIL EPSILON/K         !
CARGU ! W2...6(NCELET! TR !  -  ! TABLEAU DE TRAVAIL PRODUCTION        !
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
      INCLUDE "pointe.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
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
      INTEGER          IPHAS  , IVAR   , ISOU   , IPP
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION PRODUC(6,NCELET)
      DOUBLE PRECISION SMBR(NCELET)
      DOUBLE PRECISION COEFAX(NFABOR), COEFBX(NFABOR)
      DOUBLE PRECISION PRODUK(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION W4(NCELET),EPSK(NCELET),W6(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFACPT, IEL   , II    , JJ    , KK    , MM
      INTEGER          IRKM  , IRKI  , IRKJ  , ISKM  , ISKI  , ISKJ
      INTEGER          IR11IP, IR22IP, IR33IP, IR12IP, IR13IP, IR23IP
      INTEGER          IEIPH , IPCROM, IPCROO
      INTEGER          IFAC  , IDIMTE, ITENSO
      INTEGER          INC   , ICCOCG, IPHYDP, IVAR0 , IITYPH
C
      DOUBLE PRECISION CMU075, DISTXN, D2S3  , TRRIJ , XK
      DOUBLE PRECISION UNSSUR, VNK   , VNM   , VNI   , VNJ
      DOUBLE PRECISION DELTKI, DELTKJ, DELTKM, DELTIJ
      DOUBLE PRECISION AA    , BB    , SURFBN, XNORME
C
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IR11IP = IR11(IPHAS)
      IR22IP = IR22(IPHAS)
      IR33IP = IR33(IPHAS)
      IR12IP = IR12(IPHAS)
      IR13IP = IR13(IPHAS)
      IR23IP = IR23(IPHAS)
      IEIPH  = IEP (IPHAS)
C
      IPCROM = IPPROC(IROM  (IPHAS))
      IPCROO = IPCROM
      IF(ISTO2T(IPHAS).GT.0.AND.IROEXT(IPHAS).GT.0) THEN
        IPCROO = IPPROC(IROMA(IPHAS))
      ENDIF
C
      DELTIJ = 1.0D0
      IF(ISOU.GT.3) THEN
        DELTIJ = 0.0D0
      ENDIF
C
      CMU075 = CMU**0.75D0
      D2S3   = 2.D0/3.D0
C
C=======================================================================
C 2. CALCUL AUX CELLULES DES NORMALES EN PAROI CORRESPONDANTES
C=======================================================================
C
C     On stocke les composantes dans les tableaux de travail W2, W3 et W4
C
      IF(ABS(ICDPAR).EQ.2) THEN
C
C     On connait la face de paroi correspondante
C
          DO IEL = 1, NCEL
            IFACPT = IA(IIFAPA(IPHAS)-1+IEL)
            SURFBN = RA(ISRFBN-1+IFACPT)
            UNSSUR = 1.D0/SURFBN
            W2(IEL)= SURFBO(1,IFACPT)*UNSSUR
            W3(IEL)= SURFBO(2,IFACPT)*UNSSUR
            W4(IEL)= SURFBO(3,IFACPT)*UNSSUR
          ENDDO
C
      ELSEIF(ABS(ICDPAR).EQ.1) THEN
C
C     La normale est definie comme - gradient de la distance
C       a la paroi
C
C       La distance a la paroi vaut 0 en paroi
C         par definition et obeit a un flux nul ailleurs
C
        IITYPH = IITYPF+(IPHAS-1)*NFABOR
        DO IFAC = 1, NFABOR
          IF(IA(IITYPH+IFAC-1).EQ.IPAROI .OR.
     &       IA(IITYPH+IFAC-1).EQ.IPARUG) THEN
            COEFAX(IFAC) = 0.0D0
            COEFBX(IFAC) = 0.0D0
          ELSE
            COEFAX(IFAC) = 0.0D0
            COEFBX(IFAC) = 1.0D0
          ENDIF
        ENDDO
C
C       Calcul du gradient
C
        IF(IRANGP.GE.0) CALL PARCOM (RA(IDIPAR))
C                       ===========
        IF(IPERIO.EQ.1) THEN
          IDIMTE = 0
          ITENSO = 0
          CALL PERCOM
C         ===========
     &   ( IDIMTE , ITENSO ,
     &     RA(IDIPAR), RA(IDIPAR), RA(IDIPAR),
     &     RA(IDIPAR), RA(IDIPAR), RA(IDIPAR),
     &     RA(IDIPAR), RA(IDIPAR), RA(IDIPAR))
        ENDIF
C
        INC    = 1
        ICCOCG = 1
        IPHYDP = 0
        IVAR0  = 0
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGY , IMLIGY , IPHYDP ,
     &   IWARNY , NFECRA ,
     &   EPSRGY , CLIMGY , EXTRAY ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA     , RA     , RA     ,
     &   RA(IDIPAR) , COEFAX , COEFBX ,
     &   W2     , W3     , W4     ,
C        ------   ------   ------
     &   PRODUK , EPSK   , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C     Normalisation (attention, le gradient peut etre nul, parfois)
C
        DO IEL = 1 ,NCEL
          XNORME = MAX(SQRT(W2(IEL)**2+W3(IEL)**2+W4(IEL)**2),EPZERO)
          W2(IEL) = -W2(IEL)/XNORME
          W3(IEL) = -W3(IEL)/XNORME
          W4(IEL) = -W4(IEL)/XNORME
        ENDDO

      ENDIF
C
C=======================================================================
C 2. CALCUL DE VARIABLES DE TRAVAIL
C=======================================================================
C
C
C
C ---> Production et k
C
      DO IEL = 1 , NCEL
        PRODUK(IEL) = 0.5D0 * (PRODUC(1,IEL)  + PRODUC(2,IEL)  +
     &                         PRODUC(3,IEL))
        XK          = 0.5D0 * (RTPA(IEL,IR11IP) + RTPA(IEL,IR22IP) +
     &                         RTPA(IEL,IR33IP))
        EPSK(IEL)   = RTPA(IEL,IEIPH)/XK
      ENDDO
C
C
C
C ---> Indices des tensions
C
      IF     ((ISOU.EQ.1).OR.(ISOU.EQ.4).OR.(ISOU.EQ.5)) THEN
        II = 1
      ELSEIF ((ISOU.EQ.2).OR.(ISOU.EQ.6)) THEN
        II = 2
      ELSEIF  (ISOU.EQ.3) THEN
        II = 3
      ENDIF
C
      IF     ((ISOU.EQ.3).OR.(ISOU.EQ.5).OR.(ISOU.EQ.6)) THEN
        JJ = 3
      ELSEIF ((ISOU.EQ.2).OR.(ISOU.EQ.4)) THEN
        JJ = 2
      ELSEIF  (ISOU.EQ.1) THEN
        JJ = 1
      ENDIF
C
C ---> Boucle pour construction des termes sources
C
      DO IEL = 1, NCEL
        W6(IEL) = 0.D0
      ENDDO
C
      DO KK = 1, 3
C
C ---> Sommes sur m
C
        DO MM = 1, 3
C
C   --> Delta km
C
          IF(KK.EQ.MM) THEN
            DELTKM = 1.0D0
          ELSE
            DELTKM = 0.0D0
          ENDIF
C
C  --> R km
C
          IF     ((KK*MM).EQ.1) THEN
            IRKM = IR11IP
            ISKM = 1
          ELSEIF ((KK*MM).EQ.4) THEN
            IRKM = IR22IP
            ISKM = 2
          ELSEIF ((KK*MM).EQ.9) THEN
            IRKM = IR33IP
            ISKM = 3
          ELSEIF ((KK*MM).EQ.2) THEN
            IRKM = IR12IP
            ISKM = 4
          ELSEIF ((KK*MM).EQ.3) THEN
            IRKM = IR13IP
            ISKM = 5
          ELSEIF ((KK*MM).EQ.6) THEN
            IRKM = IR23IP
            ISKM = 6
          ENDIF
C
C  --> Termes en R km et Phi km
C
          DO IEL = 1, NCEL
C
            IF    (KK.EQ.1) THEN
              VNK    = W2(IEL)
            ELSEIF(KK.EQ.2) THEN
              VNK    = W3(IEL)
            ELSEIF(KK.EQ.3) THEN
              VNK    = W4(IEL)
            ENDIF
C
            IF    (MM.EQ.1) THEN
              VNM    = W2(IEL)
            ELSEIF(MM.EQ.2) THEN
              VNM    = W3(IEL)
            ELSEIF(MM.EQ.3) THEN
              VNM    = W4(IEL)
            ENDIF
C
            W6(IEL) = W6(IEL) + VNK*VNM*DELTIJ*(
     &             CRIJP1*RTPA(IEL,IRKM)*EPSK(IEL)
     &            -CRIJP2
     &             *CRIJ2*(PRODUC(ISKM,IEL)-D2S3*PRODUK(IEL)*DELTKM) )
          ENDDO
C
        ENDDO
C
C ---> Fin des sommes sur m
C
C
C  --> R ki
C
        IF     ((KK*II).EQ.1) THEN
          IRKI = IR11IP
          ISKI = 1
        ELSEIF ((KK*II).EQ.4) THEN
          IRKI = IR22IP
          ISKI = 2
        ELSEIF ((KK*II).EQ.9) THEN
          IRKI = IR33IP
          ISKI = 3
        ELSEIF ((KK*II).EQ.2) THEN
          IRKI = IR12IP
          ISKI = 4
        ELSEIF ((KK*II).EQ.3) THEN
          IRKI = IR13IP
          ISKI = 5
        ELSEIF ((KK*II).EQ.6) THEN
          IRKI = IR23IP
          ISKI = 6
        ENDIF
C
C  --> R kj
C
        IF     ((KK*JJ).EQ.1) THEN
          IRKJ = IR11IP
          ISKJ = 1
        ELSEIF ((KK*JJ).EQ.4) THEN
          IRKJ = IR22IP
          ISKJ = 2
        ELSEIF ((KK*JJ).EQ.9) THEN
          IRKJ = IR33IP
          ISKJ = 3
        ELSEIF ((KK*JJ).EQ.2) THEN
          IRKJ = IR12IP
          ISKJ = 4
        ELSEIF ((KK*JJ).EQ.3) THEN
          IRKJ = IR13IP
          ISKJ = 5
        ELSEIF ((KK*JJ).EQ.6) THEN
          IRKJ = IR23IP
          ISKJ = 6
        ENDIF
C
C   --> Delta ki
C
        IF (KK.EQ.II) THEN
          DELTKI = 1.D0
        ELSE
          DELTKI = 0.D0
        ENDIF
C
C   --> Delta kj
C
        IF (KK.EQ.JJ) THEN
          DELTKJ = 1.D0
        ELSE
          DELTKJ = 0.D0
        ENDIF
C
        DO IEL = 1, NCEL
C
            IF    (KK.EQ.1) THEN
              VNK    = W2(IEL)
            ELSEIF(KK.EQ.2) THEN
              VNK    = W3(IEL)
            ELSEIF(KK.EQ.3) THEN
              VNK    = W4(IEL)
            ENDIF
C
            IF    (II.EQ.1) THEN
              VNI    = W2(IEL)
            ELSEIF(II.EQ.2) THEN
              VNI    = W3(IEL)
            ELSEIF(II.EQ.3) THEN
              VNI    = W4(IEL)
            ENDIF
C
            IF    (JJ.EQ.1) THEN
              VNJ    = W2(IEL)
            ELSEIF(JJ.EQ.2) THEN
              VNJ    = W3(IEL)
            ELSEIF(JJ.EQ.3) THEN
              VNJ    = W4(IEL)
            ENDIF
C
          W6(IEL) = W6(IEL) + 1.5D0*VNK*(
     &    -CRIJP1*(RTPA(IEL,IRKI)*VNJ+RTPA(IEL,IRKJ)*VNI)*EPSK(IEL)
     &    +CRIJP2
     &     *CRIJ2*((PRODUC(ISKI,IEL)-D2S3*PRODUK(IEL)*DELTKI)*VNJ
     &            +(PRODUC(ISKJ,IEL)-D2S3*PRODUK(IEL)*DELTKJ)*VNI) )
C
        ENDDO
C
      ENDDO
C
C
C ---> Distance a la paroi et fonction d'amortissement : W3
C     Pour chaque mode de calcul : meme code, test
C       en dehors de la boucle
C
      IF(ABS(ICDPAR).EQ.2) THEN
        DO IEL = 1 , NCEL
          IFACPT = IA(IIFAPA(IPHAS)-1+IEL)
          DISTXN =
     &          (CDGFBO(1,IFACPT)-XYZCEN(1,IEL))**2
     &         +(CDGFBO(2,IFACPT)-XYZCEN(2,IEL))**2
     &         +(CDGFBO(3,IFACPT)-XYZCEN(3,IEL))**2
          DISTXN = SQRT(DISTXN)
          TRRIJ  = 0.5D0 * (RTPA(IEL,IR11IP) + RTPA(IEL,IR22IP) +
     &                         RTPA(IEL,IR33IP))
          AA = 1.D0
          BB = CMU075*TRRIJ**1.5D0/(XKAPPA*RTPA(IEL,IEIPH)*DISTXN)
          W3(IEL) = MIN(AA, BB)
        ENDDO
      ELSE
        DO IEL = 1 , NCEL
          DISTXN =  MAX(RA(IDIPAR+IEL-1),EPZERO)
          TRRIJ  = 0.5D0 * (RTPA(IEL,IR11IP) + RTPA(IEL,IR22IP) +
     &                         RTPA(IEL,IR33IP))
          AA = 1.D0
          BB = CMU075*TRRIJ**1.5D0/(XKAPPA*RTPA(IEL,IEIPH)*DISTXN)
          W3(IEL) = MIN(AA, BB)
        ENDDO
      ENDIF
C
C
C ---> Increment du terme source
C
      DO IEL = 1, NCEL
        SMBR(IEL) = SMBR(IEL)
     &            + PROPCE(IEL,IPCROO)*VOLUME(IEL)*W6(IEL)*W3(IEL)
      ENDDO
C
C
      RETURN
C
      END
c@z
