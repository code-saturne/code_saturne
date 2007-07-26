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
C     Copyright (C) 1998-2007 EDF S.A., France
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
                        SUBROUTINE PHYVAR
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
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     , W9     , W10    , XMIJ   ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC REMPLISSAGE DES GRANDEURS PHYSIQUES VARIABLES EN TEMPS
CFONC     ESSENTIELLEMENT LA VISCOSITE TURBULENTE.
CFONC     ON APPELLE UN SOUS PROGRAMME UTILISATEUR QUI PERMET DE
CFONC     SPECIFIER ROM, ROMB, VISCL, VISCLS ...
CFONC
c@fonce
C                             ARGUMENTS
c@argub
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
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! W1..10(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! XMIJ(NCELET,6! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "pointe.h"
      INCLUDE "albase.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "matiss.h"
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
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET),W4(NCELET)
      DOUBLE PRECISION W5(NCELET),W6(NCELET),W7(NCELET),W8(NCELET)
      DOUBLE PRECISION W9(NCELET),W10(NCELET),XMIJ(NCELET,6)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IVAR  , IEL   , IFAC  , ISCAL , IPHAS , IPH
      INTEGER          II    , IOK   , IOK1  , IOK2  , IISCT
      INTEGER          NPHMX , NN
      INTEGER          IBROM(NPHSMX) , IPCROM, IPBROM, IPCVST
      INTEGER          IKIPH , IEIPH , IR11IP, IR22IP, IR33IP
      INTEGER          IPCCP , IPCVIS, IPHIPH, IPCVMA
      INTEGER          IDIMTE, ITENSO, ICLIPC
      DOUBLE PRECISION XK, XE, XNU, XROM, VISMAX(NSCAMX), VISMIN(NSCAMX)
      DOUBLE PRECISION VARMN(4), VARMX(4), TT, TTMIN, TTKE, VISTOT
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IPASS = IPASS + 1
C
C=======================================================================
C 2.  PREPARATION DE LA PERIODICITE DE ROTATION
C       CALCUL DE DUDXYZ ET DRDXYZ (gradients sur les halos avec prise
C       en compte des periodicites pour exploitation dans pering, inimas)
C=======================================================================
C
      IF(IPEROT.GT.0) THEN
C
        DO IPHAS = 1,NPHAS
C
          IPH = IPHAS
          CALL PERINU
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IPH    ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER ,
     &   RA(IDUDXY) ,
     &   RA     )
C
        ENDDO
C
C
        DO IPHAS = 1, NPHAS
C
          IF(ITYTUR(IPHAS).EQ.3) THEN
C
            IPH = IPHAS
            CALL PERINR
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IPH    ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER ,
     &   RA(IDRDXY) ,
     &   RA     )
C
          ENDIF
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C 4.  ON REND LA MAIN A L'UTILISATEUR POUR LA PROGRAMMATION DES
C      GRANDEURS PHYSIQUES VARIABLES QUI LUI SONT PROPRES
C=======================================================================
C
      NPHMX = NPHSMX
      DO IPHAS = 1, NPHSMX
        IBROM(IPHAS) = 0
      ENDDO
C
C
      IF (NSCAPP.GT.0) THEN
        CALL PPPHYV
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IBROM  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C
      IF (IMATIS.EQ.1) THEN
C
C     Matisse
        CALL MTPHYV
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IBROM  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   ,
     &   PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C
      ELSE
C
C     Utilisateur standard
        CALL USPHYV
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IBROM  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   ,
     &   PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C  ROMB SUR LES BORDS : VALEUR PAR DEFAUT (CELLE DE LA CELLULE VOISINE)
C
      DO IPHAS = 1, NPHAS
        IF (IBROM(IPHAS).EQ.0) THEN
          IPCROM = IPPROC(IROM(IPHAS))
          IPBROM = IPPROB(IROM(IPHAS))
          DO IFAC = 1, NFABOR
            IEL = IFABOR(IFAC)
            PROPFB(IFAC,IPBROM) = PROPCE(IEL ,IPCROM)
          ENDDO
        ENDIF
      ENDDO
C
C  Au premier pas de temps du calcul
C     Si on a indique que rho (visc) etait constant
C       et qu'on l'a modifie dans usphyv, ca ne va pas
C     On se sert de irovar (ivivar) pour ecrire et lire
C       rho (visc) dans le fichier suite
C
      IF(NTCABS.EQ.NTPABS+1) THEN
C
C     Masse volumique aux cellules et aux faces de bord
        IOK1 = 0
        DO IPHAS = 1, NPHAS
          IF(IROVAR(IPHAS).EQ.0) THEN
            IPCROM = IPPROC(IROM(IPHAS))
            IPBROM = IPPROB(IROM(IPHAS))
            DO IEL = 1, NCEL
              IF( ABS(PROPCE(IEL ,IPCROM)-RO0   (IPHAS)).GT.EPZERO) THEN
                IOK1 = 1
              ENDIF
            ENDDO
            DO IFAC = 1, NFABOR
              IF( ABS(PROPFB(IFAC,IPBROM)-RO0   (IPHAS)).GT.EPZERO) THEN
                IOK1 = 1
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        IF(IOK1.NE.0) THEN
          WRITE(NFECRA,9001)
        ENDIF
C
C     Viscosite moleculaire aux cellules
        IOK2 = 0
        DO IPHAS = 1, NPHAS
          IF(IVIVAR(IPHAS).EQ.0) THEN
            IPCVIS = IPPROC(IVISCL(IPHAS))
            DO IEL = 1, NCEL
              IF( ABS(PROPCE(IEL ,IPCVIS)-VISCL0(IPHAS)).GT.EPZERO) THEN
                IOK2 = 1
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        IF(IOK2.NE.0) THEN
          IF ( IPPMOD(ICOMPF) .GE. 0 ) THEN
            WRITE(NFECRA,9003)
          ELSE
            WRITE(NFECRA,9002)
          ENDIF
        ENDIF
C
        IF(IOK1.NE.0.OR.IOK2.NE.0) THEN
          CALL CSEXIT(1)
        ENDIF
C
      ENDIF
C
C=======================================================================
C 3.  CALCUL DE LA VISCOSITE TURBULENTE
C=======================================================================
C
C --- Boucle sur les phases : Debut
      DO IPHAS = 1, NPHAS
C
        IF     (ITURB(IPHAS).EQ. 0) THEN
C
C 3.1 LAMINAIRE
C ==============
C
          IPCVST = IPPROC(IVISCT(IPHAS))
C
          DO IEL = 1, NCEL
            PROPCE(IEL,IPCVST) = 0.D0
          ENDDO
C
        ELSEIF (ITURB(IPHAS).EQ.10) THEN
C
C 3.2 LONGUEUR DE MELANGE
C ========================
C
          IPH = IPHAS
          CALL VISLMG
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS) , NCKPDC(IPHAS) , NCETSM(IPHAS) ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPH    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) , IA(IICESM(IPHAS)) , IA(IITPSM(IPHAS)) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)), RA(ISMACE(IPHAS)),
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C
        ELSEIF (ITYTUR(IPHAS).EQ.2) THEN
C
C 3.3 K-EPSILON
C ==============
C
          IKIPH  = IK(IPHAS)
          IEIPH  = IEP(IPHAS)
          IPCVST = IPPROC(IVISCT(IPHAS))
          IPCROM = IPPROC(IROM  (IPHAS))
C
          DO IEL = 1, NCEL
            XK = RTP(IEL,IKIPH)
            XE = RTP(IEL,IEIPH)
            PROPCE(IEL,IPCVST) = PROPCE(IEL,IPCROM)*CMU*XK**2/XE
          ENDDO
C
        ELSEIF (ITYTUR(IPHAS).EQ.3) THEN
C
C 3.4 Rij-EPSILON
C ================
C
          IR11IP = IR11(IPHAS)
          IR22IP = IR22(IPHAS)
          IR33IP = IR33(IPHAS)
          IEIPH  = IEP(IPHAS)
          IPCVST = IPPROC(IVISCT(IPHAS))
          IPCROM = IPPROC(IROM  (IPHAS))
C
          DO IEL = 1, NCEL
            XK = 0.5D0*(RTP(IEL,IR11IP)+RTP(IEL,IR22IP)+RTP(IEL,IR33IP))
            XE = RTP(IEL,IEIPH)
            PROPCE(IEL,IPCVST) = PROPCE(IEL,IPCROM)*CMU*XK**2/XE
          ENDDO
C
        ELSEIF (ITURB(IPHAS).EQ.40) THEN
C
C 3.4 LES Smagorinsky
C ===================
C
          IPH = IPHAS
C
          CALL VISSMA
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS) , NCKPDC(IPHAS) , NCETSM(IPHAS) ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPH    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) ,
     &   IA(IICESM(IPHAS)) , IA(IITPSM(IPHAS)),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)), RA(ISMACE(IPHAS)),
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C
        ELSEIF(ITURB(IPHAS).EQ.41) THEN
C
C 3.4 LES dynamique
C =================
C
          IPH = IPHAS
C
          CALL VISDYN
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS) , NCKPDC(IPHAS) , NCETSM(IPHAS) ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPH    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) ,
     &   IA(IICESM(IPHAS)) , IA(IITPSM(IPHAS)),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)), RA(ISMACE(IPHAS)),
     &   PROPCE(1,IPPROC(ISMAGO(IPHAS))) ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     , W9     , W10    , XMIJ   ,
     &   RDEVEL , RTUSER , RA     )
C
        ELSEIF (ITURB(IPHAS).EQ.50) THEN
C
C 3.5 v2f phi-model
C =================
          IKIPH  = IK(IPHAS)
          IEIPH  = IEP(IPHAS)
          IPHIPH = IPHI(IPHAS)
          IPCVIS = IPPROC(IVISCL(IPHAS))
          IPCVST = IPPROC(IVISCT(IPHAS))
          IPCROM = IPPROC(IROM  (IPHAS))
C
          DO IEL = 1, NCEL
            XK = RTP(IEL,IKIPH)
            XE = RTP(IEL,IEIPH)
            XROM = PROPCE(IEL,IPCROM)
            XNU = PROPCE(IEL,IPCVIS)/XROM
            TTKE = XK / XE
            TTMIN = CV2FCT*SQRT(XNU/XE)
            TT = MAX(TTKE,TTMIN)
            PROPCE(IEL,IPCVST) = CV2FMU*XROM*TT*RTP(IEL,IPHIPH)
     &           *RTP(IEL,IKIPH)
          ENDDO
C
        ELSEIF (ITURB(IPHAS).EQ.60) THEN
C
C 3.6 K-OMEGA SST
C ===============
C
          IPH = IPHAS
          CALL VISSST
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS) , NCKPDC(IPHAS) , NCETSM(IPHAS) ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPH    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) , IA(IICESM(IPHAS)) , IA(IITPSM(IPHAS)) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)), RA(ISMACE(IPHAS)),
     &   RA(IS2KW(IPHAS)), RA(IDVUKW(IPHAS)),
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF
C
      ENDDO
C --- Boucle sur les phases : fin
C
C=======================================================================
C 4.  MODIFICATION UTILISATEUR DE LA VISCOSITE TURBULENTE
C=======================================================================
C
      DO IPHAS = 1, NPHAS
C
        IPH = IPHAS
        CALL USVIST
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS)   , NCKPDC(IPHAS)   , NCETSM(IPHAS)   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPH    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) , IA(IICESM(IPHAS)) , IA(IITPSM(IPHAS)) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)), RA(ISMACE(IPHAS)),
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDDO
C
C=======================================================================
C 5.  CLIPPING DE LA VISCOSITE TURBULENTE EN LES DYNAMIQUE
C=======================================================================
C
      DO IPHAS = 1, NPHAS
C Pour la LES en modele dynamique on clippe la viscosite turbulente de maniere
C a ce que mu+mu_t soit positif, .e. on autorise mu_t legerement negatif
C La diffusivite turbulente des scalaires (mu_t/sigma), elle, sera clippee a 0
C dans covofi
C
        IF (ITURB(IPHAS).EQ.41) THEN
          IPCVIS = IPPROC(IVISCL(IPHAS))
          IPCVST = IPPROC(IVISCT(IPHAS))
          ICLIPC = 0
          DO IEL = 1, NCEL
            VISTOT = PROPCE(IEL,IPCVIS) + PROPCE(IEL,IPCVST)
            IF(VISTOT.LT.0.D0) THEN
              PROPCE(IEL,IPCVST) = 0.D0
              ICLIPC = ICLIPC + 1
            ENDIF
          ENDDO
          IF(IWARNI(IU(IPHAS)).GE.1) THEN
            IF(IRANGP.GE.0) THEN
              CALL PARCPT(ICLIPC)
C             ===========
            ENDIF
            WRITE(NFECRA,1000) ICLIPC
          ENDIF
        ENDIF
C
      ENDDO
C
C=======================================================================
C 6.  MODIFICATION UTILISATEUR DE LA VISCOSITE DE MAILLAGE EN ALE
C=======================================================================
C
      IF (IALE.EQ.1) THEN
C
        CALL USVIMA
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , PROPCE(1,IPPROC(IVISMA(1))) ,
     &   PROPCE(1,IPPROC(IVISMA(2))) , PROPCE(1,IPPROC(IVISMA(3))) ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C=======================================================================
C 7.  IMPRESSIONS DE CONTROLE DES VALEURS ENTREES PAR L'UTILISATEUR
C=======================================================================
C
C
C ---> Calcul des bornes des variables par phase et impressions
C
C Indicateur d'erreur
      IOK = 0
C
      DO IPHAS = 1, NPHAS
C
C Rang des variables dans PROPCE
        IPCROM = IPPROC(IROM(IPHAS))
        IPCVIS = IPPROC(IVISCL(IPHAS))
        IPCVST = IPPROC(IVISCT(IPHAS))
        IF(ICP(IPHAS).GT.0) THEN
          IPCCP  = IPPROC(ICP   (IPHAS))
          NN     = 4
        ELSE
          IPCCP = 0
          NN    = 3
        ENDIF
C
C Rang des variables dans PROPFB
        IPBROM = IPPROB(IROM(IPHAS))
C
C Min et max sur les cellules
        DO II = 1, NN
          IVAR = 0
          IF(II.EQ.1) IVAR = IPCROM
          IF(II.EQ.2) IVAR = IPCVIS
          IF(II.EQ.3) IVAR = IPCVST
          IF(II.EQ.4) IVAR = IPCCP
          IF(IVAR.GT.0) THEN
            VARMX(II) = PROPCE(1,IVAR)
            VARMN(II) = PROPCE(1,IVAR)
            DO IEL = 2, NCEL
              VARMX(II) = MAX(VARMX(II),PROPCE(IEL,IVAR))
              VARMN(II) = MIN(VARMN(II),PROPCE(IEL,IVAR))
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARMAX (VARMX(II))
C             ===========
              CALL PARMIN (VARMN(II))
C             ===========
            ENDIF
          ENDIF
        ENDDO
C
C Min et max sur les faces de bord (masse volumique uniquement)
        II   = 1
        IVAR = IPBROM
        DO IFAC = 1, NFABOR
          VARMX(II) = MAX(VARMX(II),PROPFB(IFAC,IVAR))
          VARMN(II) = MIN(VARMN(II),PROPFB(IFAC,IVAR))
        ENDDO
        IF (IRANGP.GE.0) THEN
          CALL PARMAX (VARMX(II))
C         ===========
          CALL PARMIN (VARMN(II))
C         ===========
        ENDIF
C
C Impressions
        IOK1 = 0
        DO II = 1, NN
          IF(II.EQ.1) CHAINE = NOMVAR(IPPPRO(IPPROC(IROM  (IPHAS))))
          IF(II.EQ.2) CHAINE = NOMVAR(IPPPRO(IPPROC(IVISCL(IPHAS))))
          IF(II.EQ.3) CHAINE = NOMVAR(IPPPRO(IPPROC(IVISCT(IPHAS))))
          IF(II.EQ.4) CHAINE = NOMVAR(IPPPRO(IPPROC(ICP   (IPHAS))))
          IF(IWARNI(IU(IPHAS)).GE.1.OR.IPASS.EQ.1.OR.
     &                                          VARMN(II).LT.0.D0)THEN
            IF(IOK1.EQ.0) THEN
              WRITE(NFECRA,3010)IPHAS
              IOK1 = 1
            ENDIF
            IF ((II.NE.3).OR.(ITURB(IPHAS).NE.0))
     &           WRITE(NFECRA,3011)CHAINE(1:8),VARMN(II),VARMX(II)
          ENDIF
        ENDDO
        IF(IOK1.EQ.1) WRITE(NFECRA,3012)
C
C Verifications de valeur physique
C
C Masse volumique definie
        II = 1
        CHAINE = NOMVAR(IPPPRO(IPPROC(IROM  (IPHAS))))
        IF (VARMN(II).LT.0.D0) THEN
          WRITE(NFECRA,9011)CHAINE(1:8),VARMN(II)
          IOK = IOK + 1
        ENDIF
C
C Viscosite moleculaire definie
        II = 2
        CHAINE = NOMVAR(IPPPRO(IPPROC(IVISCL(IPHAS))))
        IF (VARMN(II).LT.0.D0) THEN
          WRITE(NFECRA,9011)CHAINE(1:8),VARMN(II)
          IOK = IOK + 1
        ENDIF
C
C Viscosite turbulente definie
C on ne clippe pas mu_t en modele LES dynamique, car on a fait
C un clipping sur la viscosite totale
        II = 3
        CHAINE = NOMVAR(IPPPRO(IPPROC(IVISCT(IPHAS))))
        IF (VARMN(II).LT.0.D0.AND.ITURB(IPHAS).NE.41) THEN
          WRITE(NFECRA,9012)IPHAS,VARMN(II)
          IOK = IOK + 1
        ENDIF
C
C Chaleur specifique definie
        IF(ICP(IPHAS).GT.0) THEN
          II = 4
          CHAINE = NOMVAR(IPPPRO(IPPROC(ICP   (IPHAS))))
          IF (VARMN(II).LT.0.D0) THEN
            IISCT = 0
            DO ISCAL = 1, NSCAL
              IF (ISCSTH(ISCAL).NE.0) THEN
                IISCT = 1
              ENDIF
            ENDDO
            IF (IISCT.EQ.1) THEN
              WRITE(NFECRA,9011)CHAINE(1:8),VARMN(II)
              IOK = IOK + 1
            ENDIF
          ENDIF
        ENDIF
C
      ENDDO
C
C
C
C ---> Calcul des bornes des scalaires et impressions
C
      IF(NSCAL.GE.1) THEN
C
        IOK1 = 0
        DO ISCAL = 1, NSCAL
C
          IF(IVISLS(ISCAL).GT.0) THEN
            IPCVIS = IPPROC(IVISLS(ISCAL))
          ELSE
            IPCVIS = 0
          ENDIF
C
          VISMAX(ISCAL) = -GRAND
          VISMIN(ISCAL) =  GRAND
          IF(IPCVIS.GT.0) THEN
            DO IEL = 1, NCEL
              VISMAX(ISCAL) = MAX(VISMAX(ISCAL),PROPCE(IEL,IPCVIS))
              VISMIN(ISCAL) = MIN(VISMIN(ISCAL),PROPCE(IEL,IPCVIS))
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARMAX (VISMAX(ISCAL))
C             ===========
              CALL PARMIN (VISMIN(ISCAL))
C             ===========
            ENDIF
          ELSE
            VISMAX(ISCAL) = VISLS0(ISCAL)
            VISMIN(ISCAL) = VISLS0(ISCAL)
          ENDIF
C
          IVAR = ISCA(ISCAL)
          IF(IWARNI(IVAR).GE.1.OR.IPASS.EQ.1.OR.
     &                                        VISMIN(ISCAL).LE.0.D0)THEN
            CHAINE = NOMVAR(IPPRTP(IVAR))
            IF(IOK1.EQ.0) THEN
              WRITE(NFECRA,3110)
              IOK1 = 1
            ENDIF
            WRITE(NFECRA,3111)CHAINE(1:8),ISCAL,
     &                        VISMIN(ISCAL),VISMAX(ISCAL)
          ENDIF
C
        ENDDO
        IF(IOK1.EQ.1) WRITE(NFECRA,3112)
C
C Verifications de valeur physique
C
C IOK a deja ete initialise pour les valeurs liees aux phases
C
        DO ISCAL = 1, NSCAL
C
          IVAR = ISCA(ISCAL)
C
          IF (VISMIN(ISCAL).LT.0.D0) THEN
            CHAINE = NOMVAR(IPPRTP(IVAR))
            WRITE(NFECRA,9111)CHAINE(1:8),ISCAL,VISMIN(ISCAL)
            IOK = IOK + 1
          ENDIF
C
        ENDDO
C
      ENDIF
C
C ---> Calcul des bornes de viscosite de maillage en ALE
C
      IF(IALE.EQ.1) THEN
C
        IOK1 = 0
        NN = 1
        IF (IORTVM.EQ.1) NN = 3
        DO II = 1, NN
          IPCVMA = IPPROC(IVISMA(II))
C
C Min et max sur les cellules
          VARMX(1) = PROPCE(1,IPCVMA)
          VARMN(1) = PROPCE(1,IPCVMA)
          DO IEL = 2, NCEL
            VARMX(1) = MAX(VARMX(1),PROPCE(IEL,IPCVMA))
            VARMN(1) = MIN(VARMN(1),PROPCE(IEL,IPCVMA))
          ENDDO
          IF (IRANGP.GE.0) THEN
            CALL PARMAX (VARMX(1))
C           ===========
            CALL PARMIN (VARMN(1))
C           ===========
          ENDIF
C
C Impressions
          CHAINE = NOMVAR(IPPPRO(IPCVMA))
          IF(IWARNI(IUMA).GE.1.OR.IPASS.EQ.1.OR.
     &         VARMN(1).LT.0.D0)THEN
            IF (IOK1.EQ.0) THEN
              WRITE(NFECRA,3210)
              IOK1 = 1
            ENDIF
            WRITE(NFECRA,3211)CHAINE(1:8),VARMN(1),VARMX(1)
          ENDIF
C
C Verifications de valeur physique
C
C Viscosite de maillage definie
          CHAINE = NOMVAR(IPPPRO(IPCVMA))
          IF (VARMN(1).LE.0.D0) THEN
            WRITE(NFECRA,9011)CHAINE(1:8),VARMN(1)
            IOK = IOK + 1
          ENDIF
C
        ENDDO
C
        IF (IOK1.EQ.1) WRITE(NFECRA,3212)
C
      ENDIF
C
C --->  arret eventuel
C
      IF(IOK.NE.0) THEN
        WRITE(NFECRA,9999)IOK
        CALL CSEXIT (1)
      ENDIF
C
C=======================================================================
C 8.  ECHANGES
C=======================================================================
C
C Pour navsto et vissec on a besoin de ROM dans le halo
C
C
      DO IPHAS = 1, NPHAS
C
        IPCROM = IPPROC(IROM(IPHAS))
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (PROPCE(1,IPCROM))
C         ===========
        ENDIF
C
        IF(IPERIO.EQ.1) THEN
          IDIMTE = 0
          ITENSO = 0
          CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    PROPCE(1,IPCROM),PROPCE(1,IPCROM),PROPCE(1,IPCROM),
     &    PROPCE(1,IPCROM),PROPCE(1,IPCROM),PROPCE(1,IPCROM),
     &    PROPCE(1,IPCROM),PROPCE(1,IPCROM),PROPCE(1,IPCROM) )
        ENDIF
C
      ENDDO
C
C----
C FORMATS
C----
C
C
 1000 FORMAT(
     &' Nb de clippings de la viscosite totale (mu+mu_t>0) :',I10,/)
 3010 FORMAT(
     &' --- Phase : ',I10                                            ,/,
     &' ---------------------------------                           ',/,
     &' Propriete  Valeur min  Valeur max                           ',/,
     &' ---------------------------------                           '  )
 3011 FORMAT(
     & 2X,    A8,      E12.4,      E12.4                               )
 3012 FORMAT(
     &' ---------------------------------                           ',/)
 3110 FORMAT(
     &' --- Diffusivite :                                           ',/,
     &' ---------------------------------------                     ',/,
     &' Scalaire Numero  Valeur min  Valeur max                     ',/,
     &' ---------------------------------------                     '  )
 3111 FORMAT(
     & 1X,    A8,   I7,      E12.4,      E12.4                         )
 3112 FORMAT(
     &' ---------------------------------------                     ',/)
 3210 FORMAT(
     &' --- Viscosite de maillage (methode ALE)                     ',/,
     &' ---------------------------------                           ',/,
     &' Propriete  Valeur min  Valeur max                           ',/,
     &' ---------------------------------                           '  )
 3211 FORMAT(
     & 2X,    A8,      E12.4,      E12.4                               )
 3212 FORMAT(
     &' ---------------------------------                           ',/)
C
 9001  FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,
     &'@    *********                                               ',/,
     &'@    INCOHERENCE ENTRE USPHYV ET USINI1 POUR                 ',/,
     &'@                                          LA MASSE VOLUMIQUE',/,
     &'@                                                            ',/,
     &'@  On a indique dans usini1 que la masse volumique etait     ',/,
     &'@     constante (IROVAR(IPHAS)=0) mais on a modifie ses      ',/,
     &'@     valeurs aux cellules ou aux faces de bord dans usphyv. ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1 et usphyv.                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9002  FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,
     &'@    *********                                               ',/,
     &'@    INCOHERENCE ENTRE USPHYV ET USINI1 POUR                 ',/,
     &'@                                    LA VISCOSITE MOLECULAIRE',/,
     &'@                                                            ',/,
     &'@  On a indique dans usini1 que la viscosite moleculaire     ',/,
     &'@     etait constante (IVIVAR(IPHAS)=0) mais on a modifie ses',/,
     &'@     valeurs aux cellules dans usphyv.                      ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1 et usphyv.                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9003  FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,
     &'@    *********                                               ',/,
     &'@    MODULE COMPRESSIBLE                                     ',/,
     &'@    INCOHERENCE ENTRE USCFPV ET USCFX1 POUR                 ',/,
     &'@                                    LA VISCOSITE MOLECULAIRE',/,
     &'@                                                            ',/,
     &'@  En compressible la viscosite moleculaire est constante par',/,
     &'@     defaut (IVIVAR(IPHAS)=0) et la valeur de IVIVAR n''a   ',/,
     &'@     pas ete modifiee dans uscfx1. Pourtant, on a modifie   ',/,
     &'@     les valeurs de la viscosite moleculaire dans uscfpv.   ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier uscfx1 et uscfpv.                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9011  FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,
     &'@    *********                                               ',/,
     &'@    LA PROPRIETE PHYSIQUE ',A8  ,' N A PAS ETE              ',/,
     &'@                                       CORRECTEMENT DEFINIE.',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  La propriete physique identifiee ci-dessus est variable et',/,
     &'@    le minimum atteint est ',E12.4                           ,/,
     &'@  Verifier que cette propriete a ete definie dans usphyv et ',/,
     &'@    que la loi adoptee conduit a des valeurs correctes.     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9012  FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,
     &'@    *********                                               ',/,
     &'@    LA VISCOSITE TURBULENTE DE LA PHASE ',I10                ,/,
     &'@                           N A PAS ETE CORRECTEMENT DEFINIE.',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Le minimum atteint est ',E12.4                             ,/,
     &'@  Verifier le cas echeant la definition de la masse         ',/,
     &'@    volumique dans usphyv et la modification de la viscosite',/,
     &'@    turbulente dans usvist.                                 ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9111  FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,
     &'@    *********                                               ',/,
     &'@    LA DIFFUSIVITE DU SCALAIRE ',A8                          ,/,
     &'@       (SCALAIRE NUMERO ',I10   ,') N A PAS ETE             ',/,
     &'@                                       CORRECTEMENT DEFINIE.',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  La propriete physique identifiee ci-dessus est variable et',/,
     &'@    le minimum atteint est ',E12.4                           ,/,
     &'@  Verifier que cette propriete a ete definie dans usphyv et ',/,
     &'@    que la loi adoptee conduit a des valeurs correctes.     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9999 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,
     &'@    *********                                               ',/,
     &'@    DES GRANDEURS PHYSIQUES ONT DES VALEURS INCORRECTES     ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,
     &'@                                                            ',/,
     &'@  Se reporter aux impressions precedentes pour plus de      ',/,
     &'@    renseignements.                                         ',/,
     &'@  Verifier usphyv.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      RETURN
      END
c@z
