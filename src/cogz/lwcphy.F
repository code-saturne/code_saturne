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
                        SUBROUTINE LWCPHY
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IBROM  , IZFPPP ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC ROUTINE PHYSIQUE PARTICULIERE : FLAMME DE PREMELANGE MODELE LWC
CFONC  Calcul de RHO adiabatique ou permeatique (transport de H)
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
CARGU ! NPHMX        ! E  !  -> ! NPHSMX                               !
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
CARGU ! IBROM        ! TE !  -> ! INDICATEUR DE REMPLISSAGE DE ROMB    !
CARGU !   (NPHMX   ) !    !     !                                      !
CARGU ! IZFPPP       ! TE ! <-  ! NUMERO DE ZONE DE LA FACE DE BORD    !
CARGU ! (NFABOR)     !    !     !  POUR LE MODULE PHYS. PART.          !
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
CARGU ! W1...8(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR), IBROM(NPHMX)
      INTEGER          IZFPPP(NFABOR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IGG, IPHAS, IEL, IPCROM
      INTEGER          IZONE , IFAC, IPBROM
      DOUBLE PRECISION COEFG(NGAZGM)
      DOUBLE PRECISION NBMOL , TEMSMM
      DOUBLE PRECISION MASMG
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
C=======================================================================
C 1. INITIALISATIONS A CONSERVER
C=======================================================================
C
C --- Initialisation memoire
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C ---> Initialisation
C
      DO IGG = 1, NGAZGM
        COEFG(IGG) = ZERO
      ENDDO
C
C ---> Positions des variables, coefficients
C
      IPHAS = 1
      IPCROM = IPPROC(IROM(IPHAS))
      IPBROM = IPPROB(IROM(IPHAS))
C
C=======================================================================
C 2. DETERMINATION DES GRANDEURS THERMOCHIMIQUES MOYENNES
C=======================================================================
C
C
      IF ( (IPPMOD(ICOLWC).EQ.0) .OR. (IPPMOD(ICOLWC).EQ.1) ) THEN
C
        CALL PDFLWC
C       ===========
     &   ( NCELET        , NCEL          ,
     &     RTP(1,ISCA(IFM))    , RTP(1,ISCA(IFP2M))  ,
     &     RTP(1,ISCA(IYFM))   , RTP(1,ISCA(IYFP2M)) ,
     &     PROPCE   )
C
      ENDIF
C
       IF ( (IPPMOD(ICOLWC).EQ.2) .OR. (IPPMOD(ICOLWC).EQ.3) ) THEN
C
        CALL PDFPP3
C       ===========
     &   ( NCELET        , NCEL          ,
     &     RTP(1,ISCA(IFM))    , RTP(1,ISCA(IFP2M))  ,
     &     RTP(1,ISCA(IYFM))   , RTP(1,ISCA(IYFP2M)) ,
     &     RTP(1,ISCA(ICOYFP)) ,
     &     PROPCE   )
C
       ENDIF
C
       IF ( (IPPMOD(ICOLWC).EQ.4).OR.(IPPMOD(ICOLWC).EQ.5) ) THEN
C
        CALL PDFPP4
C       ===========
     &   ( NCELET        , NCEL          ,
     &     RTP(1,ISCA(IFM))    , RTP(1,ISCA(IFP2M))  ,
     &     RTP(1,ISCA(IYFM))   , RTP(1,ISCA(IYFP2M)) ,
     &     RTP(1,ISCA(ICOYFP)) ,
     &     PROPCE   )
C
       ENDIF
C
C=======================================================================
C 3. CALCUL DE RHO ET DES FRACTIONS MASSIQUES DES ESPECES GLOBALES
C    SUR LES BORDS
C=======================================================================
C
C --> Masse volumique au bord
C
      IPHAS = 1
      IBROM(IPHAS) = 1
C
C ---- Masse volumique au bord pour toutes les facettes
C      Les facettes d'entree seront recalculees.
C
      DO IFAC = 1, NFABOR
      IEL = IFABOR(IFAC)
        PROPFB(IFAC,IPBROM) = PROPCE(IEL,IPCROM)
      ENDDO
C
C ---- Masse volumique au bord pour les facettes d'entree UNIQUEMENT
C      Le test sur IZONE sert pour les reprises de calcul
C
      IF ( IPASS.GT.1 .OR. ISUITE.EQ.1 ) THEN
        DO IFAC = 1, NFABOR
          IZONE = IZFPPP(IFAC)
          IF(IZONE.GT.0) THEN
            IF ( IENTGB(IZONE).EQ.1 .OR. IENTGF(IZONE).EQ. 1) THEN
              COEFG(1) = FMENT(IZONE)
              COEFG(2) = 1.D0-FMENT(IZONE)
              COEFG(3) = ZERO
              IF ( IENTGB(IZONE).EQ.1 ) THEN
                COEFG(1) = MAX(ZERO,(FMENT(IZONE)-FS(1))/(1.D0-FS(1)))
                COEFG(3) = (FMENT(IZONE)-COEFG(1))/FS(1)
                COEFG(2) = 1.D0 - COEFG(1) - COEFG(3)
              ENDIF
              NBMOL = 0.D0
              DO IGG = 1, NGAZG
                NBMOL = NBMOL + COEFG(IGG)/WMOLG(IGG)
              ENDDO
             MASMG = 1.D0/NBMOL
             TEMSMM = TKENT(IZONE)/MASMG
             PROPFB(IFAC,IPBROM) = P0(IPHAS)/(RR*TEMSMM)
            ENDIF
          ENDIF
        ENDDO
      ENDIF

C----
C FIN
C----
C
      RETURN
      END
c@z
