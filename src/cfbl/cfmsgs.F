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
                        SUBROUTINE CFMSGS
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCAL  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   FLUMAS , FLUMAB , FLABGS , FLBBGS ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , W10    , W11    , W12    ,
     &   TRFLMS , TRFLMB , COEFU  , XAM    ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC CALCUL DU "FLUX DE MASSE" AUX FACES
CFONC    POUR LA RESOLUTION DE LA MASSE VOLUMIQUE
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
CARGU ! ISCAL        ! E  !  -> ! NUMERO DU SCALAIRE                   !
CARGU ! ITSPDV       ! E  !  -> ! CALCUL TERMES SOURCES PROD ET DISSIP !
CARGU !              !    !     !  (0 : NON , 1 : OUI)                 !
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
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! TSLAGR       ! TR !  -> ! TERME DE COUPLAGE RETOUR DU          !
CARGU !(NCELET,*)    !    !     !     LAGRANGIEN                       !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! FLUMAS(NFAC) ! TR ! <-  ! FLUX DE MASSE AUX FACES INTERNES     !
CARGU ! FLUMAB(NFABOR! TR ! <-  ! FLUX DE MASSE AUX FACES DE BORD      !
CARGU ! W1..12(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! TRFLMS(NFAC) ! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! TRFLMB(NFABOR! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! COEFU(NFABO,3! TR !  -  ! TABLEAU DE TRAVAIL CL DE LA QDM      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
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
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
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
      INTEGER          NCEPDP , NCKPDP , NCESMP
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          ISCAL
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICEPDC(NCEPDP)
      INTEGER          ICETSM(NCESMP), ITYPSM(NCESMP,NVAR)
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
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP), SMACEL(NCESMP,NVAR)
      DOUBLE PRECISION FLUMAS(NFAC), FLUMAB(NFABOR)
      DOUBLE PRECISION FLABGS(NFAC), FLBBGS(NFABOR)
      DOUBLE PRECISION W1(NCELET) , W2(NCELET) , W3(NCELET)
      DOUBLE PRECISION W4(NCELET) , W5(NCELET) , W6(NCELET)
      DOUBLE PRECISION W7(NCELET) , W8(NCELET) , W9(NCELET)
      DOUBLE PRECISION W10(NCELET), W11(NCELET), W12(NCELET)
      DOUBLE PRECISION TRFLMS(NFAC), TRFLMB(NFABOR)
      DOUBLE PRECISION COEFU(NFABOR,3), XAM(NFAC,2)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA, IFINIA
      INTEGER          IVAR  , IPHAS
      INTEGER          IFAC  , IEL
      INTEGER          INIT  , INC   , ICCOCG, II, JJ
      INTEGER          IPP
      INTEGER          NSWRGP, IMLIGP, IWARNP
      INTEGER          ICONVP, IDIFFP
      INTEGER          IRCFLP, ISCHCP, ISSTPP
C
      INTEGER          IIROM , IIROMB
      INTEGER          IVAR0 , IMVIS1, ICCFTH, IMODIF, ISOU
      INTEGER          IMASPE, IFLMB0, IISMPH
      INTEGER          ICLIUP, ICLIVP, ICLIWP, ICLVAR
      INTEGER          IDIMTE, ITENSO
      INTEGER          IUIPH , IVIPH , IWIPH
      INTEGER          ITSQDM, IIUN  , IEXTTS
      INTEGER          MAXELT, ILS
C
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP, BLENCP
      DOUBLE PRECISION FLUI  , FLUJ  , PFAC  , THETV
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
C --- Numero de phase associee au scalaire traite
      IPHAS  = IPHSCA(ISCAL)
C
C --- Numero des variables de calcul
C     Masse volumique
      IVAR   = ISCA(ISCAL)
C     Vitesses
      IUIPH  = IU(IPHAS)
      IVIPH  = IV(IPHAS)
      IWIPH  = IW(IPHAS)
C
C     Masse volumique dans PROPCE
      IIROM  = IPPROC(IROM  (IPHAS))
      IIROMB = IPPROB(IROM  (IPHAS))
C
C ---> Initialisation du flux de masse
C
      DO IFAC = 1, NFAC
        FLUMAS(IFAC) = 0.D0
      ENDDO
      DO IFAC = 1, NFABOR
        FLUMAB(IFAC) = 0.D0
      ENDDO
      DO IFAC = 1, NFAC
        FLABGS(IFAC) = 0.D0
      ENDDO
      DO IFAC = 1, NFABOR
        FLBBGS(IFAC) = 0.D0
      ENDDO
C
C=======================================================================
C 2. FLUX DE MASSE AUX FACES
C=======================================================================
C
C     2.1 TERMES SOURCES DE L'EQUATION DE QDM
C     =======================================
C
C     FX -> W5 , FY -> W6 , FZ -> W7
C
C     Apres premiers tests, il semble (tube à choc en double détente)
C       que ce ne soit pas une bonne idée de prendre en compte dans
C       l'équation de la masse tous les termes de l'équation de la
C       quantité de mouvement (en particulier le terme convectif, mais
C       les termes de diffusion, en gradient transposé, source de masse
C       et utilisateur n'ont pas été vraiment testés, dans la mesure ou
C       ils étaient nuls pour tous les cas unitaires considérés).
C       Bien entendu, s'agissant de tests préliminaires, il est possible
C       que le comportement insatisfaisant observé soit le fruit d'une
C       simple erreur de codage (mais on ne l'a pas trouvée...).
C     On propose donc par sécurité de ne pas prendre en compte les termes
C       source de l'équation de la quantité de mouvement, hormis le terme
C       de gravité (car il contrebalance le gradient de pression et son
C       effet est bien visible lorsqu'on est dans une situation
C       d'équilibre).
C     Cependant, pour l'avenir, on conserve ici le codage
C       préliminaire qui a permis d'effectuer les tests (version 1.1.0.h).
C       On l'encapsule dans un test qui le désactive obligatoirement
C       (ainsi, pas de question utilisateur intempestive, et toujours
C       la possibilité de poursuivre les tests si l'on a le temps  ou si
C       on rencontre des pbs...).
C     Noter que, avec ces termes, le système devient difficile à analyser
C       sur le papier (alors que sans ces termes, on est dans la
C       configuration Euler + gravité).
C
C --- Initialisation
      DO IEL = 1, NCEL
        W10(IEL) = 0.D0
      ENDDO
      DO IEL = 1, NCEL
        W11(IEL) = 0.D0
      ENDDO
      DO IEL = 1, NCEL
        W12(IEL) = 0.D0
      ENDDO
C
C
C     Test sur les termes source de qdm
      ITSQDM = 0
      IF(ITSQDM.NE.0) THEN
C
C
C --- Terme source utilisateur
C
      MAXELT = MAX(NCELET, NFAC, NFABOR)
      ILS    = IDEBIA
      IFINIA = ILS + MAXELT
      CALL IASIZE('CFMSGS',IFINIA)
C
C     Suivant X
        CALL USTSNS
C       ===========
     & ( IFINIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IUIPH  , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   W10    , W9     ,
C        ------   ------
     &   W8     , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C     Suivant Y
        CALL USTSNS
C       ===========
     & ( IFINIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVIPH  , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   W11    , W9     ,
C        ------   ------
     &   W8     , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C     Suivant Z
        CALL USTSNS
C       ===========
     & ( IFINIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IWIPH  , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   W12    , W9     ,
C        ------   ------
     &   W8     , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C --- Terme de convection de quantite de mouvement
        IF(ICONV(IUIPH).GE.1) THEN
C
          ICLIUP = ICLRTP(IUIPH ,ICOEF)
          ICLIVP = ICLRTP(IVIPH ,ICOEF)
          ICLIWP = ICLRTP(IWIPH ,ICOEF)
C
          INIT   = 1
          INC    = 1
          ICCOCG = 1
          IFLMB0 = 1
          IISMPH = IISYMP+NFABOR*(IPHAS-1)
          NSWRGP = NSWRGR(IUIPH)
          IMLIGP = IMLIGR(IUIPH)
          IWARNP = IWARNI(IUIPH)
          EPSRGP = EPSRGR(IUIPH)
          CLIMGP = CLIMGR(IUIPH)
          EXTRAP = EXTRAG(IUIPH)
C
          IMASPE = 1
C
C     Calcul du flux de masse
          CALL INIMAS
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IUIPH  , IVIPH  , IWIPH  , IMASPE , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFLMB0 , INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IISMPH) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PROPCE(1,IIROM) , PROPFB(1,IIROMB),
     &   RTPA (1,IUIPH)  , RTPA (1,IVIPH)  , RTPA (1,IWIPH)  ,
     &   COEFA(1,ICLIUP) , COEFA(1,ICLIVP) , COEFA(1,ICLIWP) ,
     &   COEFB(1,ICLIUP) , COEFB(1,ICLIVP) , COEFB(1,ICLIWP) ,
     &   FLUMAS , FLUMAB ,
C        ------   ------
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C
C     Calcul du terme convecte suivant les 3 directions
C       sans reconstruction
          DO ISOU = 1, 3
            IF(ISOU.EQ.1) IVAR0  = IUIPH
            IF(ISOU.EQ.1) ICLVAR = ICLIUP
            IF(ISOU.EQ.2) IVAR0  = IVIPH
            IF(ISOU.EQ.2) ICLVAR = ICLIVP
            IF(ISOU.EQ.3) IVAR0  = IWIPH
            IF(ISOU.EQ.3) ICLVAR = ICLIWP
C
            DO IFAC = 1, NFAC
              II = IFACEL(1,IFAC)
              JJ = IFACEL(2,IFAC)
              FLUI = 0.5D0*( FLUMAS(IFAC) +ABS(FLUMAS(IFAC)) )
              FLUJ = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
              TRFLMS(IFAC) = -(FLUI*RTPA(II,IVAR0)+FLUJ*RTPA(JJ,IVAR0))
            ENDDO
C
            DO IFAC = 1, NFABOR
              II = IFABOR(IFAC)
              FLUI = 0.5D0*( FLUMAB(IFAC) +ABS(FLUMAB(IFAC)) )
              FLUJ = 0.5D0*( FLUMAB(IFAC) -ABS(FLUMAB(IFAC)) )
              PFAC = COEFA(IFAC,ICLVAR)
     &             + COEFB(IFAC,ICLVAR)*RTPA(II,IVAR0)
              TRFLMB(IFAC) = - ( FLUI*RTPA(II,IVAR0) + FLUJ*PFAC )
            ENDDO
C
            INIT = 0
            IF(ISOU.EQ.1) THEN
              CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &             IFACEL,IFABOR,TRFLMS,TRFLMB,W10)
            ELSEIF(ISOU.EQ.2) THEN
              CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &             IFACEL,IFABOR,TRFLMS,TRFLMB,W11)
            ELSEIF(ISOU.EQ.3) THEN
              CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &             IFACEL,IFABOR,TRFLMS,TRFLMB,W12)
            ENDIF
C
          ENDDO
C
        ENDIF
C
C
C --- Terme de viscosite
C
        IF( IDIFF(IUIPH).GE.1 ) THEN
C
          DO IEL = 1, NCELET
            W8(IEL) = 1.D0
            W9(IEL) = 0.D0
          ENDDO
C
          CALL CFDIVS
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   W10    , W8     , W9     , W9     ,
C        ------
     &   W7     ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
          CALL CFDIVS
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   W11    , W9     , W8     , W9     ,
C        ------
     &   W7     ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
          CALL CFDIVS
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   W12    , W9     , W9     , W8     ,
C        ------
     &   W7     ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF
C
C
C --- Terme source de masse
C     On met tout en explicite pour l'instant ... a modifier lors des
C      tests sur ce terme (cf. remarque en debut de fichier)
        IIUN   = 1
        IEXTTS = 0
        THETV  = 0.D0
        DO ISOU = 1, 3
          IF(ISOU.EQ.1) THEN
            IVAR0  = IUIPH
            CALL CATSMA
C           ===========
     & ( NCELET, NCEL   , NCESMP , IIUN   , IEXTTS , THETV  ,
     &   ICETSM, ITYPSM(1,IVAR0) , VOLUME , RTPA(1,IVAR0)   ,
     &   SMACEL(1,IVAR0), SMACEL(1,IPR(IPHAS)),
     &   W10   , W1     , W2 )
            DO IEL = 1, NCEL
              W10(IEL) = W10(IEL) + W2(IEL)
            ENDDO
C
          ELSEIF(ISOU.EQ.2) THEN
            IVAR0  = IVIPH
            CALL CATSMA
C           ===========
     & ( NCELET, NCEL   , NCESMP , IIUN   , IEXTTS , THETV  ,
     &   ICETSM, ITYPSM(1,IVAR0) , VOLUME , RTPA(1,IVAR0)   ,
     &   SMACEL(1,IVAR0), SMACEL(1,IPR(IPHAS)),
     &   W11   , W1     , W2 )
            DO IEL = 1, NCEL
              W11(IEL) = W11(IEL) + W2(IEL)
            ENDDO
C
          ELSEIF(ISOU.EQ.3) THEN
            IVAR0  = IWIPH
            CALL CATSMA
C           ===========
     & ( NCELET, NCEL   , NCESMP , IIUN   , IEXTTS , THETV  ,
     &   ICETSM, ITYPSM(1,IVAR0) , VOLUME , RTPA(1,IVAR0)   ,
     &   SMACEL(1,IVAR0), SMACEL(1,IPR(IPHAS)),
     &   W12   , W1     , W2 )
            DO IEL = 1, NCEL
              W12(IEL) = W12(IEL) + W2(IEL)
            ENDDO
C
          ENDIF
C
        ENDDO
C
      ENDIF
C     Fin du Test sur les termes source de qdm
C
C
C --- Terme de forces volumiques (gravite)
      DO IEL = 1, NCEL
        W5(IEL) = GX + W10(IEL)/RTPA(IEL,IVAR)
        W6(IEL) = GY + W11(IEL)/RTPA(IEL,IVAR)
        W7(IEL) = GZ + W12(IEL)/RTPA(IEL,IVAR)
      ENDDO
C
C
C 2.2 CALCUL DU TERME DE DIFFUSION D'ENTROPIE
C ===========================================
C
C --- Calcul de l'entropie au centre des cellules et affectation a W1
      ICCFTH = 6
      IMODIF = 0
      CALL USCFTH
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IMODIF , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W8     , W9     , W10    ,
C        ------
     &   RDEVEL , RTUSER , RA     )
C
C --- Communication de l'entropie
      IF(IRANGP.GE.0) THEN
        CALL PARCOM (W1)
C       ===========
      ENDIF
C
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     & ( IDIMTE , ITENSO ,
     &   W1     , W1     , W1     ,
     &   W1     , W1     , W1     ,
     &   W1     , W1     , W1     )
      ENDIF
C
C --- Calcul de dt*Beta/Rho au centre des cellules et affectation a W2
      ICCFTH = 162
      IMODIF = 0
      CALL USCFTH
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IMODIF , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W2     , W8     , W9     , W10    ,
C        ------
     &   RDEVEL , RTUSER , RA     )
C
C --- Pour la condition au bord sur l'entropie
C     COEFA=COEFA(.,ITEMPK) et COEFB=COEFB(.,ITEMPK)
C
C
      DO IEL = 1, NCEL
        W2(IEL) =  DT(IEL) * W2(IEL)
      ENDDO
C
C --- Calcul de dt*Beta aux faces et affectation a TRFLMS et TRFLMB
      IMVIS1 = 1
C
      CALL VISCFA
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IMVIS1 ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W2     ,
     &   TRFLMS , TRFLMB ,
C        ------   ------
     &   RDEVEL , RTUSER , RA     )
C
C --- Calcul du flux de diffusion
C
C     Conditions aux limites de flux nul pour l'entropie
C       on utilise COEFU comme tableau de travail)
      DO IFAC = 1, NFABOR
        COEFU(IFAC,1) = 0.D0
        COEFU(IFAC,2) = 1.D0
      ENDDO
C
      INC =1
      ICCOCG = 1
      IVAR0  = 0
      ICONVP = 0
      IDIFFP = 1
      NSWRGP = NSWRGR(IVAR)
      IMLIGP = IMLIGR(IVAR)
      IRCFLP = IRCFLU(IVAR)
      ISCHCP = ISCHCV(IVAR)
      ISSTPP = ISSTPC(IVAR)
      IPP    = IPPRTP(IVAR)
      IWARNP = IWARNI(IVAR)
      BLENCP = BLENCV(IVAR)
      EPSRGP = EPSRGR(IVAR)
      CLIMGP = CLIMGR(IVAR)
      EXTRAP = EXTRAG(IVAR)
C
      CALL CFBSC3
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR0  , ICONVP , IDIFFP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , INC    , IMRGRA , ICCOCG ,
     &   IPP    , IWARNP ,
     &   BLENCP , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , COEFU(1,1)      , COEFU(1,2)      ,
     &            COEFU(1,1)      , COEFU(1,2)      ,
     &   TRFLMS , TRFLMB , TRFLMS , TRFLMB ,
     &   FLABGS , FLBBGS ,
C        ------   ------
     &   W2     , W3     , W4     , W8     , W9     , W10    ,
     &   RDEVEL , RTUSER , RA     )
C
C
C 2.3 CALCUL DU FLUX DE MASSE AUX FACES
C =====================================
C
C --- Calcul des "vitesses" de convection au centre des cellules
C
      DO IEL = 1, NCEL
        W10(IEL) = RTPA(IEL,IUIPH) + DT(IEL)*W5(IEL)
        W11(IEL) = RTPA(IEL,IVIPH) + DT(IEL)*W6(IEL)
        W12(IEL) = RTPA(IEL,IWIPH) + DT(IEL)*W7(IEL)
      ENDDO
C
C --- Calcul du flux par appel a INIMAS
C
C     Pour éviter d'appliquer une condition limite inadaptée, on impose
C       un simple flux nul. Noter que cela ne sert qu'a reconstruire des
C       gradients. La valeur au bord importe peu, dans la mesure ou le
C       flux de bord est ecrasé ensuite.
C
C     Pour éviter d'avoir recours a d'autres tableaux, on simule un flux
C       nul avec :
C       TRFLMB = 1 pour COEFB = 1 et
C       INC    = 0 pour COEFA = 0
C
C     On prend ROM = W1 = 1 et ROMB = TRFLMB = 1
C       (et ROMB sert aussi pour COEFB=1)
      DO IEL = 1, NCEL
        W1(IEL) = 1.D0
      ENDDO
      DO IFAC = 1, NFABOR
        TRFLMB(IFAC) = 1.D0
      ENDDO
C
      ICLIUP = ICLRTP(IUIPH ,ICOEF)
      ICLIVP = ICLRTP(IVIPH ,ICOEF)
      ICLIWP = ICLRTP(IWIPH ,ICOEF)
C
      INIT   = 1
      INC    = 0
C              ^ Comme indiqué ci-dessus, pour un flux nul
      ICCOCG = 1
      IVAR0  = 0
      IFLMB0 = 1
      IISMPH = IISYMP+NFABOR*(IPHAS-1)
      NSWRGP = NSWRGR(IVAR)
      IMLIGP = IMLIGR(IVAR)
      IWARNP = IWARNI(IVAR)
      EPSRGP = EPSRGR(IVAR)
      CLIMGP = CLIMGR(IVAR)
      EXTRAP = EXTRAG(IVAR)
C
      IMASPE = 1
C
      CALL INIMAS
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IVAR0  , IVAR0  , IVAR0  , IMASPE , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFLMB0 , INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IISMPH) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , TRFLMB ,
     &   W10    , W11    , W12    ,
     &   COEFA(1,ICLIUP) , COEFA(1,ICLIVP) , COEFA(1,ICLIWP) ,
     &   TRFLMB          , TRFLMB          , TRFLMB          ,
     &   FLUMAS , FLUMAB ,
C        ------   ------
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C
C--------
C FORMATS
C--------
C
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
