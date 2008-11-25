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
                        SUBROUTINE SCALAI
C                       *****************
C     --------------------------------------------------------------
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
     &   TSLAGR , COEFA  , COEFB  ,
     &   DTR    , VISCF  , VISCB  ,
     &   DAM    , XAM    ,
     &   DRTP   , SMBRS  , ROVSDT ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
CFONC    POUR LES SCALAIRES SUR UN PAS DE TEMPS
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
CARGU ! DTR(NCELET)  ! TR !  -  ! DT*CDTVAR                            !
CARGU ! VISCF(NFAC)  ! TR !  -  ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DRTP(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! SMBRS(NCELET ! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! ROVSDT(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! W1...9(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
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
      DOUBLE PRECISION TSLAGR(NCELET,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION DTR(NCELET)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET), SMBRS(NCELET)
      DOUBLE PRECISION ROVSDT(NCELET)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA, IFINIA
      INTEGER          ISCAL, IVAR, IPHAS, IEL
      INTEGER          II, IISC, ITSPDV, ICALC, IAPPEL
      INTEGER          ISPECF
      INTEGER          MAXELT, ILS
C
C NOMBRE DE PASSAGES DANS LA ROUTINE
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATIONS
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IPASS = IPASS + 1
C
C=======================================================================
C 2. TRAITEMENT DES SCALAIRES A PHYSIQUE PARTICULIERE
C    Cette section sera dediee aux traitement particuliers.
C    On donne un exemple qui n'est que la recopie (a ISCAPP pres) de
C      la section 3 sur les scalaires utilisateurs.
C=======================================================================
C
      IF (IPPMOD(IPHPAR).GE.1) THEN
C
C ---> Initialisation des RTP  a partir des CL
C
C   On initialise RTP (et pas RTPA) car RTPA ne sert pas
C      dans codits.
C
      MAXELT = MAX(NCELET, NFAC, NFABOR)
      ILS    = IDEBIA
      IFINIA = ILS + MAXELT
      CALL IASIZE('SCALAI',IFINIA)
C
        CALL PPINIV
C       ===========
     & ( IFINIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB , COEFA  , COEFB  ,
     &   RDEVEL , RTUSER , RA     )
C
C     On pourra eviter les bugs en initialisant aussi RTPA (sur NCELET)
C     (d'ailleurs, on s'en sert ensuite dans le cpflux etc)
        IF(IPASS.EQ.1.AND.ISUITE.EQ.0) THEN
          IF(NSCAPP.GT.0) THEN
            DO II = 1, NSCAPP
              ISCAL = ISCAPP(II)
              IVAR  = ISCA(ISCAL)
              DO IEL = 1, NCELET
                RTPA (IEL,IVAR) = RTP (IEL,IVAR)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
C
C ---> Calculs TS relatifs a la physique du charbon
C      GMDEV1, GMDEV2, GMHET, GMDCH
C
        IF ( IPPMOD(ICP3PL).NE.-1 ) THEN
C
          CALL CPFLUX
C         ===========
     &   ( IDEBIA , IDEBRA , NCELET , NCEL   ,
     &     RTPA   , PROPCE , VOLUME ,
     &     W1     , W2     , W3     ,
     &     RA     )
C
        ENDIF
C
C ---> Calculs TS relatifs a la physique du fuel
C      GAMEVA, GAMHTF
C
        IF ( IPPMOD(ICFUEL).NE.-1 ) THEN
C
          CALL FUFLUX
C         ===========
     &   ( IDEBIA , IDEBRA , NCELET , NCEL   ,
     &     RTPA   , PROPCE , VOLUME ,
     &     W1     , W2     , W3     ,
     &     RA     )
C
        ENDIF
C
C    ATTENTION : POUR LE CLIPPING AVEC ICLP = 1, IL FAUT
C
C                QUE LES SCALAIRES SOIENT RESOLUS AVANT
C                LEUR VARIANCE ASSOCIEE
C
C
C
C
C
C
C ---> Boucle sur les scalaires physiques particulieres.
C      On peut imaginer a la place des resolutions couplees.
C      Ici, on ne donne qu'un exemple.
C
        DO II = 1, NSCAPP
C
          ISCAL = ISCAPP(II)
          IVAR  = ISCA(ISCAL)
C
C ---> Pas de temps (avec facteur multiplicatif eventuel)
C
          IF(CDTVAR(IVAR).NE.1.D0) THEN
            DO IEL = 1, NCEL
              DTR(IEL) = DT(IEL)*CDTVAR(IVAR)
            ENDDO
          ELSE
            DO IEL = 1, NCEL
              DTR(IEL) = DT(IEL)
            ENDDO
          ENDIF
C
C     Schema compressible sans choc :
C ---> Traitement special pour la masse volumique,
C                     la temperature et l'energie
C     L'indicateur ISPECF sera non nul si l'on ne doit pas resoudre
C       le scalaire plus bas avec covofi.
C
          ISPECF = 0
C
          IF ( IPPMOD(ICOMPF).GE.0 ) THEN
C
            DO IPHAS = 1, NPHAS
              IF ( ISCAL.EQ.IRHO(IPHAS) .OR.
     &             ISCAL.EQ.ITEMPK(IPHAS) ) THEN
                ISPECF = 1
              ELSEIF ( ISCAL.EQ.IENERG(IPHAS) ) THEN
                ISPECF = 2
              ENDIF
            ENDDO
C
C ---> Masse volumique : deja resolue
C ---> Temperature     : n'est pas une variable transportee
C ---> Enthalpie       : a resoudre
C
            IF(ISPECF.EQ.2) THEN
C
              IPHAS = IPHSCA(ISCAL)
C
              CALL CFENER
C             ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS)   , NCETSM(IPHAS)   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCAL  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) , IA(IICESM(IPHAS)) , IA(IITPSM(IPHAS)) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)) , RA(ISMACE(IPHAS)) ,
     &   VISCF  , VISCB  ,
     &   DAM    , XAM    ,
     &   DRTP   , SMBRS  , ROVSDT ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
            ENDIF
C
          ENDIF
C
C     Pour le compressible, on ne resout pas celles deja resolues ou
C       qui ne doivent pas l'etre
C     Pour les autres physiques, on resout les scalaires dans l'ordre
C       (pas de tests a faire)
C
          IF(ISPECF.EQ.0) THEN
C
C ---> Variances et scalaires
C
C   iscavr dira scalaire/variance
C   itspdv dira si calcul des termes de prod et dissip supplementaires
C         ou non (1 : oui, 0 : non)
C
C         si iscavr = 0
C           scalaire
C           itspdv = 0
C         sinon
C           variance
C           si iscavr > 0 et iscavr < nscal+1
C             itspdv = 1
C           sinon
C             pour le moment, on s'arrete
C             a terme, les combustionnistes pourront donner leur propre
C               grandeur scalaire associee a la variance et
C               eventuellement reconstruite en dehors de covofi.
C               il faudra alors peut etre declarer des tableaux
C               de travail suppl
C           fin si
C         fin si
C
            IISC = ISCAL
            IF(ISCAVR(IISC).EQ.0) THEN
              ITSPDV = 0
            ELSEIF(ISCAVR(IISC).GT.0.AND.ISCAVR(IISC).LE.NSCAL) THEN
              ITSPDV = 1
            ELSE
              WRITE(NFECRA,9000)IISC,IISC,NSCAL,ISCAVR(IISC)
              CALL CSEXIT (1)
            ENDIF
C
C
C ---> Appel a covofi pour la resolution
C
            IPHAS  = IPHSCA(IISC)
C
            CALL COVOFI
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS) , NCETSM(IPHAS) ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IISC   , ITSPDV ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) , IA(IICESM(IPHAS)) , IA(IITPSM(IPHAS)) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DTR    , RTP    , RTPA   , PROPCE , PROPFA , PROPFB , TSLAGR ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)) , RA(ISMACE(IPHAS)) ,
     &   VISCF  , VISCB  ,
     &   DAM    , XAM    ,
     &   DRTP   , SMBRS  , ROVSDT ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C ---> Versions Electriques
C             Effet Joule
C             Arc Electrique
C             Conduction ionique
C
C     On calcule ici j, E, j.E reels et imagimaires
C
            IF ( IPPMOD(IELJOU).GE.1 .OR.
     &           IPPMOD(IELARC).GE.1 .OR.
     &           IPPMOD(IELION).GE.1       ) THEN
C
C
C     On utilise le  fait que les scalaires sont dans l'ordre
C       H, PotR, [PotI], [A] pour faire le calcul de j, E et j.E
C       apres  la determination de PotR [et PotI].
C
              ICALC = 0
C            On peut y aller apres PotR si on est en arc
C                                         ou en Joule sans PotI
C
              IF(IPPMOD(IELARC).GE.1.OR.IPPMOD(IELJOU).EQ.1
     &             .OR.IPPMOD(IELJOU).EQ.3) THEN
                IF(ISCAL.EQ.IPOTR) THEN
                  ICALC = 1
                ENDIF
              ENDIF
C     On y va apres PotI si on est en Joule avec PotI
              IF(IPPMOD(IELJOU).EQ.2 .OR. IPPMOD(IELJOU).EQ.4) THEN
                IF(ISCAL.EQ.IPOTI) THEN
                  ICALC = 1
                ENDIF
              ENDIF
C
              IF(ICALC.EQ.1) THEN
C
C     Calcul de j, E et j.E
                IAPPEL = 1
C
                CALL ELFLUX
C               ===========
     & ( IDEBIA , IDEBRA , IAPPEL ,
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
C
C
C     Recalage des variables electriques j, j.E (et Pot, E)
C
                IF ( IELCOR .EQ.1  .AND. NTCABS .GT. 1 ) THEN
C
                  CALL USELRC
C                 ===========
     &   (IDEBIA , IDEBRA ,
     &    NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &    NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &    NVAR   , NSCAL  , NPHAS  ,
     &    NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &    IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &    IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &    IDEVEL , ITUSER , IA     ,
     &    XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &    DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &    COEFA  , COEFB  , VISCF  , VISCB  ,
     &    W1     , W2     , W3     , W4     , W5     ,
     &    W6     , W7     , W8     , W9     ,
     &    RDEVEL , RTUSER , RA     )
C
                ENDIF
C
              ENDIF
C
            ENDIF
C
          ENDIF
C
C
C ---> Fin de la Boucle sur les scalaires physiques particulieres.
        ENDDO
      ENDIF
C
C
C     On calcule ici A, B, jxB
C
      IF ( IPPMOD(IELARC).GE.1       ) THEN
C
C     On utilise le  fait que les scalaires sont dans l'ordre
C       H, PotR, [PotI], [A] pour faire le calcul de A, B, jxB
C       apres la determination et le recalage de j
        IAPPEL = 2
C
        CALL ELFLUX
C       ===========
     & ( IDEBIA , IDEBRA , IAPPEL ,
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
C
      ENDIF
C
C=======================================================================
C 3. TRAITEMENT DES SCALAIRES UTILISATEURS STANDARD
C     On voit que meme s'ils sont numerotes en premier, on peut les
C       traiter en dernier si on veut.
C     On peut imaginer aussi de by-passer cette phase si le modele a
C       physique particuliere le demande.
C=======================================================================
C
      IF(NSCAUS.GT.0) THEN
C
C ---> Boucle sur les scalaires utilisateur.
C
        DO II = 1, NSCAUS
C
          ISCAL = II
          IVAR  = ISCA(ISCAL)
C
C ---> Pas de temps (avec facteur multiplicatif eventuel)
C
          IF(CDTVAR(IVAR).NE.1.D0) THEN
            DO IEL = 1, NCEL
              DTR(IEL) = DT(IEL)*CDTVAR(IVAR)
            ENDDO
          ELSE
            DO IEL = 1, NCEL
              DTR(IEL) = DT(IEL)
            ENDDO
          ENDIF
C
C
C ---> Variances et scalaires
C
C   iscavr dira scalaire/variance
C   itspdv dira si calcul des termes de prod et dissip supplementaires
C         ou non (1 : oui, 0 : non)
C
C         si iscavr = 0
C           scalaire
C           itspdv = 0
C         sinon
C           variance
C           si iscavr > 0 et iscavr < nscal+1
C             itspdv = 1
C           sinon
C             pour le moment, on s'arrete
C             a terme, les combustionnistes pourront donner leur propre
C               grandeur scalaire associee a la variance et
C               eventuellement reconstruite en dehors de covofi.
C               il faudra alors peut etre declarer des tableaux
C               de travail suppl
C           fin si
C         fin si
C
          IISC = ISCAL
          IF(ISCAVR(IISC).EQ.0) THEN
            ITSPDV = 0
          ELSEIF(ISCAVR(IISC).GT.0.AND.ISCAVR(IISC).LE.NSCAL) THEN
            ITSPDV = 1
          ELSE
            WRITE(NFECRA,9000)IISC,IISC,NSCAL,ISCAVR(IISC)
            CALL CSEXIT (1)
          ENDIF
C
C
C ---> Appel a covofi pour la resolution
C
          IPHAS  = IPHSCA(IISC)
C
          CALL COVOFI
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS) , NCETSM(IPHAS) ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IISC   , ITSPDV ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) , IA(IICESM(IPHAS)) , IA(IITPSM(IPHAS)) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DTR    , RTP    , RTPA   , PROPCE , PROPFA , PROPFB , TSLAGR ,
     &   COEFA  , COEFB  , RA(ICKUPD(IPHAS)) , RA(ISMACE(IPHAS)) ,
     &   VISCF  , VISCB  ,
     &   DAM    , XAM    ,
     &   DRTP   , SMBRS  , ROVSDT ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C
C
C
C
C ---> Fin de la Boucle sur les scalaires utilisateurs.
        ENDDO
C
      ENDIF
C
C=======================================================================
C 4.  FORMATS
C=======================================================================
C
#if defined(_CS_LANG_FR)
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA RESOLUTION DES SCALAIRES         ',/,
     &'@    *********                                               ',/,
     &'@    SCALAIRE NUMERO ',I10                                    ,/,
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
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#else
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT WHILE SOLVING SCALARS EQUATIONS          ',/,
     &'@    ********                                                ',/,
     &'@    SCALAR NUMBER ',I10                                      ,/,
     &'@    ISCAVR(',I10   ,') MUST BE A POSITIVE OR NULL INTEGER   ',/,
     &'@      AND LOWER OR EQUAL THAN NSCAL = ', I10                 ,/,
     &'@    ITS VALUE IS ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculaton will not be run.                           ',/,
     &'@                                                            ',/,
     &'@  If ISCAVR(I) is null, the scalar I is not a variance      ',/,
     &'@  If ISCAVR(I) is positive, the scalar I is a variance:     ',/,
     &'@    it is the variance of the fluctuations of the scalar J  ',/,
     &'@    whose number is ISCAVR(I)                               ',/,
     &'@                                                            ',/,
     &'@  Verify usini1.                                            ',/,
     &'@  Contact support.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#endif
C
C----
C FIN
C----
      RETURN
C
      END
c@z
