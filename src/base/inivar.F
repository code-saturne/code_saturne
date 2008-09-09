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
                        SUBROUTINE INIVAR
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCOFAB ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , FRCXT  ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC INITIALISATION DES VARIABLES DE CALCUL, DU PAS DE TEMPS
CFONC  ET DU TABLEAU INDICATEUR DU CALCUL DE LA DISTANCE A LA PAROI
CFONC  PAR L'UTILISATEUR (apres relecture eventuelle d'un fichier suite)
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
CARGU ! NCOFAB       ! E  !  -> ! NOMBRE DE COUPLES COEFA/B POUR LES CL!
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
CARGU ! DT(NCELET)   ! TR ! <-> ! VALEUR DU PAS DE TEMPS               !
CARGU ! RTP          ! TR ! <-> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA COEFB  ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! FRCXT(NCELET,! TR ! <-> ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU !   3,NPHAS)   !    !     !  HYDROSTATIQUE                       !
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
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "ihmpre.h"
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
      INTEGER          NVAR   , NSCAL  , NPHAS  , NCOFAB
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          MAXELT, ILS
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,NCOFAB), COEFB(NFABOR,NCOFAB)
      DOUBLE PRECISION FRCXT(NCELET,3,NPHAS)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      INTEGER          IDEBIA, IDEBRA, IFINIA
      INTEGER          IVAR  , IPHAS , ISCAL , IPHASS, IMOM
      INTEGER          IEL
      INTEGER          ICLIP , IPP  , IOK   , II
      INTEGER          IKIPH , IEIPH , IR11IP, IR22IP, IR33IP, IPHIPH
      INTEGER          IOMGIP, IDTCM , IPCMOM, IIPTOT, IPRIPH
      INTEGER          IBORMO(NBMOMX)
      DOUBLE PRECISION VALMAX, VALMIN, VFMIN , VFMAX
      DOUBLE PRECISION VDTMAX, VDTMIN
      DOUBLE PRECISION XEKMIN, XEPMIN, XOMGMN, XPHMIN, XPHMAX
      DOUBLE PRECISION X11MIN, X22MIN, X33MIN, VALMOM
      DOUBLE PRECISION VMOMAX(NBMOMX), VMOMIN(NBMOMX)
      DOUBLE PRECISION RO0IPH, P0IPH, PR0IPH, XXP0, XYP0, XZP0
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IOK = 0
C
C=======================================================================
C 2. ON REPASSE LA MAIN A L'UTILISATEUR POUR LA PROGRAMMATION DES
C    INITIALISATIONS QUI LUI SONT PROPRES
C=======================================================================
C
C Indicateur d'initialisation des scalaires par l'utilisateur
C (mis a 1 si passage dans USINIV ou PPINIV ou dans l'IHM ; a 0 sinon)
C
      IUSINI = 1
C
C   - Interface Code_Saturne
C     ======================
C
      IF (ISUITE.EQ.0 .AND. IIHMPR.EQ.1) THEN
C
        CALL UIINIV (NCELET, ISCA, RTP)
C       ===========
C
      ENDIF
C
C   - Sous-programme utilisateur
C     ==========================
C
      MAXELT = MAX(NCELET, NFAC, NFABOR)
      ILS    = IDEBIA
      IFINIA = ILS + MAXELT
      CALL IASIZE('INIVAR',IFINIA)
C
      CALL USINIV
C     ===========
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
C     Avec l'interface, il peut y avoir eu initialisation,
C       meme si usiniv n'est pas utilise.
        IF (ISUITE.EQ.0 .AND. IIHMPR.EQ.1) THEN
          IUSINI = 1
        ENDIF
C
C
C ON FAIT DE LA PHYSIQUE PARTICULIERE
C   On pourrait remonter la partie init non utilisateur de ppiniv avant lecamo
C     dans iniva0, mais il faudrait quand meme conserver ici l'appel a
C     ppiniv car il encapsule les appels aux ss pgm utilisateur similaires a
C     usiniv.
      IF (NSCAPP.GT.0) THEN
C
        IUSINI = 1
C
        CALL PPINIV
C       ===========
     & ( IDEBIA , IDEBRA ,
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
      ENDIF
C
C Si l'utilisateur a change Ptot, on change P* en consequence,
C sinon on met Ptot a P0 + rho.g.r
C A priori l'utilisateur remplira les NCEL valeurs ou rien du
C  tout, mais on ne sait jamais ...
C En compressible, Ptot n'est pas defini (correspond directement a RTP(.,IPR)
C
      IF  (IPPMOD(ICOMPF).LT.0) THEN
        DO IPHAS = 1, NPHAS
          IPRIPH = IPR(IPHAS)
          IIPTOT = IPPROC(IPRTOT(IPHAS))
          RO0IPH = RO0  (IPHAS)
          P0IPH  = P0   (IPHAS)
          PR0IPH = PRED0(IPHAS)
          XXP0   = XYZP0(1,IPHAS)
          XYP0   = XYZP0(2,IPHAS)
          XZP0   = XYZP0(3,IPHAS)
          DO IEL = 1, NCEL
            IF (PROPCE(IEL,IIPTOT).GT.-0.5D0*RINFIN) THEN
              RTP(IEL,IPRIPH) = PROPCE(IEL,IIPTOT)
     &             - RO0IPH*( GX*(XYZCEN(1,IEL)-XXP0)
     &                      + GY*(XYZCEN(2,IEL)-XYP0)
     &                      + GZ*(XYZCEN(3,IEL)-XZP0) )
     &             + PR0IPH - P0IPH
            ELSE
              PROPCE(IEL,IIPTOT) = RTP(IEL,IPRIPH)
     &             + RO0IPH*( GX*(XYZCEN(1,IEL)-XXP0)
     &                      + GY*(XYZCEN(2,IEL)-XYP0)
     &                      + GZ*(XYZCEN(3,IEL)-XZP0) )
     &             + P0IPH - PR0IPH
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
C
C
C=======================================================================
C 3.  CLIPPING DES GRANDEURS TURBULENTES (UTILISATEUR OU SUITE)
C     (pour ITYTUR=2, 3, 5 ou 6)
C     Si l'utilisateur est intervenu dans USINIV, PPINIV ou via l'interface
C         et a impose des valeurs "correctes" (au sens k, eps, Rii > 0)
C         on considere qu'il s'agit d'une initialisation admissible,
C         on la clippe pour la rendre coherente avec le clipping du code
C         et on continue le calcul
C     Si l'utilisateur est intervenu dans USINIV, PPINIV ou via l'interface
C         et a impose des valeurs visiblement erronees
C         (k, eps ou Rii < 0), on s'arrete (il s'est sans doute trompe).
C     On adopte le meme traitement en suite de calcul
C       pour assurer un comportement identique en suite entre un calcul
C       ou l'utilisateur modifie une variable avec usiniv (mais pas la
C       turbulence) et un calcul ou l'utilisateur ne modifie pas usiniv.
C     S'il n'y a ni suite ni intervention dans USINIV ou PPINIV ou via l'interface,
C       les grandeurs ont deja ete clippees par defaut, sauf si UREF n'a pas
C       (ou a mal) ete initialise. Dans ce cas on avertit aussi l'utilisateur et on
C       stoppe le calcul.
C
C     Pour resumer :
C      -en   suite  avec des valeurs positives pour k, eps, Rii : on clippe
C      -avec usiniv ou ppiniv ou interface
C                   avec des valeurs positives pour k, eps, Rii : on clippe
C      -non suite sans usiniv ni ppiniv ni interface avec UREF positif :
C                                      grandeurs par defaut (deja clippees)
C      -non suite sans usiniv ni ppiniv ni interface avec UREF negatif : stop
C      -suite ou usiniv ou ppiniv ou interface
C                   avec une valeur negative de k, eps ou Rii : stop
C                   avec une valeur hors de [0;2] pour phi : stop
C         (on souhaite indiquer a l'utilisateur que son fichier suite est
C          bizarre ou que son initialisation est fausse et qu'il a donc
C          fait au moins une erreur qui peut en cacher d'autres)
C=======================================================================
C
      IF(IUSINI.EQ.1.OR.ISUITE.EQ.1) THEN
C
        DO IPHAS = 1, NPHAS
          IF(ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50) THEN
C
            IKIPH  = IK (IPHAS)
            IEIPH  = IEP(IPHAS)
C
            XEKMIN = RTP(1,IKIPH)
            XEPMIN = RTP(1,IEIPH)
            DO IEL = 1, NCEL
              XEKMIN = MIN(XEKMIN,RTP(IEL,IKIPH) )
              XEPMIN = MIN(XEPMIN,RTP(IEL,IEIPH))
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARMIN (XEKMIN)
C             ===========
              CALL PARMIN (XEPMIN)
C             ===========
            ENDIF
C
            IF(XEKMIN.GE.0.D0.AND.XEPMIN.GE.0.D0) THEN
              ICLIP = 1
              IPHASS = IPHAS
              CALL CLIPKE( NCELET , NCEL   , NVAR   , NPHAS  ,
C             ===========
     &                     IPHASS , ICLIP  , IWARNI(IKIPH) ,
     &                     PROPCE , RTP    )
            ELSE
              WRITE(NFECRA,3020) IPHAS,XEKMIN,XEPMIN
              IOK = IOK + 1
            ENDIF
C
C     En v2-f, on verifie aussi que phi est compris entre 0 et 2
            IF (ITURB(IPHAS).EQ.50) THEN
              IPHIPH = IPHI(IPHAS)
C
              XPHMIN = RTP(1,IPHIPH)
              XPHMAX = RTP(1,IPHIPH)
              DO IEL = 1, NCEL
                XPHMIN = MIN(XPHMIN,RTP(IEL,IPHIPH) )
                XPHMAX = MAX(XPHMAX,RTP(IEL,IPHIPH) )
              ENDDO
              IF (IRANGP.GE.0) THEN
                CALL PARMIN (XPHMIN)
C               ===========
                CALL PARMAX (XPHMAX)
C               ===========
              ENDIF
C
C     Par coherence avec clpv2f, on ne clippe qu'a zero et pas a 2
C              IF(XPHMIN.LT.0.D0 .OR. XPHMAX.GT.2.D0) THEN
              IF(XPHMIN.LT.0.D0) THEN
                WRITE(NFECRA,3021) IPHAS,XPHMIN,XPHMAX
                IOK = IOK + 1
              ENDIF
C
            ENDIF
C
          ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
C
            IR11IP = IR11(IPHAS)
            IR22IP = IR22(IPHAS)
            IR33IP = IR33(IPHAS)
            IEIPH  = IEP (IPHAS)
C
            X11MIN = RTP(1,IR11IP)
            X22MIN = RTP(1,IR22IP)
            X33MIN = RTP(1,IR33IP)
            XEPMIN = RTP(1,IEIPH)
            DO IEL = 1, NCEL
              X11MIN = MIN(X11MIN,RTP(IEL,IR11IP))
              X22MIN = MIN(X22MIN,RTP(IEL,IR22IP))
              X33MIN = MIN(X33MIN,RTP(IEL,IR33IP))
              XEPMIN = MIN(XEPMIN,RTP(IEL,IEIPH) )
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARMIN (X11MIN)
C             ===========
              CALL PARMIN (X22MIN)
C             ===========
              CALL PARMIN (X33MIN)
C             ===========
              CALL PARMIN (XEPMIN)
C             ===========
            ENDIF
            IF (X11MIN.GE.0.D0.AND.X22MIN.GE.0.D0.AND.
     &          X33MIN.GE.0.D0.AND.XEPMIN.GE.0.D0 ) THEN
              ICLIP = 1
              IPHASS = IPHAS
              CALL CLPRIJ( NCELET , NCEL   , NVAR   , NPHAS  ,
C             ===========
     &                     IPHASS , ICLIP  ,
     &                     PROPCE , RTP    , RTP    )
            ELSE
              WRITE(NFECRA,3030) IPHAS,X11MIN,X22MIN,X33MIN,XEPMIN
              IOK = IOK + 1
            ENDIF
C
          ELSEIF(ITURB(IPHAS).EQ.60) THEN
C
            IKIPH   = IK  (IPHAS)
            IOMGIP  = IOMG(IPHAS)
C
            XEKMIN = RTP(1,IKIPH )
            XOMGMN = RTP(1,IOMGIP)
            DO IEL = 1, NCEL
              XEKMIN = MIN(XEKMIN,RTP(IEL,IKIPH ))
              XOMGMN = MIN(XOMGMN,RTP(IEL,IOMGIP))
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARMIN (XEKMIN)
C             ===========
              CALL PARMIN (XOMGMN)
C             ===========
            ENDIF
C
C     En k-omega on clippe seulement a 0
            IF(XEKMIN.LT.0.D0 .OR. XOMGMN.LT.0.D0) THEN
              WRITE(NFECRA,3031) IPHAS,XEKMIN,XOMGMN
              IOK = IOK + 1
            ENDIF
C
          ENDIF
        ENDDO
C
      ELSE
C
        DO IPHAS = 1, NPHAS
          IF (ITURB(IPHAS).NE.0 .AND. ITURB(IPHAS).NE.10
     &         .AND. ITYTUR(IPHAS).NE.4) THEN
            IF (UREF(IPHAS).LT.0.D0) THEN
              WRITE(NFECRA,3032) IPHAS,UREF(IPHAS)
              IOK = IOK + 1
            ENDIF
          ENDIF
        ENDDO
C
      ENDIF
C
C=======================================================================
C 4.  CLIPPING DES SCALAIRES (UTILISATEUR OU SUITE)
C     Si l'utilisateur est intervenu dans USINIV ou PPINIV et
C       a impose des valeurs "correctes" (au sens comprises dans des bornes
C         simplifiees a base de 0, scamin, scamax)
C         on considere qu'il s'agit d'une initialisation admissible,
C         on la clippe pour la rendre coherente avec le clipping du code
C         et on continue le calcul
C       si l'utilisateur a impose des valeurs visiblement erronees
C         (au sens comprises dans des bornes simplifiees a base de 0, scamin,
C          scamax), on s'arrete (il s'est sans doute trompe).
C     On adopte le meme traitement en suite de calcul
C       pour assurer un comportement identique en suite entre un calcul
C       ou l'utilisateur modifie une variable avec usiniv (mais pas un
C       scalaire) et un calcul ou l'utilisateur ne modifie pas usiniv.
C     Sinon, les grandeurs ont deja ete clippees apres les init par defaut
C
C     Pour resumer :
C      -en   suite  avec des valeurs grossierement admissibles : on clippe
C      -avec usiniv ou ppiniv
C                   avec des valeurs grossierement admissibles : on clippe
C      -non suite sans usiniv ni ppiniv :
C                                      grandeurs par defaut (deja clippees)
C      -suite ou usiniv ou ppiniv
C                   avec une valeur grossierement non admissible : stop
C         (on souhaite indiquer a l'utilisateur que son fichier suite est
C          bizarre ou que son initialisation est fausse et qu'il a donc
C          fait au moins une erreur qui peut en cacher d'autres)
C=======================================================================
C
C On traite tous les scalaires d'abord, car ils peuvent etre necessaires
C     pour clipper les variances
C
      IF(NSCAL.GT.0.AND.(IUSINI.EQ.1.OR.ISUITE.EQ.1)) THEN
C
C     Scalaires non variance
C
        DO II = 1, NSCAL
          IF(ISCAVR(II).LE.0.OR.ISCAVR(II).GT.NSCAL) THEN
C
            IF(SCAMIN(II).LE.SCAMAX(II)) THEN
              IVAR = ISCA(II)
              VALMAX = RTP(1  ,IVAR)
              VALMIN = RTP(1  ,IVAR)
              DO IEL = 1, NCEL
                VALMAX = MAX(VALMAX,RTP(IEL,IVAR))
                VALMIN = MIN(VALMIN,RTP(IEL,IVAR))
              ENDDO
              IF (IRANGP.GE.0) THEN
                CALL PARMAX (VALMAX)
C               ===========
                CALL PARMIN (VALMIN)
C               ===========
              ENDIF
C
C     Verification de la coherence pour les clippings
C                                           des scalaires non variance.
              IF(VALMIN.GE.SCAMIN(II).AND.VALMAX.LE.SCAMAX(II)) THEN
                ISCAL = II
                CALL CLPSCA
C               ===========
     &                ( NCELET , NCEL   , NVAR   , NSCAL  , ISCAL  ,
     &                  PROPCE , RA     , RTP    )
              ELSE
                CHAINE = NOMVAR(IPPRTP(ISCA(II)))
                WRITE(NFECRA,3040) II,CHAINE(1:8),
     &                             VALMIN,SCAMIN(II),VALMAX,SCAMAX(II)
                IOK = IOK + 1
              ENDIF
            ENDIF
C
          ENDIF
        ENDDO
C
C
C     Variances
C
        DO II = 1, NSCAL
          IF(ISCAVR(II).GT.0.AND.ISCAVR(II).LE.NSCAL) THEN
C
            IF(SCAMIN(II).LE.SCAMAX(II)) THEN
              IVAR = ISCA(II)
              VALMAX = RTP(1  ,IVAR)
              VALMIN = RTP(1  ,IVAR)
              DO IEL = 1, NCEL
                VALMAX = MAX(VALMAX,RTP(IEL,IVAR))
                VALMIN = MIN(VALMIN,RTP(IEL,IVAR))
              ENDDO
              IF (IRANGP.GE.0) THEN
                CALL PARMAX (VALMAX)
C               ===========
                CALL PARMIN (VALMIN)
C               ===========
              ENDIF
C
C     Verification de la coherence pour les clippings de variance.
C     Pour iclvfl = 1 on ne verifie que > 0 sinon ca va devenir difficile
C     de faire une initialisation correcte.
C
              IF(ICLVFL(II).EQ.0) THEN
C       On pourrait clipper dans le cas ou VALMIN.GE.0, mais ca
C       n'apporterait rien, par definition
                IF(VALMIN.LT.0.D0) THEN
                  CHAINE = NOMVAR(IPPRTP(ISCA(II)))
                  WRITE(NFECRA,3050)II,CHAINE(1:8),
     &                              VALMIN,SCAMIN(II),VALMAX,SCAMAX(II)
                  IOK = IOK + 1
                ENDIF
              ELSEIF(ICLVFL(II).EQ.1) THEN
C Ici on clippe pour etre coherent avec la valeur du scalaire
                IF(VALMIN.GE.0.D0) THEN
                  ISCAL = II
                  CALL CLPSCA
C                 ===========
     &            ( NCELET , NCEL   , NVAR   , NSCAL  , ISCAL  ,
     &              PROPCE , RTP(1,ISCA(ISCAVR(II))) , RTP      )
                ELSE
                  CHAINE = NOMVAR(IPPRTP(ISCA(II)))
                  WRITE(NFECRA,3050)II,CHAINE(1:8),
     &                              VALMIN,SCAMIN(II),VALMAX,SCAMAX(II)
                  IOK = IOK + 1
                ENDIF
              ELSEIF(ICLVFL(II).EQ.2) THEN
                VFMIN = 0.D0
                VFMIN = MAX(SCAMIN(ISCAL),VFMIN)
                VFMAX = SCAMAX(ISCAL)
C On pourrait clipper dans le cas ou VALMIN.GE.VFMIN.AND.VALMAX.LE.VFMAX
C     mais ca n'apporterait rien, par definition
                IF(VALMIN.LT.VFMIN.OR.VALMAX.GT.VFMAX) THEN
                  CHAINE = NOMVAR(IPPRTP(ISCA(II)))
                  WRITE(NFECRA,3051)II,CHAINE(1:8),
     &                              VALMIN,SCAMIN(II),VALMAX,SCAMAX(II),
     &                              II,ICLVFL(II)
                  IOK = IOK + 1
                ENDIF
              ENDIF
            ENDIF
C
          ENDIF
        ENDDO
C
      ENDIF
C
C
C=======================================================================
C 6.  IMPRESSIONS DE CONTROLE POUR LES INCONNUES, LE PAS DE TEMPS
C        LE CUMUL DES DUREE POUR LES MOYENNES
C=======================================================================
C
      WRITE(NFECRA,2000)
C
C     Inconnues de calcul : on affiche les bornes
      DO IPP  = 2, NVPPMX
        IF(ITRSVR(IPP ).GE.1) THEN
          IVAR = ITRSVR(IPP )
          VALMAX = -GRAND
          VALMIN =  GRAND
          DO IEL = 1, NCEL
            VALMAX = MAX(VALMAX,RTP(IEL,IVAR))
            VALMIN = MIN(VALMIN,RTP(IEL,IVAR))
          ENDDO
          IF (IRANGP.GE.0) THEN
            CALL PARMAX (VALMAX)
C           ===========
            CALL PARMIN (VALMIN)
C           ===========
          ENDIF
          CHAINE = NOMVAR(IPP )
          WRITE(NFECRA,2010)CHAINE(1:8),VALMIN,VALMAX
        ENDIF
      ENDDO
      WRITE(NFECRA,2020)
C
C     Moyennes  : on affiche les bornes
      IF(NBMOMT.GT.0) THEN
        DO IMOM = 1, NBMOMT
          IPCMOM = IPPROC(ICMOME(IMOM))
C
C       Si on ne (re)initialise pas
          IF(IMOOLD(IMOM).NE.-1) THEN
            VALMAX = -GRAND
            VALMIN =  GRAND
C         Si le cumul en temps est variable en espace
            IF(IDTMOM(IMOM).GT.0) THEN
              IDTCM  = IPPROC(ICDTMO(IDTMOM(IMOM)))
              DO IEL = 1, NCEL
                VALMOM = PROPCE(IEL,IPCMOM)/
     &                     MAX(PROPCE(IEL,IDTCM),EPZERO)
                VALMAX = MAX(VALMAX,VALMOM)
                VALMIN = MIN(VALMIN,VALMOM)
              ENDDO
C         Si le cumul en temps est uniforme
            ELSE
              IDTCM  =-IDTMOM(IMOM)
              DO IEL = 1, NCEL
                VALMOM = PROPCE(IEL,IPCMOM)/
     &                     MAX(DTCMOM(IDTCM),EPZERO)
                VALMAX = MAX(VALMAX,VALMOM)
                VALMIN = MIN(VALMIN,VALMOM)
              ENDDO
            ENDIF
            IF (IRANGP.GE.0) THEN
              CALL PARMAX (VALMAX)
C             ===========
              CALL PARMIN (VALMIN)
C             ===========
            ENDIF
C       Si on  (re)initialise
          ELSE
            VALMAX = 0.D0
            VALMIN = 0.D0
          ENDIF
C
          CHAINE = NOMVAR(IPPPRO(IPCMOM))
          WRITE(NFECRA,2010)CHAINE(1:8),VALMIN,VALMAX
C
        ENDDO
        WRITE(NFECRA,2020)
      ENDIF
C
      IF (IDTVAR.GE.0) THEN
C     Pas de temps : on affiche les bornes
C                    si < 0 on s'arrete
        VDTMAX = -GRAND
        VDTMIN =  GRAND
        DO IEL = 1, NCEL
          VDTMAX = MAX(VDTMAX,DT    (IEL))
          VDTMIN = MIN(VDTMIN,DT    (IEL))
        ENDDO
        IF (IRANGP.GE.0) THEN
          CALL PARMAX (VDTMAX)
C         ===========
          CALL PARMIN (VDTMIN)
C         ===========
        ENDIF
        WRITE(NFECRA,2010)'PasDeTmp',VDTMIN,VDTMAX
        WRITE(NFECRA,2020)
C
        IF (VDTMIN.LE.ZERO) THEN
          WRITE(NFECRA,3010) VDTMIN
          IOK = IOK + 1
        ENDIF
C
      ENDIF
C
C     Cumul du temps associe aux moments : on affiche les bornes
C                                          si < 0 on s'arrete
C
      IF(NBMOMT.GT.0) THEN
C
C     Indicateur de calcul des bornes pour les cumuls non uniformes
        DO IMOM = 1, NBMOMT
          IF(IDTMOM(IMOM).GT.0) THEN
            IBORMO(ICDTMO(IDTMOM(IMOM))) =  0
            VMOMAX(ICDTMO(IDTMOM(IMOM))) = -GRAND
            VMOMIN(ICDTMO(IDTMOM(IMOM))) =  GRAND
          ENDIF
        ENDDO
C
C     Calcul des bornes des cumuls non uniformes
        DO IMOM = 1, NBMOMT
          IF(IDTMOM(IMOM).GT.0) THEN
            IF(IBORMO(ICDTMO(IDTMOM(IMOM))).EQ.0) THEN
              IDTCM  = IPPROC(ICDTMO(IDTMOM(IMOM)))
              VDTMAX = -GRAND
              VDTMIN =  GRAND
              DO IEL = 1, NCEL
                VDTMAX = MAX(VDTMAX,PROPCE(IEL,IDTCM))
                VDTMIN = MIN(VDTMIN,PROPCE(IEL,IDTCM))
              ENDDO
              IF (IRANGP.GE.0) THEN
                CALL PARMAX (VDTMAX)
C               ===========
                CALL PARMIN (VDTMIN)
C               ===========
              ENDIF
              VMOMAX(ICDTMO(IDTMOM(IMOM))) = VDTMAX
              VMOMIN(ICDTMO(IDTMOM(IMOM))) = VDTMIN
              IBORMO(ICDTMO(IDTMOM(IMOM))) = 1
            ENDIF
          ENDIF
        ENDDO
C
C     Impression des bornes
        WRITE(NFECRA,2030)
        DO IMOM = 1, NBMOMT
          IF(IDTMOM(IMOM).GT.0) THEN
            WRITE(NFECRA,2040) IMOM,VMOMIN(ICDTMO(IDTMOM(IMOM))),
     &                              VMOMAX(ICDTMO(IDTMOM(IMOM))),
     &                              'Variable'
          ELSEIF(IDTMOM(IMOM).LT.0) THEN
            WRITE(NFECRA,2040) IMOM,DTCMOM(-IDTMOM(IMOM))       ,
     &                              DTCMOM(-IDTMOM(IMOM))       ,
     &                              'Uniforme'
          ENDIF
        ENDDO
        WRITE(NFECRA,2050)
C
C     On s'arrete si des cumuls sont negatifs
        DO IMOM = 1, NBMOMT
          IF(IDTMOM(IMOM).GT.0) THEN
            IF(VMOMIN(ICDTMO(IDTMOM(IMOM))).LT.ZERO) THEN
              WRITE(NFECRA,3011) IMOM,VMOMIN(ICDTMO(IDTMOM(IMOM)))
              IOK = IOK + 1
            ENDIF
          ELSEIF(IDTMOM(IMOM).LT.0) THEN
            IF(DTCMOM(-IDTMOM(IMOM)).LT.ZERO) THEN
              WRITE(NFECRA,3011) IMOM,DTCMOM(-IDTMOM(IMOM))
              IOK = IOK + 1
            ENDIF
          ENDIF
        ENDDO
C
      ENDIF
C
C=======================================================================
C 7.  ARRET GENERAL SI PB
C=======================================================================
C
      IF (IOK.GT.0) THEN
        WRITE(NFECRA,3090) IOK
        CALL CSEXIT (1)
      ENDIF
C
      WRITE(NFECRA,3000)
C
C----
C FORMATS
C----
C
C
#if defined(_CS_LANG_FR)
C
 2000 FORMAT(
     &'                                                             ',/,
     &' ----------------------------------------------------------- ',/,
     &'                                                             ',/,
     &'                                                             ',/,
     &' ** INITIALISATION DES VARIABLES                             ',/,
     &'    ----------------------------                             ',/,
     &'                                                             ',/,
     &' ---------------------------------                           ',/,
     &'  Variable  Valeur min  Valeur max                           ',/,
     &' ---------------------------------                           '  )
 2010 FORMAT(
     & 2X,     A8,      E12.4,      E12.4                              )
 2020 FORMAT(
     &' ---------------------------------                           ',/)
 2030 FORMAT(
     &' Duree cumulee :                                             ',/,
     &' ------------------------------------------------------------',/,
     &'   Moyenne  Valeur min  Valeur max Uniforme/Variable en espac',/,
     &' ------------------------------------------------------------'  )
 2040 FORMAT(
     &        I10,      E12.4,      E12.4,1X,   A8                     )
 2050 FORMAT(
     &' ------------------------------------------------------------',/)
C
 3000 FORMAT(/,/,
     &'-------------------------------------------------------------',/)
 3010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES     ',/,
     &'@    *********                                               ',/,
     &'@    PAS DE TEMPS NEGATIF OU NUL                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  La valeur minimale du pas de temps DT est ',E14.5          ,/,
     &'@  Verifier l''initialisation dans usiniv ou le fichier suite',/,
     &'@    dans le cas ou les valeurs lues dans le fichier suite   ',/,
     &'@    sont incorrectes, on peut les modifier par usiniv       ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES     ',/,
     &'@    *********                                               ',/,
     &'@    CUMUL DE DUREE POUR LES MOYENNES NEGATIVE               ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  La valeur minimale de la duree cumulee pour la moyenne    ',/,
     &'@    IMOM = ',I10   ,' est ',E14.5                            ,/,
     &'@                                                            ',/,
     &'@  Verifier l''initialisation dans usiniv ou le fichier suite',/,
     &'@    dans le cas ou les valeurs lues dans le fichier suite   ',/,
     &'@    sont incorrectes, on peut reinitialiser la moyenne et le',/,
     &'@    cumul temporel associe en imposant IMOOLD(IMOM) = -1    ',/,
     &'@    dans usini1.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES     ',/,
     &'@    *********                                               ',/,
     &'@     TURBULENCE NEGATIVE OU NULLE                           ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Phase                       = ',I4                         ,/,
     &'@   Valeur minimale de k       = ',E14.5                      ,/,
     &'@   Valeur minimale de epsilon = ',E14.5                      ,/,
     &'@                                                            ',/,
     &'@  Verifier l''initialisation (usiniv et/ou interface),      ',/,
     &'@    le fichier suite ou bien la valeur de UREF (usini1      ',/,
     &'@    et/ou interface).                                       ',/,
     &'@  Dans le cas ou les valeurs lues dans le fichier suite     ',/,
     &'@    sont incorrectes, on peut les modifier par usiniv ou    ',/,
     &'@    par l''interface).                                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3021 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES     ',/,
     &'@    *********                                               ',/,
     &'@     VARIABLE PHI DU V2F HORS DES BORNES [0;2]              ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Phase                       = ',I4                         ,/,
     &'@   Valeur minimale de phi     = ',E14.5                      ,/,
     &'@   Valeur maximale de phi     = ',E14.5                      ,/,
     &'@                                                            ',/,
     &'@  Verifier l''initialisation (usiniv et/ou interface),      ',/,
     &'@    ou le fichier suite.                                    ',/,
     &'@  Dans le cas ou les valeurs lues dans le fichier suite     ',/,
     &'@    sont incorrectes, on peut les modifier par usiniv ou    ',/,
     &'@    par l''interface.                                       ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3030 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES     ',/,
     &'@    *********                                               ',/,
     &'@     TURBULENCE NEGATIVE OU NULLE                           ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Phase                       = ',I4                         ,/,
     &'@   Valeur minimale de R11     = ',E14.5                      ,/,
     &'@   Valeur minimale de R22     = ',E14.5                      ,/,
     &'@   Valeur minimale de R33     = ',E14.5                      ,/,
     &'@   Valeur minimale de epsilon = ',E14.5                      ,/,
     &'@                                                            ',/,
     &'@  Verifier l''initialisation (usiniv et/ou interface),      ',/,
     &'@    le fichier suite ou bien la valeur de UREF (usini1      ',/,
     &'@    et/ou interface).                                       ',/,
     &'@  Dans le cas ou les valeurs lues dans le fichier suite     ',/,
     &'@    sont incorrectes, on peut les modifier par usiniv ou    ',/,
     &'@    par l''interface).                                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3031 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES     ',/,
     &'@    *********                                               ',/,
     &'@    TURBULENCE NEGATIVE OU NULLE                            ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Phase                       = ',I4                         ,/,
     &'@   Valeur minimale de k       = ',E14.5                      ,/,
     &'@   Valeur minimale de omega   = ',E14.5                      ,/,
     &'@                                                            ',/,
     &'@  Verifier l''initialisation (usiniv et/ou interface),      ',/,
     &'@    le fichier suite ou bien la valeur de UREF (usini1      ',/,
     &'@    et/ou interface).                                       ',/,
     &'@  Dans le cas ou les valeurs lues dans le fichier suite     ',/,
     &'@    sont incorrectes, on peut les modifier par usiniv ou    ',/,
     &'@    par l''interface.                                       ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3032 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES     ',/,
     &'@    *********                                               ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@    LA VITESSE DE REFERENCE UREF N''A PAS ETE INITIALISEE   ',/,
     &'@    OU A ETE MAL INITIALISEE (VALEUR NEGATIVE).             ',/,
     &'@    ELLE VAUT ICI ',E14.5                                    ,/,
     &'@                                                            ',/,
     &'@  La turbulence n''a pas pu etre initialisee                ',/,
     &'@  Corriger la valeur de UREF (usini1 ou interface)ou bien   ',/,
     &'@    initialiser directement la turbulence dans la routine   ',/,
     &'@    (usiniv ou interface).                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3040 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES     ',/,
     &'@    *********                                               ',/,
     &'@     GRANDEUR SCALAIRE HORS BORNES                          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Scalaire numero ',I10,' : ',A8                             ,/,
     &'@  Valeur minimale             = ',E14.5                      ,/,
     &'@    Clipping demande a SCAMIN = ',E14.5                      ,/,
     &'@  Valeur maximale             = ',E14.5                      ,/,
     &'@    Clipping demande a SCAMAX = ',E14.5                      ,/,
     &'@  Les valeurs extremes ne sont pas coherentes avec les      ',/,
     &'@    limites SCAMIN et SCAMAX imposees dans usini1.          ',/,
     &'@                                                            ',/,
     &'@  Verifier l''initialisation dans usiniv ou le fichier suite',/,
     &'@    dans le cas ou les valeurs lues dans le fichier suite   ',/,
     &'@    sont incorrectes, on peut les modifier par usiniv       ',/,
     &'@  Verifier les valeurs de clipping dans usini1.             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3050 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES     ',/,
     &'@    *********                                               ',/,
     &'@     VARIANCE NEGATIVE                                      ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Scalaire numero ',I10,' : ',A8                             ,/,
     &'@  Valeur minimale             = ',E14.5                      ,/,
     &'@  Le scalaire indique ci-dessus est une variance (ISCAVR est',/,
     &'@    postif dans usini1) mais l initialisation imposee       ',/,
     &'@    dans usiniv comporte des valeurs negatives.             ',/,
     &'@                                                            ',/,
     &'@  Verifier l''initialisation dans usiniv ou le fichier suite',/,
     &'@    dans le cas ou les valeurs lues dans le fichier suite   ',/,
     &'@    sont incorrectes, on peut les modifier par usiniv       ',/,
     &'@  Verifier la definition des variances dans usini1.         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3051 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES     ',/,
     &'@    *********                                               ',/,
     &'@     VARIANCE HORS BORNES                                   ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Scalaire numero ',I10,' : ',A8                             ,/,
     &'@  Valeur minimale             = ',E14.5                      ,/,
     &'@    Clipping demande a SCAMIN = ',E14.5                      ,/,
     &'@  Valeur maximale             = ',E14.5                      ,/,
     &'@    Clipping demande a SCAMAX = ',E14.5                      ,/,
     &'@  Le scalaire indique ci-dessus est une variance (ISCAVR est',/,
     &'@    postif dans usini1) mais l initialisation imposee       ',/,
     &'@    dans usiniv comporte des valeurs situees hors des bornes',/,
     &'@    SCAMIN, SCAMAX ou inferieures a 0 et le mode de clipping',/,
     &'@    demande est ICLVFL(',I10,') = ',I10                      ,/,
     &'@                                                            ',/,
     &'@  Verifier l''initialisation dans usiniv ou le fichier suite',/,
     &'@    dans le cas ou les valeurs lues dans le fichier suite   ',/,
     &'@    sont incorrectes, on peut les modifier par usiniv       ',/,
     &'@  Verifier la definition des variances et le mode de        ',/,
     &'@    clipping demande dans usini1.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3090 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''INITIALISATION DES VARIABLES     ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    L INITIALISATION DES VARIABLES EST INCOMPLETE OU        ',/,
     &'@      INCOHERENTE AVEC LES VALEURS DES PARAMETRES DE CALCUL ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,
     &'@                                                            ',/,
     &'@  Se reporter aux impressions precedentes pour plus de      ',/,
     &'@    renseignements.                                         ',/,
     &'@  Attention a l''initialisation du pas de temps             ',/,
     &'@                                de la turbulence            ',/,
     &'@                                des scalaires et variances  ',/,
     &'@                                des moyennes temporelles    ',/,
     &'@                                                            ',/,
     &'@  Verifier usiniv ou le fichier suite.                      ',/,
     &'@    dans le cas ou les valeurs lues dans le fichier suite   ',/,
     &'@    sont incorrectes, on peut les modifier par usiniv       ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#else
C
 2000 FORMAT(
     &'                                                             ',/,
     &' ----------------------------------------------------------- ',/,
     &'                                                             ',/,
     &'                                                             ',/,
     &' ** VARIABLES INITIALIZATION                                 ',/,
     &'    ------------------------                                 ',/,
     &'                                                             ',/,
     &' ---------------------------------                           ',/,
     &'  Variable  Min. value  Max. value                           ',/,
     &' ---------------------------------                           '  )
 2010 FORMAT(
     & 2X,     A8,      E12.4,      E12.4                              )
 2020 FORMAT(
     &' ---------------------------------                           ',/)
 2030 FORMAT(
     &' Time averages (sum over the time-steps)                     ',/,
     &' ------------------------------------------------------------',/,
     &'   Average  Min. value  Max. value Uniform/Variable in space ',/,
     &' ------------------------------------------------------------'  )
 2040 FORMAT(
     &        I10,      E12.4,      E12.4,1X,   A8                     )
 2050 FORMAT(
     &' ------------------------------------------------------------',/)
C
 3000 FORMAT(/,/,
     &'-------------------------------------------------------------',/)
 3010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION          ',/,
     &'@    ********                                                ',/,
     &'@    NEGATIVE OR NULL TIME STEP                              ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  The minimum value of the time-step dt is ',E14.5           ,/,
     &'@  Verify the initialization in usiniv or the restart file   ',/,
     &'@    In the case where the values read in the restart file   ',/,
     &'@    are incorrect, they may be modified with usiniv         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION          ',/,
     &'@    ********                                                ',/,
     &'@    NEGATIVE CUMULATIVE TIME FOR THE MOMENTS                ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  The minimum value of the cumulative time for the moment   ',/,
     &'@    IMOM = ',I10   ,' est ',E14.5                            ,/,
     &'@                                                            ',/,
     &'@  Verify the initialization in usiniv or the restart file   ',/,
     &'@    In the case where the values read in the restart file   ',/,
     &'@    are incorrect, the moment and the associated cumulative ',/,
     &'@    time may be re-initialized by setting IMOOLD(IMOM) = -1 ',/,
     &'@    in usini1.                                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION          ',/,
     &'@    ********                                                ',/,
     &'@     NEGATIVE OR NULL TURBULENCE                            ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Phase                     = ',I4                           ,/,
     &'@   Minimum value of k       = ',E14.5                        ,/,
     &'@   Minimum value of epsilon = ',E14.5                        ,/,
     &'@                                                            ',/,
     &'@  Verify the initialization (usiniv and/or interface),      ',/,
     &'@    the restart file or the value of UREF (usini1 and/or    ',/,
     &'@    interface).                                             ',/,
     &'@  In the case where the values read in the restart file     ',/,
     &'@    are incorrect, they may be modified with usiniv or      ',/,
     &'@    with the interface.                                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3021 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION          ',/,
     &'@    ********                                                ',/,
     &'@     PHI VARIABLE OF V2F OUT OF BOUNDS [0;2]                ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Phase                 = ',I4                               ,/,
     &'@   Minimum value of phi = ',E14.5                            ,/,
     &'@   Maximum value of phi = ',E14.5                            ,/,
     &'@                                                            ',/,
     &'@  Verify the initialization (usiniv and/or interface),      ',/,
     &'@    the restart file.                                       ',/,
     &'@  In the case where the values read in the restart file     ',/,
     &'@    are incorrect, they may be modified with usiniv or      ',/,
     &'@    with the interface.                                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3030 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION          ',/,
     &'@    ********                                                ',/,
     &'@     NEGATIVE OR NULL TURBULENCE                            ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Phase                     = ',I4                           ,/,
     &'@   Minimum value of R11     = ',E14.5                        ,/,
     &'@   Minimum value of R22     = ',E14.5                        ,/,
     &'@   Minimum value of R33     = ',E14.5                        ,/,
     &'@   Minimum value of epsilon = ',E14.5                        ,/,
     &'@                                                            ',/,
     &'@  Verify the initialization (usiniv and/or interface),      ',/,
     &'@    the restart file or the value of UREF (usini1 and/or    ',/,
     &'@    interface).                                             ',/,
     &'@  In the case where the values read in the restart file     ',/,
     &'@    are incorrect, they may be modified with usiniv or      ',/,
     &'@    with the interface.                                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3031 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION          ',/,
     &'@    ********                                                ',/,
     &'@     NEGATIVE OR NULL TURBULENCE                            ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Phase                     = ',I4                           ,/,
     &'@   Minimum value of k       = ',E14.5                        ,/,
     &'@   Minimum value of omega   = ',E14.5                        ,/,
     &'@                                                            ',/,
     &'@  Verify the initialization (usiniv and/or interface),      ',/,
     &'@    the restart file or the value of UREF (usini1 and/or    ',/,
     &'@    interface).                                             ',/,
     &'@  In the case where the values read in the restart file     ',/,
     &'@    are incorrect, they may be modified with usiniv or      ',/,
     &'@    with the interface.                                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3032 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION          ',/,
     &'@    ********                                                ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@    THE REFERENCE VELOCITY UREF HAS NOT BEEN INITIALIZED    ',/,
     &'@    OR HAS NOT BEEN CORRECTLY INITIALIZED (NEGATIVE VALUE)  ',/,
     &'@    ITS VALUE IS ',E14.5                                     ,/,
     &'@                                                            ',/,
     &'@  The turbulence cannot be initialized                      ',/,
     &'@  Correct the value of UREF (usini1 or interface) or        ',/,
     &'@    initialize directly the turbulence with usiniv or       ',/,
     &'@    with the interface.                                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3040 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION          ',/,
     &'@    ********                                                ',/,
     &'@     SCALAR QUANTITIES OUT OF BOUNDS                        ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Scalar number ',I10,': ',A8                                ,/,
     &'@  Minimum value                = ',E14.5                     ,/,
     &'@    Desired clipping at SCAMIN = ',E14.5                     ,/,
     &'@  Maximum value                = ',E14.5                     ,/,
     &'@    Desired clipping at SCAMAX = ',E14.5                     ,/,
     &'@  The bounds are not coherent with the limits SCAMIN and    ',/,
     &'@    SCAMAX set in usini1.                                   ',/,
     &'@                                                            ',/,
     &'@  Verify the initialization in usiniv or the restart file   ',/,
     &'@    In the case where the values read in the restart file   ',/,
     &'@    are incorrect, they may be modified with usiniv         ',/,
     &'@  Verify the clipping values in usini1.                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3050 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION          ',/,
     &'@    ********                                                ',/,
     &'@     NEGATIVE VARIANCE                                      ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Scalar number ',I10,': ',A8                                ,/,
     &'@  Minimum value               = ',E14.5                      ,/,
     &'@  This scalar is a variance (ISCAVR is positive in usini1)  ',/,
     &'@    but the initialization in usiniv has some negative      ',/,
     &'@    values.                                                 ',/,
     &'@                                                            ',/,
     &'@  Verify the initialization in usiniv or the restart file   ',/,
     &'@    In the case where the values read in the restart file   ',/,
     &'@    are incorrect, they may be modified with usiniv         ',/,
     &'@  Verify the variance definition in usini1.                 ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3051 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION          ',/,
     &'@    ********                                                ',/,
     &'@     VARIANCE OUT OF BOUNDS                                 ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Scalar number ',I10,': ',A8                                ,/,
     &'@  Minimum value                = ',E14.5                     ,/,
     &'@    Desired clipping at SCAMIN = ',E14.5                     ,/,
     &'@  Maximum value                = ',E14.5                     ,/,
     &'@    Desired clipping at SCAMAX = ',E14.5                     ,/,
     &'@  This scalar is a variance (ISCAVR is positive in usini1)  ',/,
     &'@    but the initialization in usiniv has some values out    ',/,
     &'@    of the bounds SCAMIN, SCAMAX or lower than 0 and the    ',/,
     &'@    desired clipping mode is ICLVFL(',I10,') = ',I10         ,/,
     &'@                                                            ',/,
     &'@  Verify the initialization in usiniv or the restart file   ',/,
     &'@    In the case where the values read in the restart file   ',/,
     &'@    are incorrect, they may be modified with usiniv         ',/,
     &'@  Verify the variance definition and the clipping mode in   ',/,
     &'@    usini1.                                                 ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3090 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION          ',/,
     &'@    ********                                                ',/,
     &'@                                                            ',/,
     &'@    THE VARIABLES INITIALIZATION IS INCOMPLETE OR           ',/,
     &'@    INCOHERENT WITH THE PARAMETERS VALUE OF THE CALCULATION ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run (',I10,' errors).         ',/,
     &'@                                                            ',/,
     &'@  Refer to the previous warnings for further information.   ',/,
     &'@  Pay attention to the initialization of                    ',/,
     &'@                                the time-step               ',/,
     &'@                                the turbulence              ',/,
     &'@                                the scalars and variances   ',/,
     &'@                                the time averages           ',/,
     &'@                                                            ',/,
     &'@  Verify usiniv or the restart file.                        ',/,
     &'@    In the case where the values read in the restart file   ',/,
     &'@    are incorrect, they may be modified with usiniv         ',/,
     &'@  Verify usini1.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#endif
C
C----
C FIN
C----
C
      RETURN
      END
c@z
