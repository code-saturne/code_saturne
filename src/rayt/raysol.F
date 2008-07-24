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
                        SUBROUTINE RAYSOL
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
C
     &   COFRUA , COFRUB ,
     &   FLURDS , FLURDB ,
C
     &   DTR    , VISCF  , VISCB  ,
     &   DAM    , XAM    ,
     &   DRTP   , SMBRS  , ROVSDT ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     , W10    ,
     &   RDEVEL , RTUSER ,
C
     &   RU     , RUA    ,
     &   SA     ,
     &   QX     , QY     , QZ     ,
     &   QINCID , SNPLUS ,
C
     &   RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
CFONC    -----------------------------------------
CFONC
CFONC    CALCUL DES FLUX ET DU TERME SOURCE RADIATIFS
CFONC
CFONC   1/ DONNEES DES LUMINANCES ENTRANTES AUX LIMITES DU DOMAINE
CFONC         (C.L : REFLEXION ET EMISSION ISOTROPE)
CFONC
CFONC                                ->  ->           ->
CFONC   2/ CALCUL DE LA LUMINANCE L( X , S ) AU POINT X
CFONC
CFONC                                     D L
CFONC      PAR RESOLUTION DE L'EQUATION : --- = -TK.L +TS
CFONC                                     D S
CFONC                         ->                o
CFONC      OU ENCORE : DIV (L.S ) = -TK.L + TK.L
CFONC
CFONC                                   ->   /    ->  ->  ->
CFONC   3/ CALCUL DES DENSITES DE FLUX  Q = /  L( X , S ).S DOMEGA
CFONC                                      /4.PI
CFONC
CFONC                                        /    ->  ->
CFONC          ET DE L'ABSORPTION       SA= /  L( X , S ).  DOMEGA
CFONC                                      /4.PI
CFONC
CFONC      PAR INTEGRATION DES LUMINANCES SUR LES ANGLES SOLIDES.
CFONC
CFONC      N . B : CA SERT A CALCULER LE TAUX D'ECHAUFFEMENT
CFONC      -----
CFONC                                        /    ->  ->  ->  ->
CFONC   4/ CALCUL DU FLUX INCIDENT QINCID = /  L( X , S ).S . N DOMEGA
CFONC                                      /->->
CFONC         ->                          / S.N >0
CFONC         N NORMALE FLUIDE VERS PAROI
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
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE PHASES                     !
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
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS      )!    !     !                                      !
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
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COFRUA,COFRUB! TR !  -  ! CONDITIONS AUX LIMITES AUX           !
CARGU !(NFABOR)      !    !     !    FACES DE BORD POUR LA LUMINANCE   !
CARGU ! FLURDS,FLURDB! TR !  -  ! PSEUDO FLUX DE MASSE (FACES INTERNES !
CARGU !(NFAC)(NFABOR)!    !     !    ET FACES DE BORD )                !
CARGU ! DTR(NCELET)  ! TR !  -  ! DT*CDTVAR                            !
CARGU ! VISCF(NFAC)  ! TR !  -  ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DRTP(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! SMBRS(NCELET ! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! ROVSDT(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! W1...9(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! RU  (NCELET) ! TR !  -  ! LUMINANCE                            !
CARGU ! RUA (NCELET) ! TR !  -  ! LUMINANCE PDT PRECEDENT (NULLE)      !
CARGU ! SA (NCELET)  ! TR ! <-  ! PART D'ABSORPTION DU TERME SOURCE RAD!
CARGU ! QXQYQZ(NCELET! TR ! <-  ! COMPOSANTE DU VECTEUR DENSITE DE FLUX!
CARGU !              !    !     ! RADIATIF EXPLICITE                   !
CARGU ! QINCID(NFABOR! TR ! <-  ! DENSITE DE FLUX RADIATIF AUX BORDS   !
CARGU ! SNPLUS(NFABOR! TR !  -  ! INTEGRATION DU DEMI-ESPACE EGALE A PI!
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
      INCLUDE "radiat.h"
C
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS  , IPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML) , ITYPFB(NFABOR,NPHAS)
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
C
      DOUBLE PRECISION COFRUA(NFABOR), COFRUB(NFABOR)
      DOUBLE PRECISION FLURDS(NFAC), FLURDB(NFABOR)
C
      DOUBLE PRECISION DTR(NCELET)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET), SMBRS(NCELET)
      DOUBLE PRECISION ROVSDT(NCELET)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION W10(NCELET)
C
      DOUBLE PRECISION RU(NCELET), RUA(NCELET)
      DOUBLE PRECISION SA(NCELET)
      DOUBLE PRECISION QX(NCELET), QY(NCELET), QZ(NCELET)
      DOUBLE PRECISION QINCID(NFABOR), SNPLUS(NFABOR)
C
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C
C VARIABLES LOCALES
C
      CHARACTER*80     CNOM
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC  , IEL
      INTEGER          ICONV1, IDIFF1, NDIRC1, IRESO1
      INTEGER          NITMAP, NSWRSP, NSWRGP, IWARNP
      INTEGER          IMGR1 , IMLIGP, IRCFLP, ISCHCP, ISSTPP, IESCAP
      INTEGER          NCYMAP, NITMGP
      INTEGER          IDIR  , NDIRS , KDIR  , IPP   , INUM
      INTEGER          II, JJ, KK, IDTVA0, IVAR0
      DOUBLE PRECISION EPSRGP, BLENCP, CLIMGP, EPSILP, EXTRAP, EPSRSP
      DOUBLE PRECISION SX, SY, SZ, DOMEGA
      DOUBLE PRECISION SXT(NDIRS8), SYT(NDIRS8), SZT(NDIRS8)
      DOUBLE PRECISION AA, SURFBN
      DOUBLE PRECISION RELAXP, THETAP
C
C***********************************************************************
C
C=======================================================================
C 0. GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IVAR0   = 0
      NITMAP  = 1000
c     IMRGRA  = 0
      NSWRSP  = 2
      NSWRGP  = 100
      IMLIGP  = -1
      IRCFLP  = 1
      ISCHCP  = 1
      ISSTPP  = 0
      IESCAP  = 0
      IMGR1   = 0
      NCYMAP  = 100
      NITMGP  = 10
      IWARNP  = IIMLUM
      BLENCP  = ZERO
      EPSILP  = 1.D-8
      EPSRSP  = 1.D-8
      EPSRGP  = 1.D-5
      CLIMGP  = 1.5D0
      EXTRAP  = ZERO
      RELAXP  = 1.D0
C
C--> Il y a des dirichlets
C
      NDIRC1 = 1
C
C--> Convection pure
C
      ICONV1 = 1
C
C--> JACOBI
C
      IRESO1 = 1
C
C=======================================================================
C 2. Choix du nombre de directions et calcul des coordonnees SX,SY,SZ
C=======================================================================
C
C
      NDIRS = NDIREC/8
C
      DOMEGA = PI /(2 *NDIRS)
C
      CALL RAYDIR
C     ===========
     &     (SXT, SYT, SZT, NDIRS)
C
C=======================================================================
C                                              / -> ->
C 3. CORRECTION DES C.L. POUR RESPECTER : PI= /  S. N DOMEGA
C                                            /2PI
C=======================================================================
C
      DO IFAC = 1,NFABOR
        SNPLUS(IFAC) = ZERO
      ENDDO
C
      DO II = -1,1,2
        DO JJ = -1,1,2
          DO KK = -1,1,2
            DO IDIR = 1,NDIRS
C
              SX = II *SXT (IDIR)
              SY = JJ *SYT (IDIR)
              SZ = KK *SZT (IDIR)
C
              DO IFAC = 1,NFABOR
                SURFBN = RA(ISRFBN-1+IFAC)
                AA = SX * SURFBO(1,IFAC)
     &             + SY * SURFBO(2,IFAC)
     &             + SZ * SURFBO(3,IFAC)
                AA = AA / SURFBN
                SNPLUS(IFAC) =SNPLUS(IFAC) +0.5D0 *(-AA+ABS(AA)) *DOMEGA
              ENDDO
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
      DO IFAC = 1,NFABOR
        COFRUA(IFAC) = COFRUA(IFAC) *(PI /SNPLUS(IFAC))
      ENDDO
C
C=======================================================================
C 4. INITIALISATION POUR INTEGRATION DANS LES BOUCLES SUIVANTES
C=======================================================================
C
      DO IFAC = 1, NFABOR
        QINCID(IFAC) = ZERO
        SNPLUS(IFAC) = ZERO
      ENDDO
C
      DO IEL = 1, NCELET
        SA(IEL) = ZERO
        QX(IEL) = ZERO
        QY(IEL) = ZERO
        QZ(IEL) = ZERO
      ENDDO
C
C--> Stockage du SMBRS dans tableau tampon, il sont recharges
C    a chaque changement de direction
C
      DO IEL = 1, NCEL
        W10(IEL) =  SMBRS(IEL)
      ENDDO
C
C--> ROVSDT charge une seule fois
      DO IEL = 1, NCEL
        ROVSDT(IEL) = MAX(ROVSDT(IEL),ZERO)
      ENDDO
C
C--> ON SAUVEGARDE LE NOM DE LA PREMIERE VARIABLE
C    Attention, le passage du nom de variable pour les impressions
C      est fait par common dans NOMVAR (pas tres pratique).
C      Des infos sont egalement remplies en common (residu, nbiter...)
C      Les directions n'etant pas des variables qu'on se donne
C      la possibilite d'imprimer dans le listing, on utilise
C      la position numero 1 qui est une poubelle.
C
      IPP  = 1
      NOMVAR(IPP) = 'RayonXXX'
C
C=======================================================================
C 5. RESOLUTION DE L'EQUATION DES TRANSFERTS RADIATIFS
C=======================================================================
C
C=======================================================================
C 5.1 DISCRETISATION ANGULAIRE
C=======================================================================
C
      KDIR = 0
C
      DO II = -1,1,2
        DO JJ = -1,1,2
          DO KK = -1,1,2
            DO IDIR = 1,NDIRS
C
              SX = II * SXT(IDIR)
              SY = JJ * SYT(IDIR)
              SZ = KK * SZT(IDIR)
C
              KDIR = KDIR + 1
C
              CNOM = ' '
              WRITE(CNOM,'(A5,I3.3)')'Rayon',KDIR
              INUM = IPP
              NOMVAR(INUM) = CNOM
C
C=======================================================================
C 5.2 DISCRETISATION SPATIALE
C=======================================================================
C
C=======================================================================
C 5.1.1 PREPARATION ET PARAMETRAGE DE LA RESOLUTION
C=======================================================================
C
C--> Terme source explicite
C
              DO IEL = 1, NCEL
                SMBRS(IEL) = W10(IEL)
              ENDDO
C
C--> Terme source implicite (ROVSDT vu plus haut)
C
C--> Pas de diffusion facette
C
              IDIFF1 = 0
              DO IFAC = 1,NFAC
                VISCF(IFAC) = ZERO
              ENDDO
              DO IFAC = 1,NFABOR
                VISCB(IFAC) = ZERO
              ENDDO
C
              DO IEL = 1,NCELET
                DRTP(IEL) = ZERO
                RU(IEL)   = ZERO
                RUA(IEL)  = ZERO
              ENDDO
C
              DO IFAC = 1,NFAC
                FLURDS(IFAC) =
     &               + SX*SURFAC(1,IFAC)
     &               + SY*SURFAC(2,IFAC)
     &               + SZ*SURFAC(3,IFAC)
              ENDDO
C
              DO IFAC = 1,NFABOR
                FLURDB(IFAC) =
     &               + SX*SURFBO(1,IFAC)
     &               + SY*SURFBO(2,IFAC)
     &               + SZ*SURFBO(3,IFAC)
              ENDDO
C
C=======================================================================
C 5.1.2 RESOLUTION
C=======================================================================
C
C Dans le cas d'un theta-schema on met theta = 1
C Pas de relaxation en stationnaire non plus
C
              THETAP = 1.0D0
              IDTVA0 = 0
C
              CALL CODITS
C             ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTVA0 , IVAR0  , ICONV1 , IDIFF1 , IRESO1 , NDIRC1 ,  NITMAP ,
     &   IMRGRA , NSWRSP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , IESCAP ,
     &   IMGR1  , NCYMAP , NITMGP , INUM   , IWARNP ,
     &   BLENCP , EPSILP , EPSRSP , EPSRGP , CLIMGP , EXTRAP ,
     &   RELAXP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RUA    , RU     ,
     &   COFRUA , COFRUB , COFRUA , COFRUB , FLURDS , FLURDB ,
     &   VISCF  , VISCB  , VISCF  , VISCB  ,
     &   ROVSDT , SMBRS  , RU    ,
     &   DAM    , XAM    , DRTP   ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 5.2 INTEGRATION DES FLUX ET TERME SOURCE
C=======================================================================
C
              DO IEL = 1,NCEL
                AA = RU(IEL) *DOMEGA
                SA(IEL) = SA(IEL) + AA
                QX(IEL) = QX(IEL) + AA * SX
                QY(IEL) = QY(IEL) + AA * SY
                QZ(IEL) = QZ(IEL) + AA * SZ
              ENDDO
C
C=======================================================================
C 5.3 FLUX INCIDENT A LA PAROI
C=======================================================================
C
              DO IFAC = 1,NFABOR
C
                SURFBN = RA(ISRFBN-1+IFAC)
                AA = SX * SURFBO(1,IFAC)
     &             + SY * SURFBO(2,IFAC)
     &             + SZ * SURFBO(3,IFAC)
                AA = AA / SURFBN
C
                AA = 0.5D0 *(AA+ABS(AA)) *DOMEGA
C
                SNPLUS(IFAC) = SNPLUS(IFAC) + AA
C
                QINCID(IFAC) = QINCID(IFAC) + AA*RU(IFABOR(IFAC))
C
              ENDDO
C
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C
C--------
C FORMATS
C--------
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
