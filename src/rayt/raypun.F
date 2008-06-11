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
                        SUBROUTINE RAYPUN
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
     &   COFRUA , COFRUB ,
     &   FLURDS , FLURDB ,
     &   DTR    , VISCF  , VISCB  ,
     &   DAM    , XAM    ,
     &   DRTP   , SMBRS  , ROVSDT ,
     &   THETA4 , THETAA , SA     ,
     &   QX     , QY     , QZ     ,
     &   QINCID , EPS    , TPAROI ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     , CKMEL  ,
     &   RDEVEL , RTUSER , RA     )
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
CFONC    AVEC L'APPROXIMATION P-1
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
CARGU ! THETA4(NCELET! TR !  -  ! PSEUDO TEMPERATURE RADIATIVE         !
CARGU ! THETAA(NCELET! TR !  -  ! PSEUDO TEMP RAR PDT PRECEDENT (NULLE)!
CARGU ! SA (NCELET)  ! TR ! <-  ! PART D'ABSORPTION DU TERME SOURCE RAD!
CARGU ! QXQYQZ(NCELET! TR ! <-  ! COMPOSANTE DU VECTEUR DENSITE DE FLUX!
CARGU !              !    !     ! RADIATIF EXPLICITE                   !
CARGU ! QINCID(NFABOR! TR ! <-  ! DENSITE DE FLUX RADIATIF AUX BORDS   !
CARGU ! EPS (NFABOR) ! TR !  -> ! EMISSIVITE DES FACETTES DE BORD      !
CARGU ! TPAROI(NFABOR! TR !  -> ! TEMPERATURE DE PAROI EN KELVIN       !
CARGU ! W1...9(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! CKMEL(NCELET)! TR !  -> ! COEFF D'ABSORPTION DU MELANGE        !
CARGU !              !    !     !   GAZ-PARTICULES DE CHARBON          !
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
      INCLUDE "parall.h"
      INCLUDE "period.h"
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
C
      DOUBLE PRECISION THETA4(NCELET), THETAA(NCELET)
      DOUBLE PRECISION SA(NCELET)
      DOUBLE PRECISION QX(NCELET), QY(NCELET), QZ(NCELET)
      DOUBLE PRECISION QINCID(NFABOR), TPAROI(NFABOR), EPS(NFABOR)
C
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION CKMEL(NCELET)
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
      INTEGER          INUM
      INTEGER          IDTVA0, IVAR0
      INTEGER          INC, ICCOCG, IPHYDP
      INTEGER          IDIMTE , ITENSO
      DOUBLE PRECISION EPSRGP, BLENCP, CLIMGP, EPSILP, EXTRAP
      DOUBLE PRECISION AA, AAA, AAAA, RELAXP, THETAP
C
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
C 1. PARAMETRAGE DU SOLVEUR ET INITIALISATION
C=======================================================================
C
C--> Gradient Conjugue
C
      IRESO1 = 0
C
C--> Parametrage de CODITS
C
C IVAR0= 0  LA VARIABLE N'EST ICI NI RIJ NI VITESSE
      IVAR0   = 0
      NITMAP  = 1000
c     IMRGRA  = 0
      NSWRSP  = 1
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
      EPSRGP  = 1.D-5
      CLIMGP  = 1.5D0
      EXTRAP  = ZERO
      RELAXP  = 1.D0
C
C--> Il y a des dirichlets
C
      NDIRC1 = 1
C
C--> Pas de convection pour le modele P1
C
      ICONV1 = 0
C
C--> Equation de diffusion
C
      IDIFF1 = 1
C
C--> Remise a zero des tableaux avant resolution
C
      DO IEL = 1,NCEL
        DRTP(IEL)   = ZERO
        THETA4(IEL) = ZERO
        THETAA(IEL) = ZERO
      ENDDO
C
      DO IFAC = 1,NFAC
        FLURDS(IFAC) = ZERO
      ENDDO
C
      DO IFAC = 1,NFABOR
        FLURDB(IFAC) = ZERO
      ENDDO
C
C=======================================================================
C 2. COEFFICIENT DE DIFFUSION AUX FACES
C=======================================================================
C
      DO IEL = 1,NCEL
        CKMEL(IEL) = 1.D0 / CKMEL(IEL)
      ENDDO
C
      CALL VISCFA
C     ===========
     &   ( IDEBIA , IDEBRA ,
     &     NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &     NPRFML , NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &     NIDEVE , NRDEVE , NITUSE , NRTUSE , IMVISF ,
     &     IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &     IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &     IDEVEL , ITUSER , IA     ,
     &     XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD ,
     &     VOLUME , CKMEL  , VISCF  , VISCB  ,
     &     RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 3.  RESOLUTION
C=======================================================================
C
C     Parametre pour schemas en temps et stationnaire
      THETAP = 1.D0
      IDTVA0 = 0
C
      CNOM = ' '
      WRITE(CNOM,'(A)') 'Rayon P1'
      INUM = 1
      NOMVAR(INUM) = CNOM
C
      CALL CODITS
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTVA0 , IVAR0  , ICONV1 , IDIFF1 , IRESO1 , NDIRC1 , NITMAP ,
     &   IMRGRA , NSWRSP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , IESCAP ,
     &   IMGR1  , NCYMAP , NITMGP , INUM   , IWARNP ,
     &   BLENCP , EPSILP , EPSRGP , CLIMGP , EXTRAP , RELAXP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD ,
     &   VOLUME ,
     &   THETAA , THETAA , COFRUA , COFRUB , COFRUA , COFRUB ,
     &   FLURDS , FLURDB ,
     &   VISCF  , VISCB  , VISCF  , VISCB  ,
     &   ROVSDT , SMBRS  , THETA4 ,
     &   DAM    , XAM    , DRTP   ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 4. Vecteur densite de flux radiatif
C=======================================================================
C
C    En periodique et parallele, echange avant calcul du gradient
C
C    Parallele
      IF (IRANGP.GE.0) THEN
        CALL PARCOM (THETA4)
C       ===========
      ENDIF
C
C    Periodique
      IF (IPERIO.EQ.1) THEN
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     &  ( IDIMTE , ITENSO ,
     &    THETA4 , THETA4 , THETA4 ,
     &    THETA4 , THETA4 , THETA4 ,
     &    THETA4 , THETA4 , THETA4)
      ENDIF
C
C     Calcul de la densite du flux radiatif QX, QY, QZ
C
      INC     = 1
      ICCOCG  = 1
      IMLIGP  = -1
      IWARNP  = IIMLUM
      EPSRGP  = 1.D-8
      CLIMGP  = 1.5D0
      EXTRAP  = 0.D0
      NSWRGP  = 100
      IVAR0   = 0
      IPHYDP  = 0
C
      CALL GRDCEL
C     ===========
     &   ( IDEBIA , IDEBRA ,
     &     NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML,
     &     NPRFML ,
     &     NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &     NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &     IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP,
     &     IPHYDP ,
     &     IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &     IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &     IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &     IDEVEL , ITUSER , IA     ,
     &     XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD ,
     &     VOLUME ,
     &     W7     , W7     , W7     ,
     &     THETA4 , COFRUA , COFRUB ,
     &     W1     , W2     , W3     ,
     &     W4     , W5     , W6     ,
     &     RDEVEL , RTUSER , RA     )
C
      AA = - STEPHN * 4.D0 / 3.D0
C
      DO IEL = 1,NCEL
        AAA = AA * CKMEL(IEL)
        QX(IEL) = W1(IEL) * AAA
        QY(IEL) = W2(IEL) * AAA
        QZ(IEL) = W3(IEL) * AAA
      ENDDO
C
C=======================================================================
C 5. Terme Source Radiatif d'absorption et densite de flux incident
C=======================================================================
C
C     Calcul de la part d'absorption du terme Source Radiatif
C
      AA = 4.D0 * STEPHN
      DO IEL = 1,NCEL
        SA(IEL) = AA * THETA4(IEL)
      ENDDO
C
C     Calcul du flux incident Qincid
C
      DO IFAC = 1, NFABOR
C
        IF (ITYPFB(IFAC,IPHAS).EQ.IPAROI .OR.
     &      ITYPFB(IFAC,IPHAS).EQ.IPARUG ) THEN
C
C--> Premiere version plus chere et legerement plus precise
C
          AAAA = TPAROI(IFAC)**4
C
          AAA  = 1.5D0 * RA(IDISTB-1+IFAC) / CKMEL(IFABOR(IFAC))
     &           * ( 2.D0 /(2.D0-EPS(IFAC)) -1.D0 )
          AA   = ( AAA * AAAA + THETA4(IFABOR(IFAC)) )
     &         / (1.D0 + AAA)
C
          QINCID(IFAC) = STEPHN * (2.D0 * AA - EPS(IFAC) * AAAA)
     &                       / (2.D0 - EPS(IFAC))
C
C--> Deuxieme version plus cheap mais moins precise
C
c         QINCID(IFAC) = STEPHN *
c    &    (2.D0 * THETA4(IFABOR(IFAC)) - EPS(IFAC) * TPAROI(IFAC)**4)
c    &  / (2.D0 - EPS(IFAC))
C
        ELSE
          QINCID(IFAC) = STEPHN * THETA4(IFABOR(IFAC))
     &               + ( QX(IEL) * SURFBO(1,IFAC) +
     &                   QY(IEL) * SURFBO(2,IFAC) +
     &                   QZ(IEL) * SURFBO(3,IFAC) ) /
     &                   (0.5D0 * RA(ISRFBN-1+IFAC) )
        ENDIF
C
      ENDDO
C
C=======================================================================
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
