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
                        SUBROUTINE DVVPST
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , NUMMAI ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NVLSTA , NVISBR ,
     &   NCELPS , NFACPS , NFBRPS ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITYPPS , IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   LSTCEL , LSTFAC , LSTFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , STATCE , STATIV , STATFB ,
     &   TRACEL , TRAFAC , TRAFBR , RDEVEL , RTUSER , RA     )
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  ROUTINE DEVELOPPEUR POUR LA SORTIE STANDARD
CFONC  DES VALEURS SUR LES MAILLAGES DE POST TRAITEMENT
CFONC    (appelee apres usproj)
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NUMMAI       ! EC !  -> ! NUMERO DU MAILLAGE POST              !
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
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VARIABLES STAT. LAGRANGIEN !
CARGU ! NVISBR       ! E  !  -> ! NOMBRE DE STATISTIQUES AUX FRONTIERES!
CARGU ! NCELPS       ! E  !  -> ! NOMBRE DE CELLULES DU MAILLAGE POST  !
CARGU ! NFACPS       ! E  !  -> ! NOMBRE DE FACES INTERIEUR POST       !
CARGU ! NFBRPS       ! E  !  -> ! NOMBRE DE FACES DE BORD POST         !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! ITYPPS(3)    ! TE !  -> ! INDICATEUR DE PRESENCE (0 OU 1) DE   !
CARGU !              !    !     ! CELLULES (1), FACES (2), OU FACES DE !
CARGU !              !    !     ! DE BORD (3) DANS LE MAILLAGE POST    !
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
CARGU ! LSTCEL(NCELPS! TE !  -> ! LISTE DES CELLULES DU MAILLAGE POST  !
CARGU ! LSTFAC(NFACPS! TE !  -> ! LISTE DES FACES INTERIEURES POST     !
CARGU ! LSTFBR(NFBRPS! TE !  -> ! LISTE DES FACES DE BORD POST         !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TE !  -  ! MACRO TABLEAU ENTIER                 !
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
CARGU ! (NCELET)     !    !     !                                      !
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
CARGU ! STATCE       ! TR !  -> ! STATISTIQUES CELLULES (LAGRANGIEN)   !
CARGU !(NCELET,NVLSTA!    !     !                                      !
CARGU ! STATIV       ! TR !  -> ! STATISTIQUES VARIANCE CELLULES       !
CARGU !(NCELET,NVLSTA!    !     !                          (LAGRANGIEN)!
CARGU ! STATFB       ! TR !  -> ! STATISTIQUES FACES BORD (LAGRANGIEN) !
CARGU !(NFABOR,NVISBR!    !     !                                      !
CARGU ! TRACEL(*)    ! TR ! <-> ! TAB REEL VALEURS CELLULES POST       !
CARGU ! TRAFAC(*)    ! TR ! <-> ! TAB REEL VALEURS FACES INT. POST     !
CARGU ! TRAFBR(*)    ! TR ! <-> ! TAB REEL VALEURS FACES BORD POST     !
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
      INCLUDE "dimfbr.h"
      INCLUDE "paramx.h"
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "optcal.h"
      INCLUDE "numvar.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "radiat.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NUMMAI
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS  , NVLSTA , NVISBR
      INTEGER          NCELPS , NFACPS , NFBRPS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          ITYPPS(3)
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          LSTCEL(NCELPS), LSTFAC(NFACPS), LSTFBR(NFBRPS)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTPA(NCELET,*), RTP(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION STATCE(NCELET,NVLSTA), STATFB(NFABOR,NVISBR)
      DOUBLE PRECISION STATIV(NCELET,NVLSTA)
      DOUBLE PRECISION TRACEL(NCELPS*3)
      DOUBLE PRECISION TRAFAC(NFACPS*3), TRAFBR(NFBRPS*3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*32     NAMEVR, NAMEV1, NAMEV2
      CHARACTER*80     NAME80
C
      INTEGER          IDEBIA, IDEBRA, IFINIA, IFINRA
      INTEGER          IGRADX, IGRADY, IGRADZ
      INTEGER          ITRAVX, ITRAVY, ITRAVZ, ITRECO
      INTEGER          IW1   , IW2
      INTEGER          INC   , ICCOCG, NSWRGP, IMLIGP, IWARNP, IPHYDP
      INTEGER          ISORVA, ISAUT
      INTEGER          IFAC  , ILOC  , IPHAS , IVAR , ICLVAR
      INTEGER          IRA   , IDIVDT, INEEYP
      INTEGER          IPP   , IDIMT , II    , KK   , IEL
      INTEGER          IVARL , IIP   , IPH
      INTEGER          III, IVARL1 , IVARLM , IFLU   , ILPD1  , ICLA
      INTEGER          ISCAL , IPCVSL, IPCVST, IFLMAB
      INTEGER          IDIMTE, ITENSO, IENTLA, IVARPR
      INTEGER          IYPLBP
      INTEGER          IPCCP
      DOUBLE PRECISION CP0IPH, XCP   , XVSL  , SURFBN, DISTBR
      DOUBLE PRECISION VISCT , FLUMAB, DIIPBX, DIIPBY, DIIPBZ
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP
      DOUBLE PRECISION RBID(1)
C
C
C***********************************************************************
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C     1.1. TRAITEMENT POUR LE MAILLAGE FLUIDE
C=======================================================================
C
      IF (NUMMAI .EQ. -1) THEN
C
C
C       1.1.1 TRAITEMENT DES VARIABLES POST TRAITABLES
C       ----------------------------------------------
C
        DO IPP = 2, NVPPMX
C
C         -> si chrono demande sur la variable
          IF(ICHRVR(IPP).EQ.1) THEN

C           -> pointeur de la variable dans RA (en absolu)
C              (si negatif, c'est un vecteur)
            IRA = IPP2RA(IPP)
C           -> si c'est un moment cumule, il faut diviser par le temps
C             (si   0 ce n'est pas un moment,
C              si > 0 c'est le pointeur dans RA,
C              si < 0 c'est le rang dans DTCMOM)
            IDIVDT = IPPMOM(IPP)
C           -> dimension de la variable a ecrire
            IDIMT = 1
            IF(IRA.LT.0) THEN
              IDIMT = 3
              IRA = -IRA
            ENDIF
C
C           -> nom de la variable
            NAME80 = NOMVAR(IPP)
            NAMEVR = NAME80(1:32)
C
C           -> si l'on a une variable vectorielle, on supprime
C              la partie X, Y, ou Z en derniere ou avant-derniere
C              position dans le nom
C
            IF (IDIMT.EQ.3) THEN
              NAME80 = NOMVAR(IPP+1)
              NAMEV1 = NAME80(1:32)
              NAME80 = NOMVAR(IPP+2)
              NAMEV2 = NAME80(1:32)
              CALL PSTSNV ( NAMEVR , NAMEV1 , NAMEV2 )
C             ===========
            ENDIF
C
C           -> si c'est un moment cumule, il faut diviser par le temps
            IF(IDIVDT.GT.0) THEN
              DO ILOC = 1, NCELPS
                IEL = LSTCEL(ILOC)
                TRACEL(ILOC) = RA(IRA+IEL-1)/
     &               MAX(RA(IDIVDT+IEL-1),EPZERO)
              ENDDO
            ELSEIF(IDIVDT.LT.0) THEN
              DO ILOC = 1, NCELPS
                IEL = LSTCEL(ILOC)
                TRACEL(ILOC) = RA(IRA+IEL-1)/
     &               MAX(DTCMOM(-IDIVDT),EPZERO)
              ENDDO
            ENDIF
C
C           Ecriture effective des valeurs calculees
C
C           Valeurs non entrelacées
            IENTLA = 0
C
            IF(IDIVDT.EQ.0) THEN
              IVARPR = 1
              CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
C             ===========
     &                    NTCABS, TTCABS, RA(IRA), RBID, RBID)
C
            ELSE
              IVARPR = 0
              CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
C             ===========
     &                    NTCABS, TTCABS, TRACEL, RBID, RBID)
C
            ENDIF
C
          ENDIF
C         Fin du traitement en cas de sortie de la variable
C
        ENDDO
C       Fin de la boucle sur les variables
C
C       1.1.2 VARIABLES SUPPLEMENTAIRES AUTOMATIQUES
C       --------------------------------------------
C
C       Distance a la paroi (si LES+VanDriest ou Rij+Echo ou K-w SST)
C
        IF (INEEDY.EQ.1 .AND. ABS(ICDPAR).EQ.1) THEN
C
          NAMEVR = 'DistParoi'
          IDIMT = 1
          IENTLA = 0
          IVARPR = 1
C
          CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
C         ===========
     &                NTCABS, TTCABS, RA(IDIPAR), RBID, RBID)
C
        ENDIF
C
C       Yplus (si LES+VanDriest)
C
        IF (INEEDY.EQ.1 .AND. ABS(ICDPAR).EQ.1) THEN
C
          INEEYP = 0
          DO IPHAS = 1, NPHAS
            IF(ITYTUR(IPHAS).EQ.4.AND.IDRIES(IPHAS).EQ.1) THEN
              INEEYP = 1
            ENDIF
          ENDDO
C
          IF (INEEYP.EQ.1) THEN
C
            NAMEVR = 'Yplus'
            IDIMT = 1
            IENTLA = 0
            IVARPR = 1
C
            CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
C           ===========
     &                  NTCABS, TTCABS, RA(IYPPAR), RBID, RBID)
C
          ENDIF
C
        ENDIF
C
C
C=======================================================================
C     1.2. TRAITEMENT POUR LE MAILLAGE DE BORD
C=======================================================================
C
      ELSE IF  (NUMMAI .EQ. -2) THEN
C
C
C --    1.2.1 TRAITEMENT DE YPLUS AU BORD
C       ----------------------------------
C
        IF(MOD(IPSTDV,IPSTYP).EQ.0) THEN
C
C       Phase
          DO IPHAS = 1, NPHAS
C
C       Initialisation
            DO II = 1, 32
              NAMEVR (II:II) = ' '
            ENDDO
C
C       Nom de la variable
            IF(NPHAS.GT.1) THEN
              WRITE(NAMEVR,1000)IPHAS
 1000         FORMAT('YplusPhase',I2.2)
            ELSE
              NAMEVR = 'Yplus'
            ENDIF
C
C       Dimension de la variable (3 = vecteur, 1=scalaire)
            IDIMT = 1
C
C       Calcul des valeurs de la variable sur les faces de bord
C
            IYPLBP = IYPLBR+(IPHAS-1)*NFABOR
            DO ILOC = 1, NFBRPS
              IFAC = LSTFBR(ILOC)
              TRAFBR(1 + (ILOC-1)*IDIMT) = RA(IYPLBP+IFAC-1)
            ENDDO
C
C           Valeurs non entrelacées, définies sur tableau de travail
            IENTLA = 0
            IVARPR = 0
C
            CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
C           ===========
     &                  NTCABS, TTCABS, RBID, RBID, TRAFBR)
C
          ENDDO
C     fin du test sur les phases
C
        ENDIF
C fin du test sur sortie de yplus
C
C
C --    1.2.2 TRAITEMENT DES VARIABLES AU BORD SANS RECONSTRUCTION
C       -----------------------------------------------------------
C
        IF(MOD(IPSTDV,IPSTCL).EQ.0) THEN
C
C       Le codage ci-dessous est relativement avance :
C         il accede aux variables directement au travers du macro tableau RA
C         il comprend un artifice "ISAUT" permettant de reperer les vecteurs
C           a posttraiter
C
C
C       Boucle sur les variables usuelles
C       ---------------------------------
C
          ISAUT = 0
C
          DO IVAR = 1, NVAR
C
C
C         Variables post-traitables
C         (ISAUT est utilise pour ne pas postraiter plusieurs fois les
C          composantes 2 et 3 d'un vecteur, initialise a 0 avant la boucle
C          sur IVAR)
C
            IPP = IPPRTP(IVAR)
C
            ISORVA = 0
            IF (ISAUT .GT. 0) THEN
              ISAUT = ISAUT - 1
              ISORVA = 0
            ELSE IF (ICHRVR(IPP).EQ.1) THEN
              ISORVA = 1
              NAME80 = NOMVAR(IPP)
              NAMEVR = NAME80(1:32)
            ENDIF
C
C         Traitement des variables definies aux centres cellules à sortir
C         ---------------------------------------------------------------
            IF (ISORVA .EQ. 1) THEN
C
C
C         -> on verifie le signe du pointeur de la variable dans RA
C         (si negatif, c'est un vecteur)
C
C         -> dimension de la variable a ecrire
              IDIMT = 1
              IF(IPP2RA(IPP).LT.0) THEN
                IDIMT = 3
                ISAUT = 2
              ENDIF
C
C           -> si l'on a une variable vectorielle, on supprime
C              la partie X, Y, ou Z en derniere ou avant-derniere
C              position dans le nom
C
            IF (IDIMT.EQ.3) THEN
              NAME80 = NOMVAR(IPP+1)
              NAMEV1 = NAME80(1:32)
              NAME80 = NOMVAR(IPP+2)
              NAMEV2 = NAME80(1:32)
              CALL PSTSNV ( NAMEVR , NAMEV1 , NAMEV2 )
C             ===========
            ENDIF
C
C         Calcul des valeurs (non-reconstruites) de la variable
C         sur les faces de bord
C
              DO KK = 0, IDIMT-1
C
                ICLVAR = ICLRTP(IVAR+KK,ICOEF)
                DO ILOC = 1, NFBRPS
C
                  IFAC = LSTFBR(ILOC)
                  IEL = IFABOR(IFAC)
C
                  TRAFBR(KK + (ILOC-1)*IDIMT + 1)
     &                 =   COEFA(IFAC,ICLVAR)
     &                   + COEFB(IFAC,ICLVAR)*RTP(IEL,IVAR+KK)
C
                ENDDO
C
              ENDDO
C
C             Valeurs entrelacées, définies sur tableau de travail
              IENTLA = 1
              IVARPR = 0
C
              CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
C             ===========
     &                    NTCABS, TTCABS, RBID, RBID, TRAFBR)
C
            ENDIF
C         Fin du traitement en cas de sortie de la variable
C
          ENDDO
C       Fin de la boucle sur les variables
C
        ENDIF
C     Fin du test sur sortie des variables
C
C
C
C       1.2.3 TRAITEMENT FLUX THERMIQUE AU BORD
C       ----------------------------------------
C           Si on travaille en enthalpie, on calcule un flux d'enthalpie
C
        IF(MOD(IPSTDV,IPSTFT).EQ.0) THEN
C
C       Phase
          DO IPHAS = 1, NPHAS
C
            IF(ISCALT(IPHAS).GT.0 .AND. NSCAL.GT.0 .AND.
     &           ISCALT(IPHAS).LE.NSCAL) THEN
C
C       Initialisation
              DO II = 1, 32
                NAMEVR (II:II) = ' '
              ENDDO
C
C       Nom de la variable
              IF(NPHAS.GT.1) THEN
                WRITE(NAMEVR,2000)IPHAS
 2000           FORMAT('Flux th. entrant W.m-2 Phase',I2.2)
              ELSE
                NAMEVR = 'Flux thermique entrant W.m-2'
              ENDIF
C
C       Dimension de la variable (3 = vecteur, 1=scalaire)
              IDIMT = 1
C
C       Numero de la variable
C
              ISCAL  = ISCALT(IPHAS)
              IVAR   = ISCA(ISCAL)
              ICLVAR = ICLRTP(IVAR,ICOEF)
C
C       Calcul des valeurs de la variable sur les faces de bord
C
C          Reservation de la memoire pour reconstruction
C
              IFINIA = IDEBIA
C
              IGRADX = IDEBRA
              IGRADY = IGRADX+NCELET
              IGRADZ = IGRADY+NCELET
              ITRAVX = IGRADZ+NCELET
              ITRAVY = ITRAVX+NCELET
              ITRAVZ = ITRAVY+NCELET
              ITRECO = ITRAVZ+NCELET
              IFINRA = ITRECO+NFABOR
C
C          Verification de la disponibilite de la memoire
C
              CALL IASIZE('DVVPST',IFINIA)
              CALL RASIZE('DVVPST',IFINRA)
C
C
C          Calcul du gradient de la temperature / enthalpie
C
C
C      Pour calculer le gradient de Temperature
C        - dans les calculs paralleles, il est necessaire que
C          les cellules situees sur un bord de sous-domaine connaissent
C          la valeur de temperature dans les cellules situees en
C          vis-a-vis sur le sous-domaine voisin.
C        - dans les calculs periodiques, il est necessaire que
C          les cellules periodiques aient acces a la valeur de la
C          temperature des cellules periodiques correspondantes
C
C      Pour cela, il est necessaire d'appeler les routines de
C        communication PARCOM (parallelisme) et PERCOM (periodicite)
C        pour echanger les valeurs de temperature avant de calculer le
C        gradient. L'appel a ces routines doit etre fait dans cet ordre
C        PARCOM puis PERCOM (pour les cas ou parallelisme et periodicite
C        coexistent).
C      En effet, on se situe ici a la fin du pas de temps n. Or,
C        les variables RTP ne seront echangees qu'en debut du pas de
C        temps n+1. Ici, seules les variables RTPA (obtenues a la fin
C        du pas de temps n-1) ont deja ete echangees.
C
C      Si le calcul n'est ni periodique, ni parallele, on peut conserver
C        appels (les tests sur IPERIO et IRANGP assurent la generalite)
C
C
C          Echange pour le parallelisme
C
              IF(IRANGP.GE.0) THEN
C
                CALL PARCOM (RTP(1,IVAR))
C               ===========
C
              ENDIF
C
C          Echange pour la periodicite
C
              IF(IPERIO.EQ.1) THEN
C
                IDIMTE = 0
                ITENSO = 0
                CALL PERCOM
C               ===========
     &      ( IDIMTE , ITENSO ,
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR))
C
              ENDIF
C
C
C          Calcul du gradient
C
              INC = 1
              ICCOCG = 1
              NSWRGP = NSWRGR(IVAR)
              IMLIGP = IMLIGR(IVAR)
              IWARNP = IWARNI(IVAR)
              EPSRGP = EPSRGR(IVAR)
              CLIMGP = CLIMGR(IVAR)
              EXTRAP = EXTRAG(IVAR)
              IPHYDP = 0
C
              CALL GRDCEL
C             ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA(ITRAVX) , RA(ITRAVX) , RA(ITRAVX) ,
     &   RTP(1,IVAR) , COEFA(1,ICLVAR) , COEFB(1,ICLVAR) ,
     &   RA(IGRADX) , RA(IGRADY) , RA(IGRADZ) ,
C        ----------   ----------   ----------
     &   RA(ITRAVX) , RA(ITRAVY) , RA(ITRAVZ) ,
     &   RDEVEL , RTUSER , RA     )
C
C
C          Calcul de la valeur reconstruite dans les cellules de bord
C
              DO IFAC = 1, NFABOR
                IEL = IFABOR(IFAC)
                III = IDIIPB-1+3*(IFAC-1)
                DIIPBX = RA(III+1)
                DIIPBY = RA(III+2)
                DIIPBZ = RA(III+3)
                RA(ITRECO+IFAC-1) = RTP(IEL,IVAR)
     &               + DIIPBX*RA(IGRADX+IEL-1)
     &               + DIIPBY*RA(IGRADY+IEL-1)
     &               + DIIPBZ*RA(IGRADZ+IEL-1)
              ENDDO
C
C          Calcul du flux (ouf!) convectif et diffusif
C
              IF(IVISLS(ISCAL).GT.0) THEN
                IPCVSL = IPPROC(IVISLS(ISCAL))
              ELSE
                IPCVSL = 0
              ENDIF
              IPCVST = IPPROC(IVISCT(IPHAS))
              IFLMAB = IPPROB(IFLUMA(IVAR))
C
              DO ILOC = 1, NFBRPS
                IFAC = LSTFBR(ILOC)
                IEL = IFABOR(IFAC)
C
                IF(IPCVSL.GT.0) THEN
                  XVSL = PROPCE(IEL,IPCVSL)
                ELSE
                  XVSL = VISLS0(ISCAL)
                ENDIF
                SURFBN = RA(ISRFBN-1+IFAC)
                DISTBR = RA(IDISTB-1+IFAC)
                VISCT  = PROPCE(IEL,IPCVST)
                FLUMAB = PROPFB(IFAC,IFLMAB)
C
                TRAFBR(1 + (ILOC-1)*IDIMT) =
     &               (XVSL+VISCT/SIGMAS(ISCAL))/MAX(DISTBR,EPZERO)*
     &               (COEFA(IFAC,ICLVAR)+(COEFB(IFAC,ICLVAR)-1.D0)*
     &               RTP(IEL,IVAR))
     &               - FLUMAB/MAX(SURFBN,EPZERO**2)*
     &               (COEFA(IFAC,ICLVAR)+ COEFB(IFAC,ICLVAR)*
     &               RTP(IEL,IVAR))
C
              ENDDO

C          Pour la temperature, on multiplie par CP
              IF(ABS(ISCSTH(ISCAL)).EQ.1) THEN
                IF(ICP(IPHAS).GT.0) THEN
                  IPCCP  = IPPROC(ICP   (IPHAS))
                ELSE
                  IPCCP  = 0
                  CP0IPH = CP0(IPHAS)
                ENDIF
                DO ILOC = 1, NFBRPS
                  IFAC = LSTFBR(ILOC)
                  IEL = IFABOR(IFAC)
                  IF(IPCCP.GT.0) THEN
                    XCP = PROPCE(IEL,IPCCP)
                  ELSE
                    XCP    = CP0IPH
                  ENDIF
                  TRAFBR(1 + (ILOC-1)*IDIMT)
     &                 = XCP*TRAFBR(1 + (ILOC-1)*IDIMT)
                ENDDO
              ENDIF
C
C             Valeurs entrelacées, définies sur tableau de travail
              IENTLA = 1
              IVARPR = 0
C
              CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
C             ===========
     &                    NTCABS, TTCABS, RBID, RBID, TRAFBR)
C
            ENDIF
C         Fin du test sur variable thermique
C
          ENDDO
C       Fin de boucle sur les phases
C
        ENDIF
C      Fin du test sur sortie des flux thermiques
C
C --    1.2.4 TRAITEMENT DES EFFORTS AUX BORDS
C       --------------------------------------
C
        IF(MOD(IPSTDV,IPSTFO).EQ.0) THEN
C
C       Initialisation
          DO II = 1, 32
            NAMEVR (II:II) = ' '
          ENDDO
C
C       Nom de la variable
          NAMEVR = 'Efforts'
C
C       Dimension de la variable (3 = vecteur, 1=scalaire)
          IDIMT = 3
C
C       Calcul des valeurs de la variable sur les faces de bord
C
          DO ILOC = 1, NFBRPS
            IFAC = LSTFBR(ILOC)
            SURFBN = RA(ISRFBN-1+IFAC)
            TRAFBR(1 + (ILOC-1)*IDIMT ) =
     &             RA(IFORBR+(IFAC-1)*IDIMT  )/SURFBN
            TRAFBR(2 + (ILOC-1)*IDIMT ) =
     &             RA(IFORBR+(IFAC-1)*IDIMT+1)/SURFBN
            TRAFBR(3 + (ILOC-1)*IDIMT ) =
     &             RA(IFORBR+(IFAC-1)*IDIMT+2)/SURFBN
          ENDDO
C
C           Valeurs entrelacées, définies sur tableau de travail
            IENTLA = 1
            IVARPR = 0
C
            CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
C           ===========
     &                  NTCABS, TTCABS, RBID, RBID, TRAFBR)
C
        ENDIF
C fin du test sur sortie des efforts
C
      ENDIF
C     Fin du test sur le numero de maillage post.
C
C=======================================================================
C     2.1. VARIABLES LAGRANGIENNES
C=======================================================================
C
      IF (NUMMAI .EQ. -1) THEN
C
        IF (IILAGR.GT.0 .AND. ISTALA.GE.1) THEN
C
C         Toutes les statistiques standard sont de dimension 1,
C         et sont definies ou calculees sur tableau de travail
C         de maniere non entrelacee (sans importance de toutes
C         manieres pour une variable scalaire)
C
          IDIMT  = 1
          IENTLA = 0
C
          III = NVLSTA-NVLSTS
C
          DO ICLA  = 0, NBCLST
C
C     -> si IPAS = 0 : statistiques globales
C        si 0 < IPAS =< NBCLST : statistiques par groupe
C
            DO IVARL = 1, NVLSTA
C
              IVARL1 = ICLA*NVLSTA +IVARL
              IVARLM = IVARL1
              ILPD1  = ICLA*NVLSTA +ILPD
              IFLU   = 0
C
              IF (IVARL.LE.III) THEN
                IF (IVARL.EQ.IVARL1) THEN
                  NAME80 = NOMLAG(IVARL)
                ELSE
                  WRITE(NAME80,'(A8,A4,I3)') NOMLAG(IVARL),'_grp',ICLA
                ENDIF
              ELSE IF (NVLSTS.GT.0) THEN
                IF (IVARL.EQ.IVARL1) THEN
                  NAME80 = NOMLAG(ILVU(IVARL-III))
                ELSE
                  WRITE(NAME80,'(A8,A4,I3)')
     &                 NOMLAG(ILVU(IVARL-III)),'_grp',ICLA
                ENDIF
              ENDIF
C
              NAMEVR = NAME80(1:32)
C
              CALL USLAEN
C             ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NVLSTA ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVARL  , IVARL1 , IVARLM , IFLU   , ILPD1  , ICLA   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , STATCE , STATIV , TRACEL ,
     &   RDEVEL , RTUSER , RA     )
C
C           La variable est deja definie sur le maillage volumique
C           global ; on utilise donc l'indirection  (donc IVARPR = 1)
              IVARPR = 1
C
              CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
     &                    NTCABS, TTCABS, TRACEL, RBID, RBID)
            ENDDO
C
            DO IVARL = 1, NVLSTA-1
C
              IVARL1 = ICLA*(NVLSTA-1)+IVARL
              IVARLM = ICLA*NVLSTA+IVARL
              ILPD1  = ICLA*NVLSTA +ILPD
              IFLU   = 1
C
              IF (IVARL.LE.III) THEN
                IF (IVARL.EQ.IVARL1) THEN
                  NAME80 = NOMLAV(IVARL)
                ELSE
                  WRITE(NAME80,'(A8,A4,I3)') NOMLAV(IVARL),'_grp',ICLA
                ENDIF
              ELSE IF (NVLSTS.GT.0) THEN
                IF (IVARL.EQ.IVARL1) THEN
                  NAME80 = NOMLAV(ILVU(IVARL-III))
                ELSE
                  WRITE(NAME80,'(A8,A4,I3)')
     &                 NOMLAV(ILVU(IVARL-III)),'_grp',ICLA
                ENDIF
              ENDIF
C
              NAMEVR = NAME80(1:32)
C
              CALL USLAEN
C             ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NVLSTA ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVARL  , IVARL1 , IVARLM , IFLU   , ILPD1  , ICLA   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , STATCE , STATIV , TRACEL ,
     &   RDEVEL , RTUSER , RA     )
C
C           La variable est deja definie sur le maillage volumique
C           global ; on utilise donc l'indirection  (donc IVARPR = 1)
              IVARPR = 1
C
              CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
     &                    NTCABS, TTCABS, TRACEL, RBID, RBID)
            ENDDO
C
          ENDDO
C
        ENDIF
C
      ENDIF
C
      IF (NUMMAI.EQ.-2) THEN
C
        IF (IILAGR.GT.0 .AND. IENSI3.EQ.1) THEN
C
          III = NVISBR-NUSBOR
C
          DO IVARL = 1,NVISBR
C
            IF (IVARL.LE.III) THEN
              NAME80 = NOMBRD(IVARL)
            ELSE IF (NUSBOR.GT.0) THEN
              NAME80 = NOMBRD(IUSB(IVARL-III))
            ENDIF
            NAMEVR = NAME80(1:32)
C
            IF (IMOYBR(IVARL).EQ.2) THEN
C
              DO ILOC = 1, NFBRPS
                IFAC = LSTFBR(ILOC)
                IF (STATFB(IFAC,INBR).GT.SEUILF) THEN
                  TRAFBR(ILOC) = STATFB(IFAC,IVARL)/STATFB(IFAC,INBR)
                ELSE
                  TRAFBR(ILOC) = 0.D0
                ENDIF
              ENDDO
C
            ELSE IF (IMOYBR(IVARL).EQ.1) THEN
C
              DO ILOC = 1, NFBRPS
                IFAC = LSTFBR(ILOC)
                IF (STATFB(IFAC,INBR).GT.SEUILF) THEN
                  TRAFBR(ILOC) = STATFB(IFAC,IVARL) / TSTATP
                ELSE
                  TRAFBR(ILOC) = 0.D0
                ENDIF
              ENDDO
C
            ELSE
C
              DO ILOC = 1, NFBRPS
                IFAC = LSTFBR(ILOC)
                IF (STATFB(IFAC,INBR).GT.SEUILF) THEN
                  TRAFBR(ILOC) = STATFB(IFAC,IVARL)
                ELSE
                  TRAFBR(ILOC) = 0.D0
                ENDIF
              ENDDO
C
            ENDIF
C
            IDIMT  = 1
            IENTLA = 0
            IVARPR = 0
C
            CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
     &                  NTCABS, TTCABS, RBID, RBID, TRAFBR)
C
          ENDDO
C
          NAME80 = 'lagrangian_boundary_zones'
          NAMEVR = NAME80(1:32)
C
          DO ILOC = 1, NFBRPS
            IFAC = LSTFBR(ILOC)
            TRAFBR(ILOC) = IA(IIFRLA+IFAC-1)
          ENDDO
C
          IDIMT  = 1
          IENTLA = 0
          IVARPR = 0
C
          CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
     &         NTCABS, TTCABS, RBID, RBID, TRAFBR)
C
        ENDIF
      ENDIF
C     Fin du test sur le numero de maillage post.
C
C=======================================================================
C     2.2. VARIABLES RAYONNEMENT
C=======================================================================
C
      IF (IIRAYO.EQ.1) THEN
C
        DO IIP = 1, NPHAST
          IPH   = IIP
          IPHAS = IRAPHA(IPH)
          DO IVARL = 1,NBRAYP
C
            IF (IRAYVP(IVARL,IPHAS).EQ.1) THEN
C
              CALL RAYENS
C             ===========
     & ( IFINIA , IFINRA , NUMMAI ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCELPS , NFACPS , NFBRPS ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITYPPS , IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   LSTCEL , LSTFAC , LSTFBR ,
     &   IDEVEL , ITUSER ,
     &   IVARL  , IPH    ,
     &   IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   TRACEL , TRAFAC , TRAFBR , RDEVEL , RTUSER , RA     )
C
            ENDIF
C
          ENDDO
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C     2.3. VARIABLES ELECTRIQUES
C=======================================================================
C
      IF (     IPPMOD(IELJOU).GE.1
     &    .OR. IPPMOD(IELARC).GE.1
     &    .OR. IPPMOD(IELION).GE.1) THEN
C
        IFINIA = IDEBIA
C
        IW1    = IDEBRA
        IW2    = IW1    + NCELET*3
        IFINRA = IW2    + NCELET*3
C
        CALL RASIZE ('DVVPST', IFINRA)
C       ===========
C
        CALL USELEN
C       ===========
     & ( IFINIA , IFINRA , NUMMAI ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCELPS , NFACPS , NFBRPS ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   LSTCEL , LSTFAC , LSTFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RA(IW1) , RA(IW2) ,
     &   TRACEL , TRAFAC , TRAFBR , RDEVEL , RTUSER , RA     )
C
      ENDIF
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
      END
c@z
