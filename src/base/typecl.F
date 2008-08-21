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
                        SUBROUTINE TYPECL
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , ITRIFB , ICODCL , ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RCODCL , FRCXT  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC TRAITEMENT DES CODES DE CONDITION LIMITE ITYPFB
CFONC
c@fonce
C-----------------------------------------------------------------------
C                            ARGUMENTS
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
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS)      !    !     !                                      !
CARGU ! ITRIFB(NFABOR! TE ! <-  ! TAB D'INDIRECTION POUR TRI DES FACES !
CARGU !  NPHAS)      !    !     !                                      !
CARGU ! ICODCL       ! TE ! <-> ! CODE DE CONDITION LIMITES AUX FACES  !
CARGU !  (NFABOR,NVAR!    !     !  DE BORD                             !
CARGU !              !    !     ! = 1   -> DIRICHLET                   !
CARGU !              !    !     ! = 3   -> DENSITE DE FLUX             !
CARGU !              !    !     ! = 4   -> GLISSEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 5   -> FROTTEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 6   -> RUGOSITE ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 9   -> ENTREE/SORTIE LIBRE (VITESSE!
CARGU !              !    !     !  ENTRANTE EVENTUELLE     BLOQUEE     !
CARGU ! ISOSTD       ! TE ! <-  ! INDICATEUR DE SORTIE STANDARD        !
CARGU !    (NFABOR+1)!    !     !  +NUMERO DE LA FACE DE REFERENCE     !
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
CARGU ! RCODCL       ! TR ! <-  ! VALEUR DES CONDITIONS AUX LIMITES    !
CARGU !  (NFABOR,NVAR!    !     !  AUX FACES DE BORD                   !
CARGU !              !    !     ! RCODCL(1) = VALEUR DU DIRICHLET      !
CARGU !              !    !     ! RCODCL(2) = VALEUR DU COEF. D'ECHANGE!
CARGU !              !    !     !  EXT. (INFINIE SI PAS D'ECHANGE)     !
CARGU !              !    !     ! RCODCL(3) = VALEUR DE LA DENSITE DE  !
CARGU !              !    !     !  FLUX (NEGATIF SI GAIN) W/m2 OU      !
CARGU !              !    !     !  HAUTEUR DE RUGOSITE (m) si ICODCL=6 !
CARGU !              !    !     ! POUR LES VITESSES (VISTL+VISCT)*GRADU!
CARGU !              !    !     ! POUR LA PRESSION             DT*GRADP!
CARGU !              !    !     ! POUR LES SCALAIRES                   !
CARGU !              !    !     !        CP*(VISCLS+VISCT/SIGMAS)*GRADT!
CARGU ! FRCXT(NCELET,! TR !  -> ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU !   3,NPHAS)   !    !     !  HYDROSTATIQUE                       !
CARGU ! W1,2,3,4,5,6 ! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU !  (NCELET     !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
CARGU ! RIJIPB       ! TR !  -  ! TAB DE TRAV POUR VALEURS EN IPRIME   !
CARGU ! (NFABOR,6   )!    !     !  DES RIJ AU BORD                     !
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
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "pointe.h"
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
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICODCL(NFABOR,NVAR)
      INTEGER          ITYPFB(NFABOR,NPHAS) , ITRIFB(NFABOR,NPHAS)
      INTEGER          ISOSTD(NFABOR+1,NPHAS)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NDIMFB,*)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION RCODCL(NFABOR,NVAR,3)
      DOUBLE PRECISION FRCXT(NCELET,3,NPHAS)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION W4(NCELET),W5(NCELET),W6(NCELET)
      DOUBLE PRECISION COEFU(NFABOR,3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER        CHAINE*80
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC, IVAR, IEL
      INTEGER          IOK, INC, ICCOCG, IDEB, IFIN, INB, ISUM, IWRNP
      INTEGER          IFRSLB(NPHSMX), ITBSLB(NPHSMX)
      INTEGER          ITYP, II, JJ, IPHAS, IWARU, IFLMAB
      INTEGER          NSWRGP, IMLIGP, IWARNP
      INTEGER          IPRIPH, IUIPH, IVIPH, IWIPH
      INTEGER          IR11IP, IR22IP, IR33IP, IR12IP, IR13IP,IR23IP
      INTEGER          IKIPH , IEPIPH, IPHIPH, IFBIPH, IOMGIP
      INTEGER          IPRNEW, KPHAS, III
      INTEGER          IRANGD, ICLIPR, IIPTOT
      INTEGER          ITROIS, IFADIR, NFADIR
      DOUBLE PRECISION PREF, EPSRGP, CLIMGP, EXTRAP, COEFUP
      DOUBLE PRECISION DIIPBX, DIIPBY, DIIPBZ
      DOUBLE PRECISION FLUMBF, FLUMTY(NTYPMX)
      DOUBLE PRECISION RO0IPH, P0IPH, PR0IPH, XXP0, XYP0, XZP0
      DOUBLE PRECISION XYZREF(3)
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
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
C=======================================================================
C 2.  VERIFICATION DE LA CONSISTANCE DES TYPES DONNES DANS USCLIM
C=======================================================================
C
      IOK = 0
C
      DO IPHAS = 1, NPHAS
        DO IFAC = 1, NFABOR
          ITYP = ITYPFB(IFAC,IPHAS)
          IF(ITYP.LE.0.OR.ITYP.GT.NTYPMX) THEN
            WRITE(NFECRA,1010)IPHAS, IFAC, IPRFML(IFMFBR(IFAC),1), ITYP
            IOK = IOK + 1
          ENDIF
        ENDDO
      ENDDO
C
      IF (IRANGP.GE.0) CALL PARCMX(IOK)
      IF(IOK.NE.0) THEN
        WRITE(NFECRA,1099) NTYPMX
        CALL CSEXIT (1)
      ENDIF
C
C
C=======================================================================
C 3.  TRI DES FACES DE BORD
C=======================================================================
C
C
C ---> On compte les faces de chaque type (provisoirement dans IFINTY)
C
      DO IPHAS = 1, NPHAS
        DO II = 1, NTYPMX
         IFINTY(II,IPHAS) = 0
       ENDDO
      ENDDO
C
      DO IPHAS = 1, NPHAS
        DO IFAC = 1, NFABOR
          ITYP = ITYPFB(IFAC,IPHAS)
          IFINTY(ITYP,IPHAS) = IFINTY(ITYP,IPHAS) + 1
        ENDDO
      ENDDO
C
C
C ---> On positionne le debut de chaque groupe de face
C       dans ITRIFB (classees par type) : IDEBTY
C
      DO IPHAS = 1, NPHAS
        DO II = 1, NTYPMX
          IDEBTY(II,IPHAS) = 1
        ENDDO
      ENDDO
C
      DO IPHAS = 1, NPHAS
        DO II = 1, NTYPMX-1
          DO JJ = II+1, NTYPMX
            IDEBTY(JJ,IPHAS) = IDEBTY(JJ,IPHAS) + IFINTY(II,IPHAS)
          ENDDO
        ENDDO
      ENDDO
C
C ---> On classe les faces dans ITRIFB
C       et on en profite pour remplir correctement IFINTY
C
      DO IPHAS = 1, NPHAS
        DO II = 1, NTYPMX
          IFINTY(II,IPHAS) = IDEBTY(II,IPHAS)-1
        ENDDO
      ENDDO
C
      DO IPHAS = 1, NPHAS
        DO IFAC = 1, NFABOR
          ITYP = ITYPFB(IFAC,IPHAS)
          IFIN = IFINTY(ITYP,IPHAS)+1
          ITRIFB(IFIN,IPHAS) = IFAC
          IFINTY(ITYP,IPHAS) = IFIN
        ENDDO
      ENDDO
C
C ---> On verifie sommairement
C
      IOK = 0
      DO IPHAS = 1, NPHAS
        DO II = 1, NTYPMX-1
          IF(IFINTY(II,IPHAS).GE.IDEBTY(II+1,IPHAS)) THEN
            IF (IOK.EQ.0) IOK = II
          ENDIF
        ENDDO
        IF (IRANGP.GE.0) CALL PARCMX(IOK)
        IF (IOK.GT.0) THEN
          II = IOK
          WRITE(NFECRA,2010) IPHAS
          IF(IFINTY(II,IPHAS).GE.IDEBTY(II+1,IPHAS)) THEN
            WRITE(NFECRA,2020) (IFINTY(JJ,IPHAS),JJ=1,NTYPMX)
            WRITE(NFECRA,2030) (IDEBTY(JJ,IPHAS),JJ=1,NTYPMX)
            WRITE(NFECRA,2040) (ITYPFB(JJ,IPHAS),JJ=1,NFABOR)
            WRITE(NFECRA,2098)
     &                       II,IFINTY(II,IPHAS),II+1,IDEBTY(II+1,IPHAS)
          ELSE
            WRITE(NFECRA,2099) II,II+1
          ENDIF
          CALL CSEXIT (1)
        ENDIF
      ENDDO
C
      IOK = 0
      DO IPHAS = 1, NPHAS
        ISUM = 0
        DO II = 1, NTYPMX
          ISUM = ISUM + IFINTY(II,IPHAS) - IDEBTY(II,IPHAS) + 1
        ENDDO
        IF (IRANGP.GE.0) CALL PARCPT (ISUM)
        IF(ISUM.NE.NFBRGB) THEN
          WRITE(NFECRA,3099) IPHAS, ISUM, NFBRGB
          IOK = IOK + 1
        ENDIF
      ENDDO
      IF (IOK.NE.0) THEN
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C
C ---> On ecrit les types de faces avec la borne inf et sup et le nb
C       pour chaque type de face trouve (tjrs pour les types par defaut)
C
      IF(IPASS.EQ.0.OR.IWARNI(IU(1)).GE.2) THEN
C
        IPASS = 1
C
        WRITE(NFECRA,6010)
C
        DO IPHAS = 1, NPHAS
C
          WRITE(NFECRA,6011)IPHAS
C
          IF ( IPPMOD(ICOMPF).LT.0 ) THEN
C
            II = IENTRE
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Entree           ',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
            II = IPAROI
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Paroi lisse      ',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
            II = IPARUG
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Paroi rugueuse   ',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
            II = ISYMET
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Symetrie         ',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
C
             II = ISOLIB
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Sortie libre     ',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
C
            II = IINDEF
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Indefini         ',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
C
            DO II = 1, NTYPMX
              IF( II.NE.IENTRE .AND.
     &             II.NE.IPAROI .AND.
     &             II.NE.IPARUG .AND.
     &             II.NE.ISYMET .AND.
     &             II.NE.ISOLIB .AND.
     &             II.NE.IINDEF ) THEN
                INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
                IF (IRANGP.GE.0) CALL PARCPT (INB)
                IF(INB.GT.0) THEN
                  WRITE(NFECRA,6020) 'Type utilisateur ',II,
     &                 INB,
     &                 IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &                 IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
                ENDIF
              ENDIF
            ENDDO
C
          ELSE
C
            II = IEQHCF
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Entree sub. enth.',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
C
            II = IERUCF
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Entree subsonique',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
C
            II = IESICF
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Entree/Sortie imp',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
C
            II = ISOPCF
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Sortie subsonique',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
C
            II = ISSPCF
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Sortie supersoniq',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
C
            II = IPAROI
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Paroi lisse      ',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
C
            II = IPARUG
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Paroi rugueuse   ',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
C
            II = ISYMET
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Symetrie         ',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
C
            II = IINDEF
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) CALL PARCPT (INB)
            WRITE(NFECRA,6020) 'Indefini         ',II,
     &           INB,
     &           IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &           IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
C
            DO II = 1, NTYPMX
              IF( II.NE.IESICF .AND.
     &             II.NE.ISSPCF .AND.
     &             II.NE.IEQHCF .AND.
     &             II.NE.IERUCF .AND.
     &             II.NE.ISOPCF .AND.
     &             II.NE.IPAROI .AND.
     &             II.NE.IPARUG .AND.
     &             II.NE.ISYMET .AND.
     &             II.NE.IINDEF ) THEN
                INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
                IF (IRANGP.GE.0) CALL PARCPT (INB)
                IF(INB.GT.0) THEN
                  WRITE(NFECRA,6020) 'Type utilisateur ',II,
     &                 INB,
     &                 IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1,
     &                 IDEBTY(II,IPHAS),IFINTY(II,IPHAS)
                ENDIF
              ENDIF
            ENDDO
C
          ENDIF
C
          WRITE(NFECRA,6030)
C
        ENDDO
C
      ENDIF
C
C========================================================================
C 4.  RCODCL(.,.,1) A ETE INITIALISE A RINFIN POUR VERIFIER CEUX QUE
C     L'UTILISATEUR A MODIFIE. ON REINITIALISE ICI A ZERO CEUX QUI N'ONT
C     PAS ETE MODIFIES.
C     ISOLIB et IENTRE  SONT  TRAITES PLUS TARD.
C
C
C========================================================================
C
C
      DO IPHAS = 1, NPHAS
        DO IVAR=1, NVAR
           DO IFAC = 1, NFABOR
              IF((ITYPFB(IFAC,IPHAS).NE.ISOLIB) .AND.
     &             (ITYPFB(IFAC,IPHAS).NE. IENTRE) .AND.
     &          (RCODCL(IFAC,IVAR,1).GT.RINFIN*0.5D0)) THEN
                     RCODCL(IFAC,IVAR,1) = 0.D0
             ENDIF
            ENDDO
         ENDDO
      ENDDO
C
C
C
C
C=======================================================================
C 5.  CALCUL DE LA PRESSION AU BORD (DANS COEFU(*,1))
C   (SI ON EN A BESOIN, IE S'IL Y A DES FACES DE BORD SORTIE)
C
C    LA BOUCLE SUR LES PHASES COMMENCEE ICI SE TERMINE A LA FIN DU BLOC
C    SUIVANT
C=======================================================================
C
C --- Boucle sur les phases : debut
      DO IPHAS = 1, NPHAS
C
        RO0IPH = RO0  (IPHAS)
        P0IPH  = P0   (IPHAS)
        PR0IPH = PRED0(IPHAS)
        XXP0   = XYZP0(1,IPHAS)
        XYP0   = XYZP0(2,IPHAS)
        XZP0   = XYZP0(3,IPHAS)
        IPRIPH = IPR (IPHAS)
        IUIPH  = IU  (IPHAS)
        IVIPH  = IV  (IPHAS)
        IWIPH  = IW  (IPHAS)
        IF(ITYTUR(IPHAS).EQ.2) THEN
          IKIPH  = IK(IPHAS)
          IEPIPH = IEP(IPHAS)
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
          IR11IP = IR11(IPHAS)
          IR22IP = IR22(IPHAS)
          IR33IP = IR33(IPHAS)
          IR12IP = IR12(IPHAS)
          IR13IP = IR13(IPHAS)
          IR23IP = IR23(IPHAS)
          IEPIPH = IEP(IPHAS)
        ELSEIF(ITURB(IPHAS).EQ.50) THEN
          IKIPH  = IK(IPHAS)
          IEPIPH = IEP(IPHAS)
          IPHIPH = IPHI(IPHAS)
          IFBIPH = IFB(IPHAS)
        ELSEIF(ITURB(IPHAS).EQ.60) THEN
          IKIPH  = IK (IPHAS)
          IOMGIP = IOMG(IPHAS)
        ENDIF
C
C   -- On regarde si la pression (unique) n'a pas ete traitee
C         auparavant
C
        IPRNEW = 1
        IF(IPHAS.GT.1) THEN
          DO KPHAS = 1, IPHAS-1
            IF(IPR(IPHAS).EQ.IPR(KPHAS)) THEN
              IPRNEW = 0
            ENDIF
          ENDDO
        ENDIF
C
C     IFRSLB = premiere face de sortie libre standard (ICODCL non modifie)
C     ITBSLB = max de IFRSLB sur tous les processeurs, indique s'il y a des
C              faces de sortie standard
        IFRSLB(IPHAS) = 0
        DO II = IFINTY(ISOLIB,IPHAS), IDEBTY(ISOLIB,IPHAS), -1
          IFAC = ITRIFB(II,IPHAS)
          IF (ICODCL(IFAC,IPRIPH).EQ.0)
     &         IFRSLB(IPHAS) = IFAC
        ENDDO
        ITBSLB(IPHAS) = IFRSLB(IPHAS)
        IF (IRANGP.GE.0) THEN
           CALL PARCMX (ITBSLB(IPHAS))
        ENDIF
C
        IF ((ITBSLB(IPHAS).GT.0) .AND. (IPRNEW.EQ.1)) THEN
C
          INC = 1
          ICCOCG = 1
          NSWRGP = NSWRGR(IPRIPH)
          IMLIGP = IMLIGR(IPRIPH)
          IWARNP = IWARNI(IPRIPH)
          EPSRGP = EPSRGR(IPRIPH)
          CLIMGP = CLIMGR(IPRIPH)
          EXTRAP = EXTRAG(IPRIPH)
          ICLIPR = ICLRTP(IPRIPH,ICOEF)
C
          CALL GRDCEL
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPRIPH , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDR ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   FRCXT(1,1,IPHAS), FRCXT(1,2,IPHAS), FRCXT(1,3,IPHAS),
     &   RTPA(1,IPRIPH)  , COEFA(1,ICLIPR) , COEFB(1,ICLIPR) ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C On met dans COEFU la valeur en I' ou F (suivant IPHYDR) de la
C  pression totale, calculee a partir de P*
C
          IF (IPHYDR.EQ.0) THEN
            DO IFAC = 1, NFABOR
              II = IFABOR(IFAC)
              III = IDIIPB-1+3*(IFAC-1)
              DIIPBX = RA(III+1)
              DIIPBY = RA(III+2)
              DIIPBZ = RA(III+3)
              COEFU(IFAC,1) = RTPA(II,IPRIPH)
     &             + DIIPBX*W1(II)+ DIIPBY*W2(II)  + DIIPBZ*W3(II)
     &             + RO0IPH*( GX*(CDGFBO(1,IFAC)-XXP0)
     &                      + GY*(CDGFBO(2,IFAC)-XYP0)
     &                      + GZ*(CDGFBO(3,IFAC)-XZP0))
     &             + P0IPH - PR0IPH
            ENDDO
          ELSE
            DO IFAC = 1, NFABOR
              II = IFABOR(IFAC)
              COEFU(IFAC,1) = RTPA(II,IPRIPH)
     &             + (CDGFBO(1,IFAC)-XYZCEN(1,II))*W1(II)
     &             + (CDGFBO(2,IFAC)-XYZCEN(2,II))*W2(II)
     &             + (CDGFBO(3,IFAC)-XYZCEN(3,II))*W3(II)
     &             + RO0IPH*( GX*(CDGFBO(1,IFAC)-XXP0)
     &                      + GY*(CDGFBO(2,IFAC)-XYP0)
     &                      + GZ*(CDGFBO(3,IFAC)-XZP0))
     &             + P0IPH - PR0IPH
            ENDDO
          ENDIF
C
        ENDIF


C=======================================================================
C 6.  CONVERSION EN RCODCL ICODCL
C   (SI CE DERNIER N'A PAS DEJA ETE RENSEIGNE PAR L'UTILISATEUR)
C
C      TOUT D'ABORD POUR LES VARIABLES POUR LESQUELLES IL EXISTE
C    TRAITEMENT PARTICULIER (PRESSION, VITESSE,
C=======================================================================
C
C 6.1 ENTREE
C ===========
C
C ---> La pression a un traitement Neumann, le reste Dirichlet
C                                           sera traite plus tard.
C
        IDEB = IDEBTY(IENTRE,IPHAS)
        IFIN = IFINTY(IENTRE,IPHAS)
C
        DO IVAR = 1, NVAR
          IF (IVAR.EQ.IPRIPH) THEN
            IF(IPRNEW.EQ.1) THEN
              DO II = IDEB, IFIN
                IFAC = ITRIFB(II,IPHAS)
                IF(ICODCL(IFAC,IVAR).EQ.0) THEN
                  ICODCL(IFAC,IVAR)   = 3
                  RCODCL(IFAC,IVAR,1) = 0.D0
                  RCODCL(IFAC,IVAR,2) = RINFIN
                  RCODCL(IFAC,IVAR,3) = 0.D0
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
C
C
C 6.2 SORTIE (entrée-sortie libre) (ISOLIB)
C ===================
C
C ---> La pression a un traitement Dirichlet, les vitesses 9
C        (le reste Neumann, ou Dirichlet si donnée utilisateur,
C        sera traite plus tard)
C
        IF (IPHYDR.EQ.1) THEN

C     En cas de prise en compte de la pression hydrostatique,
C     on remplit le tableau ISOSTD
C     0 -> pas une face de sortie standard (i.e. pas sortie ou sortie avec CL
C                                                de pression modifiee)
C     1 -> face de sortie libre avec CL de pression automatique.
C     le numero de la face de reference est stocke dans ISOSTD(NFABOR+1,IPHAS)
C     qui est d'abord initialise a -1 (i.e. pas de face de sortie std)
          ISOSTD(NFABOR+1,IPHAS) = -1
          DO IFAC = 1,NFABOR
            ISOSTD(IFAC,IPHAS) = 0
            IF ((ITYPFB(IFAC,IPHAS).EQ.ISOLIB).AND.
     &           (ICODCL(IFAC,IPRIPH).EQ.0)) THEN
              ISOSTD(IFAC,IPHAS) = 1
            ENDIF
          ENDDO
        ENDIF
C
C ---> Pression de recalage (unique, meme s'il y a plusieurs sorties)
C     En cas de prise en compte de la pression hydrostatique, on cherche
C     la face de reference
C
C   Determination de la pression P I' unique en parallele
C     S'il y a des faces de sortie libre, on cherche le premier
C     proc sur lequel il y en a.
C
C     On recupere aussi les coordonnees du point de reference pour les
C     sorties, pour calculer PREF ensuite
C
        IF    (ITBSLB(IPHAS).GT.0) THEN
C
          IF (IRANGP.GE.0) THEN
C
C     Indicateur de sorties locales et pointeur sur la premiere face
C       S'il y a des sorties libres standards quelque part
C         et s'il n'y en a pas en local, on affecte a IRANGD la valeur -1
            IF(IFRSLB(IPHAS).LE.0) THEN
              IRANGD = -1
C         et s'il y en a en local, on affecte a IRANGD la valeur IRANGP
            ELSE
              IRANGD = IRANGP
            ENDIF
C
C     Valeur de P I'
C         on prend le numero du dernier proc sur lequel il y en a
            CALL PARCMX(IRANGD)
C         si c'est le proc courant, on affecte P I' a COEFUP
            IF (IRANGP.EQ.IRANGD) THEN
              COEFUP = COEFU(IFRSLB(IPHAS),1)
              IF (IPHYDR.EQ.1) ISOSTD(NFABOR+1,IPHAS) = IFRSLB(IPHAS)
C         sinon on affecte 0 a COEFUP
            ELSE
              COEFUP = 0.D0
            ENDIF
C         la somme sur les procs de COEFUP donne donc P I', disponible
C            pour tous les procs (c'est plus simple qu'un bcast pour
C            lequel il faudrait au prealable que tous les proc sachent
C            quel proc envoie).
            CALL PARSOM(COEFUP)
C
C     Reference de pression pour les sorties
            IF (IRANGP.EQ.IRANGD) THEN
              XYZREF(1) = CDGFBO(1,IFRSLB(IPHAS))
              XYZREF(2) = CDGFBO(2,IFRSLB(IPHAS))
              XYZREF(3) = CDGFBO(3,IFRSLB(IPHAS))
            ELSE
              XYZREF(1) = 0.D0
              XYZREF(2) = 0.D0
              XYZREF(3) = 0.D0
            ENDIF
            ITROIS = 3
            CALL PARRSM(ITROIS,XYZREF)
C
C   Determination de la pression P I' unique en sequentiel
C     on repere la premiere face de sortie libre standard
C                       et on affecte PI' a COEFUP
C
          ELSE
C
            COEFUP = COEFU(IFRSLB(IPHAS),1)
            IF (IPHYDR.EQ.1) ISOSTD(NFABOR+1,IPHAS) = IFRSLB(IPHAS)
            XYZREF(1) = CDGFBO(1,IFRSLB(IPHAS))
            XYZREF(2) = CDGFBO(2,IFRSLB(IPHAS))
            XYZREF(3) = CDGFBO(3,IFRSLB(IPHAS))
          ENDIF
C
C     Si l'utilisateur n'a rien specifie, on met IXYZP0 a 2 pour mettre
C     a jour le point de reference
          IF (IXYZP0(IPHAS).EQ.-1) IXYZP0(IPHAS) = 2
C
        ELSEIF (IXYZP0(IPHAS).LT.0) THEN
C     S'il n'y a pas de faces de sortie, on cherche des Dirichlets
C     eventuels specifies par l'utilisateur pour y localiser
C     le point de reference.
          IFADIR = -1
          IRANGD = -1
          DO IFAC = 1, NFABOR
            IF (ICODCL(IFAC,IPRIPH).EQ.1) THEN
              IFADIR = IFAC
              IRANGD = IRANGP
            ENDIF
          ENDDO
          NFADIR = IFADIR
          IF (IRANGP.GE.0) CALL PARCMX(NFADIR)
          IF (NFADIR.GT.0) THEN
C     on met IXYZP0 a 2 pour mettre a jour le point de reference
            IXYZP0(IPHAS) = 2
C     en parallele on prend le numero du dernier proc sur lequel il y en a
            IF (IRANGP.GE.0) THEN
              CALL PARCMX(IRANGD)
              IF (IRANGP.EQ.IRANGD) THEN
                XYZREF(1) = CDGFBO(1,IFADIR)
                XYZREF(2) = CDGFBO(2,IFADIR)
                XYZREF(3) = CDGFBO(3,IFADIR)
              ELSE
                XYZREF(1) = 0.D0
                XYZREF(2) = 0.D0
                XYZREF(3) = 0.D0
              ENDIF
              ITROIS = 3
              CALL PARRSM(ITROIS,XYZREF)
            ELSE
              XYZREF(1) = CDGFBO(1,IFADIR)
              XYZREF(2) = CDGFBO(2,IFADIR)
              XYZREF(3) = CDGFBO(3,IFADIR)
            ENDIF
          ENDIF
C
        ENDIF
C
C
C   Si le point de reference n'a pas ete specifie par l'utilisateur
C   on le change et on decale alors COEFU s'il y a des sorties.
C   La pression totale dans PROPCE est aussi decalee (c'est a priori
C   inutile sauf si l'utilisateur l'utilise dans ustsns par exemple)
C
        IF (IXYZP0(IPHAS).EQ.2) THEN
          IXYZP0(IPHAS) = 1
          IIPTOT = IPPROC(IPRTOT(IPHAS))
          XXP0 = XYZREF(1) - XYZP0(1,IPHAS)
          XYP0 = XYZREF(2) - XYZP0(2,IPHAS)
          XZP0 = XYZREF(3) - XYZP0(3,IPHAS)
          XYZP0(1,IPHAS) = XYZREF(1)
          XYZP0(2,IPHAS) = XYZREF(2)
          XYZP0(3,IPHAS) = XYZREF(3)
          DO IEL = 1, NCELET
            PROPCE(IEL,IIPTOT) = PROPCE(IEL,IIPTOT)
     &           - RO0IPH*( GX*XXP0 + GY*XYP0 + GZ*XZP0 )
          ENDDO
          IF (ITBSLB(IPHAS).GT.0) THEN
            WRITE(NFECRA,8000)IPHAS,XXP0,XYP0,XZP0
            DO IFAC = 1, NFABOR
              COEFU(IFAC,1) = COEFU(IFAC,1)
     &             - RO0IPH*( GX*XXP0 + GY*XYP0 + GZ*XZP0 )
            ENDDO
            COEFUP = COEFUP - RO0IPH*( GX*XXP0 + GY*XYP0 + GZ*XZP0 )
          ELSE
            WRITE(NFECRA,8001)IPHAS,XXP0,XYP0,XZP0
          ENDIF
        ELSEIF (IXYZP0(IPHAS).EQ.-1) THEN
C     Il n'y a pas de sorties ni de Dirichlet et l'utilisateur n'a
C     rien specifie -> on met IXYZP0 a 0 pour ne plus y toucher, tout
C     en differenciant du cas =1 qui necessitera une ecriture en suite
          IXYZP0(IPHAS) = 0
        ENDIF
C
C     La pression totale doit etre recalee en Xref a la valeur
C     Po + rho_0*g.(Xref-X0)
        IF (ITBSLB(IPHAS).GT.0) THEN
          XXP0 = XYZP0(1,IPHAS)
          XYP0 = XYZP0(2,IPHAS)
          XZP0 = XYZP0(3,IPHAS)
          PREF = P0(IPHAS)
     &         + RO0IPH*( GX*(XYZREF(1)-XXP0)
     &                  + GY*(XYZREF(2)-XYP0)
     &                  + GZ*(XYZREF(3)-XZP0) )
     &         - COEFUP
        ENDIF
C
C
C ---> Entree/Sortie libre
C
        IDEB = IDEBTY(ISOLIB,IPHAS)
        IFIN = IFINTY(ISOLIB,IPHAS)
C
        DO IVAR = 1, NVAR
          IF (IVAR.EQ.IPRIPH) THEN
            IF(IPRNEW.EQ.1) THEN
              DO II = IDEB, IFIN
                IFAC = ITRIFB(II,IPHAS)
                IF(ICODCL(IFAC,IVAR).EQ.0) THEN
                  ICODCL(IFAC,IVAR)   = 1
                  RCODCL(IFAC,IVAR,1) = COEFU(IFAC,1) + PREF
                  RCODCL(IFAC,IVAR,2) = RINFIN
                  RCODCL(IFAC,IVAR,3) = 0.D0
                ENDIF
              ENDDO
            ENDIF
          ELSEIF(IVAR.EQ.IUIPH.OR.IVAR.EQ.IVIPH.OR.IVAR.EQ.IWIPH) THEN
            DO II = IDEB, IFIN
              IFAC = ITRIFB(II,IPHAS)
              IF(ICODCL(IFAC,IVAR).EQ.0) THEN
                ICODCL(IFAC,IVAR)   = 9
                RCODCL(IFAC,IVAR,1) = 0.D0
                RCODCL(IFAC,IVAR,2) = RINFIN
                RCODCL(IFAC,IVAR,3) = 0.D0
              ENDIF
            ENDDO
          ENDIF
        ENDDO
C
C
C 6.3 SYMETRIE
C =============
C
C ---> Les vecteurs et tenseurs ont un traitement particulier
C        le reste Neumann sera traite plus tard
C
        IDEB = IDEBTY(ISYMET,IPHAS)
        IFIN = IFINTY(ISYMET,IPHAS)
C
        DO IVAR = 1, NVAR
          IF ( IVAR.EQ.IUIPH.OR.IVAR.EQ.IVIPH.OR.IVAR.EQ.IWIPH.OR.
     &       ( ITYTUR(IPHAS).EQ.3.AND.
     &          (IVAR.EQ.IR11IP.OR.IVAR.EQ.IR22IP.OR.IVAR.EQ.IR33IP.OR.
     &           IVAR.EQ.IR12IP.OR.IVAR.EQ.IR13IP.OR.IVAR.EQ.IR23IP)
     &                                                          ) ) THEN
            DO II = IDEB, IFIN
              IFAC = ITRIFB(II,IPHAS)
              IF(ICODCL(IFAC,IVAR).EQ.0) THEN
                ICODCL(IFAC,IVAR)   = 4
C                RCODCL(IFAC,IVAR,1) = Modifie eventuellement par l'ALE
                RCODCL(IFAC,IVAR,2) = RINFIN
                RCODCL(IFAC,IVAR,3) = 0.D0
              ENDIF
            ENDDO
          ELSEIF(IVAR.EQ.IPRIPH) THEN
            IF(IPRNEW.EQ.1) THEN
              DO II = IDEB, IFIN
                IFAC = ITRIFB(II,IPHAS)
                IF(ICODCL(IFAC,IVAR).EQ.0) THEN
                  ICODCL(IFAC,IVAR)   = 3
                  RCODCL(IFAC,IVAR,1) = 0.D0
                  RCODCL(IFAC,IVAR,2) = RINFIN
                  RCODCL(IFAC,IVAR,3) = 0.D0
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
C
C 6.4 PAROI LISSE
C ===============
C
C ---> La vitesse et les grandeurs turbulentes ont le code 5
C        le reste Neumann sera traite plus tard
C
        IDEB = IDEBTY(IPAROI,IPHAS)
        IFIN = IFINTY(IPAROI,IPHAS)
C
        DO IVAR = 1, NVAR
          IF ( IVAR.EQ.IUIPH.OR.IVAR.EQ.IVIPH.OR.IVAR.EQ.IWIPH) THEN
            DO II = IDEB, IFIN
              IFAC = ITRIFB(II,IPHAS)
              IF(ICODCL(IFAC,IVAR).EQ.0) THEN
                ICODCL(IFAC,IVAR)   = 5
C               RCODCL(IFAC,IVAR,1) = Utilisateur
                RCODCL(IFAC,IVAR,2) = RINFIN
                RCODCL(IFAC,IVAR,3) = 0.D0
              ENDIF
            ENDDO
          ELSEIF (
     &       ( ITYTUR(IPHAS).EQ.2.AND.
     &          (IVAR.EQ.IKIPH  .OR.IVAR.EQ.IEPIPH) ).OR.
     &       ( ITYTUR(IPHAS).EQ.3.AND.
     &          (IVAR.EQ.IR11IP.OR.IVAR.EQ.IR22IP.OR.IVAR.EQ.IR33IP.OR.
     &           IVAR.EQ.IR12IP.OR.IVAR.EQ.IR13IP.OR.IVAR.EQ.IR23IP.OR.
     &           IVAR.EQ.IEPIPH)                    ).OR.
     &       ( ITURB(IPHAS).EQ.50.AND.
     &          (IVAR.EQ.IKIPH.OR.IVAR.EQ.IEPIPH.OR.IVAR.EQ.IPHIPH.OR.
     &           IVAR.EQ.IFBIPH)                    ).OR.
     &       ( ITURB(IPHAS).EQ.60.AND.
     &          (IVAR.EQ.IKIPH.OR.IVAR.EQ.IOMGIP)   ) ) THEN
            DO II = IDEB, IFIN
              IFAC = ITRIFB(II,IPHAS)
              IF(ICODCL(IFAC,IVAR).EQ.0) THEN
                ICODCL(IFAC,IVAR)   = 5
                RCODCL(IFAC,IVAR,1) = 0.D0
                RCODCL(IFAC,IVAR,2) = RINFIN
                RCODCL(IFAC,IVAR,3) = 0.D0
              ENDIF
            ENDDO
          ELSEIF(IVAR.EQ.IPRIPH) THEN
            IF(IPRNEW.EQ.1) THEN
              DO II = IDEB, IFIN
                IFAC = ITRIFB(II,IPHAS)
                IF(ICODCL(IFAC,IVAR).EQ.0) THEN
                  ICODCL(IFAC,IVAR)   = 3
                  RCODCL(IFAC,IVAR,1) = 0.D0
                  RCODCL(IFAC,IVAR,2) = RINFIN
                  RCODCL(IFAC,IVAR,3) = 0.D0
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
C
C 6.5 PAROI RUGUEUSE
C ==================
C
C ---> La vitesse et les grandeurs turbulentes ont le code 6
C      la rugosite est stockee dans rcodcl(..,..,3)
C      le reste Neumann sera traite plus tard (idem paroi lisse)
C
        IDEB = IDEBTY(IPARUG,IPHAS)
        IFIN = IFINTY(IPARUG,IPHAS)
C
        DO IVAR = 1, NVAR
          IF ( IVAR.EQ.IUIPH.OR.IVAR.EQ.IVIPH.OR.IVAR.EQ.IWIPH) THEN
            DO II = IDEB, IFIN
              IFAC = ITRIFB(II,IPHAS)
              IF(ICODCL(IFAC,IVAR).EQ.0) THEN
                ICODCL(IFAC,IVAR)   = 6
C               RCODCL(IFAC,IVAR,1) = Utilisateur
                RCODCL(IFAC,IVAR,2) = RINFIN
C               RCODCL(IFAC,IVAR,3) = Utilisateur
              ENDIF
            ENDDO
          ELSEIF (
     &       ( ITYTUR(IPHAS).EQ.2.AND.
     &          (IVAR.EQ.IKIPH  .OR.IVAR.EQ.IEPIPH) ).OR.
     &       ( ITYTUR(IPHAS).EQ.3.AND.
     &          (IVAR.EQ.IR11IP.OR.IVAR.EQ.IR22IP.OR.IVAR.EQ.IR33IP.OR.
     &           IVAR.EQ.IR12IP.OR.IVAR.EQ.IR13IP.OR.IVAR.EQ.IR23IP.OR.
     &           IVAR.EQ.IEPIPH)                    ).OR.
     &       ( ITURB(IPHAS).EQ.50.AND.
     &          (IVAR.EQ.IKIPH.OR.IVAR.EQ.IEPIPH.OR.IVAR.EQ.IPHIPH.OR.
     &           IVAR.EQ.IFBIPH)                    ).OR.
     &       ( ITURB(IPHAS).EQ.60.AND.
     &          (IVAR.EQ.IKIPH.OR.IVAR.EQ.IOMGIP)   ) ) THEN
            DO II = IDEB, IFIN
              IFAC = ITRIFB(II,IPHAS)
              IF(ICODCL(IFAC,IVAR).EQ.0) THEN
                ICODCL(IFAC,IVAR)   = 6
                RCODCL(IFAC,IVAR,1) = 0.D0
                RCODCL(IFAC,IVAR,2) = RINFIN
                RCODCL(IFAC,IVAR,3) = 0.D0
              ENDIF
            ENDDO
          ELSEIF(IVAR.EQ.IPRIPH) THEN
            IF(IPRNEW.EQ.1) THEN
              DO II = IDEB, IFIN
                IFAC = ITRIFB(II,IPHAS)
                IF(ICODCL(IFAC,IVAR).EQ.0) THEN
                  ICODCL(IFAC,IVAR)   = 3
                  RCODCL(IFAC,IVAR,1) = 0.D0
                  RCODCL(IFAC,IVAR,2) = RINFIN
                  RCODCL(IFAC,IVAR,3) = 0.D0
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
C

      ENDDO
C --- Boucle sur les phases : fin
C
C
C=======================================================================
C 6.bis  CONVERSION EN RCODCL ICODCL
C   (SI CE DERNIER N'A PAS DEJA ETE RENSEIGNE PAR L'UTILISATEUR)
C
C     MAINTENANT LES VARIABLES POUR LESQUELLES IL N'EXISTE PAS DE
C    TRAITEMENT PARTICULIER (HORS PRESSION, VITESSE ...)
C=======================================================================
C
C --- Boucle sur les phases : debut
      DO IPHAS = 1, NPHAS
C
C
C 6.1 ENTREE bis
C ===========
C
C ---> La pression a un traitement Neumann (deja traitee plus haut),
C      La vitesse  Dirichlet. Les scalaires ont un traitement
C     Dirichlet si l'utilisateur fournit une valeur, sinon on utilise
C     Neumann homogene si le flux de masse est sortant (erreur sinon).
C
        IDEB = IDEBTY(IENTRE,IPHAS)
        IFIN = IFINTY(IENTRE,IPHAS)
C
        IOK = 0
        DO IVAR = 1, NVAR
          DO II = IDEB, IFIN
            IFAC = ITRIFB(II,IPHAS)
            IF(ICODCL(IFAC,IVAR).EQ.0) THEN
C
              IF (IVAR.EQ.IUIPH.OR.IVAR.EQ.IVIPH.OR.IVAR.EQ.IWIPH)
     &             THEN
                IF (RCODCL(IFAC,IVAR,1).GT.RINFIN*0.5D0) THEN
                  IF (IOK.EQ.0) THEN
                    IOK = 1
                    CHAINE = NOMVAR(IPPRTP(IVAR))
                    WRITE(NFECRA,6040)IFAC,CHAINE(1:8)
                  ENDIF
                ELSE
                  ICODCL(IFAC,IVAR) = 1
C                 RCODCL(IFAC,IVAR,1) = Utilisateur
                  RCODCL(IFAC,IVAR,2) = RINFIN
                  RCODCL(IFAC,IVAR,3) = 0.D0
                ENDIF
C
              ELSEIF (RCODCL(IFAC,IVAR,1).GT.RINFIN*0.5D0) THEN
C
                FLUMBF = PROPFB(IFAC,IPPROB(IFLUMA(IUIPH)))
                IF( FLUMBF.GE.-EPZERO) THEN
                  ICODCL(IFAC,IVAR)   = 3
                  RCODCL(IFAC,IVAR,1) = 0.D0
                  RCODCL(IFAC,IVAR,2) = RINFIN
                  RCODCL(IFAC,IVAR,3) = 0.D0
                ELSE
                  IF (IOK.LT.2) THEN
                    IOK = IOK + 2
                    CHAINE = NOMVAR(IPPRTP(IVAR))
                    WRITE(NFECRA,6050)IFAC,CHAINE(1:8)
                  ENDIF
                ENDIF
              ELSE
                ICODCL(IFAC,IVAR) = 1
C               RCODCL(IFAC,IVAR,1) = Utilisateur
                RCODCL(IFAC,IVAR,2) = RINFIN
                RCODCL(IFAC,IVAR,3) = 0.D0
              ENDIF
C
            ENDIF
          ENDDO
        ENDDO
C
        IF (IRANGP.GE.0) CALL PARCMX(IOK)
        IF (IOK.GT.0) THEN
          IF (IOK.EQ.1 .OR. IOK.EQ.3) WRITE(NFECRA,6060)
          IF (IOK.EQ.2 .OR. IOK.EQ.3) WRITE(NFECRA,6070)
          CALL CSEXIT (1)
        ENDIF

C
C
C 6.2 SORTIE (entree sortie libre)
C ===================
C
C ---> La pression a un traitement Dirichlet, les vitesses 9 ont ete  traites plus haut.
C      Le reste Dirichlet si l'utilisateur fournit une donnee (flux de masse entrant ou sortant).
C      S'il n'y a pas de donnee utilisateur, on utilise un Neumann homogene (flux entrant et sortant)
C
C
C ---> Sortie ISOLIB
C
        IDEB = IDEBTY(ISOLIB,IPHAS)
        IFIN = IFINTY(ISOLIB,IPHAS)
C
        DO IVAR = 1, NVAR
          DO II = IDEB, IFIN
            IFAC = ITRIFB(II,IPHAS)
            IF(ICODCL(IFAC,IVAR).EQ.0) THEN
C
               IF (RCODCL(IFAC,IVAR,1).GT.RINFIN*0.5D0) THEN
                    ICODCL(IFAC,IVAR) = 3
                    RCODCL(IFAC,IVAR,1) = 0.D0
                    RCODCL(IFAC,IVAR,2) = RINFIN
                    RCODCL(IFAC,IVAR,3) = 0.D0
                ELSE
                    ICODCL(IFAC,IVAR) = 1
C                   RCODCL(IFAC,IVAR,1) = Utilisateur
                    RCODCL(IFAC,IVAR,2) = RINFIN
                    RCODCL(IFAC,IVAR,3) = 0.D0
               ENDIF
            ENDIF
          ENDDO
        ENDDO
C
C
C
C 6.3 SYMETRIE bis
C =============
C
C ---> Les vecteurs et tenseurs ont un traitement particulier
C        traite plus haut
C        le reste Neumann
C
        IDEB = IDEBTY(ISYMET,IPHAS)
        IFIN = IFINTY(ISYMET,IPHAS)
C
        DO IVAR = 1, NVAR
          DO II = IDEB, IFIN
            IFAC = ITRIFB(II,IPHAS)
            IF(ICODCL(IFAC,IVAR).EQ.0) THEN
              ICODCL(IFAC,IVAR)   = 3
              RCODCL(IFAC,IVAR,1) = 0.D0
              RCODCL(IFAC,IVAR,2) = RINFIN
              RCODCL(IFAC,IVAR,3) = 0.D0
            ENDIF
          ENDDO
        ENDDO
C
C 6.4 PAROI LISSE bis
C ===============
C
C ---> La vitesse et les grandeurs turbulentes ont le code 5
C        traite plus haut
C        le reste Neumann
C
        IDEB = IDEBTY(IPAROI,IPHAS)
        IFIN = IFINTY(IPAROI,IPHAS)
C
        DO IVAR = 1, NVAR
          DO II = IDEB, IFIN
            IFAC = ITRIFB(II,IPHAS)
            IF(ICODCL(IFAC,IVAR).EQ.0) THEN
              ICODCL(IFAC,IVAR)   = 3
              RCODCL(IFAC,IVAR,1) = 0.D0
              RCODCL(IFAC,IVAR,2) = RINFIN
              RCODCL(IFAC,IVAR,3) = 0.D0
            ENDIF
          ENDDO
        ENDDO
C
C 6.5 PAROI RUGUEUSE bis
C ==================
C
C ---> La vitesse et les grandeurs turbulentes ont le code 6
C        traite plus haut
C        le reste Neumann
C
        IDEB = IDEBTY(IPARUG,IPHAS)
        IFIN = IFINTY(IPARUG,IPHAS)
C
        DO IVAR = 1, NVAR
          DO II = IDEB, IFIN
            IFAC = ITRIFB(II,IPHAS)
            IF(ICODCL(IFAC,IVAR).EQ.0) THEN
              ICODCL(IFAC,IVAR)   = 3
              RCODCL(IFAC,IVAR,1) = 0.D0
              RCODCL(IFAC,IVAR,2) = RINFIN
              RCODCL(IFAC,IVAR,3) = 0.D0
            ENDIF
          ENDDO
        ENDDO
C
      ENDDO
C --- Boucle sur les phases : fin
C=======================================================================
C 7.  RENFORCEMENT DIAGONALE DE LA MATRICE SI AUCUN POINTS DIRICHLET
C=======================================================================
C On renforce si ISTAT=0 et si l'option est activee (IDIRCL=1)
C Si une de ces conditions est fausse, on force NDIRCL a valoir
C au moins 1 pour ne pas declaer la diagonale.
C
      DO IVAR = 1, NVAR
        NDIRCL(IVAR) = 0
        IF ( ISTAT(IVAR).GT.0 .OR. IDIRCL(IVAR).EQ.0 ) NDIRCL(IVAR) = 1
      ENDDO
C
      DO IVAR = 1, NVAR
        DO IFAC = 1, NFABOR
          IF( ICODCL(IFAC,IVAR).EQ.1 .OR. ICODCL(IFAC,IVAR).EQ.5 ) THEN
            NDIRCL(IVAR) = NDIRCL(IVAR) +1
          ENDIF
        ENDDO
        IF (IRANGP.GE.0) CALL PARCPT (NDIRCL(IVAR))
      ENDDO
C
C=======================================================================
C 8.  ON CALCULE LE FLUX DE MASSE AUX DIFFERENTS TYPES DE FACES
C       ET ON IMPRIME.
C
C     Ca serait utile de faire l'impression dans ECRLIS, mais attention,
C       on imprime le flux de masse du pas de temps precedent
C       or dans ECRLIS, on imprime a la fin du pas de temps
C       d'ou une petite incoherence possible.
C     D'autre part, ca serait utile de sortir d'autres grandeurs
C       (flux de chaleur par exemple, bilan de scalaires ...)
C
C=======================================================================
C
      IWARU = -1
      DO IPHAS = 1, NPHAS
        IUIPH  = IU(IPHAS)
        IWARU = MAX(IWARNI(IUIPH),IWARU)
      ENDDO
      IF (IRANGP.GE.0) CALL PARCMX(IWARU)
C
      IF(IWARU.GE.1 .OR. MOD(NTCABS,NTLIST).EQ.0
     &       .OR.(NTCABS.LE.NTPABS+2).OR.(NTCABS.GE.NTMABS-1)) THEN
        WRITE(NFECRA,7010)
      ENDIF
C
      DO IPHAS = 1, NPHAS
C
        IUIPH  = IU(IPHAS)
        IFLMAB = IPPROB(IFLUMA(IUIPH))
C
        IWRNP = IWARNI(IU(IPHAS))
        IF (IRANGP.GE.0) CALL PARCMX (IWRNP)
C                        ===========
C
C     On ecrit le flux de masse si IWARNI>0, a la periodicite NTLIST
C     et au deux premiers et deux derniers pas de temps.
        IF(IWRNP.GE.1 .OR. MOD(NTCABS,NTLIST).EQ.0
     &       .OR.(NTCABS.LE.NTPABS+2).OR.(NTCABS.GE.NTMABS-1)) THEN
C
          DO II = 1, NTYPMX
            FLUMTY(II) = 0.D0
          ENDDO
C
          DO II = 1, NTYPMX
            IDEB = IDEBTY(II,IPHAS)
            IFIN = IFINTY(II,IPHAS)
            DO JJ = IDEB, IFIN
              IFAC = ITRIFB(JJ,IPHAS)
              FLUMTY(II) = FLUMTY(II) + PROPFB(IFAC,IFLMAB)
            ENDDO
          ENDDO
C
C
          WRITE(NFECRA,7011) IPHAS
C
          IF ( IPPMOD(ICOMPF).LT.0 ) THEN
C
            II = IENTRE
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (INB)
              CALL PARSOM (FLUMTY(II))
            ENDIF
            WRITE(NFECRA,7020) 'Entree           ',II,INB,FLUMTY(II)
            II = IPAROI
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (INB)
              CALL PARSOM (FLUMTY(II))
            ENDIF
            WRITE(NFECRA,7020) 'Paroi lisse      ',II,INB,FLUMTY(II)
            II = IPARUG
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (INB)
              CALL PARSOM (FLUMTY(II))
            ENDIF
            WRITE(NFECRA,7020) 'Paroi rugueuse   ',II,INB,FLUMTY(II)
            II = ISYMET
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (INB)
              CALL PARSOM (FLUMTY(II))
            ENDIF
            WRITE(NFECRA,7020) 'Symetrie         ',II,INB,FLUMTY(II)

            II = ISOLIB
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (INB)
              CALL PARSOM (FLUMTY(II))
            ENDIF
            WRITE(NFECRA,7020) 'Sortie libre     ',II,INB,FLUMTY(II)


            II = IINDEF
           INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (INB)
              CALL PARSOM (FLUMTY(II))
            ENDIF
            WRITE(NFECRA,7020) 'Indefini         ',II,INB,FLUMTY(II)
C
            DO II = 1, NTYPMX
              IF( II.NE.IENTRE .AND.
     &            II.NE.IPAROI .AND.
     &            II.NE.IPARUG .AND.
     &            II.NE.ISYMET .AND.
     &            II.NE.ISOLIB .AND.
     &            II.NE.IINDEF ) THEN
                INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
                IF (IRANGP.GE.0) THEN
                  CALL PARCPT (INB)
                  CALL PARSOM (FLUMTY(II))
                ENDIF
                IF(INB.GT.0) THEN
                  WRITE(NFECRA,7020)
     &                 'Type utilisateur ',II,INB,FLUMTY(II)
                ENDIF
              ENDIF
            ENDDO
C
          ELSE
C
            II = IEQHCF
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (INB)
              CALL PARSOM (FLUMTY(II))
            ENDIF
            WRITE(NFECRA,7020) 'Entree sub. enth.',II,INB,FLUMTY(II)
C
            II = IERUCF
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (INB)
              CALL PARSOM (FLUMTY(II))
            ENDIF
            WRITE(NFECRA,7020) 'Entree subsonique',II,INB,FLUMTY(II)
C
            II = IESICF
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (INB)
              CALL PARSOM (FLUMTY(II))
            ENDIF
            WRITE(NFECRA,7020) 'Entree/Sortie imp',II,INB,FLUMTY(II)
C
            II = ISOPCF
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (INB)
              CALL PARSOM (FLUMTY(II))
            ENDIF
            WRITE(NFECRA,7020) 'Sortie subsonique',II,INB,FLUMTY(II)
C
            II = ISSPCF
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (INB)
              CALL PARSOM (FLUMTY(II))
            ENDIF
            WRITE(NFECRA,7020) 'Sortie supersoniq',II,INB,FLUMTY(II)
C
            II = IPAROI
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (INB)
              CALL PARSOM (FLUMTY(II))
            ENDIF
            WRITE(NFECRA,7020) 'Paroi            ',II,INB,FLUMTY(II)
C
            II = ISYMET
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (INB)
              CALL PARSOM (FLUMTY(II))
            ENDIF
            WRITE(NFECRA,7020) 'Symetrie         ',II,INB,FLUMTY(II)
C
            II = IINDEF
            INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
            IF (IRANGP.GE.0) THEN
              CALL PARCPT (INB)
              CALL PARSOM (FLUMTY(II))
            ENDIF
            WRITE(NFECRA,7020) 'Indefini         ',II,INB,FLUMTY(II)
C
            DO II = 1, NTYPMX
              IF( II.NE.IESICF .AND.
     &             II.NE.ISSPCF .AND.
     &             II.NE.IEQHCF .AND.
     &             II.NE.IERUCF .AND.
     &             II.NE.ISOPCF .AND.
     &             II.NE.IPAROI .AND.
     &             II.NE.ISYMET .AND.
     &             II.NE.IINDEF ) THEN
                INB = IFINTY(II,IPHAS)-IDEBTY(II,IPHAS)+1
                IF (IRANGP.GE.0) THEN
                  CALL PARCPT (INB)
                  CALL PARSOM (FLUMTY(II))
                ENDIF
                IF(INB.GT.0) THEN
                  WRITE(NFECRA,7020)
     &                 'Type utilisateur ',II,INB,FLUMTY(II)
                ENDIF
              ENDIF
            ENDDO
C
          ENDIF
C
          WRITE(NFECRA,7030)
C
        ENDIF
C
      ENDDO
C
C
C=======================================================================
C FORMATS
C=======================================================================
C
#if defined(_CS_LANG_FR)
C
 1010 FORMAT(
     &'@ @@ ATTENTION PHASE ',I3,
     &                  ', FACE ',I10,', PROP.No1 ',I5,', TYPE CL ',I4)
 1099 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,
     &'@    *********                                               ',/,
     &'@    CONDITIONS AUX LIMITES INCORRECTES OU INCOMPLETES       ',/,
     &'@                                                            ',/,
     &'@    La valeur du type ITYPFB des conditions aux limites doit',/,
     &'@      etre superieure ou egale a          1                 ',/,
     &'@        et inferieure ou egale a ',I10                       ,/,
     &'@    Une ou plusieurs erreurs sont listees ci-dessus.        ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Verifier le codage du sous-programme de definition des  ',/,
     &'@    aux limites.                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2010 FORMAT(/,' PHASE  : ',I10)
 2020 FORMAT(/,'   IFINTY : ',I10)
 2030 FORMAT(/,'   IDEBTY : ',I10)
 2040 FORMAT(/,'   ITYPFB : ',I10)
 2098 FORMAT(/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,
     &'@    *********                                               ',/,
     &'@    PROBLEME DE TRI DES FACES DE BORD                       ',/,
     &'@                                                            ',/,
     &'@    IFINTY(',I10   ,') = ',I10                               ,/,
     &'@      est superieur a                                       ',/,
     &'@    IDEBTY(',I10   ,') = ',I10                               ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Contacter l''assistance.                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2099 FORMAT(/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,
     &'@    *********                                               ',/,
     &'@    PROBLEME DE TRI DES FACES DE BORD SUR UN RANG DISTANT   ',/,
     &'@                                                            ',/,
     &'@    IFINTY(',I10   ,')                                      ',/,
     &'@      est superieur a                                       ',/,
     &'@    IDEBTY(',I10   ,')                                      ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Contacter l''assistance.                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3099 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,
     &'@    *********                                               ',/,
     &'@    PROBLEME DE TRI DES FACES DE BORD                       ',/,
     &'@                                                            ',/,
     &'@    Phase          ',I10                                     ,/,
     &'@      nombre de faces classees par type = ',I10              ,/,
     &'@      nombre de faces de bord  NFABOR   = ',I10              ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Contacter l''assistance.                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
 6010 FORMAT ( /,/,
     & '   ** INFORMATIONS SUR LE TYPE DE FACES DE BORD',/,
     & '      -----------------------------------------',/)
 6011 FORMAT (
     &'   Phase : ',I4,                                               /,
     &'---------------------------------------------------------------',
     &'----------',
     &                                                                /,
     &'Type de bord           Code    Nb faces (tot/loc)  ',
     &' Index min   Index max',
     &                                                                /,
     &'---------------------------------------------------------------',
     &'----------')
 6020 FORMAT (
     & A17,I10,I10,3I12)
 6030 FORMAT(
     &'---------------------------------------------------------------',
     &'----------'/)
C
 6040 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,
     &'@    *********                                               ',/,
     &'@    CONDITIONS AUX LIMITES INCORRECTES OU INCOMPLETES       ',/,
     &'@                                                            ',/,
     &'@    Face de bord ',I10                                       ,/,
     &'@    Face declaree en entree (ou sortie) a vitesse imposee,  ',/,
     &'@      mais la valeur de la vitesse n''a pas ete fournie pour',/,
     &'@      la composante ',A8                                     ,/,
     &'@    (ce message ne s''affiche que pour la premiere face ou  ',/,
     &'@      le probleme est detecte, mais d''autres faces peuvent ',/,
     &'@      etre concernees)                                      ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Verifier les conditions aux limites dans l''Interface   ',/,
     &'@    ou dans le sous-programme utilisateur correspondant.    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 6050 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,
     &'@    *********                                               ',/,
     &'@    CONDITIONS AUX LIMITES INCORRECTES OU INCOMPLETES       ',/,
     &'@                                                            ',/,
     &'@    Face de bord ',I10                                       ,/,
     &'@    Face declaree en entree (ou sortie) a vitesse imposee   ',/,
     &'@      avec un flux rentrant.                                ',/,
     &'@    La valeur de la variable ',A8  ,' n''a pas ete specifiee',/,
     &'@      (condition de Dirichlet).                             ',/,
     &'@    (ce message ne s''affiche que pour la premiere face ou  ',/,
     &'@      le probleme est detecte, mais d''autres faces peuvent ',/,
     &'@      etre concernees)                                      ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute                           ',/,
     &'@                                                            ',/,
     &'@    Verifier les conditions aux limites dans l''Interface   ',/,
     &'@    ou dans le sous-programme utilisateur correspondant.    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 6060 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,
     &'@    *********                                               ',/,
     &'@    CONDITIONS AUX LIMITES INCORRECTES OU INCOMPLETES       ',/,
     &'@                                                            ',/,
     &'@    Au moins une face de bord declaree en entree            ',/,
     &'@      (ou sortie) a vitesse imposee pour laquelle la valeur ',/,
     &'@      de la vitesse n''a pas ete fournie pour toutes les    ',/,
     &'@      composantes.                                          ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Verifier les conditions aux limites dans l''Interface   ',/,
     &'@    ou dans le sous-programme utilisateur correspondant.    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 6070 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE LA VERIFICATION DES COND. LIM.',/,
     &'@    *********                                               ',/,
     &'@    CONDITIONS AUX LIMITES INCORRECTES OU INCOMPLETES       ',/,
     &'@                                                            ',/,
     &'@    Au moins une face de bord declaree en entree            ',/,
     &'@      (ou sortie) a vitesse imposee avec un flux rentrant   ',/,
     &'@      pour laquelle la valeur d''une variable n''a pas ete  ',/,
     &'@      specifiee (condition de Dirichlet).                   ',/,
     &'@    Le calcul ne sera pas execute                           ',/,
     &'@                                                            ',/,
     &'@    Verifier les conditions aux limites dans l''Interface   ',/,
     &'@    ou dans le sous-programme utilisateur correspondant.    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
 7010 FORMAT ( /,/,
     & '   ** INFORMATIONS SUR LE FLUX DE MASSE AU BORD',/,
     & '      -----------------------------------------',/)
 7011 FORMAT (
     &'   Phase : ',I4,                                               /,
     &'---------------------------------------------------------------',
     &                                                                /,
     &'Type de bord           Code    Nb faces           Flux de masse',
     &                                                                /,
     &'---------------------------------------------------------------')
 7020 FORMAT (
     & A17,I10,I12,6X,E18.9)
 7030 FORMAT(
     &'---------------------------------------------------------------',
     &                                                                /)
C
 8000 FORMAT(/,
     &'PHASE ',I4,' :                                               ',/,
     &'Faces de bord d''entree/sortie libre detectees               ',/,
     &'Mise a jour du point de reference pour la pression totale    ',/,
     &' XYZP0 = ',E14.5,E14.5,E14.5                  ,/)
 8001 FORMAT(/,
     &'PHASE ',I4,' :                                               ',/,
     &'Faces de bord a Dirichlet de pression impose detectees       ',/,
     &'Mise a jour du point de reference pour la pression totale    ',/,
     &' XYZP0 = ',E14.5,E14.5,E14.5                  ,/)
C
C-----------------------------------------------------------------------
C
#else
C
 1010 FORMAT(
     &'@ @@ WARNING  PHASE ',I3,
     &                  ', FACE ',I10,', PROP.No1 ',I5,', BC TYPE ',I4)
 1099 FORMAT(
     &'@'                                                            ,/,
     &'@'                                                            ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@'                                                            ,/,
     &'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,
     &'@    ********'                                                ,/,
     &'@    INCORRECT OR INCOMPLETE BOUNDARY CONDITIONS'             ,/,
     &'@'                                                            ,/,
     &'@    Value of type ITYPFB for boundary conditions must be'    ,/,
     &'@          greater or equal to          1'                    ,/,
     &'@      and less than or equal to ', I10                       ,/,
     &'@    One or more errors are listed above.'                    ,/,
     &'@'                                                            ,/,
     &'@    The calculation will not be run.'                        ,/,
     &'@'                                                            ,/,
     &'@    Verify the boundary condition definitions in the'        ,/,
     &'@    appropriate user subroutine.'                            ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@'                                                            ,/)
C
 2010 FORMAT(/,' PHASE  : ',I10)
 2020 FORMAT(/,'   IFINTY : ',I10)
 2030 FORMAT(/,'   IDEBTY : ',I10)
 2040 FORMAT(/,'   ITYPFB : ',I10)
 2098 FORMAT(/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@'                                                            ,/,
     &'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,
     &'@    ********'                                                ,/,
     &'@    PROBLEM WITH ORDERING OF BOUNDARY FACES'                 ,/,
     &'@'                                                            ,/,
     &'@    IFINTY(',I10   ,') = ',I10                               ,/,
     &'@      is greater than'                                       ,/,
     &'@    IDEBTY(',I10   ,') = ',I10                               ,/,
     &'@'                                                            ,/,
     &'@    The calculation will not be run.'                        ,/,
     &'@'                                                            ,/,
     &'@    Contact support.'                                        ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2099 FORMAT(/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@'                                                            ,/,
     &'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,
     &'@    ********'                                                ,/,
     &'@    PROBLEM WITH ORDERING OF BOUNDARY FACES'                 ,/,
     &'@    ON A DISTANT RANK.'                                      ,/,
     &'@'                                                            ,/,
     &'@    IFINTY(',I10   ,')                                      ',/,
     &'@      is greater than'                                       ,/,
     &'@    IDEBTY(',I10   ,')                                      ',/,
     &'@'                                                            ,/,
     &'@    The calculation will not be run.'                        ,/,
     &'@'                                                            ,/,
     &'@    Contact support.'                                        ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3099 FORMAT(
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@'                                                            ,/,
     &'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,
     &'@    ********'                                                ,/,
     &'@    PROBLEM WITH ORDERING OF BOUNDARY FACES'                 ,/,
     &'@'                                                            ,/,
     &'@    Phase          ',I10                                     ,/,
     &'@      number of faces classified by type = ',I10             ,/,
     &'@      number of boundary faces (NFABOR)  = ',I10             ,/,
     &'@'                                                            ,/,
     &'@    The calculation will not be run.'                        ,/,
     &'@'                                                            ,/,
     &'@    Contact support.'                                        ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
 6010 FORMAT ( /,/,
     & '   ** INFORMATION ON BOUNDARY FACES TYPE',/,
     & '      ----------------------------------',/)
 6011 FORMAT (
     &'   Phase : ',I4,                                               /,
     &'---------------------------------------------------------------',
     &'----------',
     &                                                                /,
     &'Boundary type          Code    Nb faces (tot/loc)  ',
     &' Min index   Max index',
     &                                                                /,
     &'---------------------------------------------------------------',
     &'----------')
 6020 FORMAT (
     & A17,I10,I10,3I12)
 6030 FORMAT(
     &'---------------------------------------------------------------',
     &'----------'/)
C
 6040 FORMAT(
     &'@'                                                            ,/,
     &'@'                                                            ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@'                                                            ,/,
     &'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,
     &'@    ********'                                                ,/,
     &'@    INCORRECT OR INCOMPLETE BOUNDARY CONDITIONS'             ,/,
     &'@'                                                            ,/,
     &'@    Boundary face ',I10                                      ,/,
     &'@    Face declared as inlet (or outlet) with prescribed'      ,/,
     &'@      velocity, but velocity value has not been given for'   ,/,
     &'@      component ',A8                                         ,/,
     &'@    (this message is printed only for the first face on'     ,/,
     &'@      which this problem is detected, but other faces may'   ,/,
     &'@      be concerned)'                                         ,/,
     &'@'                                                            ,/,
     &'@    The calculation will not be run.                        ',/,
     &'@'                                                            ,/,
     &'@    Verify the boundary condition definitions in the GUI'    ,/,
     &'@    or in the appropriate user subroutine.'                  ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 6050 FORMAT(
     &'@'                                                            ,/,
     &'@'                                                            ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@'                                                            ,/,
     &'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,
     &'@    ********'                                                ,/,
     &'@    INCORRECT OR INCOMPLETE BOUNDARY CONDITIONS'             ,/,
     &'@'                                                            ,/,
     &'@    Boundary face ',I10                                      ,/,
     &'@    Face declared as inlet (or outlet) with prescribed'      ,/,
     &'@      velocity with entering flow.'                          ,/,
     &'@    Value of variable ',A8  ,' has not been specified'       ,/,
     &'@      (Dirichlet condition).'                                ,/,
     &'@    (this message is printed only for the first face on'     ,/,
     &'@      which this problem is detected, but other faces may'   ,/,
     &'@      be concerned)'                                         ,/,
     &'@'                                                            ,/,
     &'@    The calculation will not be run.                        ',/,
     &'@'                                                            ,/,
     &'@    Verify the boundary condition definitions in the GUI'    ,/,
     &'@    or in the appropriate user subroutine.'                  ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 6060 FORMAT(
     &'@'                                                            ,/,
     &'@'                                                            ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@'                                                            ,/,
     &'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,
     &'@    ********'                                                ,/,
     &'@    INCORRECT OR INCOMPLETE BOUNDARY CONDITIONS'             ,/,
     &'@'                                                            ,/,
     &'@    At least one boundary face declared as inlet (or'        ,/,
     &'@      outlet) with prescribed velocity for which the'        ,/,
     &'@      velocity value has not been assigned for all'          ,/,
     &'@      components.'                                           ,/,
     &'@    The calculation will not be run.                        ',/,
     &'@'                                                            ,/,
     &'@    Verify the boundary condition definitions in the GUI'    ,/,
     &'@    or in the appropriate user subroutine.'                  ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 6070 FORMAT(
     &'@'                                                            ,/,
     &'@'                                                            ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@'                                                            ,/,
     &'@ @@ WARNING: ABORT BY BOUNDARY CONDITION CHECK'              ,/,
     &'@    ********'                                                ,/,
     &'@    INCORRECT OR INCOMPLETE BOUNDARY CONDITIONS'             ,/,
     &'@'                                                            ,/,
     &'@    At least one boundary face declared as inlet (or'        ,/,
     &'@      outlet) with prescribed velocity with an entering'     ,/,
     &'@      flow for which the value of a variable has not been'   ,/,
     &'@      specified (Dirichlet condition).'                      ,/,
     &'@    The calculation will not be run.                        ',/,
     &'@'                                                            ,/,
     &'@    Verify the boundary condition definitions in the GUI'    ,/,
     &'@    or in the appropriate user subroutine.'                  ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
 7010 FORMAT ( /,/,
     & '   ** BOUNDARY MASS FLOW INFORMATION',/,
     & '      ------------------------------',/)
 7011 FORMAT (
     &'   Phase : ',I4,                                               /,
     &'---------------------------------------------------------------',
     &                                                                /,
     &'Boundary type          Code    Nb faces           Mass flow'   ,
     &                                                                /,
     &'---------------------------------------------------------------')
 7020 FORMAT (
     & A17,I10,I12,6X,E18.9)
 7030 FORMAT(
     &'---------------------------------------------------------------',
     &                                                                /)
C
 8000 FORMAT(/,
     &'PHASE ',I4,' :'                                               ,/,
     &'Boundary faces with free inlet/outlet detected'               ,/,
     &'Update of reference point for total pressure'                 ,/,
     &' XYZP0 = ',E14.5,E14.5,E14.5                  ,/)
 8001 FORMAT(/,
     &'PHASE ',I4,' :'                                               ,/,
     &'Boundary faces with pressure Dirichlet condition detected'    ,/,
     &'Update of reference point for total pressure'                 ,/,
     &' XYZP0 = ',E14.5,E14.5,E14.5                  ,/)
C
#endif
C
C
      RETURN
      END
c@z
