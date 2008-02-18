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
                        SUBROUTINE ALTYCL
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , IALTYB , ICODCL , IMPALE ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   RCODCL , XYZNO0 , DEPALE ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC TRAITEMENT DES CODES DE CONDITIONS AUX LIMITES IALTYB POUR L'ALE
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
CARGU !   (NFAC+1)   !    !     !                                      !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !                                      !
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS)      !    !     !                                      !
CARGU ! IALTYB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD POUR L'ALE    !
CARGU ! ICODCL       ! TE ! <-> ! CODE DE CONDITION LIMITES AUX FACES  !
CARGU !  (NFABOR,NVAR!    !     !  DE BORD                             !
CARGU !              !    !     ! = 1   -> DIRICHLET                   !
CARGU !              !    !     ! = 3   -> DENSITE DE FLUX             !
CARGU ! IMPALE(NNOD) ! TE !  -> ! INDICATEUR DE DELACEMENT IMPOSE      !
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
CARGU ! RCODCL       ! TR ! <-> ! VALEUR DES CONDITIONS AUX LIMITES    !
CARGU !  (NFABOR,NVAR!    !     !  AUX FACES DE BORD                   !
CARGU !              !    !     ! RCODCL(1) = VALEUR DU DIRICHLET      !
CARGU !              !    !     ! RCODCL(2) = VALEUR DU COEF. D'ECHANGE!
CARGU !              !    !     !  EXT. (INFINIE SI PAS D'ECHANGE)     !
CARGU !              !    !     ! RCODCL(3) = VALEUR DE LA DENSITE DE  !
CARGU !              !    !     !  FLUX (NEGATIF SI GAIN) W/m2         !
CARGU !              !    !     ! POUR LES VITESSES (VISTL+VISCT)*GRADU!
CARGU !              !    !     ! POUR LA PRESSION             DT*GRADP!
CARGU !              !    !     ! POUR LES SCALAIRES                   !
CARGU !              !    !     !        CP*(VISCLS+VISCT/SIGMAS)*GRADT!
CARGU ! DEPALE(NNOD,3! TR !  -> ! DEPLACEMENT AUX NOEUDS               !
CARGU ! XYZNO0(3,NNOD! TR !  -> ! COORDONNEES NOEUDS MAILLAGE INITIAL  !
CARGU ! W1,2,3,4,5,6 ! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU !  (NCELET     !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
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
      INTEGER          ITYPFB(NFABOR,NPHAS)
      INTEGER          IALTYB(NFABOR), ICODCL(NFABOR,NVAR)
      INTEGER          IMPALE(NNOD)
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
      DOUBLE PRECISION RCODCL(NFABOR,NVAR,3)
      DOUBLE PRECISION DEPALE(NNOD,3), XYZNO0(3,NNOD)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION W4(NCELET),W5(NCELET),W6(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC, IEL, IPHAS, IUIPH, IVIPH, IWIPH
      INTEGER          II, INOD, IECRW, ICPT, IERROR
      DOUBLE PRECISION DDEPX, DDEPY, DDEPZ
      DOUBLE PRECISION SRFBNF, RNX, RNY, RNZ
      DOUBLE PRECISION RCODCX, RCODCY, RCODCZ, RCODSN
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
C Mise a zero des RCODCL non specifies
      DO IFAC = 1, NFABOR
        IF (RCODCL(IFAC,IUMA,1).GT.RINFIN*0.5D0)
     &       RCODCL(IFAC,IUMA,1) = 0.D0
        IF (RCODCL(IFAC,IVMA,1).GT.RINFIN*0.5D0)
     &       RCODCL(IFAC,IVMA,1) = 0.D0
        IF (RCODCL(IFAC,IWMA,1).GT.RINFIN*0.5D0)
     &       RCODCL(IFAC,IWMA,1) = 0.D0
      ENDDO
C
C=======================================================================
C 2.  VERIFICATION DE LA CONSISTANCE DES IALTYB DONNES DANS USALCL
C=======================================================================
C  (valeur 0 autorisee)
C
      IERROR = 0
      DO IFAC = 1, NFABOR
        IF (IALTYB(IFAC).NE.0      .AND.
     &      IALTYB(IFAC).NE.IBFIXE .AND.
     &      IALTYB(IFAC).NE.IGLISS .AND.
     &      IALTYB(IFAC).NE.IVIMPO ) THEN
          WRITE(NFECRA,1000)IFAC,IPRFML(IFMFBR(IFAC),1),
     &         IALTYB(IFAC)
          IERROR = IERROR + 1
        ENDIF
      ENDDO
C
C
C=======================================================================
C 3.  CONVERSION EN RCODCL ICODCL
C=======================================================================
C
C Calcul ou ecrasement de RCODCL a partir du deplacement impose,
C   si tous les noeuds d'une face sont a deplacement impose
C
      DO IFAC = 1, NFABOR
        IECRW = 0
        DDEPX = 0.D0
        DDEPY = 0.D0
        DDEPZ = 0.D0
        ICPT  = 0
        DO II = IPNFBR(IFAC), IPNFBR(IFAC+1)-1
          INOD = NODFBR(II)
          IF (IMPALE(INOD).EQ.0) IECRW = IECRW + 1
          ICPT = ICPT + 1
          DDEPX = DDEPX + DEPALE(INOD,1)+XYZNO0(1,INOD)-XYZNOD(1,INOD)
          DDEPY = DDEPY + DEPALE(INOD,2)+XYZNO0(2,INOD)-XYZNOD(2,INOD)
          DDEPZ = DDEPZ + DEPALE(INOD,3)+XYZNO0(3,INOD)-XYZNOD(3,INOD)
        ENDDO
        IF (IECRW.EQ.0) THEN
          IEL = IFABOR(IFAC)
          IALTYB(IFAC) = IVIMPO
          RCODCL(IFAC,IUMA,1) = DDEPX/DT(IEL)/ICPT
          RCODCL(IFAC,IVMA,1) = DDEPY/DT(IEL)/ICPT
          RCODCL(IFAC,IWMA,1) = DDEPZ/DT(IEL)/ICPT
        ENDIF
      ENDDO
C
C Remplissage des autres RCODCL a partir des ITYALB
C
      DO IFAC = 1, NFABOR
C
        IEL = IFABOR(IFAC)
C
C --> Faces fixes
C     On force alors les noeuds en question a etre fixes, pour eviter
C       des problemes eventuels aux coins
C
        IF (IALTYB(IFAC).EQ.IBFIXE) THEN
          ICPT = 0
          IF (ICODCL(IFAC,IUMA).EQ.0) THEN
            ICPT = ICPT + 1
            ICODCL(IFAC,IUMA) = 1
            RCODCL(IFAC,IUMA,1) = 0.D0
          ENDIF
          IF (ICODCL(IFAC,IVMA).EQ.0) THEN
            ICPT = ICPT + 1
            ICODCL(IFAC,IVMA) = 1
            RCODCL(IFAC,IVMA,1) = 0.D0
          ENDIF
          IF (ICODCL(IFAC,IWMA).EQ.0) THEN
            ICPT = ICPT + 1
            ICODCL(IFAC,IWMA) = 1
            RCODCL(IFAC,IWMA,1) = 0.D0
          ENDIF
C      Si on a fixe les trois composantes, alors on fixe les noeuds
C        correspondants. Sinon c'est que l'utilisateur a modifie quelque
C        chose ... on le laisse seul maitre
          IF (ICPT.EQ.3) THEN
            DO II = IPNFBR(IFAC), IPNFBR(IFAC+1)-1
              INOD = NODFBR(II)
              IF (IMPALE(INOD).EQ.0) THEN
                DEPALE(INOD,1) = 0.D0
                DEPALE(INOD,2) = 0.D0
                DEPALE(INOD,3) = 0.D0
                IMPALE(INOD) = 1
              ENDIF
            ENDDO
          ENDIF
C
C
C --> Faces de glissement
        ELSEIF (IALTYB(IFAC).EQ.IGLISS) THEN
C
          IF (ICODCL(IFAC,IUMA).EQ.0) ICODCL(IFAC,IUMA) = 4
          IF (ICODCL(IFAC,IVMA).EQ.0) ICODCL(IFAC,IVMA) = 4
          IF (ICODCL(IFAC,IWMA).EQ.0) ICODCL(IFAC,IWMA) = 4
C
C --> Faces a vitesse imposee
        ELSEIF (IALTYB(IFAC).EQ.IVIMPO) THEN
C
          IF (ICODCL(IFAC,IUMA).EQ.0) ICODCL(IFAC,IUMA) = 1
          IF (ICODCL(IFAC,IVMA).EQ.0) ICODCL(IFAC,IVMA) = 1
          IF (ICODCL(IFAC,IWMA).EQ.0) ICODCL(IFAC,IWMA) = 1
C
        ENDIF
C
      ENDDO
C
C=======================================================================
C 4.  VERIFICATION DE LA CONSISTANCE DES ICODCL
C=======================================================================
C
C     IERROR a ete initialise plus haut
      DO IFAC = 1, NFABOR
C
        IF (ICODCL(IFAC,IUMA).NE.1 .AND.
     &      ICODCL(IFAC,IUMA).NE.3 .AND.
     &      ICODCL(IFAC,IUMA).NE.4 ) THEN
          WRITE(NFECRA,2000) IFAC,IPRFML(IFMFBR(IFAC),1),
     &         ICODCL(IFAC,IUMA)
          IERROR = IERROR + 1
        ENDIF
        IF (ICODCL(IFAC,IVMA).NE.1 .AND.
     &      ICODCL(IFAC,IVMA).NE.3 .AND.
     &      ICODCL(IFAC,IVMA).NE.4 ) THEN
          WRITE(NFECRA,2001) IFAC,IPRFML(IFMFBR(IFAC),1),
     &         ICODCL(IFAC,IVMA)
          IERROR = IERROR + 1
        ENDIF
        IF (ICODCL(IFAC,IWMA).NE.1 .AND.
     &      ICODCL(IFAC,IWMA).NE.3 .AND.
     &      ICODCL(IFAC,IWMA).NE.4 ) THEN
          WRITE(NFECRA,2002) IFAC,IPRFML(IFMFBR(IFAC),1),
     &         ICODCL(IFAC,IWMA)
          IERROR = IERROR + 1
        ENDIF
C
        IF ( ( ICODCL(IFAC,IUMA).EQ.4 .OR.
     &         ICODCL(IFAC,IVMA).EQ.4 .OR.
     &         ICODCL(IFAC,IWMA).EQ.4    ) .AND.
     &       ( ICODCL(IFAC,IUMA).NE.4 .OR.
     &         ICODCL(IFAC,IVMA).NE.4 .OR.
     &         ICODCL(IFAC,IWMA).NE.4    ) ) THEN
          WRITE(NFECRA,3000) IFAC,IPRFML(IFMFBR(IFAC),1),
     &         ICODCL(IFAC,IUMA),ICODCL(IFAC,IVMA),ICODCL(IFAC,IWMA)
          IERROR = IERROR + 1
        ENDIF
C
      ENDDO
C
      IF (IERROR.GT.0) THEN
        WRITE(NFECRA,4000)
        CALL CSEXIT(1)
C       ===========
      ENDIF
C
C=======================================================================
C 5.  VITESSE DE DEFILEMENT POUR LES PAROIS FLUIDES ET SYMETRIES
C=======================================================================
C
C Pour les symetries on rajoute toujours la vitesse de maillage, car on
C   ne conserve que la vitesse normale
C Pour les parois, on prend la vitesse de maillage si l'utilisateur n'a
C   pas specifie RCODCL, sinon on laisse RCODCL pour la vitesse tangente
C   et on prend la vitesse de maillage pour la composante normale.
C On se base uniquement sur ITYPFB, a l'utilisateur de gere les choses
C   s'il rentre en CL non standards.
      IPHAS = 1
      IUIPH = IU(IPHAS)
      IVIPH = IV(IPHAS)
      IWIPH = IW(IPHAS)
C
      DO IFAC = 1, NFABOR
C
        IF (IALTYB(IFAC).EQ.IVIMPO) THEN
C
          IF ( ITYPFB(IFAC,IPHAS).EQ.ISYMET ) THEN
            RCODCL(IFAC,IUIPH,1) = RCODCL(IFAC,IUMA,1)
            RCODCL(IFAC,IVIPH,1) = RCODCL(IFAC,IVMA,1)
            RCODCL(IFAC,IWIPH,1) = RCODCL(IFAC,IWMA,1)
          ENDIF
C
          IF ( ITYPFB(IFAC,IPHAS).EQ.IPAROI .OR.
     &         ITYPFB(IFAC,IPHAS).EQ.IPARUG ) THEN
C Si une des composantes de vitesse de glissement a ete
C    modifiee par l'utilisateur, on ne fixe que la vitesse
C    normale
            IF (RCODCL(IFAC,IUIPH,1).GT.RINFIN*0.5D0 .AND.
     &           RCODCL(IFAC,IVIPH,1).GT.RINFIN*0.5D0 .AND.
     &           RCODCL(IFAC,IWIPH,1).GT.RINFIN*0.5D0) THEN
              RCODCL(IFAC,IUIPH,1) = RCODCL(IFAC,IUMA,1)
              RCODCL(IFAC,IVIPH,1) = RCODCL(IFAC,IVMA,1)
              RCODCL(IFAC,IWIPH,1) = RCODCL(IFAC,IWMA,1)
            ELSE
C On met a 0 les composantes de RCODCL non specifiees
              IF (RCODCL(IFAC,IUIPH,1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IUIPH,1) = 0.D0
              IF (RCODCL(IFAC,IVIPH,1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IVIPH,1) = 0.D0
              IF (RCODCL(IFAC,IWIPH,1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IWIPH,1) = 0.D0
C
              SRFBNF = RA(ISRFBN-1+IFAC)
              RNX = SURFBO(1,IFAC)/SRFBNF
              RNY = SURFBO(2,IFAC)/SRFBNF
              RNZ = SURFBO(3,IFAC)/SRFBNF
              RCODCX = RCODCL(IFAC,IUIPH,1)
              RCODCY = RCODCL(IFAC,IVIPH,1)
              RCODCZ = RCODCL(IFAC,IWIPH,1)
              RCODSN = (RCODCL(IFAC,IUMA,1)-RCODCX)*RNX
     &             +   (RCODCL(IFAC,IVMA,1)-RCODCY)*RNY
     &             +   (RCODCL(IFAC,IVMA,1)-RCODCZ)*RNZ
              RCODCL(IFAC,IUIPH,1) = RCODCX + RCODSN*RNX
              RCODCL(IFAC,IVIPH,1) = RCODCY + RCODSN*RNY
              RCODCL(IFAC,IWIPH,1) = RCODCZ + RCODSN*RNZ
            ENDIF
          ENDIF
C
        ENDIF
C
      ENDDO
C
C=======================================================================
C FORMATS
C=======================================================================
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@ METHODE ALE                                                ',/,
     &'@ TYPE DE BORD NON RECONNU                                   ',/,
     &'@   FACE ',I10   ,' ; PROPRIETE 1 :',I10                      ,/,
     &'@     ALTYB FACE : ', I10                                     ,/,
     &'@                                                            '  )
 2000 FORMAT(
     &'@                                                            ',/,
     &'@ METHODE ALE                                                ',/,
     &'@ TYPE DE CONDITION A LA LIMITE NON RECONNU POUR LA VITESSE  ',/,
     &'@   DE MAILLAGE SELON X (IUMA)                               ',/,
     &'@   FACE ',I10   ,' ; PROPRIETE 1 :',I10                      ,/,
     &'@     ICODCL : ', I10                                         ,/,
     &'@                                                            ',/,
     &'@ Les seules valeurs autorisees pour ICODCL sont             ',/,
     &'@   1 : Dirichlet                                            ',/,
     &'@   3 : Neumann                                              ',/,
     &'@   4 : Glissement                                           ',/,
     &'@                                                            '  )
 2001 FORMAT(
     &'@                                                            ',/,
     &'@ METHODE ALE                                                ',/,
     &'@ TYPE DE CONDITION A LA LIMITE NON RECONNU POUR LA VITESSE  ',/,
     &'@   DE MAILLAGE SELON Y (IVMA)                               ',/,
     &'@   FACE ',I10   ,' ; PROPRIETE 1 :',I10                      ,/,
     &'@     ICODCL : ', I10                                         ,/,
     &'@                                                            ',/,
     &'@ Les seules valeurs autorisees pour ICODCL sont             ',/,
     &'@   1 : Dirichlet                                            ',/,
     &'@   3 : Neumann                                              ',/,
     &'@   4 : Glissement                                           ',/,
     &'@                                                            '  )
 2002 FORMAT(
     &'@                                                            ',/,
     &'@ METHODE ALE                                                ',/,
     &'@ TYPE DE CONDITION A LA LIMITE NON RECONNU POUR LA VITESSE  ',/,
     &'@   DE MAILLAGE SELON Z (IWMA)                               ',/,
     &'@   FACE ',I10   ,' ; PROPRIETE 1 :',I10                      ,/,
     &'@     ICODCL : ', I10                                         ,/,
     &'@                                                            ',/,
     &'@ Les seules valeurs autorisees pour ICODCL sont             ',/,
     &'@   1 : Dirichlet                                            ',/,
     &'@   3 : Neumann                                              ',/,
     &'@   4 : Glissement                                           ',/,
     &'@                                                            '  )
 3000 FORMAT(
     &'@                                                            ',/,
     &'@ METHODE ALE                                                ',/,
     &'@ INCOHERENCE DANS LES TYPES DE CONDITIONS A LA LIMITE       ',/,
     &'@   POUR LA VITESSE DE MAILLAGE                              ',/,
     &'@   FACE ',I10   ,' ; PROPRIETE 1 :',I10                      ,/,
     &'@     ICODCL(.,IUMA) : ', I10                                 ,/,
     &'@     ICODCL(.,IVMA) : ', I10                                 ,/,
     &'@     ICODCL(.,IWMA) : ', I10                                 ,/,
     &'@                                                            ',/,
     &'@ Si une composante est traitee en glissement (ICODCL=4),    ',/,
     &'@   toutes les composantes doivent etre traitees en          ',/,
     &'@   glissement.                                              ',/,
     &'@                                                            '  )
 4000 FORMAT(
     &'@                                                            ',/,
     &'@ METHODE ALE                                                ',/,
     &'@ INCOHERENCE DANS LES CONDITIONS AUX LIMITES POUR LA VITESSE',/,
     &'@   DE MAILLAGE                                              ',/,
     &'@   (cf. message(s) ci-dessus)                               ',/,
     &'@                                                            ',/,
     &'@ Verifier usalcl.F                                          ',/,
     &'@                                                            '  )
C
C
      RETURN
      END
c@z
