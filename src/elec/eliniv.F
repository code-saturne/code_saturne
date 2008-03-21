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
                        SUBROUTINE ELINIV
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB , COEFA  , COEFB  ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  INITIALISATION DES VARIABLES DE CALCUL
CFONC     POUR LA PHYSIQUE PARTICULIERE : VERSIONS ELECTRIQUES
CFONC     PENDANT DE USINIV.F
CFONC
CFONC  Cette routine est appelee en debut de calcul (suite ou non)
CFONC      avant le debut de la boucle en temps
CFONC
CFONC  Elle permet d'INITIALISER ou de MODIFIER (pour les calculs suite)
CFONC      les variables de calcul,
CFONC      les valeurs du pas de temps
CFONC
CFONC
CFONC  On dispose ici de ROM et VISCL initialises par RO0 et VISCL0
CFONC      ou relues d'un fichier suite
CFONC  On ne dispose des variables VISCLS, CP (quand elles sont
CFONC      definies) que si elles ont pu etre relues dans un fichier
CFONC      suite de calcul
CFONC
CFONC  Les proprietes physiaues sont accessibles dans le tableau
CFONC      PROPCE (prop au centre), PROPFA (aux faces internes),
CFONC      PROPFB (prop aux faces de bord)
CFONC      Ainsi,
CFONC       PROPCE(IEL,IPPROC(IROM  (IPHAS))) designe ROM   (IEL ,IPHAS)
CFONC       PROPCE(IEL,IPPROC(IVISCL(IPHAS))) designe VISCL (IEL ,IPHAS)
CFONC       PROPCE(IEL,IPPROC(ICP   (IPHAS))) designe CP    (IEL ,IPHAS)
CFONC       PROPCE(IEL,IPPROC(IVISLS(ISCAL))) designe VISLS (IEL ,ISCAL)
CFONC
CFONC       PROPFA(IFAC,IPPROF(IFLUMA(IVAR ))) designe FLUMAS(IFAC,IVAR)
CFONC
CFONC       PROPFB(IFAC,IPPROB(IROM  (IPHAS))) designe ROMB  (IFAC,IPHAS)
CFONC       PROPFB(IFAC,IPPROB(IFLUMA(IVAR ))) designe FLUMAB(IFAC,IVAR)
CFONC
CFONC  LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
CFONC      SE FERA EN STANDARD DANS LE SOUS PROGRAMME PPPHYV
CFONC      ET PAS ICI
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
CARGU ! MAXELT       !  E !  -> ! NB MAX D'ELEMENTS (CELL,FAC,FBR)     !
CARGU ! LSTELT(MAXELT) TE !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "ppincl.h"
C
C***********************************************************************
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
      INTEGER          IPRFML(NFML,NPRFML), MAXELT, LSTELT(MAXELT)
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
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IEL, MODE, IPHAS , IDIMVE , IESP
C
      DOUBLE PRECISION TINIT, HINIT, COEFE(NGAZEM)
      DOUBLE PRECISION XKENT, XEENT, D2S3
C
C NOMBRE DE PASSAGES DANS LA ROUTINE
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C***********************************************************************
C=======================================================================
C 1.  INITIALISATION VARIABLES LOCALES
C=======================================================================
C
      IPASS = IPASS + 1
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      D2S3 = 2.D0/3.D0
C
C=======================================================================
C 2. INITIALISATION DES INCONNUES :
C      UNIQUEMENT SI ON NE FAIT PAS UNE SUITE
C=======================================================================
C
C RQ IMPORTANTE : pour module electrique, 1 seul passage suffit
C
      IF ( ISUITE.EQ.0 .AND. IPASS.EQ.1 ) THEN
C
        IPHAS = 1
C
C
C --> Initialisation de k et epsilon pas standard
C
C
C ---- Initialisation de k et epsilon
C
        XKENT = 1.D-10
        XEENT = 1.D-10
C
        DO IEL = 1, NCEL
C
C ---- TURBULENCE
C
          IF (ITYTUR(IPHAS).EQ.2) THEN
C
            RTP(IEL,IK(IPHAS))  = XKENT
            RTP(IEL,IEP(IPHAS)) = XEENT
C
          ELSEIF (ITYTUR(IPHAS).EQ.3) THEN
C
            RTP(IEL,IR11(IPHAS)) = D2S3*XKENT
            RTP(IEL,IR22(IPHAS)) = D2S3*XKENT
            RTP(IEL,IR33(IPHAS)) = D2S3*XKENT
            RTP(IEL,IR12(IPHAS)) = 0.D0
            RTP(IEL,IR13(IPHAS)) = 0.D0
            RTP(IEL,IR23(IPHAS)) = 0.D0
            RTP(IEL,IEP(IPHAS))  = XEENT
C
          ELSEIF (ITURB(IPHAS).EQ.50) THEN
C
            RTP(IEL,IK(IPHAS))   = XKENT
            RTP(IEL,IEP(IPHAS))  = XEENT
            RTP(IEL,IPHI(IPHAS)) = D2S3
            RTP(IEL,IFB(IPHAS))  = 0.D0
C
          ELSEIF (ITURB(IPHAS).EQ.60) THEN
C
            RTP(IEL,IK(IPHAS))   = XKENT
            RTP(IEL,IOMG(IPHAS)) = XEENT/CMU/XKENT
C
          ENDIF
C
        ENDDO
C
C --> Enthalpie = H(T0) ou 0
C
C     En arc electrique, on initialise tout le domaine de calcul
C       a T0 avec la 1ere espece
C     En Joule, on initialise l'enthalpie a zero, et il faudra
C       que l'utilisateur intervienne, avec sa loi T->H ou une tabulation.
C
C   -- Calculs de HINIT
C
        IF ( IPPMOD(IELARC).GE.1 ) THEN
          MODE = -1
          TINIT = T0(IPHAS)
          COEFE(1) = 1.D0
          IF ( NGAZG .GT. 1 ) THEN
            DO IESP = 2, NGAZG
              COEFE(IESP) = 0.D0
            ENDDO
          ENDIF
          CALL ELTHHT(MODE,NGAZG,COEFE,HINIT,TINIT)
        ELSE
          HINIT = 0.D0
        ENDIF
C
C    -- Valeurs de l'enthalpie
C
        DO IEL = 1, NCEL
          RTP(IEL,ISCA(IHM)) = HINIT
        ENDDO
C
C
C --> Fractions massiques = 1 ou 0
C
        IF ( NGAZG .GT. 1 ) THEN
          DO IEL = 1, NCEL
            RTP(IEL,ISCA(IYCOEL(1))) = 1.D0
          ENDDO
          DO IESP = 2, NGAZG-1
            DO IEL = 1, NCEL
              RTP(IEL,ISCA(IYCOEL(IESP))) = 0.D0
            ENDDO
          ENDDO
        ENDIF
C
C
C --> Potentiels Electrique = 0
C
C     -- Potentiel Reel
        DO IEL = 1, NCEL
          RTP(IEL,ISCA(IPOTR)) = 0.D0
        ENDDO
C
C     -- Potentiel Imaginaire (Joule)
        IF(IPPMOD(IELJOU).EQ.2 .OR. IPPMOD(IELJOU).EQ.4) THEN
          DO IEL = 1, NCEL
            RTP(IEL,ISCA(IPOTI)) = 0.D0
          ENDDO
        ENDIF
C
C     -- Potentiel vecteur (arc elec. 3D)
        IF ( IPPMOD(IELARC).GE.2 ) THEN
          DO IDIMVE = 1, NDIMVE
            DO IEL = 1, NCEL
              RTP(IEL,ISCA(IPOTVA(IDIMVE))) = 0.D0
            ENDDO
          ENDDO
        ENDIF
C
C
C --> Termes sources = 0
C
C     -- Puissance volumique dissipee par effet Joule W/m3
        DO IEL = 1, NCEL
          PROPCE(IEL,IPPROC(IEFJOU)) = 0.D0
        ENDDO
C
C     -- Force de Laplace (arc elec.)
        IF ( IPPMOD(IELARC).GE.1 ) THEN
          DO IDIMVE = 1, NDIMVE
            DO IEL = 1, NCEL
              PROPCE(IEL,IPPROC(ILAPLA(IDIMVE))) = 0.D0
            ENDDO
          ENDDO
        ENDIF
C
      ENDIF
C
C=======================================================================
C 2.  ON PASSE LA MAIN A L'UTILISATEUR
C=======================================================================
C
      IF (IPASS.EQ.1) THEN
        CALL USELIV
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB , COEFA  , COEFB  ,
     &   RDEVEL , RTUSER , RA     )
      ENDIF
C
C
C----
C FORMATS
C----
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
