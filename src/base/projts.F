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
                        SUBROUTINE PROJTS
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGU , IMLIGU ,
     &   IWARNU , NFECRA ,
     &   EPSRGU , CLIMGU ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   FEXTX  , FEXTY  , FEXTZ  ,
     &   COEFBP ,
     &   FLUMAS , FLUMAB , VISCF  , VISCB  ,
     &   VISELX , VISELY , VISELZ ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  PROJECTION SUR LES FACES DES TERMES DE FORCE EXTERIEURE
CFONC  GENERANT UNE PRESSION HYDROSTATIQUE
CFONC  EN FAIT, LE TERME CALCULE EST : DTij FEXTij.Sij
CFONC                                       ----   -
CFONC  ET IL EST AJOUTE AU FLUX DE MASSE.
CFONC  LE CALCUL EST FAIT DE MANIERE COMPATIBLE AVEC ITRMAS (POUR LES
CFONC  FACES INTERNES) ET DE MANIERE A CORRIGER L'ERREUR SUR LA CL
CFONC  DE PRESSION EN PAROI (dP/dn=0 n'est pas adapte en fait)
CFONC
c@fonce
C-----------------------------------------------------------------------
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
CARGU ! INIT         ! E  !  -> ! > 0 : INITIALISATION DU FLUX DE MASSE!
CARGU ! INC          ! E  !  -> ! INDICATEUR = 0 RESOL SUR INCREMENT   !
CARGU !              !    !     !              1 SINON                 !
CARGU ! IMRGRA       ! E  !  -> ! INDICATEUR = 0 GRADRC 97             !
CARGU !              ! E  !  -> !            = 1 GRADMC 99             !
CARGU ! ICCOCG       ! E  !  -> ! INDICATEUR = 1 POUR RECALCUL DE COCG !
CARGU !              !    !     !              0 SINON                 !
CARGU ! NSWRGU       ! E  !  -> ! NOMBRE DE SWEEP POUR RECONSTRUCTION  !
CARGU !              !    !     !             DES GRADIENTS            !
CARGU ! IMLIGU       ! E  !  -> ! METHODE DE LIMITATION DU GRADIENT    !
CARGU !              !    !     !  < 0 PAS DE LIMITATION               !
CARGU !              !    !     !  = 0 A PARTIR DES GRADIENTS VOISINS  !
CARGU !              !    !     !  = 1 A PARTIR DU GRADIENT MOYEN      !
CARGU ! IWARNU       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! NFECRA       ! E  !  -> ! UNITE DU FICHIER SORTIE STD          !
CARGU ! EPSRGU       ! R  !  -> ! PRECISION RELATIVE POUR LA           !
CARGU !              !    !     !  RECONSTRUCTION DES GRADIENTS 97     !
CARGU ! CLIMGU       ! R  !  -> ! COEF GRADIENT*DISTANCE/ECART         !
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
CARGU ! COEFBP(NFABOR! TR !  -> ! TABLEAUX DES COND LIM DE PRESSION    !
CARGU ! FLUMAS(NFAC) ! TR ! <-> ! FLUX DE MASSE AUX FACES INTERNES     !
CARGU ! FLUMAB(NFABOR! TR ! <-> ! FLUX DE MASSE AUX FACES DE BORD      !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
CARGU ! RA(*)        ! TR !  -  ! MACRO TABLEAU REEL                   !
CARGU !______________!____!_____!______________________________________!
c@argue
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
      INCLUDE "pointe.h"
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
      INTEGER          INIT   , INC    , IMRGRA , ICCOCG
      INTEGER          NSWRGU , IMLIGU
      INTEGER          IWARNU , NFECRA
      DOUBLE PRECISION EPSRGU , CLIMGU
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION POND
      DOUBLE PRECISION FEXTX(NCELET),FEXTY(NCELET),FEXTZ(NCELET)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION VISELX(NCELET), VISELY(NCELET), VISELZ(NCELET)
      DOUBLE PRECISION COEFBP(NFABOR)
      DOUBLE PRECISION FLUMAS(NFAC), FLUMAB(NFABOR)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC, II, JJ, III
      DOUBLE PRECISION DIJPFX,DIJPFY,DIJPFZ
      DOUBLE PRECISION DIIPX,DIIPY,DIIPZ
      DOUBLE PRECISION DJJPX,DJJPY,DJJPZ
      DOUBLE PRECISION DIST,SURFN
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
C
C
      IF( INIT.EQ.1 ) THEN
        DO IFAC = 1, NFAC
          FLUMAS(IFAC) = 0.D0
        ENDDO
        DO IFAC = 1, NFABOR
          FLUMAB(IFAC) = 0.D0
        ENDDO
C
      ELSEIF(INIT.NE.0) THEN
        WRITE(NFECRA,1000) INIT
        CALL CSEXIT(1)
      ENDIF
C
C=======================================================================
C 2.  CALCUL DU FLUX DE MASSE SANS TECHNIQUE DE RECONSTRUCTION
C=======================================================================
C
      IF( NSWRGU.LE.1 ) THEN
C
C     FLUX DE MASSE SUR LES FACETTES FLUIDES
C
        DO IFAC = 1, NFAC
C
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
C
          FLUMAS(IFAC) =  FLUMAS(IFAC)
     &         + VISCF(IFAC)*(
     &           (CDGFAC(1,IFAC)-XYZCEN(1,II))*FEXTX(II)
     &          +(CDGFAC(2,IFAC)-XYZCEN(2,II))*FEXTY(II)
     &          +(CDGFAC(3,IFAC)-XYZCEN(3,II))*FEXTZ(II)
     &          -(CDGFAC(1,IFAC)-XYZCEN(1,JJ))*FEXTX(JJ)
     &          -(CDGFAC(2,IFAC)-XYZCEN(2,JJ))*FEXTY(JJ)
     &          -(CDGFAC(3,IFAC)-XYZCEN(3,JJ))*FEXTZ(JJ) )
C
        ENDDO
C
C
C     FLUX DE MASSE SUR LES FACETTES DE BORD
C
        DO IFAC = 1, NFABOR
C
          II = IFABOR(IFAC)
          SURFN = RA(ISRFBN-1+IFAC)
          DIST  = RA(IDISTB-1+IFAC)
C
          FLUMAB(IFAC) = FLUMAB(IFAC)+VISCB(IFAC)*DIST/SURFN
     &         *(1.D0-COEFBP(IFAC))*(FEXTX(II)*SURFBO(1,IFAC)
     &         +FEXTY(II)*SURFBO(2,IFAC)+FEXTZ(II)*SURFBO(3,IFAC) )
C
        ENDDO
C
C
      ELSE
C
C
C     FLUX DE MASSE SUR LES FACETTES FLUIDES
C
        DO IFAC = 1, NFAC
C
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
C
          POND = RA(IPOND-1+IFAC)
C
C     recuperation de I'J'
          III = IDIJPF-1+3*(IFAC-1)
          DIJPFX = RA(III+1)
          DIJPFY = RA(III+2)
          DIJPFZ = RA(III+3)
          SURFN = RA(ISRFAN-1+IFAC)
          DIST  = RA(IDIST-1+IFAC)
C
C     calcul de II' et JJ'
          DIIPX = CDGFAC(1,IFAC)-XYZCEN(1,II)-(1.D0-POND)*DIJPFX
          DIIPY = CDGFAC(2,IFAC)-XYZCEN(2,II)-(1.D0-POND)*DIJPFY
          DIIPZ = CDGFAC(3,IFAC)-XYZCEN(3,II)-(1.D0-POND)*DIJPFZ
          DJJPX = CDGFAC(1,IFAC)-XYZCEN(1,JJ)+POND*DIJPFX
          DJJPY = CDGFAC(2,IFAC)-XYZCEN(2,JJ)+POND*DIJPFY
          DJJPZ = CDGFAC(3,IFAC)-XYZCEN(3,JJ)+POND*DIJPFZ
C
          FLUMAS(IFAC) =  FLUMAS(IFAC)
     &         + VISCF(IFAC)*(
     &           (CDGFAC(1,IFAC)-XYZCEN(1,II))*FEXTX(II)
     &          +(CDGFAC(2,IFAC)-XYZCEN(2,II))*FEXTY(II)
     &          +(CDGFAC(3,IFAC)-XYZCEN(3,II))*FEXTZ(II)
     &          -(CDGFAC(1,IFAC)-XYZCEN(1,JJ))*FEXTX(JJ)
     &          -(CDGFAC(2,IFAC)-XYZCEN(2,JJ))*FEXTY(JJ)
     &          -(CDGFAC(3,IFAC)-XYZCEN(3,JJ))*FEXTZ(JJ) )
     &         +SURFN/DIST*0.5D0*(
     &       (DJJPX-DIIPX)*(VISELX(II)*FEXTX(II)+VISELX(JJ)*FEXTX(JJ))
     &      +(DJJPY-DIIPY)*(VISELY(II)*FEXTY(II)+VISELY(JJ)*FEXTY(JJ))
     &      +(DJJPZ-DIIPZ)*(VISELZ(II)*FEXTZ(II)+VISELZ(JJ)*FEXTZ(JJ)))
C
        ENDDO
C
C
C     FLUX DE MASSE SUR LES FACETTES DE BORD
C
        DO IFAC = 1, NFABOR
C
          II = IFABOR(IFAC)
          SURFN = RA(ISRFBN-1+IFAC)
          DIST  = RA(IDISTB-1+IFAC)
C
          FLUMAB(IFAC) = FLUMAB(IFAC)+VISCB(IFAC)*DIST/SURFN
     &         *(1.D0-COEFBP(IFAC))*(FEXTX(II)*SURFBO(1,IFAC)
     &         +FEXTY(II)*SURFBO(2,IFAC)+FEXTZ(II)*SURFBO(3,IFAC) )
C
        ENDDO
      ENDIF
C
C
C
C--------
C FORMATS
C--------
C
 1000 FORMAT('PROJTS APPELE AVEC INIT =',I10)
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
