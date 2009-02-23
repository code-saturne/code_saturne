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
                        SUBROUTINE RECVMC
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
     &   ROM    , FLUMAS , FLUMAB ,
     &   UX     , UY     , UZ     ,
     &   BX     , BY     , BZ     , COCG   ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RECONSTRUCTION DE LA VITESSE A PARTIR DU FLUX DE MASSE
CFONC      PAR MOINDRES CARRES (VITESSE CONSTANTE PAR ELEMENT)
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
CARGU ! ROM(NCELET   ! TR !  -> ! MASSE VOLUMIQUE AUX CELLULES         !
CARGU ! FLUMAS(NFAC) ! TR !  -> ! FLUX DE MASSE AUX FACES INTERNES     !
CARGU ! FLUMAB(NFABOR! TR !  -> ! FLUX DE MASSE AUX FACES DE BORD      !
CARGU ! UX   UY      ! TR ! <-  ! VITESSE RECONSTRUITE                 !
CARGU ! UZ   (NCELET ! TR !     !                                      !
CARGU ! BX,Y,Z(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! COCG         ! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU !   (NCELET,3,3!    !     !                                      !
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
      INCLUDE "vector.h"
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
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION ROM(NCELET)
      DOUBLE PRECISION FLUMAS(NFAC), FLUMAB(NFABOR)
      DOUBLE PRECISION UX  (NCELET), UY  (NCELET), UZ  (NCELET)
      DOUBLE PRECISION BX(NCELET),   BY(NCELET),   BZ(NCELET)
      DOUBLE PRECISION COCG(NCELET,3,3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          LBLOC
      PARAMETER       (LBLOC = 1024)
C
      INTEGER          IDEBIA, IDEBRA, II, JJ, IEL, IFAC
      INTEGER          IBLOC, NBLOC, IREL, IDIM1, IDIM2
      DOUBLE PRECISION AA(LBLOC,3,3)
      DOUBLE PRECISION A11, A22, A33, A12, A13, A23, UNSDET
      DOUBLE PRECISION COCG11, COCG12, COCG13, COCG21, COCG22, COCG23
      DOUBLE PRECISION COCG31, COCG32, COCG33
      DOUBLE PRECISION SMBX, SMBY, SMBZ, UNSRHO
      DOUBLE PRECISION VECFAC, PFACX, PFACY, PFACZ
C
C***********************************************************************
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. CALCUL DE LA MATRICE
C=======================================================================
C
C   INITIALISATION
C
      DO II = 1, 3
        DO JJ = 1, 3
          DO IEL = 1, NCELET
            COCG(IEL,II,JJ) = 0.D0
          ENDDO
        ENDDO
      ENDDO
C
C   ASSEMBLAGE A PARTIR DES FACETTES FLUIDES
C
      DO IDIM1 = 1, 3
        DO IDIM2 = IDIM1, 3
C
          IF (IVECTI.EQ.1) THEN
C
!CDIR NODEP
            DO IFAC = 1, NFAC
              II = IFACEL(1,IFAC)
              JJ = IFACEL(2,IFAC)
              VECFAC = SURFAC(IDIM1,IFAC)*SURFAC(IDIM2,IFAC)
              COCG(II,IDIM1,IDIM2) = COCG(II,IDIM1,IDIM2) + VECFAC
              COCG(JJ,IDIM1,IDIM2) = COCG(JJ,IDIM1,IDIM2) + VECFAC
            ENDDO
C
          ELSE
C
C VECTORISATION NON FORCEE
            DO IFAC = 1, NFAC
              II = IFACEL(1,IFAC)
              JJ = IFACEL(2,IFAC)
              VECFAC = SURFAC(IDIM1,IFAC)*SURFAC(IDIM2,IFAC)
              COCG(II,IDIM1,IDIM2) = COCG(II,IDIM1,IDIM2) + VECFAC
              COCG(JJ,IDIM1,IDIM2) = COCG(JJ,IDIM1,IDIM2) + VECFAC
            ENDDO
C
          ENDIF
C
          IF (IVECTB.EQ.1) THEN
C
!CDIR NODEP
            DO IFAC = 1, NFABOR
              II = IFABOR(IFAC)
              COCG(II,IDIM1,IDIM2) = COCG(II,IDIM1,IDIM2)
     &                         + SURFBO(IDIM1,IFAC)*SURFBO(IDIM2,IFAC)
            ENDDO
C
          ELSE
C
C VECTORISATION NON FORCEE
            DO IFAC = 1, NFABOR
              II = IFABOR(IFAC)
              COCG(II,IDIM1,IDIM2) = COCG(II,IDIM1,IDIM2)
     &                         + SURFBO(IDIM1,IFAC)*SURFBO(IDIM2,IFAC)
            ENDDO
C
          ENDIF
C
        ENDDO
      ENDDO
C
C
C   SYMETRISATION
C
      DO IEL = 1, NCEL
        COCG(IEL,2,1) = COCG(IEL,1,2)
        COCG(IEL,3,1) = COCG(IEL,1,3)
        COCG(IEL,3,2) = COCG(IEL,2,3)
      ENDDO
C
C=======================================================================
C 2. INVERSION DE LA MATRICE
C=======================================================================
C
C
      NBLOC = NCEL/LBLOC
      IF (NBLOC.GT.0) THEN
        DO IBLOC = 1, NBLOC
          DO II = 1, LBLOC
            IEL = (IBLOC-1)*LBLOC+II
C
            COCG11 = COCG(IEL,1,1)
            COCG12 = COCG(IEL,1,2)
            COCG13 = COCG(IEL,1,3)
            COCG21 = COCG(IEL,2,1)
            COCG22 = COCG(IEL,2,2)
            COCG23 = COCG(IEL,2,3)
            COCG31 = COCG(IEL,3,1)
            COCG32 = COCG(IEL,3,2)
            COCG33 = COCG(IEL,3,3)
C
            A11=COCG22*COCG33-COCG32*COCG23
            A12=COCG32*COCG13-COCG12*COCG33
            A13=COCG12*COCG23-COCG22*COCG13
            A22=COCG11*COCG33-COCG31*COCG13
            A23=COCG21*COCG13-COCG11*COCG23
            A33=COCG11*COCG22-COCG21*COCG12
C
            UNSDET = 1.D0/(COCG11*A11+COCG21*A12+COCG31*A13)
C
            AA(II,1,1) = A11 *UNSDET
            AA(II,1,2) = A12 *UNSDET
            AA(II,1,3) = A13 *UNSDET
            AA(II,2,2) = A22 *UNSDET
            AA(II,2,3) = A23 *UNSDET
            AA(II,3,3) = A33 *UNSDET
C
          ENDDO
C
          DO II = 1, LBLOC
            IEL = (IBLOC-1)*LBLOC+II
            COCG(IEL,1,1) = AA(II,1,1)
            COCG(IEL,1,2) = AA(II,1,2)
            COCG(IEL,1,3) = AA(II,1,3)
            COCG(IEL,2,2) = AA(II,2,2)
            COCG(IEL,2,3) = AA(II,2,3)
            COCG(IEL,3,3) = AA(II,3,3)
          ENDDO
C
        ENDDO
C
      ENDIF
C
      IREL = MOD(NCEL,LBLOC)
      IF (IREL.GT.0) THEN
        IBLOC = NBLOC + 1
        DO II = 1, IREL
          IEL = (IBLOC-1)*LBLOC+II
C
          COCG11 = COCG(IEL,1,1)
          COCG12 = COCG(IEL,1,2)
          COCG13 = COCG(IEL,1,3)
          COCG21 = COCG(IEL,2,1)
          COCG22 = COCG(IEL,2,2)
          COCG23 = COCG(IEL,2,3)
          COCG31 = COCG(IEL,3,1)
          COCG32 = COCG(IEL,3,2)
          COCG33 = COCG(IEL,3,3)
C
          A11=COCG22*COCG33-COCG32*COCG23
          A12=COCG32*COCG13-COCG12*COCG33
          A13=COCG12*COCG23-COCG22*COCG13
          A22=COCG11*COCG33-COCG31*COCG13
          A23=COCG21*COCG13-COCG11*COCG23
          A33=COCG11*COCG22-COCG21*COCG12
C
          UNSDET = 1.D0/(COCG11*A11+COCG21*A12+COCG31*A13)
C
          AA(II,1,1) = A11 *UNSDET
          AA(II,1,2) = A12 *UNSDET
          AA(II,1,3) = A13 *UNSDET
          AA(II,2,2) = A22 *UNSDET
          AA(II,2,3) = A23 *UNSDET
          AA(II,3,3) = A33 *UNSDET
C
        ENDDO
C
        DO II = 1, IREL
          IEL = (IBLOC-1)*LBLOC+II
          COCG(IEL,1,1) = AA(II,1,1)
          COCG(IEL,1,2) = AA(II,1,2)
          COCG(IEL,1,3) = AA(II,1,3)
          COCG(IEL,2,2) = AA(II,2,2)
          COCG(IEL,2,3) = AA(II,2,3)
          COCG(IEL,3,3) = AA(II,3,3)
        ENDDO
      ENDIF
C
C
C         MATRICE SYMETRIQUE
C
      DO IEL = 1, NCEL
        COCG(IEL,2,1) = COCG(IEL,1,2)
        COCG(IEL,3,1) = COCG(IEL,1,3)
        COCG(IEL,3,2) = COCG(IEL,2,3)
      ENDDO
C
C
C=======================================================================
C 3. CALCUL DU SECOND MEMBRE
C=======================================================================
C
      DO IEL = 1, NCELET
        BX(IEL) = 0.D0
        BY(IEL) = 0.D0
        BZ(IEL) = 0.D0
      ENDDO
C
C
C     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES
C
      IF (IVECTI.EQ.1) THEN
C
!CDIR NODEP
        DO IFAC = 1,NFAC
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
          PFACX = FLUMAS(IFAC)*SURFAC(1,IFAC)
          PFACY = FLUMAS(IFAC)*SURFAC(2,IFAC)
          PFACZ = FLUMAS(IFAC)*SURFAC(3,IFAC)
          BX(II) = BX(II) + PFACX
          BY(II) = BY(II) + PFACY
          BZ(II) = BZ(II) + PFACZ
          BX(JJ) = BX(JJ) + PFACX
          BY(JJ) = BY(JJ) + PFACY
          BZ(JJ) = BZ(JJ) + PFACZ
        ENDDO
C
      ELSE
C
C VECTORISATION NON FORCEE
        DO IFAC = 1,NFAC
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
          PFACX = FLUMAS(IFAC)*SURFAC(1,IFAC)
          PFACY = FLUMAS(IFAC)*SURFAC(2,IFAC)
          PFACZ = FLUMAS(IFAC)*SURFAC(3,IFAC)
          BX(II) = BX(II) + PFACX
          BY(II) = BY(II) + PFACY
          BZ(II) = BZ(II) + PFACZ
          BX(JJ) = BX(JJ) + PFACX
          BY(JJ) = BY(JJ) + PFACY
          BZ(JJ) = BZ(JJ) + PFACZ
        ENDDO
C
      ENDIF
C
C
C     ASSEMBLAGE A PARTIR DES FACETTES DE BORD
C
      IF (IVECTB.EQ.1) THEN
C
!CDIR NODEP
        DO IFAC = 1,NFABOR
          II = IFABOR(IFAC)
          BX(II) = BX(II) + FLUMAB(IFAC)*SURFBO(1,IFAC)
          BY(II) = BY(II) + FLUMAB(IFAC)*SURFBO(2,IFAC)
          BZ(II) = BZ(II) + FLUMAB(IFAC)*SURFBO(3,IFAC)
        ENDDO
C
      ELSE
C
C VECTORISATION NON FORCEE
        DO IFAC = 1,NFABOR
          II = IFABOR(IFAC)
          BX(II) = BX(II) + FLUMAB(IFAC)*SURFBO(1,IFAC)
          BY(II) = BY(II) + FLUMAB(IFAC)*SURFBO(2,IFAC)
          BZ(II) = BZ(II) + FLUMAB(IFAC)*SURFBO(3,IFAC)
        ENDDO
C
      ENDIF
C
C=======================================================================
C 4. RESOLUTION
C=======================================================================
C
C
      DO IEL = 1, NCEL
        UNSRHO = 1.D0/ROM(IEL)
        SMBX = BX(IEL)
        SMBY = BY(IEL)
        SMBZ = BZ(IEL)
        UX  (IEL) = (COCG(IEL,1,1)*SMBX+COCG(IEL,1,2)*SMBY
     &              +COCG(IEL,1,3)*SMBZ)*UNSRHO
        UY  (IEL) = (COCG(IEL,2,1)*SMBX+COCG(IEL,2,2)*SMBY
     &              +COCG(IEL,2,3)*SMBZ)*UNSRHO
        UZ  (IEL) = (COCG(IEL,3,1)*SMBX+COCG(IEL,3,2)*SMBY
     &              +COCG(IEL,3,3)*SMBZ)*UNSRHO
      ENDDO
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
