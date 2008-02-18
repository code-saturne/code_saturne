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
                        SUBROUTINE CLDIJP
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NFAC   , NFABOR , NCELET , NCEL   ,
     &   IFACEL , IFABOR ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   SURFAN , SURFBN ,
     &   POND   ,
     &   DIJPF  , DIIPB  , DOFIJ  ,
     &   IA     , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C
C  FONCTION  :
C  --------
c@foncb
CFONC
CFONC     CALCUL DE VECTEURS POUR TERMES DE NON ORTHOGONALITE :
CFONC
CFONC        SOIT UNE FACE A,B,C ET I,J LES CENTRES DES VOLUMES VOISINS.
CFONC          (SEUL I EST DEFINI POUR UNE FACE DE BORD)
CFONC
CFONC        LA FACE EST ORIENTEE DE I VERS J, DE NORMALE NIJ.
CFONC          (LES FACES DE BORD SONT ORIENTEES VERS L'EXTERIEUR)
CFONC        LA NORME DE NIJ VAUT 1.
CFONC        SIJ EST LA SURFACE DE LA FACE ABC.
CFONC        ON DEFINIT I' ET J' LES PROJETES ORTHOGONAUX DE
CFONC          I ET J SUR LA DROITE ORTHOGONALE A ABC PASSANT PAR SON CDG F
CFONC          (SEUL I' EST DEFINI POUR LES FACES DE BORD)
CFONC
CFONC        ON CALCULE ICI LE VECTEUR I'J' AUX FACES INTERNES (DIJPF)
CFONC                       LE VECTEUR II'  AUX FACES DE BORD  (DIIPB)
CFONC                       LE VECTEUR OF   AUX FACES INTERNES (DOFIJ)
CFONC
CFONC        NOTER LES RELATIONS
CFONC                  II' = IG - (IG.NIJ)NIJ
CFONC                  JJ' = JG - (JG.NIJ)NIJ
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0/IDBRA0! E  !  -> ! POINTEUR DE LA PREMIERE CAS LIBRE DES!
CARGU !              !    !     !  TABLEAUX IA/RA                      !
CARGU ! NFAC         ! E  !  -> ! NOMBRE TOTAL DE FACES INTERNES       !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE TOTAL DE FACES         DE BORD!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! IFACEL       ! TE !  -> ! IFACEL(IEL,IFAC) NUM GLOB DE L'ELEMNT!
CARGU ! (2,NFAC)     !    !     !  VOIS IEL (1 OU 2) DE LA FAC INT IFAC!
CARGU ! IFABOR       ! TE !  -> ! IFABOR(IFAC    ) NUM GLOB DE L'ELT   !
CARGU ! NFABOR  )    !    !     !  VOISIN IEL (1) DE LA FAC DE BRD IFAC!
CARGU ! XYZCEN       ! TR !  -> ! COORDS DU "CENTRE" DES NELEM ELEMENTS!
CARGU ! (3,NCELET    !    !     !                                      !
CARGU ! SURFAC       ! TR !  -> ! COORDS DU VECTEUR SURFACE DES NFAC   !
CARGU ! (3,NFAC  )   !    !     ! FACES INTERNES; DIRIGE DU VOIS 1 VERS!
CARGU !              !    !     !  LE VOIS 2 (IFACEL) ; NON UNITAIRE   !
CARGU ! SURFBO       ! TR !  -> ! COORDS DU VECTEUR SURFACE DES NFABOR !
CARGU ! (3,NFABOR)   !    !     !  FACES DE BORD ; DIRIGE VERS         !
CARGU !              !    !     !  L'EXTERIEUR DU DOMAINE; NON UNITAIRE!
CARGU ! SURFAN       ! TR !  -> ! NORME DE SURFAC (SURFACE DES FACES   !
CARGU ! (NFAC    )   !    !     ! INTERNES)                            !
CARGU ! SURFBN       ! TR !  -> ! NORME DE SURFBO (SURFACE DES FACES   !
CARGU ! (NFABOR  )   !    !     !  DE BORD)                            !
CARGU ! CDGFAC       ! TR !  -> ! COORDS DU CENTRE DE GRAVITE DES FACES!
CARGU ! 3,NFAC  )    !    !     !           INTERNES                   !
CARGU ! CDGFBO       ! TR !  -> ! COORDS DU CENTRE DE GRAVITE DES FACES!
CARGU ! 3,NFABOR)    !    !     !           DE  BORD                   !
CARGU ! POND         ! TR ! <-  ! PONDERATION POUR FACE INTERNE        !
CARGU ! (NFAC  )     !    !     !  = D2/(D1+D2) OU D1 ET D2 SONT LES   !
CARGU !              !    !     !  PROJETES SUR LA NORMALE A LA FACE   !
CARGU !              !    !     !  DES VECTEURS DEFINIS RESP. PAR :    !
CARGU !              !    !     !D1(ORIG: VOISIN 1, EXTREMITE: CDG FACE!
CARGU !              !    !     !D2(ORIG: CDG FACE, EXTREMITE: VOISIN 2!
CARGU ! DIJPF        ! TR ! <-  ! VECTEUR I'J' POUR LES FACES INTERNES !
CARGU ! (NDIM,NFAC  )!    !     !                                      !
CARGU ! DIIPB        ! TR ! <-  ! VECTEUR II' POUR LES FACES DE BORD   !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! DOFIJ        ! TR ! <-  ! VECTEUR OF POUR LES FACES INTERNES   !
CARGU ! (NDIM,NFAC  )!    !     ! O : INTERSECTION DE IJ ET LA FACE    !
CARGU !              !    !     ! F : CENTRE DE LA FACE                !
CARGU ! IA           ! TE !  -  ! TABLEAU DE TRAVAIL ENTIER            !
CARGU ! RA           ! TR !  -  ! TABLEAU DE TRAVAIL REEL              !
CARGU .______________.____._____.______________________________________.
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
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER IDBIA0,IDBRA0
      INTEGER NFAC,NFABOR,NCELET,NCEL
      INTEGER IFACEL(2,NFAC)
      INTEGER IFABOR(NFABOR)
      INTEGER IA(*)
      DOUBLE PRECISION XYZCEN(3,NCELET)
      DOUBLE PRECISION SURFAC(3,NFAC),SURFBO(3,NFABOR)
      DOUBLE PRECISION SURFAN(NFAC),SURFBN(NFABOR)
      DOUBLE PRECISION CDGFAC(3,NFAC),CDGFBO(3,NFABOR)
      DOUBLE PRECISION POND(NFAC)
      DOUBLE PRECISION DIJPF(3,NFAC)
      DOUBLE PRECISION DIIPB(3,NFABOR)
      DOUBLE PRECISION DOFIJ(3,NFAC)
      DOUBLE PRECISION RA(*)
C
      INTEGER IFAC,IVOIS1,IVOIS2
      DOUBLE PRECISION SURFNX,SURFNY,SURFNZ
      DOUBLE PRECISION VECIGX,VECIGY,VECIGZ
      DOUBLE PRECISION VECIJX,VECIJY,VECIJZ
      DOUBLE PRECISION DIPJP
      DOUBLE PRECISION PSI
C
C=======================================================================
C 1. FACES INTERNES
C=======================================================================
C
      DO IFAC = 1, NFAC
C
C--->  NUMERO DES VOISINS
C
         IVOIS1 = IFACEL(1,IFAC)
         IVOIS2 = IFACEL(2,IFAC)
C
C---> NORMALE NORMEE
C
         SURFNX = SURFAC(1,IFAC)/SURFAN(IFAC)
         SURFNY = SURFAC(2,IFAC)/SURFAN(IFAC)
         SURFNZ = SURFAC(3,IFAC)/SURFAN(IFAC)
C
C---> IJ
C
         VECIJX = XYZCEN(1,IVOIS2)-XYZCEN(1,IVOIS1)
         VECIJY = XYZCEN(2,IVOIS2)-XYZCEN(2,IVOIS1)
         VECIJZ = XYZCEN(3,IVOIS2)-XYZCEN(3,IVOIS1)
C
C---> DIJPP = IJ.NIJ
C
         DIPJP  = VECIJX*SURFNX + VECIJY*SURFNY +
     &            VECIJZ*SURFNZ
C
C---> DIJPF = (IJ.NIJ).NIJ
C
         DIJPF(1,IFAC) = DIPJP * SURFNX
         DIJPF(2,IFAC) = DIPJP * SURFNY
         DIJPF(3,IFAC) = DIPJP * SURFNZ
C
C---> DOFIJ = OF
C
        DOFIJ(1,IFAC) = CDGFAC(1,IFAC) -
     &                  (POND(IFAC) * XYZCEN(1,IVOIS1) +
     &                  (1.D0-POND(IFAC))*XYZCEN(1,IVOIS2))
        DOFIJ(2,IFAC) = CDGFAC(2,IFAC) -
     &                  (POND(IFAC)*XYZCEN(2,IVOIS1) +
     &                  (1.D0-POND(IFAC))*XYZCEN(2,IVOIS2))
        DOFIJ(3,IFAC) = CDGFAC(3,IFAC) -
     &                  (POND(IFAC)*XYZCEN(3,IVOIS1) +
     &                  (1.D0-POND(IFAC))*XYZCEN(3,IVOIS2))
C
      ENDDO
C
C=======================================================================
C 2. FACES DE BORD
C=======================================================================
C
      DO IFAC = 1, NFABOR
C
C--->  NUMERO DU VOISIN
C
         IVOIS1 = IFABOR(IFAC)
C
C---> NORMALE NORMEE
C
         SURFNX = SURFBO(1,IFAC)/SURFBN(IFAC)
         SURFNY = SURFBO(2,IFAC)/SURFBN(IFAC)
         SURFNZ = SURFBO(3,IFAC)/SURFBN(IFAC)
C
C---> IG
C
         VECIGX = CDGFBO(1,IFAC)-XYZCEN(1,IVOIS1)
         VECIGY = CDGFBO(2,IFAC)-XYZCEN(2,IVOIS1)
         VECIGZ = CDGFBO(3,IFAC)-XYZCEN(3,IVOIS1)
C
C---> PSI = IG.NIJ
C
         PSI = VECIGX*SURFNX+VECIGY*SURFNY+VECIGZ*SURFNZ
C
C---> DIIPB = IG - (IG.NIJ)NIJ
C
         DIIPB(1,IFAC) = VECIGX - PSI*SURFNX
         DIIPB(2,IFAC) = VECIGY - PSI*SURFNY
         DIIPB(3,IFAC) = VECIGZ - PSI*SURFNZ
C
      ENDDO
C
C=======================================================================
C 3. FIN
C=======================================================================
C
      RETURN
      END
c@z
