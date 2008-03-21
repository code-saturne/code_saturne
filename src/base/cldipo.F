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
                        SUBROUTINE CLDIPO
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NFAC   , NFABOR , NCELET , NCEL   ,
     &   IFACEL , IFABOR ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , SURFAN , SURFBN ,
     &   DIST   , DISTBO , POND   ,
     &   IA     , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C
C  FONCTION  :
C  --------
c@foncb
CFONC          CALCUL DES DISTANCES RELATIVES AUX FACES
CFONC                      ET DES PONDERATIONS ASSOCIEES
CFONC
c@fonce
C
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0/IDBRA0! E  !  -> ! POINTEUR DE LA PREMIERE CAS LIBRE DES!
CARGU !              !    !     !  TABLEAUX IA/RA                      !
CARGU ! NFAC  /NFABOR! E  !  -> ! NOMBRE TOTAL DE FACES INTERNES/DE BOR!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! IFACEL       ! TE !  -> ! IFACEL(IEL,IFAC) NUM GLOB DE L'ELEMNT!
CARGU ! (2,NFAC)     !    !     !  VOIS IEL (1 OU 2) DE LA FAC INT IFAC!
CARGU ! IFABOR       ! TE !  -> ! IFABOR(IFAC) NUM GLOB DE L'ELEMENT   !
CARGU ! NFAC  ,2)    !    !     !  VOISIN IEL (1) DE LA FAC DE BRD IFAC!
CARGU ! XYZCEN       ! TR !  -> ! COORDS DU "CENTRE" DES NELEM ELEMENTS!
CARGU ! (3,NCELET )  !    !     !                                      !
CARGU ! SURFAC       ! TR !  -> ! COORDS DU VECTEUR SURFACE DES NFAGLO !
CARGU ! (3,NFAC  )   !    !     ! FACES INTERNES; DIRIGE DU VOIS 1 VERS!
CARGU !              !    !     !  LE VOISIN 2 (IFACEL) ; NON UNITAIRE !
CARGU ! SURFBO       ! TR !  -> ! COORDS DU VECTEUR SURFACE DES NFAGBR !
CARGU ! (3,NFABOR)   !    !     !  FACES DE BRD;DIRIGE VERS L'EXTERIEUR!
CARGU !              !    !     !  DU DOMAINE ; NON UNITAIRE           !
CARGU ! CDGFAC       ! TR !  -> ! COORDS DU CENTRE DE GRAVITE DES FACES!
CARGU ! (NDIM,NFAC  )!    !     !                              INTERNES!
CARGU ! CDGFBO       ! TR !  -> ! COORDS DU CENTRE DE GRAVITE DES FACES!
CARGU ! (NDIM,NFABOR)!    !     !                              DE BORD !
CARGU ! DIST         ! TR ! <-  ! PROD.SCAL. DE LA NORMALE NORMEE A UNE!
CARGU ! (NFAC  )     !    !     !  FACE INTERNE ET DU VECTEUR DEFINI   !
CARGU !              !    !     !  PAR LES VOISINS (VOISIN 1 : ORIGINE,!
CARGU !              !    !     !  VOISIN 2 : EXTREMITE)               !
CARGU ! POND         ! TR ! <-  ! PONDERATION POUR FACE INTERNE        !
CARGU ! (NFAC  )     !    !     !  = D2/(D1+D2) OU D1 ET D2 SONT LES   !
CARGU !              !    !     !  PROJETES SUR LA NORMALE A LA FACE   !
CARGU !              !    !     !  DES VECTEURS DEFINIS RESP. PAR :    !
CARGU !              !    !     !D1(ORIG: VOISIN 1, EXTREMITE: CDG FACE!
CARGU !              !    !     !D2(ORIG: CDG FACE, EXTREMITE: VOISIN 2!
CARGU ! DISTBO       ! TR ! <-  ! PROD. SCAL. DE LA NORMALE NORMEE     !
CARGU ! (NFABOR)     !    !     !  UNE FACE DE BORD ET DU VECTEUR      !
CARGU !              !    !     !  VOISIN CDG FACE ((VOISIN : ORIGINE, !
CARGU !              !    !     !  CDG FACE : EXTREMITE)               !
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
      INCLUDE "paramx.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
C
      INTEGER           IDBIA0, IDBRA0
      INTEGER           NFAC  , NFABOR, NCELET, NCEL
      INTEGER           IFACEL(2,NFAC), IFABOR(NFABOR)
      INTEGER           IA(*)
      DOUBLE PRECISION  XYZCEN(3,NCELET)
      DOUBLE PRECISION  SURFAC(3,NFAC  ),SURFBO(3,NFABOR)
      DOUBLE PRECISION  CDGFAC(3,NFAC  ),CDGFBO(3,NFABOR)
      DOUBLE PRECISION  SURFAN(NFAC    ),SURFBN(NFABOR)
      DOUBLE PRECISION  DIST(NFAC)
      DOUBLE PRECISION  DISTBO(NFABOR)
      DOUBLE PRECISION  POND(NFAC)
      DOUBLE PRECISION  RA(*)
C
      INTEGER           IFAC, IVOIS1,IVOIS2,IPOND0
      DOUBLE PRECISION  SURFX,SURFY,SURFZ, SURFN,DIST2F
      DOUBLE PRECISION  XVN,YVN,ZVN,XVV,YVV,ZVV
C
C=======================================================================
C 1.  TRAITEMENT DES FACES INTERNES
C=======================================================================
C
      IPOND0 = 0
      DO IFAC = 1, NFAC
C
         SURFX = SURFAC(1,IFAC)
         SURFY = SURFAC(2,IFAC)
         SURFZ = SURFAC(3,IFAC)
         SURFN = SURFAN(IFAC)
C
         IVOIS1 = IFACEL(1,IFAC)
         IVOIS2 = IFACEL(2,IFAC)
C
C (CDG FACE,VOISIN2)
         XVN = XYZCEN(1,IVOIS2) - CDGFAC(1,IFAC)
         YVN = XYZCEN(2,IVOIS2) - CDGFAC(2,IFAC)
         ZVN = XYZCEN(3,IVOIS2) - CDGFAC(3,IFAC)
C
C (PRODUIT SCALAIRE AVEC LA NORMALE NORMEE)
         DIST2F       =    (XVN*SURFX + YVN*SURFY + ZVN*SURFZ)/SURFN
C
C (VOISIN1, VOISIN2)
         XVV = XYZCEN(1,IVOIS2) - XYZCEN(1,IVOIS1)
         YVV = XYZCEN(2,IVOIS2) - XYZCEN(2,IVOIS1)
         ZVV = XYZCEN(3,IVOIS2) - XYZCEN(3,IVOIS1)
C
C (PRODUIT SCALAIRE AVEC LA NORMALE NORMEE)
         DIST(IFAC) =    (XVV*SURFX + YVV*SURFY + ZVV*SURFZ)/SURFN
C
         IF (ABS(DIST(IFAC)).GE.EPZERO) THEN
            POND(IFAC) = DIST2F / DIST(IFAC)
         ELSE
            IPOND0 = IPOND0 + 1
            POND(IFAC) = 0.5D0
         ENDIF
C
      ENDDO
C
C=======================================================================
C 2.  TRAITEMENT DES FACES DE BORD
C=======================================================================
C
      DO IFAC = 1, NFABOR
C
         SURFX = SURFBO(1,IFAC)
         SURFY = SURFBO(2,IFAC)
         SURFZ = SURFBO(3,IFAC)
         SURFN = SURFBN(IFAC)
C
         IVOIS1 = IFABOR(IFAC)
C
C (VOISIN 1, CDG FACE)
         XVN = CDGFBO(1,IFAC) - XYZCEN(1,IVOIS1)
         YVN = CDGFBO(2,IFAC) - XYZCEN(2,IVOIS1)
         ZVN = CDGFBO(3,IFAC) - XYZCEN(3,IVOIS1)
C
         DISTBO(IFAC) =    (XVN*SURFX + YVN*SURFY + ZVN*SURFZ)/SURFN
C
      ENDDO
C
C=======================================================================
C 3.  FIN
C=======================================================================
C
      IF(IPOND0.NE.0) THEN
        WRITE(NFECRA,1000) IPOND0
      ENDIF
C
 1000 FORMAT(' CLDIPO : ',I10,' DISTANCES NULLES ENTRE CENTRES ',/,
     &       '          ON MET POND A 0.5 ')
      RETURN
      END
c@z
