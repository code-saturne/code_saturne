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
                        SUBROUTINE RIJTHE
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPHAS  , IVAR   , ISOU   , IPP    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , GRAROX , GRAROY , GRAROZ , SMBR   ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC TERMES DE GRAVITE
CFONC    POUR Rij et EPSILON
CFONC VAR  = R11 R22 R33 R12 R13 R23 EP
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
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE PHASE                      !
CARGU ! IVAR         ! E  !  -> ! NUMERO DE VARIABLE                   !
CARGU ! ISOU         ! E  !  -> ! NUMERO DE PASSAGE                    !
CARGU ! IPP          ! E  !  -> ! NUMERO DE VARIABLE POUR SORTIES POST !
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
CARGU ! GRAROX,Y,Z   ! TR !  -> ! TABLEAU DE TRAVAIL POUR GRAD ROM     !
CARGU !  (NCELET)    !    !     !                                      !
CARGU ! SMBR(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! W1...6(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INTEGER          IPHAS  , IVAR   , ISOU   , IPP
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
      DOUBLE PRECISION RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION GRAROX(NCELET), GRAROY(NCELET), GRAROZ(NCELET)
      DOUBLE PRECISION SMBR(NCELET)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IEL
      INTEGER          IR11IP, IR22IP, IR33IP, IR12IP, IR13IP, IR23IP
      INTEGER          IEIPH
C
      DOUBLE PRECISION UNS3, CONST, KSEPS
      DOUBLE PRECISION PRDTUR, R1T, R2T, R3T
      DOUBLE PRECISION G11, G22, G33, G12, G13, G23, GKKS3
      DOUBLE PRECISION G11P, G22P, G33P
      DOUBLE PRECISION PHIT11, PHIT22, PHIT33, PHIT12, PHIT13, PHIT23
      DOUBLE PRECISION AA, BB
C
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IR11IP = IR11(IPHAS)
      IR22IP = IR22(IPHAS)
      IR33IP = IR33(IPHAS)
      IR12IP = IR12(IPHAS)
      IR13IP = IR13(IPHAS)
      IR23IP = IR23(IPHAS)
      IEIPH  = IEP (IPHAS)
C
      IF(ISCALT(IPHAS).GT.0.AND.NSCAL.GE.ISCALT(IPHAS)) THEN
        PRDTUR = SIGMAS(ISCALT(IPHAS))
      ELSE
        PRDTUR = 1.D0
      ENDIF
C
      CONST = -1.5D0*CMU/PRDTUR
      UNS3  = 1.D0/3.D0
C
C=======================================================================
C 2. TERMES POUR RIJ :
C      ROM*VOLUME*dRij/dt =
C                     ... + (Gij - CRIJ3*(Gij-Delta ij Gkk/3))*VOLUME
C            Avec Gij = -(1.5 CMU/PRDTUR) (K/EPS) (Rit Gj + Rjt Gi)
C                 Rit = Rik dROM/dxk (somme sur k)
C=======================================================================
C
C
      IF     (IVAR.EQ.IR11IP) THEN
C
        DO IEL = 1, NCEL
C
          R1T = RTPA(IEL,IR11IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR12IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR13IP)*GRAROZ(IEL)
          R2T = RTPA(IEL,IR12IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR22IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR23IP)*GRAROZ(IEL)
          R3T = RTPA(IEL,IR13IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR23IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR33IP)*GRAROZ(IEL)
C
          KSEPS = (RTPA(IEL,IR11IP)+RTPA(IEL,IR22IP)+RTPA(IEL,IR33IP))
     &           /(2.D0*RTPA(IEL,IEIPH))
C
          G11 = CONST*KSEPS*2.D0*(R1T*GX       )
          G22 = CONST*KSEPS*2.D0*(R2T*GY       )
          G33 = CONST*KSEPS*2.D0*(R3T*GZ       )
          GKKS3 = UNS3*(G11+G22+G33)
C
          PHIT11 = -CRIJ3*(G11-GKKS3)
C
          SMBR(IEL) = SMBR(IEL) + (G11+PHIT11)*VOLUME(IEL)
C
        ENDDO
C
      ELSEIF (IVAR.EQ.IR22IP) THEN
C
        DO IEL = 1, NCEL
C
          R1T = RTPA(IEL,IR11IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR12IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR13IP)*GRAROZ(IEL)
          R2T = RTPA(IEL,IR12IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR22IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR23IP)*GRAROZ(IEL)
          R3T = RTPA(IEL,IR13IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR23IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR33IP)*GRAROZ(IEL)
C
          KSEPS = (RTPA(IEL,IR11IP)+RTPA(IEL,IR22IP)+RTPA(IEL,IR33IP))
     &           /(2.D0*RTPA(IEL,IEIPH))
C
          G11 = CONST*KSEPS*2.D0*(R1T*GX       )
          G22 = CONST*KSEPS*2.D0*(R2T*GY       )
          G33 = CONST*KSEPS*2.D0*(R3T*GZ       )
          GKKS3 = UNS3*(G11+G22+G33)
C
          PHIT22 = -CRIJ3*(G22-GKKS3)
C
          SMBR(IEL) = SMBR(IEL) + (G22+PHIT22)*VOLUME(IEL)
C
        ENDDO
C
      ELSEIF (IVAR.EQ.IR33IP) THEN
C
        DO IEL = 1, NCEL
C
          R1T = RTPA(IEL,IR11IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR12IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR13IP)*GRAROZ(IEL)
          R2T = RTPA(IEL,IR12IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR22IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR23IP)*GRAROZ(IEL)
          R3T = RTPA(IEL,IR13IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR23IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR33IP)*GRAROZ(IEL)
C
          KSEPS = (RTPA(IEL,IR11IP)+RTPA(IEL,IR22IP)+RTPA(IEL,IR33IP))
     &           /(2.D0*RTPA(IEL,IEIPH))
C
          G11 = CONST*KSEPS*2.D0*(R1T*GX       )
          G22 = CONST*KSEPS*2.D0*(R2T*GY       )
          G33 = CONST*KSEPS*2.D0*(R3T*GZ       )
          GKKS3 = UNS3*(G11+G22+G33)
C
          PHIT33 = -CRIJ3*(G33-GKKS3)
C
          SMBR(IEL) = SMBR(IEL) + (G33+PHIT33)*VOLUME(IEL)
C
        ENDDO
C
      ELSEIF (IVAR.EQ.IR12IP) THEN
C
        DO IEL = 1, NCEL
C
          R1T = RTPA(IEL,IR11IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR12IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR13IP)*GRAROZ(IEL)
          R2T = RTPA(IEL,IR12IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR22IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR23IP)*GRAROZ(IEL)
C
          KSEPS = (RTPA(IEL,IR11IP)+RTPA(IEL,IR22IP)+RTPA(IEL,IR33IP))
     &           /(2.D0*RTPA(IEL,IEIPH))
C
          G12 = CONST*KSEPS*     (R1T*GY+R2T*GX)
C
          PHIT12 = -CRIJ3* G12
C
          SMBR(IEL) = SMBR(IEL) + (G12+PHIT12)*VOLUME(IEL)
C
        ENDDO
C
      ELSEIF (IVAR.EQ.IR13IP) THEN
C
        DO IEL = 1, NCEL
C
          R1T = RTPA(IEL,IR11IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR12IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR13IP)*GRAROZ(IEL)
          R3T = RTPA(IEL,IR13IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR23IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR33IP)*GRAROZ(IEL)
C
          KSEPS = (RTPA(IEL,IR11IP)+RTPA(IEL,IR22IP)+RTPA(IEL,IR33IP))
     &           /(2.D0*RTPA(IEL,IEIPH))
C
          G13 = CONST*KSEPS*     (R1T*GZ+R3T*GX)
C
          PHIT13 = -CRIJ3* G13
C
          SMBR(IEL) = SMBR(IEL) + (G13+PHIT13)*VOLUME(IEL)
C
        ENDDO
C
      ELSEIF (IVAR.EQ.IR23IP) THEN
C
        DO IEL = 1, NCEL
C
          R2T = RTPA(IEL,IR12IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR22IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR23IP)*GRAROZ(IEL)
          R3T = RTPA(IEL,IR13IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR23IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR33IP)*GRAROZ(IEL)
C
          KSEPS = (RTPA(IEL,IR11IP)+RTPA(IEL,IR22IP)+RTPA(IEL,IR33IP))
     &           /(2.D0*RTPA(IEL,IEIPH))
C
          G23 = CONST*KSEPS*     (R2T*GZ+R3T*GY)
C
          PHIT23 = -CRIJ3* G23
C
          SMBR(IEL) = SMBR(IEL) + (G23+PHIT23)*VOLUME(IEL)
C
        ENDDO
C
C=======================================================================
C 3. TERMES POUR EPSILON :
C      ROM*VOLUME*dEps/dt =
C                     ... + CEPS1*(EPS/K)*MAX(0,(Gkk/2))*VOLUME
C            Avec Gij = -(1.5 CMU/PRDTUR) (K/EPS) (Rit Gj + Rjt Gi)
C                 Rit = Rik dROM/dxk (somme sur k)
C            On simplifie (EPS/K) en notant
C                GijP = -(1.5 CMU/PRDTUR)         (Rit Gj + Rjt Gi)
C      ROM*VOLUME*dEps/dt =
C                     ... + CEPS1*        MAX(0,(GkkP/2))*VOLUME
C=======================================================================
C
C
      ELSEIF (IVAR.EQ.IEIPH ) THEN
C
        DO IEL = 1, NCEL
C
          R1T = RTPA(IEL,IR11IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR12IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR13IP)*GRAROZ(IEL)
          R2T = RTPA(IEL,IR12IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR22IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR23IP)*GRAROZ(IEL)
          R3T = RTPA(IEL,IR13IP)*GRAROX(IEL)
     &        + RTPA(IEL,IR23IP)*GRAROY(IEL)
     &        + RTPA(IEL,IR33IP)*GRAROZ(IEL)
C
          G11P = CONST*      2.D0*(R1T*GX       )
          G22P = CONST*      2.D0*(R2T*GY       )
          G33P = CONST*      2.D0*(R3T*GZ       )
C
          AA = 0.D0
          BB = 0.5D0*(G11P+G22P+G33P)
          SMBR(IEL) = SMBR(IEL) + CE1*MAX(AA,BB)*VOLUME(IEL)
C
        ENDDO
C
      ENDIF
C
      RETURN
C
      END
c@z
