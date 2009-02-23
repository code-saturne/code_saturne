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
                        SUBROUTINE MATRDT
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ICONVP , IDIFFP , ISYM   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   COEFBP , FLUMAS , FLUMAB , VISCF  , VISCB  ,
     &   DA     ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC CONSTRUCTION DE LA DIAGONALE DE LA
CFONC    MATRICE DE CONVECTION UPWIND/DIFFUSION/TS
CFONC    POUR DETERMINATION DU PAS DE TEMPS VARIABLE, COURANT, FOURIER
CFONC
c@fonce
C
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
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! ICONVP       ! E  !  -> ! INDICATEUR = 1 CONVECTION, 0 SINON   !
CARGU ! IDIFFP       ! E  !  -> ! INDICATEUR = 1 DIFFUSION , 0 SINON   !
CARGU ! ISYM         ! E  !  -> ! INDICATEUR = 1 MATRICE SYMETRIQUE    !
CARGU !              !    !     !              2 MATRICE NON SYMETRIQUE!
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
CARGU ! COEFBP(NFABOR! TR !  -> ! TAB B DES CL POUR LE PDT CONSIDERE   !
CARGU ! FLUMAS(NFAC) ! TR !  -> ! FLUX DE MASSE AUX FACES INTERNES     !
CARGU ! FLUMAB(NFABOR! TR !  -> ! FLUX DE MASSE AUX FACES DE BORD      !
CARGU ! VISCF(NFAC)  ! TR !  -> ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -> ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! DA (NCELET   ! TR ! <-  ! PARTIE DIAGONALE DE LA MATRICE       !
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
      INCLUDE "vector.h"
      INCLUDE "entsor.h"
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          ICONVP , IDIFFP , ISYM
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
      DOUBLE PRECISION COEFBP(NFABOR)
      DOUBLE PRECISION FLUMAS(NFAC), FLUMAB(NFABOR)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION DA(NCELET )
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IFAC,II,JJ,IEL
      DOUBLE PRECISION FLUI,FLUJ,XAIFA1,XAIFA2
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IF(ISYM.NE.1.AND.ISYM.NE.2) THEN
         WRITE(NFECRA,1000) ISYM
         CALL CSEXIT (1)
      ENDIF
C
      DO IEL = 1, NCEL
        DA(IEL) = 0.D0
      ENDDO
      IF(NCELET.GT.NCEL) THEN
        DO IEL = NCEL+1, NCELET
          DA(IEL) = 0.D0
        ENDDO
      ENDIF
C
C=======================================================================
C 2.    CALCUL DES TERMES EXTRADIAGONAUX INUTILE
C=======================================================================
C
C=======================================================================
C 3.     CONTRIBUTION DES TERMES X-TRADIAGONAUX A LA DIAGONALE
C=======================================================================
C
      IF(ISYM.EQ.2) THEN
C
        IF (IVECTI.EQ.1) THEN
C
!CDIR NODEP
          DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            FLUJ =-0.5D0*( FLUMAS(IFAC) +ABS(FLUMAS(IFAC)) )
            FLUI = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
            XAIFA2 = ICONVP*FLUJ -IDIFFP*VISCF(IFAC)
            XAIFA1 = ICONVP*FLUI -IDIFFP*VISCF(IFAC)
            DA(II) = DA(II) -XAIFA2
            DA(JJ) = DA(JJ) -XAIFA1
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            FLUJ =-0.5D0*( FLUMAS(IFAC) +ABS(FLUMAS(IFAC)) )
            FLUI = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
            XAIFA2 = ICONVP*FLUJ -IDIFFP*VISCF(IFAC)
            XAIFA1 = ICONVP*FLUI -IDIFFP*VISCF(IFAC)
            DA(II) = DA(II) -XAIFA2
            DA(JJ) = DA(JJ) -XAIFA1
          ENDDO
C
        ENDIF
C
      ELSE
C
        IF (IVECTI.EQ.1) THEN
C
!CDIR NODEP
          DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            FLUI = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
            XAIFA1 = ICONVP*FLUI -IDIFFP*VISCF(IFAC)
            DA(II) = DA(II) -XAIFA1
            DA(JJ) = DA(JJ) -XAIFA1
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            FLUI = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
            XAIFA1 = ICONVP*FLUI -IDIFFP*VISCF(IFAC)
            DA(II) = DA(II) -XAIFA1
            DA(JJ) = DA(JJ) -XAIFA1
          ENDDO
C
        ENDIF
C
      ENDIF
C
C=======================================================================
C 4.     CONTRIBUTION DES FACETTES DE BORDS A LA DIAGONALE
C=======================================================================
C
      IF (IVECTB.EQ.1) THEN
C
!CDIR NODEP
        DO IFAC = 1, NFABOR
          II = IFABOR(IFAC)
          FLUI = 0.5D0*( FLUMAB(IFAC) -ABS(FLUMAB(IFAC)) )
          FLUJ =-0.5D0*( FLUMAB(IFAC) +ABS(FLUMAB(IFAC)) )
          DA(II) = DA(II) +ICONVP*(-FLUJ + FLUI*COEFBP(IFAC) )
     &                    +IDIFFP*VISCB(IFAC)*(1.D0-COEFBP(IFAC))
        ENDDO
C
      ELSE
C
C VECTORISATION NON FORCEE
        DO IFAC = 1, NFABOR
          II = IFABOR(IFAC)
          FLUI = 0.5D0*( FLUMAB(IFAC) -ABS(FLUMAB(IFAC)) )
          FLUJ =-0.5D0*( FLUMAB(IFAC) +ABS(FLUMAB(IFAC)) )
          DA(II) = DA(II) +ICONVP*(-FLUJ + FLUI*COEFBP(IFAC) )
     &                    +IDIFFP*VISCB(IFAC)*(1.D0-COEFBP(IFAC))
        ENDDO
C
      ENDIF
C
C--------
C FORMATS
C--------
C
#if defined(_CS_LANG_FR)
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS matrdt                           ',/,
     &'@    *********                                               ',/,
     &'@     APPEL DE matrdt              AVEC ISYM   = ',I10        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#else
C
 1000 FORMAT(
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@'                                                            ,/,
     &'@ @@ WARNING: ABORT IN matrdt'                                ,/,
     &'@    ********'                                                ,/,
     &'@     matrdt CALLED                WITH ISYM   = ',I10        ,/,
     &'@'                                                            ,/,
     &'@  The calculation will not be run.'                          ,/,
     &'@'                                                            ,/,
     &'@  Contact support.'                                          ,/,
     &'@'                                                            ,/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@'                                                            ,/)
C
#endif
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
