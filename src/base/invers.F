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
                        SUBROUTINE INVERS
C                       *****************
C     -------------------------------------------------------------
     & ( CNOM   , IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ISYM   , IPOL   , IRESLP , NITMAP , IMGRP  ,
     &   NCYMXP , NITMFP ,
     &   IWARNP , NFECRA , NITERF , ICYCLE , IINVPE ,
     &   EPSILP , RNORM  , RESIDU ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DAM    , XAM    , SMBRP  , VX     ,
     &   DAG    , XAG    , BG     , DXG    ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC APPEL AUX ROUTINE D'INVERSION DE SYSTEMES LINEAIRE
CFONC   MULTIGRILLE + GRADCO OU JACOBI OU BI-CGSTAB
CFONC   GRADCO
CFONC   JACOBI
CFONC   BI-CGSTAB
CFONC   ON SUPPOSE VX INITIALISE EN ENTREE
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! CNOM         ! A  !  -> ! NOM DE LA VARIABLE                   !
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
CARGU ! ISYM         ! E  !  -> ! INDICATEUR = 1 MATRICE SYM           !
CARGU !              !    !     !            = 2 MATRICE NON SYM       !
CARGU ! IPOL         ! E  !  -> ! DEGRE DU POLYNOME POUR PRECOND       !
CARGU !              !    !     !         (0 -> DIAGONAL)              !
CARGU ! IRESLP       ! E  !  -> ! INDICATEUR = 0 GRADCO                !
CARGU !              !    !     !            = 1 JACOBI                !
CARGU !              !    !     !            = 2 CGSTAB                !
CARGU ! NITMAP       ! E  !  -> ! NOMBRE MAX D'ITER POUR RESOL ITERATIV!
CARGU ! IMGRP        ! E  !  -> ! INDICATEUR = 0 PAS DE MULTIGRILLE    !
CARGU ! NCYMXP       ! E  !  -> ! NOMBRE DE CYCLES MAX POUR MULTIGRILLE!
CARGU ! NITMFP       ! E  !  -> ! NOMBRE D ITER SUR MAILLAGE FIN       !
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! NFECRA       ! E  !  -> ! UNITE DU FICHIER SORTIE STD          !
CARGU ! NITERF       ! E  ! <-  ! NOMBRE D'ITERATIONS EFFECTUEES       !
CARGU !              !    !     !  (NON MULTIGRILLE)                   !
CARGU ! ICYCLE       ! E  ! <-  ! NOMBRE DE CYCLES MGM EFFECTUES       !
CARGU ! IINVPE       ! E  !  -> ! INDICATEUR POUR ANNULER LES INCREMENT!
CARGU !              !    !     ! EN PERIODICITE DE ROTATION (=2) OU   !
CARGU !              !    !     ! POUR LES ECHANGER NORMALEMENT DE     !
CARGU !              !    !     ! MANIERE SCALAIRE (=1)                !
CARGU ! EPSILP       ! R  !  -> ! PRECISION POUR RESOL ITER            !
CARGU ! RNORM        ! R  !  -> ! NORMALISATION DU RESIDU              !
CARGU ! RESIDU       ! R  ! <-  ! RESIDU FINAL NON NORME               !
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
CARGU ! IFACLG       ! TE !  -  ! MACRO VOISINS DES MACROS FACES INT   !
CARGU ! (2,NFAC)     !    !     !                                      !
CARGU ! IRESPR(NCELET! TE !  -  ! MACRO EL AUQUEL APPARTIENT UN MICROEL!
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
CARGU ! DAM(NCELET   ! TR !  -> ! DIAGONALE (MAILLAGE FIN SI MGM)      !
CARGU ! XAM(NFAC,ISYM! TR !  -> ! EXTRADIAGONALE (MAILLAGE FIN SI MGM) !
CARGU ! SMBRP(NCELET ! TR !  -> ! SECOND MEMBRE (MAILLAGE FIN SI MGM)  !
CARGU ! VX   (NCELET ! TR ! <-> ! SOLUTION DU SYSTEME                  !
CARGU ! DAG(NCELET)  ! TR !  -  ! DIAGONALE MAILLAGE GROSSIER (MGM)    !
CARGU ! XAG(NFAC,2   ! TR !  -  ! EXTRADIAGONALE MAILLAGE GROSSIER (MGM!
CARGU ! BG   (NCELET ! TR !  -  ! SECOND MEMBRE MAILLAGE GROSSIER (MGM)!
CARGU ! DXG  (NCELET ! TR !  -  ! TABLEAU DE TRAVAIL MAIL GROSSIER (MGM!
CARGU ! W1,2,3,4,5,6 ! TR !  -  ! AUXILIAIRES DE TRAVAIL               !
CARGU !      (NCELET !    !     !                                      !
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
      INCLUDE "optcal.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      CHARACTER*8      CNOM
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          ISYM   , IPOL   , IRESLP , NITMAP
      INTEGER          IMGRP  , NCYMXP , NITMFP
      INTEGER          IWARNP , NFECRA
      INTEGER          NITERF , ICYCLE , IINVPE
      DOUBLE PRECISION EPSILP , RNORM  , RESIDU
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IFACLG(2,NFAC), IRESPR(NCELET)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC ,2)
      DOUBLE PRECISION DAG(NCELET), XAG(NFAC ,2)
      DOUBLE PRECISION SMBRP(NCELET), BG(NCELET)
      DOUBLE PRECISION VX(NCELET), DXG(NCELET)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET),W4(NCELET)
      DOUBLE PRECISION W5(NCELET),W6(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          LNOM
      INTEGER          IRESDS, IRESAS, NITMDS, NITMAS
C
C***********************************************************************
C
C INITIALISATIONS
C
      LNOM = LEN(CNOM)
C
      ICYCLE = 0
      NITERF = 0
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C RESOLUTION
C
      IF( IMGRP.EQ.1 ) THEN
C
        IRESDS = IRESLP
        IRESAS = IRESLP
C
        NITMDS = NITMFP
        NITMAS = NITMFP
C
        CALL RESMGR
C       ===========
     & ( CNOM   , LNOM   , NCELET , NCEL   , NFAC   ,
     &   ISYM   , IRESDS , IRESAS , IRESLP , IPOL   ,
     &   NCYMXP , NITMDS , NITMAS , NITMAP , IINVPE ,
     &   IWARNP , ICYCLE , NITERF , EPSILP , RNORM  , RESIDU ,
     &   IFACEL , SMBRP  , VX     )
C
      ELSEIF(IMGRP.EQ.0) THEN
C
        IF (IRESLP.GE.0 .AND. IRESLP.LE. 2) THEN
C
          CALL RESLIN
C         ===========
     & ( CNOM   , LNOM   , NCELET , NCEL   , NFAC   ,
     &   ISYM   , IRESLP , IPOL   , NITMAP , IINVPE ,
     &   IWARNP , NITERF , EPSILP , RNORM  , RESIDU ,
C                 ------                     ------
     &   IFACEL , DAM    , XAM    , SMBRP  , VX     )
C                                   -----
C
        ELSE
          WRITE(NFECRA,1000) CNOM, IRESLP
          CALL CSEXIT (1)
        ENDIF
C
      ENDIF
C
C
 1000 FORMAT('INVERS APPELE POUR ',A8,' AVEC IRESOL = ', I10)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
