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
                        SUBROUTINE GRADMC
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NFAC   , NFABOR , NCELBR ,
     &   INC    , ICCOCG , NSWRGP , IDIMTE , ITENSO , IPHYDP , IMRGRA ,
     &   IWARNP , NFECRA , EPSRGP , EXTRAP ,
     &   IFACEL , IFABOR , ICELBR , IPCVSE , IELVSE , ISYMPA ,
     &   VOLUME , SURFAC , SURFBO , SURFBN , POND   ,
     &   DIST   , DISTBR , DIJPF  , DIIPB  ,
     &   FEXTX  , FEXTY  , FEXTZ  ,
     &   XYZCEN , CDGFAC , CDGFBO , COEFAP , COEFBP , PVAR   ,
     &   COCGB  , COCG   ,
     &   DPDX   , DPDY   , DPDZ   ,
     &   BX     , BY     , BZ     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC CALCUL DU GRADIENT CELLULE PAR RECONSTRUCTION MOINDRES CARRES 99
CFONC   SUR UN SUPPORT ETENDU
CFONC AVEC PRISE EN COMPTE EVENTUELLE D'UN TERME DE FORCE VOLUMIQUE
CFONC   GENERANT UNE COMPOSANTE DE PRESSION HYDROSTATIQUE
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU  MOINS    !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! INC          ! E  !  -> ! INDICATEUR = 0 RESOL SUR INCREMENT   !
CARGU !              !    !     !              1 SINON                 !
CARGU ! ICCOCG       ! E  !  -> ! INDICATEUR = 1 POUR RECALCUL DE COCG !
CARGU !              !    !     !              0 SINON                 !
CARGU ! NSWRGP       ! E  !  -> ! NOMBRE DE SWEEP POUR RECONSTRUCTION  !
CARGU !              !    !     !             DES GRADIENTS            !
CARGU ! IDIMTE       ! E  !  -> ! DIMENSION DE LA VARIBLE (MAXIMUM 3)  !
CARGU !              !    !     ! 0 : SCALAIRE (VAR11) OU ASSIMILE     !
CARGU !              !    !     !     SCALAIRE                         !
CARGU !              !    !     ! 1 : VECTEUR (VAR11,VAR22,VAR33)      !
CARGU !              !    !     ! 2 : TENSEUR D'ORDRE 2 (VARIJ)        !
CARGU ! ITENSO       ! E  !  -> ! POUR L'EXPLICITATION DE LA ROTATION  !
CARGU !              !    !     ! 0 : SCALAIRE (VAR11)                 !
CARGU !              !    !     ! 1 : COMPOSANTE DE VECTEUR OU DE      !
CARGU !              !    !     !     TENSEUR (VAR11) IMPLCITE POUR LA !
CARGU !              !    !     !     TRANSLATION                      !
CARGU !              !    !     !11 : REPREND LE TRAITEMENT ITENSO=1 ET!
CARGU !              !    !     !     COMPOSANTE DE VECTEUR OU DE      !
CARGU !              !    !     !     TENSEUR (VAR11) ANNULEE  POUR LA !
CARGU !              !    !     !     ROTATION                         !
CARGU !              !    !     ! 2 : VECTEUR (VAR11 ET VAR22 ET VAR33)!
CARGU !              !    !     !     IMPLICITE POUR LA TRANSLATION    !
CARGU ! IPHYDP       ! E  !  -> ! INDICATEUR DE PRISE EN COMPTE DE LA  !
CARGU !              !    !     ! PRESSION HYDROSTATIQUE               !
CARGU ! IMRGRA       ! E  !  -> ! METHODE DE RECONSTRUCTION DU GRADIENT!
CARGU !              !    !     !  0 RECONSTRUCTION 97                 !
CARGU !              !    !     !  1 MOINDRE CARRE                     !
CARGU !              !    !     !  2 MOINDRE CARRE SUPPORT ETENDU      !
CARGU !              !    !     !  3 MOINDRE CARRE SUPPORT ETENDU REDUI!
CARGU !              !    !     !  4 RECONSTR AVEC INIT MOINDRES CARRES!
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! NFECRA       ! E  !  -> ! UNITE DU FICHIER SORTIE STD          !
CARGU ! EPSRGP       ! R  !  -> ! PRECISION RELATIVE POUR LA           !
CARGU !              !    !     !  RECONSTRUCTION DES GRADIENTS 97     !
CARGU ! EXTRAP       ! R  !  -> ! COEF EXTRAP GRADIENT                 !
CARGU ! IFACEL(2,NFAC! TE !  -> ! No DES ELTS VOISINS D'UNE FACE INTERN!
CARGU ! IFABOR(NFABOR! TE !  -> ! No DE L'ELT VOISIN D'UNE FACE DE BORD!
CARGU ! ICELBR       ! TE !  -> ! NUMERO GLOBAL DES ELEMENTS AYANT AU  !
CARGU ! (NCELBR)     !    !     !  MOINS UNE FACE DE BORD              !
CARGU !              !    !     !  RANGES PAR NUMERO CROISSANT         !
CARGU ! IPCVSE       ! TE !  -> ! POSITION DANS IELVSE DES VOISINS     !
CARGU ! (NCEL  )     !    !     !  ETENDUS DES CELLULES                !
CARGU ! IELVSE (*)   ! TE !  -> ! NUMERO DES VOISINS ETENDUS DES CELLUL!
CARGU ! ISYMPA(NFABOR! TE !  -> ! NUL SUR LES SYMETRIES                !
CARGU ! VOLUME(NCELET! TR !  -> ! VOLUME DES ELEMENTS                  !
CARGU ! SURFAC(3,NFAC! TR !  -> ! SURF VECTORIELLE DES SURFACES INTERNE!
CARGU ! SURFBO       ! TR !  -> ! SURF VECTORIELLE DES SURFACES DE BORD!
CARGU !   (3,NFABOR) !    !     !                                      !
CARGU ! SURFBN       ! TR !  -> ! NORME DE LA SURFACE DES FACES DE BORD!
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! POND(NFAC)   ! TR !  -> ! PONDERATION GEOMETRIQUE (ENTRE 0 ET 1!
CARGU ! DIST(NFAC)   ! TR !  -> ! DIST ENTRE LES PROJECTIONS ORTHOGONAL!
CARGU !              !    !     !  SUR LA NORMALE A UNE FACE DES CENTRE!
CARGU !              !    !     !  VOLUMES VOISINS                     !
CARGU ! DISTBR(NFABOR! TR !  -> ! DIST DU CENTRE A LA FACE DE BORD     !
CARGU ! DIJPF(3,NFAC)! TR !  -> ! VECT I'J', I' (RESP. J') PROJECTION  !
CARGU !              !    !     !  DU CENTRE I (RESP. J) SUR LA NORMALE!
CARGU !              !    !     !  A LA FACE INTERNE                   !
CARGU ! DIIPB        ! TR !  -> ! VECT II', II PROJECTION DU CENTRE I  !
CARGU !   (3,NFABOR) !    !     !  SUR LA NORMALE A LA FACE DE BORD    !
CARGU ! XYZCEN       ! TR !  -> ! POINT ASSOCIES AUX VOLUMES DE CONTROL!
CARGU !  (3,NCELET   !    !     !                                      !
CARGU ! CDGFAC       ! TR !  -> ! POINT ASSOCIES AUX FACETTES FLUIDES  !
CARGU !  (3,NFAC)    !    !     !                                      !
CARGU ! CDGFBO       ! TR !  -> ! POINTS ASSOCIES AUX FACETTES DE BORD !
CARGU !  (3,NFABOR)  !    !     !                                      !
CARGU ! COEFAP, B    ! TR !  -> ! TABLEAUX DES COND LIM POUR PVAR      !
CARGU !   (NFABOR)   !    !     !  SUR LA NORMALE A LA FACE DE BORD    !
CARGU ! PVAR  (NCELET! TR !  -> ! VARIABLE (NCELET + V. ETENDU EVENTUEL!
CARGU ! FEXTX,Y,Z    ! TR !  -> ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU !   (NCELET)   !    !     !  HYDROSTATIQUE                       !
CARGU ! COCGB        ! TR ! <-> ! CONTRIBUTION DES FACES INTERNES A    !
CARGU !  (NCELBR,3,3)!    !     !  COCG (CELLULES DE BORD)             !
CARGU ! COCG         ! TR ! <-> ! COUPLAGE DES COMPOSANTES DU GRADIENT !
CARGU !  NCELET,3,3  !    !     ! MODIFIE EVENTUELLEMENT AUX BORDS     !
CARGU ! DPDX DPDY    ! TR ! <-  ! GRADIENT DE PVAR                     !
CARGU ! DPDZ (NCELET ! TR !     !                                      !
CARGU ! BX,Y,Z(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR LE GRAD DE P !
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
      INCLUDE "pointe.h"
      INCLUDE "vector.h"
      INCLUDE "albase.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET , NCEL   , NFAC   , NFABOR , NCELBR
      INTEGER          INC    , ICCOCG , NSWRGP
      INTEGER          IDIMTE , ITENSO , IPHYDP , IMRGRA
      INTEGER          IWARNP , NFECRA
      DOUBLE PRECISION EPSRGP , EXTRAP
C
      INTEGER          IFACEL(2,NFAC),IFABOR(NFABOR),ICELBR(NCELBR)
      INTEGER          IPCVSE(*), IELVSE(*), ISYMPA(NFABOR)
      DOUBLE PRECISION VOLUME(NCELET), SURFAC(3,NFAC)
      DOUBLE PRECISION SURFBO(3,NFABOR), SURFBN(NFABOR)
      DOUBLE PRECISION POND(NFAC), DIST(NFAC), DISTBR(NFABOR)
      DOUBLE PRECISION DIJPF(3,NFAC), DIIPB(3,NFABOR)
      DOUBLE PRECISION XYZCEN(3,*),CDGFAC(3,NFAC),CDGFBO(3,NFABOR)
      DOUBLE PRECISION COEFAP(NFABOR), COEFBP(NFABOR), PVAR(*)
      DOUBLE PRECISION FEXTX(*),FEXTY(*),FEXTZ(*)
      DOUBLE PRECISION COCGB(NCELBR,3,3), COCG(NCELET,3,3)
      DOUBLE PRECISION DPDX (NCELET),DPDY (NCELET),DPDZ (NCELET)
      DOUBLE PRECISION BX(NCELET),BY(NCELET),BZ(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          LBLOC
      PARAMETER       (LBLOC = 1024)
C
      INTEGER          II    , JJ    , IEL   , IELB  , IFAC  , IPOS
      INTEGER          IBLOC , NBLOC , IREL  , IIII  , INVEMX
      DOUBLE PRECISION PFAC  , RKIJ  , DSIJ(3), AA(LBLOC,3,3)
      DOUBLE PRECISION A11   , A22   , A33   , A12   , A13   , A23
      DOUBLE PRECISION COCG11, COCG12, COCG13, COCG21, COCG22, COCG23
      DOUBLE PRECISION COCG31, COCG32, COCG33
      DOUBLE PRECISION USDIJ2, PFSX  , PFSY  , PFSZ
      DOUBLE PRECISION UNSDIJ, UNSSBN, UMCBSD, UNSDET
      DOUBLE PRECISION PFAC1 , PFAC2 , PFAC3 , VECFAC, UNSVOL, EXTRAB
      DOUBLE PRECISION PTMIDX, PTMIDY, PTMIDZ
C
C INDICATEUR DE PREMIER PASSAGE
      INTEGER          INICOC
      DATA             INICOC /1/
      SAVE             INICOC
C
C VERIFICATION QUE NCEL > NBRE DE VOISINS ETENDUS D'UNE CELLULE
      INTEGER          ICESVE
      DATA             ICESVE /-1/
      SAVE             ICESVE
C
C***********************************************************************
C
C=======================================================================
C 0. TRAITEMENT DE EXTRAG POUR LES SYMETRIES
C=======================================================================
C
C POUR LES CALCULS 2D EN PARTICULIER, SI ON EXTRAPOLE LE GRADIENT DE
C   PRESSION, ON SE TROUVE AVEC UNE MATRICE COCG NON INVERSIBLE, A
C   CAUSE DE LA TROISIEME DIRECTION.
C
C COMME L'ON N'A EN PRATIQUE QU'UNE SEULE PHASE,
C   ET QUE EXTRAG = 0 ou 1 (1 POUR LA PRESSION UNIQUEMENT)
C   ON SE PLACE IMPLICITEMENT DANS CETTE SITUATION, AVEC UN STOP
C   DANS VERINI AU CAS OU CE NE SERAIT PAS LE CAS.
C
C LA MODIFICATION CONSISTE A MULTIPLIER EXTRAP PAR ISYMPA QUI EST
C   NUL SUR LES SYMETRIES : ON N'EXTRAPOLE DONC PAS LE GRADIENT SUR CES
C   FACES.
C
C=======================================================================
C 1. INITIALISATION DU GRADIENT
C=======================================================================
C
      IF( NSWRGP.LE.1 ) THEN
C
        DO IEL = 1, NCELET
          BX  (IEL) = 0.D0
          BY  (IEL) = 0.D0
          BZ  (IEL) = 0.D0
        ENDDO
C
C  CAS STANDARD, SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
C  ===============================================================
        IF (IPHYDP.EQ.0) THEN
C
C     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES
C
        IF (IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
            DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            PFAC = POND(IFAC)*PVAR(II) +(1.D0-POND(IFAC))*PVAR(JJ)
            PFAC1 = PFAC*SURFAC(1,IFAC)
            PFAC2 = PFAC*SURFAC(2,IFAC)
            PFAC3 = PFAC*SURFAC(3,IFAC)
            BX(II) = BX(II) +PFAC1
            BY(II) = BY(II) +PFAC2
            BZ(II) = BZ(II) +PFAC3
            BX(JJ) = BX(JJ) -PFAC1
            BY(JJ) = BY(JJ) -PFAC2
            BZ(JJ) = BZ(JJ) -PFAC3
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
            DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            PFAC = POND(IFAC)*PVAR(II) +(1.D0-POND(IFAC))*PVAR(JJ)
            PFAC1 = PFAC*SURFAC(1,IFAC)
            PFAC2 = PFAC*SURFAC(2,IFAC)
            PFAC3 = PFAC*SURFAC(3,IFAC)
            BX(II) = BX(II) +PFAC1
            BY(II) = BY(II) +PFAC2
            BZ(II) = BZ(II) +PFAC3
            BX(JJ) = BX(JJ) -PFAC1
            BY(JJ) = BY(JJ) -PFAC2
            BZ(JJ) = BZ(JJ) -PFAC3
          ENDDO
C
        ENDIF
C
C     ASSEMBLAGE A PARTIR DES FACETTES DE BORD
C
        IF (IVECTB.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
            DO IFAC = 1,NFABOR
            II = IFABOR(IFAC)
            PFAC = INC*COEFAP(IFAC) +COEFBP(IFAC)*PVAR(II)
            BX(II) = BX(II) +PFAC*SURFBO(1,IFAC)
            BY(II) = BY(II) +PFAC*SURFBO(2,IFAC)
            BZ(II) = BZ(II) +PFAC*SURFBO(3,IFAC)
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
            DO IFAC = 1,NFABOR
            II = IFABOR(IFAC)
            PFAC = INC*COEFAP(IFAC) +COEFBP(IFAC)*PVAR(II)
            BX(II) = BX(II) +PFAC*SURFBO(1,IFAC)
            BY(II) = BY(II) +PFAC*SURFBO(2,IFAC)
            BZ(II) = BZ(II) +PFAC*SURFBO(3,IFAC)
          ENDDO
C
        ENDIF
C
C  CAS AVEC PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
C  =====================================================
        ELSE
C
C     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES
C
          IF (IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
            DO IFAC = 1,NFAC
              II = IFACEL(1,IFAC)
              JJ = IFACEL(2,IFAC)
              PFAC   = POND(IFAC)*(PVAR(II)
     &             -(XYZCEN(1,II)-CDGFAC(1,IFAC))*FEXTX(II)
     &             -(XYZCEN(2,II)-CDGFAC(2,IFAC))*FEXTY(II)
     &             -(XYZCEN(3,II)-CDGFAC(3,IFAC))*FEXTZ(II))
     &             +(1.D0-POND(IFAC))*(PVAR(JJ)
     &             -(XYZCEN(1,JJ)-CDGFAC(1,IFAC))*FEXTX(JJ)
     &             -(XYZCEN(2,JJ)-CDGFAC(2,IFAC))*FEXTY(JJ)
     &             -(XYZCEN(3,JJ)-CDGFAC(3,IFAC))*FEXTZ(JJ))
              PFAC1  = PFAC*SURFAC(1,IFAC)
              PFAC2  = PFAC*SURFAC(2,IFAC)
              PFAC3  = PFAC*SURFAC(3,IFAC)
              BX(II) = BX(II) +PFAC1
              BY(II) = BY(II) +PFAC2
              BZ(II) = BZ(II) +PFAC3
              BX(JJ) = BX(JJ) -PFAC1
              BY(JJ) = BY(JJ) -PFAC2
              BZ(JJ) = BZ(JJ) -PFAC3
            ENDDO
C
          ELSE
C
C VECTORISATION NON FORCEE
            DO IFAC = 1,NFAC
              II = IFACEL(1,IFAC)
              JJ = IFACEL(2,IFAC)
              PFAC   = POND(IFAC)*(PVAR(II)
     &             -(XYZCEN(1,II)-CDGFAC(1,IFAC))*FEXTX(II)
     &             -(XYZCEN(2,II)-CDGFAC(2,IFAC))*FEXTY(II)
     &             -(XYZCEN(3,II)-CDGFAC(3,IFAC))*FEXTZ(II))
     &             +(1.D0-POND(IFAC))*(PVAR(JJ)
     &             -(XYZCEN(1,JJ)-CDGFAC(1,IFAC))*FEXTX(JJ)
     &             -(XYZCEN(2,JJ)-CDGFAC(2,IFAC))*FEXTY(JJ)
     &             -(XYZCEN(3,JJ)-CDGFAC(3,IFAC))*FEXTZ(JJ))
              PFAC1  = PFAC*SURFAC(1,IFAC)
              PFAC2  = PFAC*SURFAC(2,IFAC)
              PFAC3  = PFAC*SURFAC(3,IFAC)
              BX(II) = BX(II) +PFAC1
              BY(II) = BY(II) +PFAC2
              BZ(II) = BZ(II) +PFAC3
              BX(JJ) = BX(JJ) -PFAC1
              BY(JJ) = BY(JJ) -PFAC2
              BZ(JJ) = BZ(JJ) -PFAC3
            ENDDO
C
          ENDIF
C
C     ASSEMBLAGE A PARTIR DES FACETTES DE BORD
C
          IF (IVECTB.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
            DO IFAC = 1,NFABOR
              II = IFABOR(IFAC)
              PFAC = INC*COEFAP(IFAC) +COEFBP(IFAC)*(PVAR(II)
     &             -(XYZCEN(1,II)-CDGFBO(1,IFAC))*FEXTX(II)
     &             -(XYZCEN(2,II)-CDGFBO(2,IFAC))*FEXTY(II)
     &             -(XYZCEN(3,II)-CDGFBO(3,IFAC))*FEXTZ(II) )
              BX(II) = BX(II) +PFAC*SURFBO(1,IFAC)
              BY(II) = BY(II) +PFAC*SURFBO(2,IFAC)
              BZ(II) = BZ(II) +PFAC*SURFBO(3,IFAC)
            ENDDO
C
          ELSE
C
C  VECTORISATION NON FORCEE
            DO IFAC = 1,NFABOR
              II = IFABOR(IFAC)
              PFAC = INC*COEFAP(IFAC) +COEFBP(IFAC)*(PVAR(II)
     &             -(XYZCEN(1,II)-CDGFBO(1,IFAC))*FEXTX(II)
     &             -(XYZCEN(2,II)-CDGFBO(2,IFAC))*FEXTY(II)
     &             -(XYZCEN(3,II)-CDGFBO(3,IFAC))*FEXTZ(II) )
              BX(II) = BX(II) +PFAC*SURFBO(1,IFAC)
              BY(II) = BY(II) +PFAC*SURFBO(2,IFAC)
              BZ(II) = BZ(II) +PFAC*SURFBO(3,IFAC)
            ENDDO
C
          ENDIF
C
        ENDIF
C
C     DPDX,DPDY,DPDZ = GRADIENT
C
        DO IEL = 1, NCEL
          UNSVOL = 1.D0/VOLUME(IEL)
          DPDX(IEL) = BX(IEL)*UNSVOL
          DPDY(IEL) = BY(IEL)*UNSVOL
          DPDZ(IEL) = BZ(IEL)*UNSVOL
        ENDDO
C
C     TRAITEMENT DU PARALLELISME
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (DPDX)
C         ===========
          CALL PARCOM (DPDY)
C         ===========
          CALL PARCOM (DPDZ)
C         ===========
        ENDIF
C
C     TRAITEMENT DE LA PERIODICITE
C
        IF(IPERIO.EQ.1) THEN
          CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    DPDX   , DPDX   , DPDX  ,
     &    DPDY   , DPDY   , DPDY  ,
     &    DPDZ   , DPDZ   , DPDZ  )
        ENDIF
C
        RETURN
C
      ENDIF
C
C
C=======================================================================
C 2. CONSTRUCTION DES GRADIENTS PAR UNE METHODE DE MOINDRE
C         CARRE POUR LES MAILLAGES TORDUS
C=======================================================================
C
C
      IF( (INICOC.EQ.1.OR.IALE.EQ.1) .AND.ICCOCG.EQ.1) THEN
C
C
C --->  2.1 CALCUL COMPLET DE COCG ET
C           SAUVEGARDE DE LA CONTRIBUTION AUX CELLULES DE BORD
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
C
C   ASSEMBLAGE A PARTIR DES FACETTES FLUIDES
C
        IF (IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1, NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            UNSDIJ = 1.D0/SQRT( (XYZCEN(1,II)-XYZCEN(1,JJ))**2
     &                         +(XYZCEN(2,II)-XYZCEN(2,JJ))**2
     &                         +(XYZCEN(3,II)-XYZCEN(3,JJ))**2 )
            DSIJ(1) = (XYZCEN(1,JJ)-XYZCEN(1,II))*UNSDIJ
            DSIJ(2) = (XYZCEN(2,JJ)-XYZCEN(2,II))*UNSDIJ
            DSIJ(3) = (XYZCEN(3,JJ)-XYZCEN(3,II))*UNSDIJ
C
            COCG(II,1,1) = COCG(II,1,1) +DSIJ(1)*DSIJ(1)
            COCG(II,2,2) = COCG(II,2,2) +DSIJ(2)*DSIJ(2)
            COCG(II,3,3) = COCG(II,3,3) +DSIJ(3)*DSIJ(3)
            COCG(II,1,2) = COCG(II,1,2) +DSIJ(1)*DSIJ(2)
            COCG(II,1,3) = COCG(II,1,3) +DSIJ(1)*DSIJ(3)
            COCG(II,2,3) = COCG(II,2,3) +DSIJ(2)*DSIJ(3)
C
            COCG(JJ,1,1) = COCG(JJ,1,1) +DSIJ(1)*DSIJ(1)
            COCG(JJ,2,2) = COCG(JJ,2,2) +DSIJ(2)*DSIJ(2)
            COCG(JJ,3,3) = COCG(JJ,3,3) +DSIJ(3)*DSIJ(3)
            COCG(JJ,1,2) = COCG(JJ,1,2) +DSIJ(1)*DSIJ(2)
            COCG(JJ,1,3) = COCG(JJ,1,3) +DSIJ(1)*DSIJ(3)
            COCG(JJ,2,3) = COCG(JJ,2,3) +DSIJ(2)*DSIJ(3)
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1, NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            UNSDIJ = 1.D0/SQRT( (XYZCEN(1,II)-XYZCEN(1,JJ))**2
     &                         +(XYZCEN(2,II)-XYZCEN(2,JJ))**2
     &                         +(XYZCEN(3,II)-XYZCEN(3,JJ))**2 )
            DSIJ(1) = (XYZCEN(1,JJ)-XYZCEN(1,II))*UNSDIJ
            DSIJ(2) = (XYZCEN(2,JJ)-XYZCEN(2,II))*UNSDIJ
            DSIJ(3) = (XYZCEN(3,JJ)-XYZCEN(3,II))*UNSDIJ
C
            COCG(II,1,1) = COCG(II,1,1) +DSIJ(1)*DSIJ(1)
            COCG(II,2,2) = COCG(II,2,2) +DSIJ(2)*DSIJ(2)
            COCG(II,3,3) = COCG(II,3,3) +DSIJ(3)*DSIJ(3)
            COCG(II,1,2) = COCG(II,1,2) +DSIJ(1)*DSIJ(2)
            COCG(II,1,3) = COCG(II,1,3) +DSIJ(1)*DSIJ(3)
            COCG(II,2,3) = COCG(II,2,3) +DSIJ(2)*DSIJ(3)
C
            COCG(JJ,1,1) = COCG(JJ,1,1) +DSIJ(1)*DSIJ(1)
            COCG(JJ,2,2) = COCG(JJ,2,2) +DSIJ(2)*DSIJ(2)
            COCG(JJ,3,3) = COCG(JJ,3,3) +DSIJ(3)*DSIJ(3)
            COCG(JJ,1,2) = COCG(JJ,1,2) +DSIJ(1)*DSIJ(2)
            COCG(JJ,1,3) = COCG(JJ,1,3) +DSIJ(1)*DSIJ(3)
            COCG(JJ,2,3) = COCG(JJ,2,3) +DSIJ(2)*DSIJ(3)
          ENDDO
C
        ENDIF
C
C ET COMPLEMENT POUR LE VOISINAGE ETENDU
C     PAS DE VECTORISATION PARTICULIERE A IMPOSER A PRIORI
C
        IF (IMRGRA.EQ.2.OR.IMRGRA.EQ.3) THEN
C
          DO II = 1, NCEL
            DO IPOS = IPCVSE(II), IPCVSE(II+1)-1
              JJ = IELVSE(IPOS)
C
              UNSDIJ = 1.D0/SQRT( (XYZCEN(1,II)-XYZCEN(1,JJ))**2
     &                           +(XYZCEN(2,II)-XYZCEN(2,JJ))**2
     &                           +(XYZCEN(3,II)-XYZCEN(3,JJ))**2 )
              DSIJ(1) = (XYZCEN(1,JJ)-XYZCEN(1,II))*UNSDIJ
              DSIJ(2) = (XYZCEN(2,JJ)-XYZCEN(2,II))*UNSDIJ
              DSIJ(3) = (XYZCEN(3,JJ)-XYZCEN(3,II))*UNSDIJ
C
              COCG(II,1,1) = COCG(II,1,1) +DSIJ(1)*DSIJ(1)
              COCG(II,2,2) = COCG(II,2,2) +DSIJ(2)*DSIJ(2)
              COCG(II,3,3) = COCG(II,3,3) +DSIJ(3)*DSIJ(3)
              COCG(II,1,2) = COCG(II,1,2) +DSIJ(1)*DSIJ(2)
              COCG(II,1,3) = COCG(II,1,3) +DSIJ(1)*DSIJ(3)
              COCG(II,2,3) = COCG(II,2,3) +DSIJ(2)*DSIJ(3)
C
            ENDDO
          ENDDO
C
        ENDIF
C
C  SAUVEGARDE DE LA CONTRIBUTION DES FACES INTERNES AUX ELEMENTS DE BORD
C
!OCL NOVREC,VRL(16)
        DO II = 1, NCELBR
          IEL = ICELBR(II)
          COCGB(II,1,1) = COCG(IEL,1,1)
          COCGB(II,1,2) = COCG(IEL,1,2)
          COCGB(II,1,3) = COCG(IEL,1,3)
          COCGB(II,2,1) = COCG(IEL,1,2)
          COCGB(II,2,2) = COCG(IEL,2,2)
          COCGB(II,2,3) = COCG(IEL,2,3)
          COCGB(II,3,1) = COCG(IEL,1,3)
          COCGB(II,3,2) = COCG(IEL,2,3)
          COCGB(II,3,3) = COCG(IEL,3,3)
        ENDDO
C
        INICOC = 0
C
C   ASSEMBLAGE A PARTIR DES FACETTES DE BORD
C
        IF (IVECTB.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1, NFABOR
C
            II = IFABOR(IFAC)
C
            EXTRAB = 1.D0-ISYMPA(IFAC)*EXTRAP*COEFBP(IFAC)
            UMCBSD = EXTRAB*(1.D0-COEFBP(IFAC))/DISTBR(IFAC)
            UNSSBN = EXTRAB/SURFBN(IFAC)
C
            DSIJ(1) = SURFBO(1,IFAC)*UNSSBN +UMCBSD*DIIPB(1,IFAC)
            DSIJ(2) = SURFBO(2,IFAC)*UNSSBN +UMCBSD*DIIPB(2,IFAC)
            DSIJ(3) = SURFBO(3,IFAC)*UNSSBN +UMCBSD*DIIPB(3,IFAC)
C
            COCG(II,1,1) = COCG(II,1,1) +DSIJ(1)*DSIJ(1)
            COCG(II,2,2) = COCG(II,2,2) +DSIJ(2)*DSIJ(2)
            COCG(II,3,3) = COCG(II,3,3) +DSIJ(3)*DSIJ(3)
            COCG(II,1,2) = COCG(II,1,2) +DSIJ(1)*DSIJ(2)
            COCG(II,1,3) = COCG(II,1,3) +DSIJ(1)*DSIJ(3)
            COCG(II,2,3) = COCG(II,2,3) +DSIJ(2)*DSIJ(3)
C
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1, NFABOR
C
            II = IFABOR(IFAC)
C
            EXTRAB = 1.D0-ISYMPA(IFAC)*EXTRAP*COEFBP(IFAC)
            UMCBSD = EXTRAB*(1.D0-COEFBP(IFAC))/DISTBR(IFAC)
            UNSSBN = EXTRAB/SURFBN(IFAC)
C
            DSIJ(1) = SURFBO(1,IFAC)*UNSSBN +UMCBSD*DIIPB(1,IFAC)
            DSIJ(2) = SURFBO(2,IFAC)*UNSSBN +UMCBSD*DIIPB(2,IFAC)
            DSIJ(3) = SURFBO(3,IFAC)*UNSSBN +UMCBSD*DIIPB(3,IFAC)
C
            COCG(II,1,1) = COCG(II,1,1) +DSIJ(1)*DSIJ(1)
            COCG(II,2,2) = COCG(II,2,2) +DSIJ(2)*DSIJ(2)
            COCG(II,3,3) = COCG(II,3,3) +DSIJ(3)*DSIJ(3)
            COCG(II,1,2) = COCG(II,1,2) +DSIJ(1)*DSIJ(2)
            COCG(II,1,3) = COCG(II,1,3) +DSIJ(1)*DSIJ(3)
            COCG(II,2,3) = COCG(II,2,3) +DSIJ(2)*DSIJ(3)
C
          ENDDO
C
        ENDIF
C
C   SYMETRISATION
C
        DO IEL = 1, NCEL
          COCG(IEL,2,1) = COCG(IEL,1,2)
          COCG(IEL,3,1) = COCG(IEL,1,3)
          COCG(IEL,3,2) = COCG(IEL,2,3)
        ENDDO
C
C   INVERSION
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
      ELSEIF(ICCOCG.EQ.1) THEN
C
C --->  2.2 CALCUL AUX BORDS DE COCG
C
        DO II = 1, 3
          DO JJ = 1, 3
!OCL NOVREC,VRL(16)
            DO IELB = 1, NCELBR
              IEL = ICELBR(IELB)
              COCG(IEL,II,JJ) = COCGB(IELB,II,JJ)
            ENDDO
          ENDDO
        ENDDO
C
        IF (IVECTB.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1, NFABOR
C
            II = IFABOR(IFAC)
C
            EXTRAB = 1.D0-ISYMPA(IFAC)*EXTRAP*COEFBP(IFAC)
            UMCBSD = EXTRAB*(1.D0-COEFBP(IFAC))/DISTBR(IFAC)
            UNSSBN = EXTRAB/SURFBN(IFAC)
C
            DSIJ(1) = SURFBO(1,IFAC)*UNSSBN+UMCBSD*DIIPB(1,IFAC)
            DSIJ(2) = SURFBO(2,IFAC)*UNSSBN+UMCBSD*DIIPB(2,IFAC)
            DSIJ(3) = SURFBO(3,IFAC)*UNSSBN+UMCBSD*DIIPB(3,IFAC)
C
            COCG(II,1,1) = COCG(II,1,1) +DSIJ(1)*DSIJ(1)
            COCG(II,2,2) = COCG(II,2,2) +DSIJ(2)*DSIJ(2)
            COCG(II,3,3) = COCG(II,3,3) +DSIJ(3)*DSIJ(3)
            COCG(II,1,2) = COCG(II,1,2) +DSIJ(1)*DSIJ(2)
            COCG(II,1,3) = COCG(II,1,3) +DSIJ(1)*DSIJ(3)
            COCG(II,2,3) = COCG(II,2,3) +DSIJ(2)*DSIJ(3)
C
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1, NFABOR
C
            II = IFABOR(IFAC)
C
            EXTRAB = 1.D0-ISYMPA(IFAC)*EXTRAP*COEFBP(IFAC)
            UMCBSD = EXTRAB*(1.D0-COEFBP(IFAC))/DISTBR(IFAC)
            UNSSBN = EXTRAB/SURFBN(IFAC)
C
            DSIJ(1) = SURFBO(1,IFAC)*UNSSBN+UMCBSD*DIIPB(1,IFAC)
            DSIJ(2) = SURFBO(2,IFAC)*UNSSBN+UMCBSD*DIIPB(2,IFAC)
            DSIJ(3) = SURFBO(3,IFAC)*UNSSBN+UMCBSD*DIIPB(3,IFAC)
C
            COCG(II,1,1) = COCG(II,1,1) +DSIJ(1)*DSIJ(1)
            COCG(II,2,2) = COCG(II,2,2) +DSIJ(2)*DSIJ(2)
            COCG(II,3,3) = COCG(II,3,3) +DSIJ(3)*DSIJ(3)
            COCG(II,1,2) = COCG(II,1,2) +DSIJ(1)*DSIJ(2)
            COCG(II,1,3) = COCG(II,1,3) +DSIJ(1)*DSIJ(3)
            COCG(II,2,3) = COCG(II,2,3) +DSIJ(2)*DSIJ(3)
C
          ENDDO
        ENDIF
C
C     SYMETRISATION
C
!OCL NOVREC,VRL(16)
        DO IELB = 1, NCELBR
          IEL = ICELBR(IELB)
          COCG(IEL,2,1) = COCG(IEL,1,2)
          COCG(IEL,3,1) = COCG(IEL,1,3)
          COCG(IEL,3,2) = COCG(IEL,2,3)
        ENDDO
C
C
        NBLOC = NCELBR/LBLOC
        IF (NBLOC.GT.0) THEN
          DO IBLOC = 1, NBLOC
            DO II = 1, LBLOC
              IELB = (IBLOC-1)*LBLOC+II
              IEL = ICELBR(IELB)
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
              AA(II,1,1) = A11*UNSDET
              AA(II,1,2) = A12*UNSDET
              AA(II,1,3) = A13*UNSDET
              AA(II,2,2) = A22*UNSDET
              AA(II,2,3) = A23*UNSDET
              AA(II,3,3) = A33*UNSDET
C
            ENDDO
C
!OCL NOVREC,VRL(16)
            DO II = 1, LBLOC
              IELB = (IBLOC-1)*LBLOC + II
              IEL = ICELBR(IELB)
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
        IREL  = MOD(NCELBR,LBLOC)
        IF (IREL .GT.0) THEN
          IBLOC = NBLOC + 1
          DO II = 1, IREL
            IELB = (IBLOC-1)*LBLOC+II
            IEL = ICELBR(IELB)
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
            AA(II,1,1) = A11*UNSDET
            AA(II,1,2) = A12*UNSDET
            AA(II,1,3) = A13*UNSDET
            AA(II,2,2) = A22*UNSDET
            AA(II,2,3) = A23*UNSDET
            AA(II,3,3) = A33*UNSDET
C
          ENDDO
C
!OCL NOVREC,VRL(16)
          DO II = 1, IREL
            IELB = (IBLOC-1)*LBLOC + II
            IEL = ICELBR(IELB)
            COCG(IEL,1,1) = AA(II,1,1)
            COCG(IEL,1,2) = AA(II,1,2)
            COCG(IEL,1,3) = AA(II,1,3)
            COCG(IEL,2,2) = AA(II,2,2)
            COCG(IEL,2,3) = AA(II,2,3)
            COCG(IEL,3,3) = AA(II,3,3)
          ENDDO
C
        ENDIF
C
C
C     SYMETRISATION DE LA MATRICE INVERSEE
C
!OCL NOVREC,VRL(16)
        DO IELB = 1, NCELBR
          IEL = ICELBR(IELB)
          COCG(IEL,2,1) = COCG(IEL,1,2)
          COCG(IEL,3,1) = COCG(IEL,1,3)
          COCG(IEL,3,2) = COCG(IEL,2,3)
        ENDDO
C
      ENDIF
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
C  CAS STANDARD, SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
C  ===============================================================
      IF (IPHYDP.EQ.0) THEN
C
C     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES
C
        IF (IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1,NFAC
C
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
C
            USDIJ2 = 1.D0/( ( XYZCEN(1,JJ)-XYZCEN(1,II) )**2
     &                     +( XYZCEN(2,JJ)-XYZCEN(2,II) )**2
     &                     +( XYZCEN(3,JJ)-XYZCEN(3,II) )**2 )
            VECFAC = ( PVAR(JJ)-PVAR(II) )*USDIJ2
C
            PFSX = ( XYZCEN(1,JJ)-XYZCEN(1,II) )*VECFAC
            PFSY = ( XYZCEN(2,JJ)-XYZCEN(2,II) )*VECFAC
            PFSZ = ( XYZCEN(3,JJ)-XYZCEN(3,II) )*VECFAC
            BX(II) = BX(II) +PFSX
            BY(II) = BY(II) +PFSY
            BZ(II) = BZ(II) +PFSZ
            BX(JJ) = BX(JJ) +PFSX
            BY(JJ) = BY(JJ) +PFSY
            BZ(JJ) = BZ(JJ) +PFSZ
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1,NFAC
C
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
C
            USDIJ2 = 1.D0/( ( XYZCEN(1,JJ)-XYZCEN(1,II) )**2
     &                     +( XYZCEN(2,JJ)-XYZCEN(2,II) )**2
     &                     +( XYZCEN(3,JJ)-XYZCEN(3,II) )**2 )
            VECFAC = ( PVAR(JJ)-PVAR(II) )*USDIJ2
C
            PFSX = ( XYZCEN(1,JJ)-XYZCEN(1,II) )*VECFAC
            PFSY = ( XYZCEN(2,JJ)-XYZCEN(2,II) )*VECFAC
            PFSZ = ( XYZCEN(3,JJ)-XYZCEN(3,II) )*VECFAC
            BX(II) = BX(II) +PFSX
            BY(II) = BY(II) +PFSY
            BZ(II) = BZ(II) +PFSZ
            BX(JJ) = BX(JJ) +PFSX
            BY(JJ) = BY(JJ) +PFSY
            BZ(JJ) = BZ(JJ) +PFSZ
          ENDDO
C
        ENDIF
C
C ET COMPLEMENT POUR LE VOISINAGE ETENDU
C     LA REECRITURE COMPLEXE DE LA BOUCLE EST DESTINEE A DIVISER LE CPU
C       TOTAL SUR VPP PAR 5. POUR CELA, ON UTILISE DPDX DPDY DPDZ COMME
C       TABLEAUX DE TRAVAIL LOCAUX
C     SI ON TOMBAIT SUR UN CAS PATHOLOGIQUES DANS LEQUEL NCEL < NOMBRE
C       DE CELLULES DU VOISINAGE ETENDU D'UNE CELLULE (AVEC PERIODICITE
C       MULTIPLE ET PEU DE MAILLES ???) ON FAIT UN TEST ET ON CONSERVE
C       LA BOUCLE INITIALE. ON AJOUTE CEPENDANT UN IF QUI PERMET SUR VPP
C       DE REDUIRE LE TEMPS CALCUL D'UN FACTEUR 4 SI LA BOUCLE N'EST PAS
C       SCINDEE EN DEUX.
C
        IF (IMRGRA.EQ.2.OR.IMRGRA.EQ.3) THEN
C
C     ON REGARDE SI NCEL >= NBRE DE VOISINS ETENDUS D'UNE CELLULE
          IF (ICESVE.LT.0) THEN
            INVEMX = 0
            DO II = 1, NCEL
              INVEMX = MAX(INVEMX,IPCVSE(II+1)-IPCVSE(II))
            ENDDO
            IF (NCEL.GE.INVEMX) THEN
              ICESVE = 1
            ELSE
              ICESVE = 0
            ENDIF
          ENDIF
C
          IF (ICESVE.GT.0) THEN
C
            DO II = 1, NCEL
              DO IPOS = IPCVSE(II), IPCVSE(II+1)-1
                JJ = IELVSE(IPOS)
                USDIJ2 = 1.D0/( ( XYZCEN(1,JJ)-XYZCEN(1,II) )**2
     &                         +( XYZCEN(2,JJ)-XYZCEN(2,II) )**2
     &                         +( XYZCEN(3,JJ)-XYZCEN(3,II) )**2 )
                VECFAC = ( PVAR(JJ)-PVAR(II) )*USDIJ2
C
                IIII = IPOS-IPCVSE(II)+1
                DPDX(IIII) = ( XYZCEN(1,JJ)-XYZCEN(1,II) )*VECFAC
                DPDY(IIII) = ( XYZCEN(2,JJ)-XYZCEN(2,II) )*VECFAC
                DPDZ(IIII) = ( XYZCEN(3,JJ)-XYZCEN(3,II) )*VECFAC
              ENDDO
              DO IPOS = IPCVSE(II), IPCVSE(II+1)-1
                IIII = IPOS-IPCVSE(II)+1
                BX(II) = BX(II) +DPDX(IIII)
                BY(II) = BY(II) +DPDY(IIII)
                BZ(II) = BZ(II) +DPDZ(IIII)
              ENDDO
            ENDDO
C
          ELSE
C
            DO II = 1, NCEL
              IF(IPCVSE(II+1).GT.IPCVSE(II)) THEN
                DO IPOS = IPCVSE(II), IPCVSE(II+1)-1
                  JJ = IELVSE(IPOS)
                  USDIJ2 = 1.D0/( ( XYZCEN(1,JJ)-XYZCEN(1,II) )**2
     &                           +( XYZCEN(2,JJ)-XYZCEN(2,II) )**2
     &                           +( XYZCEN(3,JJ)-XYZCEN(3,II) )**2 )
                  VECFAC = ( PVAR(JJ)-PVAR(II) )*USDIJ2
C
                  PFSX = ( XYZCEN(1,JJ)-XYZCEN(1,II) )*VECFAC
                  PFSY = ( XYZCEN(2,JJ)-XYZCEN(2,II) )*VECFAC
                  PFSZ = ( XYZCEN(3,JJ)-XYZCEN(3,II) )*VECFAC
                  BX(II) = BX(II) +PFSX
                  BY(II) = BY(II) +PFSY
                  BZ(II) = BZ(II) +PFSZ
                ENDDO
              ENDIF
            ENDDO
C
          ENDIF
C
        ENDIF
C
C
C     ASSEMBLAGE A PARTIR DES FACETTES DE BORD
C
        IF (IVECTB.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1,NFABOR
C
            II = IFABOR(IFAC)
C
            EXTRAB = (1.D0-ISYMPA(IFAC)*EXTRAP*COEFBP(IFAC))**2
            UNSDIJ = 1.D0/DISTBR(IFAC)
            UNSSBN = 1.D0/SURFBN(IFAC)
            UMCBSD = (1.D0-COEFBP(IFAC))*UNSDIJ
C
            DSIJ(1) = SURFBO(1,IFAC)*UNSSBN + UMCBSD*DIIPB(1,IFAC)
            DSIJ(2) = SURFBO(2,IFAC)*UNSSBN + UMCBSD*DIIPB(2,IFAC)
            DSIJ(3) = SURFBO(3,IFAC)*UNSSBN + UMCBSD*DIIPB(3,IFAC)
            RKIJ = ( INC*COEFAP(IFAC)
     &               +(COEFBP(IFAC)-1.D0)*PVAR(II) )*UNSDIJ*EXTRAB
C
            BX(II) = BX(II) +DSIJ(1)*RKIJ
            BY(II) = BY(II) +DSIJ(2)*RKIJ
            BZ(II) = BZ(II) +DSIJ(3)*RKIJ
C
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1,NFABOR
C
            II = IFABOR(IFAC)
C
            EXTRAB = (1.D0-ISYMPA(IFAC)*EXTRAP*COEFBP(IFAC))**2
            UNSDIJ = 1.D0/DISTBR(IFAC)
            UNSSBN = 1.D0/SURFBN(IFAC)
            UMCBSD = (1.D0-COEFBP(IFAC))*UNSDIJ
C
            DSIJ(1) = SURFBO(1,IFAC)*UNSSBN + UMCBSD*DIIPB(1,IFAC)
            DSIJ(2) = SURFBO(2,IFAC)*UNSSBN + UMCBSD*DIIPB(2,IFAC)
            DSIJ(3) = SURFBO(3,IFAC)*UNSSBN + UMCBSD*DIIPB(3,IFAC)
            RKIJ = ( INC*COEFAP(IFAC)
     &               +(COEFBP(IFAC)-1.D0)*PVAR(II) )*UNSDIJ*EXTRAB
C
            BX(II) = BX(II) +DSIJ(1)*RKIJ
            BY(II) = BY(II) +DSIJ(2)*RKIJ
            BZ(II) = BZ(II) +DSIJ(3)*RKIJ
C
          ENDDO
C
        ENDIF
C
C  CAS AVEC PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
C  =====================================================
      ELSEIF(IPHYDP.NE.0) THEN
C
C     ASSEMBLAGE A PARTIR DES FACETTES FLUIDES
C
        IF (IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1,NFAC
C
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
C
            USDIJ2 = 1.D0/( ( XYZCEN(1,JJ)-XYZCEN(1,II) )**2
     &           +( XYZCEN(2,JJ)-XYZCEN(2,II) )**2
     &           +( XYZCEN(3,JJ)-XYZCEN(3,II) )**2 )
            VECFAC = ( PVAR(JJ)-PVAR(II)
     &           +(XYZCEN(1,II)-CDGFAC(1,IFAC))*FEXTX(II)
     &           +(XYZCEN(2,II)-CDGFAC(2,IFAC))*FEXTY(II)
     &           +(XYZCEN(3,II)-CDGFAC(3,IFAC))*FEXTZ(II)
     &           -(XYZCEN(1,JJ)-CDGFAC(1,IFAC))*FEXTX(JJ)
     &           -(XYZCEN(2,JJ)-CDGFAC(2,IFAC))*FEXTY(JJ)
     &           -(XYZCEN(3,JJ)-CDGFAC(3,IFAC))*FEXTZ(JJ) )*USDIJ2
C
            PFSX = ( XYZCEN(1,JJ)-XYZCEN(1,II) )*VECFAC
            PFSY = ( XYZCEN(2,JJ)-XYZCEN(2,II) )*VECFAC
            PFSZ = ( XYZCEN(3,JJ)-XYZCEN(3,II) )*VECFAC
            BX(II) = BX(II) +PFSX
            BY(II) = BY(II) +PFSY
            BZ(II) = BZ(II) +PFSZ
            BX(JJ) = BX(JJ) +PFSX
            BY(JJ) = BY(JJ) +PFSY
            BZ(JJ) = BZ(JJ) +PFSZ
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1,NFAC
C
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
C
            USDIJ2 = 1.D0/( ( XYZCEN(1,JJ)-XYZCEN(1,II) )**2
     &           +( XYZCEN(2,JJ)-XYZCEN(2,II) )**2
     &           +( XYZCEN(3,JJ)-XYZCEN(3,II) )**2 )
            VECFAC = ( PVAR(JJ)-PVAR(II)
     &           +(XYZCEN(1,II)-CDGFAC(1,IFAC))*FEXTX(II)
     &           +(XYZCEN(2,II)-CDGFAC(2,IFAC))*FEXTY(II)
     &           +(XYZCEN(3,II)-CDGFAC(3,IFAC))*FEXTZ(II)
     &           -(XYZCEN(1,JJ)-CDGFAC(1,IFAC))*FEXTX(JJ)
     &           -(XYZCEN(2,JJ)-CDGFAC(2,IFAC))*FEXTY(JJ)
     &           -(XYZCEN(3,JJ)-CDGFAC(3,IFAC))*FEXTZ(JJ) )*USDIJ2
C
            PFSX = ( XYZCEN(1,JJ)-XYZCEN(1,II) )*VECFAC
            PFSY = ( XYZCEN(2,JJ)-XYZCEN(2,II) )*VECFAC
            PFSZ = ( XYZCEN(3,JJ)-XYZCEN(3,II) )*VECFAC
            BX(II) = BX(II) +PFSX
            BY(II) = BY(II) +PFSY
            BZ(II) = BZ(II) +PFSZ
            BX(JJ) = BX(JJ) +PFSX
            BY(JJ) = BY(JJ) +PFSY
            BZ(JJ) = BZ(JJ) +PFSZ
          ENDDO
C
        ENDIF
C
C
C ET COMPLEMENT POUR LE VOISINAGE ETENDU
C     ON SUPPOSE QUE LE MILIEU DU SEGMENT JOIGNANT LES CENTRES
C     PEUT REMPLACER LE CENTRE DE GRAVITE DE LA FACE FICTIVE
C
C     LA REECRITURE COMPLEXE DE LA BOUCLE EST DESTINEE A DIVISER LE CPU
C       TOTAL SUR VPP PAR 5. POUR CELA, ON UTILISE DPDX DPDY DPDZ COMME
C       TABLEAUX DE TRAVAIL LOCAUX
C     SI ON TOMBAIT SUR UN CAS PATHOLOGIQUES DANS LEQUEL NCEL < NOMBRE
C       DE CELLULES DU VOISINAGE ETENDU D'UNE CELLULE (AVEC PERIODICITE
C       MULTIPLE ET PEU DE MAILLES ???) ON FAIT UN TEST ET ON CONSERVE
C       LA BOUCLE INITIALE. ON AJOUTE CEPENDAT UN IF QUI PERMET SUR VPP
C       DE REDUIRE LE TEMPS CALCUL D'UN FACTEUR 4 SI LA BOUCLE N'EST PAS
C       SCINDEE EN DEUX.
C
        IF (IMRGRA.EQ.2.OR.IMRGRA.EQ.3) THEN
C
C     ON REGARDE SI NCEL >= NBRE DE VOISINS ETENDUS D'UNE CELLULE
          IF (ICESVE.LT.0) THEN
            INVEMX = 0
            DO II = 1, NCEL
              INVEMX = MAX(INVEMX,IPCVSE(II+1)-IPCVSE(II))
            ENDDO
            IF (NCEL.GE.INVEMX) THEN
              ICESVE = 1
            ELSE
              ICESVE = 0
            ENDIF
          ENDIF
C
          IF (ICESVE.GT.0) THEN
C
            DO II = 1, NCEL
              DO IPOS = IPCVSE(II), IPCVSE(II+1)-1
                JJ = IELVSE(IPOS)
                USDIJ2 = 1.D0/( ( XYZCEN(1,JJ)-XYZCEN(1,II) )**2
     &                         +( XYZCEN(2,JJ)-XYZCEN(2,II) )**2
     &                         +( XYZCEN(3,JJ)-XYZCEN(3,II) )**2 )
                PTMIDX = 0.5D0*(XYZCEN(1,JJ)+XYZCEN(1,II))
                PTMIDY = 0.5D0*(XYZCEN(2,JJ)+XYZCEN(2,II))
                PTMIDZ = 0.5D0*(XYZCEN(3,JJ)+XYZCEN(3,II))
                VECFAC = ( PVAR(JJ)-PVAR(II)
     &               +(XYZCEN(1,II)-PTMIDX)*FEXTX(II)
     &               +(XYZCEN(2,II)-PTMIDY)*FEXTY(II)
     &               +(XYZCEN(3,II)-PTMIDZ)*FEXTZ(II)
     &               -(XYZCEN(1,JJ)-PTMIDX)*FEXTX(JJ)
     &               -(XYZCEN(2,JJ)-PTMIDY)*FEXTY(JJ)
     &               -(XYZCEN(3,JJ)-PTMIDZ)*FEXTZ(JJ) )*USDIJ2
C
                IIII = IPOS-IPCVSE(II)+1
                DPDX(IIII) = ( XYZCEN(1,JJ)-XYZCEN(1,II) )*VECFAC
                DPDY(IIII) = ( XYZCEN(2,JJ)-XYZCEN(2,II) )*VECFAC
                DPDZ(IIII) = ( XYZCEN(3,JJ)-XYZCEN(3,II) )*VECFAC
              ENDDO
              DO IPOS = IPCVSE(II), IPCVSE(II+1)-1
                IIII = IPOS-IPCVSE(II)+1
                BX(II) = BX(II) +DPDX(IIII)
                BY(II) = BY(II) +DPDY(IIII)
                BZ(II) = BZ(II) +DPDZ(IIII)
              ENDDO
            ENDDO
C
          ELSE
C
            DO II = 1, NCEL
              IF(IPCVSE(II+1).GT.IPCVSE(II)) THEN
                DO IPOS = IPCVSE(II), IPCVSE(II+1)-1
                  JJ = IELVSE(IPOS)
                  USDIJ2 = 1.D0/( ( XYZCEN(1,JJ)-XYZCEN(1,II) )**2
     &                           +( XYZCEN(2,JJ)-XYZCEN(2,II) )**2
     &                           +( XYZCEN(3,JJ)-XYZCEN(3,II) )**2 )
                  PTMIDX = 0.5D0*(XYZCEN(1,JJ)+XYZCEN(1,II))
                  PTMIDY = 0.5D0*(XYZCEN(2,JJ)+XYZCEN(2,II))
                  PTMIDZ = 0.5D0*(XYZCEN(3,JJ)+XYZCEN(3,II))
                  VECFAC = ( PVAR(JJ)-PVAR(II)
     &                 +(XYZCEN(1,II)-PTMIDX)*FEXTX(II)
     &                 +(XYZCEN(2,II)-PTMIDY)*FEXTY(II)
     &                 +(XYZCEN(3,II)-PTMIDZ)*FEXTZ(II)
     &                 -(XYZCEN(1,JJ)-PTMIDX)*FEXTX(JJ)
     &                 -(XYZCEN(2,JJ)-PTMIDY)*FEXTY(JJ)
     &                 -(XYZCEN(3,JJ)-PTMIDZ)*FEXTZ(JJ) )*USDIJ2
C
                  PFSX = ( XYZCEN(1,JJ)-XYZCEN(1,II) )*VECFAC
                  PFSY = ( XYZCEN(2,JJ)-XYZCEN(2,II) )*VECFAC
                  PFSZ = ( XYZCEN(3,JJ)-XYZCEN(3,II) )*VECFAC
                  BX(II) = BX(II) +PFSX
                  BY(II) = BY(II) +PFSY
                  BZ(II) = BZ(II) +PFSZ
                ENDDO
              ENDIF
            ENDDO
C
          ENDIF
C
        ENDIF
C
C
C     ASSEMBLAGE A PARTIR DES FACETTES DE BORD
C      NOTER QUE EXTRAB NE DOIT PRENDRE QUE LES VALEURS 0 OU 1
C
        IF (IVECTB.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1,NFABOR
C
            II = IFABOR(IFAC)
C
            EXTRAB = (1.D0-ISYMPA(IFAC)*EXTRAP*COEFBP(IFAC))**2
            UNSDIJ = 1.D0/DISTBR(IFAC)
            UNSSBN = 1.D0/SURFBN(IFAC)
            UMCBSD = (1.D0-COEFBP(IFAC))*UNSDIJ
C
            DSIJ(1) = SURFBO(1,IFAC)*UNSSBN + UMCBSD*DIIPB(1,IFAC)
            DSIJ(2) = SURFBO(2,IFAC)*UNSSBN + UMCBSD*DIIPB(2,IFAC)
            DSIJ(3) = SURFBO(3,IFAC)*UNSSBN + UMCBSD*DIIPB(3,IFAC)
            RKIJ = ( INC*COEFAP(IFAC)
     &           +(COEFBP(IFAC)-1.D0)*PVAR(II) )*UNSDIJ
     &           +(COEFBP(IFAC)-1.D0)*UNSDIJ*(
     &           (CDGFBO(1,IFAC)-XYZCEN(1,II))*FEXTX(II)
     &           +(CDGFBO(2,IFAC)-XYZCEN(2,II))*FEXTY(II)
     &           +(CDGFBO(3,IFAC)-XYZCEN(3,II))*FEXTZ(II) )
            RKIJ = RKIJ*EXTRAB
C
            BX(II) = BX(II) +DSIJ(1)*RKIJ
            BY(II) = BY(II) +DSIJ(2)*RKIJ
            BZ(II) = BZ(II) +DSIJ(3)*RKIJ
C
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1,NFABOR
C
            II = IFABOR(IFAC)
C
            EXTRAB = (1.D0-ISYMPA(IFAC)*EXTRAP*COEFBP(IFAC))**2
            UNSDIJ = 1.D0/DISTBR(IFAC)
            UNSSBN = 1.D0/SURFBN(IFAC)
            UMCBSD = (1.D0-COEFBP(IFAC))*UNSDIJ
C
            DSIJ(1) = SURFBO(1,IFAC)*UNSSBN + UMCBSD*DIIPB(1,IFAC)
            DSIJ(2) = SURFBO(2,IFAC)*UNSSBN + UMCBSD*DIIPB(2,IFAC)
            DSIJ(3) = SURFBO(3,IFAC)*UNSSBN + UMCBSD*DIIPB(3,IFAC)
            RKIJ = ( INC*COEFAP(IFAC)
     &           +(COEFBP(IFAC)-1.D0)*PVAR(II) )*UNSDIJ
     &           +(COEFBP(IFAC)-1.D0)*UNSDIJ*(
     &           (CDGFBO(1,IFAC)-XYZCEN(1,II))*FEXTX(II)
     &           +(CDGFBO(2,IFAC)-XYZCEN(2,II))*FEXTY(II)
     &           +(CDGFBO(3,IFAC)-XYZCEN(3,II))*FEXTZ(II) )
            RKIJ = RKIJ*EXTRAB
C
            BX(II) = BX(II) +DSIJ(1)*RKIJ
            BY(II) = BY(II) +DSIJ(2)*RKIJ
            BZ(II) = BZ(II) +DSIJ(3)*RKIJ
C
          ENDDO
C
        ENDIF
C
      ENDIF
C
C=======================================================================
C 4. RESOLUTION
C=======================================================================
C
C
      DO IEL = 1, NCEL
        DPDX(IEL) = COCG(IEL,1,1)*BX(IEL)+COCG(IEL,1,2)*BY(IEL)
     &             +COCG(IEL,1,3)*BZ(IEL)
        DPDY(IEL) = COCG(IEL,2,1)*BX(IEL)+COCG(IEL,2,2)*BY(IEL)
     &             +COCG(IEL,2,3)*BZ(IEL)
        DPDZ(IEL) = COCG(IEL,3,1)*BX(IEL)+COCG(IEL,3,2)*BY(IEL)
     &             +COCG(IEL,3,3)*BZ(IEL)
      ENDDO
C
C     AJOUT DE LA COMPOSANTE HYDROSTATIQUE
      IF (IPHYDP.EQ.1) THEN
        DO IEL = 1, NCEL
          DPDX(IEL) = DPDX(IEL) + FEXTX(IEL)
          DPDY(IEL) = DPDY(IEL) + FEXTY(IEL)
          DPDZ(IEL) = DPDZ(IEL) + FEXTZ(IEL)
        ENDDO
      ENDIF
C
C     TRAITEMENT DU PARALLELISME
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (DPDX)
C         ===========
          CALL PARCOM (DPDY)
C         ===========
          CALL PARCOM (DPDZ)
C         ===========
        ENDIF
C
C     TRAITEMENT DE LA PERIODICITE
C
        IF(IPERIO.EQ.1) THEN
          CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    DPDX   , DPDX   , DPDX  ,
     &    DPDY   , DPDY   , DPDY  ,
     &    DPDZ   , DPDZ   , DPDZ  )
        ENDIF
C
C
C--------
C FORMATS
C--------
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
