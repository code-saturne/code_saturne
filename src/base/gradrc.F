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
                        SUBROUTINE GRADRC
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NFAC   , NFABOR , NCELBR ,
     &   IMRGRA , INC    , ICCOCG , NSWRGP , IDIMTE , ITENSO , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , EXTRAP ,
     &   IFACEL , IFABOR , ICELBR , IVAR   ,
     &   VOLUME , SURFAC , SURFBO , POND   , XYZCEN , CDGFAC , CDGFBO ,
     &   DIJPF  , DIIPB  , DOFIJ  , FEXTX  , FEXTY  , FEXTZ  ,
     &   COEFAP , COEFBP , PVAR   ,
     &   COCGB  , COCG   ,
     &   DPDX   , DPDY   , DPDZ   ,
     &   BX     , BY     , BZ     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC CALCUL DU GRADIENT D'UNE VARIABLE AVEC UNE TECHNIQUE IERATIVE
CFONC DE RECONSTRUCTION POUR LES MAILLAGES NON ORTHOGONAUX (NSWRGP>1)
CFONC AVEC PRISE EN COMPTE EVENTUELLE D'UN TERME DE FORCE VOLUMIQUE
CFONC GENERANT UNE COMPOSANTE DE PRESSION HYDROSTATIQUE
CFONC
CFONC CALCUL DE COCG POUR PRENDRE EN COMPTE LES C.L VARIABLES (FLUX)
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
CARGU ! IMRGRA       ! E  !  -> ! METHODE DE CALCUL DU GRADIENT        !
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
CARGU ! IVAR         ! E  !  -> ! INDICATEUR DU NUMERO DE LA VARIABLE  !
CARGU !              !    !     ! (OU 0 SI VARIABLE NON RESOLUE)       !
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! IPHYDP       ! E  !  -> ! INDICATEUR DE PRISE EN COMPTE DE LA  !
CARGU !              !    !     ! PRESSION HYDROSTATIQUE               !
CARGU ! NFECRA       ! E  !  -> ! UNITE DU FICHIER SORTIE STD          !
CARGU ! EPSRGP       ! R  !  -> ! PRECISION RELATIVE POUR LA           !
CARGU !              !    !     !  RECONSTRUCTION DES GRADIENTS 97     !
CARGU ! EXTRAP       ! R  !  -> ! COEF EXTRAP GRADIENT                 !
CARGU ! IFACEL(2,NFAC! TE !  -> ! No DES ELTS VOISINS D'UNE FACE INTERN!
CARGU ! IFABOR(NFABOR! TE !  -> ! No DE L'ELT VOISIN D'UNE FACE DE BORD!
CARGU ! ICELBR       ! TE !  -> ! NUMERO GLOBAL DES ELEMENTS AYANT AU  !
CARGU ! (NCELBR)     !    !     !  MOINS UNE FACE DE BORD              !
CARGU ! VOLUME(NCELET! TR !  -> ! VOLUME DES ELEMENTS                  !
CARGU ! SURFAC(3,NFAC! TR !  -> ! SURF VECTORIELLE DES SURFACES INTERNE!
CARGU ! SURFBO       ! TR !  -> ! SURF VECTORIELLE DES SURFACES DE BORD!
CARGU !   (3,NFABOR) !    !     !                                      !
CARGU ! POND(NFAC)   ! TR !  -> ! PONDERATION GEOMETRIQUE (ENTRE 0 ET 1!
CARGU ! XYZCEN       ! TR !  -> ! POINT ASSOCIES AUX VOLUMES DE CONTROL!
CARGU ! (3,NCELET    !    !     !                                      !
CARGU ! CDGFAC       ! TR !  -> ! CENTRE DE GRAVITE DES FACES INTERNES !
CARGU ! (3,NFAC)     !    !     !                                      !
CARGU ! CDGFBO       ! TR !  -> ! CENTRE DE GRAVITE DES FACES DE BORD  !
CARGU ! (3,NFABOR)   !    !     !                                      !
CARGU ! DIJPF(3,NFAC)! TR !  -> ! VECT I'J', I' (RESP. J') PROJECTION  !
CARGU !              !    !     !  DU CENTRE I (RESP. J) SUR LA NORMALE!
CARGU !              !    !     !  A LA FACE INTERNE                   !
CARGU ! DIIPB        ! TR !  -> ! VECT II', II PROJECTION DU CENTRE I  !
CARGU !   (3,NFABOR) !    !     !  SUR LA NORMALE A LA FACE DE BORD    !
CARGU ! COEFAP, B    ! TR !  -> ! TABLEAUX DES COND LIM POUR PVAR      !
CARGU !   (NFABOR)   !    !     !  SUR LA NORMALE A LA FACE DE BORD    !
CARGU ! PVAR  (NCELET! TR !  -> ! VARIABLE                             !
CARGU ! FEXTX,Y,Z    ! TR !  -> ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU !   (NCELET)   !    !     !  HYDROSTATIQUE                       !
CARGU ! COCGB        ! TR !  -> ! CONTRIBUTION DES FACES INTERNES A    !
CARGU !  (NCELBR,3,3)!    !     !  COCG POUR LES CELLULES DE BORD      !
CARGU ! COCG(NCELET,3! TR ! <-> ! COUPLAGE DES COMPOSANTES DU GRADIENT !
CARGU !     ,3)      !    !     ! LORS DE LA RECONSTRUCTION            !
CARGU ! DPDX DPDY    ! TR ! <-> ! GRADIENT DE PVAR                     !
CARGU ! DPDZ (NCELET ! TR !     !   (HALO REMPLI POUR LA PERIODICITE)  !
CARGU ! BX,Y,Z(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR LE GRAD DE P !
CARGU ! DOFIJ        ! TR ! <-  ! VECTEUR OF POUR LES FACES INTERNES   !
CARGU ! (NDIM,NFAC  )!    !     ! O : INTERSECTION DE IJ ET LA FACE    !
CARGU !              !    !     ! F : CENTRE DE LA FACE                !
CARGU ! RA(*)        ! TR !  -  ! TABLEAU DE TRAVAIL POUR LES REELS    !
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
      INCLUDE "cstnum.h"
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
      INTEGER          IMRGRA , INC    , ICCOCG , NSWRGP
      INTEGER          IVAR   , IDIMTE , ITENSO , IPHYDP
      INTEGER          IWARNP , NFECRA
      DOUBLE PRECISION EPSRGP , EXTRAP
C
      INTEGER          IFACEL(2,NFAC),IFABOR(NFABOR)
      INTEGER          ICELBR(NCELBR)
      DOUBLE PRECISION VOLUME(NCELET), SURFAC(3,NFAC)
      DOUBLE PRECISION SURFBO(3,NFABOR)
      DOUBLE PRECISION POND(NFAC)
      DOUBLE PRECISION XYZCEN(3,NCELET)
      DOUBLE PRECISION CDGFAC(3,NFAC),CDGFBO(3,NFABOR)
      DOUBLE PRECISION DIJPF(3,NFAC), DIIPB(3,NFABOR)
      DOUBLE PRECISION COEFAP(NFABOR), COEFBP(NFABOR), PVAR(NCELET)
      DOUBLE PRECISION COCGB(NCELBR,3,3), COCG(NCELET,3,3)
      DOUBLE PRECISION DPDX (NCELET),DPDY (NCELET),DPDZ (NCELET)
      DOUBLE PRECISION BX   (NCELET),BY   (NCELET),BZ   (NCELET)
      DOUBLE PRECISION FEXTX(NCELET),FEXTY(NCELET),FEXTZ(NCELET)
      DOUBLE PRECISION DOFIJ(3,NFAC)
C
C VARIABLES LOCALES
C
      INTEGER          LBLOC
      PARAMETER       (LBLOC = 1024)
      INTEGER          IEL, IFAC, II, JJ, KK, LL, MM
      INTEGER          ISQRT, ISWEEP, NSWMAX
      INTEGER          IBLOC,NBLOC,IREL
      DOUBLE PRECISION PFAC,PIP,DELTPX,DELTPY,DELTPZ
      DOUBLE PRECISION RNORX,RNORY,RNORZ,RNORM,RESIDU
      DOUBLE PRECISION PFACI
      DOUBLE PRECISION DOF,FMOYX,FMOYY,FMOYZ
      DOUBLE PRECISION AA(LBLOC,3,3),UNSDET,PFSX,PFSY,PFSZ
      DOUBLE PRECISION A11,A12,A13,A21,A22,A23,A31,A32,A33
      DOUBLE PRECISION COCG11,COCG12,COCG13,COCG21,COCG22,COCG23
      DOUBLE PRECISION COCG31,COCG32,COCG33
      DOUBLE PRECISION PFAC1, PFAC2, PFAC3, UNSVOL, VECFAC
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      ISQRT = 1
      NSWMAX = NSWRGP
      ISWEEP = 1
C
C
C=======================================================================
C 2. CALCUL SANS RECONSTRUCTION
C=======================================================================
C
C     SI INITIALISATION PAR MOINDRES CARRES (IMRGRA=4), B EST DEJA REMPLI
C     SINON (IMRGRA=0) ON CALCULE UN GRADIENT SANS RECONSTRUCTION
C
      IF (IMRGRA.EQ.0) THEN
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
              PFAC   = POND(IFAC)*PVAR(II) +(1.D0-POND(IFAC))*PVAR(JJ)
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
              PFAC   = POND(IFAC)*PVAR(II) +(1.D0-POND(IFAC))*PVAR(JJ)
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
     &    ( IDIMTE , ITENSO ,
     &      DPDX   , DPDX   , DPDX  ,
     &      DPDY   , DPDY   , DPDY  ,
     &      DPDZ   , DPDZ   , DPDZ  )
        ENDIF
C
      ENDIF
C
      IF( NSWRGP.LE.1 ) RETURN
C
C
C     On incremente IPASS quand on calcule COCG pour la premiere fois
      IPASS = IPASS + 1
C
C=======================================================================
C 3. RECONSTRUCTION DES GRADIENTS POUR LES MAILLAGES TORDUS
C=======================================================================
C
C     RESOLUTION SEMI-IMPLICITE SUR TOUT LE MAILLAGE
C     DPDX,DY,DZ = GRADIENT
C
      IF(IPASS.EQ.1 .OR. IALE.EQ.1) THEN
C
C ---> CALCUL DE COCG
C
        DO II = 1, 3
          DO JJ = 1, 3
            DO IEL =1,NCELET
              COCG(IEL,II,JJ) = 0.D0
            ENDDO
          ENDDO
        ENDDO
        DO IEL=1,NCEL
          COCG(IEL,1,1) = VOLUME(IEL)
          COCG(IEL,2,2) = VOLUME(IEL)
          COCG(IEL,3,3) = VOLUME(IEL)
        ENDDO
C
C ---> AJOUT DES CONTRIBUTIONS DES FACES INTERNES
        DO LL =1,3
          DO MM =1,3
C
            IF (IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
              DO IFAC=1,NFAC
                II = IFACEL(1,IFAC)
                JJ = IFACEL(2,IFAC)
C
C---> DOF = OF
                DOF = DOFIJ(MM,IFAC)
C
                PFACI  = -DOF*0.5D0
                VECFAC = PFACI*SURFAC(LL,IFAC)
                COCG(II,LL,MM) = COCG(II,LL,MM) + VECFAC
                COCG(JJ,LL,MM) = COCG(JJ,LL,MM) - VECFAC
              ENDDO
C
            ELSE
C
C VECTORISATION NON FORCEE
              DO IFAC=1,NFAC
                II = IFACEL(1,IFAC)
                JJ = IFACEL(2,IFAC)
C
C---> DOF = OF
                DOF = DOFIJ(MM,IFAC)
C
                PFACI  = -DOF*0.5D0
                VECFAC = PFACI*SURFAC(LL,IFAC)
                COCG(II,LL,MM) = COCG(II,LL,MM) + VECFAC
                COCG(JJ,LL,MM) = COCG(JJ,LL,MM) - VECFAC
              ENDDO
C
            ENDIF
C
          ENDDO
        ENDDO
C
C ---> SAUVEGADE DU COCG PARTIEL AUX FACES INTERNES DES CELLULES DE BORD
        DO II = 1, NCELBR
          IEL = ICELBR(II)
          DO LL = 1, 3
            DO MM = 1, 3
              COCGB(II,LL,MM) = COCG(IEL,LL,MM)
            ENDDO
          ENDDO
        ENDDO
C
C ---> AJOUT DES CONTRIBUTIONS DES FACES DE BORD
        DO LL =1,3
          DO MM =1,3
C
            IF (IVECTB.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
              DO IFAC=1,NFABOR
                II = IFABOR(IFAC)
                COCG(II,LL,MM) = COCG(II,LL,MM)
     &            - COEFBP(IFAC)*DIIPB(MM,IFAC)*SURFBO(LL,IFAC)
              ENDDO
C
            ELSE
C
C VECTORISATION NON FORCEE
              DO IFAC=1,NFABOR
                II = IFABOR(IFAC)
                COCG(II,LL,MM) = COCG(II,LL,MM)
     &            - COEFBP(IFAC)*DIIPB(MM,IFAC)*SURFBO(LL,IFAC)
              ENDDO
C
            ENDIF
C
          ENDDO
        ENDDO
C
C ---> ON INVERSE POUR TOUTE LES CELLULES : LE COCG POUR LES CELLULES INTERNES
C      RESTE ENSUITE LE MEME TANT QUE LE MAILLAGE NE CHANGE PAS
        NBLOC = NCEL/LBLOC
        IF (NBLOC.GT.0) THEN
          DO IBLOC = 1, NBLOC
            DO II =1, LBLOC
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
              A21=COCG31*COCG23-COCG21*COCG33
              A22=COCG11*COCG33-COCG31*COCG13
              A23=COCG21*COCG13-COCG11*COCG23
              A31=COCG21*COCG32-COCG31*COCG22
              A32=COCG31*COCG12-COCG11*COCG32
              A33=COCG11*COCG22-COCG21*COCG12
C
              UNSDET = 1.D0/(COCG11*A11 +COCG21*A12+COCG31*A13)
C
              AA(II,1,1) = A11*UNSDET
              AA(II,1,2) = A12*UNSDET
              AA(II,1,3) = A13*UNSDET
              AA(II,2,1) = A21*UNSDET
              AA(II,2,2) = A22*UNSDET
              AA(II,2,3) = A23*UNSDET
              AA(II,3,1) = A31*UNSDET
              AA(II,3,2) = A32*UNSDET
              AA(II,3,3) = A33*UNSDET
C
            ENDDO
C
            DO II = 1, LBLOC
              IEL = (IBLOC-1)*LBLOC+II
              COCG(IEL,1,1) = AA(II,1,1)
              COCG(IEL,1,2) = AA(II,1,2)
              COCG(IEL,1,3) = AA(II,1,3)
              COCG(IEL,2,1) = AA(II,2,1)
              COCG(IEL,2,2) = AA(II,2,2)
              COCG(IEL,2,3) = AA(II,2,3)
              COCG(IEL,3,1) = AA(II,3,1)
              COCG(IEL,3,2) = AA(II,3,2)
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
            A21=COCG31*COCG23-COCG21*COCG33
            A22=COCG11*COCG33-COCG31*COCG13
            A23=COCG21*COCG13-COCG11*COCG23
            A31=COCG21*COCG32-COCG31*COCG22
            A32=COCG31*COCG12-COCG11*COCG32
            A33=COCG11*COCG22-COCG21*COCG12
C
            UNSDET = 1.D0/(COCG11*A11 +COCG21*A12+COCG31*A13)
C
            AA(II,1,1) = A11*UNSDET
            AA(II,1,2) = A12*UNSDET
            AA(II,1,3) = A13*UNSDET
            AA(II,2,1) = A21*UNSDET
            AA(II,2,2) = A22*UNSDET
            AA(II,2,3) = A23*UNSDET
            AA(II,3,1) = A31*UNSDET
            AA(II,3,2) = A32*UNSDET
            AA(II,3,3) = A33*UNSDET
C
          ENDDO
C
          DO II = 1, IREL
            IEL = (IBLOC-1)*LBLOC+II
            COCG(IEL,1,1) = AA(II,1,1)
            COCG(IEL,1,2) = AA(II,1,2)
            COCG(IEL,1,3) = AA(II,1,3)
            COCG(IEL,2,1) = AA(II,2,1)
            COCG(IEL,2,2) = AA(II,2,2)
            COCG(IEL,2,3) = AA(II,2,3)
            COCG(IEL,3,1) = AA(II,3,1)
            COCG(IEL,3,2) = AA(II,3,2)
            COCG(IEL,3,3) = AA(II,3,3)
          ENDDO
C
        ENDIF
C
      ENDIF
C
C ---> SI ON DOIT RECALCULER COCG ENSUITE, ON NE LE FAIT PLUS
C      QUE POUR LES CELLULES DE BORD, AVEC LE COCGB STOCKE
C
      IF(ICCOCG.EQ.1 .AND. IPASS.GT.1 .AND. IALE.EQ.0) THEN
C
        DO LL =1,3
          DO MM =1,3
C
            DO KK = 1, NCELBR
              IEL = ICELBR(KK)
              COCG(IEL,LL,MM) = COCGB(KK,LL,MM)
            ENDDO
C
            IF (IVECTB.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
              DO IFAC=1,NFABOR
                II = IFABOR(IFAC)
                COCG(II,LL,MM) = COCG(II,LL,MM)
     &            - COEFBP(IFAC)*DIIPB(MM,IFAC)*SURFBO(LL,IFAC)
              ENDDO
C
            ELSE
C
C VECTORISATION NON FORCEE
              DO IFAC=1,NFABOR
                II = IFABOR(IFAC)
                COCG(II,LL,MM) = COCG(II,LL,MM)
     &            - COEFBP(IFAC)*DIIPB(MM,IFAC)*SURFBO(LL,IFAC)
              ENDDO
C
            ENDIF
C
          ENDDO
        ENDDO
C
        DO II = 1, NCELBR
C
          IEL = ICELBR(II)
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
          A21=COCG31*COCG23-COCG21*COCG33
          A22=COCG11*COCG33-COCG31*COCG13
          A23=COCG21*COCG13-COCG11*COCG23
          A31=COCG21*COCG32-COCG31*COCG22
          A32=COCG31*COCG12-COCG11*COCG32
          A33=COCG11*COCG22-COCG21*COCG12
C
          UNSDET = 1.D0/(COCG11*A11 +COCG21*A12+COCG31*A13)
C
          A11 = A11*UNSDET
          A12 = A12*UNSDET
          A13 = A13*UNSDET
          A21 = A21*UNSDET
          A22 = A22*UNSDET
          A23 = A23*UNSDET
          A31 = A31*UNSDET
          A32 = A32*UNSDET
          A33 = A33*UNSDET
C
          COCG(IEL,1,1) = A11
          COCG(IEL,1,2) = A12
          COCG(IEL,1,3) = A13
          COCG(IEL,2,1) = A21
          COCG(IEL,2,2) = A22
          COCG(IEL,2,3) = A23
          COCG(IEL,3,1) = A31
          COCG(IEL,3,2) = A32
          COCG(IEL,3,3) = A33
        ENDDO
C
      ENDIF
C
C ---> CALCUL DU RESIDU DE NORMALISATION
C
      CALL PRODS3(NCELET,NCEL,ISQRT,BX,BX,BY,BY,BZ,BZ,
     &            RNORX,RNORY,RNORZ)
      RNORM = RNORX +RNORY +RNORZ
      IF( RNORM.LE.EPZERO ) RETURN
C
C  LE VECTEUR OijFij EST CALCULE DANS CLDIJP
C
C ---> DEBUT DES ITERATIONS
C
 100  CONTINUE
C
      ISWEEP = ISWEEP +1
C
C
C     CALCUL DU SECOND MEMBRE
C
      DO IEL = 1, NCEL
        BX(IEL) = -DPDX(IEL)*VOLUME(IEL)
        BY(IEL) = -DPDY(IEL)*VOLUME(IEL)
        BZ(IEL) = -DPDZ(IEL)*VOLUME(IEL)
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
            VECFAC = POND(IFAC)*PVAR(II)  +(1.D0-POND(IFAC))*PVAR(JJ)
     &           +(DOFIJ(1,IFAC)*(DPDX(II)+DPDX(JJ))
     &           + DOFIJ(2,IFAC)*(DPDY(II)+DPDY(JJ))
     &           + DOFIJ(3,IFAC)*(DPDZ(II)+DPDZ(JJ)))*0.5D0
            PFSX = VECFAC*SURFAC(1,IFAC)
            PFSY = VECFAC*SURFAC(2,IFAC)
            PFSZ = VECFAC*SURFAC(3,IFAC)
            BX(II) = BX(II) +PFSX
            BY(II) = BY(II) +PFSY
            BZ(II) = BZ(II) +PFSZ
            BX(JJ) = BX(JJ) -PFSX
            BY(JJ) = BY(JJ) -PFSY
            BZ(JJ) = BZ(JJ) -PFSZ
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            VECFAC = POND(IFAC)*PVAR(II)  +(1.D0-POND(IFAC))*PVAR(JJ)
     &           +(DOFIJ(1,IFAC)*(DPDX(II)+DPDX(JJ))
     &           + DOFIJ(2,IFAC)*(DPDY(II)+DPDY(JJ))
     &           + DOFIJ(3,IFAC)*(DPDZ(II)+DPDZ(JJ)))*0.5D0
            PFSX = VECFAC*SURFAC(1,IFAC)
            PFSY = VECFAC*SURFAC(2,IFAC)
            PFSZ = VECFAC*SURFAC(3,IFAC)
            BX(II) = BX(II) +PFSX
            BY(II) = BY(II) +PFSY
            BZ(II) = BZ(II) +PFSZ
            BX(JJ) = BX(JJ) -PFSX
            BY(JJ) = BY(JJ) -PFSY
            BZ(JJ) = BZ(JJ) -PFSZ
          ENDDO
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
            II = IFABOR(IFAC)
            PIP =  PVAR(II)
     &           +DIIPB(1,IFAC)*DPDX(II) +DIIPB(2,IFAC)*DPDY(II)
     &           +DIIPB(3,IFAC)*DPDZ(II)
            PFAC = INC*COEFAP(IFAC) +COEFBP(IFAC)*PIP
            PFAC1= PVAR(II) +(CDGFBO(1,IFAC)-XYZCEN(1,II))*DPDX(II)
     &           +(CDGFBO(2,IFAC)-XYZCEN(2,II))*DPDY(II)
     &           +(CDGFBO(3,IFAC)-XYZCEN(3,II))*DPDZ(II)
            PFAC = COEFBP(IFAC)*(EXTRAP*PFAC1 +(1.D0-EXTRAP)*PFAC)
     &           +(1.D0-COEFBP(IFAC))*PFAC
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
            PIP =  PVAR(II)
     &           +DIIPB(1,IFAC)*DPDX(II) +DIIPB(2,IFAC)*DPDY(II)
     &           +DIIPB(3,IFAC)*DPDZ(II)
            PFAC = INC*COEFAP(IFAC) +COEFBP(IFAC)*PIP
            PFAC1= PVAR(II) +(CDGFBO(1,IFAC)-XYZCEN(1,II))*DPDX(II)
     &           +(CDGFBO(2,IFAC)-XYZCEN(2,II))*DPDY(II)
     &           +(CDGFBO(3,IFAC)-XYZCEN(3,II))*DPDZ(II)
            PFAC = COEFBP(IFAC)*(EXTRAP*PFAC1 +(1.D0-EXTRAP)*PFAC)
     &           +(1.D0-COEFBP(IFAC))*PFAC
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
            FMOYX=0.5D0*(FEXTX(II)+FEXTX(JJ))
            FMOYY=0.5D0*(FEXTY(II)+FEXTY(JJ))
            FMOYZ=0.5D0*(FEXTZ(II)+FEXTZ(JJ))
            VECFAC   = POND(IFAC)*(PVAR(II)
     &           -(XYZCEN(1,II)-CDGFAC(1,IFAC))*(FEXTX(II)-FMOYX)
     &           -(XYZCEN(2,II)-CDGFAC(2,IFAC))*(FEXTY(II)-FMOYY)
     &           -(XYZCEN(3,II)-CDGFAC(3,IFAC))*(FEXTZ(II)-FMOYZ))
     &           +(1.D0-POND(IFAC))*(PVAR(JJ)
     &           -(XYZCEN(1,JJ)-CDGFAC(1,IFAC))*(FEXTX(JJ)-FMOYX)
     &           -(XYZCEN(2,JJ)-CDGFAC(2,IFAC))*(FEXTY(JJ)-FMOYY)
     &           -(XYZCEN(3,JJ)-CDGFAC(3,IFAC))*(FEXTZ(JJ)-FMOYZ))
     &           +(DOFIJ(1,IFAC)*(DPDX(II)+DPDX(JJ))
     &           + DOFIJ(2,IFAC)*(DPDY(II)+DPDY(JJ))
     &           + DOFIJ(3,IFAC)*(DPDZ(II)+DPDZ(JJ)))*0.5D0
            PFSX = VECFAC*SURFAC(1,IFAC)
            PFSY = VECFAC*SURFAC(2,IFAC)
            PFSZ = VECFAC*SURFAC(3,IFAC)
            BX(II) = BX(II) +PFSX
            BY(II) = BY(II) +PFSY
            BZ(II) = BZ(II) +PFSZ
            BX(JJ) = BX(JJ) -PFSX
            BY(JJ) = BY(JJ) -PFSY
            BZ(JJ) = BZ(JJ) -PFSZ
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            FMOYX=0.5D0*(FEXTX(II)+FEXTX(JJ))
            FMOYY=0.5D0*(FEXTY(II)+FEXTY(JJ))
            FMOYZ=0.5D0*(FEXTZ(II)+FEXTZ(JJ))
            VECFAC   = POND(IFAC)*(PVAR(II)
     &           -(XYZCEN(1,II)-CDGFAC(1,IFAC))*(FEXTX(II)-FMOYX)
     &           -(XYZCEN(2,II)-CDGFAC(2,IFAC))*(FEXTY(II)-FMOYY)
     &           -(XYZCEN(3,II)-CDGFAC(3,IFAC))*(FEXTZ(II)-FMOYZ))
     &           +(1.D0-POND(IFAC))*(PVAR(JJ)
     &           -(XYZCEN(1,JJ)-CDGFAC(1,IFAC))*(FEXTX(JJ)-FMOYX)
     &           -(XYZCEN(2,JJ)-CDGFAC(2,IFAC))*(FEXTY(JJ)-FMOYY)
     &           -(XYZCEN(3,JJ)-CDGFAC(3,IFAC))*(FEXTZ(JJ)-FMOYZ))
     &           +(DOFIJ(1,IFAC)*(DPDX(II)+DPDX(JJ))
     &           + DOFIJ(2,IFAC)*(DPDY(II)+DPDY(JJ))
     &           + DOFIJ(3,IFAC)*(DPDZ(II)+DPDZ(JJ)))*0.5D0
            PFSX = VECFAC*SURFAC(1,IFAC)
            PFSY = VECFAC*SURFAC(2,IFAC)
            PFSZ = VECFAC*SURFAC(3,IFAC)
            BX(II) = BX(II) +PFSX
            BY(II) = BY(II) +PFSY
            BZ(II) = BZ(II) +PFSZ
            BX(JJ) = BX(JJ) -PFSX
            BY(JJ) = BY(JJ) -PFSY
            BZ(JJ) = BZ(JJ) -PFSZ
          ENDDO
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
            II = IFABOR(IFAC)
            PIP =  PVAR(II)
     &           +DIIPB(1,IFAC)*DPDX(II) +DIIPB(2,IFAC)*DPDY(II)
     &           +DIIPB(3,IFAC)*DPDZ(II)
            PFAC = INC*COEFAP(IFAC) +COEFBP(IFAC)*(PIP
     &           -(XYZCEN(1,II)-CDGFBO(1,IFAC)+DIIPB(1,IFAC))*FEXTX(II)
     &           -(XYZCEN(2,II)-CDGFBO(2,IFAC)+DIIPB(2,IFAC))*FEXTY(II)
     &           -(XYZCEN(3,II)-CDGFBO(3,IFAC)+DIIPB(3,IFAC))*FEXTZ(II))
            PFAC1= PVAR(II) +(CDGFBO(1,IFAC)-XYZCEN(1,II))*DPDX(II)
     &                    +(CDGFBO(2,IFAC)-XYZCEN(2,II))*DPDY(II)
     &                    +(CDGFBO(3,IFAC)-XYZCEN(3,II))*DPDZ(II)
            PFAC = COEFBP(IFAC)*(EXTRAP*PFAC1 +(1.D0-EXTRAP)*PFAC)
     &           +(1.D0-COEFBP(IFAC))*PFAC
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
            PIP =  PVAR(II)
     &           +DIIPB(1,IFAC)*DPDX(II) +DIIPB(2,IFAC)*DPDY(II)
     &           +DIIPB(3,IFAC)*DPDZ(II)
            PFAC = INC*COEFAP(IFAC) +COEFBP(IFAC)*(PIP
     &           -(XYZCEN(1,II)-CDGFBO(1,IFAC)+DIIPB(1,IFAC))*FEXTX(II)
     &           -(XYZCEN(2,II)-CDGFBO(2,IFAC)+DIIPB(2,IFAC))*FEXTY(II)
     &           -(XYZCEN(3,II)-CDGFBO(3,IFAC)+DIIPB(3,IFAC))*FEXTZ(II))
            PFAC1= PVAR(II) +(CDGFBO(1,IFAC)-XYZCEN(1,II))*DPDX(II)
     &                    +(CDGFBO(2,IFAC)-XYZCEN(2,II))*DPDY(II)
     &                    +(CDGFBO(3,IFAC)-XYZCEN(3,II))*DPDZ(II)
            PFAC = COEFBP(IFAC)*(EXTRAP*PFAC1 +(1.D0-EXTRAP)*PFAC)
     &           +(1.D0-COEFBP(IFAC))*PFAC
            BX(II) = BX(II) +PFAC*SURFBO(1,IFAC)
            BY(II) = BY(II) +PFAC*SURFBO(2,IFAC)
            BZ(II) = BZ(II) +PFAC*SURFBO(3,IFAC)
          ENDDO
C
        ENDIF
C
      ENDIF
C
C     INCREMENTATION DU GRADIENT
C
      DO IEL =1,NCEL
C
        DELTPX =
     & COCG(IEL,1,1)*BX(IEL)+COCG(IEL,1,2)*BY(IEL)+COCG(IEL,1,3)*BZ(IEL)
        DELTPY =
     & COCG(IEL,2,1)*BX(IEL)+COCG(IEL,2,2)*BY(IEL)+COCG(IEL,2,3)*BZ(IEL)
        DELTPZ =
     & COCG(IEL,3,1)*BX(IEL)+COCG(IEL,3,2)*BY(IEL)+COCG(IEL,3,3)*BZ(IEL)
C
        DPDX(IEL) = DPDX(IEL) +DELTPX
        DPDY(IEL) = DPDY(IEL) +DELTPY
        DPDZ(IEL) = DPDZ(IEL) +DELTPZ
C
      ENDDO
C
C
C     TRAITEMENT DU PARALLELISME
C
      IF(IRANGP.GE.0) THEN
        CALL PARCOM (DPDX)
C       ===========
        CALL PARCOM (DPDY)
C       ===========
        CALL PARCOM (DPDZ)
C       ===========
      ENDIF
C
C     TRAITEMENT DE LA PERIODICITE
C
      IF(IPERIO.EQ.1) THEN
        CALL PERCOM
C       ===========
     &  ( IDIMTE , ITENSO ,
     &    DPDX   , DPDX   , DPDX  ,
     &    DPDY   , DPDY   , DPDY  ,
     &    DPDZ   , DPDZ   , DPDZ  )
      ENDIF
C
C
C ---> TEST DE CONVERGENCE
C
      CALL PRODS3(NCELET,NCEL,ISQRT,BX,BX,BY,BY,BZ,BZ,
     &            RNORX,RNORY,RNORZ)
      RESIDU = RNORX +RNORY +RNORZ
C
      IF( RESIDU.LE.EPSRGP*RNORM) THEN
        IF( IWARNP.GE.2 ) THEN
          WRITE (NFECRA,1000) ISWEEP,RESIDU/RNORM,RNORM,IVAR
        ENDIF
        GOTO 101
      ELSEIF( ISWEEP.GE.NSWMAX ) THEN
        IF( IWARNP.GE.0) THEN
           WRITE (NFECRA,1000)ISWEEP,RESIDU/RNORM,RNORM,IVAR
           WRITE (NFECRA,1100)
        ENDIF
        GOTO 101
      ELSE
        GOTO 100
      ENDIF
C
 101  CONTINUE
C
C
C
C--------
C FORMATS
C--------
#if defined(_CS_LANG_FR)
C
 1000 FORMAT(1X,'GRADRC ISWEEP = ',I4,' RESIDU NORME: ',E11.4,
     &         ' NORME: ',E11.4,/,1X,'PARAMETRE IVAR = ',I4 )
 1100 FORMAT(
     &'@                                                            ',/,
     &'@ @@ ATTENTION :         NON CONVERGENCE DE GRADRC           ',/,
     &'@    *********                                               ',/,
     &'@                                                            '  )
C
#else
C
 1000 FORMAT(1X,'GRADRC ISWEEP = ',I4,' NORMED RESIDUAL: ',E11.4,
     &         ' NORM: ',E11.4,/,1X,'PARAMETER IVAR = ',I4 )
 1100 FORMAT(
     &'@'                                                            ,/,
     &'@ @@ WARNING:            NON CONVERGENCE OF GRADRC'           ,/,
     &'@    ********'                                                ,/,
     &'@'                                                              )
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
