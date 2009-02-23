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
                        SUBROUTINE CFBSC2
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , ICONVP , IDIFFP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , INC    , IMRGRA , ICCOCG , IIFBRU ,
     &   IPP    , IWARNP ,
     &   BLENCP , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , IFRUSB ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PVAR   , COEFAP , COEFBP , COFAFP , COFBFP ,
     &   FLUMAS , FLUMAB , VISCF  , VISCB  ,
     &   SMBRP  ,
     &   DPDX   , DPDY   , DPDZ   , DPDXA  , DPDYA  , DPDZA  ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ---------
c@foncb
CFONC
CFONC CALCUL DU BILAN EXPLICITE DE LA VARIABLE PVAR (VITESSE,SCALAIRES)
CFONC
CFONC                    -- .                  ----->        -->
CFONC  SMBRP = SMBRP -  \   m   PVAR +Visc   ( grad PVAR )  . n
CFONC                   /__  ij    ij     ij             ij    ij
CFONC                    j
CFONC
CFONC ATTENTION : SMBRP DEJA INITIALISE AVANT L'APPEL A BILSCA
CFONC             IL CONTIENT LES TERMES SOURCES EXPLICITES, ETC....
CFONC
CFONC BLENCP = 1 : PAS D'UPWIND EN DEHORS DU TEST DE PENTE
CFONC BLENCP = 0 : UPWIND
CFONC ISCHCP = 1 : CENTRE
CFONC ISCHCP = 0 : SECOND ORDER
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                 ARGUMENTS
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
CARGU ! IVAR         ! E  !  -> ! NUMERO DE LA VARIABLE                !
CARGU ! ICONVP       ! E  !  -> ! INDICATEUR = 1 CONVECTION, 0 SINON   !
CARGU ! IDIFFP       ! E  !  -> ! INDICATEUR = 1 DIFFUSION , 0 SINON   !
CARGU ! NSWRGP       ! E  !  -> ! NOMBRE DE SWEEP POUR RECONSTRUCTION  !
CARGU !              !    !     !             DES GRADIENTS            !
CARGU ! IMLIGP       ! E  !  -> ! METHODE DE LIMITATION DU GRADIENT    !
CARGU !              !    !     !  < 0 PAS DE LIMITATION               !
CARGU !              !    !     !  = 0 A PARTIR DES GRADIENTS VOISINS  !
CARGU !              !    !     !  = 1 A PARTIR DU GRADIENT MOYEN      !
CARGU ! IRCFLP       ! E  !  -> ! INDICATEUR = 1 REC FLUX, 0 SINON     !
CARGU ! ISCHCP       ! E  !  -> ! INDICATEUR = 1 CENTRE , 0 2ND ORDER  !
CARGU ! ISSTPP       ! E  !  -> ! INDICATEUR = 1 SANS TEST DE PENTE    !
CARGU !              !    !     !            = 0 AVEC TEST DE PENTE    !
CARGU ! INC          ! E  !  -> ! INDICATEUR = 0 RESOL SUR INCREMENT   !
CARGU !              !    !     !              1 SINON                 !
CARGU ! IMRGRA       ! E  !  -> ! INDICATEUR = 0 GRADRC 97             !
CARGU !              ! E  !  -> !            = 1 GRADMC 99             !
CARGU ! ICCOCG       ! E  !  -> ! INDICATEUR = 1 POUR RECALCUL DE COCG !
CARGU !              !    !     !              0 SINON                 !
CARGU ! IIFBRU       ! E  !  -> ! POINTEUR FLUX DE BORD RUSANOV        !
CARGU ! IPP          ! E  !  -> ! NUMERO DE VARIABLE POUR POST         !
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! BLENCP       ! R  !  -> ! 1 - PROPORTION D'UPWIND              !
CARGU ! EPSRGP       ! R  !  -> ! PRECISION RELATIVE POUR LA           !
CARGU !              !    !     !  RECONSTRUCTION DES GRADIENTS 97     !
CARGU ! CLIMGP       ! R  !  -> ! COEF GRADIENT*DISTANCE/ECART         !
CARGU ! EXTRAP       ! R  !  -> ! COEF EXTRAP GRADIENT                 !
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
CARGU ! IFRUSB(NFABOR! TE !  -> ! INDICATEUR FLUX DE RUSANOV           !
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
CARGU ! IA(*)        ! TE !  -  ! MACRO TABLEAU ENTIER                 !
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
CARGU ! PVAR (NCELET ! TR !  -> ! VARIABLE RESOLUE (INSTANT PRECEDENT) !
CARGU ! COEFAP, B    ! TR !  -> ! TABLEAUX DES COND LIM POUR P         !
CARGU !   (NFABOR)   !    !     !  SUR LA NORMALE A LA FACE DE BORD    !
CARGU ! COFAFP, B    ! TR !  -> ! TABLEAUX DES COND LIM POUR LE FLUX DE!
CARGU !   (NFABOR)   !    !     !  DIFFUSION DE P                      !
CARGU ! FLUMAS(NFAC) ! TR !  -> ! FLUX DE MASSE AUX FACES INTERNES     !
CARGU ! FLUMAB(NFABOR! TR !  -> ! FLUX DE MASSE AUX FACES DE BORD      !
CARGU ! VISCF (NFAC) ! TR !  -> ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU !              !    !     !  POUR SECOND MEMBRE                  !
CARGU ! VISCB (NFABOR! TR !  -> ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU !              !    !     !  POUR SECOND MEMBRE                  !
CARGU ! SMBRP(NCELET ! TR ! <-> ! BILAN AU SECOND MEMBRE               !
CARGU ! DPDX,Y,Z     ! TR !  -  ! TABLEAU DE TRAVAIL POUR LE GRAD DE P !
CARGU !    (NCELET)  !    !     !                                      !
CARGU ! DPDXA,YA,ZA  ! TR !  -  ! TABLEAU DE TRAVAIL POUR LE GRAD DE P !
CARGU !    (NCELET)  !    !     !  AVEC DECENTREMENT AMONT             !
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
      INCLUDE "pointe.h"
      INCLUDE "vector.h"
      INCLUDE "entsor.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
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
      INTEGER          IVAR   , ICONVP , IDIFFP , NSWRGP , IMLIGP
      INTEGER          IRCFLP , ISCHCP , ISSTPP
      INTEGER          INC    , IMRGRA , ICCOCG , IIFBRU
      INTEGER          IWARNP , IPP
      DOUBLE PRECISION BLENCP , EPSRGP , CLIMGP, EXTRAP
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IFRUSB(NFABOR)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION PVAR (NCELET), COEFAP(NFABOR), COEFBP(NFABOR)
      DOUBLE PRECISION                COFAFP(NFABOR), COFBFP(NFABOR)
      DOUBLE PRECISION FLUMAS(NFAC), FLUMAB(NFABOR)
      DOUBLE PRECISION VISCF (NFAC), VISCB (NFABOR)
      DOUBLE PRECISION SMBRP(NCELET)
      DOUBLE PRECISION DPDX (NCELET),DPDY (NCELET),DPDZ (NCELET)
      DOUBLE PRECISION DPDXA(NCELET),DPDYA(NCELET),DPDZA(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      CHARACTER*8      CNOM
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC,II,JJ,INFAC,IEL,IUPWIN, IIJ, III, IOK
      INTEGER          ITENSO, IDIMTE, IPHYDP
      INTEGER          IIU(NPHSMX),IIV(NPHSMX),IIW(NPHSMX)
      INTEGER          IITYTU(NPHSMX)
      INTEGER          IIR11(NPHSMX),IIR22(NPHSMX),IIR33(NPHSMX)
      INTEGER          IIR12(NPHSMX),IIR13(NPHSMX),IIR23(NPHSMX)
      DOUBLE PRECISION PFAC,PFACD,PIP,PJP,FLUI,FLUJ,FLUX
      DOUBLE PRECISION DIFX,DIFY,DIFZ,DJFX,DJFY,DJFZ,PIF,PJF
      DOUBLE PRECISION TESTI,TESTJ,TESTIJ
      DOUBLE PRECISION DPXF,DPYF,DPZF
      DOUBLE PRECISION DCC, DDI, DDJ, TESQCK
      DOUBLE PRECISION DIJPFX, DIJPFY, DIJPFZ
      DOUBLE PRECISION DIIPFX, DIIPFY, DIIPFZ
      DOUBLE PRECISION DJJPFX, DJJPFY, DJJPFZ
      DOUBLE PRECISION DIIPBX, DIIPBY, DIIPBZ
      DOUBLE PRECISION POND, DIST, SURFAN
      DOUBLE PRECISION PFAC1, PFAC2, PFAC3, UNSVOL
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATION
c=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      CHAINE = NOMVAR(IPP)
      CNOM   = CHAINE(1:8)
C
      IF(IWARNP.GE.2) THEN
        IF (ISCHCP.EQ.1) THEN
          WRITE(NFECRA,1000)CNOM,'    CENTRE ',
     &                                              (1.D0-BLENCP)*100.D0
        ELSE
          WRITE(NFECRA,1000)CNOM,' 2ND ORDER ',
     &                                              (1.D0-BLENCP)*100.D0
        ENDIF
      ENDIF
C
      IUPWIN = 0
      IF(BLENCP.EQ.0.D0) IUPWIN = 1
C
C=======================================================================
C 2.  CALCUL DU BILAN AVEC TECHNIQUE DE RECONSTRUCTION
C=======================================================================
C
C ======================================================================
C ---> CALCUL DU GRADIENT DE P
C ======================================================================
C    DPDX sert a la fois pour la reconstruction des flux et pour le test
C    de pente. On doit donc le calculer :
C        - quand on a de la diffusion et qu'on reconstruit les flux
C        - quand on a de la convection SOLU
C        - quand on a de la convection, qu'on n'est pas en upwind pur
C          et qu'on reconstruit les flux
C        - quand on a de la convection, qu'on n'est pas en upwind pur
C          et qu'on n'a pas shunte le test de pente
C
      IF( (IDIFFP.NE.0 .AND. IRCFLP.EQ.1) .OR.
     &    (ICONVP.NE.0 .AND. IUPWIN.EQ.0 .AND.
     &       (ISCHCP.EQ.0 .OR. IRCFLP.EQ.1 .OR. ISSTPP.EQ.0)) ) THEN
C
        IPHYDP = 0
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP ,  IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DPDXA  , DPDXA  , DPDXA  ,
     &   PVAR   , COEFAP , COEFBP ,
     &   DPDX   , DPDY   , DPDZ   ,
C        ------   ------   ------
     &   DPDXA  , DPDYA  , DPDZA  ,
     &   RDEVEL , RTUSER , RA     )
C
      ELSE
        DO IEL = 1, NCELET
          DPDX(IEL) = 0.D0
          DPDY(IEL) = 0.D0
          DPDZ(IEL) = 0.D0
        ENDDO
      ENDIF
C
C
C ======================================================================
C ---> CALCUL DU GRADIENT DECENTRE DPDXA, DPDYA, DPDZA POUR TST DE PENTE
C ======================================================================
C
      DO IEL = 1, NCELET
        DPDXA(IEL) = 0.D0
        DPDYA(IEL) = 0.D0
        DPDZA(IEL) = 0.D0
      ENDDO
C
      IF( ICONVP.GT.0.AND.IUPWIN.EQ.0.AND.ISSTPP.EQ.0 ) THEN
C
        IF (IVECTI.EQ.1) THEN
C
!CDIR NODEP
          DO IFAC = 1, NFAC
C
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
C
            DIFX = CDGFAC(1,IFAC) - XYZCEN(1,II)
            DIFY = CDGFAC(2,IFAC) - XYZCEN(2,II)
            DIFZ = CDGFAC(3,IFAC) - XYZCEN(3,II)
            DJFX = CDGFAC(1,IFAC) - XYZCEN(1,JJ)
            DJFY = CDGFAC(2,IFAC) - XYZCEN(2,JJ)
            DJFZ = CDGFAC(3,IFAC) - XYZCEN(3,JJ)
C
            PIF = PVAR(II) +DIFX*DPDX(II)+DIFY*DPDY(II)+DIFZ*DPDZ(II)
            PJF = PVAR(JJ) +DJFX*DPDX(JJ)+DJFY*DPDY(JJ)+DJFZ*DPDZ(JJ)
C
            PFAC = PJF
            IF( FLUMAS(IFAC ).GT.0.D0 ) PFAC = PIF
C
            PFAC1 = PFAC*SURFAC(1,IFAC )
            PFAC2 = PFAC*SURFAC(2,IFAC )
            PFAC3 = PFAC*SURFAC(3,IFAC )
C
            DPDXA(II) = DPDXA(II) +PFAC1
            DPDYA(II) = DPDYA(II) +PFAC2
            DPDZA(II) = DPDZA(II) +PFAC3
C
            DPDXA(JJ) = DPDXA(JJ) -PFAC1
            DPDYA(JJ) = DPDYA(JJ) -PFAC2
            DPDZA(JJ) = DPDZA(JJ) -PFAC3
C
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1, NFAC
C
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
C
            DIFX = CDGFAC(1,IFAC) - XYZCEN(1,II)
            DIFY = CDGFAC(2,IFAC) - XYZCEN(2,II)
            DIFZ = CDGFAC(3,IFAC) - XYZCEN(3,II)
            DJFX = CDGFAC(1,IFAC) - XYZCEN(1,JJ)
            DJFY = CDGFAC(2,IFAC) - XYZCEN(2,JJ)
            DJFZ = CDGFAC(3,IFAC) - XYZCEN(3,JJ)
C
            PIF = PVAR(II) +DIFX*DPDX(II)+DIFY*DPDY(II)+DIFZ*DPDZ(II)
            PJF = PVAR(JJ) +DJFX*DPDX(JJ)+DJFY*DPDY(JJ)+DJFZ*DPDZ(JJ)
C
            PFAC = PJF
            IF( FLUMAS(IFAC ).GT.0.D0 ) PFAC = PIF
C
            PFAC1 = PFAC*SURFAC(1,IFAC )
            PFAC2 = PFAC*SURFAC(2,IFAC )
            PFAC3 = PFAC*SURFAC(3,IFAC )
C
            DPDXA(II) = DPDXA(II) +PFAC1
            DPDYA(II) = DPDYA(II) +PFAC2
            DPDZA(II) = DPDZA(II) +PFAC3
C
            DPDXA(JJ) = DPDXA(JJ) -PFAC1
            DPDYA(JJ) = DPDYA(JJ) -PFAC2
            DPDZA(JJ) = DPDZA(JJ) -PFAC3
C
          ENDDO
C
        ENDIF
C
        IF (IVECTB.EQ.1) THEN
C
!CDIR NODEP
          DO IFAC = 1, NFABOR
            II = IFABOR(IFAC )
            III = IDIIPB-1+3*(IFAC-1)
            DIIPBX = RA(III+1)
            DIIPBY = RA(III+2)
            DIIPBZ = RA(III+3)
            PFAC = INC*COEFAP(IFAC )
     &            +COEFBP(IFAC )*(PVAR(II)+DIIPBX*DPDX(II)
     &                    +DIIPBY*DPDY(II)+DIIPBZ*DPDZ(II) )
            DPDXA(II) = DPDXA(II) +PFAC*SURFBO(1,IFAC )
            DPDYA(II) = DPDYA(II) +PFAC*SURFBO(2,IFAC )
            DPDZA(II) = DPDZA(II) +PFAC*SURFBO(3,IFAC )
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC =1,NFABOR
            II = IFABOR(IFAC )
            III = IDIIPB-1+3*(IFAC-1)
            DIIPBX = RA(III+1)
            DIIPBY = RA(III+2)
            DIIPBZ = RA(III+3)
            PFAC = INC*COEFAP(IFAC )
     &            +COEFBP(IFAC )*(PVAR(II)+DIIPBX*DPDX(II)
     &                    +DIIPBY*DPDY(II)+DIIPBZ*DPDZ(II) )
            DPDXA(II) = DPDXA(II) +PFAC*SURFBO(1,IFAC )
            DPDYA(II) = DPDYA(II) +PFAC*SURFBO(2,IFAC )
            DPDZA(II) = DPDZA(II) +PFAC*SURFBO(3,IFAC )
          ENDDO
C
        ENDIF
C
        DO IEL = 1, NCEL
          UNSVOL = 1.D0/VOLUME(IEL)
          DPDXA(IEL) = DPDXA(IEL)*UNSVOL
          DPDYA(IEL) = DPDYA(IEL)*UNSVOL
          DPDZA(IEL) = DPDZA(IEL)*UNSVOL
        ENDDO
C
C     TRAITEMENT DU PARALLELISME
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (DPDXA)
C         ===========
          CALL PARCOM (DPDYA)
C         ===========
          CALL PARCOM (DPDZA)
C         ===========
        ENDIF
C
C TRAITEMENT DE LA PERIODICITE
C
C  On echange pour la translation
C   pour la rotation, on prend le gradient simple (pas de temps precedent)
C
        IF(IPERIO.EQ.1) THEN
C
C        Pour les rotations
C          avec la vitesse et les tensions de Reynolds,
C          on utilise la valeur du gradient simple (PERING) a defaut de mieux.
C        Dans les autres cas, on echange DPDXA
C
C        On recupere d'abord certains COMMON necessaires a PERING
C
          CALL PERGRA
C         ===========
     &  ( NPHSMX , NPHAS  ,
     &    IIU    , IIV    , IIW    ,
     &    IITYTU ,
     &    IIR11  , IIR22  , IIR33  , IIR12  , IIR13  , IIR23  )
C
          CALL PERING
C         ===========
     &  ( NPHAS  , IVAR   ,
     &    IDIMTE , ITENSO , IPEROT , IGUPER , IGRPER ,
     &    IIU    , IIV    , IIW    , IITYTU ,
     &    IIR11  , IIR22  , IIR33  , IIR12  , IIR13  , IIR23  ,
     &    DPDXA  , DPDYA  , DPDZA  ,
     &    RA(IDUDXY) , RA(IDRDXY)  )
C
          CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    DPDXA  , DPDXA  , DPDXA ,
     &    DPDYA  , DPDYA  , DPDYA ,
     &    DPDZA  , DPDZA  , DPDZA )
        ENDIF
C
      ENDIF
C
C
C ======================================================================
C ---> ASSEMBLAGE A PARTIR DES FACETTES FLUIDES
C ======================================================================
C
      INFAC = 0
C
      IF(NCELET.GT.NCEL) THEN
        DO IEL = NCEL+1, NCELET
          SMBRP(IEL) = 0.D0
        ENDDO
      ENDIF
C
C
C  --> FLUX UPWIND PUR
C  =====================
C
      IF(IUPWIN.EQ.1) THEN
C
        IF (IVECTI.EQ.1) THEN
C
!CDIR NODEP
          DO IFAC = 1, NFAC
C
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
C
            IIJ = IDIJPF-1+3*(IFAC-1)
            DIJPFX = RA(IIJ+1)
            DIJPFY = RA(IIJ+2)
            DIJPFZ = RA(IIJ+3)
C
            POND   = RA(IPOND-1+IFAC)
C
C ON RECALCULE A CE NIVEAU II' ET JJ'
C
            DIIPFX = CDGFAC(1,IFAC) - (XYZCEN(1,II)+
     &               (1.D0-POND) * DIJPFX)
            DIIPFY = CDGFAC(2,IFAC) - (XYZCEN(2,II)+
     &               (1.D0-POND) * DIJPFY)
            DIIPFZ = CDGFAC(3,IFAC) - (XYZCEN(3,II)+
     &               (1.D0-POND) * DIJPFZ)
            DJJPFX = CDGFAC(1,IFAC) -  XYZCEN(1,JJ)+
     &                   POND  * DIJPFX
            DJJPFY = CDGFAC(2,IFAC) -  XYZCEN(2,JJ)+
     &                   POND  * DIJPFY
            DJJPFZ = CDGFAC(3,IFAC) -  XYZCEN(3,JJ)+
     &                   POND  * DIJPFZ
C
            DPXF = 0.5D0*(DPDX(II) + DPDX(JJ))
            DPYF = 0.5D0*(DPDY(II) + DPDY(JJ))
            DPZF = 0.5D0*(DPDZ(II) + DPDZ(JJ))
C
C     reconstruction uniquement si IRCFLP = 1
            PIP = PVAR(II)
     &           + IRCFLP*(DPXF*DIIPFX+DPYF*DIIPFY+DPZF*DIIPFZ)
            PJP = PVAR(JJ)
     &           + IRCFLP*(DPXF*DJJPFX+DPYF*DJJPFY+DPZF*DJJPFZ)
C
            FLUI = 0.5D0*( FLUMAS(IFAC) +ABS(FLUMAS(IFAC)) )
            FLUJ = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
C
            PIF = PVAR(II)
            PJF = PVAR(JJ)
            INFAC = INFAC+1
C
            FLUX = ICONVP*( FLUI*PIF +FLUJ*PJF )
     &           + IDIFFP*VISCF(IFAC)*( PIP -PJP )
C
            SMBRP(II) = SMBRP(II) -FLUX
            SMBRP(JJ) = SMBRP(JJ) +FLUX
C
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1, NFAC
C
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
C
            IIJ = IDIJPF-1+3*(IFAC-1)
            DIJPFX = RA(IIJ+1)
            DIJPFY = RA(IIJ+2)
            DIJPFZ = RA(IIJ+3)
C
            POND   = RA(IPOND-1+IFAC)
C
C ON RECALCULE A CE NIVEAU II' ET JJ'
C
            DIIPFX = CDGFAC(1,IFAC) - (XYZCEN(1,II)+
     &               (1.D0-POND) * DIJPFX)
            DIIPFY = CDGFAC(2,IFAC) - (XYZCEN(2,II)+
     &               (1.D0-POND) * DIJPFY)
            DIIPFZ = CDGFAC(3,IFAC) - (XYZCEN(3,II)+
     &               (1.D0-POND) * DIJPFZ)
            DJJPFX = CDGFAC(1,IFAC) -  XYZCEN(1,JJ)+
     &                   POND  * DIJPFX
            DJJPFY = CDGFAC(2,IFAC) -  XYZCEN(2,JJ)+
     &                   POND  * DIJPFY
            DJJPFZ = CDGFAC(3,IFAC) -  XYZCEN(3,JJ)+
     &                   POND  * DIJPFZ
C
            DPXF = 0.5D0*(DPDX(II) + DPDX(JJ))
            DPYF = 0.5D0*(DPDY(II) + DPDY(JJ))
            DPZF = 0.5D0*(DPDZ(II) + DPDZ(JJ))
C
            PIP = PVAR(II)
     &           + IRCFLP*(DPXF*DIIPFX+DPYF*DIIPFY+DPZF*DIIPFZ)
            PJP = PVAR(JJ)
     &           + IRCFLP*(DPXF*DJJPFX+DPYF*DJJPFY+DPZF*DJJPFZ)
C
            FLUI = 0.5D0*( FLUMAS(IFAC) +ABS(FLUMAS(IFAC)) )
            FLUJ = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
C
            PIF = PVAR(II)
            PJF = PVAR(JJ)
            INFAC = INFAC+1
C
            FLUX = ICONVP*( FLUI*PIF +FLUJ*PJF )
     &           + IDIFFP*VISCF(IFAC)*( PIP -PJP )
C
            SMBRP(II) = SMBRP(II) -FLUX
            SMBRP(JJ) = SMBRP(JJ) +FLUX
C
          ENDDO
C
        ENDIF
C
C
C  --> FLUX SANS TEST DE PENTE
C  ============================
C
      ELSEIF(ISSTPP.EQ.1) THEN
C
        IF (IVECTI.EQ.1) THEN
C
          IOK = 0
!CDIR NODEP
          DO IFAC = 1, NFAC
C
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
C
            IIJ = IDIJPF-1+3*(IFAC-1)
C
            DIJPFX = RA(IIJ+1)
            DIJPFY = RA(IIJ+2)
            DIJPFZ = RA(IIJ+3)
C
            POND   = RA(IPOND-1+IFAC)
C
C ON RECALCULE A CE NIVEAU II' ET JJ'
C
C
            DIIPFX = CDGFAC(1,IFAC) - (XYZCEN(1,II)+
     &               (1.D0-POND) * DIJPFX)
            DIIPFY = CDGFAC(2,IFAC) - (XYZCEN(2,II)+
     &               (1.D0-POND) * DIJPFY)
            DIIPFZ = CDGFAC(3,IFAC) - (XYZCEN(3,II)+
     &               (1.D0-POND) * DIJPFZ)
            DJJPFX = CDGFAC(1,IFAC) -  XYZCEN(1,JJ)+
     &                   POND  * DIJPFX
            DJJPFY = CDGFAC(2,IFAC) -  XYZCEN(2,JJ)+
     &                   POND  * DIJPFY
            DJJPFZ = CDGFAC(3,IFAC) -  XYZCEN(3,JJ)+
     &                   POND  * DIJPFZ
C
            DPXF = 0.5D0*(DPDX(II) + DPDX(JJ))
            DPYF = 0.5D0*(DPDY(II) + DPDY(JJ))
            DPZF = 0.5D0*(DPDZ(II) + DPDZ(JJ))
C
C
            PIP = PVAR(II)
     &           + IRCFLP*(DPXF*DIIPFX+DPYF*DIIPFY+DPZF*DIIPFZ)
            PJP = PVAR(JJ)
     &           + IRCFLP*(DPXF*DJJPFX+DPYF*DJJPFY+DPZF*DJJPFZ)
C
            FLUI = 0.5D0*( FLUMAS(IFAC) +ABS(FLUMAS(IFAC)) )
            FLUJ = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
C
C
C         CENTRE
C        --------
C
            IF (ISCHCP.EQ.1) THEN
C
              PIF = POND*PIP +(1.D0-POND)*PJP
              PJF = PIF
C
C
C         SECOND ORDER
C        --------------
C
            ELSEIF(ISCHCP.EQ.0) THEN
C
              DIFX = CDGFAC(1,IFAC) - XYZCEN(1,II)
              DIFY = CDGFAC(2,IFAC) - XYZCEN(2,II)
              DIFZ = CDGFAC(3,IFAC) - XYZCEN(3,II)
              DJFX = CDGFAC(1,IFAC) - XYZCEN(1,JJ)
              DJFY = CDGFAC(2,IFAC) - XYZCEN(2,JJ)
              DJFZ = CDGFAC(3,IFAC) - XYZCEN(3,JJ)
C
C     on laisse la reconstruction de PIF et PJF meme si IRCFLP=0
C     sinon cela revient a faire de l'upwind
              PIF = PVAR(II)
     &             +DIFX*DPDX(II)+DIFY*DPDY(II)+DIFZ*DPDZ(II)
              PJF = PVAR(JJ)
     &             +DJFX*DPDX(JJ)+DJFY*DPDY(JJ)+DJFZ*DPDZ(JJ)
C
            ELSE
              WRITE(NFECRA,9000)ISCHCP
              IOK = 1
            ENDIF
C
C
C        BLENDING
C       ----------
C
            PIF = BLENCP*PIF+(1.D0-BLENCP)*PVAR(II)
            PJF = BLENCP*PJF+(1.D0-BLENCP)*PVAR(JJ)
C
C
C        FLUX
C       ------
C
            FLUX = ICONVP*( FLUI*PIF +FLUJ*PJF )
     &           + IDIFFP*VISCF(IFAC)*( PIP -PJP )
C
C
C        ASSEMBLAGE
C       ------------
C
            SMBRP(II) = SMBRP(II) -FLUX
            SMBRP(JJ) = SMBRP(JJ) +FLUX
C
          ENDDO
C        Le call csexit ne doit pas etre dans la boucle, sinon
C         le VPP la devectorise.
          IF(IOK.NE.0) THEN
            CALL CSEXIT (1)
          ENDIF
C
        ELSE
C
          IOK = 0
C VECTORISATION NON FORCEE
          DO IFAC = 1, NFAC
C
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
C
            IIJ = IDIJPF-1+3*(IFAC-1)
C
            DIJPFX = RA(IIJ+1)
            DIJPFY = RA(IIJ+2)
            DIJPFZ = RA(IIJ+3)
C
            POND   = RA(IPOND-1+IFAC)
C
C ON RECALCULE II' ET JJ'
C
            DIIPFX = CDGFAC(1,IFAC) - (XYZCEN(1,II)+
     &               (1.D0-POND) * DIJPFX)
            DIIPFY = CDGFAC(2,IFAC) - (XYZCEN(2,II)+
     &               (1.D0-POND) * DIJPFY)
            DIIPFZ = CDGFAC(3,IFAC) - (XYZCEN(3,II)+
     &               (1.D0-POND) * DIJPFZ)
            DJJPFX = CDGFAC(1,IFAC) -  XYZCEN(1,JJ)+
     &                   POND  * DIJPFX
            DJJPFY = CDGFAC(2,IFAC) -  XYZCEN(2,JJ)+
     &                   POND  * DIJPFY
            DJJPFZ = CDGFAC(3,IFAC) -  XYZCEN(3,JJ)+
     &                   POND  * DIJPFZ
C
            DPXF = 0.5D0*(DPDX(II) + DPDX(JJ))
            DPYF = 0.5D0*(DPDY(II) + DPDY(JJ))
            DPZF = 0.5D0*(DPDZ(II) + DPDZ(JJ))
C
            PIP = PVAR(II)
     &           + IRCFLP*(DPXF*DIIPFX+DPYF*DIIPFY+DPZF*DIIPFZ)
            PJP = PVAR(JJ)
     &           + IRCFLP*(DPXF*DJJPFX+DPYF*DJJPFY+DPZF*DJJPFZ)
C
            FLUI = 0.5D0*( FLUMAS(IFAC) +ABS(FLUMAS(IFAC)) )
            FLUJ = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
C
C
C         CENTRE
C        --------
C
            IF (ISCHCP.EQ.1) THEN
C
              PIF = POND*PIP +(1.D0-POND)*PJP
              PJF = PIF
C
C
C         SECOND ORDER
C        --------------
C
            ELSEIF(ISCHCP.EQ.0) THEN
C
              DIFX = CDGFAC(1,IFAC) - XYZCEN(1,II)
              DIFY = CDGFAC(2,IFAC) - XYZCEN(2,II)
              DIFZ = CDGFAC(3,IFAC) - XYZCEN(3,II)
              DJFX = CDGFAC(1,IFAC) - XYZCEN(1,JJ)
              DJFY = CDGFAC(2,IFAC) - XYZCEN(2,JJ)
              DJFZ = CDGFAC(3,IFAC) - XYZCEN(3,JJ)
C
C     on laisse la reconstruction de PIF et PJF meme si IRCFLP=0
C     sinon cela revient a faire de l'upwind
              PIF = PVAR(II)
     &             +DIFX*DPDX(II)+DIFY*DPDY(II)+DIFZ*DPDZ(II)
              PJF = PVAR(JJ)
     &             +DJFX*DPDX(JJ)+DJFY*DPDY(JJ)+DJFZ*DPDZ(JJ)
C
            ELSE
              WRITE(NFECRA,9000)ISCHCP
              IOK = 1
            ENDIF
C
C
C        BLENDING
C       ----------
C
            PIF = BLENCP*PIF+(1.D0-BLENCP)*PVAR(II)
            PJF = BLENCP*PJF+(1.D0-BLENCP)*PVAR(JJ)
C
C
C        FLUX
C       ------
C
            FLUX = ICONVP*( FLUI*PIF +FLUJ*PJF )
     &           + IDIFFP*VISCF(IFAC)*( PIP -PJP )
C
C
C        ASSEMBLAGE
C       ------------
C
            SMBRP(II) = SMBRP(II) -FLUX
            SMBRP(JJ) = SMBRP(JJ) +FLUX
C
          ENDDO
C        Le call csexit ne doit pas etre dans la boucle, sinon
C         le VPP la devectorise (pour la boucle non vectorisee forcee,
C         ce n'est pas grave, c'est juste pour recopier exactement
C         la precedente)
          IF(IOK.NE.0) THEN
            CALL CSEXIT (1)
          ENDIF
C
        ENDIF
C
C
C
C
C  --> FLUX AVEC TEST DE PENTE (separe pour vectorisation)
C  =============================
C
      ELSE
C
        IF (IVECTI.EQ.1) THEN
C
          IOK = 0
!CDIR NODEP
          DO IFAC = 1, NFAC
C
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
C
            IIJ = IDIJPF-1+3*(IFAC-1)
C
            DIJPFX = RA(IIJ+1)
            DIJPFY = RA(IIJ+2)
            DIJPFZ = RA(IIJ+3)
C
            POND   = RA(IPOND-1+IFAC)
            DIST   = RA(IDIST-1+IFAC)
            SURFAN = RA(ISRFAN-1+IFAC)
C
C ON RECALCULE II' ET JJ'
C
            DIIPFX = CDGFAC(1,IFAC) - (XYZCEN(1,II)+
     &               (1.D0-POND) * DIJPFX)
            DIIPFY = CDGFAC(2,IFAC) - (XYZCEN(2,II)+
     &               (1.D0-POND) * DIJPFY)
            DIIPFZ = CDGFAC(3,IFAC) - (XYZCEN(3,II)+
     &               (1.D0-POND) * DIJPFZ)
            DJJPFX = CDGFAC(1,IFAC) -  XYZCEN(1,JJ)+
     &                   POND  * DIJPFX
            DJJPFY = CDGFAC(2,IFAC) -  XYZCEN(2,JJ)+
     &                   POND  * DIJPFY
            DJJPFZ = CDGFAC(3,IFAC) -  XYZCEN(3,JJ)+
     &                   POND  * DIJPFZ
C
            DPXF = 0.5D0*(DPDX(II) + DPDX(JJ))
            DPYF = 0.5D0*(DPDY(II) + DPDY(JJ))
            DPZF = 0.5D0*(DPDZ(II) + DPDZ(JJ))
C
            PIP = PVAR(II)
     &           + IRCFLP*(DPXF*DIIPFX+DPYF*DIIPFY+DPZF*DIIPFZ)
            PJP = PVAR(JJ)
     &           + IRCFLP*(DPXF*DJJPFX+DPYF*DJJPFY+DPZF*DJJPFZ)
C
            FLUI = 0.5D0*( FLUMAS(IFAC) +ABS(FLUMAS(IFAC)) )
            FLUJ = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
C
C
C         TEST DE PENTE
C        ---------------
C
            TESTI = DPDXA(II)*SURFAC(1,IFAC) +DPDYA(II)*SURFAC(2,IFAC)
     &            + DPDZA(II)*SURFAC(3,IFAC)
            TESTJ = DPDXA(JJ)*SURFAC(1,IFAC) +DPDYA(JJ)*SURFAC(2,IFAC)
     &            + DPDZA(JJ)*SURFAC(3,IFAC)
            TESTIJ= DPDXA(II)*DPDXA(JJ)    +DPDYA(II)*DPDYA(JJ)
     &            + DPDZA(II)*DPDZA(JJ)
C
            IF( FLUMAS(IFAC).GT.0.D0) THEN
              DCC = DPDX(II)*SURFAC(1,IFAC) +DPDY(II)*SURFAC(2,IFAC)
     &            + DPDZ(II)*SURFAC(3,IFAC)
              DDI = TESTI
              DDJ = ( PVAR(JJ)-PVAR(II) )/DIST *SURFAN
            ELSE
              DCC = DPDX(JJ)*SURFAC(1,IFAC) +DPDY(JJ)*SURFAC(2,IFAC)
     &            + DPDZ(JJ)*SURFAC(3,IFAC)
              DDI = ( PVAR(JJ)-PVAR(II) )/DIST *SURFAN
              DDJ = TESTJ
            ENDIF
            TESQCK = DCC**2 -(DDI-DDJ)**2
C
C
C         UPWIND
C        --------
C
CMO          IF( (TESTI*TESTJ).LE.0.D0 .OR. TESTIJ.LE.0.D0 ) THEN
            IF( TESQCK.LE.0.D0 .OR. TESTIJ.LE.0.D0 ) THEN
C
              PIF = PVAR(II)
              PJF = PVAR(JJ)
              INFAC = INFAC+1
C
            ELSE
C
C
C         CENTRE
C        --------
C
              IF (ISCHCP.EQ.1) THEN
C
                PIF = POND*PIP +(1.D0-POND)*PJP
                PJF = PIF
C
C
C         SECOND ORDER
C        --------------
C
              ELSEIF(ISCHCP.EQ.0) THEN
C
                DIFX = CDGFAC(1,IFAC) - XYZCEN(1,II)
                DIFY = CDGFAC(2,IFAC) - XYZCEN(2,II)
                DIFZ = CDGFAC(3,IFAC) - XYZCEN(3,II)
                DJFX = CDGFAC(1,IFAC) - XYZCEN(1,JJ)
                DJFY = CDGFAC(2,IFAC) - XYZCEN(2,JJ)
                DJFZ = CDGFAC(3,IFAC) - XYZCEN(3,JJ)
C
C     on laisse la reconstruction de PIF et PJF meme si IRCFLP=0
C     sinon cela revient a faire de l'upwind
                PIF = PVAR(II)
     &              + DIFX*DPDX(II)+DIFY*DPDY(II)+DIFZ*DPDZ(II)
                PJF = PVAR(JJ)
     &              + DJFX*DPDX(JJ)+DJFY*DPDY(JJ)+DJFZ*DPDZ(JJ)
C
              ELSE
                WRITE(NFECRA,9000)ISCHCP
                IOK = 1
              ENDIF
C
            ENDIF
C
C
C        BLENDING
C       ----------
C
            PIF = BLENCP*PIF+(1.D0-BLENCP)*PVAR(II)
            PJF = BLENCP*PJF+(1.D0-BLENCP)*PVAR(JJ)
C
C
C        FLUX
C       ------
C
            FLUX = ICONVP*( FLUI*PIF +FLUJ*PJF )
     &           + IDIFFP*VISCF(IFAC)*( PIP -PJP )
C
C
C        ASSEMBLAGE
C       ------------
C
            SMBRP(II) = SMBRP(II) -FLUX
            SMBRP(JJ) = SMBRP(JJ) +FLUX
C
          ENDDO
C        Le call csexit ne doit pas etre dans la boucle, sinon
C         le VPP la devectorise.
          IF(IOK.NE.0) THEN
            CALL CSEXIT (1)
          ENDIF
C
        ELSE
C
          IOK = 0
C VECTORISATION NON FORCEE
          DO IFAC = 1, NFAC
C
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
C
            IIJ = IDIJPF-1+3*(IFAC-1)
C
            DIJPFX = RA(IIJ+1)
            DIJPFY = RA(IIJ+2)
            DIJPFZ = RA(IIJ+3)
C
            POND   = RA(IPOND-1+IFAC)
            DIST   = RA(IDIST-1+IFAC)
            SURFAN = RA(ISRFAN-1+IFAC)
C
C ON RECALCULE II' ET JJ'
C
            DIIPFX = CDGFAC(1,IFAC) - (XYZCEN(1,II)+
     &               (1.D0-POND) * DIJPFX)
            DIIPFY = CDGFAC(2,IFAC) - (XYZCEN(2,II)+
     &               (1.D0-POND) * DIJPFY)
            DIIPFZ = CDGFAC(3,IFAC) - (XYZCEN(3,II)+
     &               (1.D0-POND) * DIJPFZ)
            DJJPFX = CDGFAC(1,IFAC) -  XYZCEN(1,JJ)+
     &                   POND  * DIJPFX
            DJJPFY = CDGFAC(2,IFAC) -  XYZCEN(2,JJ)+
     &                   POND  * DIJPFY
            DJJPFZ = CDGFAC(3,IFAC) -  XYZCEN(3,JJ)+
     &                   POND  * DIJPFZ
C
            DPXF = 0.5D0*(DPDX(II) + DPDX(JJ))
            DPYF = 0.5D0*(DPDY(II) + DPDY(JJ))
            DPZF = 0.5D0*(DPDZ(II) + DPDZ(JJ))
C
            PIP = PVAR(II)
     &           + IRCFLP*(DPXF*DIIPFX+DPYF*DIIPFY+DPZF*DIIPFZ)
            PJP = PVAR(JJ)
     &           + IRCFLP*(DPXF*DJJPFX+DPYF*DJJPFY+DPZF*DJJPFZ)
C
            FLUI = 0.5D0*( FLUMAS(IFAC) +ABS(FLUMAS(IFAC)) )
            FLUJ = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
C
C
C         TEST DE PENTE
C        ---------------
C
            TESTI = DPDXA(II)*SURFAC(1,IFAC) +DPDYA(II)*SURFAC(2,IFAC)
     &            + DPDZA(II)*SURFAC(3,IFAC)
            TESTJ = DPDXA(JJ)*SURFAC(1,IFAC) +DPDYA(JJ)*SURFAC(2,IFAC)
     &            + DPDZA(JJ)*SURFAC(3,IFAC)
            TESTIJ= DPDXA(II)*DPDXA(JJ)    +DPDYA(II)*DPDYA(JJ)
     &            + DPDZA(II)*DPDZA(JJ)
C
            IF( FLUMAS(IFAC).GT.0.D0) THEN
              DCC = DPDX(II)*SURFAC(1,IFAC) +DPDY(II)*SURFAC(2,IFAC)
     &            + DPDZ(II)*SURFAC(3,IFAC)
              DDI = TESTI
              DDJ = ( PVAR(JJ)-PVAR(II) )/DIST *SURFAN
            ELSE
              DCC = DPDX(JJ)*SURFAC(1,IFAC) +DPDY(JJ)*SURFAC(2,IFAC)
     &            + DPDZ(JJ)*SURFAC(3,IFAC)
              DDI = ( PVAR(JJ)-PVAR(II) )/DIST *SURFAN
              DDJ = TESTJ
            ENDIF
            TESQCK = DCC**2 -(DDI-DDJ)**2
C
C
C         UPWIND
C        --------
C
CMO          IF( (TESTI*TESTJ).LE.0.D0 .OR. TESTIJ.LE.0.D0 ) THEN
            IF( TESQCK.LE.0.D0 .OR. TESTIJ.LE.0.D0 ) THEN
C
              PIF = PVAR(II)
              PJF = PVAR(JJ)
              INFAC = INFAC+1
C
            ELSE
C
C
C         CENTRE
C        --------
C
              IF (ISCHCP.EQ.1) THEN
C
                PIF = POND*PIP +(1.D0-POND)*PJP
                PJF = PIF
C
C
C         SECOND ORDER
C        --------------
C
              ELSEIF(ISCHCP.EQ.0) THEN
C
                DIFX = CDGFAC(1,IFAC) - XYZCEN(1,II)
                DIFY = CDGFAC(2,IFAC) - XYZCEN(2,II)
                DIFZ = CDGFAC(3,IFAC) - XYZCEN(3,II)
                DJFX = CDGFAC(1,IFAC) - XYZCEN(1,JJ)
                DJFY = CDGFAC(2,IFAC) - XYZCEN(2,JJ)
                DJFZ = CDGFAC(3,IFAC) - XYZCEN(3,JJ)
C
C     on laisse la reconstruction de PIF et PJF meme si IRCFLP=0
C     sinon cela revient a faire de l'upwind
                PIF = PVAR(II)
     &              + DIFX*DPDX(II)+DIFY*DPDY(II)+DIFZ*DPDZ(II)
                PJF = PVAR(JJ)
     &              + DJFX*DPDX(JJ)+DJFY*DPDY(JJ)+DJFZ*DPDZ(JJ)
C
              ELSE
                WRITE(NFECRA,9000)ISCHCP
                IOK = 1
              ENDIF
C
            ENDIF
C
C
C        BLENDING
C       ----------
C
            PIF = BLENCP*PIF+(1.D0-BLENCP)*PVAR(II)
            PJF = BLENCP*PJF+(1.D0-BLENCP)*PVAR(JJ)
C
C
C        FLUX
C       ------
C
            FLUX = ICONVP*( FLUI*PIF +FLUJ*PJF )
     &           + IDIFFP*VISCF(IFAC)*( PIP -PJP )
C
C
C        ASSEMBLAGE
C       ------------
C
            SMBRP(II) = SMBRP(II) -FLUX
            SMBRP(JJ) = SMBRP(JJ) +FLUX
C
          ENDDO
C        Le call csexit ne doit pas etre dans la boucle, sinon
C         le VPP la devectorise (pour la boucle non vectorisee forcee,
C         ce n'est pas grave, c'est juste pour recopier exactement
C         la precedente)
          IF(IOK.NE.0) THEN
            CALL CSEXIT (1)
          ENDIF
C
        ENDIF
C
      ENDIF
C
C
C
      IF(IWARNP.GE.2) THEN
        WRITE(NFECRA,1100)CNOM,INFAC,NFAC
      ENDIF
C
C
C ======================================================================
C ---> ASSEMBLAGE A PARTIR DES FACETTES DE BORD
C ======================================================================
C
C     Lorsque IIFBRU.GT.0, on ne prend pas en compte le flux convectif
C       sur les faces pour lesquelles on dispose par ailleurs d'un
C       flux de C.L.
C
C     Lorsque IIFBRU.LE.0, on adopte le schéma standard (i.e. on prend
C       en compte le flux convectitf).
C
C     Dans les deux cas, si on prend en compte le flux convectif, il
C       pourrait etre utile de vraiment utiliser la condition imposee
C       a la limite, i.e. de ne pas utiliser un schéma upwind (sinon,
C       on ne verra pas certaines C.L.). Actuellement, avec les conditions
C       aux limites proposées, il n'y a que 3 cas pour lesquels on
C       prend en compte le flux convectif ici : paroi, symetrie et
C       sortie supersonique. Dans les deux premiers, le flux est nul.
C       Pour le dernier, un traitement upwind correspond a utiliser
C       effectivement la valeur de bord. On peut donc conserver
C       le code tel quel.
C
C     On n'impose pas le flux convectif sur les faces pour lesquelles
C       il sera imposé par les C.L.
      IF(IIFBRU.GT.0) THEN
C
        IF (IVECTB.EQ.1) THEN
C
!CDIR NODEP
          DO IFAC = 1, NFABOR
C
            II = IFABOR(IFAC)
C
            III = IDIIPB-1+3*(IFAC-1)
            DIIPBX = RA(III+1)
            DIIPBY = RA(III+2)
            DIIPBZ = RA(III+3)
C
            FLUI = 0.5D0*( FLUMAB(IFAC) +ABS(FLUMAB(IFAC)) )
            FLUJ = 0.5D0*( FLUMAB(IFAC) -ABS(FLUMAB(IFAC)) )
C
            PIP = PVAR(II)
     &       +IRCFLP*(DPDX(II)*DIIPBX+DPDY(II)*DIIPBY+DPDZ(II)*DIIPBZ)
C
            PFAC  = INC*COEFAP(IFAC) +COEFBP(IFAC)*PIP
            PFACD = INC*COFAFP(IFAC) +COFBFP(IFAC)*PIP
C
c            FLUX = ICONVP*( FLUI*PVAR(II) +FLUJ*PFAC )
c     &           + IDIFFP*VISCB(IFAC)*( PIP -PFACD )
            FLUX = ICONVP*( FLUI*PVAR(II) +FLUJ*PFAC )
     &                   *DBLE(1-IFRUSB(IFAC))
     &           + IDIFFP*VISCB(IFAC)*( PIP -PFACD )
            SMBRP(II) = SMBRP(II) -FLUX
C
          ENDDO
C
        ELSE
C
          DO IFAC = 1, NFABOR
C
            II = IFABOR(IFAC)
C
            III = IDIIPB-1+3*(IFAC-1)
            DIIPBX = RA(III+1)
            DIIPBY = RA(III+2)
            DIIPBZ = RA(III+3)
C
            FLUI = 0.5D0*( FLUMAB(IFAC) +ABS(FLUMAB(IFAC)) )
            FLUJ = 0.5D0*( FLUMAB(IFAC) -ABS(FLUMAB(IFAC)) )
C
            PIP = PVAR(II)
     &       +IRCFLP*(DPDX(II)*DIIPBX+DPDY(II)*DIIPBY+DPDZ(II)*DIIPBZ)
C
            PFAC  = INC*COEFAP(IFAC) +COEFBP(IFAC)*PIP
            PFACD = INC*COFAFP(IFAC) +COFBFP(IFAC)*PIP
C
c            FLUX = ICONVP*( FLUI*PVAR(II) +FLUJ*PFAC )
c     &           + IDIFFP*VISCB(IFAC)*( PIP -PFACD )
            FLUX = ICONVP*( FLUI*PVAR(II) +FLUJ*PFAC )
     &                   *DBLE(1-IFRUSB(IFAC))
     &           + IDIFFP*VISCB(IFAC)*( PIP -PFACD )
            SMBRP(II) = SMBRP(II) -FLUX
C
          ENDDO
C
        ENDIF
C
C     On ne dispose pas de flux issu des C.L. : traitement std
      ELSE
C
        IF (IVECTB.EQ.1) THEN
C
!CDIR NODEP
          DO IFAC = 1, NFABOR
C
            II = IFABOR(IFAC)
C
            III = IDIIPB-1+3*(IFAC-1)
            DIIPBX = RA(III+1)
            DIIPBY = RA(III+2)
            DIIPBZ = RA(III+3)
C
            FLUI = 0.5D0*( FLUMAB(IFAC) +ABS(FLUMAB(IFAC)) )
            FLUJ = 0.5D0*( FLUMAB(IFAC) -ABS(FLUMAB(IFAC)) )
C
            PIP = PVAR(II)
     &       +IRCFLP*(DPDX(II)*DIIPBX+DPDY(II)*DIIPBY+DPDZ(II)*DIIPBZ)
C
            PFAC  = INC*COEFAP(IFAC) +COEFBP(IFAC)*PIP
            PFACD = INC*COFAFP(IFAC) +COFBFP(IFAC)*PIP
C
            FLUX = ICONVP*( FLUI*PVAR(II) +FLUJ*PFAC )
     &           + IDIFFP*VISCB(IFAC)*( PIP -PFACD )
            SMBRP(II) = SMBRP(II) -FLUX
C
          ENDDO
C
        ELSE
C
          DO IFAC = 1, NFABOR
C
            II = IFABOR(IFAC)
C
            III = IDIIPB-1+3*(IFAC-1)
            DIIPBX = RA(III+1)
            DIIPBY = RA(III+2)
            DIIPBZ = RA(III+3)
C
            FLUI = 0.5D0*( FLUMAB(IFAC) +ABS(FLUMAB(IFAC)) )
            FLUJ = 0.5D0*( FLUMAB(IFAC) -ABS(FLUMAB(IFAC)) )
C
            PIP = PVAR(II)
     &       +IRCFLP*(DPDX(II)*DIIPBX+DPDY(II)*DIIPBY+DPDZ(II)*DIIPBZ)
C
            PFAC  = INC*COEFAP(IFAC) +COEFBP(IFAC)*PIP
            PFACD = INC*COFAFP(IFAC) +COFBFP(IFAC)*PIP
C
            FLUX = ICONVP*( FLUI*PVAR(II) +FLUJ*PFAC )
     &           + IDIFFP*VISCB(IFAC)*( PIP -PFACD )
            SMBRP(II) = SMBRP(II) -FLUX
C
          ENDDO
C
        ENDIF
C
      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(1X,A8,' : CONVECTION EN ',A11,
     &                               ' BLENDING A ',F4.0,' % D''UPWIND')
 1100 FORMAT(1X,A8,' : ',I10,' FACES UPWIND SUR ',
     &                               I10,' FACES INTERNES ')
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS cfbsc2                           ',/,
     &'@    *********                                               ',/,
     &'@     APPEL DE cfbsc2 POUR ',A8 ,' AVEC ISCHCP = ',I10        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
