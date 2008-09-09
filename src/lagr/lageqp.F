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
                        SUBROUTINE LAGEQP
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
     &   DT     , PROPCE , PROPFA , PROPFB ,
     &   VISCF  , VISCB  ,
     &   DAM    , XAM    ,
     &   DRTP   , SMBRS  , ROVSDT ,
     &   FMALA  , FMALB  ,
     &   UL     , VL     , WL     , ALPHAL , PHIA   , PHI    ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     , W10     , W11    , W12 ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC           RESOLUTION D'UNE EQUATION DE POISSON
CFONC
CFONC             div[ALPHA grad(PHI)] = div(ALPHA <Up>)
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
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC           !
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR           !
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
CARGU !  (NFML,NPRFML!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC            !
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !                                      !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR            !
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !                                      !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! XYZCEN       ! TR !  -> ! POINT ASSOCIES AUX VOLUMES DE CONTROL!
CARGU ! (NDIM,NCELET)!    !     !                                      !
CARGU ! SURFAC       ! TR !  -> ! VECTEUR SURFACE DES FACES INTERNES   !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! SURFBO       ! TR !  -> ! VECTEUR SURFACE DES FACES DE BORD    !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! CDGFAC       ! TR !  -> ! CENTRE DE GRAVITE DES FACES INTERNES !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! CDGFBO       ! TR !  -> ! CENTRE DE GRAVITE DES FACES DE BORD  !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS                !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME(NCELET! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! VISCF(NFAC)  ! TR !  -  ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! DAM(NCELET)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DRTP(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! SMBRS(NCELET ! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! ROVSDT(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! FMALA(NFAC)  ! TR !  -  ! FLUX DE MASSE AU FACES INTERNES      !
CARGU ! FMALB(NFABOR)! TR !  -  ! FLUX DE MASSE AU FACES DE BORD       !
CARGU ! UL,VL,WL     ! TR !  -> ! VITESSE LAGRANGIEN                   !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! ALPHAL       ! TR !  -> ! TAUX DE PRESENCE                     !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! PHI , PHIA   ! TR ! <-  ! TERME DE CORRECTION EN N et N-1      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! W1..W9(NCELET! TR !  -  ! TABLEAUX DE TRAVAIL                  !
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
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
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
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION UL(NCELET), VL(NCELET), WL(NCELET)
      DOUBLE PRECISION PHIA(NCELET), PHI(NCELET), ALPHAL(NCELET)
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION FMALA(NFAC) , FMALB(NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET), SMBRS(NCELET)
      DOUBLE PRECISION ROVSDT(NCELET)
      DOUBLE PRECISION W1(NCELET),  W2(NCELET),  W3(NCELET)
      DOUBLE PRECISION W4(NCELET),  W5(NCELET),  W6(NCELET)
      DOUBLE PRECISION W7(NCELET),  W8(NCELET),  W9(NCELET)
      DOUBLE PRECISION W10(NCELET), W11(NCELET), W12(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      INTEGER          IDEBIA, IDEBRA , IFINIA , IFINRA
      INTEGER          IDTVA0, IVAR
      INTEGER          IFAC, IEL
      INTEGER          IPP
      INTEGER          NSWRGP, IMLIGP, IWARNP , IESCAP
      INTEGER          ICONVP, IDIFFP, NDIRCP, IRESLP, NITMAP
      INTEGER          NSWRSP, IRCFLP, ISCHCP, ISSTPP
      INTEGER          IMGRP, NCYMXP, NITMFP
      INTEGER          ICOEFAX,ICOEFAY,ICOEFAZ,ICEFAP
      INTEGER          ICOEFBX,ICOEFBY,ICOEFBZ,ICEFBP
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP, BLENCP, EPSILP, EPSRSP
      DOUBLE PRECISION RELAXP, THETAP
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
      CHAINE = 'Correction pression'
      WRITE(NFECRA,1000) CHAINE(1:19)
C
C=======================================================================
C 2. TERMES SOURCES
C=======================================================================
C
C --> Initialisation
C
      DO IEL = 1, NCEL
        SMBRS(IEL) = 0.D0
      ENDDO
      DO IEL = 1, NCEL
        ROVSDT(IEL) = 0.D0
      ENDDO
      DO IEL = 1, NCEL
        PHI(IEL)  = 0.D0
        PHIA(IEL) = 0.D0
        DRTP(IEL) = 0.D0
      ENDDO
C
C     "VITESSE" DE DIFFUSION FACE
C
        CALL VISCFA
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IMVISF ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   ALPHAL ,
     &   VISCF  , VISCB  ,
     &   RDEVEL , RTUSER , RA     )
C
C CALCUL  de div(Alpha Up) avant correction
C
      DO IEL = 1, NCEL
        W1(IEL) = -UL(IEL)*ALPHAL(IEL)
        W2(IEL) = -VL(IEL)*ALPHAL(IEL)
        W3(IEL) = -WL(IEL)*ALPHAL(IEL)
      ENDDO
C
C --> Calcul du gradient de W1
C     ========================
C
C       On alloue localement 6 tableaux de NFABOR pour le calcul
C         de COEFA et COEFB de W1,W2,W3
C
      ICOEFAX = IDEBRA
      ICOEFBX = ICOEFAX + NFABOR
      ICOEFAY = ICOEFBX + NFABOR
      ICOEFBY = ICOEFAY + NFABOR
      ICOEFAZ = ICOEFBY + NFABOR
      ICOEFBZ = ICOEFAZ + NFABOR
      IFINRA  = ICOEFBZ + NFABOR
      CALL RASIZE ('LAGEQP',IFINRA)
C     ===========
C
      DO IFAC = 1, NFABOR
        IEL = IFABOR(IFAC)
C
        RA(ICOEFAX+IFAC-1) = W1(IEL)
        RA(ICOEFBX+IFAC-1) = ZERO
C
        RA(ICOEFAY+IFAC-1) = W2(IEL)
        RA(ICOEFBY+IFAC-1) = ZERO
C
        RA(ICOEFAZ+IFAC-1) = W3(IEL)
        RA(ICOEFBZ+IFAC-1) = ZERO
C
      ENDDO
C
      CALL DIVERV
C     ===========
     & ( IDEBIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     ,
     &   SMBRS  , W1     , W2     , W3     ,
     &   RA(ICOEFAX) , RA(ICOEFAY) , RA(ICOEFAZ) ,
     &   RA(ICOEFBX) , RA(ICOEFBY) , RA(ICOEFBZ) ,
     &   W4     , W5     , W6     , W7     , W8     ,
     &   W9     , W10    , W11    , W12    ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C
C      On libere la place dans RA
C
      IFINRA = IDEBRA
C
C --> Conditions aux limites sur PHI
C     ==============================
C
C       On alloue localement 2 tableaux de NFABOR pour le calcul
C         de COEFA et COEFB de PHI
C
      IFINIA = IDEBIA
      ICEFAP = IDEBRA
      ICEFBP = ICEFAP + NFABOR
      IFINRA = ICEFBP + NFABOR
      CALL RASIZE ('LAGEQP',IFINRA)
C     ===========
C
      DO IFAC = 1, NFABOR
        IEL = IFABOR(IFAC)
C
        IF ( IA(IITYPF+IFAC-1) .EQ. IENTRE ) THEN
C
C      Flux Nul
C
          RA(ICEFAP+IFAC-1) = ZERO
          RA(ICEFBP+IFAC-1) = 1.D0
C
        ELSE IF ( IA(IITYPF+IFAC-1) .EQ. IPAROI) THEN
C
C      FLux nul
C
          RA(ICEFAP+IFAC-1) = ZERO
          RA(ICEFBP+IFAC-1) = 1.D0
C
        ELSE IF ( IA(IITYPF+IFAC-1) .EQ. IPARUG) THEN
C
C      FLux nul
C
          RA(ICEFAP+IFAC-1) = ZERO
          RA(ICEFBP+IFAC-1) = 1.D0
C
        ELSE IF ( IA(IITYPF+IFAC-1) .EQ. ISYMET) THEN
C
C      FLux nul
C
          RA(ICEFAP+IFAC-1) = ZERO
          RA(ICEFBP+IFAC-1) = 1.D0
C
        ELSE IF ( IA(IITYPF+IFAC-1) .EQ. ISOLIB ) THEN
C
C      Valeur Imposee
C
          RA(ICEFAP+IFAC-1) = PHIA(IEL)
          RA(ICEFBP+IFAC-1) = ZERO
C
        ELSE
          WRITE(NFECRA,1100) IA(IITYPF+IFAC-1)
          CALL CSEXIT (1)
C              ======
        ENDIF
C
      ENDDO
C
C=======================================================================
C 3. RESOLUTION
C=======================================================================
C
C Pas de stationnaire
      IDTVA0 = 0
C Pas de terme de convection
      ICONVP = 0
C Diffusion
      IDIFFP = 1
C Methode de resolution : Gradient conjugue (pas de convection)
      IRESLP = 0
C Valeur par defaut
      NDIRCP = 1
      NITMAP = 1000
      NSWRSP = 2
      NSWRGP = 10000
      IMLIGP = 1
      IRCFLP = 1
      ISCHCP = 1
      ISSTPP = 0
      IMGRP  = 0
      NCYMXP = 100
      NITMFP = 100
      IWARNP = 10
      BLENCP = 0.D0
      EPSILP = 1.D-8
      EPSRSP = 1.D-8
      EPSRGP = 1.D-5
      CLIMGP = 1.5D0
      EXTRAP = 0.D0
      RELAXP = 1.D0
      IESCAP = 0
C
      IPP  = 1
      NOMVAR(IPP) = 'PoissonL'
C
C  IVAR = 0 (indique pour la periodicite de rotation que la variable
C     n'est pas la vitesse ni Rij)
C
      IVAR = 0
C
C On annule les flux de masse
C
      DO IFAC = 1,NFAC
        FMALA(IFAC) = ZERO
      ENDDO
C
      DO IFAC = 1,NFABOR
        FMALB(IFAC) = ZERO
      ENDDO
C
C Dans le cas d'un theta-schema on met theta = 1 (ordre 1)
C
      THETAP = 1.0D0
C
      CALL CODITS
C     ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTVA0 , IVAR   , ICONVP , IDIFFP , IRESLP , NDIRCP , NITMAP ,
     &   IMRGRA , NSWRSP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , IESCAP ,
     &   IMGRP  , NCYMXP , NITMFP , IPP    , IWARNP ,
     &   BLENCP , EPSILP , EPSRSP , EPSRGP , CLIMGP , EXTRAP ,
     &   RELAXP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PHIA   , PHIA   , RA(ICEFAP) , RA(ICEFBP) ,
     &             RA(ICEFAP) , RA(ICEFBP) ,
     &             FMALA       , FMALB       ,
     &   VISCF  , VISCB  , VISCF  , VISCB  ,
     &   ROVSDT , SMBRS  , PHI    ,
     &   DAM    , XAM    , DRTP   ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(/,
     &'   ** RESOLUTION POUR LA VARIABLE ',A19                       ,/,
     &'      ---------------------------                            ',/)
C
 1100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    ERREUR A LA RESOLUTION DE L''EQUATION DE POISSON :      ',/,
     &'@      CONDITIONS AUX LIMITES SUR PHI NON PREVUES (LAGEQP).  ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Contacter l''equipe de developpement.                     ',/,
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
