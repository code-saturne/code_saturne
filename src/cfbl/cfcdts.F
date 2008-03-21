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
                        SUBROUTINE CFCDTS
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , ICONVP , IDIFFP , IRESLP , NDIRCP , NITMAP ,
     &   IMRGRA , NSWRSP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , IESCAP , IIFBRU ,
     &   IMGRP  , NCYMXP , NITMFP , IPP    , IWARNP ,
     &   BLENCP , EPSILP , EPSRGP , CLIMGP , EXTRAP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , IFRUSB ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PVARA  , COEFAP , COEFBP , COFAFP , COFBFP , FLUMAS , FLUMAB ,
     &   VISCFM , VISCBM , VISCFS , VISCBS ,
     &   ROVSDT , SMBRP  , PVAR   ,
     &   DAM    , XAM    , DAG    , XAG    , DPVAR  ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , SMBINI ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  RESOLUTION SUR UN PAS DE TEMPS D'UNE EQUATION DE CONVECTION
CFONC  /DIFFUSION/TERME SOURCE POUR LA VARIABLE PVAR EN COMPRESSIBLE
CFONC
CFONC  ROVSDT.( PVAR -PVARA )
CFONC         (    ->           --->        )
CFONC    + DIV( RO.U  PVAR -VISC GRAD( PVAR ) ).VOLUME  = SMBR
CFONC         (                             )
CFONC
CFONC ON RESOUT EN FAIT :
CFONC
CFONC  ROVSDT.DPVAR
CFONC         (    ->            --->         )
CFONC    + DIV( RO.U  DPVAR -VISC GRAD( DPVAR ) ).VOLUME  = SMBR
CFONC         (                               )
CFONC  AVEC
CFONC   SMBR = SMBR
CFONC         (    ->n           --->      n  )
CFONC    - DIV( RO.U  DPVAR -VISC GRAD( DPVAR ) ).VOLUME
CFONC         (                               )
CFONC ET DPVAR = INCREMENT VARIABLE PVAR
CFONC
CFONC  Attention, on suppose qu'on arrive ici avec PVAR initialise
CFONC   Y COMPRIS DANS LE HALO EN PERIODICITE (appel a promav pour
CFONC   le calcul de la norme)
CFONC
CFONC
CFONC
CFONC
CFONC
CFONC
c@fonce
C-----------------------------------------------------------------------
C                       ARGUMENTS
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
CARGU ! ICONVP       ! E  !  -> ! INDICATEUR = 1 CONVECTION, 0 SINON   !
CARGU ! IDIFFP       ! E  !  -> ! INDICATEUR = 1 DIFFUSION , 0 SINON   !
CARGU ! NDIRCP       ! E  !  -> ! INDICATEUR = 0 SI DECALAGE DIAGONALE !
CARGU ! IRESLP       ! E  !  -> ! INDICATEUR = 0 GRADCO                !
CARGU !              !    !     !            = 1 JACOBI                !
CARGU !              !    !     !            = 2 BI-CGSTAB             !
CARGU ! IMRGRA       ! E  !  -> ! INDICATEUR = 0 GRADRC 97             !
CARGU !              ! E  !  -> !            = 1 GRADMC 99             !
CARGU ! NSWRSP       ! E  !  -> ! NOMBRE DE SWEEP POUR RECONSTRUCTION  !
CARGU !              !    !     !             DU SECOND MEMBRE         !
CARGU ! NSWRGP       ! E  !  -> ! NOMBRE DE SWEEP POUR RECONSTRUCTION  !
CARGU !              !    !     !             DES GRADIENTS            !
CARGU ! IMLIGP       ! E  !  -> ! METHODE DE LIMITATION DU GRADIENT    !
CARGU !              !    !     !  < 0 PAS DE LIMITATION               !
CARGU !              !    !     !  = 0 A PARTIR DES GRADIENTS VOISINS  !
CARGU !              !    !     !  = 1 A PARTIR DU GRADIENT MOYEN      !
CARGU ! IRCFLP       ! E  !  -> ! INDICATEUR = 1 REC FLUX ; 0 SINON    !
CARGU ! ISCHCP       ! E  !  -> ! INDICATEUR = 1 CENTRE , 0 2ND ORDER  !
CARGU ! ISSTPP       ! E  !  -> ! INDICATEUR = 1 SANS TEST DE PENTE    !
CARGU !              !    !     !            = 0 AVEC TEST DE PENTE    !
CARGU ! IESCAP       ! E  !  -> ! =1 CALCUL DE L'INDICATEUR PREDICTION !
CARGU ! IIFBRU       ! E  !  -> ! POINTEUR FLUX DE BORD RUSANOV        !
CARGU ! IMGRP        ! E  !  -> ! INDICATEUR = 0 PAS DE MGM            !
CARGU !              !    !     !            = 1 SINON                 !
CARGU ! NITMAP       ! E  !  -> ! NOMBRE MAX D'ITER POUR RESOL ITERATIV!
CARGU ! IPP          ! E  !  -> ! NUMERO DE VARIABLE POUR POST         !
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! NFECRL       ! E  !  -> ! UNITE DU FICHIER SORTIE STD          !
CARGU ! BLENCP       ! R  !  -> ! 1 - PROPORTION D'UPWIND              !
CARGU ! EPSILP       ! R  !  -> ! PRECISION POUR RESOL ITER            !
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
CARGU ! PVARA(NCELET ! TR !  -> ! VARIABLE RESOLUE (INSTANT PRECEDENT) !
CARGU ! COEFAP, B    ! TR !  -> ! TABLEAUX DES COND LIM POUR P         !
CARGU !   (NFABOR)   !    !     !  SUR LA NORMALE A LA FACE DE BORD    !
CARGU ! COFAFP, B    ! TR !  -> ! TABLEAUX DES COND LIM POUR LE FLUX DE!
CARGU !   (NFABOR)   !    !     !  DIFFUSION DE P                      !
CARGU ! FLUMAS(NFAC) ! TR !  -> ! FLUX DE MASSE AUX FACES INTERNES     !
CARGU ! FLUMAB(NFABOR! TR !  -> ! FLUX DE MASSE AUX FACES DE BORD      !
CARGU ! VISCFM(NFAC) ! TR !  -> ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU !              !    !     !  POUR LA MATRICE                     !
CARGU ! VISCBM(NFABOR! TR !  -> ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU !              !    !     !  POUR LA MATRICE                     !
CARGU ! VISCFS(NFAC) ! TR !  -> ! IDEM VISCFM POUR SECOND MEMBRE       !
CARGU ! VISCBS(NFABOR! TR !  -> ! IDEM VISCBM POUR SECOND MEMBRE       !
CARGU ! ROVSDT(NCELET! TR !  -> ! RHO*VOLUME/DT                        !
CARGU ! SMBRP(NCELET ! TR ! <-> ! BILAN AU SECOND MEMBRE               !
CARGU ! PVAR (NCELET ! TR ! <-> ! VARIABLE RESOLUE                     !
CARGU ! DAM(NCELET   ! TR ! <-  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU !              !    !     !  ET RESULTAT ESTIMATEUR              !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DAG(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE (MGM)!
CARGU ! XAG(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE (MGM)!
CARGU ! W1...8(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! SMBINI(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "numvar.h"
      INCLUDE "cstnum.h"
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
      INTEGER          IVAR   , ICONVP , IDIFFP , NDIRCP
      INTEGER          NITMAP
      INTEGER          IMRGRA , NSWRSP , NSWRGP , IMLIGP , IRCFLP
      INTEGER          ISCHCP , ISSTPP , IESCAP , IIFBRU , IMGRP
      INTEGER          NCYMXP , NITMFP
      INTEGER          IPP    , IWARNP
      DOUBLE PRECISION BLENCP , EPSILP , EPSRGP , CLIMGP , EXTRAP
      DOUBLE PRECISION THETAP
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IFRUSB(NFABOR)
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
      DOUBLE PRECISION PVARA(NCELET), COEFAP(NFABOR), COEFBP(NFABOR)
      DOUBLE PRECISION                COFAFP(NFABOR), COFBFP(NFABOR)
      DOUBLE PRECISION FLUMAS(NFAC), FLUMAB(NFABOR)
      DOUBLE PRECISION VISCFM(NFAC), VISCBM(NFABOR)
      DOUBLE PRECISION VISCFS(NFAC), VISCBS(NFABOR)
      DOUBLE PRECISION ROVSDT(NCELET), SMBRP(NCELET)
      DOUBLE PRECISION PVAR(NCELET)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC ,2)
      DOUBLE PRECISION DAG(NCELET), XAG(NFAC ,2)
      DOUBLE PRECISION DPVAR(NCELET)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET), W4(NCELET)
      DOUBLE PRECISION W5(NCELET), W6(NCELET), W7(NCELET), W8(NCELET)
      DOUBLE PRECISION SMBINI(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      CHARACTER*8      CNOM
      INTEGER          IDEBIA, IDEBRA, IFINIA, IFINRA
      INTEGER          ISYM,IRESLP,IRESLQ,IPOL,ISQRT,NGR
      INTEGER          INC,ISWEEP,NITERF,ICCOCG,IEL,ICYCLE,NSWMOD
      INTEGER          IPHAS,IDIMTE,ITENSO,IINVPE, IINVPP
      DOUBLE PRECISION RESIDU,RNORM
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C NOMS
      CHAINE = NOMVAR(IPP)
      CNOM   = CHAINE(1:8)
C
C MATRICE A PRIORI SYMETRIQUE ( = 1)
      ISYM  = 1
      IF( ICONVP.GT.0 ) ISYM  = 2
C
C
C METHODE DE RESOLUTION ET DEGRE DU PRECOND DE NEUMANN
C     0 SI CHOIX AUTOMATIQUE GRADCO OU BICGSTAB
C     0 SI CHOIX AUTOMATIQUE JACOBI
C     DONNE PAR IRESLP/1000 SI NON AUTOMATIQUE
      IF (IRESLP.EQ.-1) THEN
        IRESLQ = 0
        IPOL   = 0
        IF( ICONVP.GT.0 ) THEN
          IRESLQ = 1
          IPOL   = 0
        ENDIF
      ELSE
        IRESLQ = MOD(IRESLP,1000)
        IPOL   = (IRESLP-IRESLQ)/1000
      ENDIF
C
C
C PRISE DE SQRT DANS PS
      ISQRT = 1
C
C PRISE EN COMPTE DE LA PERIODICITE
C
C    Initialisation pour test avant promav
      IINVPE = 0
C
      IF(IPERIO.EQ.1) THEN
C
C
C    Par defaut, toutes les periodicites seront traitees dans percom,
C      les variables etant assimilees a des scalaires (meme si ce sont
C      des composantes de vecteurs ou de tenseur)
        IDIMTE = 0
        ITENSO = 0
C
        IINVPE = 1
C
        DO IPHAS = 1, NPHAS
          IF(IVAR.EQ.IU(IPHAS).OR.IVAR.EQ.IV(IPHAS).OR.
     &                            IVAR.EQ.IW(IPHAS).OR.
     &       IVAR.EQ.IR11(IPHAS).OR.IVAR.EQ.IR12(IPHAS).OR.
     &       IVAR.EQ.IR13(IPHAS).OR.IVAR.EQ.IR22(IPHAS).OR.
     &       IVAR.EQ.IR23(IPHAS).OR.IVAR.EQ.IR33(IPHAS)) THEN
C
C    Pour la vitesse et les tensions de Reynolds
C      seules seront echangees les informations sur les faces periodiques
C      de translation dans percom ; on ne touche pas aux informations
C      relatives aux faces de periodicite de rotation.
            IDIMTE = 0
            ITENSO = 1
C
C      Lors de la resolution par increments, on echangera egalement les
C      informations relatives aux faces de periodicite de translation.
C      Pour les faces de periodicite de rotation, l'increment sera
C      annule dans percom (iinvpe=2).
            IINVPE = 2
C
          ENDIF
        ENDDO
C
      ENDIF
C
C=======================================================================
C 1.  CONSTRUCTION MATRICE "SIMPLIFIEE" DE RESOLUTION
C=======================================================================
C
      CALL MATRIX
C     ===========
     & ( NCELET , NCEL   , NFAC   , NFABOR ,
     &   ICONVP , IDIFFP , NDIRCP , ISYM   , NFECRA ,
     &   THETAP ,
     &   IFACEL , IFABOR ,
     &   COEFBP , ROVSDT , FLUMAS , FLUMAB , VISCFM , VISCBM ,
     &   DAM    , XAM    )
C
C      CREATION DE LA HIERARCHIE DE MAILLAGE SI MULTIGRILLE
C
      IF( IMGRP.GT.0 ) THEN
C
        NGR = 0
C
        CALL CLMLGA
C        ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ISYM   , IWARNP , NFECRA , IFINIA , IFINRA , NGR    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DAM    , XAM    , DAG    , XAG    ,
     &   RDEVEL , RTUSER , RA     )
C
        IDEBIA = IFINIA
        IDEBRA = IFINRA
C
      ENDIF
C
C
C=======================================================================
C 2.  BOUCLES SUR LES NON ORTHOGONALITES
C       (A PARTIR DE LA SECONDE ITERATION)
C=======================================================================
C
C
C     AVANT DE BOUCLER SUR LES SWEEP, ON STOCKE LE SECOND MEMBRE SANS
C     RECONSTRUCTION DANS LE TABLEAU AUXILIAIRE SMBINI
C
      DO IEL = 1, NCEL
         SMBINI(IEL) = SMBRP(IEL)
      ENDDO
C
C     On passe toujours dans cfbsc2 avec INC=1
      INC = 1
C     Sauf pour les matrices poids (NSWRSP=-1 ... a priori bloque dans
C                                   verini en compressible)
      IF (NSWRSP.EQ.-1) THEN
        NSWRSP = 1
        INC = 0
      ENDIF
C
C
C  Ca serait bien de le simplifier plus tard aussi ce nombre de sweeps
C  Attention, pour les matrices poids il faut pouvoir ne faire
C     qu'un seul sweep
      NSWMOD = MAX( NSWRSP, 1 )
      DO 100 ISWEEP = 1, NSWMOD
C
C ---> INCREMENTATION ET RECONSTRUCTION DU SECOND MEMBRE
C      ON NE RECALCULE COCG QU'AU PREMIER PASSAGE (PRESQUE)
C
        IF( ISWEEP.EQ.1) THEN
          ICCOCG = 1
C
        ELSE
          ICCOCG = 0
          DO IEL = 1, NCEL
C     SMBINI CONTIENT LES TERMES INSTAT, EN DIV(RHO U) ET SOURCE DE MASSE
C     DU SECOND MEMBRE  MIS A JOUR A CHAQUE SWEEP
            SMBINI(IEL) = SMBINI(IEL) - ROVSDT(IEL)*DPVAR(IEL)
            SMBRP(IEL)  = SMBINI(IEL)
          ENDDO
        ENDIF
C
        CALL CFBSC2
C       ===========
     & ( IDEBIA , IDEBRA ,
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
     &   FLUMAS , FLUMAB , VISCFS , VISCBS ,
     &   SMBRP  ,
C        ------
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
        CALL PRODSC(NCELET,NCEL,ISQRT,SMBRP,SMBRP,RESIDU)
C
C ---> RESIDU DE NORMALISATION CALCULE AU PREMIER SWEEP
C    (NORME C.L +TERMES SOURCES+ TERMES DE NON ORTHOGONALITE)
C
C       Attention, lors de l'appel a promav, ici pour une variable qui
C         n'est pas en increments et qui est supposee initialisee
C         y compris dans le halo.
C         Pour les variables vitesse et les tensions de Reynolds
C         (IINVPE=2), il ne faudra donc pas annuler le halo
C         des periodicites de rotation, mais au contraire le laisser
C         inchange.
C         Pour les autres variables (scalaires) IINVPE=1 permettra de
C         tout echanger, meme si c'est superflu.
        IF( ISWEEP.EQ.1 ) THEN
          IF(IINVPE.EQ.2) THEN
            IINVPP = 3
          ELSE
            IINVPP = IINVPE
          ENDIF
          CALL PROMAV(NCELET,NCEL,NFAC,ISYM,IINVPP,IFACEL,
     &                DAM,XAM,PVAR,W1)
          DO IEL = 1, NCEL
            W1(IEL) = W1(IEL) + SMBRP(IEL)
          ENDDO
          CALL PRODSC(NCELET,NCEL,ISQRT,W1,W1,RNORM)
          RNSMBR(IPP) = RNORM
        ENDIF
C
C ---> RESOLUTION IMPLICITE SUR L'INCREMENT DPVAR
C
        DO IEL = 1, NCEL
          DPVAR(IEL) = 0.D0
        ENDDO
C
        CALL INVERS
C       ===========
     & ( CNOM   , IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ISYM   , IPOL   , IRESLQ , NITMAP , IMGRP  , NGR    ,
     &   NCYMXP , NITMFP ,
     &   IWARNP , NFECRA , NITERF , ICYCLE , IINVPE ,
     &   EPSILP , RNORM  , RESIDU ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DAM    , XAM    , SMBRP  , DPVAR  ,
     &   DAG    , XAG    , W1     , W2     ,
     &   W3     , W4     , W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C
C
        NBIVAR(IPP) = NITERF
        IF(ABS(RNORM).GT.EPZERO) THEN
          RESVAR(IPP) = RESIDU/RNORM
        ELSE
          RESVAR(IPP) = 0.D0
        ENDIF
C
C ---> INCREMENTATION SOLUTION
C
        DO IEL = 1, NCEL
          PVAR(IEL) = PVAR(IEL)+DPVAR(IEL)
        ENDDO
C
C ---> TRAITEMENT DU PARALLELISME
C
      IF(IRANGP.GE.0) CALL PARCOM (PVAR)
C                     ===========
C
C ---> TRAITEMENT DE LA PERIODICITE : SEULE LA PERIODICITE IMPLICITE
C      EST ASSUREE (SCALAIRE ET TRANSLATION DE VECTEUR ET DE TENSEUR)
C
      IF(IPERIO.EQ.1) THEN
        CALL PERCOM
C       ===========
     &  ( IDIMTE , ITENSO ,
     &    PVAR   , PVAR   , PVAR  ,
     &    PVAR   , PVAR   , PVAR  ,
     &    PVAR   , PVAR   , PVAR  )
      ENDIF
C
C ---> TEST DE CONVERGENCE
C
      CALL PRODSC(NCELET,NCEL,ISQRT,SMBRP,SMBRP,RESIDU)
C
      IF( RESIDU.LE.EPSILP*RNORM ) THEN
         IF(IWARNP.GE.1) THEN
            WRITE( NFECRA,1000) CNOM,ISWEEP,RESIDU,RNORM
         ENDIF
         GOTO 200
      ENDIF
C
      IF(IWARNP.GE.3) THEN
         WRITE( NFECRA,1000) CNOM,ISWEEP,RESIDU,RNORM
      ENDIF
C
 100  CONTINUE
C
      IF(IWARNP.GE.2) THEN
         WRITE( NFECRA,1100)CNOM, NSWMOD
      ENDIF
C
C=======================================================================
C 3.  SORTIE OU CALCUL D'ESTIMATEURS POUR LES VITESSES
C       A L'ETAPE DE PREDICTION
C=======================================================================
C
 200  CONTINUE
C
C ---> TEST DE PASSAGE DANS LE CALCUL
C
      IF (IESCAP.GT.0) THEN
C
C ---> CALCUL DE LA CONTRIBUTION COMPOSANTE PAR COMPOSANTE. DE L ESTIMATEUR
C
C
C     SMBINI CONTIENT LES TERMES INSTAT ET EN DIV(U) DU SECOND MEMBRE
C     MIS A JOUR A CHAQUE SWEEP,DONC AU DERNIER, POUR KMAX +1, ON A:
C
        DO IEL = 1,NCEL
          SMBRP(IEL) = SMBINI(IEL) - ROVSDT(IEL)*DPVAR(IEL)
        ENDDO
C
        INC    = 1
        ICCOCG = 1
C
        CALL CFBSC2
C       ===========
     & ( IDEBIA , IDEBRA ,
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
     &   FLUMAS , FLUMAB , VISCFS , VISCBS ,
     &   SMBRP  ,
C        ------
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C     CONTRIBUTION DES NORMES L2 DES DIFFERENTES COMPOSANTES
C       DANS LE TABLEAU DAM QUI EST ICI DISPONIBLE.
C
        DO IEL = 1,NCEL
          DAM(IEL) = (SMBRP(IEL)/ VOLUME(IEL))**2
        ENDDO

      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT (
     & 1X,A8,' : CV-DIF-TS',I5,' IT - RES= ',E12.5,' NORME= ', E12.5)
 1100 FORMAT (
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ',A8 ,' CONVECTION-DIFFUSION-TERMES SOURCES ',/,
     &'@    *********                                               ',/,
     &'@  Nombre d''iterations maximal ',I10   ,' atteint           ',/,
     &'@                                                            '  )
C
C12345678 : CV-DIF-TS 2000 IT - RES= 1234567890234 NORME= 12345678901234
CATTENTION 12345678 : NON CONVERGENCE DU SYSTEME CONV-DIFF-TS
C----
C FIN
C----
C
      END
c@z
