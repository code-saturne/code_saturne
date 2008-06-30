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
                        SUBROUTINE CLPTRG
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  , ISVHB  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB , RCODCL ,
     &   COEFU  , RIJIPB , COEFA  , COEFB  , UETBOR , VISVDR ,
     &   HBORD  , THBORD ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CONDITIONS LIMITES EN PAROI RUGUEUSE POUR TOUTES LES VARIABLES
CFONC   DE LA PHASE IPHAS
CFONC
CFONC ON SUPPOSE QUE ICODCL(IU) = 6 =>
CFONC              PAROI RUGUEUSE POUR TOUTES LES VARIABLES TURBULENTES
CFONC   (A PRIORI PEU RESTRICTIF EN MONOPHASIQUE)
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
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
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE PHASE                      !
CARGU ! ISVHB        ! E  !  -> ! INDICATEUR DE SAUVEGARDE DES         !
CARGU !              !    !     !  COEFFICIENTS D'ECHANGE AUX BORDS    !
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
CARGU ! ICODCL       ! TE ! <-  ! CODE DE CONDITION LIMITES AUX FACES  !
CARGU !  (NFABOR,NVAR!    !     !  DE BORD                             !
CARGU !              !    !     ! = 1   -> DIRICHLET                   !
CARGU !              !    !     ! = 3   -> DENSITE DE FLUX             !
CARGU !              !    !     ! = 4   -> GLISSEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 5   -> FROTTEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 6   -> RUGOSITE ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 9   -> ENTREE/SORTIE LIBRE (VITESSE!
CARGU !              !    !     !  ENTRANTE EVENTUELLE     BLOQUEE     !
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
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! RCODCL       ! TR ! <-  ! VALEUR DES CONDITIONS AUX LIMITES    !
CARGU !  (NFABOR,NVAR!    !     !  AUX FACES DE BORD                   !
CARGU !              !    !     ! RCODCL(1) = VALEUR DU DIRICHLET      !
CARGU !              !    !     ! RCODCL(2) = VALEUR DU COEF. D'ECHANGE!
CARGU !              !    !     !  EXT. (INFINIE SI PAS D'ECHANGE)     !
CARGU !              !    !     ! RCODCL(3) = VALEUR DE LA DENSITE DE  !
CARGU !              !    !     !  FLUX (NEGATIF SI GAIN) W/m2 OU      !
CARGU !              !    !     !  HAUTEUR DE RUGOSITE (m) si ICODCL=6 !
CARGU !              !    !     ! POUR LES VITESSES (VISTL+VISCT)*GRADU!
CARGU !              !    !     ! POUR LA PRESSION             DT*GRADP!
CARGU !              !    !     ! POUR LES SCALAIRES                   !
CARGU !              !    !     !        CP*(VISCLS+VISCT/SIGMAS)*GRADT!
CARGU ! COEFU        ! TR !  -> ! TAB DE TRAV POUR VALEURS EN IPRIME   !
CARGU ! (NFABOR,3   )!    !     !  DES COMP DE LA VITESSE AU BORD      !
CARGU ! RIJIPB       ! TR !  -> ! TAB DE TRAV POUR VALEURS EN IPRIME   !
CARGU ! (NFABOR,6   )!    !     !  DES RIJ AU BORD                     !
CARGU ! COEFA, COEFB ! TR ! <-> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! UETBOR       ! TR ! <-  ! VITESSE DE FROTTEMENT AU BORD        !
CARGU ! (NFABOR,NPHAS!    !     !  POUR VAN DRIEST EN LES              !
CARGU ! VISVDR(NPHAS)! TR ! <-> ! VISCOSITE DYNAMIQUE DS LES CELLULES  !
CARGU ! (NCELET,NPHAS!    !     !  DE BORD APRES AMORTISST DE V DRIEST !
CARGU ! HBORD        ! TR ! <-  ! COEFFICIENTS D'ECHANGE AUX BORDS     !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! THBORD       ! TR !  -> ! TEMPERATURE AUX BORDS EN I'          !
CARGU ! (NFABOR)     !    !     !    (PLUS EXACTMT : VAR. ENERGETIQUE) !
CARGU ! W1,2,3,4,5,6 ! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU !  (NCELET     !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
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
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "albase.h"
      INCLUDE "radiat.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
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
      INTEGER          IPHAS  , ISVHB
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICODCL(NFABOR,NVAR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION RCODCL(NFABOR,NVAR,3)
      DOUBLE PRECISION COEFU(NFABOR,NDIM), RIJIPB(NFABOR,6)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION UETBOR(NFABOR,NPHAS), VISVDR(NCELET,NPHAS)
      DOUBLE PRECISION HBORD(NFABOR),THBORD(NFABOR)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION W4(NCELET),W5(NCELET),W6(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC, IEL, IVAR, ISOU, II, JJ, KK, LL, ISVHBL
      INTEGER          IHCP, ISCAL
      INTEGER          MODNTL
      INTEGER          IUIPTN
      INTEGER          IUIPH , IVIPH , IWIPH
      INTEGER          IKIPH , IEPIPH, IPHIPH, IFBIPH, IOMGIP
      INTEGER          IR11IP, IR22IP, IR33IP, IR12IP, IR13IP, IR23IP
      INTEGER          ICLU  , ICLV  , ICLW  , ICLK  , ICLEP
      INTEGER          ICL11 , ICL22 , ICL33 , ICL12 , ICL13 , ICL23
      INTEGER          ICLUF , ICLVF , ICLWF , ICLPHI, ICLFB , ICLOMG
      INTEGER          IPCROM, IPCVIS, IPCVST, IPCCP , IPCCV
      INTEGER          ICLVAR, IPCVSL, ICLVAF
      INTEGER          IPH   , IIPH  , IYPLBP
      DOUBLE PRECISION RNX, RNY, RNZ, RXNN
      DOUBLE PRECISION TX, TY, TZ, TXN, TXN0, T2X, T2Y, T2Z
      DOUBLE PRECISION UTAU, UPX, UPY, UPZ, USN
      DOUBLE PRECISION UIPTN, UIPTNF, UIPTMN, UIPTMX
      DOUBLE PRECISION UETMAX, UETMIN, UKMAX, UKMIN, YPLUMX, YPLUMN
      DOUBLE PRECISION UK, UET, YPLUS
      DOUBLE PRECISION SQRCMU, CLSYME, EK
      DOUBLE PRECISION XMUTLM
      DOUBLE PRECISION RCPROD, RCFLUX
      DOUBLE PRECISION CPP, RKL,  PRDTL
      DOUBLE PRECISION HFLUI, HREDUI, HINT, HEXT, PIMP
      DOUBLE PRECISION UND0, DEUXD0
      DOUBLE PRECISION ELOGLO(3,3), ALPHA(6,6)
      DOUBLE PRECISION RCODCX, RCODCY, RCODCZ, RCODSN
      DOUBLE PRECISION VISCLC, VISCTC, ROMC  , DISTBF, SRFBNF, CPSCV
      DOUBLE PRECISION DISTBF0,RUGD,RUGT,YDEP,ACT,UIPTFN,UIPTFX
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
C --- Memoire
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C --- Constantes
      UET = 1.D0
      UTAU = 1.D0
      SQRCMU = SQRT(CMU)
C
      UND0   = 1.D0
      DEUXD0 = 2.D0
C
C --- Variables
      IUIPH  = IU  (IPHAS)
      IVIPH  = IV  (IPHAS)
      IWIPH  = IW  (IPHAS)
      IF(ITYTUR(IPHAS).EQ.2) THEN
        IKIPH  = IK  (IPHAS)
        IEPIPH = IEP (IPHAS)
      ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
        IR11IP = IR11(IPHAS)
        IR22IP = IR22(IPHAS)
        IR33IP = IR33(IPHAS)
        IR12IP = IR12(IPHAS)
        IR13IP = IR13(IPHAS)
        IR23IP = IR23(IPHAS)
        IEPIPH = IEP (IPHAS)
      ELSEIF(ITURB(IPHAS).EQ.50) THEN
        IKIPH  = IK  (IPHAS)
        IEPIPH = IEP (IPHAS)
        IPHIPH = IPHI(IPHAS)
        IFBIPH = IFB (IPHAS)
      ELSEIF(ITURB(IPHAS).EQ.60) THEN
        IKIPH  = IK  (IPHAS)
        IOMGIP = IOMG(IPHAS)
      ENDIF
C
C --- Conditions aux limites
      ICLU   = ICLRTP(IUIPH ,ICOEF)
      ICLV   = ICLRTP(IVIPH ,ICOEF)
      ICLW   = ICLRTP(IWIPH ,ICOEF)
      IF(ITYTUR(IPHAS).EQ.2) THEN
        ICLK   = ICLRTP(IKIPH ,ICOEF)
        ICLEP  = ICLRTP(IEPIPH,ICOEF)
      ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
        ICL11  = ICLRTP(IR11IP,ICOEF)
        ICL22  = ICLRTP(IR22IP,ICOEF)
        ICL33  = ICLRTP(IR33IP,ICOEF)
        ICL12  = ICLRTP(IR12IP,ICOEF)
        ICL13  = ICLRTP(IR13IP,ICOEF)
        ICL23  = ICLRTP(IR23IP,ICOEF)
        ICLEP  = ICLRTP(IEPIPH,ICOEF)
      ELSEIF(ITURB(IPHAS).EQ.50) THEN
        ICLK   = ICLRTP(IKIPH ,ICOEF)
        ICLEP  = ICLRTP(IEPIPH,ICOEF)
        ICLPHI = ICLRTP(IPHIPH,ICOEF)
        ICLFB  = ICLRTP(IFBIPH,ICOEF)
      ELSEIF(ITURB(IPHAS).EQ.60) THEN
        ICLK   = ICLRTP(IKIPH ,ICOEF)
        ICLOMG = ICLRTP(IOMGIP,ICOEF)
      ENDIF
C
      ICLUF  = ICLRTP(IUIPH ,ICOEFF)
      ICLVF  = ICLRTP(IVIPH ,ICOEFF)
      ICLWF  = ICLRTP(IWIPH ,ICOEFF)
C
C --- Grandeurs physiques
      IPCROM = IPPROC(IROM  (IPHAS))
      IPCVIS = IPPROC(IVISCL(IPHAS))
      IPCVST = IPPROC(IVISCT(IPHAS))
      IF(ICP(IPHAS).GT.0) THEN
        IPCCP  = IPPROC(ICP   (IPHAS))
      ELSE
        IPCCP = 0
      ENDIF
C
C --- Compressible
C
      IF ( IPPMOD(ICOMPF) .GE. 0 ) THEN
        IF(ICV(IPHAS).GT.0) THEN
          IPCCV  = IPPROC(ICV   (IPHAS))
        ELSE
          IPCCV = 0
        ENDIF
      ENDIF
C
C --- Rayonnement
      IPH = 0
      IF (IRAYON(IPHAS).GE.1) THEN
        DO IIPH = 1, NPHAST
          IF (IRAPHA(IIPH).EQ.IPHAS) THEN
            IPH = IIPH
          ENDIF
        ENDDO
        IF(IPH.EQ.0) THEN
          WRITE(NFECRA,9000)
          CALL CSEXIT (1)
C         ===========
        ENDIF
      ENDIF
C
C --- Post traitement de Yplus
      IYPLBP = IYPLBR+(IPHAS-1)*NFABOR
C
C
C MIN ET MAX DE LA VITESSE TANGENTIELLE EN PAROI
      UIPTMX = -GRAND
      UIPTMN =  GRAND
C
C MIN ET MAX DE LA VITESSE DE FROTTEMENT EN PAROI
      UETMAX = -GRAND
      UETMIN =  GRAND
      UKMAX  = -GRAND
      UKMIN  =  GRAND
C
C MIN ET MAX DE YPLUS
      YPLUMX = -GRAND
      YPLUMN =  GRAND
C
C COMPTEUR RETOURNEMENT
C
      IUIPTN = 0
C
C
C     En v2f on met directement u=0 donc UIPTMX et UIPTMN vaudront
C     forcement 0
      IF (ITURB(IPHAS).EQ.50) THEN
        UIPTMX = 0.D0
        UIPTMN = 0.D0
      ENDIF
C
C --- Boucle sur les faces : debut
      DO IFAC = 1, NFABOR
C
C --- Test sur la presence d'une condition de paroi rugueuse
        IF( ICODCL(IFAC,IUIPH).EQ.6 ) THEN
C
          IEL = IFABOR(IFAC)
C
C --- Proprietes physiques
          VISCLC = PROPCE(IEL,IPCVIS)
          VISCTC = PROPCE(IEL,IPCVST)
          ROMC   = PROPCE(IEL,IPCROM)
C
C --- Grandeurs geometriques
          DISTBF = RA(IDISTB-1+IFAC)
          SRFBNF = RA(ISRFBN-1+IFAC)
C
C=======================================================================
C 1. REPERE LOCAL
C=======================================================================
C
C
C ---> NORMALE UNITAIRE
C
          RNX = SURFBO(1,IFAC)/SRFBNF
          RNY = SURFBO(2,IFAC)/SRFBNF
          RNZ = SURFBO(3,IFAC)/SRFBNF
C
C ---> PRISE EN COMPTE DE LA VITESSE DE DEFILEMENT
C
          RCODCX = RCODCL(IFAC,IUIPH,1)
          RCODCY = RCODCL(IFAC,IVIPH,1)
          RCODCZ = RCODCL(IFAC,IWIPH,1)
C
C     Si on n'est pas en ALE, on force la vitesse de deplacement
C       de la face a etre tangentielle (et on met a jour RCODCL
C       pour un passage eventuel dans usruet)
          IF (IALE.EQ.0) THEN
            RCODSN = RCODCX*RNX+RCODCY*RNY+RCODCZ*RNZ
            RCODCX = RCODCX -RCODSN*RNX
            RCODCY = RCODCY -RCODSN*RNY
            RCODCZ = RCODCZ -RCODSN*RNZ
            RCODCL(IFAC,IUIPH,1) = RCODCX
            RCODCL(IFAC,IVIPH,1) = RCODCY
            RCODCL(IFAC,IWIPH,1) = RCODCZ
          ENDIF
C
C
C ---> VITESSE TANGENTIELLE RELATIVE
C
          UPX = COEFU(IFAC,1) - RCODCX
          UPY = COEFU(IFAC,2) - RCODCY
          UPZ = COEFU(IFAC,3) - RCODCZ
C
          USN = UPX*RNX+UPY*RNY+UPZ*RNZ
          TX  = UPX -USN*RNX
          TY  = UPY -USN*RNY
          TZ  = UPZ -USN*RNZ
          TXN = SQRT( TX**2 +TY**2 +TZ**2 )
          UTAU= TXN
C
C ---> TANGENTE UNITAIRE
C
          IF( TXN.GE.EPZERO) THEN
C
            TXN0 = 1.D0
C
            TX  = TX/TXN
            TY  = TY/TXN
            TZ  = TZ/TXN
C
          ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
C
C      SI LA VITESSE EST NULLE, LE VECTEUR T EST NORMAL ET QCQUE
C        ON EN A BESOIN POUR LE CHGT DE REPERE DE RIJ
C        ET ON ANNULERA LA VITESSE
C
            TXN0 = 0.D0
C
            IF(ABS(RNY).GE.EPZERO.OR.ABS(RNZ).GE.EPZERO)THEN
              RXNN = SQRT(RNY**2+RNZ**2)
              TX  =  0.D0
              TY  =  RNZ/RXNN
              TZ  = -RNY/RXNN
            ELSEIF(ABS(RNX).GE.EPZERO.OR.ABS(RNZ).GE.EPZERO)THEN
              RXNN = SQRT(RNX**2+RNZ**2)
              TX  =  RNZ/RXNN
              TY  =  0.D0
              TZ  = -RNX/RXNN
            ELSE
              WRITE(NFECRA,1000)IFAC,RNX,RNY,RNZ
              CALL CSEXIT (1)
            ENDIF
C
          ELSE
C
C       SI LA VITESSE EST NULLE ET QU'ON N'EST PAS EN RIJ
C         TX, TY, TZ NE SERT PAS (ON ANNULE LA VITESSE)
C         ET ON LUI DONNE UNE VALEUR BIDON (NULLE PAR EXEMPLE)
C
            TXN0 = 0.D0
C
            TX  = 0.D0
            TY  = 0.D0
            TZ  = 0.D0
C
          ENDIF
C
C ---> ON COMPLETE EVENTUELLEMENT POUR LE RIJ-EPSILON
C
          IF (ITYTUR(IPHAS).EQ.3) THEN
C
C     --> T2 = RN X T (OU X EST LE PRODUIT VECTORIEL)
C
C
            T2X = RNY*TZ - RNZ*TY
            T2Y = RNZ*TX - RNX*TZ
            T2Z = RNX*TY - RNY*TX
C
C     --> MATRICE ORTHOGONALE DE CHANGEMENT DE BASE ELOGLOij
C         (DE LA BASE LOCALE VERS LA BASE GLOBALE)
C
C                            |TX  -RNX  T2X|
C                   ELOGLO = |TY  -RNY  T2Y|
C                            |TZ  -RNZ  T2Z|
C
C         SA TRANSPOSEE ELOGLOt EST SON INVERSE
C
C
            ELOGLO(1,1) =  TX
            ELOGLO(1,2) = -RNX
            ELOGLO(1,3) =  T2X
            ELOGLO(2,1) =  TY
            ELOGLO(2,2) = -RNY
            ELOGLO(2,3) =  T2Y
            ELOGLO(3,1) =  TZ
            ELOGLO(3,2) = -RNZ
            ELOGLO(3,3) =  T2Z
C
C     --> ON CALCULE ALPHA(6,6)
C
C       SOIT f LE CENTRE DE LA FACE DE BORD ET
C            I LE CENTRE DE LA CELLULE CORRESPONDANTE
C
C       EN NOTE RG (RESP RL) INDICE PAR f OU PAR I
C          LE TENSEUR DE REYNOLDS DANS LA BASE GLOBALE (RESP LOCALE)
C
C       LA MATRICE ALPHA APPLIQUEE AU VECTEUR GLOBAL EN I'
C         (RG11,I'|RG22,I'|RG33,I'|RG12,I'|RG13,I'|RG23,I')t
C         DOIT DONNER LES VALEURS A IMPOSER A LA FACE
C         (RG11,f |RG22,f |RG33,f |RG12,f |RG13,f |RG23,f )t
C         AUX CONDITIONS LIMITES DE DIRICHLET PRES (AJOUTEES ENSUITE)
C
C       ON LA DEFINIT EN CALCULANT RG,f EN FONCTION DE RG,I' COMME SUIT
C
C         RG,f = ELOGLO.RL,f.ELOGLOt (PRODUITS MATRICIELS)
C
C                          | RL,I'(1,1)     B*U*.Uk     C*RL,I'(1,3) |
C           AVEC    RL,f = | B*U*.Uk       RL,I'(2,2)       0        |
C                          | C*RL,I'(1,3)     0         RL,I'(3,3)   |
C
C                  AVEC    RL,I = ELOGLOt.RG,I'.ELOGLO
C                          B = 0
C                    ET    C = 0 EN PAROI (1 EN SYMETRIE)
C
C
C
C          ON CALCULE EN FAIT   ELOGLO.PROJECTEUR.ELOGLOt
C
C
            CLSYME=0.D0
            CALL CLCA66 ( CLSYME , ELOGLO , ALPHA )
C           ===========

          ENDIF
C
C=======================================================================
C 2. VITESSES DE FROTTEMENT
C==================================================================
C
C RUGD : rugosite de paroi pour les variables dynamiques
C        seule la valeur stockee pour IU est utilisee
         RUGD=RCODCL(IFAC,IUIPH,3)
C NB: en rugueux la signification de yplus change !
         YPLUS=DISTBF/RUGD
C
C PSEUDO DECALAGE DE LA PAROI de RUGD ((DISTBF+RUGD)/RUGD):
         UET = UTAU/LOG(YPLUS+1.D0)*XKAPPA
C
          IF (IDEUCH(IPHAS).EQ.0) THEN
            UK = UET
          ELSE
C Si IDEUCH=1 ou 2 on calcule uk et uet
C
            IF(ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50
     &           .OR. ITURB(IPHAS).EQ.60) THEN
              EK = RTP(IEL,IKIPH)
            ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
              EK = 0.5D0*
     &               (RTP(IEL,IR11IP)+RTP(IEL,IR22IP)+RTP(IEL,IR33IP))
            ENDIF
C
            UK = CMU025*SQRT(EK)
C
          ENDIF
C
          CALL USRUET
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  , IFAC   , IEL    ,
     &   UK     , UTAU   , YPLUS  ,
     &   UET    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICODCL ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    ,          PROPCE , PROPFA , PROPFB , RCODCL ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
          IF(IDEUCH(IPHAS).EQ.0) THEN
            UK = UET
          ENDIF
C
          UETMAX = MAX(UET,UETMAX)
          UETMIN = MIN(UET,UETMIN)
          UKMAX  = MAX(UK,UKMAX)
          UKMIN  = MIN(UK,UKMIN)
          YPLUMX = MAX(YPLUS,YPLUMX)
          YPLUMN = MIN(YPLUS,YPLUMN)
C
C Sauvegarde de la vitesse de frottement et de la viscosite turbulente
C apres amortissement de van Driest pour la LES
C On n'amortit pas mu_t une seconde fois si on l'a deja fait
C (car une cellule peut avoir plusieurs faces de paroi)
C
          IF(ITYTUR(IPHAS).EQ.4.AND.IDRIES(IPHAS).EQ.1) THEN
            UETBOR(IFAC,IPHAS) = UET
            IF (VISVDR(IEL,IPHAS).LT.-900.D0) THEN
              PROPCE(IEL,IPCVST) = PROPCE(IEL,IPCVST)
C NB amortissement de van Driest a revoir en rugueux :
C    &             *(1.D0-EXP(-YPLUS/CDRIES(IPHAS)))**2
              VISVDR(IEL,IPHAS) = PROPCE(IEL,IPCVST)
              VISCTC = PROPCE(IEL,IPCVST)
            ENDIF
          ENDIF
C
C
C Sauvegarde de yplus si post traite
C
          IF(MOD(IPSTDV,IPSTYP).EQ.0) THEN
            RA(IYPLBP+IFAC-1) = YPLUS
          ENDIF
C
C
C=======================================================================
C 3. CONDITIONS AUX LIMITES SUR LA VITESSE
C=======================================================================
C
C              UIPTN  respecte la production de k
C              de facon conditionnelle    --> Coef RCPROD
C              UIPTNF respecte le flux
C               de facon conditionnelle   --> Coef RCFLUX
C
C On traite tous les modeles de la meme maniere (a revisiter)
          IF (ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.60
     &        .OR.ITURB(IPHAS).EQ.0 .OR.ITURB(IPHAS).EQ.10.OR.
     &           ITYTUR(IPHAS).EQ.3.OR.ITYTUR(IPHAS).EQ.4) THEN
C
C Pseudo decalage de la paroi de la distance RUGD :
            DISTBF0=DISTBF+RUGD
            XMUTLM = XKAPPA*UK*DISTBF0*ROMC
C
            RCPROD = MAX(DISTBF/DISTBF0,
     &              DEUXD0*DISTBF*SQRT(XMUTLM/VISCTC/DISTBF0**2)
     &              -1./(2.+RUGD/DISTBF0))

            RCFLUX = MAX(XMUTLM,VISCTC)/(VISCLC+VISCTC)*DISTBF/DISTBF0

            UIPTN  = MIN(UTAU,MAX(UTAU - UK/XKAPPA*RCPROD,0.D0))
            UIPTNF = UTAU - UET/XKAPPA*RCFLUX
C
C     Clipping :
C     On borne U_f,grad entre 0 et Utau (il y a surement mieux...)
C     - 0    : on interdit le retournement en face de bord, qui est en
C              contradiction avec l'hypothèse de loi log.
C     - Utau : la production turbulente ne peut etre nulle
C     On empeche U_f,flux d'etre negatif
C
            IF (UIPTNF.LT.0) UIPTNF=0

            UIPTMX = MAX(UIPTN,UIPTMX)
            UIPTMN = MIN(UIPTN,UIPTMN)

            IF(UIPTN.LT.-EPZERO) IUIPTN = IUIPTN + 1

            UIPTFX = MAX(UIPTNF,UIPTFX)
            UIPTFN = MIN(UIPTNF,UIPTFN)
            UETMAX = MAX(UET,UETMAX)
            UETMIN = MIN(UET,UETMIN)
            UKMAX = MAX(UK,UKMAX)
            UKMIN = MIN(UK,UKMIN)
C
            COEFA(IFAC,ICLU)   = UIPTN *TX *TXN0
            COEFA(IFAC,ICLV)   = UIPTN *TY *TXN0
            COEFA(IFAC,ICLW)   = UIPTN *TZ *TXN0
            COEFA(IFAC,ICLUF)  = UIPTNF*TX *TXN0
            COEFA(IFAC,ICLVF)  = UIPTNF*TY *TXN0
            COEFA(IFAC,ICLWF)  = UIPTNF*TZ *TXN0
C
            COEFA(IFAC,ICLU)   = COEFA(IFAC,ICLU)  +RCODCX
            COEFA(IFAC,ICLV)   = COEFA(IFAC,ICLV)  +RCODCY
            COEFA(IFAC,ICLW)   = COEFA(IFAC,ICLW)  +RCODCZ
            COEFA(IFAC,ICLUF)  = COEFA(IFAC,ICLUF) +RCODCX
            COEFA(IFAC,ICLVF)  = COEFA(IFAC,ICLVF) +RCODCY
            COEFA(IFAC,ICLWF)  = COEFA(IFAC,ICLWF) +RCODCZ
C
            COEFB(IFAC,ICLU)   = 0.D0
            COEFB(IFAC,ICLV)   = 0.D0
            COEFB(IFAC,ICLW)   = 0.D0
            COEFB(IFAC,ICLUF)  = 0.D0
            COEFB(IFAC,ICLVF)  = 0.D0
            COEFB(IFAC,ICLWF)  = 0.D0
C
C
          ELSEIF(ITURB(IPHAS).EQ.50) THEN
C
C     Avec ces conditions, pas besoin de calculer UIPTMX, UIPTMN
C     et IUIPTN qui sont nuls (valeur d'initialisation)
C
            COEFA(IFAC,ICLU)   = RCODCX
            COEFA(IFAC,ICLV)   = RCODCY
            COEFA(IFAC,ICLW)   = RCODCZ
C
            COEFB(IFAC,ICLU)   = 0.D0
            COEFB(IFAC,ICLV)   = 0.D0
            COEFB(IFAC,ICLW)   = 0.D0
C
          ENDIF
C
C=======================================================================
C 4. CONDITIONS AUX LIMITES SUR K ET EPSILON
C=======================================================================
C
          YDEP = DISTBF*0.5D0+RUGD
          IF (ITYTUR(IPHAS).EQ.2) THEN
C
            COEFA(IFAC,ICLK)   = UK**2/SQRCMU
            COEFB(IFAC,ICLK)   = 0.D0
C
            COEFA(IFAC,ICLEP)  = UK**3/(XKAPPA*YDEP**2)*DISTBF
            COEFB(IFAC,ICLEP)  = 1.D0
C
C=======================================================================
C 5. CONDITIONS AUX LIMITES SUR RIJ ET EPSILON
C=======================================================================
C
          ELSEIF (ITYTUR(IPHAS).EQ.3) THEN
C
C ---> TENSEUR RIJ (PARTIELLEMENT IMPLICITE)
C
            DO ISOU = 1, 6
C
              IF(ISOU.EQ.1) ICLVAR = ICL11
              IF(ISOU.EQ.2) ICLVAR = ICL22
              IF(ISOU.EQ.3) ICLVAR = ICL33
              IF(ISOU.EQ.4) ICLVAR = ICL12
              IF(ISOU.EQ.5) ICLVAR = ICL13
              IF(ISOU.EQ.6) ICLVAR = ICL23
C
              COEFA(IFAC,ICLVAR) = 0.0D0
              COEFB(IFAC,ICLVAR) = 0.0D0
C
            ENDDO
C
            DO ISOU = 1,6
C
              IF(ISOU.EQ.1) THEN
                ICLVAR = ICL11
                JJ = 1
                KK = 1
              ELSE IF(ISOU.EQ.2) THEN
                ICLVAR = ICL22
                JJ = 2
                KK = 2
              ELSE IF(ISOU.EQ.3) THEN
                ICLVAR = ICL33
                JJ = 3
                KK = 3
              ELSE IF(ISOU.EQ.4) THEN
                ICLVAR = ICL12
                JJ = 1
                KK = 2
              ELSE IF(ISOU.EQ.5) THEN
                ICLVAR = ICL13
                JJ = 1
                KK = 3
              ELSE IF(ISOU.EQ.6) THEN
                ICLVAR = ICL23
                JJ = 2
                KK = 3
              ENDIF
C
              IF (ICLPTR(IPHAS).EQ.1) THEN
                DO II = 1, 6
                  IF(II.NE.ISOU) THEN
                    COEFA(IFAC,ICLVAR) = COEFA(IFAC,ICLVAR) +
     &                   ALPHA(ISOU,II) * RIJIPB(IFAC,II)
                  ENDIF
                ENDDO
                COEFB(IFAC,ICLVAR) = ALPHA(ISOU,ISOU)
              ELSE
                DO II = 1, 6
                  COEFA(IFAC,ICLVAR) = COEFA(IFAC,ICLVAR) +
     &                 ALPHA(ISOU,II) * RIJIPB(IFAC,II)
                ENDDO
                COEFB(IFAC,ICLVAR) = 0.D0
              ENDIF
C
              COEFA(IFAC,ICLVAR) = COEFA(IFAC,ICLVAR)  -
     &                           (ELOGLO(JJ,1)*ELOGLO(KK,2)+
     &                            ELOGLO(JJ,2)*ELOGLO(KK,1))*UET*UK
C
            ENDDO
C
C ---> SCALAIRE EPSILON
C      on traite comme en k-eps
C
C
            COEFA(IFAC,ICLEP)  = UK**3/(XKAPPA*YDEP**2)*DISTBF
            COEFB(IFAC,ICLEP) = 1.D0
C
C=======================================================================
C 6. CONDITIONS AUX LIMITES SUR K, EPSILON, F_BARRE ET PHI
C=======================================================================
C
          ELSEIF (ITURB(IPHAS).EQ.50) THEN
C
            COEFA(IFAC,ICLK) = 0.D0
            COEFB(IFAC,ICLK) = 0.D0
            COEFA(IFAC,ICLEP) =
     &           2.0D0*PROPCE(IEL,IPCVIS)/PROPCE(IEL,IPCROM)
     &           *RTP(IEL,IKIPH)/DISTBF**2
            COEFB(IFAC,ICLEP) = 0.D0
            COEFA(IFAC,ICLPHI) = 0.0D0
            COEFB(IFAC,ICLPHI) = 0.0D0
            COEFA(IFAC,ICLFB) = 0.0D0
            COEFB(IFAC,ICLFB) = 0.0D0
C=======================================================================
C 7. CONDITIONS AUX LIMITES SUR K ET OMEGA
C=======================================================================
C
          ELSEIF (ITURB(IPHAS).EQ.60) THEN
C
C     on est toujours hors de la sous couche visqeuse
              COEFA(IFAC,ICLK)   = UK**2/SQRCMU
              COEFB(IFAC,ICLK)   = 0.D0
              COEFA(IFAC,ICLOMG) = DISTBF*4.D0*UK**3*ROMC**2/
     &              (SQRCMU*XKAPPA*VISCLC**2*YPLUS**2)
              COEFB(IFAC,ICLOMG) = 1.D0
C
          ENDIF
C
C=======================================================================
C 8. CONDITIONS AUX LIMITES SUR LES SCALAIRES
C              (AUTRES QUE PRESSION, K, EPSILON, RIJ, VARIANCES)
C    Pour les variances, pas de traitement specifique en paroi : voir
C      condli.
C=======================================================================
C
          IF(NSCAL.GE.1) THEN
C
            DO LL = 1, NSCAL
C
              IF(IPHSCA(LL).EQ.IPHAS.AND.ISCAVR(LL).LE.0) THEN
C
                IVAR = ISCA(LL)
                ICLVAR = ICLRTP(IVAR,ICOEF)
                ICLVAF = ICLRTP(IVAR,ICOEFF)
C
                ISVHBL = 0
                IF(LL.EQ.ISVHB) THEN
                  ISVHBL = ISVHB
                ENDIF
C
                IHCP = 0
                ISCAL = LL
                IF(ISCSTH(ISCAL).EQ.0.OR.ISCSTH(ISCAL).EQ.2
     &                               .OR.ISCSTH(ISCAL).EQ.3) THEN
                  IHCP = 0
                ELSEIF(ABS(ISCSTH(ISCAL)).EQ.1) THEN
                  IF(IPCCP.GT.0) THEN
                    IHCP = 2
                  ELSE
                    IHCP = 1
                  ENDIF
                ENDIF
C
                CPP = 1.D0
                IF(IHCP.EQ.0) THEN
                  CPP = 1.D0
                ELSEIF(IHCP.EQ.2) THEN
                  CPP = PROPCE(IEL,IPCCP )
                ELSEIF(IHCP.EQ.1) THEN
                  CPP = CP0(IPHAS)
                ENDIF
                HINT = CPP
C
                IF(IVISLS(LL).GT.0) THEN
                  IPCVSL = IPPROC(IVISLS(LL))
                ELSE
                  IPCVSL = 0
                ENDIF
                IF (IPCVSL.LE.0) THEN
                  RKL = VISLS0(LL)
                  PRDTL = VISCLC/RKL
                ELSE
                  RKL = PROPCE(IEL,IPCVSL)
                  PRDTL = VISCLC/RKL
                ENDIF
C
C  Compressible : On suppose que le nombre de Pr doit etre
C               defini de la meme façon que l'on resolve
C               en enthalpie ou en energie, soit Mu*Cp/Lambda.
C               Si l'on resout en energie, on a calcule ci-dessus
C               Mu*Cv/Lambda.
C
                IF ( IPPMOD(ICOMPF).GE.0 ) THEN
                  IF(ISCSTH(ISCAL).EQ.3) THEN
                    IF(IPCCP.GT.0) THEN
                      PRDTL = PRDTL*PROPCE(IEL,IPCCP )
                    ELSE
                      PRDTL = PRDTL*CP0(IPHAS)
                    ENDIF
                    IF(IPCCV.GT.0) THEN
                      PRDTL = PRDTL/PROPCE(IEL,IPCCV )
                    ELSE
                      PRDTL = PRDTL/CV0(IPHAS)
                    ENDIF
                  ENDIF
                ENDIF
C
C          CAS TURBULENT
                IF (ITURB(IPHAS).NE.0) THEN
                  IF ( IPPMOD(ICOMPF) .GE. 0 ) THEN
C                 En compressible, pour l'energie LAMBDA/CV+CP/CV*(MUT/SIGMAS)
                    IF(IPCCP.GT.0) THEN
                      CPSCV = PROPCE(IEL,IPPROC(ICP(IPHAS)))
                    ELSE
                      CPSCV = CP0(IPHAS)
                    ENDIF
                    IF(IPCCV.GT.0) THEN
                      CPSCV = CPSCV/PROPCE(IEL,IPPROC(ICV(IPHAS)))
                    ELSE
                      CPSCV = CPSCV/CV0(IPHAS)
                    ENDIF
                    HINT = HINT*(RKL+CPSCV*VISCTC/SIGMAS(LL))/DISTBF
                  ELSE
                    HINT = HINT*(RKL+VISCTC/SIGMAS(LL))/DISTBF
                  ENDIF
C          CAS LAMINAIRE
                ELSE
                  HINT  = HINT*RKL/DISTBF
                ENDIF
C
                IF(ITURB(IPHAS).NE.0.AND.ICODCL(IFAC,IVAR).EQ.6)THEN
C
C Loi rugueuse, on recalcule le coefficient d'echange fluide - paroi
                  RUGT=RCODCL(IFAC,IVAR,3)
                  ACT = XKAPPA/LOG((DISTBF+RUGT)/RUGT)
                  HFLUI = ROMC*CPP*UET*ACT
                ELSE
                  HFLUI = HINT
                ENDIF
C
                IF (ISVHBL .GT. 0) HBORD(IFAC) = HFLUI
C
C
C
C --->  C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE
C
C     Si on a deux types de conditions aux limites (ICLVAR, ICLVAF)
C       il faut que le flux soit traduit par ICLVAF.
C     Si on n'a qu'un type de condition, peu importe (ICLVAF=ICLVAR)
C     Pour le moment, dans cette version compressible, on impose un
C       flux nul pour ICLVAR, lorsqu'il est différent de ICLVAF (cette
C       condition ne sert qu'à la reconstruction des gradients et
C       s'applique à l'energie totale qui inclut l'energie cinétique :

C
                IF( ICODCL(IFAC,IVAR).EQ.6 ) THEN
                  HEXT = RCODCL(IFAC,IVAR,2)
                  PIMP = RCODCL(IFAC,IVAR,1)
                  HREDUI = HINT/HFLUI
                  COEFA(IFAC,ICLVAF) = HEXT*PIMP/(HINT+HEXT*HREDUI)
                  COEFB(IFAC,ICLVAF) = (HINT-(1.D0-HREDUI)*HEXT)/
     &                                 (HINT+HEXT*HREDUI)
                  IF(ICLVAR.NE.ICLVAF) THEN
                    COEFA(IFAC,ICLVAR) = 0.D0
                    COEFB(IFAC,ICLVAR) = 1.D0
                  ENDIF
C
C--> Rayonnement :
C
C      On stocke le coefficient d'echange lambda/distance
C      (ou son equivalent en turbulent) quelle que soit la
C      variable thermique transportee (temperature ou enthalpie)
C      car on l'utilise pour realiser des bilans aux parois qui
C      sont faits en temperature (on cherche la temperature de
C      paroi quelle que soit la variable thermique transportee pour
C      ecrire des eps sigma T4.
C
C     donc :
C
C       lorsque la variable transportee est la temperature
C         ABS(ISCSTH(II)).EQ.1 : RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) = HINT
C         puisque HINT = VISLS * CP / DISTBR
C                      = lambda/distance en W/(m2 K)
C
C       lorsque la variable transportee est l'enthalpie
C         ISCSTH(II).EQ.2 : RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) = HINT*CPR
C         avec
C            IF(IPCCP.GT.0) THEN
C              CPR = PROPCE(IEL,IPCCP )
C            ELSE
C              CPR = CP0(IPHAS)
C            ENDIF
C         puisque HINT = VISLS / DISTBR
C                      = lambda/(CP * distance)
C
C       lorsque la variable transportee est l'energie (compressible)
C         ISCSTH(II).EQ.3 :
C         on procede comme pour l'enthalpie avec CV au lieu de CP
C         (rq : il n'y a pas d'hypothèse, sf en non orthogonal :
C               le flux est le bon et le coef d'echange aussi)
C
C      De meme dans condli.
C
C
C
C               Si on rayonne sur la phase et que
C                  le scalaire est la variable energetique
C
                  IF (IRAYON(IPHAS).GE.1            .AND.
     &                LL           .EQ.ISCALT(IPHAS)     ) THEN
C
C                On calcule le coefficient d'echange en W/(m2 K)
C
C                Si on resout en enthalpie
                    IF(ISCSTH(LL).EQ.2) THEN
C                  Si Cp variable
                      IF(IPCCP.GT.0) THEN
                        RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                       HFLUI*PROPCE(IEL,IPCCP )
                      ELSE
                        RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                       HFLUI*CP0(IPHAS)
                      ENDIF
C
C                  Si on resout en energie (compressible)
                    ELSEIF(ISCSTH(LL).EQ.3) THEN
C                    Si Cv variable
                      IF(IPCCV.GT.0) THEN
                        RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                       HFLUI*PROPCE(IEL,IPCCV )
                      ELSE
                        RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                       HFLUI*CV0(IPHAS)
                      ENDIF
C
C                Si on resout en temperature
                    ELSEIF(ABS(ISCSTH(LL)).EQ.1) THEN
                      RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) = HFLUI
                    ENDIF
C
C                On recupere le flux h(Ti'-Tp) (sortant ou
C                             negatif si gain pour le fluide) en W/m2
C
                    RA(IFCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                   HINT*( (1.D0-COEFB(IFAC,ICLVAF))*THBORD(IFAC)
     &                         - COEFA(IFAC,ICLVAF))
                  ENDIF
C
                ENDIF
C
C --->  C.L DE TYPE FLUX : VOIR CONDLI
C
              ENDIF
C
            ENDDO
C
          ENDIF
C
C
        ENDIF
C --- Test sur la presence d'une condition de paroi vitesse : fin
C
C
C
      ENDDO
C --- Boucle sur les faces : fin
C
      IF (IRANGP.GE.0) THEN
        CALL PARMIN (UIPTMN)
C       ===========
        CALL PARMAX (UIPTMX)
C       ===========
        CALL PARMIN (UETMIN)
C       ===========
        CALL PARMAX (UETMAX)
C       ===========
        CALL PARMIN (UKMIN)
C       ===========
        CALL PARMAX (UKMAX)
C       ===========
        CALL PARMIN (YPLUMN)
C       ===========
        CALL PARMAX (YPLUMX)
C       ===========
        CALL PARCPT (IUIPTN)
C       ===========
      ENDIF
C
C=======================================================================
C 9.  IMPRESSIONS
C=======================================================================
C
C     Remarque : afin de ne pas surcharger les listings dans le cas ou
C       quelques yplus ne sont pas corrects, on ne produit le message
C       qu'aux deux premiers pas de temps ou le message apparait et
C       aux deux derniers pas de temps du calcul, ou si IWARNI est >= 2
C       On indique aussi le numero du dernier pas de temps auquel on
C       a rencontre des yplus hors bornes admissibles
C
      IF(IWARNI(IUIPH).GE.0) THEN
        IF(NTLIST.GT.0) THEN
          MODNTL = MOD(NTCABS,NTLIST)
        ELSEIF(NTLIST.EQ.-1.AND.NTCABS.EQ.NTMABS) THEN
          MODNTL = 0
        ELSE
          MODNTL = 1
        ENDIF
C
        IF (MODNTL.EQ.0 .OR. IWARNI(IUIPH).GE.2) THEN
          WRITE(NFECRA,2010) IPHAS,
     &         UIPTMN,UIPTMX,UETMIN,UETMAX,UKMIN,UKMAX,YPLUMN,YPLUMX,
     &         IUIPTN
        ENDIF
C
      ENDIF
C
C=======================================================================
C 10.  FORMATS
C=======================================================================
C
C
 1000 FORMAT(/,' LA NORMALE A LA FACE DE BORD DE PAROI ',I10,/,
     &         ' EST NULLE ; COORDONNEES : ',3E12.5)
C
 2010 FORMAT(/,
     & 3X,'** CONDITIONS AUX LIMITES EN PAROI RUGUEUSE',/,
     & '   -------------------------------------------',/,
     & '------------------------------------------------------------',/,
     & ' Phase ',I6,'                            Minimum     Maximum',/,
     & '------------------------------------------------------------',/,
     & '   Vitesse rel. en paroi    uiptn : ',2E12.5                 ,/,
     & '   Vitesse de frottement    uet   : ',2E12.5                 ,/,
     & '   Vitesse de frottement    uk    : ',2E12.5                 ,/,
     & '   Distance adim. rugueuse  yplus : ',2E12.5                 ,/,
     & '   ------------------------------------------------------   ',/,
     & '   Nbre de retournements de la vitesse en paroi : ',I10      ,/,
     & '------------------------------------------------------------',
     & /,/)
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT                                 ',/,
     &'@    *********                                               ',/,
     &'@      ERREUR SUR LE NOMBRE DE PHASES DANS CONDLI            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      RETURN
      END
c@z
