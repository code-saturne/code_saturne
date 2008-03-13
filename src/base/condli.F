c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS                  Code_Saturne version 1.3
C                      ------------------------
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
                        SUBROUTINE CONDLI
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISVHB  , ISVTB  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB , RCODCL ,
     &   COEFA  , COEFB  , UETBOR , VISVDR , HBORD  , THBORD , FRCXT  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   COEFU  , RIJIPB ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC TRADUCTION DES CONDITIONS AUX LIMITES FOURNIES PAR USCLIM.F
CFONC SOUS UNE FORME "SIMPLEMENT" ADMISSIBLE PAR LE SOLVEUR
CFONC
CFONC CETTE TRADUCTION SE PRESENTE SOUS LA FORME D'UNE VALEUR PFAC DE
CFONC  LA VARIABLE P CONSIDEREE A LA FACETTE :
CFONC         PFAC = COEFA +COEFB.P(I)
CFONC P(I) : VALEUR DE LA VARIABLE DANS LA CELLULE FLUIDE ADJACENTE
CFONC
CFONC ATTENTION : SI ON CONSIDERE L'INCREMENT DE LA VARIABLE, LA C.L SE
CFONC REDUIT A : d(PFAC) = COEFB.d(P(I))
CFONC
CFONC CAS PARTICULIER DES VITESSES :
CFONC  --> C.L PEUVENT COUPLER LES 3 COMPOSANTES DE VITESSES
CFONC          (POUR L'INSTANT CE N'EST PAS LE CAS)
CFONC
CFONC   UXFAC = COEFAX +COEFBX  *UX(I) +COEFU(1)*UY(I) +COEFU(2)*UZ(I)
CFONC   UYFAC = COEFAY +COEFU(1)*UX(I) +COEFBY  *UY(I) +COEFU(3)*UZ(I)
CFONC   UZFAC = COEFAZ +COEFU(2)*UX(I) +COEFU(3)*UY(I) +COEFBZ  *UZ(I)
CFONC
CFONC  On dispose du tableau de tri des faces de bord du
CFONC    pas de temps precedent (sauf au premier pas de temps, ou
CFONC    ITRIFB n'a pas ete renseigne)
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
CARGU ! ISVHB        ! E  !  -> ! INDICATEUR DE SAUVEGARDE DES         !
CARGU !              !    !     !  COEFFICIENTS D'ECHANGE AUX BORDS    !
CARGU ! ISVTB        ! E  !  -> ! INDICATEUR DE SAUVEGARDE DES         !
CARGU !              !    !     !  TEMPERATURES AUX BORDS              !
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
CARGU ! ISOSTD       ! TE ! <-  ! INDICATEUR DE SORTIE STANDARD        !
CARGU !    (NFABOR+1)!    !     !  +NUMERO DE LA FACE DE REFERENCE     !
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
CARGU ! COEFA, COEFB ! TR ! <-> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! UETBOR       ! TR ! <-  ! VITESSE DE FROTTEMENT AU BORD        !
CARGU ! (NFABOR,NPHAS!    !     !  POUR VAN DRIEST EN LES              !
CARGU ! VISVDR(NPHAS)! TR ! <-  ! VISCOSITE DYNAMIQUE DS LES CELLULES  !
CARGU ! (NCELET,NPHAS!    !     !  DE BORD APRES AMORTISST DE V DRIEST !
CARGU ! HBORD        ! TR ! <-  ! COEFFICIENTS D'ECHANGE AUX BORDS     !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! THBORD       ! TR ! <-  ! TEMPERATURE AUX BORDS EN I'          !
CARGU ! (NFABOR)     !    !     !    (PLUS EXACTMT : VAR. ENERGETIQUE) !
CARGU ! FRCXT(NCELET,! TR !  -> ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU !   3,NPHAS)   !    !     !  HYDROSTATIQUE                       !
CARGU ! W1,2,3,4,5,6 ! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU !  (NCELET     !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
CARGU ! COEFU        ! TR !  -  ! TAB DE TRAV POUR VALEURS EN IPRIME   !
CARGU ! (NFABOR,3   )!    !     !  DES COMP DE LA VITESSE AU BORD      !
CARGU ! RIJIPB       ! TR !  -  ! TAB DE TRAV POUR VALEURS EN IPRIME   !
CARGU ! (NFABOR,6   )!    !     !  DES RIJ AU BORD                     !
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
      INCLUDE "dimfbr.h"
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "albase.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "radiat.h"
      INCLUDE "parall.h"
      INCLUDE "matiss.h"
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
      INTEGER          ISVHB  , ISVTB
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICODCL(NFABOR,NVAR)
      INTEGER          ISOSTD(NFABOR+1,NPHAS)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NDIMFB,*)
      DOUBLE PRECISION RCODCL(NFABOR,NVAR,3)
      DOUBLE PRECISION FRCXT(NCELET,3,NPHAS)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION UETBOR(NFABOR,NPHAS), VISVDR(NCELET,NPHAS)
      DOUBLE PRECISION HBORD(NFABOR),THBORD(NFABOR)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION W4(NCELET),W5(NCELET),W6(NCELET)
      DOUBLE PRECISION COEFU(NFABOR,NDIM), RIJIPB(NFABOR,6)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC  , IEL   , IVAR  , IPHAS
      INTEGER          ISOU  , II    , III   , IIPH
      INTEGER          IHCP  , ISCAL , ISCAT
      INTEGER          INC   , ICCOCG, IPHYDP
      INTEGER          IOK   , IOK1
      INTEGER          ICODCU
      INTEGER          ISOENT, ISORTI
      INTEGER          ICLSYM, IPATUR, IPATRG, ISVHBL
      INTEGER          IUIPH , IVIPH , IWIPH , IPRIPH
      INTEGER          IKIPH , IEPIPH, IPHIPH, IFBIPH, IOMGIP
      INTEGER          IR11IP, IR22IP, IR33IP, IR12IP, IR13IP, IR23IP
      INTEGER          IPCVIS, IPCVST, IPCCP , IPCVSL, IPCCV
      INTEGER          ICLPR , ICLU  , ICLV  , ICLW  , ICLK  , ICLEP
      INTEGER          ICL11 , ICL22 , ICL33 , ICL12 , ICL13 , ICL23
      INTEGER          ICLUF , ICLVF , ICLWF , ICLPHI, ICLFB , ICLOMG
      INTEGER          ICLVAR, ICLVAF, ICLUMA, ICLVMA, ICLWMA
      INTEGER          IISMPH, IYPLBP
      INTEGER          NSWRGP, IMLIGP, IWARNP, ICLIVA
      INTEGER          IPH
      DOUBLE PRECISION SIGMA , CPP   , RKL
      DOUBLE PRECISION HINT  , HEXT  , PIMP  , XDIS
      DOUBLE PRECISION FLUMBF, VISCLC, VISCTC, DISTBF, SURFBN
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP
      DOUBLE PRECISION RO0IPH, P0IPH , PR0IPH, XXP0, XYP0, XZP0
      DOUBLE PRECISION SRFBNF, RNX   , RNY   , RNZ
      DOUBLE PRECISION UPX   , UPY   , UPZ   , VISTOT
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
C  On a besoin des COEFA et COEFB pour le calcul des gradients
C     pour les cond lim turb en paroi
C   Leur valeur en entree n'est donc pas ecrasee (au premier pas de
C     temps ils sont initialises dans INIVAR a flux nul)
C
C  COEFU sert a stocker la vitesse en I'
C    On l'utilise aussi pour stocker la pression en I' (dans TYPECL), etc
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C  Initialisation du tableau pour stockage de yplus
C     On le remplit dans clptur
C
      IF(MOD(IPSTDV,IPSTYP).EQ.0) THEN
        DO IPHAS = 1, NPHAS
          IYPLBP = IYPLBR+(IPHAS-1)*NFABOR
          DO IFAC = 1, NFABOR
            RA(IYPLBP+IFAC-1) = 0.D0
          ENDDO
        ENDDO
      ENDIF
C
C
C=======================================================================
C 2.  TRAITEMENT DES CONDITIONS DONNES PAR ITYPFB
C=======================================================================
C
C
      IF(NSCAPP.GT.0) THEN
        CALL PPTYCL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , IA(IITRIF) , IA(IITYPF)  , IA(IIZFPP) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RCODCL ,
     &   W1     , W2     , W3     , W4     , W5     , W6     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
      ENDIF
C
      IF (IMATIS.EQ.1) THEN
        CALL MTTYCL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IITYPF)      , IA(IITRIF)      , ICODCL , ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RCODCL , FRCXT  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
      ENDIF
C
      IF (IALE.EQ.1) THEN
        CALL ALTYCL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IITYPF)      , IA(IIALTY)      , ICODCL , IA(IIMPAL)      ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   RCODCL , RA(IXYZN0)      , RA(IDEPAL)      ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
      ENDIF
C
      CALL TYPECL
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IITYPF)      , IA(IITRIF)      , ICODCL , ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RCODCL , FRCXT  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 2.  VERIFICATION DE LA CONSISTANCE DES CL
C=======================================================================
C
      CALL VERICL
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RCODCL ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 4. DISTANCE A LA PAROI ANCIEN MODELE
C=======================================================================
C attention, si on s'en sert pour autre chose, disons le module X
C   bien faire attention dans verini avec quelles options le module
C   X doit alors etre incompatible (perio, parall).
C
      IOK1 = 0
C
      IF(INEEDY.EQ.1.AND.ABS(ICDPAR).EQ.2) THEN
C
        DO IPHAS = 1, NPHAS
C
          IF(IA(IIFAPA(IPHAS)).LE.0) THEN
C
C ON FERA ATTENTION EN PARALLELISME OU PERIODICITE
C    (UNE PAROI PEUT ETRE PLUS PROCHE EN TRAVERSANT UN BORD ...)
C
            DO IEL = 1, NCEL
              W1(IEL) = GRAND
            ENDDO
C
            IUIPH = IU(IPHAS)
C
            DO IFAC = 1, NFABOR
              ICODCU = ICODCL(IFAC,IUIPH)
              IF( ICODCU.EQ.5 .OR. ICODCU.EQ.6 ) THEN
                DO IEL = 1, NCEL
                  XDIS =
     &              (CDGFBO(1,IFAC)-XYZCEN(1,IEL))**2
     &             +(CDGFBO(2,IFAC)-XYZCEN(2,IEL))**2
     &             +(CDGFBO(3,IFAC)-XYZCEN(3,IEL))**2
                  IF(W1(IEL).GT.XDIS) THEN
                    W1(IEL) = XDIS
                    IA(IIFAPA(IPHAS)-1+IEL) = IFAC
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
C
          ENDIF
C
          IOK = 0
          DO IEL = 1, NCEL
            IF(IA(IIFAPA(IPHAS)-1+IEL).LE.0)THEN
              IOK = IOK + 1
            ENDIF
          ENDDO
          IF(IOK.GT.0) THEN
            WRITE(NFECRA,1000) IPHAS,IPHAS,IRIJEC(IPHAS),
     &                               IPHAS,IDRIES(IPHAS)
            IOK1 = 1
          ENDIF
C
        ENDDO
C
      ENDIF
C
C     Normalement, on ne passe pas en parallele ici,
C       mais au cas ou ...
      IF(IRANGP.GE.0) THEN
        CALL PARCPT(IOK1)
      ENDIF
C
      IF(IOK1.NE.0) THEN
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C
C=======================================================================
C 6.  DEBUT DE LA BOUCLE SUR LES PHASES
C        ET REPERAGE DES VARIABLES
C=======================================================================
C
C --- Boucle sur les phases : debut
      DO IPHAS = 1, NPHAS
C
C --- Variables
        RO0IPH = RO0  (IPHAS)
        P0IPH  = P0   (IPHAS)
        PR0IPH = PRED0(IPHAS)
        XXP0   = XYZP0(1,IPHAS)
        XYP0   = XYZP0(2,IPHAS)
        XZP0   = XYZP0(3,IPHAS)
        IPRIPH = IPR (IPHAS)
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
        ICLPR  = ICLRTP(IPRIPH,ICOEF)
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
        IPCVIS = IPPROC(IVISCL(IPHAS))
        IPCVST = IPPROC(IVISCT(IPHAS))
        IF(ICP(IPHAS).GT.0) THEN
          IPCCP  = IPPROC(ICP   (IPHAS))
        ELSE
          IPCCP = 0
        ENDIF
C --- Compressible
        IF ( IPPMOD(ICOMPF).GE.0 ) THEN
          IF(ICV(IPHAS).GT.0) THEN
            IPCCV  = IPPROC(ICV   (IPHAS))
          ELSE
            IPCCV = 0
          ENDIF
        ENDIF
C
C     On determine le numero de phase reduit pour le rayonnement
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
C           ===========
          ENDIF
        ENDIF
C
C
C=======================================================================
C 6.  CONSTRUCTION DE LA TEMPERATURE OU ENTHALPIE
C        AU CENTRE DES FACES DE BORD (OBTENUS PAR Fi + II'.GRAD(Fi))
C
C          POUR LE COUPLAGE SYRTHES
C            THBORD EST UTILISE PAR COUPBO EN SORTIE DE CONDLI
C          POUR LE COUPLAGE AVEC LE MODULE THERMIQUE 1D DE PAROI
C            THBORD EST UTILISE PAR COU1DO EN SORTIE DE CONDLI
C          POUR LE RAYONNEMENT
C            THBORD EST DANS LA BOUCLE POUR CONSTUIRE LE FLUX QUI SERT
C            DANS RAYPAR.
C
C        LE CAS PLUSIEURS PHASES COUPLEES AVEC SYRTHES N'EST PAS PREVU.
C        LE CAS PLUSIEURS PHASES COUPLEES AVEC LE MODULE 1D  N'EST PAS PREVU.
C
C
C        CECI POURRAIT EN PRATIQUE ETRE HORS DE LA BOUCLE.
C
C=======================================================================
C
C  Pour le couplage SYRTHES ou module thermique 1D
C  -----------------------------------------------
C  Ici, on fait une boucle "inutile"  (on ne fait quelque chose
C    que pour ICPSYR(ISCAL) = 1). C'est pour preparer le traitement
C    eventuel de plusieurs temperatures (ie plusieurs couplages
C    SYRTHES a la fois ; noter cependant que meme dans ce cas,
C    une seule temperature sera recue de chaque couplage. En polyph,
C    il faudrait ensuite reconstruire les enthalpies des phases ...
C    plus tard si necessaire).
C  Ici, il ne peut y avoir qu'un seul scalaire avec ICPSYR = 1 et
C    ce uniquement s'il y a effectivement couplage avec SYRTHES
C    (sinon, on s'est arrete dans verini)
C Dans le cas du couplage avec le module 1D, on utilise le scalaire
C    couple avec Syrthes s'il y a couplage, sinon ISCALT(1).
C  La valeur de ISVTB a ete initialisee dans tridim
C    au numero du scalaire couple.
C
C
C  Pour le rayonnement
C  -------------------
C  On calcule la valeur en I' s'il y a une variable
C    thermique sur la phase
C
C
C  On recherche l'unique scalaire qui convient pour la phase courante
C     (ce peut etre T, H, ou E (en compressible))
C
        ISCAT = 0
C
C     Si un scalaire est couple a SYRTHES ou au module 1D
        IF(ISVTB.NE.0) THEN
C       et qu'il est relatif a la phase courante
          IF(IPHSCA(ISVTB).EQ.IPHAS) THEN
C         si ce n'est pas la variable thermique, ca ne va pas.
            IF(ISVTB.NE.ISCALT(IPHAS)) THEN
              WRITE(NFECRA,8000)ISVTB,IPHAS,ISCALT(IPHAS)
              CALL CSEXIT (1)
C             ===========
C         sinon, on calcule le gradient.
            ELSE
              ISCAT = ISVTB
            ENDIF
          ENDIF
       ENDIF
C
C
C     S'il y a du rayonnement sur la phase
C       (il y a forcement une variable energetique)
C       on en calcule le gradient
        IF(IRAYON(IPHAS).GE.1) THEN
          ISCAT = ISCALT(IPHAS)
        ENDIF
C
C     S'il y a un scalaire dont il faut calculer le gradient
C       ... on le calcule.
        IF (ISCAT .GT. 0) THEN
C
          IVAR   = ISCA(ISCAT)
C
          IF (NTCABS.GT.1 .AND. ITBRRB.EQ.1) THEN
C
            INC = 1
            ICCOCG = 1
            NSWRGP = NSWRGR(IVAR)
            IMLIGP = IMLIGR(IVAR)
            IWARNP = IWARNI(IVAR)
            EPSRGP = EPSRGR(IVAR)
            CLIMGP = CLIMGR(IVAR)
            EXTRAP = EXTRAG(IVAR)
            ICLIVA = ICLRTP(IVAR,ICOEF)
            IPHYDP = 0
C
            CALL GRDCEL
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP ,  IPHYDP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   RTPA(1,IVAR)    , COEFA(1,ICLIVA) , COEFB(1,ICLIVA) ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
            DO IFAC = 1 , NFABOR
              IEL = IFABOR(IFAC)
              III = IDIIPB-1+3*(IFAC-1)
              THBORD(IFAC)       =
     &          W1(IEL)*RA(III+1)+W2(IEL)*RA(III+2)+W3(IEL)*RA(III+3)
     &        + RTPA(IEL,IVAR)
            ENDDO
C
          ELSE
C
            DO IFAC = 1 , NFABOR
              IEL = IFABOR(IFAC)
              THBORD(IFAC) = RTPA(IEL,IVAR)
            ENDDO
C
          ENDIF
C
        ENDIF
C
C --- La boucle sur les phases continue
C=======================================================================
C 6.  CONSTRUCTION DE LA VITESSE ET DU TENSEUR DE REYNOLDS
C        AU CENTRE DES FACES DE BORD (OBTENUS PAR Fi + II'.GRAD(Fi))
C        S'IL Y A DES SYMETRIES OU DES PAROIS TURBULENTES
C=======================================================================
C
C ---> INDICATEUR SYMETRIES OU PAROIS TURBULENTES
C
        ICLSYM = 0
        IPATUR = 0
        IPATRG = 0
        DO IFAC = 1, NFABOR
          IF ( ICODCL(IFAC,IUIPH).EQ.4 ) THEN
            ICLSYM = 1
          ELSEIF ( ICODCL(IFAC,IUIPH).EQ.5 ) THEN
            IPATUR = 1
          ELSEIF ( ICODCL(IFAC,IUIPH).EQ.6 ) THEN
            IPATRG = 1
          ENDIF
          IF (ICLSYM.NE.0.AND.IPATUR.NE.0.AND.IPATRG.NE.0 ) GOTO 100
        ENDDO
 100    CONTINUE
C
        IF (IRANGP.GE.0) THEN
           CALL PARCMX(ICLSYM)
           CALL PARCMX(IPATUR)
           CALL PARCMX(IPATRG)
        ENDIF
C
C
C ---> CONSTRUCTION DE LA VITESSE AU CENTRE DES FACES DE BORD
C
        IF (ICLSYM.NE.0.OR.IPATUR.NE.0.OR.IPATRG.NE.0) THEN
C
C
          DO ISOU = 1, 3
C
            IF(ISOU.EQ.1) IVAR = IUIPH
            IF(ISOU.EQ.2) IVAR = IVIPH
            IF(ISOU.EQ.3) IVAR = IWIPH
C
            IF(NTCABS.GT.1) THEN
C
              ICCOCG = 1
              INC    = 1
              NSWRGP = NSWRGR(IVAR)
              IMLIGP = IMLIGR(IVAR)
              IWARNP = IWARNI(IVAR)
              EPSRGP = EPSRGR(IVAR)
              CLIMGP = CLIMGR(IVAR)
              EXTRAP = EXTRAG(IVAR)
              ICLIVA = ICLRTP(IVAR,ICOEF)
              IPHYDP = 0
C
              CALL GRDCEL
C             ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   RTPA(1,IVAR)    , COEFA(1,ICLIVA) , COEFB(1,ICLIVA) ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
              DO IFAC = 1, NFABOR
                IEL = IFABOR(IFAC)
                III = IDIIPB-1+3*(IFAC-1)
                COEFU(IFAC,ISOU) =
     &            W1(IEL)*RA(III+1)+W2(IEL)*RA(III+2)+W3(IEL)*RA(III+3)
     &          + RTPA(IEL,IVAR)
              ENDDO
C
            ELSE
C
              DO IFAC = 1, NFABOR
                IEL = IFABOR(IFAC)
                COEFU(IFAC,ISOU) = RTPA(IEL,IVAR)
              ENDDO
C
            ENDIF
C
          ENDDO
C
        ENDIF
C
C
C ---> CONSTRUCTION DU TENSEUR DE REYNOLDS AU CENTRE DES FACES DE BORD
C
        IF ((ICLSYM.NE.0.OR.IPATUR.NE.0.OR.IPATRG.NE.0)
     &                         .AND.ITYTUR(IPHAS).EQ.3) THEN
C
C
          DO ISOU = 1 , 6
C
            IF(ISOU.EQ.1) IVAR = IR11IP
            IF(ISOU.EQ.2) IVAR = IR22IP
            IF(ISOU.EQ.3) IVAR = IR33IP
            IF(ISOU.EQ.4) IVAR = IR12IP
            IF(ISOU.EQ.5) IVAR = IR13IP
            IF(ISOU.EQ.6) IVAR = IR23IP
C
C
            IF(NTCABS.GT.1.AND.IRIJRB(IPHAS).EQ.1) THEN
C
C CALCUL DU GRADIENT CELLULE DE Rij EN I
C
              INC = 1
              ICCOCG = 1
              NSWRGP = NSWRGR(IVAR)
              IMLIGP = IMLIGR(IVAR)
              IWARNP = IWARNI(IVAR)
              EPSRGP = EPSRGR(IVAR)
              CLIMGP = CLIMGR(IVAR)
              EXTRAP = EXTRAG(IVAR)
              ICLIVA = ICLRTP(IVAR,ICOEF)
              IPHYDP = 0
C
              CALL GRDCEL
C             ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   RTPA(1,IVAR)    , COEFA(1,ICLIVA) , COEFB(1,ICLIVA) ,
     &   W1     , W2     , W3     ,
C        ------   ------   ------
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C CALCUL DE LA VALEUR EN I' DE Rij
C
              DO IFAC = 1 , NFABOR
                IEL = IFABOR(IFAC)
                III = IDIIPB-1+3*(IFAC-1)
                RIJIPB(IFAC,ISOU) =
     &            W1(IEL)*RA(III+1)+W2(IEL)*RA(III+2)+W3(IEL)*RA(III+3)
     &          + RTPA(IEL,IVAR)
              ENDDO
C
C
C   AU PREMIER PAS DE TEMPS, ON NE CONNAIT PAS COEFA ET COEFB
C   (ILS SONT ANNULES DANS CONDLI), LE CALCUL DE RI' EST SIMPLIFIE
C
            ELSE
C
              DO IFAC = 1 , NFABOR
                IEL = IFABOR(IFAC)
                RIJIPB(IFAC,ISOU) = RTPA(IEL,IVAR)
              ENDDO
C
            ENDIF
C
          ENDDO

        ENDIF
C
C --- La boucle sur les phases continues
C=======================================================================
C 7.  TURBULENCE EN PAROI : TOUTES LES VARIABLES CONCERNEES PAR PHASE
C       (U,V,W,K,EPSILON,RIJ,TEMPERATURE)
C=======================================================================
C --- On a besoin de COEFU et de RIJIPB (et THBORD pour le rayonnement)
C --- On suppose que tout scalaire qui dispose de cl de paroi est
C       associe a une phase (qui permet entre autre de calculer yplus)
C
C     On initialise VISVDR a -999.D0.
C     Dans clptur, on amortit la viscosite turbulente sur les cellules
C     de paroi si on a active van Driest. La valeur finale est aussi
C     stockee dans VISVDR.
C     Plus loin, dans vandri, la viscosite sur les cellules
C     de paroi sera amortie une seconde fois. On se sert alors de
C     VISVDR pour lui redonner une valeur correcte.
        IF(ITYTUR(IPHAS).EQ.4.AND.IDRIES(IPHAS).EQ.1) THEN
          DO IEL=1,NCEL
            VISVDR(IEL,IPHAS) = -999.D0
          ENDDO
        ENDIF
C
        IF (IPATUR.NE.0) THEN
C
C Smooth wall laws
          CALL CLPTUR
C         ===========
     & ( IDEBIA , IDEBRA ,
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
C
        ENDIF
C
        IF (IPATRG.NE.0) THEN
C
C Rough wall laws
          CALL CLPTRG
C         ===========
     & ( IDEBIA , IDEBRA ,
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
C
        ENDIF
C
C --- La boucle sur les phases continue
C=======================================================================
C 7.  SYMETRIES POUR LES VECTEURS ET TENSEURS
C       (U,V,W,RIJ)
C=======================================================================
C   On a besoin de COEFU et de RIJIPB
C
        IISMPH = IISYMP     +NFABOR*(IPHAS-1)
        DO IFAC = 1, NFABOR
          IA(IISMPH+IFAC-1) = 1
        ENDDO
C
        IF (ICLSYM.NE.0) THEN
C
          CALL CLSYVT
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , IA(IISMPH) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB , RCODCL ,
     &   COEFU  , RIJIPB , COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF
C
C --- La boucle sur les phases continue
C=======================================================================
C 8.  VITESSE : SORTIE, DIRICHLET, NEUMANN
C=======================================================================
C
C ---> SORTIE : SI FLUX ENTRANT, ON "BLOQUE" A L'INFINI AVAL
C
        ISOENT = 0
        ISORTI = 0
        DO IFAC = 1, NFABOR
C
          FLUMBF = PROPFB(IFAC,IPPROB(IFLUMA(IUIPH)))
C
          IF( ICODCL(IFAC,IUIPH).EQ.9 ) THEN
C
            ISORTI = ISORTI + 1
            IF( FLUMBF.LT.-EPZERO) THEN
              COEFA(IFAC,ICLU) = 0.D0
              COEFB(IFAC,ICLU) = 0.D0
              COEFA(IFAC,ICLV) = 0.D0
              COEFB(IFAC,ICLV) = 0.D0
              COEFA(IFAC,ICLW) = 0.D0
              COEFB(IFAC,ICLW) = 0.D0
              ISOENT = ISOENT + 1
            ELSE
              COEFA(IFAC,ICLU) = 0.D0
              COEFB(IFAC,ICLU) = 1.D0
              COEFA(IFAC,ICLV) = 0.D0
              COEFB(IFAC,ICLV) = 1.D0
              COEFA(IFAC,ICLW) = 0.D0
              COEFB(IFAC,ICLW) = 1.D0
            ENDIF
C
          ENDIF
C
        ENDDO
C
        IF ( MOD(NTCABS,NTLIST).EQ.0 .OR. IWARNI(IU(1)).GE. 0 ) THEN
          IF(ISORTI.GT.0.AND.(IWARNI(IUIPH).GE.2.OR.ISOENT.GT.0)) THEN
            WRITE(NFECRA,3010)IPHAS, ISOENT,ISORTI
          ENDIF
        ENDIF
C
C ---> DIRICHLET ET FLUX
C
        DO II = 1, 3
C
          IF(II.EQ.1) THEN
            IVAR   = IUIPH
            ICLVAR = ICLU
          ELSEIF(II.EQ.2) THEN
            IVAR   = IVIPH
            ICLVAR = ICLV
          ELSEIF(II.EQ.3) THEN
            IVAR   = IWIPH
            ICLVAR = ICLW
          ENDIF
C
          DO IFAC = 1, NFABOR
C
            IEL = IFABOR(IFAC)
C
C --- Proprietes physiques
            VISCLC = PROPCE(IEL,IPCVIS)
            VISCTC = PROPCE(IEL,IPCVST)
C
C --- Grandeurs geometriques
            DISTBF = RA(IDISTB-1+IFAC)
C
            IF (ITYTUR(IPHAS).EQ.3) THEN
              HINT =   VISCLC         /DISTBF
            ELSE
              HINT = ( VISCLC+VISCTC )/DISTBF
            ENDIF
C
C      C.L DE TYPE DIRICHLET
            IF( ICODCL(IFAC,IVAR).EQ.1 ) THEN
              HEXT = RCODCL(IFAC,IVAR,2)
              IF(ABS(HEXT).GT.RINFIN*0.5D0) THEN
                PIMP = RCODCL(IFAC,IVAR,1)
                COEFA(IFAC,ICLVAR) = PIMP
                COEFB(IFAC,ICLVAR) = 0.D0
              ELSE
                PIMP = RCODCL(IFAC,IVAR,1)
                COEFA(IFAC,ICLVAR) = HEXT*PIMP/(HINT +HEXT)
                COEFB(IFAC,ICLVAR) = HINT     /(HINT +HEXT)
              ENDIF
C
C      C.L DE TYPE FLUX
            ELSEIF( ICODCL(IFAC,IVAR).EQ.3 ) THEN
              COEFA(IFAC,ICLVAR) = -RCODCL(IFAC,IVAR,3)/HINT
              COEFB(IFAC,ICLVAR) = 1.D0
            ENDIF
C
          ENDDO
C
        ENDDO
C
C ---> COEFAF ET COEFBF
C       POUR TOUS LES CODES SAUF 4, 5 ET 6 TRAITES SEPAREMENT
C
        DO II = 1, 3
C
          IF(II.EQ.1) THEN
            IVAR   = IUIPH
            ICLVAR = ICLU
            ICLVAF = ICLUF
          ELSEIF(II.EQ.2) THEN
            IVAR   = IVIPH
            ICLVAR = ICLV
            ICLVAF = ICLVF
          ELSEIF(II.EQ.3) THEN
            IVAR   = IWIPH
            ICLVAR = ICLW
            ICLVAF = ICLWF
          ENDIF
C
          IF(ICLVAF.NE.ICLVAR) THEN
            DO IFAC = 1, NFABOR
              IF( ICODCL(IFAC,IVAR).EQ.1.OR.ICODCL(IFAC,IVAR).EQ.3.OR.
     &            ICODCL(IFAC,IVAR).EQ.9                          ) THEN
                COEFA(IFAC,ICLVAF) = COEFA(IFAC,ICLVAR)
                COEFB(IFAC,ICLVAF) = COEFB(IFAC,ICLVAR)
              ENDIF
            ENDDO
          ENDIF
C
        ENDDO
C
C
C --- La boucle sur les phases continue
C=======================================================================
C 9.  PRESSION : DIRICHLET, NEUMANN
C=======================================================================
C
        DO IFAC = 1, NFABOR
C
          IEL = IFABOR(IFAC)
C
C --- Grandeurs geometriques
          DISTBF = RA(IDISTB-1+IFAC)
C
C ON MET UN FLUX EN DT.GRAD P (W/m2) DANS USCLIM
          HINT = DT(IEL)/DISTBF
C
C On doit remodifier la valeur du  Dirichlet de pression de manière
C  à retrouver P*. Car dans typecl.F on a travaillé avec la pression
C totale fournie par l'utilisateur :  Ptotale= P*+ rho.g.r
C En compressible, on laisse RCODCL tel quel
C
C      C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE
C
          IF( ICODCL(IFAC,IPRIPH).EQ.1 ) THEN
            HEXT = RCODCL(IFAC,IPRIPH,2)
            IF ( IPPMOD(ICOMPF).GE.0 ) THEN
              PIMP = RCODCL(IFAC,IPRIPH,1)
            ELSE
              PIMP = RCODCL(IFAC,IPRIPH,1)
     &             - RO0IPH*( GX*(CDGFBO(1,IFAC)-XXP0)
     &                      + GY*(CDGFBO(2,IFAC)-XYP0)
     &                      + GZ*(CDGFBO(3,IFAC)-XZP0) )
     &             + PR0IPH - P0IPH
            ENDIF
            IF( ABS(HEXT).GT.RINFIN*0.5D0 ) THEN
              COEFA(IFAC,ICLPR) = PIMP
              COEFB(IFAC,ICLPR) = 0.D0
            ELSE
              COEFA(IFAC,ICLPR) = HEXT*PIMP/(HINT +HEXT)
              COEFB(IFAC,ICLPR) = HINT     /(HINT +HEXT)
            ENDIF
          ENDIF
C
C      C.L DE TYPE FLUX
          IF( ICODCL(IFAC,IPRIPH).EQ.3 ) THEN
            COEFA(IFAC,ICLPR) = -RCODCL(IFAC,IPRIPH,3)/HINT
            COEFB(IFAC,ICLPR) = 1.D0
          ENDIF
C
        ENDDO
C
C
C --- La boucle sur les phases continue
C=======================================================================
C 10.  K, EPSILON, RIJ, V2F, OMEGA : DIRICHLET, NEUMANN
C=======================================================================
C
C ---> K-EPSILON ET K-OMEGA
C
        IF(ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.60) THEN
C
          DO II = 1, 2
C
C     Pour le k-omega, on met les valeurs sigma_k2 et sigma_w2 car ce terme
C     ne concerne en pratique que les entrees (pas de pb en paroi ou en flux
C     nul)
            IF(II.EQ.1 .AND. ITYTUR(IPHAS).EQ.2) THEN
              IVAR   = IKIPH
              ICLVAR = ICLK
              SIGMA  = SIGMAK
            ELSEIF(II.EQ.1 .AND. ITURB(IPHAS).EQ.60) THEN
              IVAR   = IKIPH
              ICLVAR = ICLK
              SIGMA  = CKWSK2
            ELSEIF (ITYTUR(IPHAS).EQ.2) THEN
              IVAR   = IEPIPH
              ICLVAR = ICLEP
              SIGMA  = SIGMAE
            ELSE
              IVAR   = IOMGIP
              ICLVAR = ICLOMG
              SIGMA  = CKWSW2
            ENDIF
C
            DO IFAC = 1, NFABOR
C
              IEL = IFABOR(IFAC)
C
C --- Proprietes physiques
              VISCLC = PROPCE(IEL,IPCVIS)
              VISCTC = PROPCE(IEL,IPCVST)
              FLUMBF = PROPFB(IFAC,IPPROB(IFLUMA(IKIPH)))
C
C --- Grandeurs geometriques
              DISTBF = RA(IDISTB-1+IFAC)
C
              HINT = (VISCLC+VISCTC/SIGMA)/DISTBF
C
C      C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE
              IF(ICODCL(IFAC,IVAR).EQ.1) THEN
                HEXT = RCODCL(IFAC,IVAR,2)
                PIMP = RCODCL(IFAC,IVAR,1)
                COEFA(IFAC,ICLVAR) = HEXT*PIMP/(HINT +HEXT)
                COEFB(IFAC,ICLVAR) = HINT     /(HINT +HEXT)
C      C.L DE TYPE FLUX
              ELSEIF(ICODCL(IFAC,IVAR).EQ.3)THEN
                COEFA(IFAC,ICLVAR) = -RCODCL(IFAC,IVAR,3)/HINT
                COEFB(IFAC,ICLVAR) = 1.D0
              ENDIF
            ENDDO
C
          ENDDO
C
C ---> RIJ-EPSILON
C         (ATTENTION, PAS DE VISCT)
C
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
C
C   --> RIJ
C
          DO ISOU = 1, 6
C
            IF(ISOU.EQ.1) THEN
              IVAR   = IR11IP
              ICLVAR = ICL11
            ELSEIF(ISOU.EQ.2) THEN
              IVAR   = IR22IP
              ICLVAR = ICL22
            ELSEIF(ISOU.EQ.3) THEN
              IVAR   = IR33IP
              ICLVAR = ICL33
            ELSEIF(ISOU.EQ.4) THEN
              IVAR   = IR12IP
              ICLVAR = ICL12
            ELSEIF(ISOU.EQ.5) THEN
              IVAR   = IR13IP
              ICLVAR = ICL13
            ELSEIF(ISOU.EQ.6) THEN
              IVAR   = IR23IP
              ICLVAR = ICL23
            ENDIF
C
            DO IFAC = 1, NFABOR
C
              IEL = IFABOR(IFAC)
C
C --- Proprietes physiques
              VISCLC = PROPCE(IEL,IPCVIS)
              FLUMBF = PROPFB(IFAC,IPPROB(IFLUMA(IR11IP)))
C
C --- Grandeurs geometriques
              DISTBF = RA(IDISTB-1+IFAC)
C
              IF(ICODCL(IFAC,IVAR).EQ.1) THEN
                HINT = VISCLC/DISTBF
C
                HEXT = RCODCL(IFAC,IVAR,2)
                PIMP = RCODCL(IFAC,IVAR,1)
                COEFA(IFAC,ICLVAR) = HEXT*PIMP/(HINT +HEXT)
                COEFB(IFAC,ICLVAR) = HINT     /(HINT +HEXT)
C
              ELSEIF(ICODCL(IFAC,IVAR).EQ.3)THEN
C
                HINT = VISCLC/DISTBF
C
                COEFA(IFAC,ICLVAR) = -RCODCL(IFAC,IVAR,3)/HINT
                COEFB(IFAC,ICLVAR) = 1.D0
C
              ENDIF
C
            ENDDO
C
          ENDDO
C
C
C   --> EPSILON
C
          IVAR   = IEPIPH
          ICLVAR = ICLEP
C
          DO IFAC = 1, NFABOR
C
            IEL = IFABOR(IFAC)
C
C --- Proprietes physiques
            VISCLC = PROPCE(IEL,IPCVIS)
            FLUMBF = PROPFB(IFAC,IPPROB(IFLUMA(IEPIPH)))
C
C --- Grandeurs geometriques
            DISTBF = RA(IDISTB-1+IFAC)
C
            HINT = VISCLC/DISTBF
C
C      C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE
            IF( ICODCL(IFAC,IVAR).EQ.1) THEN
              HEXT = RCODCL(IFAC,IVAR,2)
              PIMP = RCODCL(IFAC,IVAR,1)
              COEFA(IFAC,ICLVAR) = HEXT*PIMP/(HINT +HEXT)
              COEFB(IFAC,ICLVAR) = HINT     /(HINT +HEXT)
C      C.L DE TYPE FLUX
            ELSEIF(
     &          ICODCL(IFAC,IVAR).EQ.3)THEN
              COEFA(IFAC,ICLVAR) = -RCODCL(IFAC,IVAR,3)/HINT
              COEFB(IFAC,ICLVAR) = 1.D0
            ENDIF
C
          ENDDO
C
C ---> V2F
C
        ELSEIF(ITURB(IPHAS).EQ.50) THEN
C
C   --> K, EPSILON ET PHI
          DO II = 1, 3
C
            IF(II.EQ.1) THEN
              IVAR   = IKIPH
              ICLVAR = ICLK
              SIGMA  = SIGMAK
            ELSEIF(II.EQ.2) THEN
              IVAR   = IEPIPH
              ICLVAR = ICLEP
              SIGMA  = SIGMAE
            ELSE
              IVAR   = IPHIPH
              ICLVAR = ICLPHI
              SIGMA  = SIGMAK
            ENDIF
C
            DO IFAC = 1, NFABOR
C
              IEL = IFABOR(IFAC)
C
C --- Proprietes physiques
              VISCLC = PROPCE(IEL,IPCVIS)
              VISCTC = PROPCE(IEL,IPCVST)
              FLUMBF = PROPFB(IFAC,IPPROB(IFLUMA(IKIPH)))
C
C --- Grandeurs geometriques
              DISTBF = RA(IDISTB-1+IFAC)
C
              HINT = (VISCLC+VISCTC/SIGMA)/DISTBF
C
C      C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE
              IF(ICODCL(IFAC,IVAR).EQ.1) THEN
                HEXT = RCODCL(IFAC,IVAR,2)
                PIMP = RCODCL(IFAC,IVAR,1)
                COEFA(IFAC,ICLVAR) = HEXT*PIMP/(HINT +HEXT)
                COEFB(IFAC,ICLVAR) = HINT     /(HINT +HEXT)
C      C.L DE TYPE FLUX

              ELSEIF(ICODCL(IFAC,IVAR).EQ.3)THEN
                COEFA(IFAC,ICLVAR) = -RCODCL(IFAC,IVAR,3)/HINT
                COEFB(IFAC,ICLVAR) = 1.D0
              ENDIF
            ENDDO
C
          ENDDO
C
C   --> FB
C
          IVAR   = IFBIPH
          ICLVAR = ICLFB
C
          DO IFAC = 1, NFABOR
C
C --- Proprietes physiques
            VISCLC = 1.D0
            FLUMBF = PROPFB(IFAC,IPPROB(IFLUMA(IFBIPH)))
C
C --- Grandeurs geometriques
            DISTBF = RA(IDISTB-1+IFAC)
C
            HINT = VISCLC/DISTBF
C
C      C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE
            IF( ICODCL(IFAC,IVAR).EQ.1) THEN
              HEXT = RCODCL(IFAC,IVAR,2)
              PIMP = RCODCL(IFAC,IVAR,1)
              COEFA(IFAC,ICLVAR) = HEXT*PIMP/(HINT +HEXT)
              COEFB(IFAC,ICLVAR) = HINT     /(HINT +HEXT)
C      C.L DE TYPE FLUX
            ELSEIF(ICODCL(IFAC,IVAR).EQ.3)THEN
              COEFA(IFAC,ICLVAR) = -RCODCL(IFAC,IVAR,3)/HINT
              COEFB(IFAC,ICLVAR) = 1.D0
            ENDIF
C
          ENDDO
C
        ENDIF
C
C --- La boucle sur les phases continue
C=======================================================================
C 11. SCALAIRES (AUTRES QUE PRESSION, K, EPSILON, RIJ, OMEGA, VARIANCES)
C                     : DIRICHLET, NEUMANN
C=======================================================================
C
        IF(NSCAL.GE.1) THEN
C
          DO II = 1, NSCAL
C
            IF(IPHSCA(II).EQ.IPHAS) THEN
C
              IVAR   = ISCA(II)
              ICLVAR = ICLRTP(IVAR,ICOEF)
              ICLVAF = ICLRTP(IVAR,ICOEFF)
C
              ISVHBL = 0
              IF(II.EQ.ISVHB) THEN
                ISVHBL = ISVHB
              ENDIF
C
              IF(IVISLS(II).GT.0) THEN
                IPCVSL = IPPROC(IVISLS(II))
              ELSE
                IPCVSL = 0
              ENDIF
C
C --- Indicateur de prise en compte de Cp ou non
C       (selon si le scalaire (scalaire associe pour une fluctuation)
C        doit etre ou non traite comme une temperature)
C      Si le scalaire est une variance et que le
C        scalaire associe n'est pas resolu, on suppose alors qu'il
C        doit etre traite comme un scalaire passif (defaut IHCP = 0)
              IHCP = 0
              IF(ISCAVR(II).LE.NSCAL) THEN
                IF(ISCAVR(II).GT.0) THEN
                  ISCAL = ISCAVR(II)
                ELSE
                  ISCAL = II
                ENDIF
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
              ENDIF
C
C --- Boucle sur les faces
              DO IFAC = 1, NFABOR
C
                IEL = IFABOR(IFAC)
C
C --- Proprietes physiques
                VISCTC = PROPCE(IEL,IPCVST)
                FLUMBF = PROPFB(IFAC,IPPROB(IFLUMA(IVAR)))
C
C --- Grandeurs geometriques
                DISTBF = RA(IDISTB-1+IFAC)
C
C --- Prise en compte de Cp ou CV
C      (dans le Cas compressible IHCP=0)
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
C --- Viscosite variable ou non
                IF (IPCVSL.LE.0) THEN
                  RKL = VISLS0(II)
                ELSE
                  RKL = PROPCE(IEL,IPCVSL)
                ENDIF
C
C --- Cas turbulent
                IF (ITURB(IPHAS).NE.0) THEN
                  HINT = HINT*(RKL+VISCTC/SIGMAS(II))/DISTBF
C     Cas laminaire
                ELSE
                  HINT  = HINT*RKL/DISTBF
                ENDIF
C
C --->  C.L DE TYPE DIRICHLET AVEC OU SANS COEFFICIENT D'ECHANGE
C
                IF( ICODCL(IFAC,IVAR).EQ.1) THEN
                  HEXT = RCODCL(IFAC,IVAR,2)
                  IF(ABS(HEXT).GE.RINFIN*0.5D0) THEN
                    PIMP = RCODCL(IFAC,IVAR,1)
                    COEFA(IFAC,ICLVAR) = PIMP
                    COEFB(IFAC,ICLVAR) = 0.D0
                  ELSE
                    PIMP = RCODCL(IFAC,IVAR,1)
                    COEFA(IFAC,ICLVAR) = HEXT*PIMP/(HINT+HEXT)
                    COEFB(IFAC,ICLVAR) = HINT     /(HINT+HEXT)
                  ENDIF
C
C     On utilise le Dirichlet pour les calculs de gradients
C       et pour les flux de bord.
C
                  IF(ICLVAF.NE.ICLVAR) THEN
                    COEFA(IFAC,ICLVAF) = COEFA(IFAC,ICLVAR)
                    COEFB(IFAC,ICLVAF) = COEFB(IFAC,ICLVAR)
                  ENDIF
C
C ---> COUPLAGE : on stocke le hint (lambda/d      en temperature,
C                                    lambda/(cp d) en enthalpie,
C                                    lambda/(cv d) en energie)
C
                  IF (ISVHBL .GT. 0) THEN
                    HBORD(IFAC) = HINT
                  ENDIF
C
C
C--> Rayonnement :
C
C      On stocke le coefficient d'echange lambda/distance
C      (ou son equivalent en turbulent) quelle que soit la
C      variable thermique transportee (temperature ou enthalpie)
C      car on l'utilise pour realiser des bilans aux parois qui
C      sont faits en temperature (on cherche la temperature de
C      paroi quelle que soit la variable thermique transportee pour
C      ecrire des eps sigma T4).
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
C       lorsque la variable transportee est l'energie
C         ISCSTH(II).EQ.3 :
C         on procede comme pour l'enthalpie avec CV au lieu de CP
C         (rq : il n'y a pas d'hypothèse, sf en non orthogonal :
C               le flux est le bon et le coef d'echange aussi)
C
C      De meme plus bas et de meme dans clptur.
C
C               Si on rayonne sur la phase et que
C                  le scalaire est la variable energetique
C
                  IF (IRAYON(IPHAS).GE.1            .AND.
     &                II           .EQ.ISCALT(IPHAS)     ) THEN
C
C                On calcule le coefficient d'echange en W/(m2 K)
C
C                  Si on resout en enthalpie
                    IF(ISCSTH(II).EQ.2) THEN
C                    Si Cp variable
                      IF(IPCCP.GT.0) THEN
                        RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                       HINT*PROPCE(IEL,IPCCP )
                      ELSE
                        RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                       HINT*CP0(IPHAS)
                      ENDIF
C                  Si on resout en energie (compressible)
                    ELSEIF(ISCSTH(II).EQ.3) THEN
C                    Si Cv variable
                      IF(IPCCV.GT.0) THEN
                        RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                       HINT*PROPCE(IEL,IPCCV )
                      ELSE
                        RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                       HINT*CV0(IPHAS)
                      ENDIF
C                  Si on resout en temperature
                    ELSEIF(ABS(ISCSTH(II)).EQ.1) THEN
                      RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) = HINT
                    ENDIF
C
C                On recupere le flux h(Ti'-Tp) (sortant ou
C                             negatif si gain pour le fluide) en W/m2
C
                    RA(IFCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                   HINT*( (1.D0-COEFB(IFAC,ICLVAF))*THBORD(IFAC)
     &                         - COEFA(IFAC,ICLVAF))
C
                  ENDIF
C
C --->  C.L DE TYPE FLUX
C
                ELSEIF(ICODCL(IFAC,IVAR).EQ.3)THEN
                  COEFA(IFAC,ICLVAF)  = -RCODCL(IFAC,IVAR,3)/HINT
                  COEFB(IFAC,ICLVAF)  = 1.D0
                  IF(ICLVAR.NE.ICLVAF) THEN
                    COEFA(IFAC,ICLVAR)  = 0.D0
                    COEFB(IFAC,ICLVAR)  = 1.D0
                  ENDIF
                  IF (ISVHBL .GT. 0) HBORD(IFAC) = HINT
C
C
C--> Rayonnement :
C
                  IF (IRAYON(IPHAS).GE.1            .AND.
     &                II           .EQ.ISCALT(IPHAS)     ) THEN
C
C                On calcule le coefficient d'echange en W/(m2 K)
C
C                Si on resout en enthalpie
                    IF(ISCSTH(II).EQ.2) THEN
C                  Si Cp variable
                      IF(IPCCP.GT.0) THEN
                        RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                       HINT*PROPCE(IEL,IPCCP )
                      ELSE
                        RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                       HINT*CP0(IPHAS)
                      ENDIF
                    ELSEIF(ISCSTH(II).EQ.3) THEN
C                    Si Cv variable
                      IF(IPCCV.GT.0) THEN
                        RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                       HINT*PROPCE(IEL,IPCCV )
                      ELSE
                        RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                       HINT*CV0(IPHAS)
                      ENDIF
C                Si on resout en temperature
                    ELSEIF(ABS(ISCSTH(II)).EQ.1) THEN
                      RA(IHCONV-1+IFAC+NFABOR*(IPH-1)) = HINT
                    ENDIF
C
C              On recupere le flux h(Ti'-Tp) (sortant ou
C                             negatif si gain pour le fluide)
C
                    RA(IFCONV-1+IFAC+NFABOR*(IPH-1)) =
     &                   RCODCL(IFAC,IVAR,3)
                  ENDIF
C
                ENDIF
C
              ENDDO
C
            ENDIF
C
          ENDDO
C
        ENDIF
C
      ENDDO
C --- Boucle sur les phases : fin
C=======================================================================
C 13.  VITESSE DE MAILLAGE EN ALE : DIRICHLET, NEUMANN
C=======================================================================
C (les conditions de glissement on ete traitees dans ALTYCL
C
      IF (IALE.EQ.1) THEN
C
        ICLUMA = ICLRTP(IUMA ,ICOEF)
        ICLVMA = ICLRTP(IVMA ,ICOEF)
        ICLWMA = ICLRTP(IWMA ,ICOEF)
C
        DO IFAC = 1, NFABOR
C
          IEL = IFABOR(IFAC)
          DISTBF = RA(IDISTB-1+IFAC)
          SURFBN = RA(ISRFBN-1+IFAC)**2
          IF (IORTVM.EQ.1) THEN
            HINT = PROPCE(IEL,IPPROC(IVISMA(1)))/DISTBF
          ELSE
            HINT = ( PROPCE(IEL,IPPROC(IVISMA(1)))*SURFBO(1,IFAC)**2
     &             + PROPCE(IEL,IPPROC(IVISMA(2)))*SURFBO(2,IFAC)**2
     &             + PROPCE(IEL,IPPROC(IVISMA(3)))*SURFBO(3,IFAC)**2 )
     &           /DISTBF/SURFBN
          ENDIF
C
          DO II = 1, 3
            IF (II.EQ.1) IVAR = IUMA
            IF (II.EQ.2) IVAR = IVMA
            IF (II.EQ.3) IVAR = IWMA
            ICLVAR = ICLRTP(IVAR,ICOEF)
C
C      C.L DE TYPE DIRICHLET
            IF( ICODCL(IFAC,IVAR).EQ.1 ) THEN
              HEXT = RCODCL(IFAC,IVAR,2)
              IF(ABS(HEXT).GT.RINFIN*0.5D0) THEN
                PIMP = RCODCL(IFAC,IVAR,1)
                COEFA(IFAC,ICLVAR) = PIMP
                COEFB(IFAC,ICLVAR) = 0.D0
              ELSE
                PIMP = RCODCL(IFAC,IVAR,1)
                COEFA(IFAC,ICLVAR) = HEXT*PIMP/(HINT +HEXT)
                COEFB(IFAC,ICLVAR) = HINT     /(HINT +HEXT)
              ENDIF
C
C      C.L DE TYPE FLUX
            ELSEIF( ICODCL(IFAC,IVAR).EQ.3 ) THEN
              COEFA(IFAC,ICLVAR) = -RCODCL(IFAC,IVAR,3)/HINT
              COEFB(IFAC,ICLVAR) = 1.D0
            ENDIF
C
          ENDDO
C
          IF (ICODCL(IFAC,IUMA).EQ.4) THEN
C     Face de glissement (si IUMA est de type 4, les autres aussi)
C     On force la vitesse de maillage normale a etre nulle, CL de Neumann
C       sur les autres composantes (calque sur CLSYVT). On prend directement
C       la valeur au centre de la cellule et pas reconstruite a la face (on
C       ne cherche pas une precision particuliere sur w)
            SRFBNF = RA(ISRFBN-1+IFAC)
            RNX = SURFBO(1,IFAC)/SRFBNF
            RNY = SURFBO(2,IFAC)/SRFBNF
            RNZ = SURFBO(3,IFAC)/SRFBNF
            UPX = RTPA(IEL,IUMA)
            UPY = RTPA(IEL,IVMA)
            UPZ = RTPA(IEL,IWMA)
            COEFA(IFAC,IUMA) = - RNX*(RNY*UPY+RNZ*UPZ)
            COEFB(IFAC,IUMA) = 1.D0-RNX**2
            COEFA(IFAC,IVMA) = - RNY*(RNZ*UPZ+RNX*UPX)
            COEFB(IFAC,IVMA) = 1.D0-RNY**2
            COEFA(IFAC,IWMA) = - RNZ*(RNX*UPX+RNY*UPY)
            COEFB(IFAC,IWMA) = 1.D0-RNZ**2
          ENDIF
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C 14.  CALCUL DES EFFORTS AUX BORDS (partie 1/5)
C=======================================================================
C
      IF (INEEDF.EQ.1) THEN
        IPHAS = 1
        DO IFAC = 1, NFABOR
          IEL = IFABOR(IFAC)
          VISCLC = PROPCE(IEL,IPPROC(IVISCL(IPHAS)))
          VISCTC = PROPCE(IEL,IPPROC(IVISCT(IPHAS)))
          IF (ITYTUR(IPHAS).EQ.3) THEN
            VISTOT = VISCLC
          ELSE
            VISTOT = VISCLC + VISCTC
          ENDIF
          DISTBF = RA(IDISTB-1+IFAC)
          SRFBNF = RA(ISRFBN-1+IFAC)
          RA(IFORBR+(IFAC-1)*NDIM)     = -VISTOT * ( COEFA(IFAC,ICLUF)
     &         + (COEFB(IFAC,ICLUF)-1.D0)*COEFU(IFAC,1) )/DISTBF*SRFBNF
          RA(IFORBR+(IFAC-1)*NDIM + 1) = -VISTOT * ( COEFA(IFAC,ICLVF)
     &         + (COEFB(IFAC,ICLVF)-1.D0)*COEFU(IFAC,2) )/DISTBF*SRFBNF
          RA(IFORBR+(IFAC-1)*NDIM + 2) = -VISTOT * ( COEFA(IFAC,ICLWF)
     &         + (COEFB(IFAC,ICLWF)-1.D0)*COEFU(IFAC,3) )/DISTBF*SRFBNF
        ENDDO
      ENDIF
C
C=======================================================================
C 15.  FORMATS
C=======================================================================
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,
     &'@    *********                                               ',/,
     &'@      INCOHERENCE ENTRE OPTIONS DE CALCUL ET COND. LIM.     ',/,
     &'@                                                            ',/,
     &'@    Pour la phase ',I10   ,'                                ',/,
     &'@      la prise en compte des termes d''echo de paroi        ',/,
     &'@      du modele de turbulence Rij-epsilon est activee       ',/,
     &'@      IRIJEC(',I10   ,') = ',I10,'                          ',/,
     &'@      Ou bien l amortissement de la viscosite turbulente    ',/,
     &'@      est active IDRIES(',I10   ,') = ',I10,'en LES         ',/,
     &'@    mais aucune face de bord de type paroi n''est detectee. ',/,
     &'@    L''incoherence indiquee ci-dessus n''est pas bloquante  ',/,
     &'@      mais peut resulter d''une erreur lors de la           ',/,
     &'@      specification des conditions aux limites.             ',/,
     &'@                                                            ',/,
     &'@    Par securite, le calcul ne sera pas execute.            ',/,
     &'@                                                            ',/,
     &'@    Verifier les conditions aux limites dans usclim si le   ',/,
     &'@      domaine comporte des parois.                          ',/,
     &'@    Eliminer l''option IRIJEC de usini1 si le domaine ne    ',/,
     &'@      comporte pas de paroi (ou conditions ICODCL = 5 en    ',/,
     &'@      vitesse).                                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3010 FORMAT(
     & 'Phase ',I4,' debit entrant retenu en ',I10   ,
     &                                      ' faces de sortie sur ',I10)
 8000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : COND. LIM.                                  ',/,
     &'@    *********                                               ',/,
     &'@     Le scalaire ',I10   ,' est couple a SYRTHES            ',/,
     &'@     Il est relatif a la phase ',I10                         ,/,
     &'@      mais n''est pas la variable energetique               ',/,
     &'@         ISCALT(IPHAS) = ',I10                               ,/,
     &'@                                                            ',/,
     &'@     Le calcul ne sera pas execute.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
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
