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
                        SUBROUTINE MTTYCL
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , ITRIFB , ICODCL , ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RCODCL , FRCXT  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC TRAITEMENT DES CODES DE CONDITIONS AUX LIMITES DE MATISSE
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
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS)      !    !     !                                      !
CARGU ! ITRIFB(NFABOR! TE ! <-  ! TAB D'INDIRECTION POUR TRI DES FACES !
CARGU !  NPHAS)      !    !     !                                      !
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
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
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
CARGU ! FRCXT(NCELET,! TR !  -> ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU !   3,NPHAS)   !    !     !  HYDROSTATIQUE                       !
CARGU ! W1,2,3,4,5,6 ! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU !  (NCELET     !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
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
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "pointe.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
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
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICODCL(NFABOR,NVAR)
      INTEGER          ITYPFB(NFABOR,NPHAS) , ITRIFB(NFABOR,NPHAS)
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
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION RCODCL(NFABOR,NVAR,3)
      DOUBLE PRECISION FRCXT(NCELET,3,NPHAS)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION W4(NCELET),W5(NCELET),W6(NCELET)
      DOUBLE PRECISION COEFU(NFABOR,3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IPHAS , IFAC  , IFML  , ICOUL
      DOUBLE PRECISION DBM   , ROE
C
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
C --- Gestion memoire
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C --- Une seule phase
C
      IPHAS = 1
C
C=======================================================================
C 2.  BOUCLE SUR LES FACES DE BORD
C=======================================================================
C
      DO IFAC = 1, NFABOR
C
C --- Couleur de la face de bord
C
        IFML  = IFMFBR(IFAC  )
        ICOUL = IPRFML(IFML,1)
C
C
C=======================================================================
C 3.  ENTREE : classique en convection forcee, a pression imposee sinon
C=======================================================================
C
        IF(ICOUL.EQ.ICMTFI) THEN
C
C
C --- Convection forcee
          IF(ICOFOR.EQ.1)THEN
C
C     Type de base : IENTRE
            ITYPFB(IFAC,IPHAS) = IENTRE
C
C     Temperature exterieure prise a TINIT
            RCODCL(IFAC,ISCA(ITAAMT),1) = TINIT
            RCODCL(IFAC,ISCA(ITPCMT),1) = TINIT
            RCODCL(IFAC,ISCA(ITPPMT),1) = TINIT
C
C     Vitesses transverses nulles
            RCODCL(IFAC,IU(IPHAS),1) = 0.D0
            RCODCL(IFAC,IV(IPHAS),1) = 0.D0
C
C     Vitesse debitante W est calculee a partir du debit massique par
C       W = - debit massique / (masse volumique * surface)
C     . On suppose la vitesse verticale descendante (d'ou le signe).
C     . ROE est la masse volumique exterieure evaluee a partir de
C       la temperature TINIT (degres C) en utilisant la loi d'etat.
C     . La surface de l'entree du maillage est calculee par la formule
C       NPTRAN*PTRRES*EPCHEM/RCONVE. En effet, la surface au bas de la
C       cheminee est NPTRAN*PTRRES*EPCHEM et on divise par le rapport
C       du convergent represente sur le maillage (>=1) pour obtenir
C       la surface d'entree. (en 2D, il faut RCONVE=1, sinon, on impose
C       un debit plus grand que le debit souhaite).
C     . La correction du debit massique par FRDTRA correspond au
C       rapport d'echelle transverse eventuel entre le maillage et la
C       realite.
C
            ROE = RRFMAT*(TRFMAT+TKELVI)/(TINIT+TKELVI)
            DBM = DEBMAS/FRDTRA
            RCODCL(IFAC,IW(IPHAS),1) =
     &             -DBM*RCONVE/(ROE*NPTRAN*PTRRES*EPCHEM)
C
C
C --- Convection naturelle (ICOFOR=0)
          ELSE
C
C     Type de base indefini
C       (dans une version ulterieure de Code_Saturne, on disposera de
C        sorties a pression imposee)
            ITYPFB(IFAC,IPHAS) = IINDEF
C
C     Dirichlet sur les temperatures
            ICODCL(IFAC,ISCA(ITAAMT)  ) = 1
            RCODCL(IFAC,ISCA(ITAAMT),1) = TINIT
            ICODCL(IFAC,ISCA(ITPCMT)  ) = 1
            RCODCL(IFAC,ISCA(ITPCMT),1) = TINIT
            ICODCL(IFAC,ISCA(ITPPMT)  ) = 1
            RCODCL(IFAC,ISCA(ITPPMT),1) = TINIT
C
C     Dirichlet nul sur les vitesses transverses (nul par symetrie)
            ICODCL(IFAC,IU(IPHAS)  ) = 1
            RCODCL(IFAC,IU(IPHAS),1) = 0.D0
            ICODCL(IFAC,IV(IPHAS)  ) = 1
            RCODCL(IFAC,IV(IPHAS),1) = 0.D0
C
C     Neumann homogene sur la vitesse debitante (on impose la pression)
            ICODCL(IFAC,IW(IPHAS)  ) = 3
C
C     Pression nulle
            ICODCL(IFAC,IPR(IPHAS)  ) = 1
            RCODCL(IFAC,IPR(IPHAS),1) = P0(IPHAS)
C
          ENDIF
C
C
C=======================================================================
C 4.  SORTIE (pression imposee)
C=======================================================================
C
        ELSEIF(ICOUL.EQ.ICMTFO) THEN
C
C     Type de base : sortie libre
          ITYPFB(IFAC,IPHAS) = ISOLIB
C
C     Temperature imposee en cas de reentree
          RCODCL(IFAC,ISCA(ITAAMT),1) = TCRIT
          RCODCL(IFAC,ISCA(ITPCMT),1) = TINIT
          RCODCL(IFAC,ISCA(ITPPMT),1) = TINIT
C
C     Pression imposee (ecart hydrostatique par rapport a l'entree
C                       + decalage eventuel DPVENT)
          ICODCL(IFAC,IPR(IPHAS)  ) = 1
          RCODCL(IFAC,IPR(IPHAS),1) = P0(IPHAS) + DPVENT
     &         - DABS(GZ) * RO0(IPHAS) * (HCHEVA - HCHALI)
C
C=======================================================================
C 5.  LE RESTE : glissement et Neumann homogene
C=======================================================================
C
C --- Sol
        ELSEIF(ICOUL.EQ.ICMTFG) THEN
          ITYPFB(IFAC,IPHAS)   = ISYMET
C
C --- Plafond
        ELSEIF(ICOUL.EQ.ICMTFC) THEN
          ITYPFB(IFAC,IPHAS)   = ISYMET
C
C --- Symetries
        ELSEIF(ICOUL.EQ.ICMTFS) THEN
          ITYPFB(IFAC,IPHAS)   = ISYMET
C
C --- Parois
        ELSEIF(ICOUL.EQ.ICMTFW) THEN
          ITYPFB(IFAC,IPHAS)   = ISYMET
C
        ENDIF
C
      ENDDO
C
C
C----
C FORMATS
C----
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
