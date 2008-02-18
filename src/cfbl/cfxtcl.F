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
                        SUBROUTINE CFXTCL
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , ITRIFB , ITYPFB , IZFPPP ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RCODCL ,
     &   W1     , W2     , W3     , W4     , W5     , W6     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC     CONDITIONS AUX LIMITES AUTOMATIQUES
CFONC
CFONC            COMPRESSIBLE SANS CHOC
CFONC
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
CARGU ! ICODCL       ! TE ! <-  ! CODE DE CONDITION LIMITES AUX FACES  !
CARGU !  (NFABOR,NVAR!    !     !  DE BORD                             !
CARGU !              !    !     ! = 1   -> DIRICHLET                   !
CARGU !              !    !     ! = 3   -> DENSITE DE FLUX             !
CARGU !              !    !     ! = 4   -> GLISSEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 5   -> FROTTEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 6   -> RUGOSITE ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 9   -> ENTREE/SORTIE LIBRE (VITESSE!
CARGU ! ITRIFB(NFABOR! TE !  -> ! INDIRECTION POUR TRI DES FACES DE BRD!
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! IZFPPP       ! TE !  -> ! NUMERO DE ZONE DE LA FACE DE BORD    !
CARGU ! (NFABOR)     !    !     !  POUR LE MODULE PHYS. PART.          !
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
CARGU ! W1,2,3,4,5,6 ! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU !  (NCELET     !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
CARGU ! COEFU        ! TR !  -  ! TAB DE TRAV                          !
CARGU !  (NFABOR,3)  !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
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
C ARGUMENTS
C
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "cfpoin.h"
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
      INTEGER          ITRIFB(NFABOR,NPHAS), ITYPFB(NFABOR,NPHAS)
      INTEGER          IZFPPP(NFABOR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION RCODCL(NFABOR,NVAR,3)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION W4(NCELET),W5(NCELET),W6(NCELET)
      DOUBLE PRECISION COEFU(NFABOR,NDIM)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IPHAS , IVAR  , IFAC  , IEL
      INTEGER          II    , III   , IMODIF, ICCFTH
      INTEGER          ICALEP, ICALGM
      INTEGER          IFLMAB
      INTEGER          IPRIPH, IUIPH , IVIPH , IWIPH
      INTEGER          IRHIPH, IENIPH, ITKIPH
      INTEGER          ICLP  , ICLR
      INTEGER          ICLU  , ICLV  , ICLW
      INTEGER          NVARCF
C
      INTEGER          NVCFMX
      PARAMETER       (NVCFMX=7)
      INTEGER          IVARCF(NVCFMX)
C
      DOUBLE PRECISION HINT  , GAMMAG
C
C***********************************************************************
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
      DO IPHAS = 1, NPHAS
C
        IPRIPH = IPR   (IPHAS)
        IUIPH  = IU    (IPHAS)
        IVIPH  = IV    (IPHAS)
        IWIPH  = IW    (IPHAS)
        IRHIPH = ISCA(IRHO  (IPHAS))
        IENIPH = ISCA(IENERG(IPHAS))
        ITKIPH = ISCA(ITEMPK(IPHAS))
        ICLP   = ICLRTP(IPRIPH,ICOEF)
        ICLR   = ICLRTP(IRHIPH,ICOEF)
        ICLU   = ICLRTP(IUIPH ,ICOEF)
        ICLV   = ICLRTP(IVIPH ,ICOEF)
        ICLW   = ICLRTP(IWIPH ,ICOEF)
C
        IFLMAB = IPPROB(IFLUMA(IENIPH))
C
C     Liste des variables compressible :
        IVARCF(1) = IPRIPH
        IVARCF(2) = IUIPH
        IVARCF(3) = IVIPH
        IVARCF(4) = IWIPH
        IVARCF(5) = IRHIPH
        IVARCF(6) = IENIPH
        IVARCF(7) = ITKIPH
        NVARCF    = 7
C
C     Calcul de epsilon_sup = e - CvT
C       On en a besoin si on a des parois a temperature imposee.
C       Il est calculé aux cellules W5 et aux faces de bord COEFU.
C       On n'en a besoin ici qu'aux cellules de bord : s'il est
C         nécessaire de gagner de la mémoire, on pourra modifier
C         uscfth.
C
        ICALEP = 0
        DO IFAC = 1, NFABOR
          IF(ICODCL(IFAC,ITKIPH).EQ.5) THEN
            ICALEP = 1
          ENDIF
        ENDDO
        IF(ICALEP.NE.0) THEN
          ICCFTH = 7
          IMODIF = 0
          CALL USCFTH
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IMODIF , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W5     , COEFU(1,1) , W3 , W4     ,
C        ------   ---------
     &   RDEVEL , RTUSER , RA     )
        ENDIF
C
C
C     Calcul de gamma (constant ou variable ; pour le moment : cst)
C       On en a besoin pour les entrees sorties avec rusanov
C
        ICALGM = 0
        DO IFAC = 1, NFABOR
          IF ( ( ITYPFB(IFAC,IPHAS).EQ.IESICF ) .OR.
     &         ( ITYPFB(IFAC,IPHAS).EQ.ISOPCF ) .OR.
     &         ( ITYPFB(IFAC,IPHAS).EQ.IERUCF ) .OR.
     &         ( ITYPFB(IFAC,IPHAS).EQ.IEQHCF ) ) THEN
            ICALGM = 1
          ENDIF
        ENDDO
        IF(ICALGM.NE.0) THEN
          ICCFTH = 1
          IMODIF = 0
          CALL USCFTH
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IMODIF , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W6     , W4     ,
     &   RDEVEL , RTUSER , RA     )
C
          IF(IEOS(IPHAS).EQ.1) THEN
            GAMMAG = W6(1)
          ELSE
C     Gamma doit etre passe a cfrusb ; s'il est variable
C       il est dans le tableau W6 et il faut ajouter
C           GAMMAG = W6(IFABOR(IFAC)) selon IEOS
C       dans la boucle sur les faces.
C     En attendant que IEOS different de 1 soit code, on stoppe
            WRITE(NFECRA,7000)
            CALL CSEXIT (1)
          ENDIF
C
        ENDIF
C
C
C
C     Boucle sur les faces
C
        DO IFAC = 1, NFABOR
          IEL = IFABOR(IFAC)
C
C=======================================================================
C 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
C       ON BOUCLE SUR TOUTES LES FACES DE PAROI
C=======================================================================
C
          IF ( ITYPFB(IFAC,IPHAS).EQ.IPAROI) THEN
C
C     Les RCODCL ont ete initialises a -RINFIN pour permettre de
C       verifier ceux que l'utilisateur a modifies. On les remet a zero
C       si l'utilisateur ne les a pas modifies.
C       En paroi, on traite toutes les variables.
            DO IVAR = 1, NVAR
              IF(RCODCL(IFAC,IVAR,1).LE.-RINFIN*0.5D0) THEN
                RCODCL(IFAC,IVAR,1) = 0.D0
              ENDIF
            ENDDO
C
C     Le flux de masse est nul
C
            PROPFB(IFAC,IFLMAB) = 0.D0
C
C     Pression :
C
C       Si la gravite est predominante : pression hydrostatique
C         (approximatif et surtout explicite en rho)
C
            IF(ICFGRP(IPHAS).EQ.1) THEN
C
              ICODCL(IFAC,IPRIPH) = 3
              HINT = DT(IEL)/RA(IDISTB-1+IFAC)
              RCODCL(IFAC,IPRIPH,3) = -HINT
     &         * ( GX*(CDGFBO(1,IFAC)-XYZCEN(1,IEL))
     &           + GY*(CDGFBO(2,IFAC)-XYZCEN(2,IEL))
     &           + GZ*(CDGFBO(3,IFAC)-XYZCEN(3,IEL)) )
     &         * RTP(IEL,IRHIPH)
C
            ELSE

C       En général : proportionnelle a la valeur interne
C         (Pbord = COEFB*Pi)
C       Si on détend trop : Dirichlet homogene
C
              ICCFTH = 91
C
              CALL USCFTH
C             ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IFAC   , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   RDEVEL , RTUSER , RA     )
C
C       En outre, il faut appliquer une pre-correction pour compenser
C        le traitement fait dans condli... Si on pouvait remplir COEFA
C        et COEFB directement, on gagnerait en simplicite, mais cela
C        demanderait un test sur IPPMOD dans condli : à voir)
C
              ICODCL(IFAC,IPRIPH) = 1
              IF(COEFB(IFAC,ICLP).LT.RINFIN*0.5D0.AND.
     &           COEFB(IFAC,ICLP).GT.0.D0  ) THEN
                HINT = DT(IEL)/RA(IDISTB-1+IFAC)
                RCODCL(IFAC,IPRIPH,1) = 0.D0
                RCODCL(IFAC,IPRIPH,2) =
     &               HINT*(1.D0/COEFB(IFAC,ICLP)-1.D0)
              ELSE
                RCODCL(IFAC,IPRIPH,1) = 0.D0
              ENDIF
C
            ENDIF
C
C
C       La vitesse et la turbulence sont traitées de manière standard,
C         dans condli.
C
C       Pour la thermique, on doit effectuer ici un prétraitement,
C         la variable résolue étant l'energie
C         (energie interne+epsilon sup+energie cinétique). En particulier
C         lorsque la paroi est à température imposée, on prépare le
C         travail de clptur. Hormis l'énergie résolue, toutes les
C         variables rho et s prendront arbitrairement une condition de
C         flux nul (leurs conditions aux limites ne servent qu'à la
C         reconstruction des gradients et il parait délicat d'imposer
C         autre chose qu'un flux nul sans risque de créer des valeurs
C         aberrantes au voisinage de la couche limite)
C
C       Par défaut : adiabatique
            IF(  ICODCL(IFAC,ITKIPH).EQ.0.AND.
     &           ICODCL(IFAC,IENIPH).EQ.0) THEN
              ICODCL(IFAC,ITKIPH) = 3
              RCODCL(IFAC,ITKIPH,3) = 0.D0
            ENDIF
C
C       Temperature imposee
            IF(ICODCL(IFAC,ITKIPH).EQ.5) THEN
C
C           On impose la valeur de l'energie qui conduit au bon flux.
C             On notera cependant qu'il s'agit de la condition à la
C               limite pour le flux diffusif. Pour la reconstruction
C               des gradients, il faudra utiliser autre chose.
C               Par exemple un flux nul ou encore toute autre
C               condition respectant un profil : on se calquera sur
C               ce qui sera fait pour la température si c'est possible,
C               sachant que l'energie contient l'energie cinetique,
C               ce qui rend le choix du profil délicat.
C
              ICODCL(IFAC,IENIPH) = 5
              IF(ICV(IPHAS).EQ.0) THEN
                RCODCL(IFAC,IENIPH,1) =
     &               CV0(IPHAS)*RCODCL(IFAC,ITKIPH,1)
              ELSE
                RCODCL(IFAC,IENIPH,1) = PROPCE(IEL,IPPROC(ICV(IPHAS)))
     &               *RCODCL(IFAC,ITKIPH,1)
              ENDIF
              RCODCL(IFAC,IENIPH,1) = RCODCL(IFAC,IENIPH,1)
     &           + 0.5D0*(RTP(IEL,IUIPH)**2+
     &                    RTP(IEL,IVIPH)**2+RTP(IEL,IWIPH)**2)
     &           + W5(IEL)
C                   ^epsilon sup (cf USCFTH)
C
C           Les flux en grad epsilon sup et énergie cinétique doivent
C             être nuls puisque tout est pris par le terme de
C             diffusion d'energie.
              IA(IIFBET+IFAC-1+(IPHAS-1)*NFABOR) = 1
C
C           Flux nul pour la reconstruction éventuelle de température
              ICODCL(IFAC,ITKIPH) = 3
              RCODCL(IFAC,ITKIPH,3) = 0.D0
C
C       Flux impose
            ELSEIF(ICODCL(IFAC,ITKIPH).EQ.3) THEN
C
C           On impose le flux sur l'energie
              ICODCL(IFAC,IENIPH) = 3
              RCODCL(IFAC,IENIPH,3) = RCODCL(IFAC,ITKIPH,3)
C
C           Les flux en grad epsilon sup et énergie cinétique doivent
C             être nuls puisque tout est pris par le terme de
C             diffusion d'energie.
              IA(IIFBET+IFAC-1+(IPHAS-1)*NFABOR) = 1
C
C           Flux nul pour la reconstruction éventuelle de température
              ICODCL(IFAC,ITKIPH) = 3
              RCODCL(IFAC,ITKIPH,3) = 0.D0
C
            ENDIF
C
C
C     Scalaires : flux nul (par defaut dans typecl pour iparoi)
C
C
C=======================================================================
C 3.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
C       ON BOUCLE SUR TOUTES LES FACES DE SYMETRIE
C=======================================================================
C
          ELSEIF ( ITYPFB(IFAC,IPHAS).EQ.ISYMET ) THEN
C
C     Les RCODCL ont ete initialises a -RINFIN pour permettre de
C       verifier ceux que l'utilisateur a modifies. On les remet a zero
C       si l'utilisateur ne les a pas modifies.
C       En symetrie, on traite toutes les variables.
            DO IVAR = 1, NVAR
              IF(RCODCL(IFAC,IVAR,1).LE.-RINFIN*0.5D0) THEN
                RCODCL(IFAC,IVAR,1) = 0.D0
              ENDIF
            ENDDO
C
C     Le flux de masse est nul
C
            PROPFB(IFAC,IFLMAB) = 0.D0
C
C     Condition de Pression
C
            ICCFTH = 90
C
            CALL USCFTH
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IFAC   , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C     Pression :
C       En général : proportionnelle a la valeur interne
C         (Pbord = COEFB*Pi)
C       Si on détend trop : Dirichlet homogene
C
C       En outre, il faut appliquer une pre-correction pour compenser le
C        traitement fait dans condli... Si on pouvait remplir COEFA
C        et COEFB directement, on gagnerait en simplicite, mais cela
C        demanderait un test sur IPPMOD dans condli : à voir)
C
            ICODCL(IFAC,IPRIPH) = 3
            RCODCL(IFAC,IPRIPH,1) = 0.D0
            RCODCL(IFAC,IPRIPH,2) = RINFIN
            RCODCL(IFAC,IPRIPH,3) = 0.D0
C
C       Toutes les autres variables prennent un flux nul (sauf la vitesse
C         normale, qui est nulle) : par defaut dans typecl pour isymet.
C
C=======================================================================
C 4.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
C       ON BOUCLE SUR TOUTES LES FACES D'ENTREE/SORTIE
C       ETAPE DE THERMO
C=======================================================================
C
C
C=======================================================================
C     4.1 Entree/sortie imposée (par exemple : entree supersonique)
C=======================================================================
C
          ELSEIF ( ITYPFB(IFAC,IPHAS).EQ.IESICF ) THEN
C
C     On a
C       - la vitesse,
C       - 2 variables parmi P, rho, T, E (mais pas (T,E)),
C       - la turbulence
C       - les scalaires
C
C     On recherche la variable a initialiser
C       (si on a donne une valeur nulle, c'est pas adapte : on supposera
C        qu'on n'a pas initialise et on sort en erreur)
            ICCFTH = 10000
            IF(RCODCL(IFAC,IPRIPH,1).GT.0.D0) ICCFTH = 2*ICCFTH
            IF(RCODCL(IFAC,IRHIPH,1).GT.0.D0) ICCFTH = 3*ICCFTH
            IF(RCODCL(IFAC,ITKIPH,1).GT.0.D0) ICCFTH = 5*ICCFTH
            IF(RCODCL(IFAC,IENIPH,1).GT.0.D0) ICCFTH = 7*ICCFTH
            IF((ICCFTH.LE.70000.AND.ICCFTH.NE.60000).OR.
     &         (ICCFTH.EQ.350000)) THEN
              WRITE(NFECRA,1000)ICCFTH
              CALL CSEXIT (1)
            ENDIF
            ICCFTH = ICCFTH + 900
C
C     Les RCODCL ont ete initialises a -RINFIN pour permettre de
C       verifier ceux que l'utilisateur a modifies. On les remet a zero
C       si l'utilisateur ne les a pas modifies.
C       On traite d'abord les variables autres que la turbulence et les
C       scalaires passifs : celles-ci sont traitees plus bas.
            DO III = 1, NVARCF
              IVAR = IVARCF(III)
              IF(RCODCL(IFAC,IVAR,1).LE.-RINFIN*0.5D0) THEN
                RCODCL(IFAC,IVAR,1) = 0.D0
              ENDIF
            ENDDO
C
C     On calcule les variables manquantes parmi P,rho,T,E
C     COEFA sert de tableau de transfert dans USCFTH
C
            DO IVAR = 1, NVAR
              COEFA(IFAC,ICLRTP(IVAR,ICOEF)) = RCODCL(IFAC,IVAR,1)
            ENDDO
C
            CALL USCFTH
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IFAC   , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C     Rusanov, flux de masse et type de conditions aux limites :
C       voir plus bas
C
C
C=======================================================================
C     4.2 Sortie supersonique
C=======================================================================
C
          ELSEIF ( ITYPFB(IFAC,IPHAS).EQ.ISSPCF ) THEN
C
C     On impose un Dirichlet égal à la valeur interne pour rho u E
C       (on impose des Dirichlet déduit pour les autres variables).
C       Il est inutile de passer dans Rusanov.
C     Il serait nécessaire de reconstruire ces valeurs en utilisant
C       leur gradient dans la cellule de bord : dans un premier temps,
C       on utilise des valeurs non reconstruites (non consistant mais
C       potentiellement plus stable).
C     On pourrait imposer des flux nuls (a tester), ce qui éviterait
C       la nécessité de reconstruire les valeurs.
C
C     Les RCODCL ont ete initialises a -RINFIN pour permettre de
C       verifier ceux que l'utilisateur a modifies. On les remet a zero
C       si l'utilisateur ne les a pas modifies.
C       On traite d'abord les variables autres que la turbulence et les
C       scalaires passifs : celles-ci sont traitees plus bas.
            DO III = 1, NVARCF
              IVAR = IVARCF(III)
              IF(RCODCL(IFAC,IVAR,1).LE.-RINFIN*0.5D0) THEN
                RCODCL(IFAC,IVAR,1) = 0.D0
              ENDIF
            ENDDO
C
C     Valeurs de rho u E
            RCODCL(IFAC,IRHIPH,1) = RTP(IEL,IRHIPH)
            RCODCL(IFAC,IUIPH ,1) = RTP(IEL,IUIPH)
            RCODCL(IFAC,IVIPH ,1) = RTP(IEL,IVIPH)
            RCODCL(IFAC,IWIPH ,1) = RTP(IEL,IWIPH)
            RCODCL(IFAC,IENIPH,1) = RTP(IEL,IENIPH)
C
C     Valeurs de P et s déduites
            ICCFTH = 924
C
            DO IVAR = 1, NVAR
              COEFA(IFAC,ICLRTP(IVAR,ICOEF)) = RCODCL(IFAC,IVAR,1)
            ENDDO
C
            CALL USCFTH
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IFAC   , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   RDEVEL , RTUSER , RA     )
C
C               flux de masse et type de conditions aux limites :
C       voir plus bas
C
C
C=======================================================================
C     4.3 Sortie a pression imposee
C=======================================================================
C
          ELSEIF ( ITYPFB(IFAC,IPHAS).EQ.ISOPCF ) THEN
C
C       Sortie subsonique a priori (si c'est supersonique dans le
C         domaine, ce n'est pas pour autant que c'est supersonique
C         à la sortie, selon la pression que l'on a imposée)
C
C     On utilise un scenario dans lequel on a une 1-détente et un
C       2-contact entrant dans le domaine. On détermine les conditions
C       sur l'interface selon la thermo et on passe dans Rusanov
C       ensuite pour lisser.
C
C     Si P n'est pas donné, erreur ; on sort aussi en erreur si P
C       négatif, même si c'est possible, dans la plupart des cas ce
C       sera une erreur
            IF(RCODCL(IFAC,IPRIPH,1).LT.-RINFIN*0.5D0) THEN
              WRITE(NFECRA,1100)
              CALL CSEXIT (1)
            ENDIF
C
C     Les RCODCL ont ete initialises a -RINFIN pour permettre de
C       verifier ceux que l'utilisateur a modifies. On les remet a zero
C       si l'utilisateur ne les a pas modifies.
C       On traite d'abord les variables autres que la turbulence et les
C       scalaires passifs : celles-ci sont traitees plus bas.
            DO III = 1, NVARCF
              IVAR = IVARCF(III)
              IF(RCODCL(IFAC,IVAR,1).LE.-RINFIN*0.5D0) THEN
                RCODCL(IFAC,IVAR,1) = 0.D0
              ENDIF
            ENDDO
C
C     Valeurs de rho, u, E, s
            ICCFTH = 93
C
            DO IVAR = 1, NVAR
              COEFA(IFAC,ICLRTP(IVAR,ICOEF)) = RCODCL(IFAC,IVAR,1)
            ENDDO
C
            CALL USCFTH
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IFAC   , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   RDEVEL , RTUSER , RA     )
C
C     Rusanov, flux de masse et type de conditions aux limites :
C       voir plus bas
C
C
C=======================================================================
C     4.4 Entree à rho et U imposes
C=======================================================================
C
          ELSEIF ( ITYPFB(IFAC,IPHAS).EQ.IERUCF ) THEN
C
C       Entree subsonique a priori (si c'est supersonique dans le
C         domaine, ce n'est pas pour autant que c'est supersonique
C         à l'entree, selon les valeurs que l'on a imposées)
C
C     On utilise un scenario détente ou choc.
C       On détermine les conditions sur l'interface
C       selon la thermo et on passe dans Rusanov ensuite pour lisser.
C
C     Si rho et u ne sont pas donnés, erreur
            IF(RCODCL(IFAC,IRHIPH,1).LT.-RINFIN*0.5D0.OR.
     &         RCODCL(IFAC,IUIPH ,1).LT.-RINFIN*0.5D0.OR.
     &         RCODCL(IFAC,IVIPH ,1).LT.-RINFIN*0.5D0.OR.
     &         RCODCL(IFAC,IWIPH ,1).LT.-RINFIN*0.5D0) THEN
              WRITE(NFECRA,1200)
              CALL CSEXIT (1)
            ENDIF
C
C     Les RCODCL ont ete initialises a -RINFIN pour permettre de
C       verifier ceux que l'utilisateur a modifies. On les remet a zero
C       si l'utilisateur ne les a pas modifies.
C       On traite d'abord les variables autres que la turbulence et les
C       scalaires passifs : celles-ci sont traitees plus bas.
            DO III = 1, NVARCF
              IVAR = IVARCF(III)
              IF(RCODCL(IFAC,IVAR,1).LE.-RINFIN*0.5D0) THEN
                RCODCL(IFAC,IVAR,1) = 0.D0
              ENDIF
            ENDDO
C
C     Valeurs de P, E, s
            ICCFTH = 92
C
            DO IVAR = 1, NVAR
              COEFA(IFAC,ICLRTP(IVAR,ICOEF)) = RCODCL(IFAC,IVAR,1)
            ENDDO
C
            CALL USCFTH
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IFAC   , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   RDEVEL , RTUSER , RA     )
C
C     Rusanov, flux de masse et type de conditions aux limites :
C       voir plus bas
C
C
C=======================================================================
C     4.5 Entree à rho*U et rho*U*H imposes
C=======================================================================
C
          ELSEIF ( ITYPFB(IFAC,IPHAS).EQ.IEQHCF ) THEN
C
C       Entree subsonique a priori (si c'est supersonique dans le
C         domaine, ce n'est pas pour autant que c'est supersonique
C         à l'entree, selon les valeurs que l'on a imposées)
C
C     On utilise un scenario dans lequel on a un 2-contact et une
C       3-détente entrant dans le domaine. On détermine les conditions
C       sur l'interface selon la thermo et on passe dans Rusanov
C       ensuite pour lisser.
C
C     Si rho et u ne sont pas donnés, erreur
            IF(RCODCL(IFAC,IRUN (IPHAS),1).LT.-RINFIN*0.5D0.OR.
     &         RCODCL(IFAC,IRUNH(IPHAS),1).LT.-RINFIN*0.5D0) THEN
              WRITE(NFECRA,1300)
              CALL CSEXIT (1)
            ENDIF
C
C     Les RCODCL ont ete initialises a -RINFIN pour permettre de
C       verifier ceux que l'utilisateur a modifies. On les remet a zero
C       si l'utilisateur ne les a pas modifies.
C       On traite d'abord les variables autres que la turbulence et les
C       scalaires passifs : celles-ci sont traitees plus bas.
            DO III = 1, NVARCF
              IVAR = IVARCF(III)
              IF(RCODCL(IFAC,IVAR,1).LE.-RINFIN*0.5D0) THEN
                RCODCL(IFAC,IVAR,1) = 0.D0
              ENDIF
            ENDDO
C
C     A coder
C
C     Noter que IRUN(IPHAS)  = ISCA(IRHO (IPHAS))
C            et IRUNH(IPHAS) = ISCA(IENER(IPHAS))
C     (aliases pour simplifier uscfcl)
C
            WRITE(NFECRA,1301)
            CALL CSEXIT (1)
C
C=======================================================================
C 5. CONDITION NON PREVUE
C=======================================================================
C     Stop
          ELSE
C
            WRITE(NFECRA,1400)
            CALL CSEXIT (1)
C
C --- Fin de test sur les types de faces
          ENDIF
C
C
C=======================================================================
C 6. FIN DU TRAITEMENT DES ENTREE/SORTIES
C     CALCUL DU FLUX DE MASSE,
C     CALCUL DES FLUX DE BORD AVEC RUSANOV (SI BESOIN)
C     TYPE DE C.L. (DIRICHLET NEUMANN)
C=======================================================================
C
          IF ( ( ITYPFB(IFAC,IPHAS).EQ.IESICF ) .OR.
     &         ( ITYPFB(IFAC,IPHAS).EQ.ISSPCF ) .OR.
     &         ( ITYPFB(IFAC,IPHAS).EQ.ISOPCF ) .OR.
     &         ( ITYPFB(IFAC,IPHAS).EQ.IERUCF ) .OR.
     &         ( ITYPFB(IFAC,IPHAS).EQ.IEQHCF ) ) THEN
C
C=======================================================================
C     6.1 Flux de bord Rusanov ou simplement flux de masse
C         Attention a bien avoir calcule gamma pour Rusanov
C=======================================================================
C
C     Sortie supersonique :
            IF ( ITYPFB(IFAC,IPHAS).EQ.ISSPCF ) THEN
C
C     Seul le flux de masse est calcule (on n'appelle pas Rusanov)
C       (toutes les variables sont connues)
C
              PROPFB(IFAC,IFLMAB) = COEFA(IFAC,ICLR)*
     &             ( COEFA(IFAC,ICLU)*SURFBO(1,IFAC)
     &             + COEFA(IFAC,ICLV)*SURFBO(2,IFAC)
     &             + COEFA(IFAC,ICLW)*SURFBO(3,IFAC) )
C
C     Entree subsonique
C
            ELSE IF ( ITYPFB(IFAC,IPHAS).EQ.IERUCF ) THEN
C
C     Seul le flux de masse est calcule (on n'appelle pas Rusanov)
C
              PROPFB(IFAC,IFLMAB) = COEFA(IFAC,ICLR)*
     &             ( COEFA(IFAC,ICLU)*SURFBO(1,IFAC)
     &             + COEFA(IFAC,ICLV)*SURFBO(2,IFAC)
     &             + COEFA(IFAC,ICLW)*SURFBO(3,IFAC) )
C
C
C
C     Autres entrees/sorties :
            ELSE
C
C     On calcule des flux par Rusanov (PROPFB)
C       (en particulier, le flux de masse est complete)
C
              CALL CFRUSB
C             ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IFAC   , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   GAMMAG ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   RDEVEL , RTUSER , RA     )
C
            ENDIF
C
C=======================================================================
C     6.2 Recuperation de COEFA
C=======================================================================
C
C     On rétablit COEFA dans RCODCL
            DO IVAR = 1, NVAR
              RCODCL(IFAC,IVAR,1) = COEFA(IFAC,ICLRTP(IVAR,ICOEF))
            ENDDO
C
C=======================================================================
C     6.3 Types de C.L.
C=======================================================================
C
C     P               : Dirichlet sauf IESICF : Neumann (choix arbitraire)
C     rho, U, E, T    : Dirichlet
C     k, R, eps, scal : Dirichlet/Neumann selon flux de masse
C
C     Pour P, le Neumann est censé etre moins genant pour les
C       reconstructions de gradient si la valeur de P fournie par
C       l'utilisateur est tres differente de la valeur interne.
C       Le choix est cependant arbitraire.
C
C     On suppose que par defaut,
C            RCODCL(IFAC,X,1) = utilisateur ou calcule ci-dessus
C            RCODCL(IFAC,X,2) = RINFIN
C            RCODCL(IFAC,X,3) = 0.D0
C       et si ICODCL(IFAC,X) = 3, seul RCODCL(IFAC,X,3) est utilisé
C
C
C-----------------------------------------------------------------------
C     Pression : Dirichlet ou Neumann homogene
C-----------------------------------------------------------------------
C
C       Entree sortie imposee : Neumann
            IF ( ITYPFB(IFAC,IPHAS).EQ.IESICF ) THEN
              ICODCL(IFAC,IPRIPH)   = 3
C       Entree subsonique
            ELSE IF ( ITYPFB(IFAC,IPHAS).EQ.IERUCF ) THEN
              ICODCL(IFAC,IPRIPH)   = 3
              RCODCL(IFAC,IPRIPH,3) = 0.D0
C       Autres entrees/sorties : Dirichlet
            ELSE
              ICODCL(IFAC,IPRIPH)   = 1
            ENDIF
C
C-----------------------------------------------------------------------
C     rho U E T : Dirichlet
C-----------------------------------------------------------------------
C
C     Masse volumique
            ICODCL(IFAC,IRHIPH)   = 1
C     Vitesse
            ICODCL(IFAC,IUIPH)    = 1
            ICODCL(IFAC,IVIPH)    = 1
            ICODCL(IFAC,IWIPH)    = 1
C     Energie totale
            ICODCL(IFAC,IENIPH)   = 1
C     Temperature
            ICODCL(IFAC,ITKIPH)   = 1
C
C-----------------------------------------------------------------------
C     turbulence et scalaires passifs : Dirichlet/Neumann selon flux
C-----------------------------------------------------------------------
C
C       Dirichlet ou Neumann homogène
C       On choisit un Dirichlet si le flux de masse est entrant et
C       que l'utilisateur a donné une valeur dans RCODCL
C
            IF(PROPFB(IFAC,IFLMAB).GE.0.D0) THEN
              IF(ITYTUR(IPHAS).EQ.2) THEN
                ICODCL(IFAC,IK (IPHAS)) = 3
                ICODCL(IFAC,IEP(IPHAS)) = 3
              ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
                ICODCL(IFAC,IR11(IPHAS)) = 3
                ICODCL(IFAC,IR22(IPHAS)) = 3
                ICODCL(IFAC,IR33(IPHAS)) = 3
                ICODCL(IFAC,IR12(IPHAS)) = 3
                ICODCL(IFAC,IR13(IPHAS)) = 3
                ICODCL(IFAC,IR23(IPHAS)) = 3
                ICODCL(IFAC,IEP (IPHAS)) = 3
              ELSEIF(ITURB(IPHAS).EQ.50) THEN
                ICODCL(IFAC,IK  (IPHAS)) = 3
                ICODCL(IFAC,IEP (IPHAS)) = 3
                ICODCL(IFAC,IPHI(IPHAS)) = 3
                ICODCL(IFAC,IFB (IPHAS)) = 3
              ELSEIF(ITURB(IPHAS).EQ.60) THEN
                ICODCL(IFAC,IK  (IPHAS)) = 3
                ICODCL(IFAC,IOMG(IPHAS)) = 3
              ENDIF
              IF(NSCAUS.GT.0) THEN
                DO II = 1, NSCAUS
                  ICODCL(IFAC,ISCA(II)) = 3
                ENDDO
              ENDIF
            ELSE
              IF(ITYTUR(IPHAS).EQ.2) THEN
                IF(RCODCL(IFAC,IK (IPHAS),1).GT.0.D0.AND.
     &             RCODCL(IFAC,IEP(IPHAS),1).GT.0.D0) THEN
                  ICODCL(IFAC,IK (IPHAS)) = 1
                  ICODCL(IFAC,IEP(IPHAS)) = 1
                ELSE
                  ICODCL(IFAC,IK (IPHAS)) = 3
                  ICODCL(IFAC,IEP(IPHAS)) = 3
                ENDIF
              ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
                IF(RCODCL(IFAC,IR11(IPHAS),1).GT.0.D0.AND.
     &             RCODCL(IFAC,IR22(IPHAS),1).GT.0.D0.AND.
     &             RCODCL(IFAC,IR33(IPHAS),1).GT.0.D0.AND.
     &             RCODCL(IFAC,IR12(IPHAS),1).GT.-RINFIN*0.5D0.AND.
     &             RCODCL(IFAC,IR13(IPHAS),1).GT.-RINFIN*0.5D0.AND.
     &             RCODCL(IFAC,IR23(IPHAS),1).GT.-RINFIN*0.5D0.AND.
     &             RCODCL(IFAC,IEP (IPHAS),1).GT.0.D0) THEN
                  ICODCL(IFAC,IR11(IPHAS)) = 1
                  ICODCL(IFAC,IR22(IPHAS)) = 1
                  ICODCL(IFAC,IR33(IPHAS)) = 1
                  ICODCL(IFAC,IR12(IPHAS)) = 1
                  ICODCL(IFAC,IR13(IPHAS)) = 1
                  ICODCL(IFAC,IR23(IPHAS)) = 1
                  ICODCL(IFAC,IEP (IPHAS)) = 1
                ELSE
                  ICODCL(IFAC,IR11(IPHAS)) = 3
                  ICODCL(IFAC,IR22(IPHAS)) = 3
                  ICODCL(IFAC,IR33(IPHAS)) = 3
                  ICODCL(IFAC,IR12(IPHAS)) = 3
                  ICODCL(IFAC,IR13(IPHAS)) = 3
                  ICODCL(IFAC,IR23(IPHAS)) = 3
                  ICODCL(IFAC,IEP (IPHAS)) = 3
                ENDIF
              ELSEIF(ITURB(IPHAS).EQ.50) THEN
                IF(RCODCL(IFAC,IK  (IPHAS),1).GT.0.D0.AND.
     &             RCODCL(IFAC,IEP (IPHAS),1).GT.0.D0.AND.
     &             RCODCL(IFAC,IPHI(IPHAS),1).GT.0.D0.AND.
     &             RCODCL(IFAC,IFB (IPHAS),1).GT.-RINFIN*0.5D0 ) THEN
                  ICODCL(IFAC,IK  (IPHAS)) = 1
                  ICODCL(IFAC,IEP (IPHAS)) = 1
                  ICODCL(IFAC,IPHI(IPHAS)) = 1
                  ICODCL(IFAC,IFB (IPHAS)) = 1
                ELSE
                  ICODCL(IFAC,IK  (IPHAS)) = 3
                  ICODCL(IFAC,IEP (IPHAS)) = 3
                  ICODCL(IFAC,IPHI(IPHAS)) = 3
                  ICODCL(IFAC,IFB (IPHAS)) = 3
                ENDIF
              ELSEIF(ITURB(IPHAS).EQ.60) THEN
               IF(RCODCL(IFAC,IK  (IPHAS),1).GT.0.D0.AND.
     &            RCODCL(IFAC,IOMG(IPHAS),1).GT.0.D0 ) THEN
                  ICODCL(IFAC,IK  (IPHAS)) = 1
                  ICODCL(IFAC,IOMG(IPHAS)) = 1
                ELSE
                  ICODCL(IFAC,IK  (IPHAS)) = 3
                  ICODCL(IFAC,IOMG(IPHAS)) = 3
                ENDIF
              ENDIF
              IF(NSCAUS.GT.0) THEN
                DO II = 1, NSCAUS
                  IF(RCODCL(IFAC,ISCA(II),1).GT.-RINFIN*0.5D0) THEN
                    ICODCL(IFAC,ISCA(II)) = 1
                  ELSE
                    ICODCL(IFAC,ISCA(II)) = 3
                  ENDIF
                ENDDO
              ENDIF
            ENDIF
C
C
C     Les RCODCL ont ete initialises a -RINFIN pour permettre de
C       verifier ceux que l'utilisateur a modifies. On les remet a zero
C       si l'utilisateur ne les a pas modifies.
C       On traite la turbulence et les scalaires passifs (pour
C       simplifier la boucle, on traite toutes les variables : les
C       variables du compressible sont donc vues deux fois, mais ce
C       n'est pas grave).
            DO IVAR = 1, NVAR
              IF(RCODCL(IFAC,IVAR,1).LE.-RINFIN*0.5D0) THEN
                RCODCL(IFAC,IVAR,1) = 0.D0
              ENDIF
            ENDDO
C
C
C --- Fin de test sur les faces d'entree sortie
          ENDIF
C
C --- Fin de boucle sur les faces de bord
        ENDDO
C
C --- Fin de boucle sur les phases
      ENDDO
C
C----
C FORMATS
C----
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION                        ',/,
     &'@    *********                                               ',/,
     &'@    Deux variables independantes et deux seulement parmi    ',/,
     &'@    P, rho, T et E doivent etre imposees aux bords de type  ',/,
     &'@    IESICF dans uscfcl (ICCFTH = ',I10,').                  ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les conditions aux limites dans uscfcl.          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1100 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION                        ',/,
     &'@    *********                                               ',/,
     &'@    La pression n''a pas ete fournie en sortie a pression   ',/,
     &'@    imposée.                                                ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les conditions aux limites dans uscfcl.          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1200 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION                        ',/,
     &'@    *********                                               ',/,
     &'@    La masse volumique ou la vitesse n''a pas été fournie   ',/,
     &'@    en entree a masse volumique et vitesse imposee.         ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les conditions aux limites dans uscfcl.          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1300 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION                        ',/,
     &'@    *********                                               ',/,
     &'@    Le debit massique ou le debit enthalpique n''a pas été  ',/,
     &'@    fourni en entree a debit massique et enthalpique imposé.',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les conditions aux limites dans uscfcl.          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1301 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION                        ',/,
     &'@    *********                                               ',/,
     &'@    Entree à debit massique et debit enthalpique non prevue ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Contacter l''equipe de developpement pour uscfcl.         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1400 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION                        ',/,
     &'@    *********                                               ',/,
     &'@    Une condition a la limite ne fait pas partie des        ',/,
     &'@      conditions aux limites predefinies en compressible.   ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier les conditions aux limites dans uscfcl.          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 7000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION                        ',/,
     &'@    *********                                               ',/,
     &'@    cfxtcl doit etre modifie pour prendre en compte une loi ',/,
     &'@      d''etat a gamma variable. Seul est pris en compte le  ',/,
     &'@      cas IEOS = 1                                          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier IEOS dans uscfth.                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C----
C FIN
C----
C
      RETURN
      END
c@z
