c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS
C
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
                        SUBROUTINE USCFTH
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
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
     &   SORTI1 , SORTI2 , GAMAGR , XMASMR ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ---------
c@foncb
CFONC
CFONC ROUTINE UTILISATEUR POUR ENTREE DES PARAMETRES THERMODYNAMIQUES
CFONC
CFONC
CFONC     CE SOUS PROGRAMME UTILISATEUR EST OBLIGATOIRE
CFONC     =============================================
CFONC
CFONC
CFONC
CFONC CALCUL DE LA PROPRIETE THERMODYNAMIQUE VOULUE
CFONC EN FONCTION DES PARAMETRES DONNES ET DE LA THERMO CHOISIE
CFONC
CFONC
CFONC LOIS THERMODYNAMIQUES IMPLEMENTEES :
CFONC ====================================
CFONC
CFONC   1. GAZ PARFAIT (fournir la masse molaire XMASML)
CFONC
CFONC   2. GAZ PARFAIT A GAMMA VARIABLE (à adapter)
CFONC
CFONC   3. VAN DER WAALS (pas encore implemente)
CFONC
CFONC
CFONC CALCULS IMPLEMENTES :
CFONC =====================
CFONC
CFONC   VARIABLE     P   rho   T    e    h    s  epsilon-CvT  C.L.
CFONC       CODE     1    2    3    4    5    6      7         9
CFONC
CFONC
CFONC   CODE DU CALCUL : *VARIABLE i : i
CFONC
CFONC                                       di|
CFONC                    *DERIVEE PARTIELLE --| : ijk
CFONC                                       dj|k
CFONC
CFONC                    *CALCUL DES VARIABLES A PARTIR DE i ET j : ij
CFONC
CFONC   Et pour un reperage automatique, on utilise les nombres premiers
CFONC            P   rho   T   e   s
CFONC            2     3   5   7   13  (* 10000 pour les calculs aux cellules)
CFONC
CFONC
CFONC
CFONC   -> OPTIONS DE CALCUL                           : ICCFTH = -1
CFONC
CFONC   -> INITIALISATIONS PAR DEFAUT                  : ICCFTH = 0
CFONC
CFONC   -> CALCUL DE GAMMA                             : ICCFTH = 1
CFONC
CFONC   -> VERIFICATION DE rho                         : ICCFTH = -2
CFONC
CFONC   -> VERIFICATION DE E                           : ICCFTH = -4
CFONC
CFONC   -> CALCUL DE T ET e EN FONCTION DE P ET rho    : ICCFTH = 12 ou  60000
CFONC
CFONC   -> CALCUL DE rho ET e EN FONCTION DE P ET T    : ICCFTH = 13 ou 100000
CFONC
CFONC   -> CALCUL DE rho ET T EN FONCTION DE P ET e    : ICCFTH = 14 ou 140000
CFONC
CFONC   -> CALCUL DE P ET e EN FONCTION DE rho ET T    : ICCFTH = 23 ou 150000
CFONC
CFONC   -> CALCUL DE P ET T EN FONCTION DE rho ET e    : ICCFTH = 24 ou 210000
CFONC
CFONC                 2    dP |
CFONC   -> CALCUL DE c  = ----|                        : ICCFTH = 126
CFONC                     drho|s
CFONC
CFONC                       dP|
CFONC   -> CALCUL DE beta = --|                        : ICCFTH = 162
CFONC                       ds|rho
CFONC
CFONC                     de|
CFONC   -> CALCUL DE Cv = --|                          : ICCFTH = 432
CFONC                     dT|rho
CFONC
CFONC   -> CALCUL DE L'ENTROPIE                        : ICCFTH = 6
CFONC
CFONC
CFONC   -> CALCUL DE epsilon - Cv.T                    : ICCFTH = 7
CFONC
CFONC
CFONC   -> CALCUL DES CONDITIONS AUX LIMITES
CFONC               - SYMETRIE                         : ICCFTH = 90
CFONC               - PAROI                            : ICCFTH = 91
CFONC               - ENTREE                           : ICCFTH = 92
CFONC               - SORTIE                           : ICCFTH = 93
CFONC               - T ET e EN FONCTION DE P ET rho   : ICCFTH = 912 ou  60900
CFONC               - rho ET e EN FONCTION DE P ET T   : ICCFTH = 913 ou 100900
CFONC               - rho ET T EN FONCTION DE P ET e   : ICCFTH = 914 ou 140900
CFONC               - P ET e EN FONCTION DE rho ET T   : ICCFTH = 923 ou 150900
CFONC               - P ET T EN FONCTION DE rho ET e   : ICCFTH = 924 ou 210900
CFONC
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
CARGU ! ICCFTH       ! E  !  -> ! NUMERO DU TYPE DE CALCUL DEMANDE     !
CARGU ! IMODIF       ! E  !  -> ! SI CALCUL AUX CELLULES : MODIF DE RTP!
CARGU !              !    !     !   QUAND IMODIF > 0 (EN PARTICULIER   !
CARGU !              !    !     !   AUTOMATISE L'INIT (VOIR USCFXI)    !
CARGU !              !    !     ! SI CALCUL AUX FACES : NUMERO DE FACE !
CARGU !              !    !     !   AUTOMATISE L'INIT (VOIR USCFXI)    !
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE LA PHASE COURANTE          !
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
CARGU ! DT(NCELET)   ! TR !  -> ! VALEUR DU PAS DE TEMPS               !
CARGU ! RTP,RTPA     ! TR ! <-> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA COEFB  ! TR ! <-> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! SORTI1,2(*)  ! TR ! <-  ! VARIABLES DE SORTIE                  !
CARGU !              !    !     !   (INUTILISE SI ICCFTH.LT.0)         !
CARGU ! GAMAGR(*)    ! TR ! <-  ! CONSTANTE GAMMA EQUIVALENT DU GAZ    !
CARGU !              !    !     !   (INUTILISE SI ICCFTH.LT.0)         !
CARGU !              !    !     !   (PREMIERE CASE UTILISEE EN G.P.)   !
CARGU ! XMASMR(*)    ! TR ! <-  ! MASSE MOLAIRE DES CONSTITUANTS DU GAZ!
CARGU !              !    !     !   (INUTILISE SI ICCFTH.LT.0)         !
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
      INTEGER          ICCFTH   , IMODIF , IPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*),PROPFA(NFAC,*),PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
C
      DOUBLE PRECISION SORTI1(*), SORTI2(*), GAMAGR(*), XMASMR(*)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IIPH   , IFAC0
      INTEGER          IERR
      INTEGER          IEL    , IFAC   , IVAR
      INTEGER          IPRIPH , IRHIPH , ITKIPH , IENIPH
      INTEGER          IUIPH  , IVIPH  , IWIPH
      INTEGER          ICLP   , ICLR   , ICLT   , ICLE
      INTEGER          ICLU   , ICLV   , ICLW
      INTEGER          IUTILE
      DOUBLE PRECISION GAMAGP , XMASML , ENINT
      DOUBLE PRECISION XMACH  , XMACHI , XMACHE , DXMACH
C
      INTEGER          NPMAX
      PARAMETER (NPMAX = 1000)
      DOUBLE PRECISION CSTGR(NPMAX)
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C
      IF(1.EQ.1) RETURN
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C 0. INITIALISATION
C     Pas d'intervention utilisateur requise.
C=======================================================================
C
C     Pointeurs mémoire
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C     Indicateur d'erreur (stop si non nul)
      IERR   = 0
C
C     Position des variables
      IF(ICCFTH.GE.0.OR.ICCFTH.LE.-2) THEN
        IPRIPH = IPR(IPHAS)
        IRHIPH = ISCA(IRHO  (IPHAS))
        ITKIPH = ISCA(ITEMPK(IPHAS))
        IENIPH = ISCA(IENERG(IPHAS))
        IUIPH = IU(IPHAS)
        IVIPH = IV(IPHAS)
        IWIPH = IW(IPHAS)
        ICLP = ICLRTP(IPRIPH,ICOEF)
        ICLR = ICLRTP(IRHIPH,ICOEF)
        ICLT = ICLRTP(ITKIPH,ICOEF)
        ICLE = ICLRTP(IENIPH,ICOEF)
        ICLU = ICLRTP(IUIPH,ICOEF)
        ICLV = ICLRTP(IVIPH,ICOEF)
        ICLW = ICLRTP(IWIPH,ICOEF)
      ENDIF
C
C
C     Si calcul aux cellules,
C       indique, si > 0, que RTP doit etre modifié
C     Si calcul aux faces,
C       indique le numéro de la face considérée
      IFAC0 = IMODIF
C
C=======================================================================
C 1. CHOIX DE LA THERMODYNAMIQUE POUR CHAQUE PHASE
C     Intervention utilisateur requise.
C=======================================================================
C
C --> IEOS = 1 : Gaz Parfait avec Gamma constant
C --> IEOS = 2 : Gaz Parfait avec Gamma variable
C --> IEOS = 3 : Equation d'etat de Van Der Waals
C
      DO IIPH = 1, NPHAS
        IEOS(IIPH) = 1
      ENDDO
C
C
C            NE PAS OUBLIER DE RENSEIGNER LES PARAMETRES
C                  DANS LA THERMODYNAMIQUE CHOISIE
C            ===========================================
C
C
C=======================================================================
C 2.  GAZ PARFAIT
C=======================================================================
C
      IF(IEOS(IPHAS).EQ.1) THEN
C
C=======================================================================
C 2.1. PARAMETRES A RENSEIGNER PAR L'UTILISATEUR
C     Intervention utilisateur requise.
C=======================================================================
C
C     Pour chaque phase
C
C     Masse molaire du gaz (en kg/mol)
C
        IF(ICCFTH.GE.0) THEN
C
          IF(IPHAS.EQ.1) THEN
C
            XMASML = 28.8D-3
C
          ENDIF
C
        ENDIF
C
C=======================================================================
C 2.2. LOIS IMPLEMENTEES PAR DEFAUT
C     Pas d'intervention utilisateur requise.
C=======================================================================
C
C     CALCUL DE LA CONSTANTE GAMAGP
C     =============================
C
C     On la suppose supérieure ou égale à 1.
C     On la calcule à chaque appel, même si cela peut paraître cher, dans
C       par cohérence avec le cas gamma variable, pour lequel, on n'a pas
C       prévu de la conserver. Avec un save et un test, on pourrait
C       s'affranchir du calcul si nécessaire.
C
        IF(ICCFTH.GT.0) THEN
C
          GAMAGP = 1.D0 + RR/(XMASML*CP0(IPHAS)-RR)
C
          IF(GAMAGP.LT.1.D0) THEN
            WRITE(NFECRA,1010) GAMAGP
            CALL CSEXIT (1)
          ENDIF
C
C     On renvoie gamma si demande
          IF(ICCFTH.EQ.1) THEN
            GAMAGR(1) = GAMAGP
          ENDIF
C
        ENDIF
C
C
C     OPTIONS DE CALCUL  : Cp ET Cv CONSTANTS (hypothèse gaz parfait)
C     =================
C
C     La valeur de CP0 est à fournir dans usini1. La valeur de CV0
C       est calculée plus bas (à partir de CP0).
C
C
        IF(ICCFTH.EQ.-1) THEN
C
          ICP(IPHAS) = 0
          ICV(IPHAS) = 0
C
C
C     INITIALISATIONS PAR DEFAUT (avant uscfxi)
C     ==========================
C     T0 est forcément positif (vérifié dans verini)
C
        ELSEIF(ICCFTH.EQ.0) THEN
C
          CV0(IPHAS) = CP0(IPHAS) - RR/XMASML
C
          IF ( ISUITE .EQ. 0 ) THEN
            DO IEL = 1, NCEL
              RTP(IEL,IRHIPH) = P0(IPHAS)*XMASML/(RR*T0(IPHAS))
              RTP(IEL,IENIPH) = CV0(IPHAS)*T0(IPHAS)
            ENDDO
          ENDIF
C
C
C     VERIFICATION DE rho
C     ===================
C
        ELSEIF(ICCFTH.EQ.-2) THEN
C
C
C --- Clipping si rho est exactement nul + write + stop
C       En effet, si c'est le cas, on va nécessairement se planter
C       dans la thermo (bon, d'accord, pas dans tous les cas, mais
C       pour le moment, je présume que les fluides à traiter seront
C       assez classiques, avec des pbs à rho = 0)
C       Appelé en fin de resolution de la masse volumique (apres
C       clipping classique et avant communication parallele).
C
          IERR = 0
          DO IEL = 1, NCEL
            IF(RTP(IEL,IRHIPH).LE.0.D0) THEN
              RTP(IEL,IRHIPH) = EPZERO
              IERR = IERR + 1
            ENDIF
          ENDDO
          IF(IRANGP.GE.0) THEN
            CALL PARCPT (IERR)
          ENDIF
          IF(IERR.GT.0) THEN
            NTMABS = NTCABS
            WRITE(NFECRA,8000)IERR, EPZERO
          ENDIF
C
C
C     VERIFICATION DE e
C     ===================
C
        ELSEIF(ICCFTH.EQ.-4) THEN
C
C
C --- Clipping si e est exactement nul + write + stop
C       En effet, si c'est le cas, on va nécessairement se planter
C       dans la thermo (bon, d'accord, pas dans tous les cas, mais
C       pour le moment, je présume que les fluides à traiter seront
C       assez classiques, avec des pbs à e = 0)
C
          IERR = 0
          DO IEL = 1, NCEL
            ENINT = RTP(IEL,IENIPH)
     &               - 0.5D0*( RTP(IEL,IUIPH)**2
     &                       + RTP(IEL,IVIPH)**2
     &                       + RTP(IEL,IWIPH)**2 )
            IF(ENINT.LE.0.D0) THEN
              RTP(IEL,IENIPH) = EPZERO
     &               + 0.5D0*( RTP(IEL,IUIPH)**2
     &                       + RTP(IEL,IVIPH)**2
     &                       + RTP(IEL,IWIPH)**2 )
              IERR = IERR + 1
            ENDIF
          ENDDO
          IF(IRANGP.GE.0) THEN
            CALL PARCPT (IERR)
          ENDIF
          IF(IERR.GT.0) THEN
            NTMABS = NTCABS
            WRITE(NFECRA,8100)IERR, EPZERO
          ENDIF
C
C
C     CALCUL DE T ET e EN FONCTION DE P ET rho :
C     ========================================
C
        ELSEIF(ICCFTH.EQ.12.OR.ICCFTH.EQ.60000) THEN
C
C     Test des valeurs de rho
          IERR = 0
          DO IEL = 1, NCEL
            IF(RTP(IEL,IRHIPH).LE.0.D0) THEN
              WRITE(NFECRA,3010)RTP(IEL,IRHIPH),IEL
            ENDIF
          ENDDO
C     Si pb on s'arrete, car rho a été fourni par l'utilisateur
C       (initialisation erronée sans doute)
          IF(IERR.EQ.1) THEN
            CALL CSEXIT (1)
          ENDIF
C
          DO IEL = 1, NCEL
C     Temperature
            SORTI1(IEL) = XMASML*RTP(IEL,IPRIPH)/(RR*RTP(IEL,IRHIPH))
C     Energie totale
            SORTI2(IEL) = CV0(IPHAS)*SORTI1(IEL)
     &           + 0.5D0*( RTP(IEL,IUIPH)**2 + RTP(IEL,IVIPH)**2
     &                                       + RTP(IEL,IWIPH)**2 )
          ENDDO
C
C     Affectation a RTP
          IF(IMODIF.GT.0) THEN
            DO IEL = 1, NCEL
              RTP(IEL,ITKIPH) = SORTI1(IEL)
              RTP(IEL,IENIPH) = SORTI2(IEL)
            ENDDO
          ENDIF
C
C
C     CALCUL DE rho ET e EN FONCTION DE P ET T :
C     ========================================
C
        ELSEIF(ICCFTH.EQ.13.OR.ICCFTH.EQ.100000) THEN
C
C     Test des valeurs de T
          IERR = 0
          DO IEL = 1, NCEL
            IF(RTP(IEL,ITKIPH).LE.0.D0) THEN
              WRITE(NFECRA,2010)RTP(IEL,ITKIPH),IEL
            ENDIF
          ENDDO
C     Si pb on s'arrete, car T a été fourni par l'utilisateur
C       (initialisation erronée par exemple, non en K ?)
          IF(IERR.EQ.1) THEN
            CALL CSEXIT (1)
          ENDIF
C
          DO IEL = 1, NCEL
C     Masse volumique
            SORTI1(IEL) = XMASML*RTP(IEL,IPRIPH)/(RR*RTP(IEL,ITKIPH))
C     Energie totale
            SORTI2(IEL) = CV0(IPHAS)*RTP(IEL,ITKIPH)
     &           + 0.5D0*( RTP(IEL,IUIPH)**2 + RTP(IEL,IVIPH)**2
     &                                       + RTP(IEL,IWIPH)**2 )
          ENDDO
C
C     Affectation a RTP
          IF(IMODIF.GT.0) THEN
            DO IEL = 1, NCEL
              RTP(IEL,IRHIPH) = SORTI1(IEL)
              RTP(IEL,IENIPH) = SORTI2(IEL)
            ENDDO
          ENDIF
C
C
C     CALCUL DE rho ET T EN FONCTION DE P ET e :
C     ========================================
C
        ELSEIF(ICCFTH.EQ.14.OR.ICCFTH.EQ.140000) THEN
C
          DO IEL = 1, NCEL
C     Energie interne (évite de diviser par T)
            ENINT = RTP(IEL,IENIPH)
     &               - 0.5D0*( RTP(IEL,IUIPH)**2
     &                       + RTP(IEL,IVIPH)**2
     &                       + RTP(IEL,IWIPH)**2 )
C     Masse volumique
            SORTI1(IEL) = RTP(IEL,IPRIPH) / ( (GAMAGP-1.D0) * ENINT )
C     Temperature
            SORTI2(IEL) = XMASML * (GAMAGP-1.D0) * ENINT / RR
          ENDDO
C
C     Affectation a RTP
          IF(IMODIF.GT.0) THEN
            DO IEL = 1, NCEL
              RTP(IEL,IRHIPH) = SORTI1(IEL)
              RTP(IEL,ITKIPH) = SORTI2(IEL)
            ENDDO
          ENDIF
C
C
C     CALCUL DE P ET e EN FONCTION DE rho ET T :
C     ========================================
C
        ELSEIF(ICCFTH.EQ.23.OR.ICCFTH.EQ.150000) THEN
C
          DO IEL = 1, NCEL
C     Pression
            SORTI1(IEL) = RTP(IEL,IRHIPH)*RTP(IEL,ITKIPH)*RR/XMASML
C     Energie totale
            SORTI2(IEL) = CV0(IPHAS)*RTP(IEL,ITKIPH)
     &           + 0.5D0*( RTP(IEL,IUIPH)**2 + RTP(IEL,IVIPH)**2
     &                                       + RTP(IEL,IWIPH)**2 )
          ENDDO
C
C     Affectation a RTP
          IF(IMODIF.GT.0) THEN
            DO IEL = 1, NCEL
              RTP(IEL,IPRIPH) = SORTI1(IEL)
              RTP(IEL,IENIPH) = SORTI2(IEL)
            ENDDO
          ENDIF
C
C
C     CALCUL DE P ET T EN FONCTION DE rho ET e :
C     ========================================
C
        ELSEIF(ICCFTH.EQ.24.OR.ICCFTH.EQ.210000) THEN
C
          DO IEL = 1, NCEL
C     Energie interne (évite de diviser par T)
            ENINT = RTP(IEL,IENIPH)
     &               - 0.5D0*( RTP(IEL,IUIPH)**2
     &                       + RTP(IEL,IVIPH)**2
     &                       + RTP(IEL,IWIPH)**2 )
C     Pression
            SORTI1(IEL) = (GAMAGP-1.D0) * RTP(IEL,IRHIPH) * ENINT
C     Temperature
            SORTI2(IEL) = XMASML * (GAMAGP-1.D0) * ENINT / RR
          ENDDO
C
C       Affectation a RTP
          IF(IMODIF.GT.0) THEN
            DO IEL = 1, NCEL
              RTP(IEL,IPRIPH) = SORTI1(IEL)
              RTP(IEL,ITKIPH) = SORTI2(IEL)
            ENDDO
          ENDIF
C
C                2                            2          P
C     CALCUL DE c  EN FONCTION DE P ET rho : c  = gamma*---
C     ====================================              rho
C
        ELSEIF(ICCFTH.EQ.126) THEN
C
C     Test des valeurs de rho (on peut éliminer ce test pour optimiser)
          IUTILE = 0
          IF(IUTILE.EQ.1) THEN
            IERR = 0
            DO IEL = 1, NCEL
              IF(RTP(IEL,IRHIPH).LE.0.D0) THEN
                WRITE(NFECRA,4010)RTP(IEL,IRHIPH),IEL
              ENDIF
            ENDDO
            IF(IERR.EQ.1) THEN
              CALL CSEXIT (1)
            ENDIF
          ENDIF
C
          DO IEL = 1, NCEL
            SORTI1(IEL) = GAMAGP * RTP(IEL,IPRIPH) / RTP(IEL,IRHIPH)
          ENDDO
C
C                                                    gamma
C CALCUL DE beta EN FONCTION DE P ET rho : beta = rho
C ======================================
C
        ELSEIF(ICCFTH.EQ.162) THEN
C
C     Test des valeurs de rho (on peut éliminer ce test pour optimiser)
          IUTILE = 0
          IF(IUTILE.EQ.1) THEN
            IERR = 0
            DO IEL = 1, NCEL
              IF(RTP(IEL,IRHIPH).LT.0.D0) THEN
                WRITE(NFECRA,4020)RTP(IEL,IRHIPH),IEL
              ENDIF
            ENDDO
            IF(IERR.EQ.1) THEN
              CALL CSEXIT (1)
            ENDIF
          ENDIF
C
          DO IEL = 1, NCEL
            SORTI1(IEL) = RTP(IEL,IRHIPH)**GAMAGP
          ENDDO
C
C
C     CALCUL DE LA CHALEUR MASSIQUE A VOLUME CONSTANT
C     ===============================================
C
C     elle est constante !
C
C                                                           P
C     CALCUL DE L'ENTROPIE EN FONCTION DE P ET rho : s = --------
C     ============================================          gamma
C                                                        rho
        ELSEIF(ICCFTH.EQ.6) THEN
C
C     Test des valeurs de rho (on peut éliminer ce test pour optimiser)
          IERR = 0
          DO IEL = 1, NCEL
            IF(RTP(IEL,IRHIPH).LE.0.D0) THEN
              WRITE(NFECRA,4030)RTP(IEL,IRHIPH),IEL
            ENDIF
          ENDDO
          IF(IERR.EQ.1) THEN
            CALL CSEXIT (1)
          ENDIF
C
          DO IEL = 1, NCEL
            SORTI1(IEL) = RTP(IEL,IPRIPH) / (RTP(IEL,IRHIPH)**GAMAGP)
          ENDDO
C
C
C     CALCUL DE epsilon - Cv.T : epsilon - Cv.T = 0
C     ========================
C
        ELSEIF(ICCFTH.EQ.7) THEN
C
C --- A l'interieur du domaine
C
          DO IEL = 1, NCEL
            SORTI1(IEL) = 0.D0
          ENDDO
C
C --- Sur les bords
C
          DO IFAC = 1, NFABOR
            SORTI2(IFAC) = 0.D0
          ENDDO
C
C
C     CALCUL DES CONDITIONS AUX LIMITES : sur la face IFAC = IFAC0
C     =================================
C
C     PAROI
C     -----
C
        ELSEIF(ICCFTH.EQ.91) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C --- Calcul du nombre de Mach normal a la frontiere
C
          XMACH =
     &         ( RTP(IEL,IUIPH)*SURFBO(1,IFAC)
     &         + RTP(IEL,IVIPH)*SURFBO(2,IFAC)
     &         + RTP(IEL,IWIPH)*SURFBO(3,IFAC) ) / RA(ISRFBN+IFAC-1)
     &         / SQRT( GAMAGP * RTP(IEL,IPRIPH) / RTP(IEL,IRHIPH) )
C
C --- Pression
C
C     On utilise un Neumann en pression. Bien que cela ne permette
C       pas d'utiliser Rusanov, on espere que cela stabilise.
C     On ajoute un test pour eviter de basculer
C       de detente a choc d'un pas de temps à l'autre
C
C     Detente
          IF(XMACH.LT.0.D0.AND.COEFB(IFAC,ICLP).LE.1.D0) THEN
C
            IF(XMACH.GT.2.D0/(1.D0-GAMAGP)) THEN
              COEFB(IFAC,ICLP) = (1.D0 + (GAMAGP-1.D0)/2.D0 * XMACH)
     &             ** (2.D0*GAMAGP/(GAMAGP-1.D0))
C         Si on detend trop fort : Dirichlet nul en pression
C           (la valeur de COEFB ici est utilisee comme indicateur
C            et sera retraitée plus tard dans cfxtcl)
            ELSE
              COEFB(IFAC,ICLP) = RINFIN
            ENDIF
C
C     Choc
          ELSEIF(XMACH.GT.0.D0.AND.COEFB(IFAC,ICLP).GE.1.D0) THEN
C
            COEFB(IFAC,ICLP) = 1.D0 + GAMAGP*XMACH
     &            *( (GAMAGP+1.D0)/4.D0*XMACH
     &                + SQRT(1.D0 + (GAMAGP+1.D0)**2/16.D0*XMACH**2) )
C
C     Oscillation ou Mach nul
          ELSE
            COEFB(IFAC,ICLP) = 1.D0
          ENDIF

C
C     SYMETRIE
C     --------
C
        ELSEIF(ICCFTH.EQ.90) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C --- On impose un Neumann Homogene en Pression par defaut
C         Pas d'intervention requise
C
C
C --- ENTREE SUBSONIQUE, RHO ET U DONNES (subsonique par hypothèse)
C     ----------------------------------
C
C     Contrairement à la conception initiale, on impose un Dirichlet
C       explicite sur la pression et non pas un flux (on se base cpdt
C       sur la meme valeur physique).
C     Ceci a l'avantage de permettre l'utilisation de Rusanov pour
C       lisser les conditions d'entree.
C     En outre, ceci permet de ne pas avoir à gérer le remplissage de
C       coefb ici (c'est secondaire, car il faut le faire pour la paroi).
C     Si on a des oscillations en temps, on pourra essayer d'eviter
C       le basculement de détente à choc d'un pas de temps à l'autre.
C     La pertinence de cette démarche reste à démontrer.
C
        ELSEIF(ICCFTH.EQ.92) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C --- Calcul du nombre de Mach entrant normalement a la frontiere
C
          XMACHI =
     &         ( RTP(IEL,IUIPH)*SURFBO(1,IFAC)
     &         + RTP(IEL,IVIPH)*SURFBO(2,IFAC)
     &         + RTP(IEL,IWIPH)*SURFBO(3,IFAC) ) / RA(ISRFBN+IFAC-1)
     &         / SQRT( GAMAGP * RTP(IEL,IPRIPH) / RTP(IEL,IRHIPH) )
          XMACHE =
     &         ( COEFA(IFAC,ICLU)*SURFBO(1,IFAC)
     &         + COEFA(IFAC,ICLV)*SURFBO(2,IFAC)
     &         + COEFA(IFAC,ICLW)*SURFBO(3,IFAC) ) /RA(ISRFBN+IFAC-1)
     &         / SQRT( GAMAGP * RTP(IEL,IPRIPH) / RTP(IEL,IRHIPH) )
          DXMACH = XMACHI - XMACHE
C
C --- Pression : Detente
          IF(DXMACH.LE.0.D0) THEN
C
            IF(DXMACH.GT.2.D0/(1.D0-GAMAGP)) THEN
              COEFA(IFAC,ICLP) = RTP(IEL,IPRIPH)*
     &             ( (1.D0 + (GAMAGP-1.D0)*0.50D0*DXMACH)
     &               ** (2.D0*GAMAGP/(GAMAGP-1.D0))    )
            ELSEIF(DXMACH.LE.2.D0/(1.D0-GAMAGP) ) THEN
              COEFA(IFAC,ICLP) = 0.D0
            ENDIF
C
C --- Pression : Choc
          ELSE
            COEFA(IFAC,ICLP) = RTP(IEL,IPRIPH)*
     &           (  1.D0 + GAMAGP*DXMACH
     &           *( (GAMAGP+1.D0)*0.25D0*DXMACH
     &           + SQRT(1.D0 + (GAMAGP+1.D0)**2/16.D0*DXMACH**2) )  )
          ENDIF
C
C --- Energie totale
          COEFA(IFAC,ICLP) = RTP(IEL,IPRIPH)
          COEFA(IFAC,ICLE) =
     &         COEFA(IFAC,ICLP)/((GAMAGP-1.D0)*COEFA(IFAC,ICLR))
     &         + 0.5D0*(COEFA(IFAC,ICLU)**2
     &                + COEFA(IFAC,ICLV)**2 + COEFA(IFAC,ICLW)**2)
C
C
C --- ENTREE SUBSONIQUE, RHO*U ET RHO*U*H DONNES (subsonique par hypothèse)
C     ------------------------------------------
C
C     Ceci reste a completer :
        ELSEIF(ICCFTH.EQ.94) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
          WRITE(NFECRA,7000)
C
          CALL CSEXIT (1)
C         ===========
C
C     P   : a calculer avec un newton
C     u et rho en fonction de P
C     E en fonction de H
C     (je l'ai ecrit, reste à le coder ...)
C
C
C     SORTIE SUBSONIQUE (subsonique par hypothèse)
C     -----------------
C
        ELSEIF(ICCFTH.EQ.93) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C --- Detente :
          IF(COEFA(IFAC,ICLP).LE.RTP(IEL,IPRIPH)) THEN
C
C       Masse volumique
            COEFA(IFAC,ICLR) = RTP(IEL,IRHIPH)
     &           * (COEFA(IFAC,ICLP)/RTP(IEL,IPRIPH))**(1.D0/GAMAGP)
C
C       Vitesse
            COEFA(IFAC,ICLU) = RTP(IEL,IUIPH)
     &           + 2.D0/(GAMAGP-1.D0)
     &           * SQRT(GAMAGP*RTP(IEL,IPRIPH)/RTP(IEL,IRHIPH))
     &           * (1.D0-(COEFA(IFAC,ICLP)/RTP(IEL,IPRIPH)
     &                        )**((GAMAGP-1.D0)/(2.D0*GAMAGP)))
     &           * SURFBO(1,IFAC)/RA(ISRFBN+IFAC-1)
C
            COEFA(IFAC,ICLV) = RTP(IEL,IVIPH)
     &           + 2.D0/(GAMAGP-1.D0)
     &           * SQRT( GAMAGP*RTP(IEL,IPRIPH)/RTP(IEL,IRHIPH))
     &           * (1.D0-(COEFA(IFAC,ICLP)/RTP(IEL,IPRIPH)
     &                        )**((GAMAGP-1.D0)/(2.D0*GAMAGP)))
     &           * SURFBO(2,IFAC)/RA(ISRFBN+IFAC-1)
C
            COEFA(IFAC,ICLW) = RTP(IEL,IWIPH)
     &           + 2.D0/(GAMAGP-1.D0)
     &           * SQRT( GAMAGP*RTP(IEL,IPRIPH)/RTP(IEL,IRHIPH))
     &           * (1.D0-(COEFA(IFAC,ICLP)/RTP(IEL,IPRIPH)
     &                        )**((GAMAGP-1.D0)/(2.D0/GAMAGP)))
     &           * SURFBO(3,IFAC)/RA(ISRFBN+IFAC-1)
C
C       Energie totale
            COEFA(IFAC,ICLE) =
     &           COEFA(IFAC,ICLP)/((GAMAGP-1.D0)*COEFA(IFAC,ICLR))
     &           + 0.5D0*(COEFA(IFAC,ICLU)**2
     &                  + COEFA(IFAC,ICLV)**2 + COEFA(IFAC,ICLW)**2)
C
C --- Choc :
          ELSE
C
C       Masse volumique
            COEFA(IFAC,ICLR) = RTP(IEL,IRHIPH)
     &           * ( (GAMAGP+1.D0)*COEFA(IFAC,ICLP)
     &             + (GAMAGP-1.D0)*RTP(IEL,IPRIPH) )
     &           / ( (GAMAGP-1.D0)*COEFA(IFAC,ICLP)
     &             + (GAMAGP+1.D0)*RTP(IEL,IPRIPH) )
C
C       Vitesse
            COEFA(IFAC,ICLU) = RTP(IEL,IUIPH)
     &           - (COEFA(IFAC,ICLP)-RTP(IEL,IPRIPH))
     &           * SQRT(2.D0/
     &                  (RTP(IEL,IRHIPH)
     &                   *((GAMAGP+1.D0)*COEFA(IFAC,ICLP)
     &                    +(GAMAGP-1.D0)*RTP(IEL,IPRIPH) )))
     &           * SURFBO(1,IFAC)/RA(ISRFBN+IFAC-1)
C
            COEFA(IFAC,ICLV) = RTP(IEL,IVIPH)
     &           - (COEFA(IFAC,ICLP)-RTP(IEL,IPRIPH))
     &           * SQRT(2.D0/
     &                  (RTP(IEL,IRHIPH)
     &                   *((GAMAGP+1.D0)*COEFA(IFAC,ICLP)
     &                    +(GAMAGP-1.D0)*RTP(IEL,IPRIPH) )))
     &           * SURFBO(2,IFAC)/RA(ISRFBN+IFAC-1)
C
            COEFA(IFAC,ICLW) = RTP(IEL,IWIPH)
     &           - (COEFA(IFAC,ICLP)-RTP(IEL,IPRIPH))
     &           * SQRT(2.D0/
     &                  (RTP(IEL,IRHIPH)
     &                   *((GAMAGP+1.D0)*COEFA(IFAC,ICLP)
     &                    +(GAMAGP-1.D0)*RTP(IEL,IPRIPH) )))
     &           * SURFBO(3,IFAC)/RA(ISRFBN+IFAC-1)
C
C       Energie totale
            COEFA(IFAC,ICLE) =
     &           COEFA(IFAC,ICLP)/((GAMAGP-1.D0)*COEFA(IFAC,ICLR))
     &           + 0.5D0*(COEFA(IFAC,ICLU)**2
     &                  + COEFA(IFAC,ICLV)**2 + COEFA(IFAC,ICLW)**2)
C
          ENDIF
C
C
C     T ET e EN FONCTION DE P ET rho
C     ------------------------------
C      on a donné des valeurs positives strictement (hypothèse)
C
        ELSEIF(ICCFTH.EQ.912.OR.ICCFTH.EQ.60900) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C     Temperature
          COEFA(IFAC,ICLT) =
     &         XMASML*COEFA(IFAC,ICLP)/(RR*COEFA(IFAC,ICLR))
C
C     Energie totale
          COEFA(IFAC,ICLE) =
     &         CV0(IPHAS)*COEFA(IFAC,ICLT)
     &         + 0.5D0*( COEFA(IFAC,ICLU)**2
     &                 + COEFA(IFAC,ICLV)**2 + COEFA(IFAC,ICLW)**2 )
C
C
C     rho ET e EN FONCTION DE P ET T
C     ------------------------------
C
        ELSEIF(ICCFTH.EQ.913.OR.ICCFTH.EQ.100900) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C     Masse volumique
          COEFA(IFAC,ICLR) =
     &         XMASML*COEFA(IFAC,ICLP)/(RR*COEFA(IFAC,ICLT))
C
C     Energie totale
          COEFA(IFAC,ICLE) =
     &         CV0(IPHAS)*COEFA(IFAC,ICLT)
     &         + 0.5D0*( COEFA(IFAC,ICLU)**2
     &                 + COEFA(IFAC,ICLV)**2 + COEFA(IFAC,ICLW)**2 )
C
C
C     rho ET T EN FONCTION DE P ET e
C     ------------------------------
C
        ELSEIF(ICCFTH.EQ.914.OR.ICCFTH.EQ.140900) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C     Masse volumique
          COEFA(IFAC,ICLR) = COEFA(IFAC,ICLP)/( (GAMAGP-1.D0)*
     &         (COEFA(IFAC,ICLE)
     &         - 0.5D0*( COEFA(IFAC,ICLU)**2
     &                 + COEFA(IFAC,ICLV)**2
     &                 + COEFA(IFAC,ICLW)**2 ) ) )
C
C     Temperature
          COEFA(IFAC,ICLT)=
     &         XMASML*COEFA(IFAC,ICLP)/(RR*COEFA(IFAC,ICLR))
C
C
C     P ET e EN FONCTION DE rho ET T
C     ------------------------------
C
        ELSEIF(ICCFTH.EQ.923.OR.ICCFTH.EQ.150900) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C     Pression
          COEFA(IFAC,ICLP) = COEFA(IFAC,ICLR)*RR/XMASML
     &                                       *COEFA(IFAC,ICLT)
C
C     Energie totale
          COEFA(IFAC,ICLE) = CV0(IPHAS) * COEFA(IFAC,ICLT)
     &         + 0.5D0*( COEFA(IFAC,ICLU)**2
     &                 + COEFA(IFAC,ICLV)**2 + COEFA(IFAC,ICLW)**2 )
C
C
C     P ET T EN FONCTION DE rho ET e
C     ------------------------------
C
        ELSEIF(ICCFTH.EQ.924.OR.ICCFTH.EQ.210900) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C     Pression
          COEFA(IFAC,ICLP) = (GAMAGP-1.D0)*COEFA(IFAC,ICLR)
     &          *( COEFA(IFAC,ICLE)
     &            - 0.5D0*( COEFA(IFAC,ICLU)**2
     &                    + COEFA(IFAC,ICLV)**2
     &                    + COEFA(IFAC,ICLW)**2 ) )

C
C     Temperature
          COEFA(IFAC,ICLT)=
     &         XMASML*COEFA(IFAC,ICLP)/(RR*COEFA(IFAC,ICLR))
C
C
C --- Fin de test sur GAZ PARFAIT
        ENDIF
C
C
C=======================================================================
C 2.  GAZ PARFAIT AVEC GAMMA VARIABLE
C=======================================================================
C
      ELSEIF(IEOS(IPHAS).EQ.2) THEN
C
C***********************************************************************
C
C     PARAMETRES
C     ==========
C
C-----------------------------------------------------------------------
C
C     EXEMPLES  (A recopier et a adapter apres 'PARAMETRES')
        IF(0.EQ.1) THEN
C-----------------------------------------------------------------------
C
C --> EXEMPLE 1 : GAZ PARFAIT A 3 CONSTITUANTS
C
C     Phase
        IF(IPHAS.EQ.1) THEN
C
C     Masses molaires des constituants (en kg/mol)
          CSTGR(1)  = 18.D-3
          CSTGR(2)  = 32.D-3
          CSTGR(3)  = 28.D-3
C
C     Calcul de la masse molaire du melange au centre des cellules
          IF(ICCFTH.GT.0) THEN
            DO IEL = 1, NCEL
              XMASMR(IEL) = 1.D0 / ( RTP(IEL,ISCA(1))/CSTGR(1)
     &                             + RTP(IEL,ISCA(2))/CSTGR(2)
     &                             + RTP(IEL,ISCA(3))/CSTGR(3) )
            ENDDO

C
C     Calcul du GAMMA equivalent au centre des cellules
            DO IEL = 1, NCEL
              GAMAGR(IEL) = PROPCE(IEL,IPPROC(ICP(IPHAS)))
     &           / ( PROPCE(IEL,IPPROC(ICP(IPHAS))) - RR/XMASMR(IEL) )
            ENDDO
          ENDIF
C
        ENDIF
C
C-----------------------------------------------------------------------
C
C     FIN DES EXEMPLES
        ENDIF
C-----------------------------------------------------------------------
C
C***********************************************************************
C
C
C TEST DE VALIDITE DE GAMAGR : GAMAGR >= 1
C ==========================
C
        IERR = 0
C
        DO IEL = 1, NCEL
          IF(ICCFTH.GT.0 .AND. GAMAGR(IEL).LT.1.D0) THEN
            IERR = 1
            WRITE(NFECRA,1020) IEL, GAMAGR(IEL)
          ENDIF
        ENDDO
C
        IF(IERR.EQ.1) THEN
          CALL CSEXIT (1)
        ENDIF
C
C
C OPTIONS DE CALCUL  : Cp ET Cv VARIABLES
C =================
C
        IF(ICCFTH.EQ.-1) THEN
C
          ICP(IPHAS) = 1
          CP0(IPHAS) = EPZERO
          ICV(IPHAS) = 1
          CV0(IPHAS) = EPZERO
C
C
C INITIALISATIONS PAR DEFAUT :
C ==========================
C
        ELSEIF(ICCFTH.EQ.0) THEN
C
          DO IEL = 1, NCEL
            PROPCE(IEL,IPPROC(ICP(IPHAS))) = CP0(IPHAS)
            PROPCE(IEL,IPPROC(ICV(IPHAS))) =
     &           CP0(IPHAS) - RR/XMASMR(IEL)
            RTP(IEL,IRHIPH) = P0(IPHAS)*XMASMR(IEL)/RR/T0(IPHAS)
            RTP(IEL,IENIPH) =
     &           PROPCE(IEL,IPPROC(ICV(IPHAS)))*T0(IPHAS)
          ENDDO
C
C
C CALCUL DE T ET e EN FONCTION DE P ET rho :
C ========================================
C
        ELSEIF(ICCFTH.EQ.12) THEN
C
          DO IEL = 1, NCEL
C
C     Temperature
            SORTI1(IEL) =
     &           XMASMR(IEL)/RR*RTP(IEL,IPRIPH)/RTP(IEL,IRHIPH)
C
C     Energie totale
            SORTI2(IEL) = PROPCE(IEL,IPPROC(ICV(IPHAS)))*SORTI1(IEL)
     &    + 0.5D0*( RTP(IEL,IUIPH)**2
     &           + RTP(IEL,IVIPH)**2 + RTP(IEL,IWIPH)**2 )
C
          ENDDO
C
C
          IF(IMODIF.GT.0) THEN
C       Affectation directe a RTP
            DO IEL = 1, NCEL
              RTP(IEL,ITKIPH) = SORTI1(IEL)
              RTP(IEL,IENIPH) = SORTI2(IEL)
            ENDDO
          ENDIF
C
C
C CALCUL DE rho ET e EN FONCTION DE P ET T :
C ========================================
C
        ELSEIF(ICCFTH.EQ.13) THEN
C
          DO IEL = 1, NCEL
C
C     Masse volumique
            SORTI1(IEL) =
     &           XMASMR(IEL)/RR*RTP(IEL,IPRIPH)/RTP(IEL,ITKIPH)
C
C     Energie totale
            SORTI2(IEL) =
     &           PROPCE(IEL,IPPROC(ICV(IPHAS)))*RTP(IEL,ITKIPH)
     &    + 0.5D0*( RTP(IEL,IUIPH)**2
     &           + RTP(IEL,IVIPH)**2 + RTP(IEL,IWIPH)**2 )
C
          ENDDO
C
C
          IF(IMODIF.GT.0) THEN
C       Affectation directe a RTP
            DO IEL = 1, NCEL
              RTP(IEL,IRHIPH) = SORTI1(IEL)
              RTP(IEL,IENIPH) = SORTI2(IEL)
            ENDDO
          ENDIF
C
C
C CALCUL DE rho ET T EN FONCTION DE P ET e :
C ========================================
C
        ELSEIF(ICCFTH.EQ.14) THEN
C
          DO IEL = 1, NCEL
C
C     Masse volumique
            SORTI1(IEL) =
     &           RTP(IEL,IPRIPH)/(GAMAGR(IEL)-1.D0)/( RTP(IEL,IENIPH)
     &  - 0.5D0*( RTP(IEL,IUIPH)**2
     &           + RTP(IEL,IVIPH)**2 + RTP(IEL,IWIPH)**2 ) )
C
C     Temperature
            SORTI2(IEL) = XMASMR(IEL)/RR*RTP(IEL,IPRIPH)/SORTI1(IEL)
C
          ENDDO
C
C
          IF(IMODIF.GT.0) THEN
C       Affectation directe a RTP
            DO IEL = 1, NCEL
              RTP(IEL,IRHIPH) = SORTI1(IEL)
              RTP(IEL,ITKIPH) = SORTI2(IEL)
            ENDDO
          ENDIF
C
C
C CALCUL DE P ET e EN FONCTION DE rho ET T :
C ========================================
C
        ELSEIF(ICCFTH.EQ.23) THEN
C
          DO IEL = 1, NCEL
C
C     Pression
            SORTI1(IEL) =
     &           RTP(IEL,IRHIPH)*RR/XMASMR(IEL)*RTP(IEL,ITKIPH)
C
C     Energie totale
            SORTI2(IEL) =
     &           PROPCE(IEL,IPPROC(ICV(IPHAS)))*RTP(IEL,ITKIPH)
     &    + 0.5D0*( RTP(IEL,IUIPH)**2
     &           + RTP(IEL,IVIPH)**2 + RTP(IEL,IWIPH)**2 )
C
          ENDDO
C
C
          IF(IMODIF.GT.0) THEN
C       Affectation directe a RTP
            DO IEL = 1, NCEL
              RTP(IEL,IPRIPH) = SORTI1(IEL)
              RTP(IEL,IENIPH) = SORTI2(IEL)
            ENDDO
          ENDIF
C
C
C CALCUL DE P ET T EN FONCTION DE rho ET e :
C ========================================
C
        ELSEIF(ICCFTH.EQ.24) THEN
C
          DO IEL = 1, NCEL
C
C     Pression
            SORTI1(IEL) =
     &           (GAMAGR(IEL)-1.D0)*RTP(IEL,IRHIPH)*( RTP(IEL,IENIPH)
     &  - 0.5D0*( RTP(IEL,IUIPH)**2
     &           + RTP(IEL,IVIPH)**2 + RTP(IEL,IWIPH)**2 ) )
C
C     Temperature
            SORTI2(IEL) = XMASMR(IEL)/RR*SORTI1(IEL)/RTP(IEL,IRHIPH)
C
          ENDDO
C
          IF(IMODIF.GT.0) THEN
C       Affectation directe a RTP
            DO IEL = 1, NCEL
              RTP(IEL,IPRIPH) = SORTI1(IEL)
              RTP(IEL,ITKIPH) = SORTI2(IEL)
            ENDDO
          ENDIF
C
C            2                            2          P
C CALCUL DE c  EN FONCTION DE P ET rho : c  = gamma*---
C ====================================              rho
C
        ELSEIF(ICCFTH.EQ.126) THEN
C
          DO IEL = 1, NCEL
C
C --- Test de positivite de P
            IF(RTP(IEL,IPRIPH).LT.0.D0) THEN
              WRITE(NFECRA,1110) IEL , RTP(IEL,IPRIPH)
              IERR = 1
C
C --- Test de positivite de rho
            ELSEIF(RTP(IEL,IRHIPH).LE.0.D0) THEN
              WRITE(NFECRA,1120) IEL , RTP(IEL,IRHIPH)
              IERR = 1
C
            ELSE
C
C --- Calcul
              SORTI1(IEL) =
     &             GAMAGR(IEL) * RTP(IEL,IPRIPH) / RTP(IEL,IRHIPH)
C
            ENDIF
C
          ENDDO
C
          IF(IERR.EQ.1) CALL CSEXIT (1)
C
C                                                    gamma
C CALCUL DE beta EN FONCTION DE P ET rho : beta = rho
C ======================================
C
        ELSEIF(ICCFTH.EQ.162) THEN
C
          DO IEL = 1, NCEL
C
C --- Test de positivite de rho
            IF(RTP(IEL,IRHIPH).LT.0.D0) THEN
              WRITE(NFECRA,1220) IEL , RTP(IEL,IRHIPH)
              IERR = 1
C
            ELSE
C
C --- Calcul
              SORTI1(IEL) = RTP(IEL,IRHIPH)**GAMAGR(IEL)
C
            ENDIF
C
          ENDDO
C
          IF(IERR.EQ.1) CALL CSEXIT (1)
C
C                           R
C CALCUL DE Cv : Cv = Cp - ---
C ============              M
C
        ELSEIF(ICCFTH.EQ.432) THEN
C
          DO IEL = 1, NCEL
C
            SORTI1(IEL) = PROPCE(IEL,IPPROC(ICP(IPHAS)))-RR/XMASMR(IEL)
C
          ENDDO
C
          IF(IERR.EQ.1) CALL CSEXIT (1)
C
C                                                       P
C CALCUL DE L'ENTROPIE EN FONCTION DE P ET rho : s = --------
C ============================================          gamma
C                                                    rho
        ELSEIF(ICCFTH.EQ.6) THEN
C
          DO IEL = 1, NCEL
C
C --- Test de positivite de P
            IF(RTP(IEL,IPRIPH).LT.0.D0) THEN
              WRITE(NFECRA,1310) IEL , RTP(IEL,IPRIPH)
              IERR = 1
C
C --- Test de positivite de rho
            ELSEIF(RTP(IEL,IRHIPH).LE.0.D0) THEN
              WRITE(NFECRA,1320) IEL , RTP(IEL,IRHIPH)
              IERR = 1
C
            ELSE
C
C --- Calcul
              SORTI1(IEL) =
     &             RTP(IEL,IPRIPH) / (RTP(IEL,IRHIPH)**GAMAGR(IEL))
C
            ENDIF
C
          ENDDO
C
          IF(IERR.EQ.1) CALL CSEXIT (1)
C
C
C CALCUL DE epsilon - Cv.T : epsilon - Cv.T = 0
C ========================
C
        ELSEIF(ICCFTH.EQ.7) THEN
C
C --- A l'interieur du domaine
C
          DO IEL = 1, NCEL
C
            SORTI1(IEL) = 0.D0
C
          ENDDO
C
C --- Sur les bords
C
          DO IFAC = 1, NFABOR
C
            SORTI2(IFAC) = 0.D0
C
          ENDDO
C
          IF(IERR.EQ.1) CALL CSEXIT (1)
C
C
C CALCUL DES CONDITIONS AUX LIMITES :
C =================================
C
C     PAROI/SYMETRIE
C     --------------
        ELSEIF(ICCFTH.EQ.91) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C --- Calcul du nombre de Mach normal a la frontiere
C
          XMACH = ( RTP(IEL,IUIPH)*SURFBO(1,IFAC)
     &           + RTP(IEL,IVIPH)*SURFBO(2,IFAC)
     &           + RTP(IEL,IWIPH)*SURFBO(3,IFAC) ) / RA(ISRFBN+IFAC-1)
     &         / SQRT( GAMAGR(IEL)*RTP(IEL,IPRIPH)/RTP(IEL,IRHIPH) )
C
C --- Pression / Entropie : Detente
C
          COEFA(IFAC,ICLP) = 0.D0
C
          IF(XMACH.LE.0.D0 .AND. XMACH.GT.2.D0/(1.D0-GAMAGR(IEL))) THEN
            COEFB(IFAC,ICLP) = (1.D0 + (GAMAGR(IEL)-1.D0)/2.D0 * XMACH)
     &           ** (2.D0*GAMAGR(IEL)/(GAMAGR(IEL)-1.D0))
            COEFB(IFAC,ICLT) = 1.D0
C
          ELSEIF(XMACH.LE.2.D0/(1.D0-GAMAGR(IEL)) ) THEN
            COEFB(IFAC,ICLP) = 0.D0
            COEFB(IFAC,ICLT) = 1.D0
C
C --- Pression / Entropie : Choc
C
          ELSE
            COEFB(IFAC,ICLP) = 1.D0 + GAMAGR(IEL)*XMACH
     &            *( (GAMAGR(IEL)+1.D0)/4.D0*XMACH
     &           + SQRT(1.D0 + (GAMAGR(IEL)+1.D0)**2/16.D0*XMACH**2) )
            COEFB(IFAC,ICLT) = COEFB(IFAC,ICLP)/(1.D0-COEFB(IFAC,ICLP))
     &          / RTP(IEL,IPRIPH) * ( RTP(IEL,IRHIPH)
     &              * (RTP(IEL,IUIPH)**2
     &                +RTP(IEL,IVIPH)**2+RTP(IEL,IWIPH)**2)
     &              + RTP(IEL,IPRIPH) *(1.D0-COEFB(IFAC,ICLP)) )
          ENDIF
C
C --- Energie totale : Epsilon sup
C
          COEFA(IFAC,ICLE) = 0.D0
C
          IF(IERR.EQ.1) CALL CSEXIT (1)
C
C     ENTREE
C     ------
        ELSEIF(ICCFTH.EQ.92) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C --- Calcul du nombre de Mach entrant normalement a la frontiere
C
          XMACHI = ( RTP(IEL,IUIPH)*SURFBO(1,IFAC)
     &         + RTP(IEL,IVIPH)*SURFBO(2,IFAC)
     &         + RTP(IEL,IWIPH)*SURFBO(3,IFAC) )/RA(ISRFBN+IFAC-1)
     &         / SQRT(GAMAGR(IEL)*RTP(IEL,IPRIPH)/RTP(IEL,IRHIPH))
          XMACHE = ( COEFA(IFAC,ICLU)*SURFBO(1,IFAC)
     &         + COEFA(IFAC,ICLV)*SURFBO(2,IFAC)
     &         + COEFA(IFAC,ICLW)*SURFBO(3,IFAC) )/RA(ISRFBN+IFAC-1)
     &         / SQRT(GAMAGR(IEL)*RTP(IEL,IPRIPH)/RTP(IEL,IRHIPH))
          DXMACH = XMACHI - XMACHE
C
C --- Pression : Detente
          IF(DXMACH.LE.0.D0) THEN
C
            IF(DXMACH.GT.2.D0/(1.D0-GAMAGR(IEL))) THEN
              COEFA(IFAC,ICLP) = RTP(IEL,IPRIPH)*
     &             ( (1.D0 + (GAMAGR(IEL)-1.D0)*0.50D0*DXMACH)
     &               ** (2.D0*GAMAGR(IEL)/(GAMAGR(IEL)-1.D0))  )
            ELSEIF(DXMACH.LE.2.D0/(1.D0-GAMAGR(IEL)) ) THEN
              COEFA(IFAC,ICLP) = 0.D0
            ENDIF
C
C --- Pression : Choc
          ELSE
            COEFA(IFAC,ICLP) = RTP(IEL,IPRIPH)*
     &           (  1.D0 + GAMAGR(IEL)*DXMACH
     &           *( (GAMAGR(IEL)+1.D0)*0.25D0*DXMACH
     &           + SQRT(1.D0 + (GAMAGR(IEL)+1.D0)**2/16.D0
     &                                           *DXMACH**2) )  )
          ENDIF
C
C --- Energie totale
C
          COEFA(IFAC,ICLP) = RTP(IEL,IPRIPH)
C
          COEFA(IFAC,ICLE) =
     &         COEFA(IFAC,ICLP)/((GAMAGR(IEL)-1.D0)*COEFA(IFAC,ICLR))
     &         + 0.5D0*(COEFA(IFAC,ICLU)**2
     &                + COEFA(IFAC,ICLV)**2 + COEFA(IFAC,ICLW)**2)
C
C     SORTIE
C     ------
        ELSEIF(ICCFTH.EQ.93) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C --- Calcul du nombre de Mach sortant normalement a la frontiere
C
          XMACH = ( RTP(IEL,IUIPH)*SURFBO(1,IFAC)
     &           + RTP(IEL,IVIPH)*SURFBO(2,IFAC)
     &           + RTP(IEL,IWIPH)*SURFBO(3,IFAC) ) / RA(ISRFBN+IFAC-1)
     &         / SQRT(GAMAGR(IEL)*RTP(IEL,IPRIPH)/RTP(IEL,IRHIPH))
C
C --- Sortie supersonique : Dirichlet sur toutes les variables
C     -------------------
          IF(XMACH.GE.1.D0) THEN
            DO IVAR = 1, NVAR
              COEFA(IFAC,ICLRTP(IVAR,ICOEF)) = RTP(IEL,IVAR)
            ENDDO
C
C           Entropie
            COEFA(IFAC,ICLT) =
     &           RTP(IEL,IPRIPH)/RTP(IEL,IRHIPH)**GAMAGR(IEL)
C
C --- Sortie subsonique
C     -----------------
          ELSEIF(XMACH.LT.1.D0 .AND. XMACH.GE.0.D0) THEN
C
C --- Detente :
            IF(COEFA(IFAC,ICLP).LE.RTP(IEL,IPRIPH)) THEN
C
C       Masse volumique
              COEFA(IFAC,ICLR) = RTP(IEL,IRHIPH)
     &             * (COEFA(IFAC,ICLP)/RTP(IEL,IPRIPH))
     &                **(1.D0/GAMAGR(IEL))
C
C       Vitesse
              COEFA(IFAC,ICLU) = RTP(IEL,IUIPH)
     &             + 2.D0/(GAMAGR(IEL)-1.D0)
     & * SQRT( GAMAGR(IEL) * RTP(IEL,IPRIPH) / RTP(IEL,IRHIPH) )
     &             * ( 1.D0
     & - (COEFA(IFAC,ICLP)/RTP(IEL,IPRIPH))
     &               **((GAMAGR(IEL)-1.D0)/2.D0/GAMAGR(IEL)) )
     & * SURFBO(1,IFAC) / RA(ISRFBN+IFAC-1)
C
              COEFA(IFAC,ICLV) = RTP(IEL,IVIPH)
     &             + 2.D0/(GAMAGR(IEL)-1.D0)
     & * SQRT( GAMAGR(IEL) * RTP(IEL,IPRIPH) / RTP(IEL,IRHIPH) )
     &             * ( 1.D0
     & - (COEFA(IFAC,ICLP)/RTP(IEL,IPRIPH))
     &               **((GAMAGR(IEL)-1.D0)/2.D0/GAMAGR(IEL)) )
     & * SURFBO(2,IFAC) / RA(ISRFBN+IFAC-1)
C
              COEFA(IFAC,ICLW) = RTP(IEL,IWIPH)
     &             + 2.D0/(GAMAGR(IEL)-1.D0)
     & * SQRT( GAMAGR(IEL) * RTP(IEL,IPRIPH) / RTP(IEL,IRHIPH) )
     &             * ( 1.D0
     & - (COEFA(IFAC,ICLP)/RTP(IEL,IPRIPH))
     &               **((GAMAGR(IEL)-1.D0)/2.D0/GAMAGR(IEL)) )
     & * SURFBO(3,IFAC) / RA(ISRFBN+IFAC-1)
C
C       Energie totale
              COEFA(IFAC,ICLE) = COEFA(IFAC,ICLP)
     & /( (GAMAGR(IEL)-1.D0)*COEFA(IFAC,ICLR) )
     &              + 0.5D0*(COEFA(IFAC,ICLU)**2
     &                     + COEFA(IFAC,ICLV)**2
     &                     + COEFA(IFAC,ICLW)**2)
C
C       Entropie
              COEFA(IFAC,ICLT) = COEFA(IFAC,ICLP)
     &                             /COEFA(IFAC,ICLR)**GAMAGR(IEL)
C
C --- Choc :
            ELSE
C
C       Masse volumique
              COEFA(IFAC,ICLR) = RTP(IEL,IRHIPH)
     & * ( (GAMAGR(IEL)+1.D0)*COEFA(IFAC,ICLP)
     &   + (GAMAGR(IEL)-1.D0)*RTP(IEL,IPRIPH) )
     & / ( (GAMAGR(IEL)-1.D0)*COEFA(IFAC,ICLP)
     &   + (GAMAGR(IEL)+1.D0)*RTP(IEL,IPRIPH) )
C
C       Vitesse
              COEFA(IFAC,ICLU) = RTP(IEL,IUIPH)
     & - (COEFA(IFAC,ICLP)-RTP(IEL,IPRIPH))*SQRT(2.D0/RTP(IEL,IRHIPH)
     & / ( (GAMAGR(IEL)+1.D0)*COEFA(IFAC,ICLP)
     &   + (GAMAGR(IEL)-1.D0)*RTP(IEL,IPRIPH) ))
     & * SURFBO(1,IFAC) / RA(ISRFBN+IFAC-1)
C
              COEFA(IFAC,ICLV) = RTP(IEL,IVIPH)
     & - (COEFA(IFAC,ICLP)-RTP(IEL,IPRIPH))*SQRT(2.D0/RTP(IEL,IRHIPH)
     & / ( (GAMAGR(IEL)+1.D0)*COEFA(IFAC,ICLP)
     &   + (GAMAGR(IEL)-1.D0)*RTP(IEL,IPRIPH) ))
     & * SURFBO(2,IFAC) / RA(ISRFBN+IFAC-1)
C
              COEFA(IFAC,ICLW) = RTP(IEL,IWIPH)
     & - (COEFA(IFAC,ICLP)-RTP(IEL,IPRIPH))*SQRT(2.D0/RTP(IEL,IRHIPH)
     & / ( (GAMAGR(IEL)+1.D0)*COEFA(IFAC,ICLP)
     &   + (GAMAGR(IEL)-1.D0)*RTP(IEL,IPRIPH) ))
     & * SURFBO(3,IFAC) / RA(ISRFBN+IFAC-1)
C
C       Energie totale
              COEFA(IFAC,ICLE) = COEFA(IFAC,ICLP)
     & /( (GAMAGR(IEL)-1.D0)*COEFA(IFAC,ICLR) )
     &     + 0.5D0*(COEFA(IFAC,ICLU)**2
     &            + COEFA(IFAC,ICLV)**2 + COEFA(IFAC,ICLW)**2)
C
C       Entropie
              COEFA(IFAC,ICLT) = COEFA(IFAC,ICLP)
     &                             /COEFA(IFAC,ICLR)**GAMAGR(IEL)
C
            ENDIF
C
          ELSE
            WRITE(NFECRA,*) 'ICCFTH = ',ICCFTH,'  MACH = ',XMACH
            IERR = 1
          ENDIF
C
          IF(IERR.EQ.1) CALL CSEXIT (1)
C
C
C     T ET e EN FONCTION DE P ET rho
C     ------------------------------
C
        ELSEIF(ICCFTH.EQ.912.OR.ICCFTH.EQ.60900) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C     Temperature
          COEFA(IFAC,ICLT) = XMASMR(IEL)/RR*COEFA(IFAC,ICLP)
     &                                        /COEFA(IFAC,ICLR)
C
C     Energie totale
          COEFA(IFAC,ICLE) = PROPCE(IEL,IPPROC(ICV(IPHAS)))
     &               * COEFA(IFAC,ICLT) + 0.5D0*( COEFA(IFAC,ICLU)**2
     &                    + COEFA(IFAC,ICLV)**2 + COEFA(IFAC,ICLW)**2)
C
C
C     rho ET e EN FONCTION DE P ET T
C     ------------------------------
C
        ELSEIF(ICCFTH.EQ.913.OR.ICCFTH.EQ.100900) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C     Masse volumique
          COEFA(IFAC,ICLR) = XMASMR(IEL)/RR*COEFA(IFAC,ICLP)
     &                                       /COEFA(IFAC,ICLT)
C
C     Energie totale
          COEFA(IFAC,ICLE) = PROPCE(IEL,IPPROC(ICV(IPHAS)))
     &               * COEFA(IFAC,ICLT) + 0.5D0*( COEFA(IFAC,ICLU)**2
     &                    + COEFA(IFAC,ICLV)**2 + COEFA(IFAC,ICLW)**2)
C
C
C     rho ET T EN FONCTION DE P ET e
C     ------------------------------
C
        ELSEIF(ICCFTH.EQ.914.OR.ICCFTH.EQ.140900) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C     Masse volumique
          COEFA(IFAC,ICLR) = COEFA(IFAC,ICLP)/(GAMAGR(IEL)-1.D0)
     &           / (COEFA(IFAC,ICLE) - 0.5D0*( COEFA(IFAC,ICLU)**2
     &                 + COEFA(IFAC,ICLV)**2 + COEFA(IFAC,ICLW)**2 ))
C
C     Temperature
          COEFA(IFAC,ICLT)= XMASMR(IEL)/RR*COEFA(IFAC,ICLP)
     &                                       /COEFA(IFAC,ICLR)
C
C
C     P ET e EN FONCTION DE rho ET T
C     ------------------------------
C
        ELSEIF(ICCFTH.EQ.923.OR.ICCFTH.EQ.150900) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C     Pression
          COEFA(IFAC,ICLP) = COEFA(IFAC,ICLR)*RR/XMASMR(IEL)
     &                                       *COEFA(IFAC,ICLT)
C
C     Energie totale
          COEFA(IFAC,ICLE) = PROPCE(IEL,IPPROC(ICV(IPHAS)))
     &               * COEFA(IFAC,ICLT) + 0.5D0*( COEFA(IFAC,ICLU)**2
     &                    + COEFA(IFAC,ICLV)**2 + COEFA(IFAC,ICLW)**2)
C
C
C     P ET T EN FONCTION DE rho ET e
C     ------------------------------
C
        ELSEIF(ICCFTH.EQ.924.OR.ICCFTH.EQ.210900) THEN
C
          IFAC = IFAC0
          IEL  = IFABOR(IFAC)
C
C     Pression
          COEFA(IFAC,ICLP) = (GAMAGR(IEL)-1.D0)*COEFA(IFAC,ICLR)
     &          *( COEFA(IFAC,ICLE) - 0.5D0*( COEFA(IFAC,ICLU)**2
     &                + COEFA(IFAC,ICLV)**2 + COEFA(IFAC,ICLW)**2 ) )

C
C     Temperature
          COEFA(IFAC,ICLT)= XMASMR(IEL)/RR*COEFA(IFAC,ICLP)
     &                                       /COEFA(IFAC,ICLR)
C
C
C --- Fin de test sur GAZ PARFAIT A GAMMA VARIABLE
        ENDIF
C
C --- Fin de test sur les thermos
      ENDIF
C
C
C--------
C FORMATS
C--------
C
 1010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : PROBLEME DANS USCFTH                        ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    GAMMA DOIT ETRE UN REEL SUPERIEUR OU EGAL A 1           ',/,
     &'@    IL A POUR VALEUR ',E12.4                                 ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : PROBLEME DANS USCFTH                        ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    GAMMA DOIT ETRE UN REEL SUPERIEUR OU EGAL A 1           ',/,
     &'@    IL EST NEGATIF OU NUL DANS LA CELLULE ',I10              ,/,
     &'@    IL A POUR VALEUR ',E12.4                                 ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : PROBLEME DANS USCFTH                        ',/,
     &'@    *********   CALCUL DE LA MASSE VOLUMIQUE                ',/,
     &'@                                                            ',/,
     &'@    LA TEMPERATURE DOIT ETRE UN REEL POSITIF STRICTEMENT    ',/,
     &'@    ELLE VAUT ',E12.4   ,' DANS LA CELLULE ',I10             ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : PROBLEME DANS USCFTH                        ',/,
     &'@    *********   CALCUL DE LA TEMPERATURE                    ',/,
     &'@                                                            ',/,
     &'@    LA MASSE VOLUMIQUE DOIT ETRE UN REEL POSITIF STRICTEMENT',/,
     &'@    ELLE VAUT ',E12.4   ,' DANS LA CELLULE ',I10             ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 4010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : PROBLEME DANS USCFTH                        ',/,
     &'@    *********   CALCUL DE LA VARIABLE THERMODYNAMIQUE ''C2''',/,
     &'@                                                            ',/,
     &'@    LA MASSE VOLUMIQUE DOIT ETRE UN REEL POSITIF STRICTEMENT',/,
     &'@    ELLE VAUT ',E12.4   ,' DANS LA CELLULE ',I10             ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : PROBLEME DANS USCFTH                        ',/,
     &'@    *********   CALCUL DE LA VARIABLE THERMO        ''BETA''',/,
     &'@                                                            ',/,
     &'@    LA MASSE VOLUMIQUE DOIT ETRE UN REEL POSITIF            ',/,
     &'@    ELLE VAUT ',E12.4   ,' DANS LA CELLULE ',I10             ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4030 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : PROBLEME DANS USCFTH                        ',/,
     &'@    *********   CALCUL DE LA VARIABLE THERMODYNAMIQUE ''S'' ',/,
     &'@                                                            ',/,
     &'@    LA MASSE VOLUMIQUE DOIT ETRE UN REEL POSITIF STRICTEMENT',/,
     &'@    ELLE VAUT ',E12.4   ,' DANS LA CELLULE ',I10             ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 7000 FORMAT (
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : CONDITION A LA LIMITE NON DISPONIBLE USCFTH ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    La condition a la limite de type debit et debit         ',/,
     &'@      enthalpique imposes n''est pas disponible dans la     ',/,
     &'@      version courante.                                     ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne peut pas etre execute.                     ',/,
     &'@                                                            ',/,
     &'@    Modifier uscfcl.                                        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8000 FORMAT (
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : MASSE VOLUMIQUE NEGATIVE OU NULLE (USCFTH)  ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    rho est negatif ou nul en ',I10   ,' cellules           ',/,
     &'@    On le limite a ',E12.4                                   ,/,
     &'@    et on met fin au calcul.                                ',/,
     &'@                                                            ',/,
     &'@    Si on souhaite que le calcul se poursuive, on peut      ',/,
     &'@      forcer un clipping standard par valeur inferieure en  ',/,
     &'@      renseignant la variable SCAMIN associee au scalaire   ',/,
     &'@      representant la masse volumique : ISCA(IRHO(IPHAS))   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 8100 FORMAT (
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ENERGIE INTERNE NEGATIVE OU NULLE (USCFTH)  ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    e-0.5*u*u est negatif ou nul en ',I10   ,' cellules     ',/,
     &'@    On le limite a ',E12.4                                   ,/,
     &'@    et on met fin au calcul.                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
C
C formats a eliminer apres mise au point de gama variable
C
C

 1110 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : PROBLEME DANS USTHMO                        ',/,
     &'@    *********   CALCUL DE LA VARIABLE THERMODYNAMIQUE ''C2''',/,
     &'@                                                            ',/,
     &'@    P DOIT ETRE UN REEL POSITIF                             ',/,
     &'@    IL EST NEGATIF DANS LA CELLULE ',I10                     ,/,
     &'@    IL A POUR VALEUR ',E12.4                                 ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1120 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : PROBLEME DANS USTHMO                        ',/,
     &'@    *********   CALCUL DE LA VARIABLE THERMODYNAMIQUE ''C2''',/,
     &'@                                                            ',/,
     &'@    RHO DOIT ETRE UN REEL POSITIF STRICTEMENT               ',/,
     &'@    IL EST NEGATIF OU NUL DANS LA CELLULE ',I10              ,/,
     &'@    IL A POUR VALEUR ',E12.4                                 ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1220 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : PROBLEME DANS USTHMO                        ',/,
     &'@    *********   CALCUL DE LA VARIABLE THERMO ''BETA''       ',/,
     &'@                                                            ',/,
     &'@    RHO DOIT ETRE UN REEL POSITIF                           ',/,
     &'@    IL EST NEGATIF DANS LA CELLULE ',I10                     ,/,
     &'@    IL A POUR VALEUR ',E12.4                                 ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1310 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : PROBLEME DANS USTHMO                        ',/,
     &'@    *********   CALCUL DE LA VARIABLE THERMODYNAMIQUE ''S'' ',/,
     &'@                                                            ',/,
     &'@    P DOIT ETRE UN REEL POSITIF                             ',/,
     &'@    IL EST NEGATIF DANS LA CELLULE ',I10                     ,/,
     &'@    IL A POUR VALEUR ',E12.4                                 ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1320 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : PROBLEME DANS USTHMO                        ',/,
     &'@    *********   CALCUL DE LA VARIABLE THERMODYNAMIQUE ''S'' ',/,
     &'@                                                            ',/,
     &'@    RHO DOIT ETRE UN REEL POSITIF STRICTEMENT               ',/,
     &'@    IL EST NEGATIF OU NUL DANS LA CELLULE ',I10              ,/,
     &'@    IL A POUR VALEUR ',E12.4                                 ,/,
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
