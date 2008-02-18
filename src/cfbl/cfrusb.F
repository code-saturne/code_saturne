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
                        SUBROUTINE CFRUSB
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IMODIF , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   GAMMAG ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   SORTI1 , SORTI2 , GAMAGR , MASMOR ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ---------
c@foncb
CFONC
CFONC Flux de rusanov au bord pour euler + energie
CFONC
CFONC  d rho   /dt + div rho u             = 0
CFONC  d rho u /dt + div rho u u + grad  P = 0
CFONC  d E     /dt + div rho u E + div u P = 0
CFONC
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
CARGU ! IMODIF       ! E  !  -> ! MODIFICATION DIRECTE DE RTP (IMODIF=1!
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
CARGU ! GAMMAG       ! R  !  -> ! GAMMA DU GAZ                         !
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
CARGU ! GAMAGR(*)    ! TR !  -  ! CONSTANTE GAMMA EQUIVALENT DU GAZ    !
CARGU ! MASMOR(*)    ! TR !  -  ! MASSE MOLAIRE DES CONSTITUANTS DU GAZ!
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
      INTEGER          IMODIF , IPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION GAMMAG
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*),PROPFA(NFAC,*),PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION SORTI1(*), SORTI2(*), GAMAGR(*), MASMOR(*)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IFAC0
      INTEGER          IEL    , IFAC
      INTEGER          IPRIPH , IRHIPH , IENIPH
      INTEGER          IUIPH  , IVIPH  , IWIPH
      INTEGER          ICLP   , ICLR   , ICLE
      INTEGER          ICLU   , ICLV   , ICLW
      INTEGER          IFLMAB
      DOUBLE PRECISION UND    , UNI    , RUND   , RUNI   , CD     , CI
      DOUBLE PRECISION RRUS   , RUNB

C
C***********************************************************************
C
C=======================================================================
C 0. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IPRIPH = IPR(IPHAS)
      IRHIPH = ISCA(IRHO  (IPHAS))
      IENIPH = ISCA(IENERG(IPHAS))
      IUIPH  = IU(IPHAS)
      IVIPH  = IV(IPHAS)
      IWIPH  = IW(IPHAS)
      ICLP = ICLRTP(IPRIPH,ICOEF)
      ICLR = ICLRTP(IRHIPH,ICOEF)
      ICLE = ICLRTP(IENIPH,ICOEF)
      ICLU = ICLRTP(IUIPH,ICOEF)
      ICLV = ICLRTP(IVIPH,ICOEF)
      ICLW = ICLRTP(IWIPH,ICOEF)
C
      IFLMAB = IPPROB(IFLUMA(IENIPH))
C
      IFAC0 = IMODIF
      IFAC  = IFAC0
      IEL   = IFABOR(IFAC)
C
C=======================================================================
C 1. GRANDEURS LIEES A RUSANOV
C=======================================================================
C
C
      UND   = (COEFA(IFAC,ICLU)*SURFBO(1,IFAC)
     &       + COEFA(IFAC,ICLV)*SURFBO(2,IFAC)
     &       + COEFA(IFAC,ICLW)*SURFBO(3,IFAC))/RA(ISRFBN+IFAC-1)
      UNI   = (RTP(IEL,IUIPH)*SURFBO(1,IFAC)
     &       + RTP(IEL,IVIPH)*SURFBO(2,IFAC)
     &       + RTP(IEL,IWIPH)*SURFBO(3,IFAC))/RA(ISRFBN+IFAC-1)
      RUND  = COEFA(IFAC,ICLR)*UND
      RUNI  = RTP(IEL,IRHIPH)     *UNI
      CD    = SQRT(GAMMAG*COEFA(IFAC,ICLP)/COEFA(IFAC,ICLR))
      CI    = SQRT(GAMMAG*RTP(IEL,IPRIPH)/RTP(IEL,IRHIPH))
      RRUS  = MAX(ABS(UND)+CD,ABS(UNI)+CI)
C
      RUNB  = 0.5D0*(COEFA(IFAC,ICLR)*UND+RTP(IEL,IRHIPH)*UNI)
     &      - 0.5D0*RRUS*(COEFA(IFAC,ICLR)-RTP(IEL,IRHIPH))
C
C=======================================================================
C 2. FLUX CONVECTIFS DE RUSANOV
C=======================================================================
C
C
C     Reperage de la face pour annuler les flux convectifs
C       calcules au bord par bilsc2 ou cfbsc2 pour la qdm (div(rho u u))
C       et l'energie (div(rho u E)) ainsi que les termes en
C       grad(P) et div(u P)
C
      IA(IIFBRU+IFAC-1+(IPHAS-1)*NFABOR) = 1
C
C     Flux de masse
      PROPFB(IFAC,IFLMAB) = RUNB*RA(ISRFBN+IFAC-1)
C
C     Flux de Qdm (la partie centree en pression pourrait etre prise dans
C       la condition à la limite de pression, ce qui eviterait de retoucher
C       le gradient de pression de la qdm, mais qui donne moins de
C       flexibilité quant à la condition à la limite de pression utilisee
C       pour la reconstruction du gradient, si le maillage est non
C       orthogonal en entree)
      PROPFB(IFAC,IPPROB(IFBRHU(IPHAS))) = RA(ISRFBN+IFAC-1)*
     &             0.5D0*(
     &             (RUND*COEFA(IFAC,ICLU)+RUNI*RTP(IEL,IUIPH))
     &             -RRUS*( COEFA(IFAC,ICLR)*COEFA(IFAC,ICLU)
     &                    -RTP(IEL,IRHIPH)     *RTP(IEL,IUIPH)     ))
     &                        + SURFBO(1,IFAC)*
     &             0.5D0*(COEFA(IFAC,ICLP)+RTP(IEL,IPRIPH))
      PROPFB(IFAC,IPPROB(IFBRHV(IPHAS))) = RA(ISRFBN+IFAC-1)*
     &             0.5D0*(
     &             (RUND*COEFA(IFAC,ICLV)+RUNI*RTP(IEL,IVIPH))
     &             -RRUS*( COEFA(IFAC,ICLR)*COEFA(IFAC,ICLV)
     &                    -RTP(IEL,IRHIPH)     *RTP(IEL,IVIPH)     ))
     &                        + SURFBO(2,IFAC)*
     &             0.5D0*(COEFA(IFAC,ICLP)+RTP(IEL,IPRIPH))
      PROPFB(IFAC,IPPROB(IFBRHW(IPHAS))) = RA(ISRFBN+IFAC-1)*
     &             0.5D0*(
     &             (RUND*COEFA(IFAC,ICLW)+RUNI*RTP(IEL,IWIPH))
     &             -RRUS*( COEFA(IFAC,ICLR)*COEFA(IFAC,ICLW)
     &                    -RTP(IEL,IRHIPH)     *RTP(IEL,IWIPH)     ))
     &                        + SURFBO(3,IFAC)*
     &             0.5D0*(COEFA(IFAC,ICLP)+RTP(IEL,IPRIPH))
C     Flux de E
      PROPFB(IFAC,IPPROB(IFBENE(IPHAS))) = RA(ISRFBN+IFAC-1)*
     &             0.5D0*(
     &              RUND*COEFA(IFAC,ICLE)+RUNI*RTP(IEL,IENIPH)
     &              +UND*COEFA(IFAC,ICLP)+ UNI*RTP(IEL,IPRIPH)
     &             -RRUS*( COEFA(IFAC,ICLR)*COEFA(IFAC,ICLE)
     &                    -RTP(IEL,IRHIPH)     *RTP(IEL,IENIPH)     ))
C
C
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
