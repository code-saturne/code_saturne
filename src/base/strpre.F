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
C     Copyright (C) 1998-2007 EDF S.A., France
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
                        SUBROUTINE STRPRE
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , ITRALE , ITALIM , INEEFL ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IMPALE ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   FLMALF , FLMALB , XPRALE , COFALE , DEPALE , RDEVEL , RTUSER ,
     &   RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  PREDICTION DU DEPLACEMENT DES STRUCTURES EN ALE
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
CARGU ! ITRALE       ! E  !  -> ! NUMERO D'ITERATION POUR L'ALE        !
CARGU ! ITALIM       ! E  !  -> ! NUMERO D'ITERATION COUPLAGE IMPLICITE!
CARGU ! INEEDF       ! E  !  -> ! INDICATEUR DE SAUVEGARDE DES FLUX    !
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
CARGU !   (NFAC+1)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (LNDFAC)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (NFABOR+1) !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (LNDFBR)   !    !     !  (OPTIONNEL)                         !
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
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR ! <-> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! FLMALF(NFAC) ! TR ! <-  ! SAUVEGARDE DU FLUX DE MASSE FACES INT!
CARGU ! FLMALB(NFABOR! TR ! <-  ! SAUVEGARDE DU FLUX DE MASSE FACES BRD!
CARGU ! COFALE       ! TR ! <-  ! SAUVEGARDE DES CL DE P ET U          !
CARGU !    (NFABOR,8)!    !     !                                      !
CARGU ! XPRALE(NCELET! TR ! <-  ! SAUVEGARDE DE LA PRESSION, SI NTERUP !
CARGU !              !    !     !    EST >1                            !
CARGU ! DEPALE(NNOD,3! TR !  -> ! DEPLACEMENT AUX NOEUDS               !
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
C
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
      INCLUDE "optcal.h"
      INCLUDE "numvar.h"
      INCLUDE "pointe.h"
      INCLUDE "albase.h"
      INCLUDE "alstru.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "entsor.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0 , ITRALE , ITALIM , INEEFL
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IMPALE(NNOD)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION FLMALF(NFAC), FLMALB(NFABOR), XPRALE(NCELET)
      DOUBLE PRECISION COFALE(NFABOR,8)
      DOUBLE PRECISION DEPALE(NNOD,3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          ISTR, II, IFAC, INOD, IEL
      INTEGER          IFLMAS, IFLMAB,ICLP,ICLU,ICLV,ICLW
C
C***********************************************************************
C
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IFLMAS = IPPROF(IFLUMA(IU(1)))
      IFLMAB = IPPROB(IFLUMA(IU(1)))
      ICLP = ICLRTP(IPR(1),ICOEF)
      ICLU = ICLRTP(IU(1),ICOEF)
      ICLV = ICLRTP(IV(1),ICOEF)
      ICLW = ICLRTP(IW(1),ICOEF)
C
C=======================================================================
C 2. PREDICTION DU DEPLACEMENT
C=======================================================================
C
C Lors de la phase d'initialisation de l'ALE (ITRALE=0), XSTP contient
C    - la valeur du deplacement initial des structures si l'utilisateur
C        les a touchees (suite ou pas)
C    - 0 si debut de calcul avec les structures
C    - le deplacement utilise pour le calcul precedent si suite sans
C        modification utilisateur
C   Il faut cependant transferer sa valeur dans XSTR (qui est utilise
C   par Newmark)
C Lors des iterations suivantes (ITRALE>0) on utilise les schemas standard
C   de calcul de XSTP
C
      IF (ITRALE.EQ.0) THEN
C
        DO ISTR = 1, NBSTRU
          DO II = 1, NDIM
            XSTR(II,ISTR) = XSTP(II,ISTR)
          ENDDO
        ENDDO
      ELSE
C
C 2.1 : SCHEMA DE COUPLAGE EXPLICITE
C-----------------------------------
        IF (NALIMX.EQ.1) THEN
          DO ISTR = 1, NBSTRU
            DO II = 1, 3
              XSTP(II,ISTR) = XSTR(II,ISTR)
     &             + AEXXST*DTSTR(ISTR)*XPSTR(II,ISTR)
     &             + BEXXST*DTSTR(ISTR)*(XPSTR(II,ISTR)-XPSTA(II,ISTR))
            ENDDO
          ENDDO
C
C 2.2 : SCHEMA DE COUPLAGE IMPLICITE
C-----------------------------------
        ELSE
          DO ISTR = 1, NBSTRU
            DO II = 1, 3
              XSTP(II,ISTR)  = XSTR(II,ISTR)
            ENDDO
          ENDDO
        ENDIF
C
      ENDIF
C
      DO IFAC = 1, NFABOR
        ISTR = IA(IIDFST+IFAC-1)
        IF (ISTR.GT.0) THEN
          DO II = IPNFBR(IFAC), IPNFBR(IFAC+1)-1
            INOD = NODFBR(II)
            IMPALE(INOD) = 1
            DEPALE(INOD,1) = XSTP(1,ISTR)
            DEPALE(INOD,2) = XSTP(2,ISTR)
            DEPALE(INOD,3) = XSTP(3,ISTR)
          ENDDO
        ENDIF
      ENDDO
C
C=======================================================================
C 3. DEPLACEMENT AU PAS DE TEMPS PRECEDENT ET SAUVEGARDE FLUX ET DE P
C=======================================================================
      IF (ITALIM.EQ.1) THEN
        DO ISTR = 1, NBSTRU
          DO II = 1, 3
            XSTA(II,ISTR)   = XSTR(II,ISTR)
            XPSTA(II,ISTR)  = XPSTR(II,ISTR)
            XPPSTA(II,ISTR) = XPPSTR(II,ISTR)
          ENDDO
        ENDDO
        IF (INEEFL.EQ.1) THEN
          DO IFAC = 1, NFAC
            FLMALF(IFAC) = PROPFA(IFAC,IFLMAS)
          ENDDO
          DO IFAC = 1, NFABOR
            FLMALB(IFAC) = PROPFB(IFAC,IFLMAB)
            COFALE(IFAC,1) = COEFA(IFAC,ICLP)
            COFALE(IFAC,2) = COEFA(IFAC,ICLU)
            COFALE(IFAC,3) = COEFA(IFAC,ICLV)
            COFALE(IFAC,4) = COEFA(IFAC,ICLW)
            COFALE(IFAC,5) = COEFB(IFAC,ICLP)
            COFALE(IFAC,6) = COEFB(IFAC,ICLU)
            COFALE(IFAC,7) = COEFB(IFAC,ICLV)
            COFALE(IFAC,8) = COEFB(IFAC,ICLW)
          ENDDO
          IF (NTERUP.GT.1) THEN
            DO IEL = 1, NCELET
              XPRALE(IEL) = RTPA(IEL,IPR(1))
            ENDDO
          ENDIF
        ENDIF
      ENDIF
C
C----
C FORMATS
C----
C

C
C----
C FIN
C----
C
      END
c@z
