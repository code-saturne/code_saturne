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
                        SUBROUTINE VORPRE
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IREPVO ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PROPCE , PROPFA , PROPFB ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC     ROUTINE DE PREPATATION DE LA METHODE DES VORTEX
CFONC     Gestion memoire, connectivites, ...
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! IFINIA       ! E  ! <-  ! NUMERO DE LA DERNIERE CASE LIBRE DANS!
CARGU !              !    !     ! IA APRES LA RESERVATION MEMOIRE      !
CARGU ! IFINRA       ! E  ! <-  ! NUMERO DE LA DERNIERE CASE LIBRE DANS!
CARGU !              !    !     ! RA APRES LA RESERVATION MEMOIRE      !
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
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! IREPVO       ! TE ! <-> ! TAB ENTIER POUR REPERAGE DES FACES DE!
CARGU !              !    !     ! BORD POUR LA METHODE DES VORTEX      !
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
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "vortex.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0 , IFINIA , IFINRA
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
      INTEGER          IREPVO(NFABOR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC, IEL, II, IPHAS
      INTEGER          IENT, IPCVIS, IPCROM
      INTEGER          IAPPEL
      INTEGER          ISURF(NENTMX)
      DOUBLE PRECISION XX, YY, ZZ
      DOUBLE PRECISION XXV, YYV, ZZV
C
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      NVOMAX = 0
      DO IENT = 1, NNENT
        NVOMAX = MAX(NVORT(IENT),NVOMAX)
      ENDDO
C
C NVOMAX = nombre max de vortex (utilise plus tard)
C
      DO IENT = 1, NNENT
        ICVOR2(IENT) = 0
      ENDDO
C
      DO IFAC = 1, NFABOR
        IENT = IREPVO(IFAC)
        IF(IENT.NE.0) THEN
          ICVOR2(IENT) = ICVOR2(IENT) + 1
        ENDIF
      ENDDO
C
C ICVOR2 = compteur du nombre local de faces
C   utilisant des vortex a l'entree IENT
C
      ICVMAX = 0
      IF(IRANGP.GE.0) THEN
        DO IENT = 1, NNENT
          ICVOR(IENT) = ICVOR2(IENT)
          CALL PARCPT(ICVOR(IENT))
C         ===========
          ICVMAX = MAX(ICVMAX,ICVOR(IENT))
        ENDDO
      ELSE
        DO IENT = 1, NNENT
          ICVOR(IENT) = ICVOR2(IENT)
          ICVMAX = MAX(ICVMAX,ICVOR(IENT))
        ENDDO
      ENDIF
C
C ICVOR = nombre global de faces utilisant des vortex a l'entree IENT
C
C ICVMAX = max du nombre global de faces utilisant des vortex
C (toutes entrees confondues).
C
      IAPPEL = 2
      CALL MEMVOR
C     ===========
     & ( IDEBIA , IDEBRA , IAPPEL , NFABOR , IFINIA , IFINRA )
C
      IDEBIA = IFINIA
      IDEBRA = IFINRA
C
C=======================================================================
C 2. CONSTRUCTION DE LA " GEOMETRIE GOBALE "
C=======================================================================
C
      DO IENT = 1, NNENT
        ICVOR2(IENT) = 0
        XSURFV(IENT) = 0.D0
        ISURF(IENT)  = 0
      ENDDO
C
C Chaque processeur stocke dans les tableaux RA(IW1X),...
C les coordonnees des faces ou il doit ensuite utiliser des vortex
C
      IPHAS  = 1
      IPCVIS = IPPROC(IVISCL(IPHAS))
      IPCROM = IPPROC(IROM(IPHAS))
      DO IFAC = 1, NFABOR
        IENT = IREPVO(IFAC)
        IF(IENT.NE.0) THEN
          IEL = IFABOR(IFAC)
          ICVOR2(IENT) = ICVOR2(IENT) + 1
          RA(IW1X+(IENT-1)*ICVMAX+ICVOR2(IENT)-1)= CDGFBO(1,IFAC)
          RA(IW1Y+(IENT-1)*ICVMAX+ICVOR2(IENT)-1)= CDGFBO(2,IFAC)
          RA(IW1Z+(IENT-1)*ICVMAX+ICVOR2(IENT)-1)= CDGFBO(3,IFAC)
          RA(IW1V+(IENT-1)*ICVMAX+ICVOR2(IENT)-1) =
     &      PROPCE(IEL,IPCVIS)/PROPCE(IEL,IPCROM)
          XSURFV(IENT) = XSURFV(IENT) + SQRT(SURFBO(1,IFAC)**2
     &      + SURFBO(2,IFAC)**2 + SURFBO(3,IFAC)**2)
C         Vecteur surface d'une face de l'entree
          IF (ISURF(IENT).EQ.0) THEN
            SURF(1,IENT) = SURFBO(1,IFAC)
            SURF(2,IENT) = SURFBO(2,IFAC)
            SURF(3,IENT) = SURFBO(3,IFAC)
            ISURF(IENT)  = 1
          ENDIF
        ENDIF
      ENDDO
C
      IF(IRANGP.GE.0) THEN
        DO IENT = 1, NNENT
          CALL PARSOM(XSURFV(IENT))
C         ===========
        ENDDO
      ENDIF
C
C -------------
C En parallele
C -------------
      IF(IRANGP.GE.0) THEN
        DO IENT = 1, NNENT
          CALL PARAGV
C         ===========
     & ( ICVOR2(IENT), ICVOR(IENT),
     &   RA(IW1X  + (IENT-1)*ICVMAX)   ,
     &   RA(IXYZV + (IENT-1)*3*ICVMAX) )
          CALL PARAGV
C         ===========
     & ( ICVOR2(IENT), ICVOR(IENT),
     &   RA(IW1Y  + (IENT-1)*ICVMAX)              ,
     &   RA(IXYZV + (IENT-1)*3*ICVMAX +   ICVMAX) )
          CALL PARAGV
C         ===========
     & ( ICVOR2(IENT), ICVOR(IENT),
     &   RA(IW1Z  + (IENT-1)*ICVMAX)              ,
     &   RA(IXYZV + (IENT-1)*3*ICVMAX + 2*ICVMAX) )
          CALL PARAGV
C         ===========
     & ( ICVOR2(IENT), ICVOR(IENT),
     &   RA(IW1V  + (IENT-1)*ICVMAX)   ,
     &   RA(IVISV + (IENT-1)*3*ICVMAX) )
        ENDDO
C
C  -> A la fin de cette etape, tous les processeurs connaissent
C     les coordonees des faces d'entree
C
      ELSE
C ----------------------
C Sur 1 seul processeur
C ----------------------
        DO IENT = 1,NNENT
          DO II = 1, ICVOR(IENT)
            RA(IXYZV+(IENT-1)*3*ICVMAX+II-1)=
     &           RA(IW1X+(IENT-1)*ICVMAX + II -1)
            RA(IXYZV+(IENT-1)*3*ICVMAX+ICVMAX + II-1) =
     &           RA(IW1Y+(IENT-1)*ICVMAX + II -1)
            RA(IXYZV+(IENT-1)*3*ICVMAX+2*ICVMAX + II-1) =
     &           RA(IW1Z+(IENT-1)*ICVMAX + II -1)
            RA(IVISV+(IENT-1)*ICVMAX+II-1) =
     &           RA(IW1V+(IENT-1)*ICVMAX+II-1)
          ENDDO
        ENDDO
      ENDIF
C
C=======================================================================
C 3. CONSTRUCTION DE LA CONNECTIVITE
C=======================================================================
C
      DO IENT = 1, NNENT
        ICVOR2(IENT) = 0
        DO IFAC = 1, ICVMAX
          IA(IIFAGL+(IENT-1)*ICVMAX+IFAC-1) = 0
        ENDDO
      ENDDO
C
C On cherche ensuite le numero de la ligne du tableau RA(IXYZV) qui est
C associe a la Ieme face d'entree utilisant des vortex (dans la
C numerotation chronologique que suit ICVOR2).
C
      DO IFAC = 1, NFABOR
        IENT = IREPVO(IFAC)
        IF(IENT.NE.0) THEN
          ICVOR2(IENT) = ICVOR2(IENT) + 1
          DO II = 1, ICVOR(IENT)
            XX = CDGFBO(1,IFAC)
            YY = CDGFBO(2,IFAC)
            ZZ = CDGFBO(3,IFAC)
            XXV = RA(IXYZV+(IENT-1)*3*ICVMAX+II-1)
            YYV = RA(IXYZV+(IENT-1)*3*ICVMAX+ICVMAX+II-1)
            ZZV = RA(IXYZV+(IENT-1)*3*ICVMAX+2*ICVMAX+II-1)
            IF(ABS(XXV-XX).LT.EPZERO.AND.ABS(YYV-YY).LT.EPZERO.AND.
     &           ABS(ZZV-ZZ).LT.EPZERO) THEN
              IA(IIFAGL+(IENT-1)*ICVMAX+ICVOR2(IENT)-1) = II
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
C La methode de vortex va generer un tableau de vitesse RA(IUVOR)
C qui aura la meme structure que RA(IXYZV).
C Le tableau RA(IXYZV) sera envoyee a tous les processeurs
C la vitesse a imposer à la Ieme face se trouvera à la ligne IA(IIFAGL+I)
C
      IAPPEL = 3
      CALL MEMVOR
C     ===========
     & ( IDEBIA , IDEBRA , IAPPEL , NFABOR , IFINIA , IFINRA )
C
C ---
C FIN
C ---
C
      RETURN
      END
c@z
