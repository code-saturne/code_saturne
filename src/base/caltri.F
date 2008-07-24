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
                        SUBROUTINE CALTRI
C                       *****************
C     -------------------------------------------------------------
     & ( IVERIF ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD ,
     &   VOLUME ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ----------
c@foncb
CFONC
CFONC  GESTION DU PROGRAMME (LECTURE, RESOLUTION, ECRITURE)
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IVERIF       ! E  !  -> ! INDICATEUR DES TESTS ELEMENTAIRES    !
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
CARGU ! (NFML,NPRFML)!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (NFAC+1)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (LNDFAC)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (NFABOR+1) !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (LNDFBR)   !    !     !  (OPTIONNEL)                         !
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
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! RA(*)        ! TR !  -  ! MACRO TABLEAU REEL                   !
CARGU .______________.____._____.______________________________________.
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
      INCLUDE "dimens.h"
      INCLUDE "pointe.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "albase.h"
      INCLUDE "radiat.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagdim.h"
      INCLUDE "lagran.h"
      INCLUDE "vortex.h"
      INCLUDE "ihmpre.h"
      INCLUDE "matiss.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IVERIF
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

      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IPROPC , IPROPF , IPROPB
      INTEGER          ICOEFA , ICOEFB
      INTEGER          IRTP   , IRTPA
      INTEGER          IDT
      INTEGER          IISSTD , IFRCX
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IFINIA , IFINRA , IDBIA1 , IDBRA1, IDBIA2
      INTEGER          IDITIA, IDITRA
      INTEGER          IFNIA1 , IFNRA1 , IFNIA2 , IFNIA3, IFNRA2
      INTEGER          JCELBR
      INTEGER          IIII
C
      INTEGER          MODHIS, IAPPEL, MODNTL, IISUIT, IWARN0
      INTEGER          NTSDEF, NTHDEF, NTCREL
      INTEGER          IMSFIN, IPHAS , IVAR
C
      INTEGER          IICOCE , IITYCE
      INTEGER          IIITEP , IITEPA , ISTATC , ISTATF
      INTEGER          IETTP  , IETTPA , IAUXL  , ITSLAG , ISTATV
      INTEGER          ITAUP  , IITLAG , IPIIL  , IINDEP , IIBORD
      INTEGER          IVAGAU , ITSUF  , ITSUP  , IBX    , IAUXL2
      INTEGER          IBRGAU , ITEBRU
      INTEGER          IGRADP , IGRADV , ICROUL
      INTEGER          ITEPCT , ITSFEX , ITSVAR
      INTEGER          ICPGD1 , ICPGD2 , ICPGHT
      INTEGER          ILAGIA , ILAGRA , IIWORK
      INTEGER          IW1    , IW2    , IW3
      INTEGER          INOD   , IDIM
      INTEGER          ITRALE , INDACT , INDWRI
      INTEGER          MAXELT , ILS
C
      DOUBLE PRECISION TITER1, TITER2
      DOUBLE PRECISION TECRF1, TECRF2
C
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATIONS
C=======================================================================
C
C     Initialisation de la premiere case libre des tableaux IA et RA
      IDEBIA = 1
      IDEBRA = 1
C
C     Initialisation du generateur de nombres aleatoires
C     (pas toujours necessaire mais ne coute rien)
      CALL ZUFALLI(0)
C     ============
C
C---> MISE A ZERO DES TABLEAUX UTILISATEUR ET DEVELOPPEUR
C
      IF(NITUSE.GT.0) THEN
        DO IIII = 1, NITUSE
          ITUSER(IIII) = 0
        ENDDO
      ENDIF
      IF(NRTUSE.GT.0) THEN
        DO IIII = 1, NRTUSE
          RTUSER(IIII) = 0.D0
        ENDDO
      ENDIF
C
      IF(NIDEVE.GT.0) THEN
        DO IIII = 1, NIDEVE
          IDEVEL(IIII) = 0
        ENDDO
      ENDIF
      IF(NRDEVE.GT.0) THEN
        DO IIII = 1, NRDEVE
          RDEVEL(IIII) = 0.D0
        ENDDO
      ENDIF
C
C---> Test d'arret mis a 1 si le rayonnement P-1 voit trop de cellules
C       a epaisseur optique superieure a l'unite (voir ppcabs)
      ISTPP1 = 0
C
C---> Nombre max d'elements pour le selector
      MAXELT = MAX(NCELET,NFAC,NFABOR)
C
C=======================================================================
C 2. GEOMETRIE
C=======================================================================
C
C---> CALCULS GEOMETRIQUES
C
C     (MEMCLG remplit directement pointe.h)
      CALL MEMCLG
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFINIA , IFINRA )
C
C---> La memoire sera conservee jusqu'a la fin.
      IDEBIA = IFINIA
      IDEBRA = IFINRA
C
C
      CALL CREGEO
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   VOLTOT ,
     &   RDEVEL , RTUSER , RA     )
C
C
C
C---> CALCUL DE JCELBR (=NCELBR) et REMPLISSAGE DE IA(IICELB)
C         directement en passant par le pointeur et IA.
C         (c'est pas bien, mais c'est mieux qu'avant)
C
      IICELB = IDEBIA
      CALL MEMCBR
C     ===========
     & ( IICELB , NCELET , NCEL   , NFABOR ,
     &   JCELBR , IFINIA ,
     &   IFABOR ,
     &   IA     )
C
C---> La memoire sera conservee jusqu'a la fin.
      IDEBIA = IFINIA
      IDEBRA = IFINRA
C
C=======================================================================
C 3. FIN INITIALISATION DES COMMONS
C=======================================================================
C
      CALL INITI2
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   JCELBR ,
     &   IA     , RA     )
C
      IF (IILAGR.GT.0) THEN
C
C--> Calcul de LNDNOD (lagran.h)
C
C   Tableau NCELET de travail entier
        IIWORK = IDEBIA
        IFINIA = IIWORK + NCELET
        CALL IASIZE ('LAGINI',IFINIA)
C
        IFINRA = IDEBRA
C
        CALL LAGINI
C       ===========
     & ( IFINIA , IFINRA ,
     &   NCELET , NCEL   , NFAC   , NFABOR ,
     &   LNDNOD ,
     &   IFACEL , IFABOR ,
     &   IA(IIWORK) ,
     &   IA     , RA     )
C
      ENDIF
C
C
C=======================================================================
C 4. AUTRES TABLEAUX
C=======================================================================
C
C---> GESTION MEMOIRE
C
      CALL MEMTRI
C     ===========
     & ( IDEBIA , IDEBRA , IVERIF ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCOFAB , NPROCE , NPROFA , NPROFB ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IISSTD , IFRCX  ,
     &   IDT    , IRTP   , IRTPA  , IPROPC , IPROPF , IPROPB ,
     &   ICOEFA , ICOEFB ,
     &   IFINIA , IFINRA )
C
C     Reservations complementaires pour Matisse
C       On remplit un pointeur dans matiss.h, pour eviter de charger
C       les arguments par des tableaux specifiques a Matisse.
C
      IF (IMATIS.EQ.1) THEN
C
        IDBIA1 = IFINIA
        IDBRA1 = IFINRA
        CALL MEMMAT
C       ===========
     & ( IDBIA1 , IDBRA1 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCOFAB , NPROCE , NPROFA , NPROFB ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFINIA , IFINRA )
C
      ENDIF
C
C=======================================================================
C 4.1 RESERVATION DE LA MEMOIRE POUR LE RAYONNEMENT SEMI-TRANSPARENT
C=======================================================================
C
      IF (IIRAYO.GT.0) THEN
C
        IDBIA1 = IFINIA
        IDBRA1 = IFINRA
C
        CALL MEMRA1
C       ===========
     & ( IDBIA1 , IDBRA1 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IFINIA , IFINRA )
C
      ENDIF
C
C=======================================================================
C 4.2 RESERVATION DE LA MEMOIRE POUR LE LAGRANGIEN
C=======================================================================
C
C     Si on ne fait pas de Lagrangien, on initialise
C       quand meme les "pointeurs".
C
      IDBIA1 = IFINIA
      IDBRA1 = IFINRA
C
      CALL MEMLA1
C     ===========
     &  ( IDBIA1 , IDBRA1 ,
     &    NDIM   , NCELET , NCEL   , NFAC   , NFABOR ,
     &    LNDNOD ,
     &    NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &    NTERSL , NVLSTA , NVISBR ,
     &    IIITEP , IICOCE , IITYCE ,
     &    IETTP  , IETTPA , IITEPA , ISTATC , ISTATV, ITSLAG , ISTATF ,
     &    IFINIA , IFINRA )
C
C=======================================================================
C 4.3 TESTS ELEMENTAIRES : APPEL A TESTEL.F
C=======================================================================
C
      IF (IVERIF.GT.0) THEN
C
        WRITE(NFECRA, 1000) IVERIF
C
        CALL TESTEL
C       ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  , NVAR   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME , RA(IRTP) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RDEVEL , RTUSER , RA     )
C
        GOTO 200
C
      ENDIF
C
C=======================================================================
C 5. INITIALISATIONS PAR DEFAUT
C=======================================================================
C
      CALL INIVA0
C     ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCOFAB ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTP)   , RA(IPROPC) , RA(IPROPF) , RA(IPROPB),
     &   RA(ICOEFA) , RA(ICOEFB) , RA(IFRCX ) ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 6. CALCUL SUITE EVENTUEL
C=======================================================================
C
      IF (ISUITE.EQ.1) THEN
C
        CALL LECAMO
C       ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NNOD   ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDEVEL , ITUSER , IA     ,
     &   RA(IDT)    , RA(IRTP)   , RA(IPROPC) , RA(IPROPF) , RA(IPROPB),
     &   RA(ICOEFA) , RA(ICOEFB) , RA(IFRCX ) ,
     &   RDEVEL , RTUSER , RA     )
C
C     En ALE, il faut recalculer les parametres geometriques
        IF (IALE.EQ.1) THEN
C
          DO INOD = 1, NNOD
            DO IDIM = 1, NDIM
              XYZNOD(IDIM,INOD) = RA(IXYZN0+(INOD-1)*NDIM+IDIM-1)
     &             + RA(IDEPAL+(IDIM-1)*NNOD+INOD-1)
            ENDDO
          ENDDO
C
          CALL ALGRMA
C         ===========
          CALL CALGEO
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   VOLTOT ,
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF
C
      ENDIF
C
C=======================================================================
C 7. INITIALISATIONS (Utilisateur et complementaires)
C    RTP DT ROM ROMB VISCL VISCT VISCLS
C    (TPUCOU en PERIODICITE)
C=======================================================================
C
      CALL INIVAR
C     ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCOFAB ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTP)   , RA(IPROPC) , RA(IPROPF) , RA(IPROPB),
     &   RA(ICOEFA) , RA(ICOEFB) , RA(IFRCX ) ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 8.1 MODULE DE RAYONNEMENT : CALCUL SUITE EVENTUEL
C=======================================================================
C
      IF (IIRAYO.GT.0 .AND. ISUIRD.EQ.1) THEN
C
        CALL RAYLEC
C       ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDEVEL , ITUSER , IA     ,
     &   RA(ITSRE)  , RA(ITSRI)  ,
     &   RA(ITPARO) , RA(IQINCI) ,  RA(IFCONV) , RA(IHCONV) ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C=======================================================================
C 8.2 INITIALISATIONS DES PARTICULES POUR LE LAGRANGIEN
C=======================================================================
C
      IF (IILAGR.GT.0) THEN
C
        CALL LAGLEC
C       ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   IA(IIITEP) , IA ,
     &   RA(IRTPA)  , RA(IPROPC) ,
     &   RA(IETTP)  , RA(IITEPA) , RA(ISTATC) , RA(ISTATV) ,
     &   RA(ISTATF) , RA(ITSLAG) , RA         )
C
      ENDIF
C
C=======================================================================
C 8.3 INITIALISATIONS POUR LE MODULE THERMIQUE 1D EN PAROI
C=======================================================================
C On suppose que toutes les phases voient la meme temperature de paroi
C USPT1D a un fonctionnement similaire a USKPDC et USTSMA, mais comme
C on ecrit des infos dans un fichier suite, on a besoin d'une partie de
C la memoire meme apres la boucle en temps -> IFPT1D et TPPT1D
C                                            (IFNIA1 et IFNRA1)
C
      IDBIA1 = IFINIA
      IDBRA1 = IFINRA
C
      ILS    = IDBIA1
      IFNIA2 = ILS + MAXELT
      CALL IASIZE('CALTRI',IFNIA2)
C
      IPHAS = 1
C
C     Premier appel : definition de NFPT1D et ISUIT1
      IAPPEL = 1
      CALL USPT1D
C     ===========
     & ( IFNIA2 , IDBRA1 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NFPT1D , IPHAS  , IAPPEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IDBIA1) , IA(IDBIA1) , IA(IDBIA1) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDBRA1) , RA(IDBRA1) , RA(IDBRA1) ,
     &   RA(IDBRA1) , RA(IDBRA1) , RA(IDBRA1) ,
     &   RA(IDBRA1) , RA(IDBRA1) , RA(IDBRA1) ,
     &   RA(IDT)    , RA(IRTPA)  ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RDEVEL , RTUSER , RA     )
C
      IAPPEL = 1
      CALL VERT1D
C     ===========
     & (IDBIA1     , IDBRA1     ,
     &  NFABOR     , NFPT1D     , IAPPEL     ,
     &  IA(IDBIA1) , IA(IDBIA1) , IA(IDBIA1) , IA     ,
     &  RA(IDBRA1) , RA(IDBRA1) ,
     &  RA(IDBRA1) , RA(IDBRA1) , RA(IDBRA1) , RA     )
C
      CALL MEMT1D
C     ===========
     & ( IDBIA1 , IDBRA1 , NFABOR , IFNIA1 , IFNRA1 ,IFNIA2 , IFNRA2 ,
     &   IFINIA , IFINRA , IA     , RA     )
C
C On appelle uspt1d lorqu'il y a sur un processeur au moins des faces de
C     bord avec module thermique 1D.
C
      IF (NFPT1T.GT.0) THEN
C Deuxieme appel : remplissage des tableaux de definition de la geometrie
C            et de l'initialisation (IFPT1D,NPPT1D,EPPT1D,RGPT1D,TPPT1D)
        ILS    = IFINIA
        IFNIA3 = ILS + MAXELT
        CALL IASIZE('CALTRI',IFNIA3)
C
        IAPPEL = 2
        CALL  USPT1D
C       ============
     & ( IFNIA3 , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NFPT1D , IPHAS  , IAPPEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IIFPT1) , IA(INPPT1) , IA(IICLT1) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(ITPPT1) , RA(IRGPT1) , RA(IEPPT1) ,
     &   RA(ITEPT1) , RA(IHEPT1) , RA(IFEPT1) ,
     &   RA(IXLMT1) , RA(IRCPT1) , RA(IDTPT1) ,
     &   RA(IDT)    , RA(IRTPA)  ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RDEVEL , RTUSER , RA     )
C
        IAPPEL = 2
        CALL VERT1D
C       ===========
     & (IFINIA     , IFINRA     ,
     &  NFABOR     , NFPT1D     , IAPPEL     ,
     &  IA(IIFPT1) , IA(INPPT1) , IA(IICLT1) , IA     ,
     &  RA(IRGPT1) , RA(IEPPT1) ,
     &  RA(IXLMT1) , RA(IRCPT1) , RA(IDTPT1) , RA     )
C
C     Calcul du max des NPPT1D (pour les fichiers suite)
        NMXT1D = 0
        DO IIII = 1, NFPT1D
          NMXT1D = MAX(IA(INPPT1+IIII-1),NMXT1D)
        ENDDO
        IF (IRANGP.GE.0) CALL PARCMX(NMXT1D)
C                        ===========
C
        IF (ISUIT1.EQ.1) THEN
C
          CALL LECT1D
C         ===========
     & ( FICMT1     , LEN(FICMT1), IFOVT1     , NFPT1D     , NFPT1T    ,
     &   NMXT1D     , NFABOR     , IA(INPPT1) , IA(IIFPT1) , RA(IEPPT1),
     &   RA(IRGPT1) , RA(ITPPT1))
C
        ELSE
C     Creation du maillage, initialisation de la temperature.
C
          CALL MAIT1D
C         ===========
     & ( NFPT1D, IA(INPPT1), RA(IEPPT1), RA(IRGPT1),RA(ITPPT1))
C
        ENDIF
C
      ENDIF
C     Les infos sur l'epaisseur de la paroi, le nombre de points de
C     discretisation et la raison geometrique ont ete transmises a
C     la structure C. Elles sont maintenant inutiles dans le Fortran.
C     -> on libere la memoire.
      IFINIA = IFNIA2
      IFINRA = IFNRA2
C
C=======================================================================
C 9. TABLEAUX POUR BLC EN TEMPS MAIS A OUBLIER ENSUITE
C=======================================================================
C
C  En fin de bloc en temps on doit retrouver IFNIA1 et IFNRA1
      IDITIA = IFNIA1
      IDITRA = IFNRA1
C
      IDBIA1 = IFINIA
      IDBRA1 = IFINRA
C
C
      DO IPHAS = 1, NPHAS
C
        IAPPEL = 1
C
        IF (IMATIS.EQ.1) THEN
C
C     Noter que uskpdc n'est pas permis avec Matisse
C       (si necessaire, on pourrait regrouper toutes les pertes de
C        charge de Matisse dans un unique sous-programme et reactiver
C        uskpdc)
C
          CALL  MTKPDC
C         ============
     & ( IDBIA1 , IDBRA1 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   NCEPDC(IPHAS) , NCKPDC(IPHAS)     , IPHAS  , IAPPEL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IDBIA1),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) , RA(IDBRA1) ,
     &   RDEVEL , RTUSER , RA     )
C
        ELSE
C
          ILS    = IDBIA1
          IDBIA2 = ILS + MAXELT
          CALL IASIZE('CALTRI',IDBIA2)
C
          CALL  USKPDC
C         ============
     & ( IDBIA2 , IDBRA1 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   NCEPDC(IPHAS) , NCKPDC(IPHAS)     , IPHAS  , IAPPEL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IDBIA1),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) , RA(IDBRA1) ,
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF
C
      ENDDO
C
      CALL MEMPDC
C     ===========
     & ( IDBIA1, IDBRA1, NCELET, NCEL,  NPHAS, NDIM, IFINIA, IFINRA)
C
C
C On appelle uskpdc lorqu'il y a sur un processeur au moins des cellules
C     avec terme source de masse.
C     On ne fait que remplir le tableau d'indirection des cellules
C     On appelle cependant uskpdc avec tous les processeurs, au cas ou
C     l'utilisateur aurait mis en oeuvre des operations globales.
C
      DO IPHAS = 1, NPHAS
C
        IF(NCPDCT(IPHAS).GT.0) THEN
C
          IAPPEL = 2
C
          IF (IMATIS.EQ.1) THEN
C
            CALL  MTKPDC
C           ============
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   NCEPDC(IPHAS) , NCKPDC(IPHAS)     , IPHAS  , IAPPEL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) , RA(ICKUPD(IPHAS)) ,
     &   RDEVEL , RTUSER , RA     )
C
          ELSE
C
            ILS    = IFINIA
            IFNIA2 = ILS + MAXELT
            CALL IASIZE('CALTRI',IFNIA2)
C
            CALL  USKPDC
C           ============
     & ( IFNIA2 , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   NCEPDC(IPHAS) , NCKPDC(IPHAS)     , IPHAS  , IAPPEL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) , RA(ICKUPD(IPHAS)) ,
     &   RDEVEL , RTUSER , RA     )
C
          ENDIF
C
        ENDIF
C
      ENDDO
C
      IDBIA1 = IFINIA
      IDBRA1 = IFINRA
C
      ILS    = IDBIA1
      IDBIA2 = ILS + MAXELT
      CALL IASIZE('CALTRI',IDBIA2)
C
      DO IPHAS = 1, NPHAS
C
        IAPPEL = 1
        CALL  USTSMA
C       ============
     & ( IDBIA2 , IDBRA1 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDC(IPHAS)   , NCKPDC(IPHAS)   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   NCETSM(IPHAS) ,   IPHAS  , IAPPEL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) ,
     &   IA(IDBIA1) , IA(IDBIA1),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) , RA(ICKUPD(IPHAS))       , RA(IDBRA1),
     &   RDEVEL , RTUSER , RA     )
C
      ENDDO
C
      CALL MEMTSM
C     ===========
     &     ( IDBIA1 , IDBRA1 ,
     &       NCELET , NCEL   , NVAR   , NPHAS  ,
     &       IFINIA , IFINRA )
C
C On appelle ustsma lorqu'il y a sur un processeur au moins des cellules
C     avec terme source de masse.
C     On ne fait que remplir le tableau d'indirection des cellules
C     On appelle cependant ustsma avec tous les processeurs, au cas ou
C     l'utilisateur aurait mis en oeuvre des operations globales.
C
      DO IPHAS = 1, NPHAS
C
        IF(NCTSMT(IPHAS).GT.0) THEN
C
          ILS    = IFINIA
          IFNIA2 = ILS + MAXELT
          CALL IASIZE('CALTRI',IFNIA2)
C
          IAPPEL = 2
          CALL  USTSMA
C         ============
     & ( IFNIA2 , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDC(IPHAS)   , NCKPDC(IPHAS)   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   NCETSM(IPHAS) ,   IPHAS  , IAPPEL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS)) ,
     &   IA(IICESM(IPHAS)) , IA(IITPSM(IPHAS)),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) , RA(ICKUPD(IPHAS)), RA(ISMACE(IPHAS)),
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF
C
      ENDDO
C
C
C -- Methode des vortex pour la L.E.S.
C    (dans verini on s'est deja assure que ITYTUR=4 si IVRTEX=1)
C
      IF (IVRTEX.EQ.1) THEN
C
        IDBIA1 = IFINIA
        IDBRA1 = IFINRA
C
        IPHAS  = 1
        IAPPEL = 1
C
C  On met une valeur factice a certains parametres non utilise en IAPPEL=1
C
        CALL MEMVOR
C       ===========
     & ( IDBIA1 , IDBRA1 , IAPPEL , NFABOR , IFINIA , IFINRA )
C
        CALL VORIN0( NFABOR , IA(IIREPV) )
C       ===========
C
        ILS    = IFINIA
        IFNIA2 = ILS + MAXELT
        CALL IASIZE('CALTRI',IFNIA2)
C
        CALL USVORT
C       ===========
     & ( IFNIA2 , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPHAS  , IAPPEL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IIREPV)      ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RDEVEL , RTUSER , RA     )
C
        CALL VORVER ( NFABOR , IA(IIREPV)  , IAPPEL )
C       ===========
C
        IDBIA1 = IFINIA
        IDBRA1 = IFINRA
C
C Attention, vorpre reserve de la memoire qu'il faut garder ensuite
C           (-> on utilise IFINIA/IFINRA ensuite)
C
        CALL VORPRE
C       ===========
     & ( IDBIA1 , IDBRA1 , IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IIREPV),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD ,
     &   VOLUME , RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C -- Fin de zone Methode des vortex pour la L.E.S.
C
C -- Structures mobiles en ALE
C
      IF (IALE.EQ.1) THEN
C
        IDBIA1 = IFINIA
        IDBRA1 = IFINRA
C
C Attention, strini reserve de la memoire qu'il faut garder ensuite
C           (-> on utilise IFINIA/IFINRA ensuite)
        CALL STRINI
C       ===========
     & ( IDBIA1 , IDBRA1 , IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME , RA(IDT),
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C -- Fin de zone Structures mobiles en ALE
C
C
C     EN CAS DE COUPLAGE, ON LIT DES MAINTENANT L'ENTETE DU PREMIER
C     MESSAGE, AU CAS OU IL S'AGISSE D'UN MESSAGE DE TERMINAISON
C
      CALL TSTSYR (IMSFIN, NTMABS, NTCABS)
C     ===========
C
      IF (IMSFIN.EQ.1) THEN
        CALL CSEXIT (0)
C       ===========
      ENDIF
C
C=======================================================================
C 10. DEBUT DE LA BOUCLE EN TEMPS
C=======================================================================
C
      WRITE(NFECRA,2000)
C
      NTCABS = NTPABS
      TTCABS = TTPABS
C
      IWARN0 = 1
      DO IVAR = 1, NVAR
        IWARN0 = MAX(IWARN0,IWARNI(IVAR))
      ENDDO
C
      IF(IWARN0.GT.0) THEN
        WRITE(NFECRA,3000)
      ENDIF
C
      IF(INPDT0.EQ.1) THEN
        NTMABS = NTCABS
      ENDIF
C
C     Nb d'iter ALE (nb relatif a l'execution en cours)
C     Si ITALIN=1, on fait une iteration d'initialisation
C     (si ITALIN=-999, c'est qu'on a fait une suite de calculs
C      sans relire lecamx -> comme si on faisait une suite
C      d'un calcul sans ALE)
      IF (ITALIN.EQ.-999) ITALIN = 1
      ITRALE = 1
      IF (ITALIN.EQ.1) THEN
        ITRALE = 0
        WRITE(NFECRA,3002) TTCABS
      ENDIF
C
 100  CONTINUE
C
      IF(INPDT0.EQ.0 .AND. ITRALE.GT.0) THEN
        NTCABS = NTCABS + 1
        IF(IDTVAR.EQ.0.OR.IDTVAR.EQ.1) THEN
          TTCABS = TTCABS + RA(IDT)
        ELSE
          TTCABS = TTCABS + DTREF
        ENDIF
        IF(IWARN0.GT.0) THEN
          WRITE(NFECRA,3001) TTCABS,NTCABS
        ENDIF
      ENDIF
C
C
C=======================================================================
C 11. AVANCEE EN TEMPS
C=======================================================================
C
C
C     On teste la presence de ficstp pour modifier NTMABS le cas echeant
      CALL MODPAR(NTCABS,NTMABS)
C     ===========
C
      CALL DMTMPS(TITER1)
C     ===========
C
C     Synchronisation Syrthes, si ITRALE>0
      IF (ITRALE.GT.0) THEN
        CALL ITDSYR(NTCABS,NTMABS)
C       ===========
      ENDIF
C
      CALL TRIDIM
C     ===========
     & ( IFINIA , IFINRA , ITRALE ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IISSTD),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ITSLAG) , RA(ICOEFA) , RA(ICOEFB) ,
     &   RA(IFRCX)  ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 12. CALCUL DES MOYENNES TEMPORELLES (ACCUMULATION)
C=======================================================================
C
C
      IF(INPDT0.EQ.0 .AND. ITRALE.GT.0) THEN
      CALL CALMOM
C     ===========
     &     ( IFINIA , IFINRA , NCEL   , NCELET ,
     &       NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &       IDEVEL , ITUSER , IA     ,
     &       RA(IRTP  ) , RA(IDT   ) , RA(IPROPC) ,
     &       RDEVEL , RTUSER , RA     )
      ENDIF
C
C
C=======================================================================
C 13. APPEL DU MODULE LAGRANGIEN
C=======================================================================
C
      IF (IILAGR.GT.0 .AND. INPDT0.EQ.0 .AND. ITRALE.GT.0) THEN
C
        CALL MEMLA2
C       ===========
     & ( IFINIA , IFINRA ,
     &   NFABOR , NCELET , NFAC   ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   IINDEP , IIBORD , IETTPA , IAUXL  , IAUXL2 ,
     &   ITAUP  , IITLAG , IPIIL  ,
     &   IVAGAU , ITSUF  , ITSUP  , IBX    ,
     &   IGRADP , IGRADV , ICROUL ,
     &   ITEPCT , ITSFEX , ITSVAR ,
     &   ICPGD1 , ICPGD2 , ICPGHT ,
     &   IBRGAU , ITEBRU ,
     &   IW1    , IW2    , IW3    ,
     &   ILAGIA , ILAGRA  )
C
        CALL LAGUNE
C       ===========
     & ( ILAGIA , ILAGRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICOCE) , IA(IITYCE) , IA(IIFRLA) , IA(IIITEP) ,
     &   IA(IINDEP) , IA(IIBORD) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RA(IETTP)  , RA(IETTPA) , RA(IITEPA) , RA(ISTATC) , RA(ISTATV),
     &   RA(ITSLAG),
     &   RA(ISTATF) , RA(ITAUP)  , RA(IITLAG) , RA(IPIIL)  , RA(IBX  ) ,
     &   RA(IVAGAU) , RA(ITSUF ) , RA(ITSUP ) , RA(ITSVAR) ,
     &   RA(ITEPCT) , RA(ITSFEX) ,
     &   RA(ICPGD1) , RA(ICPGD2) , RA(ICPGHT) ,
     &   RA(IGRADP) , RA(IGRADV) , RA(ICROUL) ,
     &   RA(IBRGAU) , RA(ITEBRU) ,
     &   RA(IW1   ) , RA(IW2   ) , RA(IW3   ) , RA(IAUXL)  , RA(IAUXL2),
     &   RDEVEL , RTUSER , RA     )
C
C--> Ici on libere la memoire reserve par MEMLA2
C      (i.e. on oublie ILAGIA et ILAGRA)
C
      ENDIF
C
C=======================================================================
C 14. BRANCHEMENT UTILISATEUR POUR MODIF DES VARIABLES EVENTUELLES
C=======================================================================
C
C     Appel pour Matisse d'une routine de bilans, d'impression, ...
      IF (IMATIS.EQ.1 .AND. ITRALE.GT.0) THEN
C
        CALL MTPROJ
C       ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVEP   , NIVEP  , NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IIITEP),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RA(IETTP)  , RA(IETTPA) , RA(IITEPA) , RA(ISTATC) , RA(ITSLAG),
     &   RA(ISTATF) ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
      IF (ITRALE.GT.0) THEN
C
        ILS    = IFINIA
        IFNIA2 = ILS + MAXELT
        CALL IASIZE('CALTRI',IFNIA2)
C
        CALL USPROJ
C       ===========
     & ( IFNIA2 , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVEP   , NIVEP  , NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IIITEP),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RA(IETTP)  , RA(IETTPA) , RA(IITEPA) , RA(ISTATC) , RA(ISTATV),
     &   RA(ITSLAG) , RA(ISTATF) ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C=======================================================================
C 15. MISE A JOUR DU MAILLAGE (ALE)
C=======================================================================
C
      IF (IALE.EQ.1 .AND. INPDT0.EQ.0) THEN
C
        IF (ITRALE.EQ.0 .OR. ITRALE.GT.NALINF) THEN
C
          CALL ALEMAJ
C         ===========
     & ( IFINIA , IFINRA , ITRALE ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IIMPAL)      ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RA(IDEPAL) , RA(IXYZN0) ,
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF
C
      ENDIF
C
C=======================================================================
C 16. TEST D'ARRET PAR MANQUE DE TEMPS
C=======================================================================
C
      CALL ARMTPS(NTCABS,NTMABS)
C     ===========
C
C=======================================================================
C 17. TEST D'ARRET ISSU DU MODULE DE RAYONNEMENT P-1
C=======================================================================
      IF (ISTPP1.EQ.1) THEN
        NTMABS = NTCABS
      ENDIF
C
C=======================================================================
C 18. TEST D'ARRET PAR DEMANDE DE SYRTHES
C=======================================================================
C
C     EN CAS DE COUPLAGE, ON LIT DES MAINTENANT L'ENTETE DU PREMIER
C     MESSAGE, AU CAS OU IL S'AGISSE D'UN MESSAGE DE TERMINAISON
C
      IF (NTCABS.LT.NTMABS) THEN
        CALL TSTSYR (IMSFIN, NTMABS, NTCABS)
C       ===========
      ENDIF
C
C=======================================================================
C 19. SORTIE EVENTUELLE DU FICHIER SUITE
C      (SAUF SI ON EST AU DERNIER PAS DE TEMPS :
C       ON LE FERA APRES)
C=======================================================================
C
      IISUIT = 0
      IF(NTCABS.LT.NTMABS) THEN
        IF(NTSUIT.EQ.0) THEN
          NTSDEF = MAX((NTMABS-NTPABS)/4,10)
          IF(NTSDEF.GT.0) THEN
            IF(MOD(NTCABS-NTPABS,NTSDEF).EQ.0) THEN
              IISUIT = 1
            ENDIF
          ENDIF
        ELSEIF(NTSUIT.GT.0) THEN
          IF(MOD(NTCABS,NTSUIT).EQ.0) THEN
            IISUIT = 1
          ENDIF
        ENDIF
      ENDIF
      IF (ITRALE.EQ.0) IISUIT = 0
C
      IF(IISUIT.EQ.1) THEN
        CALL ECRAVA
C       ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NNOD   ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN     , SURFAC     , SURFBO     , CDGFAC     , CDGFBO    ,
     &   RA(IDT)    , RA(IRTP)   , RA(IPROPC) , RA(IPROPF) , RA(IPROPB),
     &   RA(ICOEFA) , RA(ICOEFB) , RA(IFRCX)  ,
     &   RDEVEL , RTUSER , RA     )
C
        IF (NFPT1T.GT.0) THEN
          CALL ECRT1D
C         ===========
     & ( FICVT1   , LEN(FICVT1), IFOVT1    , NFPT1D   ,  NMXT1D  ,
     &   NFABOR   , RA(ITPPT1) , IA(IIFPT1))
        ENDIF
C
        IF (IILAGR.GT.0) THEN
C
          CALL LAGOUT
C         ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICOCE) , IA(IITYCE) , IA(IIITEP) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RA(IETTP)  , RA(IITEPA) , RA(ISTATF) ,
     &   RA(ISTATC) , RA(ISTATV) , RA(ITSLAG) ,
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF
C
        IF (IIRAYO.GT.0) THEN
          CALL RAYOUT
C         ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RA(ITSRE)  , RA(ITSRI)  ,
     &   RA(ITPARO) , RA(IQINCI) , RA(IXLAM)  , RA(IEPA)   , RA(IEPS)  ,
     &   RA(IFNET)  , RA(IFCONV) , RA(IHCONV) ,
     &   RDEVEL , RTUSER , RA     )
        ENDIF
C
        IF(IWARN0.GT.0) THEN
          WRITE(NFECRA,3020)NTCABS,TTCABS
        ENDIF
C
      ENDIF
C
C=======================================================================
C 20. TEST POUR SAVOIR SI ON SORT UN FICHIER POST OU NON
C=======================================================================
C
      CALL USNPST
C     ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NVLSTA ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) , RA(ISTATC) ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 21. SORTIE DES FICHIERS POST STANDARDS
C=======================================================================
C
C     Si ITRALE=0 on desactive tous les writers (car la geometrie n'a pas ete
C       ecrite)
      IF (ITRALE.EQ.0) THEN
        INDWRI = 0
        INDACT = 0
        CALL PSTACT(INDWRI, INDACT)
C       ===========
      ENDIF
C
      CALL PSTVAR
C     ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NTCABS , NCELET , NCEL   , NFAC   , NFABOR ,
     &   NFML   , NPRFML , NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   TTCABS , XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RA(ISTATC) , RA(ISTATV) , RA(ISTATF) ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 22. HISTORIQUES
C=======================================================================
C
C ON STOCKE SUR TMP SI ON A DECIDE DE LE FAIRE ET QUE C'EST LE MOMENT
C ON SAUVE SI ON A DECIDE DE LE FAIRE, QUE C'EST LE MOMENT ET
C   QUE CE N'EST PAS LE DERNIER PAS DE TEMPS
C
      IF(NTHIST.GT.0 .AND. ITRALE.GT.0) THEN
        IF(MOD(NTCABS,NTHIST).EQ.0) THEN
          MODHIS = 0
          IF(NTHSAV.GT.0) THEN
            IF(MOD(NTCABS,NTHSAV).LT.NTHIST ) MODHIS = 1
          ELSEIF(NTHSAV.EQ.0) THEN
            NTCREL =  NTCABS-NTPABS
            NTHDEF = (NTMABS-NTPABS)/4
            IF( (NTCREL-10).GE.0 .AND. (NTCREL-10).LT.NTHIST) THEN
              MODHIS = 1
            ELSEIF(NTHDEF.GT.0) THEN
              IF(MOD(NTCREL,NTHDEF).LT.NTHIST) THEN
                MODHIS = 1
              ENDIF
            ENDIF
          ENDIF
          IF (NTCABS.EQ.NTMABS) MODHIS=0

          CALL ECRHIS
C         ===========
     &   (IFINIA , IFINRA , NDIM   , NCELET , NCEL,
     &    NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &    MODHIS ,
     &    IDEVEL , ITUSER , IA     ,
     &    XYZCEN ,
     &    RDEVEL , RTUSER , RA     )
C
          IF (IILAGR.GT.0) THEN
            CALL LAGHIS
C           ===========
     &     (IFINIA , IFINRA , NDIM   , NCELET , NCEL,
     &      NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &      MODHIS , NVLSTA ,
     &      IDEVEL , ITUSER , IA     ,
     &      XYZCEN , VOLUME ,
     &      RA(ISTATC) , RA(ISTATV) ,
     &      RDEVEL , RTUSER , RA     )
          ENDIF
C
          IF (IHISTR.EQ.1) THEN
            CALL STRHIS
C           ===========
     &     (IFINIA , IFINRA , NCELET , NCEL,
     &      NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &      MODHIS ,
     &      IDEVEL , ITUSER , IA     ,
     &      RDEVEL , RTUSER , RA     )
          ENDIF
C
        ENDIF
      ENDIF
C
      CALL USHIST
C     ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 23. ECRITURE LISTING TOUTES LES NTLIST ITERATIONS
C=======================================================================
C
C
      IF(NTLIST.GT.0) THEN
        MODNTL = MOD(NTCABS,NTLIST)
      ELSEIF(NTLIST.EQ.-1.AND.NTCABS.EQ.NTMABS) THEN
        MODNTL = 0
      ELSE
        MODNTL = 1
      ENDIF
      IF(MODNTL.EQ.0) THEN
        CALL ECRLIS
C       ===========
     &  ( IFINIA , IFINRA ,
     &    NVAR   , NPHAS  , NDIM   , NCELET , NCEL   ,
     &    NIDEVE , NRDEVE , NITUSE , NRTUSE , IRTP   ,
     &    IDEVEL , ITUSER , IA     ,
     &    RA(IRTP  ) , RA(IRTPA ) , RA(IDT ) , VOLUME , XYZCEN,
     &    RDEVEL , RTUSER , RA     )
C
        IF (IILAGR.GT.0) THEN
C
          CALL LAGLIS
C         ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IIITEP),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RA(IETTP)  , RA(IITEPA) , RA(ISTATC) ,  RA(ISTATV) ,
     &   RA(ITSLAG) , RA(ISTATF),
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF
C
      ENDIF
C
      CALL DMTMPS(TITER2)
C     ===========
C
      IF(IWARN0.GT.0) THEN
        IF (ITRALE.GT.0) THEN
          WRITE(NFECRA,3010)NTCABS,TITER2-TITER1
        ELSE
          WRITE(NFECRA,3012)TITER2-TITER1
        ENDIF
      ENDIF
C
C
C=======================================================================
C 24. FIN DE LA BOUCLE EN TEMPS
C=======================================================================
C
      ITRALE = ITRALE + 1
C
      IF(NTCABS.LT.NTMABS) GOTO 100
C
C
C LIBERATION DES TABLEAUX INTERMEDIAIRES (PDC+TSM)
C
      IFINIA = IDITIA
      IFINRA = IDITRA
C
C
C=======================================================================
C 25. ECRITURE DES SUITES + FINALISATION HISTORIQUES
C=======================================================================
C
      CALL DMTMPS(TECRF1)
C     ===========
C
      IF(IWARN0.GT.0) THEN
        WRITE(NFECRA,4000)
      ENDIF
C
      CALL ECRAVA
C     ===========
     &( IFINIA , IFINRA ,
     &  NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NNOD   ,
     &  NVAR   , NSCAL  , NPHAS  ,
     &  NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &  IDEVEL , ITUSER , IA     ,
     &  XYZCEN     , SURFAC     , SURFBO     , CDGFAC     , CDGFBO    ,
     &  RA(IDT)    , RA(IRTP)   , RA(IPROPC) , RA(IPROPF) , RA(IPROPB),
     &  RA(ICOEFA) , RA(ICOEFB) , RA(IFRCX)  ,
     &  RDEVEL , RTUSER , RA     )
C
      IF (NFPT1T.GT.0) THEN

        CALL ECRT1D
C       ===========
     & ( FICVT1   , LEN(FICVT1), IFOVT1    , NFPT1D   ,  NMXT1D  ,
     &   NFABOR   , RA(ITPPT1) , IA(IIFPT1))
      ENDIF
C
      IF (IILAGR.GT.0) THEN
C
        CALL LAGOUT
C       ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICOCE) , IA(IITYCE) , IA(IIITEP) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RA(IETTP)  , RA(IITEPA) , RA(ISTATF) ,
     &   RA(ISTATC) , RA(ISTATV) , RA(ITSLAG) ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C Ecriture du fichier suite
      IF (IIRAYO.GT.0) THEN
C
        CALL RAYOUT
C       ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   RA(IDT)    , RA(IRTPA)  , RA(IRTP)   ,
     &   RA(IPROPC) , RA(IPROPF) , RA(IPROPB) ,
     &   RA(ICOEFA) , RA(ICOEFB) ,
     &   RA(ITSRE)  , RA(ITSRI)  ,
     &   RA(ITPARO) , RA(IQINCI) , RA(IXLAM)  , RA(IEPA)   , RA(IEPS)  ,
     &   RA(IFNET)  , RA(IFCONV) , RA(IHCONV) ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C ICI ON SAUVE LES HISTORIQUES (SI ON EN A STOCKE)
C
      MODHIS = 2
      CALL ECRHIS
C     ===========
     &(  IFINIA , IFINRA , NDIM   , NCELET , NCEL,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   MODHIS ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN ,
     &   RDEVEL , RTUSER , RA     )
C
      IF (IILAGR.GT.0) THEN
          CALL LAGHIS
C         ===========
     &   (IFINIA , IFINRA , NDIM   , NCELET , NCEL,
     &    NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &    MODHIS , NVLSTA ,
     &   IDEVEL , ITUSER , IA     ,
     &    XYZCEN , VOLUME ,
     &    RA(ISTATC) , RA(ISTATV) ,
     &   RDEVEL , RTUSER , RA     )
      ENDIF
      IF (IHISTR.EQ.1) THEN
        CALL STRHIS
C       ===========
     & (IFINIA , IFINRA , NCELET , NCEL,
     &  NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &  MODHIS ,
     &  IDEVEL , ITUSER , IA     ,
     &  RDEVEL , RTUSER , RA     )
      ENDIF
C
C     LE CAS ECHEANT, ON LIBERE LES STRUCTURES C DU MODULE THERMIQUE 1D
C     ET/OU ON FERME LE LISTING LAGRANGIEN
C
      IF (NFPT1D.GT.0) THEN
        CALL LBRT1D
C       ===========
      ENDIF
C
      IF (IALE.GT.0) THEN
        CALL LBRALE
C       ===========
      ENDIF
C
      IF (IILAGR.GT.0) CLOSE(IMPLAL)
C
      CALL DMTMPS(TECRF2)
C     ===========
C
      IF(IWARN0.GT.0) THEN
        WRITE(NFECRA,4010)TECRF2-TECRF1
      ENDIF
C
 200  CONTINUE
C
C=======================================================================
C 26. MEMOIRE UTILISEE
C=======================================================================
C
#if defined(_CS_HAVE_XML)
C     Libération des structures liées à la lecture du fichier xml
C
      IF (IIHMPR.EQ.1) THEN
C
        CALL MEMUI2
C       ===========
C
        CALL MEMUI1(NCHARB)
C       ===========
      ENDIF
C
#endif
C
      WRITE(NFECRA,7000)
C
C----
C FORMATS
C----
C
C
 1000 FORMAT(/,
     &'***************************************************************',
     &                                                              /,/,
     &'              FONCTIONNEMENT EN MODE VERIFICATION            ',/,
     &'              ===================================            ',/,
     &'                                                             ',/,
     &' MODE IVERIF = ',I7                                           ,/,
     &'                                                             ',/,
     &' =========================================================== ',/,
     &                                                                /,
     &                                                                /)
C
 2000 FORMAT(/,/,
     &'***************************************************************',
     &                                                              /,/,
     &                                                                /,
     &                                                                /,
     &'                       CORPS DU CALCUL                       ',/,
     &'                       ===============                       ',/,
     &                                                                /,
     &                                                                /,
     &'***************************************************************',
     &                                                              /,/,
     &                                                                /)
 3000 FORMAT(/,
     &'***************************************************************',
     & /)
 3001 FORMAT(/,' INSTANT ',E18.9,        '   ITERATION NUMERO ',I15,/,
     &' ============================================================= ',
     & /,/)
 3002 FORMAT(/,' INSTANT ',E18.9,        '   INITIALISATION ALE ',/,
     &' ============================================================= ',
     & /,/)
 3010 FORMAT(/,' TEMPS CPU POUR L''ITERATION ',I15,' :    ',E14.5,/,/,
     &'***************************************************************',
     & /)
 3012 FORMAT(/,' TEMPS CPU POUR L''INITIALISATION ALE :    ',E14.5,/,/,
     &'***************************************************************',
     & /)
 3020 FORMAT(/,/,
     & ' Sortie intermediaire d''un fichier suite ',/,
     & '   Sauvegarde a l''iteration ',    I10,
     &                                    ', Temps physique ',E14.5,/,/)
C
 4000 FORMAT(/,/,
     &'***************************************************************',
     &                                                              /,/,
     &                                                                /,
     &                                                                /,
     &'                   ETAPES FINALES DU CALCUL                  ',/,
     &'                   ========================                  ',/,
     &                                                                /,
     &                                                                /,
     &' =========================================================== ',/,
     &                                                                /,
     &                                                                /)
 4010 FORMAT(                                                         /,
     & 3X,'** TEMPS CPU POUR LES SORTIES FINALES : ',E14.5           ,/,
     & 3X,'   ----------------------------------                    ',/)
 7000 FORMAT(/,/,
     &' =========================================================== ',/,
     &                                                              /,/,
     &                                                                /,
     &                                                                /,
     &'                 FIN DE L''EXECUTION DU CALCUL               ',/,
     &'                 ============================                ',/,
     &                                                                /,
     &                                                                /,
     &'***************************************************************')
C
C=======================================================================
C 26. FIN
C=======================================================================
C
C
      RETURN
      END
c@z
