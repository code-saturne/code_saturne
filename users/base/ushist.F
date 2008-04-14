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
                        SUBROUTINE USHIST
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC SORTIE D'HISTORIQUES NON STD LIVREE A L'UTILISATEUR
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
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
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
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          II, KK, NODE, NDRANG, NVARPP, NUMCEL, LNG
      DOUBLE PRECISION XX, YY, ZZ, XYZTMP(3)
C
C   Numero des noeuds ou on sort des historiques
      INTEGER          NCAPMX
      PARAMETER       (NCAPMX=100)
      INTEGER          ICAPT(NCAPMX)
      SAVE             ICAPT
      INTEGER          IRCAPT(NCAPMX)
      SAVE             IRCAPT
C
C   Nombre de noeuds ou on sort des historiques
      INTEGER          NCAPTS
      SAVE             NCAPTS
C
C   Numero du passage actuel dans ce ss pgm
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C   Tableau de valeurs temporaires
      DOUBLE PRECISION VACAPT(NCAPMX)
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
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
C ---> Gestion memoire
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C ---> Numero du passage actuel dans ce ss pgm
C
      IPASS = IPASS + 1
C
C=======================================================================
C 2. RECHERCHE DES CAPTEURS
C=======================================================================
C Les numeros stockes dans IRCAPT donnent le rang du processeur sur
C   lequel se trouve la sonde. L'utilisateur n'a pas a s'en preoccuper
C   specialement tant qu'il utilise bien la fonction FINDPT pour reperer
C   les sondes.
C
C
C  Au premier passage : reperage des numeros de cellule dont le centre est
C    le plus proche des coordonnees XX YY ZZ.
C    En parallelisme, le numero de cellule ICAPT(II) est local au processeur
C    dont le rang est donne par IRCAPT(II) (de 0 a nombre de processeurs-1).
C    NCAPTS donne le nombre de sondes total.
C
      IF (IPASS.EQ.1) THEN
C
        II = 0
C
        XX = 0.20D0
        YY = 0.15D0
        ZZ = 0.01D0
        CALL FINDPT
C       ===========
     & ( NCELET , NCEL   , XYZCEN ,
     &   XX     , YY     , ZZ     , NODE  , NDRANG)
        II = II + 1
        ICAPT(II) = NODE
        IRCAPT(II) = NDRANG
C
        XX = 0.70D0
        YY = 0.15D0
        ZZ = 0.01D0
        CALL FINDPT
C       ===========
     & ( NCELET , NCEL   , XYZCEN ,
     &   XX     , YY     , ZZ     , NODE  , NDRANG)
        II = II + 1
        ICAPT(II) = NODE
        IRCAPT(II) = NDRANG
C
        XX = 0.20D0
        YY = 0.75D0
        ZZ = 0.01D0
        CALL FINDPT
C       ===========
     & ( NCELET , NCEL   , XYZCEN ,
     &   XX     , YY     , ZZ     , NODE  , NDRANG)
        II = II + 1
        ICAPT(II) = NODE
        IRCAPT(II) = NDRANG
C
        XX = 0.70D0
        YY = 0.75D0
        ZZ = 0.01D0
        CALL FINDPT
C       ===========
     & ( NCELET , NCEL   , XYZCEN ,
     &   XX     , YY     , ZZ     , NODE  , NDRANG)
        II = II + 1
        ICAPT(II) = NODE
        IRCAPT(II) = NDRANG
C
        NCAPTS = II
C
        IF(II.GT.NCAPMX) THEN
          WRITE(NFECRA,*) ' USHIST : NCAPMX = ',II,' AU MINIMUM '
          CALL CSEXIT (1)
        ENDIF
C
      ENDIF
C
C
C=======================================================================
C 3. OUVERTURE DES FICHIERS
C     EXEMPLE D'UNE VARIABLE PAR FICHIER
C=======================================================================
C
C ---> Nombre de variables = nombre de fichiers
C
      NVARPP = NVAR
C
C
C ---> Au premier passage, on ouvre les fichiers et on ecrit une entete
C
      IF(IPASS.EQ.1) THEN
C
C   --> Test du nombre max de fichiers
C
        IF(NVARPP.GT.NUSHMX) THEN
          WRITE(NFECRA,*)
     &  ' USHIST : PAS DROIT A PLUS DE ',NUSHMX,' FICHIERS HISTORIQUES'
          CALL CSEXIT (1)
        ENDIF
C
        DO II = 1, NVARPP
C
C   --> Ouverture des fichiers avec les unites disponibles
C
          IF (IRANGP.LE.0) THEN
            OPEN(FILE=FICUSH(II),UNIT=IMPUSH(II))
          ENDIF
C
C   --> On imprime le numero (global) de la cellule et les coordonnees
C        du centre
C
          DO KK = 1, NCAPTS
C           Numero de cellule (en parallele : local au processeur courant)
            NUMCEL    = ICAPT(KK)
            IF (IRANGP.LT.0 .OR. IRANGP.EQ.IRCAPT(KK)) THEN
C           Coordonnees de la cellule (en parallele, c'est le processeur
C             qui la contient qui travaille)
              XYZTMP(1) = XYZCEN(1,NUMCEL)
              XYZTMP(2) = XYZCEN(2,NUMCEL)
              XYZTMP(3) = XYZCEN(3,NUMCEL)
            ELSE
C           Valeurs bidons sur les autres processeurs
              XYZTMP(1) = 0.D0
              XYZTMP(2) = 0.D0
              XYZTMP(3) = 0.D0
            ENDIF
C           En parallele, le processeur qui a trouve la cellule
C           envoie son numero global et ses coordonnees aux autres.
            IF (IRANGP.GE.0) THEN
              CALL PARCEL(ICAPT(KK), IRCAPT(KK), NUMCEL)
C             ===========
              LNG = 3
              CALL PARBCR(IRCAPT(KK), LNG , XYZTMP)
C             ===========
            ENDIF
C           On ecrit les informations (seul le processeur 0
C           travaille en parallele : on n'a pasa besoin de
C           plusieurs exemplaires du fichier)
            IF (IRANGP.LE.0) THEN
              WRITE(IMPUSH(II),1000) '#',' Cellule ',NUMCEL,
     &            ' Coord ',XYZTMP(1),XYZTMP(2),XYZTMP(3)
            ENDIF
C
          ENDDO
C
        ENDDO
C
      ENDIF
C
 1000 FORMAT(A,A9,I10,A7,3E14.5)
C
C=======================================================================
C 4. ECRITURE
C     EXEMPLE D'UNE VARIABLE PAR FICHIER
C=======================================================================
C
C Ecriture du numero du pas de temps,
C          de la valeur du temps physique
C          de la variable en tous les points d'historiques
C En sequentiel, la valeur a ecrire est simplement RTP(ICAPT(KK),II)
C en parallele, la valeur a ecrire peut etre sur un autre processeur
C  et il faut la determiner dans  VACAPT(KK) avec PARHIS.
C
      DO II = 1 , NVARPP
        DO KK = 1, NCAPTS
          IF (IRANGP.LT.0) THEN
            VACAPT(KK) = RTP(ICAPT(KK),II)
          ELSE
            CALL PARHIS(ICAPT(KK), IRCAPT(KK), RTP(1,II), VACAPT(KK))
C           ===========
          ENDIF
        ENDDO
        IF (IRANGP.LE.0) THEN
          WRITE (IMPUSH(II),1010) NTCABS,TTCABS,
     &                            (VACAPT(KK),KK=1,NCAPTS)
        ENDIF
      ENDDO
C
C
C
C   ATTENTION : IL FAUT ADAPTER LE FORMAT POUR PLUS DE  9 CAPTEURS
C
C
 1010 FORMAT(I10,10E17.9)
C
C=======================================================================
C 4. FERMETURE
C=======================================================================
C
      IF(NTCABS.EQ.NTMABS .AND. IRANGP.LE.0) THEN
        DO II = 1, NVARPP
          CLOSE(IMPUSH(II))
        ENDDO
      ENDIF
C
C
C=======================================================================
C 5. SORTIE
C=======================================================================
C
      RETURN
      END
c@z
