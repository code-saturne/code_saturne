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
                        SUBROUTINE STRINI
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RDEVEL , RTUSER ,
     &   RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  INITILISATION DES DONNEES DES STRUCTURES MOBILES EN ALE
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
CARGU ! IFINIA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU ! IFINRA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              !    !     !  DANS IA EN SORTIE                   !
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
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
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
      INCLUDE "paramx.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "pointe.h"
      INCLUDE "albase.h"
      INCLUDE "alstru.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0 , IFINIA , IFINRA
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
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC  , ISTR, ICOMPT, II
      INTEGER          MBSTRU
      INTEGER          MAXELT, IFNIA2, ILS
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
      DO ISTR = 1, NSTRMX
        DTSTR(ISTR) = DT(1)
      ENDDO
C
C     NBSTRU vaut -999 si le calcul n'est pas un calcul suite ou s'il
C       est une suite d'un calcul sans ALE
      MBSTRU = NBSTRU
C
C     En ALE on met IHISTR a 1 par defaut (remis a zero sans structure)
      IHISTR = 1
C
C=======================================================================
C 2.  RESERVATION DU TABLEAU IDFSTR
C=======================================================================
C
      IIDFST = IDEBIA
      IFINIA = IIDFST + NFABOR
C
      MAXELT = MAX(NCELET,NFAC,NFABOR)
      ILS    = IFINIA
      IFNIA2 = ILS + MAXELT
C
      CALL IASIZE('STRINI',IFNIA2)
C     ===========
C
      DO IFAC = 1, NFABOR
        IA(IIDFST+IFAC-1) = 0
      ENDDO
C
C
C=======================================================================
C 2.  REMPLISSAGE DE IDFSTR PAR L'UTILISATEUR
C=======================================================================
C
C
      CALL USSTR1
C     ===========
     & ( IFNIA2 , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IIDFST),
     &   IDEVEL , ITUSER , IA     ,
     &   AEXXST , BEXXST , CFOPRE ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   XSTP   , XPSTR  , XSTREQ ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 3.  CALCUL DE NBSTRU
C=======================================================================
C
      NBSTRU = 0
      DO IFAC = 1, NFABOR
        IF (IA(IIDFST+IFAC-1).GT.NBSTRU) NBSTRU = IA(IIDFST+IFAC-1)
      ENDDO
C
      IF (IRANGP.GE.0) CALL PARCMX(NBSTRU)
C                      ===========
C     On compare NBSTRU a la valeur eventuelle anterieure
C
      IF (MBSTRU.GT.-999) THEN
        IF (NBSTRU.NE.MBSTRU) THEN
          WRITE(NFECRA,4000)MBSTRU,NBSTRU
          CALL CSEXIT(1)
        ENDIF
      ENDIF
C
C     Valeur par defaut et verifiction de IHISTR
      IF (NBSTRU.EQ.0) IHISTR = 0
C
      ICOMPT = 0
      DO II = 2, NVPPMX
        IF(IHISVR(II,1).NE.0) ICOMPT = ICOMPT+1
      ENDDO
      IF( (ICOMPT.EQ.0.OR.NCAPT.EQ.0) .AND. IHISTR.EQ.0 ) THEN
        NTHIST = -1
      ENDIF
C
      IF (IHISTR.NE.0 .AND. IHISTR.NE.1) THEN
        WRITE(NFECRA,1000)IHISTR
        CALL CSEXIT(1)
      ENDIF
C
C     Si NBSTRU=0, on desalloue IDFSTR et on passe NALIMX a 1 si necessaire
      IF (NBSTRU.EQ.0) THEN
        IFINIA = IDEBIA
        WRITE(NFECRA,2000) NBSTRU
        IF (NALIMX.GT.1) THEN
          WRITE(NFECRA,2001)
          NALIMX = 1
        ENDIF
      ELSE
        WRITE(NFECRA,2010) NBSTRU,ALPNMK,BETNMK,GAMNMK,IHISTR
        IF (NALIMX.EQ.1) THEN
          WRITE(NFECRA,2020) AEXXST, BEXXST, CFOPRE
        ELSE
          CFOPRE = 1.D0
          WRITE(NFECRA,2030) NALIMX, EPALIM
        ENDIF
      ENDIF
      WRITE(NFECRA,3000)
C
C----
C FORMATS
C----
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR D''ECRITURE DES FICHIERS HISTORIQUES DES  ',/,
     &'@      STRUCTURES MOBILES NE PEU VALOIR QUE 0 OU 1.          ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes dans usstru.               ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2000 FORMAT(
     &    /,'TTES PHASES  : MODE COUPLAGE DE STRUCTURES NON ACTIVE  ',/,
     &      '                 NBSTRU = ',I10                         ,/)
 2001 FORMAT(
     &      '               NALIMX INUTILE ET POSITIONNE A 1        ',/)
 2010 FORMAT(
     &    /,'TTES PHASES  : MODE COUPLAGE DE STRUCTURES ACTIVE      ',/,
     &      '                 AVEC NBSTRU = ',I10   ,' STRUCTURE(S) ',/,
     &      '                                                       ',/,
     &      '               COEFFICIENTS DE NEWMARK :               ',/,
     &      '                 ALPNMK = ',E12.4                       ,/,
     &      '                 BETNMK = ',E12.4                       ,/,
     &      '                 GAMNMK = ',E12.4                       ,/,
     &      '                                                       ',/,
     &      '               FICHIERS HISTORIQUES DES STRUCTURES :   ',/,
     &      '                 IHISTR = ',I4,' ( 1 : active)         ',/)
 2020 FORMAT(
     &    /,'TTES PHASES  : SCHEMA DE COUPLAGE EXPLICTE ACTIVE      ',/,
     &      '                                                       ',/,
     &      '               COEFFICIENTS DU SCHEMA :                ',/,
     &      '                 AEXXST = ',E12.4                       ,/,
     &      '                 BEXXST = ',E12.4                       ,/,
     &      '                 CFOPRE = ',E12.4                       ,/)
 2030 FORMAT(
     &    /,'TTES PHASES  : SCHEMA DE COUPLAGE IMPLICITE ACTIVE     ',/,
     &      '                                                       ',/,
     &      '               NB DE SOUS-ITERATIONS MAX. : ',I10       ,/,
     &      '               SEUIL DE CONVERGENCE       : ',E12.4     ,/)
C
 3000 FORMAT(
     &'-------------------------------------------------------------',/)
 4000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA DEFINITION DES STRUCTURES MOBILES',/,
     &'@                                               (METHODE ALE)',/,
     &'@                                                            ',/,
     &'@    Le nombre de structures definies est different de celui ',/,
     &'@      du calcul precedent :                                 ',/,
     &'@      Nombre de structures calcul precedent : ',I10          ,/,
     &'@      Nombre de structures calcul actuel    : ',I10          ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas  execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier le fichier suite auxiliaire ou la specification',/,
     &'@      des structures dans usstru.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      END
c@z
