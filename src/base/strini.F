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
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "pointe.h"
      INCLUDE "albase.h"
      INCLUDE "alstru.h"
      INCLUDE "alaste.h"
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
      INTEGER          IDEBIA, IDEBRA, ITRAV
      INTEGER          IFAC  , ISTR, ICOMPT, II
      INTEGER          MBSTRU, MBASTE
      INTEGER          MAXELT, IFNIA2, ILS
C
      INTEGER          JJ, INOD
      INTEGER          ILSTFA, INDAST, IIDFLO, IIDNLO
      INTEGER          NBPDT, NBSSIT, IHI, CHRO
C
      DOUBLE PRECISION DELTA, PDT, TT
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
      DO ISTR = 1, NASTMX
        DO II = 1, 3
          ASDDLF(II,ISTR) = 1
          ASDDLC(II,ISTR) = 1
        ENDDO
      ENDDO
C
C     NBSTRU et NBASTE valent -999 si le calcul n'est pas un calcul
C       suite ou s'il est une suite d'un calcul sans ALE
      MBSTRU = NBSTRU
      MBASTE = NBASTE
C
C     En ALE on met IHISTR et ISYNCP a 1 par defaut
C       (remis a zero sans structure)
      IHISTR = 1
      ISYNCP = 1
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
C 2.1 STRUCTURES INTERNES :
C -----------------------
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
C 2.2 STRUCTURES EXTERNES : COUPLAGE CODE_SATURNE / CODE_ASTER
C -----------------------
C
      CALL USASTE
C     ===========
     & ( IFNIA2 , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IIDFST),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RDEVEL , RTUSER , RA     )
C
C
C=======================================================================
C 3.  CALCUL DE NBSTRU ET NBASTE
C=======================================================================
C
C 3.1 STRUCTURES INTERNES :
C -----------------------
C
      NBSTRU = 0
      DO IFAC = 1, NFABOR
        IF (IA(IIDFST+IFAC-1).GT.NBSTRU) NBSTRU = IA(IIDFST+IFAC-1)
      ENDDO
C
      IF (IRANGP.GE.0) CALL PARCMX(NBSTRU)
C                      ===========
C
      IF (NBSTRU.GT.NSTRMX) THEN
        WRITE(NFECRA,4000)
        CALL CSEXIT(1)
      ENDIF
C
C     On compare NBSTRU a la valeur eventuelle anterieure
C
      IF (MBSTRU.GT.-999) THEN
        IF (NBSTRU.NE.MBSTRU) THEN
          WRITE(NFECRA,4001)MBSTRU,NBSTRU
          CALL CSEXIT(1)
        ENDIF
      ENDIF
C
C 3.2 STRUCTURES EXTERNES : COUPLAGE CODE_SATURNE / CODE_ASTER
C -----------------------
C
      NBASTE = 0
      DO IFAC = 1, NFABOR
        IF (-IA(IIDFST+IFAC-1).GT.NBASTE) NBASTE = -IA(IIDFST+IFAC-1)
      ENDDO
C
      IF (IRANGP.GE.0) CALL PARCMX(NBASTE)
C                      ===========
C
      IF (NBASTE.GT.NASTMX) THEN
        WRITE(NFECRA,4002)
        CALL CSEXIT(1)
      ENDIF
C
C     On compare NBASTE a la valeur eventuelle anterieure
C
      IF (MBASTE.GT.-999) THEN
        IF (NBASTE.NE.MBASTE) THEN
          WRITE(NFECRA,4003)MBASTE,NBASTE
          CALL CSEXIT(1)
        ENDIF
      ENDIF
C
C
C=======================================================================
C 5.  CALCUL ET ENVOI A CODE_ASTER DES PARAMETRES GEOMETRIQUES
C=======================================================================
C
      IF (NBASTE.GT.0) THEN
C
        ITRAV  = IFINIA
        IFNIA2 = ITRAV + NNOD
        CALL IASIZE('STRINI',IFNIA2)
C
        DO INOD = 1, NNOD
           IA(ITRAV + INOD-1) = 0
        ENDDO
C
        NBFAST = 0
        NBNAST = 0
C
C       Calcul du nombre de faces et noeuds couples avec Code_Aster
        DO IFAC = 1, NFABOR
          ISTR = IA(IIDFST+IFAC-1)
          IF (ISTR.LT.0) THEN
            NBFAST = NBFAST + 1
            DO II = IPNFBR(IFAC), IPNFBR(IFAC+1)-1
              INOD = NODFBR(II)
              IA(ITRAV + INOD-1) = ISTR
            ENDDO
          ENDIF
        ENDDO
        DO INOD = 1, NNOD
          IF (IA(ITRAV+INOD-1).GT.0) NBNAST = NBNAST + 1
        ENDDO
C
        ILSTFA = IFINIA
        IIDFLO = ILSTFA + NBFAST
        IIDNLO = IIDFLO + NBFAST
        IFNIA2 = IIDNLO + NBNAST
        CALL IASIZE('STRINI',IFNIA2)
C
        INDAST = 0
        DO IFAC = 1, NFABOR
          ISTR = IA(IIDFST+IFAC-1)
          IF (ISTR.LT.0) THEN
            INDAST = INDAST + 1
            IA(ILSTFA + INDAST-1) = IFAC
            IA(IIDFLO + INDAST-1) = -ISTR
          ENDIF
        ENDDO
        NBFAST = INDAST
C
        INDAST = 0
        DO INOD = 1, NNOD
          ISTR = IA(ITRAV+INOD-1)
          IF (ISTR.LT.0) THEN
            INDAST = INDAST + 1
            IA(IIDNLO + INDAST-1) = -ISTR
          ENDIF
        ENDDO
        NBNAST = INDAST
C
C       Envoi des donnees geometriques a Code_Aster
        CALL ASTGEO
C       ===========
     & ( NBFAST, NBNAST, IA(ILSTFA), IA(IIDFLO), IA(IIDNLO), ALMAX )
C
C       Recuperation des parametres commun du couplage
        CALL ASTPAR
C       ===========
     & ( NTMABS, NALIMX, EPALIM, ISYNCP, NTCHR, TTPABS, DTREF )
C
      ENDIF
C
C
C=======================================================================
C 6.  MESSAGES D'INFORMATION SUR LE COUPLAGE
C=======================================================================
C
C     Valeur par defaut et verifiction de IHISTR
      IF (NBSTRU.EQ.0) IHISTR = 0
      IF (NBASTE.EQ.0) ISYNCP = 0
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
      IF (ISYNCP.NE.0 .AND. ISYNCP.NE.1) THEN
        WRITE(NFECRA,1002)ISYNCP
        CALL CSEXIT(1)
      ENDIF
C
C     Si NBSTRU=0 et NBASTE=0, on desalloue IDFSTR et on passe NALIMX a 1
C       si necessaire
      IF (NBSTRU.GT.0) THEN
        WRITE(NFECRA,2010) NBSTRU,ALPNMK,BETNMK,GAMNMK,IHISTR
      ELSE
        WRITE(NFECRA,2000) NBSTRU
      ENDIF
      IF (NBASTE.GT.0) THEN
        WRITE(NFECRA,2012) NBASTE,ISYNCP
      ELSE
        WRITE(NFECRA,2002) NBASTE
      ENDIF
      IF (NBSTRU.EQ.0.AND.NBASTE.EQ.0) THEN
        IFINIA = IDEBIA
        IFINRA = IDEBRA
        IF (NALIMX.GT.1) THEN
          WRITE(NFECRA,2001)
          NALIMX = 1
        ENDIF
      ELSE IF (NBSTRU.GT.0) THEN
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
#if defined(_CS_LANG_FR)
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR D''ECRITURE DES FICHIERS HISTORIQUES DES  ',/,
     &'@      STRUCTURES MOBILES NE PEUT VALOIR QUE 0 OU 1.         ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les parametres donnes dans usstru.               ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR D''ECRITURE AU MEME INSTANT DES SORTIES   ',/,
     &'@      CODE_SATURNE ET CODE_ASTER NE PEUT VALOIR QUE 0 OU 1. ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@    LA PERIODE D''IMPRESSION EST DEFINIE PAR LA VARIABLE :  ',/,
     &'@      NTCHR                                                 ',/,
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
 2002 FORMAT(
     &    /,'TTES PHASES  : MODE COUPLAGE CODE_ASTER NON ACTIVE     ',/,
     &      '                 NBASTE = ',I10                         ,/)
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
 2012 FORMAT(
     &    /,'TTES PHASES  : MODE COUPLAGE CODE_ASTER ACTIVE         ',/,
     &      '                 AVEC NBASTE = ',I10   ,' STRUCTURE(S) ',/,
     &      '                                                       ',/,
     &      '               IMPRESSIONS  DES SORTIES AUX MEMES      ',/,
     &      '               INSTANTS DES DEUX CODES :               ',/,
     &      '                 ISYNCP = ',I4,' ( 1 : active)         ',/)
 2020 FORMAT(
     &    /,'TTES PHASES  : SCHEMA DE COUPLAGE EXPLICITE ACTIVE     ',/,
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
     &'@ @@ ATTENTION : ARRET A LA DEFINITION DES STRUCTURES        ',/,
     &'@                INTERNES                                    ',/,
     &'@                                                            ',/,
     &'@    Le nombre de structures definies est superieur au nombre',/,
     &'@    maximum autorise NSTRMX :                               ',/,
     &'@      Nombre de structures definies      : ',I10             ,/,
     &'@      Nombre de structures autorisees    : ',I10             ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas  execute.                         ',/,
     &'@                                                            ',/,
     &'@    Reduire le nombre de structure                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA DEFINITION DES STRUCTURES MOBILES',/,
     &'@                INTERNES                                    ',/,
     &'@                                                            ',/,
     &'@    Le nombre de structures internes definies est           ',/,
     &'@      different de celui du calcul precedent :              ',/,
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
 4002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA DEFINITION DES STRUCTURES        ',/,
     &'@                              EXTERNES (COUPLAGE CODE_ASTER)',/,
     &'@                                                            ',/,
     &'@    Le nombre de structures definies est superieur au nombre',/,
     &'@    maximum autorise NASTMX :                               ',/,
     &'@      Nombre de structures definies      : ',I10             ,/,
     &'@      Nombre de structures autorisees    : ',I10             ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas  execute.                         ',/,
     &'@                                                            ',/,
     &'@    Reduire le nombre de structure                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4003 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA DEFINITION DES STRUCTURES MOBILES',/,
     &'@                             EXTERNES  (COUPLAGE CODE_ASTER)',/,
     &'@                                                            ',/,
     &'@    Le nombre de structures externes definies est           ',/,
     &'@      different de celui du calcul precedent :              ',/,
     &'@      Nombre de structures calcul precedent : ',I10          ,/,
     &'@      Nombre de structures calcul actuel    : ',I10          ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas  execute.                         ',/,
     &'@                                                            ',/,
     &'@    Verifier le fichier suite auxiliaire ou la specification',/,
     &'@      des structures dans usaste.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#else
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    THE TIME MONITORING FILES INDICATOR FOR THE MOBILE      ',/,
     &'@      STRUCTURES CAN ONLY TAKE THE VALUES 0 OR 1.           ',/,
     &'@    ITS VALUE IS ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not run.                             ',/,
     &'@                                                            ',/,
     &'@  Verify the parameters given in usstru.                    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,
     &'@    ********                                                ',/,
     &'@    THE INDICATOR OF SYNCHRONIZED OUTPUTS FOR CODE_SATURNE  ',/,
     &'@      AND CODE_ASTER CAN ONLY TAKE THE VALUES 0 OR 1.       ',/,
     &'@    ITS VALUE IS ',I10                                       ,/,
     &'@                                                            ',/,
     &'@  The calculation will not run.                             ',/,
     &'@                                                            ',/,
     &'@  Verify the parameters given in usstru.                    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2000 FORMAT(
     &    /,'ALL PHASES: COUPLING MODE FOR STRUCTURES NOT ACTIVATED ',/,
     &      '              NBSTRU = ',I10                            ,/)
 2001 FORMAT(
     &      '            NALIMX USELESS AND SET TO 1                ',/)
 2002 FORMAT(
     &    /,'ALL PHASES: CODE_ASTER COUPLING MODE NOT ACTIVATED     ',/,
     &      '              NBASTE = ',I10                            ,/)
 2010 FORMAT(
     &    /,'ALL PHASES: COUPLING MODE FOR STRUCTURES ACTIVATED     ',/,
     &      '              WITH NBSTRU = ',I10   ,' STRUCTURE(S)    ',/,
     &      '                                                       ',/,
     &      '            NEWMARK COEFFICIENTS:                      ',/,
     &      '              ALPNMK = ',E12.4                          ,/,
     &      '              BETNMK = ',E12.4                          ,/,
     &      '              GAMNMK = ',E12.4                          ,/,
     &      '                                                       ',/,
     &      '            MONITORING FILES FOR STRUCTURES:           ',/,
     &      '                 IHISTR = ',I4,' ( 1 : activated)      ',/)
 2012 FORMAT(
     &    /,'ALL PHASES: CPDE_ASTER COUPLING MODE ACTIVATED         ',/,
     &      '              WITH NBASTE = ',I10   ,' STRUCTURE(S)    ',/,
     &      '                                                       ',/,
     &      '            SYNCHRONIZED OUTPUTS FOR BOTH CODES:       ',/,
     &      '                 ISYNCP = ',I4,' ( 1 : activated)      ',/)
 2020 FORMAT(
     &    /,'ALL PHASES: EXPLICIT SCHEME FOR COUPLING ACTIVATED     ',/,
     &      '                                                       ',/,
     &      '            SCHEME COEFFICIENTS:                       ',/,
     &      '              AEXXST = ',E12.4                          ,/,
     &      '              BEXXST = ',E12.4                          ,/,
     &      '              CFOPRE = ',E12.4                          ,/)
 2030 FORMAT(
     &    /,'ALL PHASES: IMPLICIT SCHEME FOR COUPING ACTIVATED      ',/,
     &      '                                                       ',/,
     &      '            NB OF MAX INNER ITERATIONS : ',I10          ,/,
     &      '            CONVERGENCE THRESHOLD      : ',E12.4        ,/)
C
 3000 FORMAT(
     &'-------------------------------------------------------------',/)
 4000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE INTERNAL STRUCTURES SPECIFICATION ',/,
     &'@                                                            ',/,
     &'@    The number of defined structures is greater than the    ',/,
     &'@      allowed maximum NSTRMX:                               ',/,
     &'@      Number of defined structures: ',I10                    ,/,
     &'@      Number of allowed structures: ',I10                    ,/,
     &'@                                                            ',/,
     &'@    The calculation will not be run.                        ',/,
     &'@                                                            ',/,
     &'@    Decrease the number of structures                       ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE INTERNAL MOBILE STRUCTURES        ',/,
     &'@             SPECIFICATION                                  ',/,
     &'@                                                            ',/,
     &'@    The number of defined structures is different from the  ',/,
     &'@      previous calculation:                                 ',/,
     &'@      Number of structures previous calculation: ',I10       ,/,
     &'@      Number of structures current  calculation: ',I10       ,/,
     &'@                                                            ',/,
     &'@    The calculation will not be run.                        ',/,
     &'@                                                            ',/,
     &'@    Verify the auxilliary restart file or the structures    ',/,
     &'@      specifications in usstru.                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE EXTERNAL MOBILE STRUCTURES        ',/,
     &'@             SPECIFICATION (CODE_ASTER COUPLING)            ',/,
     &'@                                                            ',/,
     &'@    The number of defined structures is greater than the    ',/,
     &'@      allowed maximum NASTMX:                               ',/,
     &'@      Number of defined structures: ',I10                    ,/,
     &'@      Number of allowed structures: ',I10                    ,/,
     &'@                                                            ',/,
     &'@    The calculation will not be run.                        ',/,
     &'@                                                            ',/,
     &'@    Decrease the number of structures                       ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4003 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN THE EXTERNAL MOBILE STRUCTURES        ',/,
     &'@             SPECIFICATION (CODE_ASTER COUPLING)            ',/,
     &'@                                                            ',/,
     &'@    The number of defined structures is different from the  ',/,
     &'@      previous calculation:                                 ',/,
     &'@      Number of structures previous calculation: ',I10       ,/,
     &'@      Number of structures current  calculation: ',I10       ,/,
     &'@                                                            ',/,
     &'@    The calculation will not be run.                        ',/,
     &'@                                                            ',/,
     &'@    Verify the auxilliary restart file or the structures    ',/,
     &'@      specifications in usstru.                             ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#endif
C
C----
C FIN
C----
C
      END
c@z
