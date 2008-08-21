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
                        SUBROUTINE STRDEP
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , ITRALE , ITALIM , ITRFIN ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NVAR   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   FLMALF , FLMALB , COFALE , XPRALE , RDEVEL , RTUSER ,
     &   RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  DEPLACEMENT DES STRUCTURES MOBILES EN ALE EN COUPLAGE INTERNE
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
CARGU ! ITRFIN       ! E  !  -> ! INDICATEUR DE DERNIERE ITERATION DE  !
CARGU !              !    !     !                    COUPLAGE IMPLICITE!
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
      INCLUDE "cstphy.h"
      INCLUDE "numvar.h"
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
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          ITRALE , ITALIM , ITRFIN
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR , NVAR
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
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION FLMALF(NFAC), FLMALB(NFABOR), XPRALE(NCELET)
      DOUBLE PRECISION COFALE(NFABOR,8)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          ISTR, II, IEL, IFAC, IPHAS, NTAB
      INTEGER          IFLMAS, IFLMAB, ICLP, ICLU, ICLV, ICLW
      DOUBLE PRECISION DELTA
C
C***********************************************************************
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
C 2.  CALCUL DES EFFORTS SUR LA STRUCTURE
C=======================================================================
C
      DO ISTR = 1, NBSTRU
        DO II = 1, NDIM
          FORSTA(II,ISTR) = FORSTR(II,ISTR)
c-a tester          FORSTA(II,ISTR) = FORSTP(II,ISTR)
          FORSTR(II,ISTR) = 0.D0
        ENDDO
      ENDDO
C
      DO IFAC = 1, NFABOR
        ISTR = IA(IIDFST+IFAC-1)
        IF (ISTR.GT.0) THEN
          DO II = 1, 3
            FORSTR(II,ISTR) = FORSTR(II,ISTR)
     &           + RA(IFORBR+(IFAC-1)*NDIM+II-1)
          ENDDO
        ENDIF
      ENDDO
C
      IF (IRANGP.GE.0) THEN
        NTAB = NDIM*NBSTRU
        CALL PARRSM(NTAB,FORSTR)
      ENDIF
C
      DO ISTR = 1, NBSTRU
        DO II = 1, NDIM
          FORSTP(II,ISTR) = CFOPRE*FORSTR(II,ISTR)+
     &         (1.D0-CFOPRE)*FORSTA(II,ISTR)
        ENDDO
      ENDDO
C
C     Si on est en phase d'initialisation du fluide
      IF (ITRALE.LE.NALINF) THEN
        ITRFIN = -1
        RETURN
      ENDIF
C
C=======================================================================
C 3.  DEFINITION UTILISATEUR DES CARACTERISTIQUES DE LA STRUCTURE
C=======================================================================
C
C
      CALL USSTR2
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NBSTRU ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , RA(IIDFST),
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     ,
     &   XMSTRU , XCSTRU , XKSTRU , XSTREQ , XSTR   , XPSTR  , FORSTP ,
     &   DTSTR  ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C
C=======================================================================
C 4.  DEPLACEMENT DES STRUCTURES
C=======================================================================
C
C
      DO ISTR = 1, NBSTRU
C
        CALL NEWMRK
C       ===========
     & ( ISTR  , ALPNMK  , BETNMK          , GAMNMK          ,
     &   XMSTRU(1,1,ISTR), XCSTRU(1,1,ISTR), XKSTRU(1,1,ISTR),
     &   XSTREQ(1,ISTR)  ,
     &   XSTR(1,ISTR)    , XPSTR(1,ISTR)   , XPPSTR(1,ISTR)  ,
     &   XSTA(1,ISTR)    , XPSTA(1,ISTR)   , XPPSTA(1,ISTR)  ,
     &   FORSTP(1,ISTR)  , FORSTA(1,ISTR)  , DTSTR(ISTR)     )
C
      ENDDO
C
C=======================================================================
C 5.  TEST DE CONVERGENCE
C=======================================================================
C
      IPHAS = 1
      DELTA = 0.D0
      DO ISTR = 1, NBSTRU
        DO II = 1, 3
          DELTA = DELTA + (XSTR(II,ISTR)-XSTP(II,ISTR))**2
        ENDDO
      ENDDO
      DELTA = SQRT(DELTA)/ALMAX(IPHAS)/NBSTRU
C
      IF (IWARNI(IUMA).GE.2) WRITE(NFECRA,1000) ITALIM, DELTA
C
C     si convergence
      IF (DELTA.LT.EPALIM) THEN
        IF (ITRFIN.EQ.1) THEN
C       si ITRFIN=1 on sort
          IF (IWARNI(IUMA).GE.1) WRITE(NFECRA,1001) ITALIM, DELTA
          ITRFIN = -1
        ELSE
C       sinon on refait une derniere iteration pour Syrthes/T1D/rayonnement
          ITRFIN = 1
        ENDIF
      ELSEIF (ITRFIN.EQ.0 .AND. ITALIM.EQ.NALIMX-1) THEN
C       ce sera la derniere iteration
        ITRFIN = 1
      ELSEIF (ITALIM.EQ.NALIMX) THEN
C       on a forcement ITRFIN=1 et on sort
        IF (NALIMX.GT.1) WRITE(NFECRA,1100) ITALIM, DELTA
        ITRFIN = -1
      ENDIF
C
C=======================================================================
C 6.  RETOUR AUX VALEURS ANTERIEURES SI NECESSAIRE
C=======================================================================
C
C     Si NTERUP.GT.1, RTPA a ete touche apres NAVSTO, on doit donc
C       revenir a une valeur anterieure
      IF (ITRFIN.NE.-1) THEN
        DO II = 1, NVAR
          IF (II.EQ.IPR(1) .AND. NTERUP.GT.1) THEN
            DO IEL = 1, NCELET
              RTPA(IEL,II) = XPRALE(IEL)
            ENDDO
          ENDIF
          DO IEL = 1, NCELET
            RTP(IEL,II) = RTPA(IEL,II)
          ENDDO
        ENDDO
        DO IFAC = 1, NFAC
           PROPFA(IFAC,IFLMAS) = FLMALF(IFAC)
        ENDDO
        DO IFAC = 1, NFABOR
           PROPFB(IFAC,IFLMAB) = FLMALB(IFAC)
           COEFA(IFAC,ICLP) = COFALE(IFAC,1)
           COEFA(IFAC,ICLU) = COFALE(IFAC,2)
           COEFA(IFAC,ICLV) = COFALE(IFAC,3)
           COEFA(IFAC,ICLW) = COFALE(IFAC,4)
           COEFB(IFAC,ICLP) = COFALE(IFAC,5)
           COEFB(IFAC,ICLU) = COFALE(IFAC,6)
           COEFB(IFAC,ICLV) = COFALE(IFAC,7)
           COEFB(IFAC,ICLW) = COFALE(IFAC,8)
        ENDDO
      ENDIF
C
C----
C FORMATS
C----
C
#if defined(_CS_LANG_FR)
C
 1000 FORMAT (
     & '            ALE IMPLICITE : ITER=',I5,' DERIVE=',E12.5     )
 1001 FORMAT (
     & 'CONVERGENCE ALE IMPLICITE : ITER=',I5,' DERIVE=',E12.5     )
 1100 FORMAT (
     &'@                                                            ',/,
     &'@ @@ ATTENTION : COUPLAGE IMPLICITE ALE                      ',/,
     &'@    *********                                               ',/,
     &'@  Nombre d''iterations maximal ',I10   ,' atteint           ',/,
     &'@  Derive normee :',E12.5                                     ,/,
     &'@                                                            '  )
C
#else
C
 1000 FORMAT (
     & '            IMPLICIT ALE: ITER=',I5,' DERIVE=',E12.5     )
 1001 FORMAT (
     & 'CONVERGENCE IMPLICIT ALE: ITER=',I5,' DERIVE=',E12.5     )
 1100 FORMAT (
     &'@                                                            ',/,
     &'@ @@ WARNING: IMPLICIT ALE                                   ',/,
     &'@    ********                                                ',/,
     &'@  Maximum number of iterations ',I10   ,' reached           ',/,
     &'@  Normed derive :',E12.5                                     ,/,
     &'@                                                            '  )
C
#endif
C
C----
C FIN
C----
C
      END
c@z
