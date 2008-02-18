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
                        SUBROUTINE ELFLUX
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , IAPPEL ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , VISCF  , VISCB  ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  ROUTINE PHYSIQUE PARTICULIERE POUR LE MODULE ELECTRIQUE
CFONC
CFONC               CALCULS DES VARIABLES PROPCE
CFONC         ELLES SONT UTILISEES POUR LE CALCUL DES TS
CFONC
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! IAPPEL       ! E  !  -> ! NUMERO D'APPEL                       !
CARGU !              !    !     ! 1 : j, E, j.E                        !
CARGU !              !    !     ! 2 : B,    jxB                        !
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
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! ITYPSM       ! TE !  -> ! TYPE DE SOURCE DE MASSE POUR LES     !
CARGU ! (NCESMP,NVAR)!    !     !  VARIABLES (cf. USTSMA)              !
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
CARGU ! SMACEL       ! TR !  -> ! VALEUR DES VARIABLES ASSOCIEE A LA   !
CARGU ! (NCESMP,*   )!    !     !  SOURCE DE MASSE                     !
CARGU !              !    !     !  POUR IVAR=IPR, SMACEL=FLUX DE MASSE !
CARGU ! VISCF(NFAC)  ! TR !  -  ! TABLEAU DE TRAVAIL    FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! TABLEAU DE TRAVAIL    FACES DE BORD  !
CARGU ! W1..9(NCELET ! TR !  -  ! TABLEAU DE TRAVAIL    CELLULES       !
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
C-----------------------------------------------------------------------
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
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "elincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0 , IAPPEL
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
      INTEGER          IDEVEL(NIDEVE)
      INTEGER          ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NDIMFB,*)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IEL
      INTEGER          IPCEFJ, IPCSIG, IPCSII
      INTEGER          IPCLA1, IPCLA2, IPCLA3
      INTEGER          IPCDC1, IPCDC2, IPCDC3
      INTEGER          IPCDI1, IPCDI2, IPCDI3
      INTEGER          INC   , ICCOCG, NSWRGP, IMLIGP, IWARNP
      INTEGER          IVAR0 , ICLIMV
      INTEGER          IPHYDP, IDIMTE, ITENSO, IVAR  , MODNTL
C
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP, VRMIN, VRMAX, VAR
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
C --- Numero des grandeurs physiques
      IPCSIG = IPPROC(IVISLS(IPOTR))
      IPCSII = IPPROC(IVISLS(IPOTI))
      IPCEFJ = IPPROC(IEFJOU)
      IPCDC1 = IPPROC(IDJR(1))
      IPCDC2 = IPPROC(IDJR(2))
      IPCDC3 = IPPROC(IDJR(3))
C
      IF ( IPPMOD(IELJOU).EQ.4 ) THEN
        IPCDI1 = IPPROC(IDJI(1))
        IPCDI2 = IPPROC(IDJI(2))
        IPCDI3 = IPPROC(IDJI(3))
      ENDIF
C
C --- Necessite d'une impression (impression si MODNTL=0)
      IF(NTLIST.GT.0) THEN
        MODNTL = MOD(NTCABS,NTLIST)
      ELSEIF(NTLIST.EQ.-1.AND.NTCABS.EQ.NTMABS) THEN
        MODNTL = 0
      ELSE
        MODNTL = 1
      ENDIF
C
C=======================================================================
C 1. PREMIER APPEL : J, E => J.E
C=======================================================================
C
      IF(IAPPEL.EQ.1) THEN
C
C=======================================================================
C 1.1 PRISE EN COMPTE DES TERMES SOURCES ET VARIABLES STANDARD ET
C   COMMUNES A TOUTES LES VERSIONS ELECTRIQUES : EFFET JOULE, E, J
C         CAS IELJOU >= 1 ou IELARC >= 1 ou IELION >= 1
C=======================================================================
C
C    Pour toutes les versions electriques :
C      on doit calculer le terme joule (partie reelle) = j . E
C      produit de la densite de courant par le champ electrique
C
C
C   2.1 Calcul du grad (potR) (W4, W5, W6)
C  ---------------------------
C
        IVAR = ISCA(IPOTR)
        ICLIMV = ICLRTP(IVAR,ICOEF)
C
        INC = 1
        ICCOCG = 1
        NSWRGP = NSWRGR(IVAR)
        IMLIGP = IMLIGR(IVAR)
        IWARNP = IWARNI(IVAR)
        EPSRGP = EPSRGR(IVAR)
        CLIMGP = CLIMGR(IVAR)
        EXTRAP = EXTRAG(IVAR)
C
C En periodique et parallele, echange avant calcul du gradient
C     C'est indispensable car on vient de calculer IVAR
C
C    Parallele
        IF(IRANGP.GE.0) THEN
          CALL PARCOM(RTP(1,IVAR))
C         ===========
        ENDIF
C
C    Periodique
        IF(IPERIO.EQ.1) THEN
          IDIMTE = 0
          ITENSO = 0
          CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR))
        ENDIF
C
C  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
C     n'est pas la vitesse ni Rij)
        IVAR0 = 0
C
C    Sans prise en compte de la pression hydrostatique
C
        IPHYDP = 0
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA     , RA     , RA     ,
     &   RTP(1,IVAR), COEFA(1,ICLIMV) , COEFB(1,ICLIMV)  ,
C       POTR
     &   W4     , W5     , W6     ,
C       d POTR /dx   d POTR /dy   d POTR /dz
     &   W7     , W8     , W9    ,
     &   RDEVEL , RTUSER , RA     )
C
C
C   2.2 Calcul du champ electrique E = - grad (potR) : (-W4, -W5, -W6)
C  -------------------------------------------------
C
C   2.3 Calcul de la densite de courant j = sig E  :
C  -------------------------------------------------
C                                   PROPCE(IEL,IPCSIG) (-W4, -W5, -W6)
C
C
        IF( IPPMOD(IELJOU).GE.1 .OR. IPPMOD(IELARC).GE.1 ) THEN
          DO IEL = 1, NCEL
            PROPCE(IEL,IPCDC1)= - PROPCE(IEL,IPCSIG) * W4(IEL)
            PROPCE(IEL,IPCDC2)= - PROPCE(IEL,IPCSIG) * W5(IEL)
            PROPCE(IEL,IPCDC3)= - PROPCE(IEL,IPCSIG) * W6(IEL)
          ENDDO
        ENDIF
C
C   2.4 Calcul de l'Effet Joule j . E : PROPCE( .,IPPROC(IEFJOU) )
C  -----------------------------------------------------------------
C                                                           sig E.E
        DO IEL = 1, NCEL
          PROPCE(IEL,IPCEFJ)=
     &         PROPCE(IEL,IPCSIG)*(W4(IEL)**2+W5(IEL)**2+W6(IEL)**2)
        ENDDO
C
C
C   2.5 On imprime les extrema de E et j
C  -------------------------------------
C
        IF(MODNTL.EQ.0) THEN
C
          WRITE(NFECRA,1000)
C
C     Grad PotR = -E
          VAR    = W4(1)
          VRMIN = VAR
          VRMAX = VAR
          DO IEL = 1, NCEL
            VAR    = W4(IEL)
            VRMIN = MIN(VRMIN,VAR)
            VRMAX = MAX(VRMAX,VAR)
          ENDDO
          IF (IRANGP.GE.0) THEN
            CALL PARMIN (VRMIN)
            CALL PARMAX (VRMAX)
          ENDIF
          WRITE(NFECRA,1010)'Gr_PotRX',VRMIN,VRMAX
C
          VAR    = W5(1)
          VRMIN = VAR
          VRMAX = VAR
          DO IEL = 1, NCEL
            VAR    = W5(IEL)
            VRMIN = MIN(VRMIN,VAR)
            VRMAX = MAX(VRMAX,VAR)
          ENDDO
          IF (IRANGP.GE.0) THEN
            CALL PARMIN (VRMIN)
            CALL PARMAX (VRMAX)
          ENDIF
          WRITE(NFECRA,1010)'Gr_PotRY',VRMIN,VRMAX
C
          VAR    = W6(1)
          VRMIN = VAR
          VRMAX = VAR
          DO IEL = 1, NCEL
            VAR    = W6(IEL)
            VRMIN = MIN(VRMIN,VAR)
            VRMAX = MAX(VRMAX,VAR)
          ENDDO
          IF (IRANGP.GE.0) THEN
            CALL PARMIN (VRMIN)
            CALL PARMAX (VRMAX)
          ENDIF
          WRITE(NFECRA,1010)'Gr_PotRZ',VRMIN,VRMAX
C
          VAR    = -PROPCE(1,IPCSIG) * W4(1)
          VRMIN = VAR
          VRMAX = VAR
          DO IEL = 1, NCEL
            VAR = -PROPCE(IEL,IPCSIG) * W4(IEL)
            VRMIN = MIN(VRMIN,VAR)
            VRMAX = MAX(VRMAX,VAR)
          ENDDO
          IF (IRANGP.GE.0) THEN
            CALL PARMIN (VRMIN)
            CALL PARMAX (VRMAX)
          ENDIF
          WRITE(NFECRA,1010)'Cour_ReX',VRMIN,VRMAX
C
          VAR    = -PROPCE(1,IPCSIG) * W5(1)
          VRMIN = VAR
          VRMAX = VAR
          DO IEL = 1, NCEL
            VAR = -PROPCE(IEL,IPCSIG) * W5(IEL)
            VRMIN = MIN(VRMIN,VAR)
            VRMAX = MAX(VRMAX,VAR)
          ENDDO
          IF (IRANGP.GE.0) THEN
            CALL PARMIN (VRMIN)
            CALL PARMAX (VRMAX)
          ENDIF
          WRITE(NFECRA,1010)'Cour_ReY',VRMIN,VRMAX
C
          VAR    = -PROPCE(1,IPCSIG) * W6(1)
          VRMIN = VAR
          VRMAX = VAR
          DO IEL = 1, NCEL
            VAR = -PROPCE(IEL,IPCSIG) * W6(IEL)
            VRMIN = MIN(VRMIN,VAR)
            VRMAX = MAX(VRMAX,VAR)
          ENDDO
          IF (IRANGP.GE.0) THEN
            CALL PARMIN (VRMIN)
            CALL PARMAX (VRMAX)
          ENDIF
          WRITE(NFECRA,1010)'Cour_ReZ',VRMIN,VRMAX
C
        ENDIF
C
C
C
C=======================================================================
C 1.2 PRISE EN COMPTE ET AJOUT DES TERMES SOURCES ET VARIABLES
C    RELATIVES A LA PRISE EN COMPTE DU POTENTIEL COMPLEXE
C                           CAS IELJOU = 2
C=======================================================================
C
        IF(IPPMOD(IELJOU).GE.2 .OR. IPPMOD(IELJOU).EQ.4) THEN
C
C
C   3.1 Calcul du grad (potI) :  (W4, W5, W6)
C  ----------------------------
C
          IVAR = ISCA(IPOTI)
          ICLIMV = ICLRTP(IVAR,ICOEF)
C
          INC = 1
          ICCOCG = 1
          NSWRGP = NSWRGR(IVAR)
          IMLIGP = IMLIGR(IVAR)
          IWARNP = IWARNI(IVAR)
          EPSRGP = EPSRGR(IVAR)
          CLIMGP = CLIMGR(IVAR)
          EXTRAP = EXTRAG(IVAR)
C
C En periodique et parallele, echange avant calcul du gradient
C     C'est indispensable car on vient de calculer IVAR
C
C    Parallele
          IF(IRANGP.GE.0) THEN
            CALL PARCOM(RTP(1,IVAR))
C           ===========
          ENDIF
C
C    Periodique
          IF(IPERIO.EQ.1) THEN
            IDIMTE = 0
            ITENSO = 0
            CALL PERCOM
C           ===========
     &  ( IDIMTE , ITENSO ,
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR))
          ENDIF
C
C  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
C     n'est pas la vitesse ni Rij)
          IVAR0 = 0
C
C
C    SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
C
          IPHYDP = 0
C
          CALL GRDCEL
C         ===========
     &  ( IDBIA0 , IDBRA0 ,
     &    NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &    NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &    NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &    IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &    IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &    IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &    IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &    IDEVEL , ITUSER , IA     ,
     &    XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &    RA     , RA     , RA     ,
     &    RTP(1,IVAR), COEFA(1,ICLIMV) , COEFB(1,ICLIMV) ,
     &    W4     , W5     , W6     ,
C       d POTI /dx   d POTI /dy   d POTI /dz
     &    W7     , W8     , W9     ,
     &    RDEVEL , RTUSER , RA     )
C
C
C   3.2 Calcul du champ electrique Ei = - grad (potI) : (-W4, -W5, -W6)
C  -------------------------------------------------
C
C   3.3 Partie imaginaire de la densite de courant :
C  ------------------------------------------------
C                                   PROPCE(IEL,IPCSIG) (-W4, -W5, -W6)
C
        IF ( IPPMOD(IELJOU).EQ.4 ) THEN
          DO IEL = 1, NCEL
            PROPCE(IEL,IPCDI1)= -PROPCE(IEL,IPCSIG)*W4(IEL)
            PROPCE(IEL,IPCDI2)= -PROPCE(IEL,IPCSIG)*W5(IEL)
            PROPCE(IEL,IPCDI3)= -PROPCE(IEL,IPCSIG)*W6(IEL)
          ENDDO
        ENDIF
C
C
C   3.4 Effet Joule total : PROPCE( .,IPPROC(IEFJOU) )
C  ----------------------------------------------------
C
          DO IEL = 1, NCEL
C
C             ajout de la partie imaginaire et ...
            PROPCE(IEL,IPCEFJ) = PROPCE(IEL,IPCEFJ)
     &         + PROPCE(IEL,IPCSII)*(W4(IEL)**2+W5(IEL)**2+W6(IEL)**2)
C             ...division par 2
            PROPCE(IEL,IPCEFJ) = 0.5D0*PROPCE(IEL,IPCEFJ)
C
          ENDDO
C
C
C
C   3.5 On imprime les extrema de E et j
C  -------------------------------------
C
          IF(MODNTL.EQ.0) THEN
C
C     Grad PotI = -Ei
            VAR    = W4(1)
            VRMIN = VAR
            VRMAX = VAR
            DO IEL = 1, NCEL
              VAR    = W4(IEL)
              VRMIN = MIN(VRMIN,VAR)
              VRMAX = MAX(VRMAX,VAR)
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARMIN (VRMIN)
              CALL PARMAX (VRMAX)
            ENDIF
            WRITE(NFECRA,1010)'Gr_PotIX',VRMIN,VRMAX
C
            VAR    = W5(1)
            VRMIN = VAR
            VRMAX = VAR
            DO IEL = 1, NCEL
              VAR    = W5(IEL)
              VRMIN = MIN(VRMIN,VAR)
              VRMAX = MAX(VRMAX,VAR)
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARMIN (VRMIN)
              CALL PARMAX (VRMAX)
            ENDIF
            WRITE(NFECRA,1010)'Gr_PotIY',VRMIN,VRMAX
C
            VAR    = W6(1)
            VRMIN = VAR
            VRMAX = VAR
            DO IEL = 1, NCEL
              VAR    = W6(IEL)
              VRMIN = MIN(VRMIN,VAR)
              VRMAX = MAX(VRMAX,VAR)
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARMIN (VRMIN)
              CALL PARMAX (VRMAX)
            ENDIF
            WRITE(NFECRA,1010)'Gr_PotIZ',VRMIN,VRMAX
C
C     j=sigma E
            VAR    = -PROPCE(1,IPCSII) * W4(1)
            VRMIN = VAR
            VRMAX = VAR
            DO IEL = 1, NCEL
              VAR = -PROPCE(IEL,IPCSII) * W4(IEL)
              VRMIN = MIN(VRMIN,VAR)
              VRMAX = MAX(VRMAX,VAR)
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARMIN (VRMIN)
              CALL PARMAX (VRMAX)
            ENDIF
            WRITE(NFECRA,1010)'Cour_ImX',VRMIN,VRMAX
C
            VAR    = -PROPCE(1,IPCSII) * W5(1)
            VRMIN = VAR
            VRMAX = VAR
            DO IEL = 1, NCEL
              VAR = -PROPCE(IEL,IPCSII) * W5(IEL)
              VRMIN = MIN(VRMIN,VAR)
              VRMAX = MAX(VRMAX,VAR)
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARMIN (VRMIN)
              CALL PARMAX (VRMAX)
            ENDIF
            WRITE(NFECRA,1010)'Cour_ImY',VRMIN,VRMAX
C
            VAR    = -PROPCE(1,IPCSII) * W6(1)
            VRMIN = VAR
            VRMAX = VAR
            DO IEL = 1, NCEL
              VAR = -PROPCE(IEL,IPCSII) * W6(IEL)
              VRMIN = MIN(VRMIN,VAR)
              VRMAX = MAX(VRMAX,VAR)
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARMIN (VRMIN)
              CALL PARMAX (VRMAX)
            ENDIF
            WRITE(NFECRA,1010)'Cour_ImZ',VRMIN,VRMAX
C
            WRITE(NFECRA,1001)
C
          ENDIF
C
        ENDIF
C
      ENDIF
C     Fin du test IAPPEL = 1
C
C=======================================================================
C 2. DEUXIEME APPEL : A, B, JXB
C=======================================================================
C
      IF (IAPPEL.EQ.2) THEN
C
C=======================================================================
C 2.1 PRISE EN COMPTE ET AJOUT DES TERMES SOURCES ET VARIABLES RELATIVES
C    A LA PRISE EN COMPTE DE L'ARC ELECTRIQUE
C                           CAS IELARC = 1 OU 2
C=======================================================================
C
C
        IF(IPPMOD(IELARC).GE.1) THEN
          IPCLA1 = IPPROC(ILAPLA(1))
          IPCLA2 = IPPROC(ILAPLA(2))
          IPCLA3 = IPPROC(ILAPLA(3))
        ENDIF
C
C   4.1 ARC ELECTRIQUE DIT 3D : IELARC = 2
C  ----------------------------------------
C         ON PASSE PAR LA RESOLUTION D'EQUATIONS DE POISSON SUR LES
C        -----------------------------------------------------------
C           COMPOSANTES DU POTENTIEL VECTEUR
C          ----------------------------------
C
        IF( IPPMOD(IELARC).GE.2 ) THEN
C
C --> Calcul des composantes du champ magnetique B = (W1, W2, W3)
C    ********************************************
C
C
C    Sur Ax
C
          IVAR = ISCA(IPOTVA(1))
          ICLIMV = ICLRTP(IVAR,ICOEF)
C
          INC = 1
          ICCOCG = 1
          NSWRGP = NSWRGR(IVAR)
          IMLIGP = IMLIGR(IVAR)
          IWARNP = IWARNI(IVAR)
          EPSRGP = EPSRGR(IVAR)
          CLIMGP = CLIMGR(IVAR)
          EXTRAP = EXTRAG(IVAR)
C
C En periodique et parallele, echange avant calcul du gradient
C     C'est indispensable car on vient de calculer IVAR
C
C    Parallele
          IF(IRANGP.GE.0) THEN
            CALL PARCOM(RTP(1,IVAR))
C           ===========
          ENDIF
C
C    Periodique
          IF(IPERIO.EQ.1) THEN
            IDIMTE = 0
            ITENSO = 0
            CALL PERCOM
C           ===========
     &  ( IDIMTE , ITENSO ,
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR))
          ENDIF
C
C  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
C     n'est pas la vitesse ni Rij)
C
          IVAR0 = 0
C
C    SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
C
          IPHYDP = 0
C
          CALL GRDCEL
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RA     , RA     , RA     ,
     &   RTP(1,IVAR), COEFA(1,ICLIMV) , COEFB(1,ICLIMV)  ,
     &   W4     , W5     , W6     ,
C       d Ax /dx   d Ax /dy   d Ax /dz
     &   W7     , W8     , W9    ,
     &   RDEVEL , RTUSER , RA     )
C
C       B = rot A
C
          DO IEL = 1, NCEL
            W1(IEL)=  ZERO
            W2(IEL)=  W6(IEL)
            W3(IEL)= -W5(IEL)
          ENDDO
C
C    Sur Ay
C
          IVAR = ISCA(IPOTVA(2))
          ICLIMV = ICLRTP(IVAR,ICOEF)
C
          INC = 1
          ICCOCG = 1
          NSWRGP = NSWRGR(IVAR)
          IMLIGP = IMLIGR(IVAR)
          IWARNP = IWARNI(IVAR)
          EPSRGP = EPSRGR(IVAR)
          CLIMGP = CLIMGR(IVAR)
          EXTRAP = EXTRAG(IVAR)
C
C En periodique et parallele, echange avant calcul du gradient
C     C'est indispensable car on vient de calculer IVAR
C
C    Parallele
          IF(IRANGP.GE.0) THEN
            CALL PARCOM(RTP(1,IVAR))
C           ===========
          ENDIF
C
C    Periodique
          IF(IPERIO.EQ.1) THEN
            IDIMTE = 0
            ITENSO = 0
            CALL PERCOM
C         ===========
     &  ( IDIMTE , ITENSO ,
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR))
          ENDIF
C
C  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
C     n'est pas la vitesse ni Rij)
C
          IVAR0 = 0
C
C    SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
C
          IPHYDP = 0
C
          CALL GRDCEL
C         ===========
     &  ( IDBIA0 , IDBRA0 ,
     &    NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &    NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &    NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &    IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &    IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &    IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &    IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &    IDEVEL , ITUSER , IA     ,
     &    XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &    RA     , RA     , RA     ,
     &    RTP(1,IVAR), COEFA(1,ICLIMV) , COEFB(1,ICLIMV) ,
     &    W4     , W5     , W6     ,
C       d Ay /dx   d Ay /dy   d Ay /dz
     &    W7     , W8     , W9     ,
     &    RDEVEL , RTUSER , RA     )
C
C       B = rot A
C
          DO IEL = 1, NCEL
            W1(IEL)= W1(IEL) - W6(IEL)
            W2(IEL)= W2(IEL) + ZERO
            W3(IEL)= W3(IEL) + W4(IEL)
          ENDDO
C
C    Sur Az
C
          IVAR = ISCA(IPOTVA(3))
          ICLIMV = ICLRTP(IVAR,ICOEF)
C
          INC = 1
          ICCOCG = 1
          NSWRGP = NSWRGR(IVAR)
          IMLIGP = IMLIGR(IVAR)
          IWARNP = IWARNI(IVAR)
          EPSRGP = EPSRGR(IVAR)
          CLIMGP = CLIMGR(IVAR)
          EXTRAP = EXTRAG(IVAR)
C
C En periodique et parallele, echange avant calcul du gradient
C     C'est indispensable car on vient de calculer IVAR
C
C    Parallele
          IF(IRANGP.GE.0) THEN
            CALL PARCOM(RTP(1,IVAR))
C           ===========
          ENDIF
C
C    Periodique
          IF(IPERIO.EQ.1) THEN
            IDIMTE = 0
            ITENSO = 0
            CALL PERCOM
C           ===========
     &  ( IDIMTE , ITENSO ,
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &    RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR))
          ENDIF
C
C  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
C     n'est pas la vitesse ni Rij)
C
          IVAR0 = 0
C
C    SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE
C
          IPHYDP = 0
C
          CALL GRDCEL
C         ===========
     &  ( IDBIA0 , IDBRA0 ,
     &    NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &    NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &    NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &    IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &    IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &    IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &    IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &    IDEVEL , ITUSER , IA     ,
     &    XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &    RA     , RA     , RA     ,
     &    RTP(1,IVAR), COEFA(1,ICLIMV) , COEFB(1,ICLIMV) ,
     &    W4     , W5     , W6     ,
C       d Az /dx   d Az /dy   d Az /dz
     &    W7     , W8     , W9     ,
     &    RDEVEL , RTUSER , RA     )
C
C       B = rot A
C
          DO IEL = 1, NCEL
            W1(IEL)= W1(IEL) + W5(IEL)
            W2(IEL)= W2(IEL) - W4(IEL)
            W3(IEL)= W3(IEL) + ZERO
          ENDDO
C
C   4.2 CAS D'ARC AXISYMETRIQUE : IELARC = 1
C  ------------------------------------------
C         ON PEUT UTILISER LE THEOREME D'AMPERE
C        ---------------------------------------
C
        ELSE IF(IPPMOD(IELARC).EQ.1) THEN
C
C      Calcul du Champ magnetique calcule :
C     ************************************
C
C       B = / j d S
C
C       Cette version n'a pas ete developpe pour l'instant :
C         - elle ne fonctionne que dans le cas d'un arc axisymetrique
C         - elle demande donc un reperage des plans perpendiculaires
C           a l'arc
C         - elle necessite la connaissance des coordonnees de chaque
C           point du maillage
C
          WRITE(NFECRA,2000)
          CALL CSEXIT (1)
C
        ENDIF
C
C   4.3 ARC ELECTRIQUE : CALCUL DES FORCES DE LAPLACE j X B
C  ---------------------------------------------------------
C         POUR TOUS LES CAS D'ARCS ELECTRIQUES
C        --------------------------------------
C
        IF( IPPMOD(IELARC) .GE. 1 ) THEN
          DO IEL = 1, NCEL
            PROPCE(IEL,IPCLA1)= PROPCE(IEL,IPCDC2) * W3(IEL)
     &                         -PROPCE(IEL,IPCDC3) * W2(IEL)
            PROPCE(IEL,IPCLA2)= PROPCE(IEL,IPCDC3) * W1(IEL)
     &                         -PROPCE(IEL,IPCDC1) * W3(IEL)
            PROPCE(IEL,IPCLA3)= PROPCE(IEL,IPCDC1) * W2(IEL)
     &                         -PROPCE(IEL,IPCDC2) * W1(IEL)
          ENDDO
        ENDIF
C
C
C   4.4 Impression de B en arc electrique
C  --------------------------------------
C
        IF( IPPMOD(IELARC) .GE. 2 ) THEN

          IF(MODNTL.EQ.0) THEN
C
            WRITE(NFECRA,1000)
C
C     B=rot A
            VAR    = W1(1)
            VRMIN = VAR
            VRMAX = VAR
            DO IEL = 1, NCEL
              VAR    = W1(IEL)
              VRMIN = MIN(VRMIN,VAR)
              VRMAX = MAX(VRMAX,VAR)
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARMIN (VRMIN)
              CALL PARMAX (VRMAX)
            ENDIF
            WRITE(NFECRA,1010)'Ch_MagX ',VRMIN,VRMAX
C
            VAR    = W2(1)
            VRMIN = VAR
            VRMAX = VAR
            DO IEL = 1, NCEL
              VAR    = W2(IEL)
              VRMIN = MIN(VRMIN,VAR)
              VRMAX = MAX(VRMAX,VAR)
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARMIN (VRMIN)
              CALL PARMAX (VRMAX)
            ENDIF
            WRITE(NFECRA,1010)'Ch_MagY ',VRMIN,VRMAX
C
            VAR    = W3(1)
            VRMIN = VAR
            VRMAX = VAR
            DO IEL = 1, NCEL
              VAR    = W3(IEL)
              VRMIN = MIN(VRMIN,VAR)
              VRMAX = MAX(VRMAX,VAR)
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARMIN (VRMIN)
              CALL PARMAX (VRMAX)
            ENDIF
            WRITE(NFECRA,1010)'Ch_MagZ ',VRMIN,VRMAX
C
            WRITE(NFECRA,1001)
C
          ENDIF
C
        ENDIF
C
      ENDIF
C     Fin du test IAPPEL = 2
C
C--------
C FORMATS
C--------
C
C
 1000 FORMAT(/,
     &'-----------------------------------------                    ',/,
     &'   Variable         Minimum       Maximum                    ',/,
     &'-----------------------------------------                    '  )
 1010 FORMAT(
     &'v  ',A8,'    ',E12.5,'  ',E12.5                                 )
 1001 FORMAT(
     &'-----------------------------------------                    '  )
C
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DANS ELTSSC            ',/,
     &'@    *********                                               ',/,
     &'@     LA VERSION ARC ELECTRIQUE AVEC THEOREME D''AMPERE      ',/,
     &'@       (IELARC = 1)    N''EST PAS DISPONIBLE                ',/,
     &'@                                                            ',/,
     &'@            VEUILLEZ UTILISER LA VERSION DITE 3D            ',/,
     &'@                      IELARC = 2                            ',/,
     &'@   AVEC EQUATION DE TRANSPORT SUR LE POTENTIEL VECTEUR A    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
