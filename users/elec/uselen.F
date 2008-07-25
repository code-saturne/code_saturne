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
                        SUBROUTINE USELEN
C                       *****************
C     ---------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , NUMMAI ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCELPS , NFACPS , NFBRPS ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   LSTCEL , LSTFAC , LSTFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     ,
     &   TRACEL , TRAFAC , TRAFBR , RDEVEL , RTUSER , RA     )
C     ---------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC POUR LA SORTIE POST-TRAITEMENT MODULE ELECTRIQUE
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
CARGU ! NUMMAI       ! EC !  -> ! NUMERO DU MAILLAGE POST              !
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
CARGU ! NCELPS       ! E  !  -> ! NOMBRE DE CELLULES DU MAILLAGE POST  !
CARGU ! NFACPS       ! E  !  -> ! NOMBRE DE FACES INTERIEUR POST       !
CARGU ! NFBRPS       ! E  !  -> ! NOMBRE DE FACES DE BORD POST         !
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
CARGU ! LSTCEL(NCELPS! TE !  -> ! LISTE DES CELLULES DU MAILLAGE POST  !
CARGU ! LSTFAC(NFACPS! TE !  -> ! LISTE DES FACES INTERIEURES POST     !
CARGU ! LSTFBR(NFBRPS! TE !  -> ! LISTE DES FACES DE BORD POST         !
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
CARGU ! TRACEL(*)    ! TR ! <-> ! TAB REEL VALEURS CELLULES POST       !
CARGU ! TRAFAC(*)    ! TR ! <-> ! TAB REEL VALEURS FACES INT. POST     !
CARGU ! TRAFBR(*)    ! TR ! <-> ! TAB REEL VALEURS FACES BORD POST     !
CARGU ! W1-W2        ! TR !  -  ! TAB REEL POUR CALCUL GRADIENT        !
CARGU ! (NCELET,3)   !    !     !                                      !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
CARGU ! RA(*)        ! TR !  -  ! MACRO TABLEAU REEL                   !
CARGU !______________!____!_____!______________________________________!

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
      INCLUDE "dimfbr.h"
      INCLUDE "paramx.h"
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
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
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NUMMAI
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NCELPS , NFACPS , NFBRPS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IDIMT
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          LSTCEL(NCELPS), LSTFAC(NFACPS), LSTFBR(NFBRPS)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NDIMFB,*)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION TRACEL(NCELPS*3)
      DOUBLE PRECISION TRAFAC(NFACPS*3), TRAFBR(NFBRPS*3)
      DOUBLE PRECISION W1(NCELET,3), W2(NCELET,3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*32     NAMEVR
      INTEGER          IDEBIA, IDEBRA, IEL   , ILOC
      INTEGER          IVAR  , IVAR0 , INC   , ICCOCG
      INTEGER          IPHYDP, NSWRGP, IMLIGP, IWARNP, ICLIMV
      INTEGER          IPCSII
      INTEGER          IENTLA, IVARPR
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP
      DOUBLE PRECISION RBID(1)
C
C***********************************************************************
C=======================================================================
C 0.  PAR DEFAUT, ON CONSIDERE QUE LE SOUS PROGRAMME CI-DESSOUS CONVIENT
C       A L'UTILISATEUR, C'EST-A-DIRE QUE LA MISE EN OEUVRE DU MODULE
C       ELECTRIQUE DECLENCHE LA PRODUCTION DE CHAMPS STANDARD DANS LE
C       POST-TRAITEMENT.
C     L'UTILISATEUR N'A PAS A MODIFIER LE PRESENT SOUS-PROGRAMME DANS
C       LES CONDITIONS D'UTILISATION STANDARD.
C     DANS LE CAS OU IL SOUHAITE PRODUIRE DES VARIABLES SUPPLEMENTAIRES
C       ILPEUT LES AJOUTER A LA FIN, VOIR LA DOCUMENTATION DE USEEVO
C=======================================================================
C
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IF(NUMMAI.EQ.-1) THEN
C
C=======================================================================
C 1.   Gradient du potentiel reel
C=======================================================================
C
        IDIMT  = 3
        NAMEVR = 'Gr_PotR'
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
C     Cela a deja ete fait puisqu'on a deja fait le calcul de cette variable
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
     &   W1(1,1) , W1(1,2) , W1(1,3) ,
C       d POTR /dx   d POTR /dy   d POTR /dz
     &   W2(1,1) , W2(1,2) , W2(1,3) ,
     &   RDEVEL , RTUSER , RA     )
C
C       Le gradient est defini sur le maillage principal tout entier ;
C       inutile de le recopier, on utilise l'indirection (IVARPR = 1),
C       et les valeurs sont non entrelacees (IENTLA = 0)
        IENTLA = 0
        IVARPR = 1
C
        CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
C       ===========
     &              NTCABS, TTCABS, W1, RBID, RBID)
C
C=======================================================================
C 2.   Gradient du potentiel imaginaire si Joule
C=======================================================================
C
        IF (IPPMOD(IELJOU).EQ.2 .OR. IPPMOD(IELJOU).EQ.4) THEN
C
          IDIMT  = 3
          NAMEVR = 'Gr_PotI'
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
C     Cela a deja ete fait puisqu'on a deja fait le calcul de cette variable
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
C       POTI
     &   W1(1,1) , W1(1,2) , W1(1,3) ,
C       d POTR /dx   d POTR /dy   d POTR /dz
     &   W2(1,1) , W2(1,2) , W2(1,3) ,
     &   RDEVEL , RTUSER , RA     )
C
C         Le gradient est defini sur le maillage principal tout entier ;
C         inutile de le recopier, on utilise l'indirection (IVARPR = 1),
C         et les valeurs sont non entrelacees (IENTLA = 0)
          IENTLA = 0
          IVARPR = 1
C
          CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
C         ===========
     &                NTCABS, TTCABS, W1, RBID, RBID)
C
        ENDIF
C
C=======================================================================
C 3.   Courant imaginaire si Joule
C=======================================================================
C
        IF(IPPMOD(IELJOU).EQ.2 .OR. IPPMOD(IELJOU).EQ.4 ) THEN
C
          IDIMT  = 3
          NAMEVR = 'Cour_Im'
C
          IVAR = ISCA(IPOTI)
          ICLIMV = ICLRTP(IVAR,ICOEF)
C
C     Commme dans elflux
          IPCSII = IPPROC(IVISLS(IPOTI))
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
C     Cela a deja ete fait puisqu'on a deja fait le calcul de cette variable
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
C       POTI
     &   W1(1,1) , W1(1,2) , W1(1,3) ,
C       d POTR /dx   d POTR /dy   d POTR /dz
     &   W2(1,1) , W2(1,2) , W2(1,3) ,
     &   RDEVEL , RTUSER , RA     )
C
          DO ILOC = 1, NCELPS
            IEL = LSTCEL(ILOC)
            TRACEL(ILOC)          = -PROPCE(IEL,IPCSII)*W1(IEL,1)
            TRACEL(ILOC+NCELPS)   = -PROPCE(IEL,IPCSII)*W1(IEL,2)
            TRACEL(ILOC+2*NCELPS) = -PROPCE(IEL,IPCSII)*W1(IEL,3)
          ENDDO
C
C         La variable est définie sur le tableau de travail. On
C         a deja utilise l'indirection via LSTCEL (donc IVARPR = 0),
C         et les valeurs sont non entrelacees (IENTLA = 0)
          IENTLA = 0
          IVARPR = 0
C
          CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
C         ===========
     &                NTCABS, TTCABS, TRACEL, RBID, RBID)
C
        ENDIF
C
C=======================================================================
C 5.   Champ Magnetique si Arc
C=======================================================================
C
        IF( IPPMOD(IELARC).GE.2 ) THEN
C
          IDIMT  = 3
          NAMEVR = 'Ch_Mag'
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
C     Cela a deja ete fait puisqu'on a deja fait le calcul de cette variable
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
     &   W1(1,1) , W1(1,2) , W1(1,3) ,
C       d Ax /dx   d Ax /dy   d Ax /dz
     &   W2(1,1) , W2(1,2) , W2(1,3) ,
     &   RDEVEL , RTUSER , RA     )
C
C       B = rot A
C
          DO ILOC = 1, NCELPS
            IEL = LSTCEL(ILOC)
            TRACEL(ILOC)          =  ZERO
            TRACEL(ILOC+NCELPS)   =  W1(IEL,3)
            TRACEL(ILOC+2*NCELPS) = -W1(IEL,2)
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
C     Cela a deja ete fait puisqu'on a deja fait le calcul de cette variable
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
     &    W1(1,1) , W1(1,2) , W1(1,3) ,
C       d Ay /dx   d Ay /dy   d Ay /dz
     &    W2(1,1) , W2(1,2) , W2(1,3) ,
     &    RDEVEL , RTUSER , RA     )
C
C       B = rot A
C
          DO ILOC = 1, NCELPS
            IEL = LSTCEL(ILOC)
            TRACEL(ILOC)          = TRACEL(ILOC)          - W1(IEL,3)
            TRACEL(ILOC+NCELPS)   = TRACEL(ILOC + NCELPS) + ZERO
            TRACEL(ILOC+2*NCELPS) = TRACEL(ILOC+2*NCELPS) + W1(IEL,1)
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
C     Cela a deja ete fait puisqu'on a deja fait le calcul de cette variable
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
     &    W1(1,1) , W1(1,2) , W1(1,3) ,
C       d Az /dx   d Az /dy   d Az /dz
     &    W2(1,1) , W2(1,2) , W2(1,3) ,
     &    RDEVEL , RTUSER , RA     )
C
C       B = rot A
C
          DO ILOC = 1, NCELPS
            IEL = LSTCEL(ILOC)
            TRACEL(ILOC)          = TRACEL(ILOC)          + W1(IEL,2)
            TRACEL(ILOC+NCELPS)   = TRACEL(ILOC+NCELPS)   - W1(IEL,1)
            TRACEL(ILOC+2*NCELPS) = TRACEL(ILOC+2*NCELPS) + ZERO
          ENDDO
C
C         La variable est définie sur le tableau de travail. On
C         a deja utilise l'indirection via LSTCEL (donc IVARPR = 0),
C         et les valeurs sont non entrelacees (IENTLA = 0)
          IENTLA = 0
          IVARPR = 0
C
          CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
C         ===========
     &                NTCABS, TTCABS, TRACEL, RBID, RBID)
C
        ENDIF
C
C
C=======================================================================
C 4.   Module et Argument du potentiel si IELJOU = 4
C=======================================================================
C
        IF (IPPMOD(IELJOU).EQ.4) THEN
C
          IDIMT  = 1
          NAMEVR = 'ModPot'
C
          IVAR = ISCA(IPOTR)
C
          DO ILOC = 1, NCELPS
            IEL = LSTCEL(ILOC)
            TRACEL(ILOC) =
     &        SQRT( RTP(IEL,ISCA(IPOTR))*RTP(IEL,ISCA(IPOTR))
     &             +RTP(IEL,ISCA(IPOTI))*RTP(IEL,ISCA(IPOTI)) )
          ENDDO
C
          IENTLA = 0
          IVARPR = 0
C
          CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
C         ===========
     &                NTCABS, TTCABS, TRACEL, RBID, RBID)
C
          IDIMT  = 1
          NAMEVR = 'ArgPot'
C
          IVAR = ISCA(IPOTR)
C
          DO ILOC = 1, NCELPS
            IEL = LSTCEL(ILOC)
            IF ( RTP(IEL,ISCA(IPOTR)) .NE. 0.D0 ) THEN
              IF ( RTP(IEL,ISCA(IPOTR)) .GE. 0.D0 ) THEN
                TRACEL(ILOC) =
     &           ATAN( RTP(IEL,ISCA(IPOTI))/RTP(IEL,ISCA(IPOTR)))
              ELSE
                IF ( RTP(IEL,ISCA(IPOTI)) .GT. 0.D0 ) THEN
                  TRACEL(ILOC) =
     &              4.D0*ATAN(1.D0)
     &             +ATAN( RTP(IEL,ISCA(IPOTI))
     &                   /RTP(IEL,ISCA(IPOTR)))
                ELSE
                  TRACEL(ILOC) =
     &             -4.D0*ATAN(1.D0)
     &             +ATAN( RTP(IEL,ISCA(IPOTI))
     &                   /RTP(IEL,ISCA(IPOTR)))
                ENDIF
              ENDIF
            ELSE
              TRACEL(ILOC) = 2.D0*ATAN(1.D0)
            ENDIF
C
            IF (TRACEL(ILOC) .LT. 0.D0) THEN
              TRACEL(ILOC) = TRACEL(ILOC) + 8.D0**ATAN(1.D0)
            ENDIF
C
          ENDDO
C
          IENTLA = 0
          IVARPR = 0
C
          CALL PSTEVA(NUMMAI, NAMEVR, IDIMT, IENTLA, IVARPR,
C         ===========
     &                NTCABS, TTCABS, TRACEL, RBID, RBID)
C
        ENDIF
C
      ENDIF
C
      RETURN
C
      END
c@z
