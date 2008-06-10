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
                        SUBROUTINE CALHYD
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR , INDHYD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   FEXTX  , FEXTY  , FEXTZ  ,
     &   DFEXTX , DFEXTY , DFEXTZ ,
     &   PHYDR  , FLUMAS , FLUMAB ,
     &   COEFA  , COEFB  ,
     &   VISCF  , VISCB  ,
     &   DAM    , XAM    , DAG    , XAG    ,
     &   DRTP   , SMBR   ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , W10    ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION D'UNE EQUATION DE POISSON SUR LA PRESSION HYDROSTATIQUE
CFONC                 DIV( GRAD(P) ) = DIV( F )
CFONC                      ----             -
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
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
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE PHASE                      !
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
CARGU ! IFACLG(2,NFAC! TE !  -  ! TAB ENTIER MULTIGRILLE               !
CARGU ! IRESPR(NCELET! TE !  -  ! TAB ENTIER MULTIGRILLE               !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! INDHYD       ! E  ! <-  ! INDICATEUR DE MISE A JOUR DE PHYDR   !
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
CARGU ! PHYDR(NCELET)! TR !  -> ! INCREMENT DE PRESSION HYDROSTATIQUE  !
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR)    !    !     !    FACES DE BORD                     !
CARGU ! FEXTX,FEXTY  ! TR !  -> ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU ! FEXTZ(NCELET)!    !     !    HYDROSTATIQUE                     !
CARGU ! DFEXTX,DFEXTY! TR !  -> ! INCREMENT DE FORCE EXTERIEURE        !
CARGU ! DFEXTZ(NCELET!    !     !    GENERANT LA PRESSION HYDROSTATIQUE!
CARGU ! VISCF(NFAC)  ! TR !  -  ! 1*SURFACE/DIST AUX FACES INTERNES    !
CARGU ! VISCB(NFABOR ! TR !  -  ! 1*SURFACE/DIST AUX FACES DE BORD     !
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DAG(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAG(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DRTP(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! SMBR  (NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! W1..10(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
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
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS
C
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IFACLG(2,NFAC), IRESPR(NCELET), INDHYD
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION FEXTX(NCELET),FEXTY(NCELET),FEXTZ(NCELET)
      DOUBLE PRECISION DFEXTX(NCELET),DFEXTY(NCELET),DFEXTZ(NCELET)
      DOUBLE PRECISION PHYDR(NCELET)
      DOUBLE PRECISION FLUMAS(NFAC), FLUMAB(NFABOR)
      DOUBLE PRECISION COEFA(NFABOR), COEFB(NFABOR)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DAG(NCELET), XAG(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET)
      DOUBLE PRECISION SMBR(NCELET)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION W10(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      INTEGER          IDEBIA, IDEBRA
      INTEGER          ICCOCG, INC   , INIT  , ISYM  , IPOL  , ISQRT
      INTEGER          IEL   , ICAL
      INTEGER          IRESLP, NSWMPR
      INTEGER          ISWEEP, NITERF, ICYCLE
      INTEGER          IPHYDP
      INTEGER          NSWRGP, IMLIGP, IWARNP
      INTEGER          IPRIPH
      INTEGER          IINVPE
      INTEGER          IDIFFP, ICONVP, NDIRCP
      INTEGER          NITMAP, IMGRP , NCYMAP, NITMGP
      DOUBLE PRECISION RESIDU, RNORM , RNRMF , RNRMDF
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP, EPSILP
      DOUBLE PRECISION PRECRE, PRECAB, THETAP
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
C --- Memoire
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
C --- Variables
      IPRIPH = IPR(IPHAS)
C
C
C
C --- Options de resolution
C     Symetrique
C     Preconditionnement diagonal par defaut
      ISYM  = 1
      IF (IRESOL(IPRIPH).EQ.-1) THEN
        IRESLP = 0
        IPOL   = 0
      ELSE
        IRESLP = MOD(IRESOL(IPRIPH),1000)
        IPOL   = (IRESOL(IPRIPH)-IRESLP)/1000
      ENDIF
C
      ISQRT = 1
C
C     TEST DE VARIATION DE LA PRESSION HYDROSTATIQUE EN SORTIE
C
C
C     on regarde si les terme source ont varie
C     on ne passe dans calhyd que si on a des faces de sortie std
C     la precision pour les tests est a peu pres arbitraire.
      PRECRE = SQRT(EPZERO)
      PRECAB = 1.D2*EPZERO
C
      ICAL = 0
      DO IEL = 1, NCEL
        RNRMF  = FEXTX(IEL)**2+FEXTY(IEL)**2+FEXTZ(IEL)**2
        RNRMDF = DFEXTX(IEL)**2+DFEXTY(IEL)**2+DFEXTZ(IEL)**2
        IF ((RNRMDF.GE.PRECRE*RNRMF).AND.(RNRMDF.GE.PRECAB)) THEN
          ICAL = 1
        ENDIF
      ENDDO
      IF (IRANGP.GE.0) THEN
        CALL PARCPT (ICAL)
      ENDIF
      IF (ICAL.EQ.0) THEN
        DO IEL = 1,NCEL
          PHYDR(IEL) = 0.D0
        ENDDO
        INDHYD = 0
        RETURN
      ENDIF
C
      IF ( MOD(NTCABS,NTLIST).EQ.0 .OR. IWARNI(IU(1)) .GE.0 )
     &     WRITE(NFECRA,1000)
C
 1000 FORMAT(
     &'  Calcul de la pression hydrostatique : ',/,
     &'         mise a jour des Dirichlets en sortie (CALHYD)',/)
C
      INDHYD = 1
C
C=======================================================================
C 2.  PREPARATION DE LA MATRICE DU SYSTEME A RESOUDRE
C=======================================================================
C
C ---> TERME INSTATIONNAIRE
C
      DO IEL = 1, NCEL
        W1(IEL) = 0.D0
      ENDDO
C
C ---> "VITESSE" DE DIFFUSION FACETTE
C
      DO IEL = 1, NCEL
        W10(IEL) = 1.D0
      ENDDO
C
      CALL VISCFA
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IMVISF ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W10    ,
     &   VISCF  , VISCB  ,
     &   RDEVEL , RTUSER , RA     )
C
C
      ICONVP = 0
      IDIFFP = 1
C  On resout avec des CL de flux nul partout
      NDIRCP = 0
C
      THETAP = 1.D0
      CALL MATRIX
C     ===========
     & ( NCELET , NCEL   , NFAC   , NFABOR ,
     &   ICONVP , IDIFFP , NDIRCP ,
     &   ISYM   , NFECRA ,
     &   THETAP ,
     &   IFACEL , IFABOR ,
     &   COEFB  , W1     ,
     &   FLUMAS , FLUMAB , VISCF  , VISCB  ,
     &   DAM    , XAM    )
C
C
C=======================================================================
C 4.  INITIALISATION DU FLUX DE MASSE
C=======================================================================
C
C
C     PROJECTION AUX FACES DES TERMES SOURCES
      INIT   = 1
      INC    = 0
      ICCOCG = 1
      NSWRGP = NSWRGR(IPRIPH)
      IMLIGP = IMLIGR(IPRIPH)
      IWARNP = IWARNI(IPRIPH)
      EPSRGP = EPSRGR(IPRIPH)
      CLIMGP = CLIMGR(IPRIPH)
C
      CALL PROJTS
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DFEXTX , DFEXTY , DFEXTZ ,
     &   COEFB     ,
     &   FLUMAS, FLUMAB ,
     &   VISCF  , VISCB  ,
     &   W10    , W10    , W10    ,
     &   RDEVEL , RTUSER , RA     )
C
      INIT = 1
      CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &            IFACEL,IFABOR,FLUMAS,FLUMAB,W7)
      CALL PRODSC(NCELET,NCEL,ISQRT,W7,W7,RNORM)
C
C=======================================================================
C 6.  BOUCLES SUR LES NON ORTHOGONALITES (RESOLUTION)
C=======================================================================
C
C --- Nombre de sweeps
      NSWMPR = NSWRSM(IPRIPH)
C
C --- Mise a zero des variables
C       RTP(.,IPR) sera l'increment de pression cumule
C       DRTP       sera l'increment d'increment a chaque sweep
C       W7         sera la divergence du flux de masse predit
      DO IEL = 1,NCEL
        PHYDR(IEL) = 0.D0
        DRTP(IEL) = 0.D0
        SMBR(IEL) = 0.D0
      ENDDO
C
C
C --- Boucle de reconstruction : debut
      DO ISWEEP = 1, NSWMPR
C
C --- Mise a jour du second membre
C     (signe "-" a cause de celui qui est implicitement dans la matrice)
        DO IEL = 1, NCEL
          SMBR(IEL) = - W7(IEL) - SMBR(IEL)
        ENDDO
C
C --- Test de convergence du calcul
C
        CALL PRODSC(NCELET,NCEL,ISQRT,SMBR,SMBR,RESIDU)
        IF (IWARNI(IPRIPH).GE.2) THEN
           CHAINE = 'PresHydr'
           WRITE(NFECRA,1400)CHAINE(1:8),ISWEEP,RESIDU
        ENDIF
C
CMO IL FAUDRA VERIFIER LA PERTINENCE DU TEST
C
        IF( RESIDU .LE. 10.D0*EPSILO(IPRIPH)*RNORM ) THEN
C     Si convergence,  sortie
C
          GOTO 101
C
        ENDIF
C
C --- Resolution implicite sur l'increment d'increment DRTP
        DO IEL = 1, NCEL
          DRTP(IEL) = 0.D0
        ENDDO
C
        CHAINE = 'PresHydr'
        NITMAP = NITMAX(IPRIPH)
        IMGRP  = IMGR  (IPRIPH)
        NCYMAP = NCYMAX(IPRIPH)
        NITMGP = NITMGF(IPRIPH)
        IWARNP = IWARNI(IPRIPH)
        EPSILP = EPSILO(IPRIPH)
        IINVPE = 1
C
        CALL INVERS
C       ===========
     & ( CHAINE(1:8)     , IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ISYM   , IPOL   , IRESLP , NITMAP , IMGRP  ,
     &   NCYMAP , NITMGP ,
     &   IWARNP , NFECRA , NITERF , ICYCLE , IINVPE ,
     &   EPSILP , RNORM  , RESIDU ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DAM    , XAM    , SMBR   , DRTP   ,
     &   DAG    , XAG    , W1     , W2     ,
     &   W3     , W4     , W5     , W6     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
C
        IF( ISWEEP.EQ.NSWMPR ) THEN
C     Mise a jour de l'increment de pression
          DO IEL = 1, NCEL
            PHYDR(IEL) = PHYDR(IEL) + DRTP(IEL)
          ENDDO
C
C
        ELSE
C
C --- Si ce n'est pas le dernier sweep
C       Mise a jour de l'increment de pression et calcul direct de la
C       partie en gradient d'increment de pression du second membre
C       (avec reconstruction)
C
          DO IEL = 1, NCEL
            PHYDR(IEL) = PHYDR(IEL) + DRTP(IEL)
          ENDDO
C
          ICCOCG = 1
          INIT = 1
          INC  = 1
          NSWRGP = NSWRGR(IPRIPH)
          IMLIGP = IMLIGR(IPRIPH)
          IWARNP = IWARNI(IPRIPH)
          EPSRGP = EPSRGR(IPRIPH)
          CLIMGP = CLIMGR(IPRIPH)
          EXTRAP = 0.D0
          IPHYDP = 1
C
          CALL ITRGRP
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DFEXTX , DFEXTY , DFEXTZ ,
     &   PHYDR  ,
     &   COEFA  , COEFB  ,
     &   VISCF  , VISCB  ,
     &   W10         , W10         , W10         ,
     &   SMBR   ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF
C
      ENDDO
C --- Boucle de reconstruction : fin
C
      IF(IWARNI(IPRIPH).GE.2) THEN
         CHAINE = 'PresHydr'
         WRITE( NFECRA,1600)CHAINE(1:8),NSWMPR
      ENDIF
C
 101  CONTINUE
C
C--------
C FORMATS
C--------
C
 1400 FORMAT(1X,A8,' : SWEEP = ',I5,' NORME SECOND MEMBRE = ',E14.6)
 1600 FORMAT(
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ',A8 ,' ETAPE DE PRESSION HYDROSTATIQUE     ',/,
     &'@    *********                                               ',/,
     &'@  Nombre d''iterations maximal ',I10   ,' atteint           ',/,
     &'@                                                            '  )
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
