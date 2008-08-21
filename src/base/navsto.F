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
                        SUBROUTINE NAVSTO
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , ITERNS , ICVRGE ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFACLG , IRESPR , ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   TSLAGR , COEFA  , COEFB  , FRCXT  ,
     &   TRAVA  , XIMPA  , UVWK   ,
     &   VISCF  , VISCB  , VISCFI , VISCBI ,
     &   DAM    , XAM    , DAG    , XAG    ,
     &   DRTP   , TRAV   , SMBR   , ROVSDT ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , W10    , DFRCXT , FRCHY  , DFRCHY ,
     &   COEFU  , ESFLUM , ESFLUB ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC RESOLUTION DES EQUATIONS N-S 1 PHASE INCOMPRESSIBLE OU RO VARIABLE
CFONC SUR UN PAS DE TEMPS (CONVECTION/DIFFUSION - PRESSION /CONTINUITE)
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
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE                !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NFML         ! E  !  -> ! NOMBRE DE FAMILLES D ENTITES         !
CARGU ! NPRFML       ! E  !  -> ! NOMBRE DE PROPRIETES DES FAMILLES    !
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE SOMMETS                    !
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! ITERNS       ! E  !  -> ! NUMERO D'ITERATION SUR NAVSTO        !
CARGU ! ICVRGE       ! E  !  -> ! INDICATEUR DE CONVERGENCE DU PNT FIX !
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
CARGU ! IFACLG(2,NFAC! TE !  -  ! TAB ENTIER MULTIGRILLE               !
CARGU ! IRESPR(NCELET! TE !  -  ! TAB ENTIER MULTIGRILLE               !
CARGU ! ISOSTD       ! TE !  -> ! INDICATEUR DE SORTIE STANDARD        !
CARGU !    (NFABOR+1)!    !     !  +NUMERO DE LA FACE DE REFERENCE     !
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
CARGU ! FRCXT(NCELET,! TR !  -> ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU !   3,NPHAS)   !    !     !  HYDROSTATIQUE                       !
CARGU ! TSLAGR       ! TR !  -> ! TERME DE COUPLAGE RETOUR DU          !
CARGU !(NCELET,*)    !    !     !     LAGRANGIEN                       !
CARGU ! TRAVA,XIMPA  ! TR ! <-> ! TABLEAU DE TRAVAIL POUR COUPLAGE     !
CARGU !NCELET,3,NPHAS!    !     ! VITESSE PRESSION PAR POINT FIXE      !
CARGU ! UVWK         ! TR ! <-> ! TABLEAU DE TRAVAIL POUR COUPLAGE U/P !
CARGU !NCELET,3,NPHAS!    !     ! SERT A STOCKER LA VITESSE DE         !
CARGU !              !    !     ! L'ITERATION PRECEDENTE               !
CARGU ! VISCF(NFAC)  ! TR !  -  ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! VISCFI(NFAC) ! TR !  -  ! IDEM VISCF POUR INCREMENTS           !
CARGU ! VISCBI(NFABOR! TR !  -  ! IDEM VISCB POUR INCREMENTS           !
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DAG(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAG(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DRTP(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! TRAV(NCELET,3! TR !  -  ! TABLEAU DE TRAVAIL POUR GRADIENT     !
CARGU ! SMBR  (NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! ROVSDT(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! W1..10(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! DFRCXT(NCELET! TR !  -  ! VARIATION DE FORCE EXTERIEURE        !
CARGU !   3,NPHAS)   !    !     !  GENERANT LA PRESSION HYDROSTATIQUE  !
CARGU ! FRCHY(NCELET ! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU !  NDIM  )     !    !     !  PRESSION HYDROSTATIQUE              !
CARGU ! DFRCHY(NCELET! TR !  -  ! TABLEAU DE TRAVAIL VARIATION DE      !
CARGU !  NDIM  )     !    !     !  PRESSION HYDROSTATIQUE              !
CARGU ! COEFU(NFAB,3)! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! ESFLUM(NFAC) ! TR !  -  ! TABLEAU DE TRAVAIL (IESTOT  )        !
CARGU ! ESFLUB(NFABOR! TR !  -  ! TABLEAU DE TRAVAIL (IESTOT  )        !
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
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "pointe.h"
      INCLUDE "albase.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS  , ITERNS , ICVRGE
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IFACLG(2,NFAC), IRESPR(NCELET)
      INTEGER          ISOSTD(NFABOR+1,NPHAS)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NDIMFB,*)
      DOUBLE PRECISION TSLAGR(NCELET,*)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION FRCXT(NCELET,3,NPHAS)
      DOUBLE PRECISION TRAVA(NCELET,NDIM,NPHAS),XIMPA(NCELET,NDIM,NPHAS)
      DOUBLE PRECISION UVWK(NCELET,NDIM,NPHAS)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION VISCFI(NFAC), VISCBI(NFABOR)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DAG(NCELET), XAG(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET), TRAV(NCELET,3)
      DOUBLE PRECISION SMBR(NCELET), ROVSDT(NCELET)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET), W10(NCELET)
      DOUBLE PRECISION DFRCXT(NCELET,3,NPHAS)
      DOUBLE PRECISION FRCHY(NCELET,NDIM), DFRCHY(NCELET,NDIM)
      DOUBLE PRECISION COEFU(NFABOR,3)
      DOUBLE PRECISION ESFLUM(NFAC), ESFLUB(NFABOR)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          ICCOCG, INC, IEL, IEL1, IEL2, IFAC, IMAX, III
      INTEGER          II    , INOD
      INTEGER          IPHAS , IPH, ISOU, IVAR, IITSM, IGAMM1
      INTEGER          IPRIPH, IUIPH , IVIPH , IWIPH , ICLIPR, ICLIPF
      INTEGER          ICLIUP, ICLIVP, ICLIWP, INIT
      INTEGER          ICLUMA, ICLVMA, ICLWMA
      INTEGER          IFLMAS, IFLMAB, IPCROM, IPBROM
      INTEGER          IFLMS1, IFLMB1, IFLMB0, IISMPH
      INTEGER          NSWRGP, IMLIGP, IWARNP, IMASPE
      INTEGER          IDIMTE, ITENSO
      INTEGER          NBRVAL, IAPPEL, IESCOP, IDTSCA
      INTEGER          IFLINT, IFLBRD, ICOCGV, IFINRA
      INTEGER          NDIRCP, ICPT  , IECRW
      DOUBLE PRECISION RNORM , RNORMA, RNORMI, VITNOR
      DOUBLE PRECISION DTSROM, UNSROM, SURF  , RHOM
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP, XYZMAX(3)
      DOUBLE PRECISION THETAP, XDU, XDV, XDW
      DOUBLE PRECISION RO0IPH, P0IPH, PR0IPH, XXP0 , XYP0 , XZP0
      DOUBLE PRECISION RHOFAC, DTFAC, DDEPX , DDEPY, DDEPZ
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATION
C=======================================================================
C
      IF(IWARNI(IU(1)).GE.1) THEN
        WRITE(NFECRA,1000)
      ENDIF
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IF(NTERUP.GT.1) THEN
C
        DO IPHAS = 1, NPHAS
C
          IUIPH   = IU (IPHAS)
          IVIPH   = IV (IPHAS)
          IWIPH   = IW (IPHAS)
          IPRIPH  = IPR(IPHAS)
          DO ISOU = 1, 3
            IF(ISOU.EQ.1) IVAR = IUIPH
            IF(ISOU.EQ.2) IVAR = IVIPH
            IF(ISOU.EQ.3) IVAR = IWIPH
C     La boucle sur NCELET est une securite au cas
C       ou on utiliserait UVWK par erreur a ITERNS = 1
            DO IEL = 1,NCELET
              UVWK(IEL,ISOU,IPHAS) = RTP(IEL,IVAR)
            ENDDO
          ENDDO
C
C Calcul de la norme L2 de la vitesse
          IF(ITERNS.EQ.1) THEN
            XNRMU0(IPHAS) = 0.D0
            IUIPH   = IU (IPHAS)
            IVIPH   = IV (IPHAS)
            IWIPH   = IW (IPHAS)
            DO IEL = 1, NCEL
              XNRMU0(IPHAS) = XNRMU0(IPHAS) +(RTPA(IEL,IUIPH)**2
     &                                      + RTPA(IEL,IVIPH)**2
     &                                      + RTPA(IEL,IWIPH)**2)
     &                                      * VOLUME(IEL)
            ENDDO
            IF(IRANGP.GE.0) THEN
              CALL PARSOM (XNRMU0(IPHAS))
C             ===========
            ENDIF
            XNRMU0(IPHAS) = SQRT(XNRMU0(IPHAS))
          ENDIF
C
C On assure la periodicite ou le parallelisme de UVWK et la pression
C (cette derniere vaut la pression a l'iteration precedente)
          IF(ITERNS.GT.1) THEN
            IF(IRANGP.GE.0) THEN
              CALL PARCOM (UVWK(1,1,IPHAS))
C             ===========
              CALL PARCOM (UVWK(1,2,IPHAS))
C             ===========
              CALL PARCOM (UVWK(1,3,IPHAS))
C             ===========
              CALL PARCOM (RTPA(1,IPRIPH))
C             ===========
            ENDIF
            IF(IPERIO.EQ.1) THEN
              IDIMTE = 1
              ITENSO = 0
              CALL PERCOM
C             ===========
     &      ( IDIMTE , ITENSO ,
     &        UVWK(1,1,IPHAS),UVWK(1,1,IPHAS),UVWK(1,1,IPHAS),
     &        UVWK(1,2,IPHAS),UVWK(1,2,IPHAS),UVWK(1,2,IPHAS),
     &        UVWK(1,3,IPHAS),UVWK(1,3,IPHAS),UVWK(1,3,IPHAS))
              IDIMTE = 0
              ITENSO = 0
              CALL PERCOM
C             ===========
     &      ( IDIMTE , ITENSO ,
     &        RTPA(1,IPRIPH),RTPA(1,IPRIPH),RTPA(1,IPRIPH),
     &        RTPA(1,IPRIPH),RTPA(1,IPRIPH),RTPA(1,IPRIPH),
     &        RTPA(1,IPRIPH),RTPA(1,IPRIPH),RTPA(1,IPRIPH))
            ENDIF
          ENDIF
C
        ENDDO
C
      ENDIF
C
C
C=======================================================================
C 2.  ETAPE DE PREDICTION DES VITESSES
C=======================================================================
C
      DO IPHAS = 1, NPHAS
C
        IAPPEL = 1
        IUIPH  = IU(IPHAS)
        IFLMAS = IPPROF(IFLUMA(IUIPH))
        IFLMAB = IPPROB(IFLUMA(IUIPH))
        IPH    = IPHAS
C
        CALL PREDUV
C       ===========
     & ( IDEBIA , IDEBRA , IAPPEL ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , ITERNS ,
     &   NCEPDC(IPHAS)   , NCKPDC(IPHAS)   , NCETSM(IPHAS)   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPH    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS))        , IA(IICESM(IPHAS))       ,
     &   IA(IITPSM(IPHAS))        , IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   TSLAGR , COEFA  , COEFB  ,
     &   RA(ICKUPD(IPHAS))        , RA(ISMACE(IPHAS))        ,  FRCXT ,
     &   TRAVA  , XIMPA  , UVWK   , DFRCXT , RA(ITPUCO)      ,  TRAV  ,
     &   VISCF  , VISCB  , VISCFI , VISCBI ,
     &   DAM    , XAM    , DAG    , XAG    ,
     &   DRTP   , SMBR   , ROVSDT ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , W10    , COEFU  ,
     &   RDEVEL , RTUSER , RA     )

      ENDDO
C
C
C --- Sortie si pas de pression continuite (on suppose que
C       la pression est unique, donc pas de iphas dans le test),
C       on met a jour les flux de masse, et on sort
C
      IF( IPRCO.LE.0 ) THEN
C
        DO IPHAS = 1, NPHAS
C
          IUIPH  = IU(IPHAS)
          IVIPH  = IV(IPHAS)
          IWIPH  = IW(IPHAS)
C
          ICLIUP = ICLRTP(IUIPH ,ICOEF)
          ICLIVP = ICLRTP(IVIPH ,ICOEF)
          ICLIWP = ICLRTP(IWIPH ,ICOEF)
C
          IFLMAS = IPPROF(IFLUMA(IUIPH))
          IFLMAB = IPPROB(IFLUMA(IUIPH))
          IPCROM = IPPROC(IROM  (IPHAS))
          IPBROM = IPPROB(IROM  (IPHAS))
C
          INIT   = 1
          INC    = 1
          ICCOCG = 1
          IFLMB0 = 1
          IF (IALE.EQ.1) IFLMB0 = 0
          IISMPH = IISYMP+NFABOR*(IPHAS-1)
          NSWRGP = NSWRGR(IUIPH)
          IMLIGP = IMLIGR(IUIPH)
          IWARNP = IWARNI(IUIPH)
          EPSRGP = EPSRGR(IUIPH)
          CLIMGP = CLIMGR(IUIPH)
          EXTRAP = EXTRAG(IUIPH)
C
          IPH    = IPHAS
C
          IMASPE = 1
C
          CALL INIMAS
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IUIPH  , IVIPH  , IWIPH  , IMASPE , IPH    ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFLMB0 , INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IISMPH) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PROPCE(1,IPCROM), PROPFB(1,IPBROM),
     &   RTP(1,IUIPH) , RTP(1,IVIPH) , RTP(1,IWIPH) ,
     &   COEFA(1,ICLIUP), COEFA(1,ICLIVP), COEFA(1,ICLIWP),
     &   COEFB(1,ICLIUP), COEFB(1,ICLIVP), COEFB(1,ICLIWP),
     &   PROPFA(1,IFLMAS), PROPFB(1,IFLMAB) ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C
        ENDDO
C
C     En ALE on doit rajouter la composante en vitesse de maillage
        IF (IALE.EQ.1) THEN
C
          ICLUMA = ICLRTP(IUMA ,ICOEF)
          ICLVMA = ICLRTP(IVMA ,ICOEF)
          ICLWMA = ICLRTP(IWMA ,ICOEF)
C
C     On change de signe car on veut l'oppose de la vitesse de maillage
C       aux faces
          DO IEL = 1, NCELET
            RTP(IEL,IUMA) = -RTP(IEL,IUMA)
            RTP(IEL,IVMA) = -RTP(IEL,IVMA)
            RTP(IEL,IWMA) = -RTP(IEL,IWMA)
          ENDDO
          DO IFAC = 1, NFABOR
            COEFA(IFAC,ICLUMA) = -COEFA(IFAC,ICLUMA)
            COEFA(IFAC,ICLVMA) = -COEFA(IFAC,ICLVMA)
            COEFA(IFAC,ICLWMA) = -COEFA(IFAC,ICLWMA)
          ENDDO
C
C     One temporary array needed for internal faces, in case some internal vertices
C       are moved directly by the user
          IFLINT = IDEBRA
          IFINRA = IFLINT + NFAC
C
          CALL RASIZE('NAVSTO',IFINRA)
C         ===========
C
          DO IPHAS = 1, NPHAS
C
            IFLMAS = IPPROF(IFLUMA(IU(IPHAS)))
            IFLMAB = IPPROB(IFLUMA(IU(IPHAS)))
            IPCROM = IPPROC(IROM  (IPHAS))
            IPBROM = IPPROB(IROM  (IPHAS))
C
            INIT   = 0
            INC    = 1
            ICCOCG = 1
            IFLMB0 = 1
            NSWRGP = NSWRGR(IUMA )
            IMLIGP = IMLIGR(IUMA )
            IWARNP = IWARNI(IUMA )
            EPSRGP = EPSRGR(IUMA )
            CLIMGP = CLIMGR(IUMA )
            EXTRAP = EXTRAG(IUMA )
C
            IPH    = IPHAS
C
            IMASPE = 1
C
            DO IFAC = 1, NFAC
              RA(IFLINT+IFAC-1) = 0.D0
            ENDDO
C
            CALL INIMAS
C           ===========
     & ( IDEBIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IUIPH  , IVIPH  , IWIPH  , IMASPE , IPH    ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFLMB0 , INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IISMPH) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PROPCE(1,IPCROM), PROPFB(1,IPBROM),
     &   RTP(1,IUMA )    , RTP(1,IVMA )    , RTP(1,IWMA)     ,
     &   COEFA(1,ICLUMA), COEFA(1,ICLVMA), COEFA(1,ICLWMA),
     &   COEFB(1,ICLUMA), COEFB(1,ICLVMA), COEFB(1,ICLWMA),
     &   RA(IFLINT)     , PROPFB(1,IFLMAB) ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C
          ENDDO
C
          DO IEL = 1, NCELET
            RTP(IEL,IUMA) = -RTP(IEL,IUMA)
            RTP(IEL,IVMA) = -RTP(IEL,IVMA)
            RTP(IEL,IWMA) = -RTP(IEL,IWMA)
          ENDDO
          DO IFAC = 1, NFABOR
            COEFA(IFAC,ICLUMA) = -COEFA(IFAC,ICLUMA)
            COEFA(IFAC,ICLVMA) = -COEFA(IFAC,ICLVMA)
            COEFA(IFAC,ICLWMA) = -COEFA(IFAC,ICLWMA)
          ENDDO
C
          DO IFAC = 1, NFAC
            IECRW = 0
            DDEPX = 0.D0
            DDEPY = 0.D0
            DDEPZ = 0.D0
            ICPT  = 0
            DO II = IPNFAC(IFAC),IPNFAC(IFAC+1)-1
              INOD = NODFAC(II)
              IF (IA(IIMPAL+INOD-1).EQ.0) IECRW = IECRW + 1
              ICPT = ICPT + 1
              DDEPX = DDEPX + RA(IDEPAL       +INOD-1)
     &               +RA(IXYZN0+(INOD-1)*NDIM  )-XYZNOD(1,INOD)
              DDEPY = DDEPY + RA(IDEPAL+NNOD  +INOD-1)
     &               +RA(IXYZN0+(INOD-1)*NDIM+1)-XYZNOD(2,INOD)
              DDEPZ = DDEPZ + RA(IDEPAL+2*NNOD+INOD-1)
     &               +RA(IXYZN0+(INOD-1)*NDIM+2)-XYZNOD(3,INOD)
            ENDDO
C     If all the face vertices have imposed displacement, w is evaluated from
C       this displacement
            IF (IECRW.EQ.0) THEN
              IEL1 = IFACEL(1,IFAC)
              IEL2 = IFACEL(2,IFAC)
              DTFAC = 0.5D0*(DT(IEL1) + DT(IEL2))
              RHOFAC = 0.5D0*(PROPCE(IEL1,IPCROM) + PROPCE(IEL2,IPCROM))
              PROPFA(IFAC,IFLMAS) = PROPFA(IFAC,IFLMAS) - RHOFAC*(
     &                              DDEPX*SURFAC(1,IFAC)
     &                             +DDEPY*SURFAC(2,IFAC)
     &                             +DDEPZ*SURFAC(3,IFAC) )/DTFAC/ICPT
C     Else w is calculated from the cell-centre mesh velocity
            ELSE
              PROPFA(IFAC,IFLMAS) = PROPFA(IFAC,IFLMAS)
     &                            + RA(IFLINT+IFAC-1)
            ENDIF
          ENDDO
        ENDIF
C
        RETURN
C
      ENDIF
C
C
C=======================================================================
C 3.  ETAPE DE PRESSION/CONTINUITE ( VITESSE/PRESSION )
C=======================================================================
C
      IF(IWARNI(IU(1)).GE.1) THEN
        WRITE(NFECRA,1200)
      ENDIF
C
C On n'appelle resolp qu'une seule fois, pour LA phase qu'il faut
C
      IPHAS = 1
      IPH   = IPHAS
C
C --- Pas de temps scalaire ou pas
      IDTSCA = 0
      IF ((IPUCOU.EQ.1).OR.(NCPDCT(IPHAS).GT.0)) IDTSCA = 1
C
      CALL RESOLP
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCEPDC(IPHAS)   , NCKPDC(IPHAS)   , NCETSM(IPHAS)   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPH    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS))        , IA(IICESM(IPHAS))       ,
     &   IA(IITPSM(IPHAS))        , ISOSTD , IFACLG , IRESPR , IDTSCA ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RA(ICKUPD(IPHAS))        , RA(ISMACE(IPHAS))        ,
     &   FRCXT  , DFRCXT , RA(ITPUCO)      , TRAV   ,
     &   VISCF  , VISCB  , VISCFI , VISCBI ,
     &   DAM    , XAM    , DAG    , XAG    ,
     &   DRTP   , SMBR   , ROVSDT , TSLAGR ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , FRCHY  , DFRCHY , COEFU  , TRAVA ,
     &   RDEVEL , RTUSER , RA     )
C
C
C Si on est en polyphasique, il faudra de toutes maniere modifier
C  tout ceci. Pour le moment, on se contente de mettre a jour les
C  flux de masse des phases 2 a N en postulant qu'ils sont tous egaux
C  au flux de masse de la phase 1 ...
C Ca ne sert qu'aux tests bien sur
C
      IF(NPHAS.GT.1) THEN
        IFLMS1 = IPPROF(IFLUMA(IU(1    )))
        IFLMB1 = IPPROB(IFLUMA(IU(1    )))
        DO IPHAS = 2, NPHAS
          IFLMAS = IPPROF(IFLUMA(IU(IPHAS)))
          IFLMAB = IPPROB(IFLUMA(IU(IPHAS)))
          DO IFAC = 1, NFAC
            PROPFA(IFAC,IFLMAS) = PROPFA(IFAC,IFLMS1)
          ENDDO
          DO IFAC = 1, NFABOR
            PROPFB(IFAC,IFLMAB) = PROPFB(IFAC,IFLMB1)
          ENDDO
        ENDDO
      ENDIF
C
C=======================================================================
C 4.  REACTUALISATION DU CHAMP DE VITESSE
C=======================================================================
C
C
      DO IPHAS = 1, NPHAS
C
        IPRIPH = IPR(IPHAS)
        IUIPH  = IU(IPHAS)
        IVIPH  = IV(IPHAS)
        IWIPH  = IW(IPHAS)
C
        ICLIPR = ICLRTP(IPRIPH,ICOEF)
        ICLIPF = ICLRTP(IPRIPH,ICOEFF)
        ICLIUP = ICLRTP(IUIPH ,ICOEF)
        ICLIVP = ICLRTP(IVIPH ,ICOEF)
        ICLIWP = ICLRTP(IWIPH ,ICOEF)
C
        IFLMAS = IPPROF(IFLUMA(IUIPH))
        IFLMAB = IPPROB(IFLUMA(IUIPH))
        IPCROM = IPPROC(IROM  (IPHAS))
        IPBROM = IPPROB(IROM  (IPHAS))
        IISMPH = IISYMP+NFABOR*(IPHAS-1)
C
C
C
C       IREVMC = 0 : Methode standard (pas par moindres carres) : on
C                      ajoute un gradient d'increment de pression standard
C                      a la vitesse predite opur obtenir la vitesse corrigee
C
C       IREVMC = 1 : On applique la methode par moindres carres a
C                      l'ecart entre le flux de masse predit et le flux
C                      de masse actualise,
C                      c'est-a-dire au gradient d'increment de pression
C                    On ajoute la grandeur obtenue aux cellules a la vitesse
C                      predite pour obtenir la vitesse actualisee
C                    Cette methode correspond a IREVMC = 0 avec
C                      gradient par moindres carres IMRGRA=1 dans la
C                      reactualisation des vitesses.
C
C       IREVMC = 2 : On applique la methode par moindres carres au
C                      flux de masse actualise
C                      pour obtenir la vitesse actualisee
C                    Cette methode correspond a la methode RT0.
C
C       La methode IREVMC = 2 semble plus "diffusive", mais semble aussi la
C         seule issue pour certains ecoulements atmospheriques de mercure.
C       La methode IREVMC = 1 semble ne pas trop "diffuser", avec un
C         gain du a l'utilisation du gradient moindres carres. Elle
C         se rapproche beaucoup de IREVMC=0.
C
C
        IF( IREVMC(IPHAS).EQ.1 ) THEN
C
C     On a besoin de trois tableaux de travail
          IFLINT = IDEBRA
          IFLBRD = IFLINT + NFAC
          ICOCGV = IFLBRD + NFABOR
          IFINRA = ICOCGV + NCELET*9
C
          CALL RASIZE('NAVSTO',IFINRA)
C         ===========
C
C     on ote la partie en u-predit dans le flux de masse final,
C     on projete des faces vers le centre, puis on rajoute u-predit.
          INIT   = 1
          INC    = 1
          ICCOCG = 1
          IFLMB0 = 1
          IF (IALE.EQ.1) IFLMB0 = 0
          NSWRGP = NSWRGR(IUIPH )
          IMLIGP = IMLIGR(IUIPH )
          IWARNP = IWARNI(IUIPH )
          EPSRGP = EPSRGR(IUIPH )
          CLIMGP = CLIMGR(IUIPH )
          EXTRAP = EXTRAG(IUIPH )
C
          IPH  = IPHAS
C
          IMASPE = 1
C
          CALL INIMAS
C         ===========
     & ( IDEBIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IUIPH  , IVIPH  , IWIPH  , IMASPE , IPH    ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFLMB0 , INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IISMPH) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PROPCE(1,IPCROM), PROPFB(1,IPBROM),
     &   RTP(1,IUIPH)    , RTP(1,IVIPH)    , RTP(1,IWIPH)    ,
     &   COEFA(1,ICLIUP), COEFA(1,ICLIVP), COEFA(1,ICLIWP),
     &   COEFB(1,ICLIUP), COEFB(1,ICLIVP), COEFB(1,ICLIWP),
     &   RA(IFLINT), RA(IFLBRD),
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C
          DO IFAC = 1, NFAC
            RA(IFLINT+IFAC-1) = PROPFA(IFAC,IFLMAS) - RA(IFLINT+IFAC-1)
          ENDDO
          DO IFAC = 1, NFABOR
            RA(IFLBRD+IFAC-1) = PROPFB(IFAC,IFLMAB) - RA(IFLBRD+IFAC-1)
          ENDDO
C
          CALL RECVMC
C         ===========
     & ( IDEBIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PROPCE(1,IPCROM), RA(IFLINT)      , RA(IFLBRD)      ,
     &   W1     , W2     , W3     ,
     &   W4     , W5     , W6     , RA(ICOCGV)      ,
     &   RDEVEL , RTUSER , RA     )
C
          DO IEL = 1, NCEL
            RTP(IEL,IUIPH) = RTP(IEL,IUIPH) + W1(IEL)
            RTP(IEL,IVIPH) = RTP(IEL,IVIPH) + W2(IEL)
            RTP(IEL,IWIPH) = RTP(IEL,IWIPH) + W3(IEL)
          ENDDO
C
        ELSEIF( IREVMC(IPHAS).EQ.2 ) THEN
C
C     On calcule la vitesse corrigee directement a partir du flux de masse
C       corrige.
C     On a besoin de trois tableaux de travail
          ICOCGV = IDEBRA
          IFINRA = ICOCGV + NCELET*9
C
          CALL RASIZE('NAVSTO',IFINRA)
C         ===========
C
          CALL RECVMC
C         ===========
     & ( IDEBIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PROPCE(1,IPCROM), PROPFA(1,IFLMAS), PROPFB(1,IFLMAB),
     &   RTP(1,IUIPH), RTP(1,IVIPH), RTP(1,IWIPH),
     &   W4     , W5     , W6     , RA(ICOCGV)   ,
     &   RDEVEL , RTUSER , RA     )
C
C
        ELSE
C
C     On corrige la vitesse predite par le gradient cellule de
C       l'increment de pression
C
C     GRADIENT DE L'INCREMENT TOTAL DE PRESSION
C
          IF (IDTVAR.LT.0) THEN
            DO IEL = 1, NCEL
              DRTP(IEL) = (RTP(IEL,IPRIPH) -RTPA(IEL,IPRIPH))
     &                   /RELAXV(IPRIPH)
            ENDDO
          ELSE
            DO IEL = 1, NCEL
              DRTP(IEL) = RTP(IEL,IPRIPH) -RTPA(IEL,IPRIPH)
            ENDDO
          ENDIF
C
C --->    TRAITEMENT DU PARALLELISME
C
          IF(IRANGP.GE.0) CALL PARCOM (DRTP)
C                         ===========
C
C ---> TRAITEMENT DE LA PERIODICITE
C
          IF(IPERIO.EQ.1) THEN
            IDIMTE = 0
            ITENSO = 0
            CALL PERCOM
C           ===========
     &      ( IDIMTE , ITENSO ,
     &        DRTP   , DRTP   , DRTP  ,
     &        DRTP   , DRTP   , DRTP  ,
     &        DRTP   , DRTP   , DRTP  )
          ENDIF
C
C
          ICCOCG = 1
          INC = 0
          IF (IPHYDR.EQ.1) INC = 1
          NSWRGP = NSWRGR(IPRIPH)
          IMLIGP = IMLIGR(IPRIPH)
          IWARNP = IWARNI(IPRIPH)
          EPSRGP = EPSRGR(IPRIPH)
          CLIMGP = CLIMGR(IPRIPH)
          EXTRAP = EXTRAG(IPRIPH)
C
          CALL GRDCEL
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPRIPH , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDR ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DFRCXT(1,1,IPHAS),DFRCXT(1,2,IPHAS),DFRCXT(1,3,IPHAS),
     &   DRTP   , COEFA(1,ICLIPF) , COEFB(1,ICLIPR)  ,
     &   TRAV(1,1)       , TRAV(1,2)       , TRAV(1,3) ,
C        ---------         ---------         ---------
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
C     REACTUALISATION DU CHAMP DE VITESSES
C
          THETAP = THETAV(IPRIPH)
          IF (IPHYDR.EQ.0) THEN
            IF (IDTSCA.EQ.0) THEN
              DO IEL = 1, NCEL
                DTSROM = -THETAP*DT(IEL)/PROPCE(IEL,IPCROM)
                RTP(IEL,IUIPH) = RTP(IEL,IUIPH)+DTSROM*TRAV(IEL,1)
                RTP(IEL,IVIPH) = RTP(IEL,IVIPH)+DTSROM*TRAV(IEL,2)
                RTP(IEL,IWIPH) = RTP(IEL,IWIPH)+DTSROM*TRAV(IEL,3)
              ENDDO
            ELSE
              DO IEL = 1, NCEL
                UNSROM = -THETAP/PROPCE(IEL,IPCROM)
                III = ITPUCO-1+IEL
                RTP(IEL,IUIPH) = RTP(IEL,IUIPH)
     &                          +UNSROM*RA(III         )*TRAV(IEL,1)
                RTP(IEL,IVIPH) = RTP(IEL,IVIPH)
     &                          +UNSROM*RA(III+NCELET  )*TRAV(IEL,2)
                RTP(IEL,IWIPH) = RTP(IEL,IWIPH)
     &                          +UNSROM*RA(III+2*NCELET)*TRAV(IEL,3)
              ENDDO
            ENDIF
          ELSE
            IF (IDTSCA.EQ.0) THEN
              DO IEL = 1, NCEL
                DTSROM = THETAP*DT(IEL)/PROPCE(IEL,IPCROM)
                RTP(IEL,IUIPH) = RTP(IEL,IUIPH)
     &               +DTSROM*(DFRCXT(IEL,1,IPHAS)-TRAV(IEL,1) )
                RTP(IEL,IVIPH) = RTP(IEL,IVIPH)
     &               +DTSROM*(DFRCXT(IEL,2,IPHAS)-TRAV(IEL,2) )
                RTP(IEL,IWIPH) = RTP(IEL,IWIPH)
     &               +DTSROM*(DFRCXT(IEL,3,IPHAS)-TRAV(IEL,3) )
              ENDDO
            ELSE
              DO IEL = 1, NCEL
                UNSROM = THETAP/PROPCE(IEL,IPCROM)
                III = ITPUCO-1+IEL
                RTP(IEL,IUIPH) = RTP(IEL,IUIPH)
     &               +UNSROM*RA(III         )
     &               *(DFRCXT(IEL,1,IPHAS)-TRAV(IEL,1) )
                RTP(IEL,IVIPH) = RTP(IEL,IVIPH)
     &               +UNSROM*RA(III+NCELET  )
     &               *(DFRCXT(IEL,2,IPHAS)-TRAV(IEL,2) )
                RTP(IEL,IWIPH) = RTP(IEL,IWIPH)
     &               +UNSROM*RA(III+2*NCELET)
     &               *(DFRCXT(IEL,3,IPHAS)-TRAV(IEL,3) )
              ENDDO
            ENDIF
C     mise a jour des forces exterieures pour le calcul des gradients
            DO IEL=1,NCEL
              FRCXT(IEL,1,IPHAS) = FRCXT(IEL,1,IPHAS)
     &             + DFRCXT(IEL,1,IPHAS)
              FRCXT(IEL,2,IPHAS) = FRCXT(IEL,2,IPHAS)
     &             + DFRCXT(IEL,2,IPHAS)
              FRCXT(IEL,3,IPHAS) = FRCXT(IEL,3,IPHAS)
     &             + DFRCXT(IEL,3,IPHAS)
            ENDDO
            IF(IRANGP.GE.0) THEN
              CALL PARCOM (FRCXT(1,1,IPHAS))
C             ===========
              CALL PARCOM (FRCXT(1,2,IPHAS))
C             ===========
              CALL PARCOM (FRCXT(1,3,IPHAS))
C             ===========
            ENDIF
            IF(IPERIO.EQ.1) THEN
              IDIMTE = 1
              ITENSO = 0
              CALL PERCOM
C             ===========
     &  ( IDIMTE , ITENSO ,
     &    FRCXT(1,1,IPHAS),FRCXT(1,1,IPHAS),FRCXT(1,1,IPHAS),
     &    FRCXT(1,2,IPHAS),FRCXT(1,2,IPHAS),FRCXT(1,2,IPHAS),
     &    FRCXT(1,3,IPHAS),FRCXT(1,3,IPHAS),FRCXT(1,3,IPHAS) )
            ENDIF
c     mise a jour des Dirichlets de pression en sortie dans COEFA
            ICLIPR = ICLRTP(IPRIPH,ICOEF)
            ICLIPF = ICLRTP(IPRIPH,ICOEFF)
            DO IFAC = 1,NFABOR
              IF (ISOSTD(IFAC,IPHAS).EQ.1)
     &             COEFA(IFAC,ICLIPR) = COEFA(IFAC,ICLIPR)
     &             + COEFA(IFAC,ICLIPF)
            ENDDO
          ENDIF
C
        ENDIF
C
      ENDDO
C
C
C     Ajout de la vitesse de maillage dans le flux convectif en ALE
      IF (IALE.EQ.1) THEN
C
        ICLUMA = ICLRTP(IUMA ,ICOEF)
        ICLVMA = ICLRTP(IVMA ,ICOEF)
        ICLWMA = ICLRTP(IWMA ,ICOEF)
C
C     On change de signe car on veut l'oppose de la vitesse de maillage
C       aux faces
        DO IEL = 1, NCELET
          RTP(IEL,IUMA) = -RTP(IEL,IUMA)
          RTP(IEL,IVMA) = -RTP(IEL,IVMA)
          RTP(IEL,IWMA) = -RTP(IEL,IWMA)
        ENDDO
        DO IFAC = 1, NFABOR
          COEFA(IFAC,ICLUMA) = -COEFA(IFAC,ICLUMA)
          COEFA(IFAC,ICLVMA) = -COEFA(IFAC,ICLVMA)
          COEFA(IFAC,ICLWMA) = -COEFA(IFAC,ICLWMA)
        ENDDO
C
C     One temporary array needed for internal faces, in case some internal vertices
C       are moved directly by the user
          IFLINT = IDEBRA
          IFINRA = IFLINT + NFAC
C
          CALL RASIZE('NAVSTO',IFINRA)
C         ===========
C
        DO IPHAS = 1, NPHAS
C
          IFLMAS = IPPROF(IFLUMA(IU(IPHAS)))
          IFLMAB = IPPROB(IFLUMA(IU(IPHAS)))
          IPCROM = IPPROC(IROM  (IPHAS))
          IPBROM = IPPROB(IROM  (IPHAS))
C
          INIT   = 0
          INC    = 1
          ICCOCG = 1
          IFLMB0 = 1
          NSWRGP = NSWRGR(IUMA )
          IMLIGP = IMLIGR(IUMA )
          IWARNP = IWARNI(IUMA )
          EPSRGP = EPSRGR(IUMA )
          CLIMGP = CLIMGR(IUMA )
          EXTRAP = EXTRAG(IUMA )
C
          IPH    = IPHAS
C
          IMASPE = 1
C
          DO IFAC = 1, NFAC
            RA(IFLINT+IFAC-1) = 0.D0
          ENDDO
C
          CALL INIMAS
C         ===========
     & ( IDEBIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IUIPH  , IVIPH  , IWIPH  , IMASPE , IPH    ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFLMB0 , INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IISMPH) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PROPCE(1,IPCROM), PROPFB(1,IPBROM),
     &   RTP(1,IUMA )    , RTP(1,IVMA )    , RTP(1,IWMA)     ,
     &   COEFA(1,ICLUMA), COEFA(1,ICLVMA), COEFA(1,ICLWMA),
     &   COEFB(1,ICLUMA), COEFB(1,ICLVMA), COEFB(1,ICLWMA),
     &   RA(IFLINT)     , PROPFB(1,IFLMAB) ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C
        ENDDO
C
        DO IEL = 1, NCELET
          RTP(IEL,IUMA) = -RTP(IEL,IUMA)
          RTP(IEL,IVMA) = -RTP(IEL,IVMA)
          RTP(IEL,IWMA) = -RTP(IEL,IWMA)
        ENDDO
        DO IFAC = 1, NFABOR
          COEFA(IFAC,ICLUMA) = -COEFA(IFAC,ICLUMA)
          COEFA(IFAC,ICLVMA) = -COEFA(IFAC,ICLVMA)
          COEFA(IFAC,ICLWMA) = -COEFA(IFAC,ICLWMA)
        ENDDO
C
        DO IFAC = 1, NFAC
          IECRW = 0
          DDEPX = 0.D0
          DDEPY = 0.D0
          DDEPZ = 0.D0
          ICPT  = 0
          DO II = IPNFAC(IFAC),IPNFAC(IFAC+1)-1
            INOD = NODFAC(II)
            IF (IA(IIMPAL+INOD-1).EQ.0) IECRW = IECRW + 1
            ICPT = ICPT + 1
            DDEPX = DDEPX + RA(IDEPAL       +INOD-1)
     &             +RA(IXYZN0+(INOD-1)*NDIM  )-XYZNOD(1,INOD)
            DDEPY = DDEPY + RA(IDEPAL+NNOD  +INOD-1)
     &             +RA(IXYZN0+(INOD-1)*NDIM+1)-XYZNOD(2,INOD)
            DDEPZ = DDEPZ + RA(IDEPAL+2*NNOD+INOD-1)
     &             +RA(IXYZN0+(INOD-1)*NDIM+2)-XYZNOD(3,INOD)
          ENDDO
C     If all the face vertices have imposed displacement, w is evaluated from
C       this displacement
          IF (IECRW.EQ.0) THEN
            IEL1 = IFACEL(1,IFAC)
            IEL2 = IFACEL(2,IFAC)
            DTFAC = 0.5D0*(DT(IEL1) + DT(IEL2))
            RHOFAC = 0.5D0*(PROPCE(IEL1,IPCROM) + PROPCE(IEL2,IPCROM))
            PROPFA(IFAC,IFLMAS) = PROPFA(IFAC,IFLMAS) - RHOFAC*(
     &                            DDEPX*SURFAC(1,IFAC)
     &                           +DDEPY*SURFAC(2,IFAC)
     &                           +DDEPZ*SURFAC(3,IFAC) )/DTFAC/ICPT
C     Else w is calculated from the cell-centre mesh velocity
          ELSE
            PROPFA(IFAC,IFLMAS) = PROPFA(IFAC,IFLMAS)
     &                          + RA(IFLINT+IFAC-1)
          ENDIF
        ENDDO
C
      ENDIF
C
C=======================================================================
C 5.  CALCUL D'UN ESTIMATEUR D'ERREUR DE L'ETAPE DE CORRECTION ET TOTAL
C=======================================================================
C
C
      DO IPHAS = 1, NPHAS
C
        IF(IESCAL(IESCOR,IPHAS).GT.0.OR.IESCAL(IESTOT,IPHAS).GT.0) THEN
C
C ---> REPERAGE DES VARIABLES
C
          IPRIPH = IPR(IPHAS)
          IUIPH  = IU(IPHAS)
          IVIPH  = IV(IPHAS)
          IWIPH  = IW(IPHAS)
C
          ICLIUP = ICLRTP(IUIPH ,ICOEF)
          ICLIVP = ICLRTP(IVIPH ,ICOEF)
          ICLIWP = ICLRTP(IWIPH ,ICOEF)
C
          IPCROM = IPPROC(IROM  (IPHAS))
          IPBROM = IPPROB(IROM  (IPHAS))
          IISMPH = IISYMP+NFABOR*(IPHAS-1)
C
C
C
C ---> ECHANGE DES VITESSES ET PRESSION EN PERIODICITE ET PARALLELISME
C
C    Pour les estimateurs IESCOR et IESTOT, la vitesse doit etre echangee.
C
C    Pour l'estimateur IESTOT, la pression doit etre echangee aussi.
C
C    Cela ne remplace pas l'echange du debut de pas de temps
C     a cause de usproj qui vient plus tard et des calculs suite)
C
C
C --- Vitesse
C
          IF(IRANGP.GE.0) THEN
            CALL PARCOM (RTP(1,IUIPH ))
C           ===========
            CALL PARCOM (RTP(1,IVIPH ))
C           ===========
            CALL PARCOM (RTP(1,IWIPH ))
C           ===========
          ENDIF
C
          IF(IPERIO.EQ.1) THEN
            IDIMTE = 1
            ITENSO = 0
            CALL PERCOM
C           ===========
     &         ( IDIMTE , ITENSO ,
     &           RTP(1,IUIPH), RTP(1,IUIPH), RTP(1,IUIPH),
     &           RTP(1,IVIPH), RTP(1,IVIPH), RTP(1,IVIPH),
     &           RTP(1,IWIPH), RTP(1,IWIPH), RTP(1,IWIPH))
          ENDIF
C
C
C  -- Pression
C
          IF(IESCAL(IESTOT,IPHAS).GT.0) THEN
C
            IF(IRANGP.GE.0) THEN
              CALL PARCOM (RTP(1,IPRIPH))
C             ===========
            ENDIF
C
            IF(IPERIO.EQ.1) THEN
              IDIMTE = 0
              ITENSO = 0
              CALL PERCOM
C             ===========
     &           ( IDIMTE , ITENSO ,
     &           RTP(1,IPRIPH), RTP(1,IPRIPH), RTP(1,IPRIPH),
     &           RTP(1,IPRIPH), RTP(1,IPRIPH), RTP(1,IPRIPH),
     &           RTP(1,IPRIPH), RTP(1,IPRIPH), RTP(1,IPRIPH))
            ENDIF
C
          ENDIF
C
C
C ---> CALCUL DU FLUX DE MASSE DEDUIT DE LA VITESSE REACTUALISEE
C
          INIT   = 1
          INC    = 1
          ICCOCG = 1
          IFLMB0 = 1
          IF (IALE.EQ.1) IFLMB0 = 0
          NSWRGP = NSWRGR(IUIPH )
          IMLIGP = IMLIGR(IUIPH )
          IWARNP = IWARNI(IUIPH )
          EPSRGP = EPSRGR(IUIPH )
          CLIMGP = CLIMGR(IUIPH )
          EXTRAP = EXTRAG(IUIPH )
C
          IPH    = IPHAS
C
          IMASPE = 1
C
          CALL INIMAS
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   IUIPH  , IVIPH  , IWIPH  , IMASPE , IPH    ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFLMB0 , INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP ,
     &   IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IA(IISMPH) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PROPCE(1,IPCROM), PROPFB(1,IPBROM),
     &   RTP(1,IUIPH)    , RTP(1,IVIPH)    , RTP(1,IWIPH)    ,
     &   COEFA(1,ICLIUP), COEFA(1,ICLIVP), COEFA(1,ICLIWP),
     &   COEFB(1,ICLIUP), COEFB(1,ICLIVP), COEFB(1,ICLIWP),
     &   ESFLUM , ESFLUB ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C
C
C ---> CALCUL DE L'ESTIMATEUR CORRECTION : DIVERGENCE DE ROM * U (N + 1)
C                                          - GAMMA
C
          IF(IESCAL(IESCOR,IPHAS).GT.0) THEN
            INIT = 1
            CALL DIVMAS(NCELET,NCEL,NFAC,NFABOR,INIT,NFECRA,
     &                                   IFACEL,IFABOR,ESFLUM,ESFLUB,W1)
C
            IF (NCETSM(IPHAS).GT.0) THEN
              IGAMM1 = ISMACE(IPHAS)+(IPRIPH-1)*NCETSM(IPHAS)-1
              DO IITSM = 1, NCETSM(IPHAS)
                IEL = IA(IICESM(IPHAS)+IITSM-1)
                W1(IEL) = W1(IEL)-VOLUME(IEL)*RA(IGAMM1+IITSM)
              ENDDO
            ENDIF
C
            IF(IESCAL(IESCOR,IPHAS).EQ.2) THEN
              IESCOP = IPPROC(IESTIM(IESCOR,IPHAS))
              DO IEL = 1, NCEL
                PROPCE(IEL,IESCOP) =  ABS(W1(IEL))
              ENDDO
            ELSEIF(IESCAL(IESCOR,IPHAS).EQ.1) THEN
              IESCOP = IPPROC(IESTIM(IESCOR,IPHAS))
              DO IEL = 1, NCEL
                PROPCE(IEL,IESCOP) =  ABS(W1(IEL)) / VOLUME(IEL)
              ENDDO
            ENDIF
          ENDIF
C
C
C ---> CALCUL DE L'ESTIMATEUR TOTAL
C
          IF(IESCAL(IESTOT,IPHAS).GT.0) THEN
C
C   INITIALISATION DE TRAV AVEC LE TERME INSTATIONNAIRE
C
            DO IEL = 1, NCEL
              TRAV(IEL,1) = PROPCE(IEL,IPCROM) * VOLUME(IEL) *
     &             ( RTPA(IEL,IUIPH)- RTP(IEL,IUIPH) )/DT(IEL)
              TRAV(IEL,2) = PROPCE(IEL,IPCROM) * VOLUME(IEL) *
     &             ( RTPA(IEL,IVIPH)- RTP(IEL,IVIPH) )/DT(IEL)
              TRAV(IEL,3) = PROPCE(IEL,IPCROM) * VOLUME(IEL) *
     &             ( RTPA(IEL,IWIPH)- RTP(IEL,IWIPH) )/DT(IEL)
            ENDDO
C
C   APPEL A PREDUV AVEC RTP ET RTP AU LIEU DE RTP ET RTPA
C                  AVEC LE FLUX DE MASSE RECALCULE
            IAPPEL = 2
            IPH    = IPHAS
            CALL PREDUV
C           ===========
     & ( IDEBIA , IDEBRA , IAPPEL ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , ITERNS ,
     &   NCEPDC(IPHAS)   , NCKPDC(IPHAS)   , NCETSM(IPHAS)   ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPH    ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IA(IICEPD(IPHAS))        , IA(IICESM(IPHAS))       ,
     &   IA(IITPSM(IPHAS))        , IFACLG , IRESPR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTP    , PROPCE , PROPFA , PROPFB ,
     &   ESFLUM , ESFLUB ,
     &   TSLAGR , COEFA  , COEFB  ,
     &   RA(ICKUPD(IPHAS))        , RA(ISMACE(IPHAS))        , FRCXT  ,
     &   TRAVA  , XIMPA  , UVWK   , DFRCXT , RA(ITPUCO)      , TRAV   ,
     &   VISCF  , VISCB  , VISCFI , VISCBI ,
     &   DAM    , XAM    , DAG    , XAG    ,
     &   DRTP   , SMBR   , ROVSDT ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   W7     , W8     , W9     , W10    , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C
          ENDIF
C
        ENDIF
C
      ENDDO
C
C=======================================================================
C 6.  TRAITEMENT DU POINT FIXE SUR LE SYSTEME VITESSE/PRESSION
C=======================================================================
C
      IF(NTERUP.GT.1) THEN
C TEST DE CONVERGENCE DE L'ALGORITHME ITERATIF
C On initialise ICVRGE a 1 et on le met a 0 si une des phases n'est
C pas convergee
C
        ICVRGE = 1
C
        DO IPHAS = 1,NPHAS
          XNRMU(IPHAS) = 0.D0
          DO IEL = 1,NCEL
            XDU = RTP(IEL,IUIPH) - UVWK(IEL,1,IPHAS)
            XDV = RTP(IEL,IVIPH) - UVWK(IEL,2,IPHAS)
            XDW = RTP(IEL,IWIPH) - UVWK(IEL,3,IPHAS)
            XNRMU(IPHAS) = XNRMU(IPHAS) +(XDU**2 + XDV**2 + XDW**2)
     &                                  * VOLUME(IEL)
          ENDDO
C --->    TRAITEMENT DU PARALLELISME
C
          IF(IRANGP.GE.0) CALL PARSOM (XNRMU(IPHAS))
C                           ===========
          XNRMU(IPHAS) = SQRT(XNRMU(IPHAS))
C
C Indicateur de convergence du point fixe
          IF(XNRMU(IPHAS).GE.EPSUP(IPHAS)*XNRMU0(IPHAS)) ICVRGE = 0
C
        ENDDO

C
      ENDIF
C
C ---> RECALAGE DE LA PRESSION SUR UNE PRESSION A MOYENNE NULLE
C  On recale si on n'a pas de Dirichlet. Or le nombre de Dirichlets
C  calcule dans typecl.F est NDIRCL si IDIRCL=1 et NDIRCL-1 si IDIRCL=0
C  (ISTAT vaut toujours 0 pour la pression)
C
      DO IPHAS = 1, NPHAS
        IPRIPH  = IPR(IPHAS)
        IF (IDIRCL(IPR(IPHAS)).EQ.1) THEN
          NDIRCP = NDIRCL(IPR(IPHAS))
        ELSE
          NDIRCP = NDIRCL(IPR(IPHAS))-1
        ENDIF
        IF(NDIRCP.LE.0) THEN
          CALL PRMOY0
C         ===========
     &    ( IDEBIA , IDEBRA ,
     &      NCELET , NCEL   , NFAC   , NFABOR ,
     &      NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &      IPHAS  , IDEVEL , ITUSER , IA     ,
     &      VOLUME , RTP(1,IPRIPH) ,
     &      RDEVEL , RTUSER , RA     )
        ENDIF
C
C Calcul de la pression totale IPRTOT : (definie comme propriete )
C En compressible, la pression resolue est deja la pression totale
C
        IF (IPPMOD(ICOMPF).LT.0) THEN
          RO0IPH = RO0  (IPHAS)
          P0IPH  = P0   (IPHAS)
          PR0IPH = PRED0(IPHAS)
          XXP0   = XYZP0(1,IPHAS)
          XYP0   = XYZP0(2,IPHAS)
          XZP0   = XYZP0(3,IPHAS)
          DO IEL=1,NCEL
            PROPCE(IEL,IPPROC(IPRTOT(IPHAS)))= RTP(IEL,IPR(IPHAS))
     &           + RO0IPH*( GX*(XYZCEN(1,IEL)-XXP0)
     &                    + GY*(XYZCEN(2,IEL)-XYP0)
     &                    + GZ*(XYZCEN(3,IEL)-XZP0) )
     &           + P0IPH - PR0IPH
          ENDDO
        ENDIF
C
C
      ENDDO
C
C=======================================================================
C 7.  IMPRESSIONS
C=======================================================================
C
      DO IPHAS = 1, NPHAS
C
        IPRIPH = IPR(IPHAS)
        IUIPH  = IU(IPHAS)
        IVIPH  = IV(IPHAS)
        IWIPH  = IW(IPHAS)
C
        IFLMAS = IPPROF(IFLUMA(IUIPH))
        IFLMAB = IPPROB(IFLUMA(IUIPH))
        IPCROM = IPPROC(IROM  (IPHAS))
        IPBROM = IPPROB(IROM  (IPHAS))
C
        IF (IWARNI(IUIPH).GE.1) THEN
C
          WRITE(NFECRA,2000)IPHAS
C
          RNORM = -1.D0
          DO IEL = 1, NCEL
            RNORM  = MAX(RNORM,ABS(RTP(IEL,IPRIPH)))
          ENDDO
          IF (IRANGP.GE.0) CALL PARMAX (RNORM)
C                          ===========
          WRITE(NFECRA,2100)RNORM
C
          RNORM = -1.D0
          DO IEL = 1, NCEL
            VITNOR =
     &       SQRT(RTP(IEL,IUIPH)**2+RTP(IEL,IVIPH)**2+RTP(IEL,IWIPH)**2)
            IF(VITNOR.GE.RNORM) THEN
              RNORM = VITNOR
              IMAX  = IEL
            ENDIF
          ENDDO
C
          XYZMAX(1) = XYZCEN(1,IMAX)
          XYZMAX(2) = XYZCEN(2,IMAX)
          XYZMAX(3) = XYZCEN(3,IMAX)
C
          IF (IRANGP.GE.0) THEN
            NBRVAL = 3
            CALL PARMXL (NBRVAL, RNORM, XYZMAX)
C           ===========
          ENDIF
C                          ===========
C
          WRITE(NFECRA,2200) RNORM,XYZMAX(1),XYZMAX(2),XYZMAX(3)
C
C
C Pour la periodicite et le parallelisme, rom est echange dans phyvar
C
C
          RNORMA = -GRAND
          RNORMI =  GRAND
          DO IFAC = 1, NFAC
            IEL1 = IFACEL(1,IFAC)
            IEL2 = IFACEL(2,IFAC)
            SURF = RA(ISRFAN-1+IFAC)
            RHOM = (PROPCE(IEL1,IPCROM)+PROPCE(IEL2,IPCROM))*0.5D0
            RNORM = PROPFA(IFAC,IFLMAS)/(SURF*RHOM)
            RNORMA = MAX(RNORMA,RNORM)
            RNORMI = MIN(RNORMI,RNORM)
          ENDDO
          IF (IRANGP.GE.0) THEN
            CALL PARMAX (RNORMA)
C           ===========
            CALL PARMIN (RNORMI)
C           ===========
          ENDIF
          WRITE(NFECRA,2300)RNORMA, RNORMI
C
          RNORMA = -GRAND
          RNORMI =  GRAND
          DO IFAC = 1, NFABOR
            RNORM = PROPFB(IFAC,IFLMAB)/
     &             (RA(ISRFBN-1+IFAC)*PROPFB(IFAC,IPBROM))
            RNORMA = MAX(RNORMA,RNORM)
            RNORMI = MIN(RNORMI,RNORM)
          ENDDO
          IF (IRANGP.GE.0) THEN
            CALL PARMAX (RNORMA)
C           ===========
            CALL PARMIN (RNORMI)
C           ===========
          ENDIF
          WRITE(NFECRA,2400)RNORMA, RNORMI
C
          RNORM = 0.D0
          DO IFAC = 1, NFABOR
            RNORM = RNORM + PROPFB(IFAC,IFLMAB)
          ENDDO
C
          IF (IRANGP.GE.0) CALL PARSOM (RNORM)
C                          ===========
C
          WRITE(NFECRA,2500)RNORM
C
          WRITE(NFECRA,2001)
C
          IF(NTERUP.GT.1) THEN
            IF(ICVRGE.EQ.0) THEN
              WRITE(NFECRA,2600) ITERNS
              WRITE(NFECRA,2601) XNRMU(IPHAS),
     &                           XNRMU0(IPHAS), EPSUP(IPHAS)
              WRITE(NFECRA,2001)
              IF(ITERNS.EQ.NTERUP) THEN
                WRITE(NFECRA,2603)
                WRITE(NFECRA,2001)
              ENDIF
            ELSE
              WRITE(NFECRA,2602) ITERNS
              WRITE(NFECRA,2601) XNRMU(IPHAS),
     &                           XNRMU0(IPHAS), EPSUP(IPHAS)
              WRITE(NFECRA,2001)
            ENDIF
          ENDIF
C
        ENDIF
C
      ENDDO
C
C--------
C FORMATS
C--------
#if defined(_CS_LANG_FR)
C
 1000 FORMAT(/,
     &'   ** RESOLUTION POUR LA VITESSE                             ',/,
     &'      --------------------------                             ',/)
 1200 FORMAT(/,
     &'   ** RESOLUTION POUR LA PRESSION CONTINUITE                 ',/,
     &'      --------------------------------------                 ',/)
 2000 FORMAT(/,' APRES PRESSION CONTINUITE',/,
     &'  -- Phase : ',I10                                            ,/,
     &'-------------------------------------------------------------'  )
 2100 FORMAT(
     &' Pression max.',E12.4   ,' (max. de la valeur absolue)       ',/)
 2200 FORMAT(
     &' Vitesse  max.',E12.4   ,' en',3E11.3                         ,/)
 2300 FORMAT(
     &' Vitesse  en face interne max.',E12.4   ,' ; min.',E12.4        )
 2400 FORMAT(
     &' Vitesse  en face de bord max.',E12.4   ,' ; min.',E12.4        )
 2500 FORMAT(
     &' Bilan de masse   au bord   ',E14.6                             )
 2600 FORMAT(
     &' Informations Point fixe a l''iteration :',I10                ,/)
 2601 FORMAT('norme = ',E12.4,' norme 0 = ',E12.4,' toler  = ',E12.4 ,/)
 2602 FORMAT(
     &' Convergence du point fixe a l''iteration ',I10               ,/)
 2603 FORMAT(
     &' Non convergence du couplage vitesse pression par point fixe  ' )
 2001 FORMAT(
     &'-------------------------------------------------------------',/)
C
#else
C
 1000 FORMAT(/,
     &'   ** SOLVING VELOCITY'                                       ,/,
     &'      ----------------'                                       ,/)
 1200 FORMAT(/,
     &'   ** SOLVING CONTINUITY PRESSURE'                            ,/,
     &'      ---------------------------'                            ,/)
 2000 FORMAT(/,' AFTER CONTINUITY PRESSURE',/,
     &'  -- Phase : ',I10                                            ,/,
     &'-------------------------------------------------------------'  )
 2100 FORMAT(
     &' Max. pressure',E12.4   ,' (max. absolute value)'             ,/)
 2200 FORMAT(
     &' Max. velocity',E12.4   ,' en',3E11.3                         ,/)
 2300 FORMAT(
     &' Max. velocity at interior face',E12.4   ,' ; min.',E12.4       )
 2400 FORMAT(
     &' Max. velocity at boundary face',E12.4   ,' ; min.',E12.4       )
 2500 FORMAT(
     &' Mass balance  at boundary  ',E14.6                             )
 2600 FORMAT(
     &' Fixed point informations at iteration:',I10                  ,/)
 2601 FORMAT('norm = ',E12.4,' norm 0 = ',E12.4,' toler  = ',E12.4   ,/)
 2602 FORMAT(
     &' Fixed point convergence at iteration ',I10                   ,/)
 2603 FORMAT(
     &' Non convergence of fixed point for velocity pressure coupling' )
 2001 FORMAT(
     &'-------------------------------------------------------------',/)
C
#endif
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
