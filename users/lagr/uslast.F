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
                        SUBROUTINE USLAST
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , ETTPA  , TEPA   , TAUP   , TLAG   , TEMPCT ,
     &   STATIS , STATIV ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC        SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC        -----------------------------------
CFONC
CFONC     SOUS-PROGRAMME UTILISATEUR (INTERVENTION NON OBLIGATOIRE)
CFONC
CFONC     MODIFICATIONS UTILSATEUR SUR LES VARIABLES EN FIN D'ITERATION
CFONC     LAGRANGIENNES ET CALCUL DES STATISTIQUES UTILISATEUR
CFONC     SUPPLEMENTAIRES SUR LES PARTICULES
CFONC
CFONC    POUR LES STATISTIQUES UTILISATEUR SUPPLEMENTAIRES,
CFONC    ON RAPPELLE QUE :
CFONC
CFONC    ISTTIO = 0 : calcul instationnaire pour le lagrangien
CFONC           = 1 : calcul stationnaire   pour le lagrangien
CFONC
CFONC    ISTALA : calcul statistiques       si  >= 1 sinon pas de stat
CFONC
CFONC    ISUIST : suite calcul statistiques si  >= 1 sinon pas de stat
CFONC
CFONC    IDSTNT : Numero du pas de temps pour debut statistque
CFONC
CFONC    NSTIST : iteration Lagrangienne du debut calcul stationnaire
CFONC
CFONC    NPST   : nombre d'iterations de calcul de stat stationnaires
CFONC
CFONC    NPSTT  : nombre d'iterations total des stats depuis le debut
CFONC             du calcul, partie instationnaire comprise
CFONC
CFONC    TSTAT  : Temps physique d'enregistrement des stats volumiques
CFONC             stationnaires
CFONC             (en instationnaire TSTAT=DTP le pas de temps Lagrangien)
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
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC           !
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR           !
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NTERSL       ! E  !  -> ! NBR TERMES SOURCES DE COUPLAGE RETOUR!
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VAR STATISTIQUES LAGRANGIEN!
CARGU ! NVISBR       ! E  !  -> ! NOMBRE DE STATISTIQUES AUX FRONTIERES!
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
CARGU !  (NFABOR+1)  !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (LNDFBR  ) !    !     !  (OPTIONNEL)                         !
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
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
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS                !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME(NCELET! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
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
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA(NBPMAX, ! TR !  -> ! CARACTERISTIQUES DES PARTICULES      !
CARGU !       NVEP)  !    !     !  AUX PARTICULES (poids, ...)         !
CARGU ! TAUP(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE DYNAMIQUE      !
CARGU ! TLAG(NBPMAX) ! TR !  -> ! TEMPS CARACTERISTIQUE FLUIDE         !
CARGU ! TEMPCT       ! TR !  -> ! TEMPS CARACTERISTIQUE THERMIQUE      !
CARGU !  (NBPMAX,2)  !    !     !                                      !
CARGU ! STATIS       ! TR ! <-> ! CUMUL POUR LES MOYENNES DES          !
CARGU !(NCELET,NVLSTA!    !     !   STATISTIQUES VOLUMIQUES            !
CARGU ! STATIV       ! TR ! <-> ! CUMUL POUR LES VARIANCES DES         !
CARGU !(NCELET,      !    !     !    STATISTIQUES VOLUMIQUES           !
CARGU !   NVLSTA-1)  !    !     !                                      !
CARGU ! W1..W3(NCELET! TR !  -  ! TABLEAUX DE TRAVAIL                  !
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
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "radiat.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "cstphy.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "cpincl.h"
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
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1) , NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1) , NODFBR(LNDFBR)
      INTEGER          ITEPA(NBPMAX,NIVEP)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD) , VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET) , RTP(NCELET,*) , RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*) , COEFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION TAUP(NBPMAX) , TLAG(NBPMAX,3) , TEMPCT(NBPMAX,2)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA)
      DOUBLE PRECISION STATIV(NCELET,NVLSTA-1)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE) , RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IFINIA, IFINRA
      INTEGER          NPT ,  IEL , IPHAS
C
      INTEGER          IVF , IVFF , ITABVR , IFLU , ICLA
C
C VARIABLES LOCALES UTILISATEUR
C
      INTEGER          NXLIST
      PARAMETER       (NXLIST=100)
C
      INTEGER          IPLAN
      INTEGER          II, IND, IL
      INTEGER          INOEUD, IRANG0, INDIC
      INTEGER          IST(6)
C
      DOUBLE PRECISION ZZ(4), ZZZ(8), XLIST(NXLIST,8), XYZPT(3)

      CHARACTER        NAME(8)*4
C
      DOUBLE PRECISION DEBM(4)
      SAVE             DEBM
C
C***********************************************************************
C
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C
      IF(ISTALA.EQ.1 .AND. IPLAS.GE.IDSTNT .AND. NVLSTS.GT.0) THEN
C
C     Si l'on passe ici, il faut que l'utilisateur complete
C       l'exemple ci-dessous et l'adapte...
C
        IF(1.EQ.1) THEN
          WRITE(NFECRA,9000)NVLSTS
          CALL CSEXIT (1)
        ENDIF
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS LE MODULE LAGRANGIEN             ',/,
     &'@    *********                                               ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR uslast DOIT ETRE COMPLETE',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Des variables statistiques supplementaires ont ete        ',/,
     &'@    demandees dans uslag1 (NVLSTS=',   I10,')               ',/,
     &'@  Le sous-programme uslast doit etre complete pour preciser ',/,
     &'@    le  calcul de leur cumul.                               ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      ELSE
C
C     On entre toujours dans ce sous programme en lagrangien,
C       si on ne souhaite rien y faire, on sort immediatement.
C
        RETURN
C
      ENDIF
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C
C=======================================================================
C 0.  GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IPHAS = ILPHAS
C
C=======================================================================
C 2 - CALCUL DES STATISTIQUES PARTICULAIRES UTILISATEURS
C=======================================================================
C
C   D'une facon generale, dans cette routine on realise les cumuls
C   de la quantite dont on souhaite faire les statistiques.
C   La moyenne et la variance sont calculees dans la routine
C   USLAEN.F. Ce calcul est le plus souvent obtenu par division
C   des cumuls soit par le temps du cumul stationnaire contenu dans
C   la variable TSTAT, soit par le nombre de particules en poids
C   statistiques. Cette division est appliquee pour chaque ecriture
C   dans le listing et pour les sorties post-processing.
C
C   Cet exemple est desactive et doit etre adapte au cas traite
C
      IF (1.EQ.0) THEN
C
       IF(ISTALA.EQ.1 .AND. IPLAS.GE.IDSTNT .AND. NVLSTS.GT.0) THEN
C
        DO NPT = 1,NBPART
C
          IF( ITEPA(NPT,JISOR).GT.0 ) THEN
C
            IEL = ITEPA(NPT,JISOR)
C
C -------------------------------------------------
C EXEMPLE 1 : Cumul pour la concentration massique
C -------------------------------------------------
C
            STATIS(IEL,ILVU(1)) = STATIS(IEL,ILVU(1))
     &        + TEPA(NPT,JRPOI) *ETTP(NPT,JMP)
C
            STATIV(IEL,ILVU(1)) = STATIV(IEL,ILVU(1))
     &        + TEPA(NPT,JRPOI) *ETTP(NPT,JMP) *ETTP(NPT,JMP)
C
          ENDIF
C
        ENDDO
C
       ENDIF
C
      ENDIF
C
C=======================================================================
C 3 - CALCUL UTILISATEUR DU DEBIT MASSIQUE DE PARTICULES SUR 4 PLANS
C=======================================================================
C
C   Cet exemple est desactive et doit etre adapte au cas traite
C
      IF (1.EQ.0) THEN
C
        ZZ(1) = 0.1D0
        ZZ(2) = 0.15D0
        ZZ(3) = 0.20D0
        ZZ(4) = 0.25D0
C
C   Si on est en instationnaire, ou si le debut des stat stationnaires
C   n'est pas encore atteint, toutes les statistiques sont remises a
C   zero a chaque pas de temps avant d'entrer dans ce sous-programme.
C
        IF(ISTTIO.EQ.0 .OR. NPSTT.LE.NSTIST) THEN
          DO IPLAN = 1,4
            DEBM(IPLAN) = 0.D0
          ENDDO
        ENDIF
C
        DO IPLAN = 1,4

          DO NPT = 1,NBPART
C
            IF(ITEPA(NPT,JISOR).GT.0) THEN
C
              IEL = ITEPA(NPT,JISOR)
C
              IF( ETTP(NPT,JXP).GT.ZZ(IPLAN) .AND.
     &            ETTPA(NPT,JXP).LE.ZZ(IPLAN)      ) THEN
                DEBM(IPLAN) = DEBM(IPLAN) +TEPA(NPT,JRPOI)*ETTP(NPT,JMP)
              ENDIF
C
            ENDIF
C
          ENDDO
        ENDDO
C
        DO IPLAN = 1,4
          WRITE(NFECRA,1001)IPLAN,DEBM(IPLAN)/TSTAT
        ENDDO
C
 1001   FORMAT(' Debit massique particulaire en Z(',I10,') : ',E14.5)
C
      ENDIF
C
C
C=======================================================================
C 4 - EXTRACTION DE STATISTIQUES VOLUMIQUES EN FIN DE CALCUL
C=======================================================================
C
C   Cet exemple est desactive et doit etre adapte au cas traite
C
      IF (1.EQ.0) THEN
C
        IF(NTCABS.EQ.NTMABS) THEN
C
          ZZZ(1) = 0.005D0
          ZZZ(2) = 0.025D0
          ZZZ(3) = 0.050D0
          ZZZ(4) = 0.075D0
          ZZZ(5) = 0.100D0
          ZZZ(6) = 0.150D0
          ZZZ(7) = 0.200D0
          ZZZ(8) = 0.250D0

          NAME(1) = 'XB01'
          NAME(2) = 'XB05'
          NAME(3) = 'XB10'
          NAME(4) = 'XB15'
          NAME(5) = 'XB20'
          NAME(6) = 'XB30'
          NAME(7) = 'XB40'
          NAME(8) = 'XB50'
C
          IST(1) = ILVX
          IST(2) = ILVZ
          IST(3) = ILFV
          IST(4) = ILPD
C
          NPTS = NXLIST
C
          IFINIA = IDEBIA
          ITABVR = IDEBRA
          IFINRA = ITABVR + NCELET
          CALL RASIZE('USLAST',IFINRA)
C         ===========
C
          DO IPLAN = 1,8
C
C     Pour le fichier ci-dessous :
C       l'utilisateur verifiera qu'il n'a pas laisse ouverte l'unite
C       IMPUSR(1), dans un autre sous-programme utilisateur
            OPEN(FILE=NAME(IPLAN),UNIT=IMPUSR(1),FORM='formatted')
C
            XYZPT(1) = ZZZ(IPLAN)
C
            DO IVF = 1,4
C
              IVFF = IST(IVF)
              ICLA = 0
              IFLU = 0
C
              CALL USLAEN
C             ===========
     & ( IFINIA , IFINRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NVLSTA ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVFF   , IVFF   , IVFF   , IFLU   , ILPD   , ICLA   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , STATIS , STATIV , RA(ITABVR) ,
     &   RDEVEL , RTUSER , RA     )
C
              IND = 0
              DO II = 1, NPTS
C
                XYZPT(2) = 0.D0
                XYZPT(3) = FLOAT(II-1)/FLOAT(NPTS-1)*150.D-3
C
                CALL FINDPT
C               ===========
     &          (NCELET, NCEL, XYZCEN,
     &           XYZPT(1), XYZPT(2), XYZPT(3), INOEUD, IRANG0)
C
                INDIC = ITUSER(INOEUD)
                ITUSER(INOEUD) = 1
                IF (INDIC.EQ.1) THEN
                  IND = IND +1
                  XLIST(IND,1) = XYZCEN(1,INOEUD)
                  XLIST(IND,2) = XYZCEN(3,INOEUD) * (1.D3 / 5.D0)
                  XLIST(IND,IVF+2) = RA(ITABVR+INOEUD-1)
                ENDIF
              ENDDO
            ENDDO
C
            DO IL = 1, IND
              WRITE (IMPUSR(1),'(8E13.5)') (XLIST(IL,II), II=1,6)
            ENDDO
C
            CLOSE(IMPUSR(1))
C
          ENDDO
C
        ENDIF
C
      ENDIF
C
C
C
C=======================================================================
C
C====
C FIN
C====
C
      RETURN
C
      END
c@z
