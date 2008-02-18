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
                        SUBROUTINE D3PPHY
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IBROM  , IZFPPP ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   DIRMIN , DIRMAX , FDEB   , FFIN   , HREC   ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA      )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC ROUTINE PHYSIQUE PARTICULIERE : FLAMME DE DIFFUSION
CFONC  Calcul de RHO mutualise pour chimie 3 points
CFONC   adiabatique ou permeatique (transport de H)
CFONC
CFONC
c@fonce
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
CARGU ! NPHMX        ! E  !  -> ! NPHSMX                               !
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
CARGU ! IBROM        ! TE !  -> ! INDICATEUR DE REMPLISSAGE DE ROMB    !
CARGU !   (NPHMX   ) !    !     !                                      !
CARGU ! IZFPPP       ! TE ! <-  ! NUMERO DE ZONE DE LA FACE DE BORD    !
CARGU ! (NFABOR)     !    !     !  POUR LE MODULE PHYS. PART.          !
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
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! DIRMIN       ! TR !  -  ! PDF : DIRAC EN FMIN                  !
CARGU ! DIRMAX       ! TR !  -  ! PDF : DIRAC EN FMAX                  !
CARGU ! FDEB         ! TR !  -  ! PDF : ABSCISSE DEBUT RECTANGLE       !
CARGU ! FFIN         ! TR !  -  ! PDF : ABSCISSE FIN RECTANGLE         !
CARGU ! HREC         ! TR !  -  ! PDF : HAUTEUR RECTANGLE              !
CARGU ! W1-W3(NCELET)! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
CARGU ! RA(*)        ! TR !  -  ! MACRO TABLEAU REEL                   !
CARGU ! RPP          ! TR !  -  ! MACRO TABLEAU REEL PP                !
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
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          IINPDF
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR), IBROM(NPHMX)
      INTEGER          IZFPPP(NFABOR)
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
      DOUBLE PRECISION DIRMIN(NCELET),DIRMAX(NCELET)
      DOUBLE PRECISION FDEB(NCELET),FFIN(NCELET)
      DOUBLE PRECISION HREC(NCELET)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA, IFINIA
      INTEGER          IF, IH, IEL, ICG
      INTEGER          IFAC, MODE, IZONE
      INTEGER          IPHAS, IPBROM, IPCROM, IPBYCG, IPCYCG
      DOUBLE PRECISION COEFG(NGAZGM), FSIR, HHLOC, TSTOEA, TIN
      DOUBLE PRECISION TEMSMM
C
      INTEGER       IPASS
      DATA          IPASS /0/
      SAVE          IPASS
C
C***********************************************************************
C=======================================================================
C 0. ON COMPTE LES PASSAGES
C=======================================================================
C
      IPASS = IPASS + 1
C
C=======================================================================
C 1. INITIALISATIONS A CONSERVER
C=======================================================================
C
C --- Initialisation memoire
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C     On reserve la memoire pour le tabeau INDPDF : passage ou non
C     par les PDF (on le garde pendant tout le sous-programme
C
      IINPDF = IDEBIA
      IFINIA = IINPDF + NCELET
      CALL IASIZE('D3PPHY',IFINIA)
C     ===========
C
C=======================================================================
C 2. DETERMINATION DES GRANDEURS THERMOCHIMIQUES
C=======================================================================
C
      IF ( IPASS.LE.2 ) THEN
C
C Rq : Il faut avoir vu usd3pc.F pour calculer correctement HH et FF
C
C ---> Calcul de TSTOEA
C
C ---- Initialisation
        DO ICG = 1, NGAZGM
          COEFG(ICG) = ZERO
        ENDDO
C
        HSTOEA = FS(1)*HINFUE + (1.D0-FS(1))*HINOXY
        COEFG(1) = ZERO
        COEFG(2) = ZERO
        COEFG(3) = 1.D0
        MODE = 1
        CALL COTHHT
C       ===========
     &  ( MODE   , NGAZG , NGAZGM  , COEFG  ,
     &    NPO    , NPOT   , TH     , EHGAZG ,
     &    HSTOEA , TSTOEA )
C
C
C ---> Construction d'une table Temperature en fonction de la richesse
C        et de l'enthalpie stoechiometrique
C        de dimension 9X9
C
        FSIR = FS(1)
C
C ---- Calcul du tableau FF(IF)
C
        DO IF = 1, (NMAXF/2+1)
          FF(IF) = FSIR * DBLE(2*IF-2)/DBLE(NMAXF-1)
        ENDDO
        DO IF = (NMAXF/2+2), NMAXF
          FF(IF) = FSIR + DBLE(2*IF-NMAXF-1)
     &                  / DBLE(NMAXF-1)*(1.D0-FSIR)
        ENDDO
C
C ----  Remplissage du tableau HH(IH)
C
        COEFG(1) = ZERO
        COEFG(2) = ZERO
        COEFG(3) = 1.D0
        TIN = MIN(TINFUE,TINOXY)
        MODE    = -1
        CALL COTHHT
C       ===========
     &  ( MODE      , NGAZG , NGAZGM  , COEFG  ,
     &    NPO       , NPOT   , TH     , EHGAZG ,
     &    HH(NMAXH) , TIN    )
        HH(1) = HSTOEA
        DO IH = 2, (NMAXH-1)
          HH(IH) = HH(1) + (HH(NMAXH)-HH(1))*
     &                    DBLE(IH-1)/DBLE(NMAXH-1)
        ENDDO
C
C ---- Remplissage du tableau TFH(IF,IH)
C
        DO IH = 1, NMAXH
          DO IF = 1, (NMAXF/2+1)
C ----- Melange pauvre
            COEFG(1) = ZERO
            COEFG(2) = (FSIR-FF(IF))/FSIR
            COEFG(3) = FF(IF)/FSIR
            HHLOC = HINOXY + DBLE(2*IF-2)/DBLE(NMAXF-1)
     &                     * (HH(IH)-HINOXY)
            MODE = 1
            CALL COTHHT
C           ===========
     &      ( MODE   , NGAZG , NGAZGM  , COEFG  ,
     &        NPO    , NPOT   , TH     , EHGAZG ,
     &        HHLOC , TFH(IF,IH) )
          ENDDO
          DO IF = (NMAXF/2+2), NMAXF
C ----- Melange riche
            COEFG(1) = (FF(IF)-FSIR)/(1.D0-FSIR)
            COEFG(2) = 0.D0
            COEFG(3) = (1.D0-FF(IF))/(1.D0-FSIR)
            HHLOC = ( DBLE(2*IF)*(HINFUE-HH(IH))
     &                + DBLE(2*NMAXF)*HH(IH)
     &                - HINFUE*DBLE(NMAXF+1) )
     &            / DBLE(NMAXF-1)
            MODE = 1
            CALL COTHHT
C           ===========
     &      ( MODE   , NGAZG , NGAZGM  , COEFG  ,
     &        NPO    , NPOT   , TH     , EHGAZG ,
     &        HHLOC , TFH(IF,IH) )
          ENDDO
C
        ENDDO
C
      ENDIF
C
C
C=======================================================================
C 3. CALCUL DE PARAMETRES DE LA FONCTION DENSITE DE PROBABILITE
C    POUR LA FRACTION DE MELANGE
C=======================================================================
C
C --- Definition des bornes min et max de la pdf
C     dans les 2 tableaux de travail W1 et W2
C
      DO IEL = 1, NCEL
        W1(IEL) = 0.D0
        W2(IEL) = 1.D0
      ENDDO
C
      CALL PPPDFR
C     ===========
     & ( NCELET,NCEL, IA(IINPDF),
     &   RTP(1,ISCA(IFM)), RTP(1,ISCA(IFP2M)),
     &   W1, W2,
     &   DIRMIN, DIRMAX, FDEB, FFIN, HREC )
C
C
C
C=======================================================================
C 4. INTEGRATION DE LA FONCTION DENSITE DE PROBABILITE
C    POUR DETERMINER LA TEMPERATURE
C                    LES FRACTIONS MASSIQUES
C                    LA MASSE VOLUMIQUE
C                    qsp RAYONNEMENT
C    Ces variables d'etat sont dans PROPCE
C=======================================================================
C
      IF ( IPPMOD(ICOD3P).EQ.1 ) THEN
C
        CALL D3PINT
C       ===========
     & ( NCELET,NCEL, IA(IINPDF),
     &   DIRMIN,DIRMAX,FDEB,FFIN,HREC,
     &   RTP(1,ISCA(IFM)),RTP(1,ISCA(IHM)),RTP(1,IPR(1)),
     &   PROPCE,
     &   W1 )
C
      ELSE
C
        CALL D3PINT
C       ===========
     & ( NCELET,NCEL, IA(IINPDF),
     &   DIRMIN,DIRMAX,FDEB,FFIN,HREC,
     &   RTP(1,ISCA(IFM)),W2,RTP(1,IPR(1)),
     &   PROPCE,
     &   W1 )
C
      ENDIF
C
C
C=======================================================================
C 4. CALCUL DE RHO ET DES FRACTIONS MASSIQUES DES ESPECES GLOBALES
C    SUR LES BORDS
C=======================================================================
C
C --> Masse volumique au bord
C
      IPHAS = 1
      IBROM(IPHAS) = 1
      IPBROM = IPPROB(IROM(IPHAS))
      IPCROM = IPPROC(IROM(IPHAS))
C
C ---- Masse volumique au bord pour toutes les facettes
C      Les facettes d'entree seront recalculees apres
C      a partir des CL (si IPASS > 1). .
C
C
C ---- Au premier passage sans suite ou si on n'a pas relu la
C      masse volumique dans le fichier suite, on n'a pas recalcule la
C      masse volumique dans d3pint, pas la peine de la reprojeter aux
C      faces.
C
      IF ( IPASS.GT.1.OR.(ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1)) THEN
C
        DO IFAC = 1, NFABOR
          IEL = IFABOR(IFAC)
          PROPFB(IFAC,IPBROM) = PROPCE(IEL,IPCROM)
        ENDDO
C
      ENDIF
C
C ---- Masse volumique au bord pour les facettes d'entree UNIQUEMENT
C      Le test sur IZONE sert pour les reprises de calcul
C      On suppose implicitement que les elements ci-dessus ont ete relus
C      dans le fichier suite (i.e. pas de suite en combustion d'un calcul
C      a froid) -> sera pris en compte eventuellement dans les versions
C      suivantes
C
      IF(IPASS.GT.1 .OR. ISUITE.EQ.1 ) THEN
        DO IFAC = 1, NFABOR
          IZONE = IZFPPP(IFAC)
          IF(IZONE.GT.0) THEN
            IF ( IENTFU(IZONE).EQ.1 .OR. IENTOX(IZONE).EQ.1 ) THEN
              TEMSMM = TINFUE/WMOLG(1)
              IF ( IENTOX(IZONE).EQ.1 ) TEMSMM = TINOXY/WMOLG(2)
              PROPFB(IFAC,IPBROM) = P0(IPHAS)/(RR*TEMSMM)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C --> Fractions massiques des especes globales au bord
C     Uniquement si rayonnement
C
      IF ( IRAYPP.GT.0 ) THEN
        DO ICG = 1, NGAZG
          IPBYCG = IPPROB(IYM(ICG))
          IPCYCG = IPPROC(IYM(ICG))
          DO IFAC = 1, NFABOR
            IEL = IFABOR(IFAC)
            PROPFB(IFAC,IPBYCG) = PROPCE(IEL,IPCYCG)
          ENDDO
        ENDDO
      ENDIF
C
C
C=======================================================================
C FORMATS
C----
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
