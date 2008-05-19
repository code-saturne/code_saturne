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
                        SUBROUTINE ELPHYV
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
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC    REMPLISSAGE DES VARIABLES PHYSIQUES : Version Electrique
CFONC
CFONC      ----> Effet Joule
CFONC      ----> Arc Electrique
CFONC      ----> Conduction Ionique
CFONC
CFONC       1) Masse Volumique
CFONC       2) Viscosite moleculaire
CFONC       3) Cp
CFONC       4) Lambda/Cp moleculaire
CFONC       4) Diffusivite moleculaire
CFONC
CFONC
CFONC
CFONC  ATTENTION :
CFONC  =========
CFONC
CFONC
CFONC  Il est INTERDIT de modifier la viscosite turbulente VISCT ici
CFONC         ========
CFONC   (une routine specifique est dediee a cela : usvist)
CFONC
CFONC
CFONC   Il FAUT AVOIR PRECISE ICP(IPHAS) = 1
CFONC      ==================
CFONC     dans usini1 si on souhaite imposer une chaleur specifique
CFONC     CP variable pour la phase IPHAS (sinon: ecrasement memoire).
CFONC
CFONC
CFONC   Il FAUT AVOIR PRECISE IVISLS(Numero de scalaire) = 1
CFONC      ==================
CFONC      dans usini1 si on souhaite une diffusivite VISCLS variable
CFONC      pour le scalaire considere (sinon: ecrasement memoire).
CFONC
CFONC
CFONC
CFONC
CFONC  Remarques :
CFONC  ---------
CFONC
CFONC  Cette routine est appelee au debut de chaque pas de temps
CFONC
CFONC     Ainsi, AU PREMIER PAS DE TEMPS (calcul non suite), les seules
CFONC     grandeurs initialisees avant appel sont celles donnees
CFONC       - dans usini1 :
CFONC          . la masse volumique (initialisee a RO0(IPHAS))
CFONC          . la viscosite       (initialisee a VISCL0(IPHAS))
CFONC       - dans usiniv :
CFONC          . les variables de calcul  (initialisees a 0 par defaut
CFONC              ou a la valeur donnee dans usiniv)
CFONC
CFONC  On peut donner ici les lois de variation aux cellules
CFONC      - de la masse volumique                      ROM    kg/m3
CFONC          (et eventuellememt aux faces de bord     ROMB   kg/m3)
CFONC      - de la viscosite moleculaire                VISCL  kg/(m s)
CFONC      - de la chaleur specifique associee          CP     J/(kg degres)
CFONC      - des "diffusivites" associees aux scalaires VISCLS kg/(m s)
CFONC
CFONC
CFONC  On dispose des types de faces de bord au pas de temps
CFONC    precedent (sauf au premier pas de temps, ou les tableaux
CFONC    ITYPFB et ITRIFB n'ont pas ete renseignes)
CFONC
CFONC
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
CARGU ! IZFPPP       ! TE !  -> ! NUMERO DE ZONE DE LA FACE DE BORD    !
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
CARGU ! W1...8(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "elincl.h"
      INCLUDE "radiat.h"
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
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET),W4(NCELET)
      DOUBLE PRECISION W5(NCELET),W6(NCELET),W7(NCELET),W8(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA, IFINIA
      INTEGER          IEL   , IPHAS
      INTEGER          IPCROM, IPCVIS, IPCCP , IPCRAY
      INTEGER          IPCVSL, ITH   , ISCAL , II
      INTEGER          IIII  , IPCSIG, IT
      INTEGER          IESP  , IESP1 , IESP2 , MODE , ISRROM
      INTEGER          MAXELT, ILS
C
      DOUBLE PRECISION TP    , DELT  , SOMPHI, VAL
      DOUBLE PRECISION ALPRO , ALPVIS, ALPCP , ALPSIG, ALPLAB , ALPKAB
      DOUBLE PRECISION RHONP1
      DOUBLE PRECISION YM    (NGAZGM),YVOL  (NGAZGM)
      DOUBLE PRECISION COEF(NGAZGM,NGAZGM)
      DOUBLE PRECISION ROESP (NGAZGM),VISESP(NGAZGM),CPESP(NGAZGM)
      DOUBLE PRECISION SIGESP(NGAZGM),XLABES(NGAZGM),XKABES(NGAZGM)
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C***********************************************************************
C=======================================================================
C 0 - INITIALISATIONS A CONSERVER
C=======================================================================
C
C --- Initialisation memoire
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IPASS = IPASS + 1
C
      IPHAS = 1
C
C     Sous relaxation de la masse volumique (pas au premier pas de temps)
      IF(NTCABS.GT.1.AND.SRROM.GT.0.D0) THEN
        ISRROM = 1
      ELSE
        ISRROM = 0
      ENDIF
C
C=======================================================================
C 1 - EFFET JOULE
C=======================================================================
C
C  -- Les lois doivent etre imposees par l'utilisateur
C       donc on ne fait rien.
C
C      IF ( IPPMOD(IELJOU).GE.1 ) THEN
C
C
C  -- Attention, dans les modules electriques, la chaleur massique, la
C       conductivite thermique et la conductivite electriques sont
C       toujours dans le tableau PROPCE
C       qu'elles soient physiquement variables ou non.
C
C       On n'utilisera donc PAS les variables
C          =====================
C                                CP0(IPHAS), VISLS0(ISCALT(IPHAS))
C                                VISLS0(IPOTR) et VISLS0(IPOTI)
C
C       Informatiquement, ceci se traduit par le fait que
C                                ICP(IPHAS)>0, IVISLS(ISCALT(IPHAS))>0,
C                                IVISLS(IPOTR)>0 et IVISLS(IPOTI)>0
C
C       Les verifications ont ete faites dans elveri
C
C  -- Si la conductivite electrique est toujours la meme pour
C       le potentiel reel et le potentiel imaginaire, on pourrait
C       n'en avoir qu'une seule (modif dans varpos pour definir
C       IVISLS(IPOTI) = IVISLS(IPOTR)) et economiser NCEL reels .
C
C      IPCROM = IPPROC(IROM(IPHAS))
C      IPCVIS = IPPROC(IVISCL(IPHAS))
C      IPCCP  = IPPROC(ICP(IPHAS))
C      IPCVSL = IPPROC(IVISLS(ISCALT(IPHAS)))
C      IPCSIR = IPPROC(IVISLS(IPOTR))
C      IPCSII = IPPROC(IVISLS(IPOTI))
C
C      PROPCE(IEL,IPPROC(ITEMP)) =
C      PROPCE(IEL,IPCROM) =
C      PROPCE(IEL,IPCVIS) =
C      PROPCE(IEL,IPCCP) =
C      PROPCE(IEL,IPCVSL) =
C      PROPCE(IEL,IPCSIR) =
C      PROPCE(IEL,IPCSII) =
C
C      ENDIF
C
C=======================================================================
C 2 - ARC ELECTRIQUE
C=======================================================================
C
      IF ( IPPMOD(IELARC).GE.1 ) THEN
C
C      Un message une fois au moins pour dire
C                                    qu'on prend les valeurs sur fichier
        IF(IPASS.EQ.1) THEN
          WRITE(NFECRA,1000)
        ENDIF
C
C      Calcul de la temperature a partir de l'enthalpie
C
        MODE = 1
C
        IF ( NGAZG .EQ. 1 ) THEN
          YM(1) = 1.D0
          MODE = 1
          DO IEL = 1, NCEL
            CALL ELTHHT(MODE,NGAZG,YM,RTP(IEL,ISCA(IHM)),
     &                                PROPCE(IEL,IPPROC(ITEMP)))
          ENDDO
        ELSE
          DO IEL = 1, NCEL
            YM(NGAZG) = 1.D0
            DO IESP = 1, NGAZG-1
              YM(IESP) = RTP(IEL,ISCA(IYCOEL(IESP)))
              YM(NGAZG) = YM(NGAZG) - YM(IESP)
            ENDDO
            CALL ELTHHT(MODE,NGAZG,YM,RTP(IEL,ISCA(IHM)),
     &                                PROPCE(IEL,IPPROC(ITEMP)))
          ENDDO
        ENDIF
C
C      Pointeurs pour les differentes variables
C
        IPCROM = IPPROC(IROM(IPHAS))
        IPCVIS = IPPROC(IVISCL(IPHAS))
        IF(ICP(IPHAS).GT.0) THEN
          IPCCP  = IPPROC(ICP(IPHAS))
        ENDIF
        IF(IVISLS(ISCALT(IPHAS)).GT.0) THEN
          IPCVSL = IPPROC(IVISLS(ISCALT(IPHAS)))
        ENDIF
        IF ( IVISLS(IPOTR).GT.0 ) THEN
          IPCSIG = IPPROC(IVISLS(IPOTR))
        ENDIF
        IF ( IXKABE .GT. 0 ) THEN
           IPCRAY = IPPROC(IDRAD)
        ENDIF
C
C       Interpolation des donnees sur le fichier de donnees
C         en fonction de la temperature
C
        DO IEL = 1, NCEL
C
C        Valeur de la temperature
C
          TP = PROPCE(IEL,IPPROC(ITEMP))
C
C        On determine  le IT ou il faut interpoler
C
          IT = 0
          IF ( TP .LE. TH(1) ) THEN
            IT = 1
          ELSE IF ( TP .GE. TH(NPO) ) THEN
            IT = NPO
          ELSE
            DO IIII = 1, NPO-1
              IF ( TP .GT. TH(IIII) .AND. TP .LE. TH(IIII+1) ) THEN
                IT = IIII
              ENDIF
            ENDDO
          ENDIF
          IF ( IT .EQ. 0 ) THEN
            WRITE(NFECRA,9900) TP
            CALL CSEXIT(1)
          ENDIF
C
C        Fraction massique
C
          YM(NGAZG) = 1.D0
          DO IESP = 1, NGAZG-1
            YM(IESP)  = RTP(IEL,ISCA(IYCOEL(IESP)))
            YM(NGAZG) = YM(NGAZG) - YM(IESP)
          ENDDO
C
C     Masse volumique, Viscosite, CP, Sigm et Lambda de chaque constituant
C
          IF ( TP . LE. TH(1) ) THEN
C
C         Extrapolation : Valeur constante = 1ere valeur de la table
C
            DO IESP = 1, NGAZG
              ROESP (IESP) = RHOEL (IESP,1)
              VISESP(IESP) = VISEL (IESP,1)
              CPESP (IESP) = CPEL  (IESP,1)
              SIGESP(IESP) = SIGEL (IESP,1)
              XLABES(IESP) = XLABEL(IESP,1)
              IF ( IXKABE .GT. 0 ) THEN
                XKABES(IESP) = XKABEL(IESP,1)
              ENDIF
            ENDDO
C
          ELSE IF ( TP .GE. TH(NPO) ) THEN
C
C         Extrapolation : valeur constante = derniere valeur de la table
C
            DO IESP = 1, NGAZG
              ROESP (IESP) = RHOEL (IESP,NPO)
              VISESP(IESP) = VISEL (IESP,NPO)
              CPESP (IESP) = CPEL  (IESP,NPO)
              SIGESP(IESP) = SIGEL (IESP,NPO)
              XLABES(IESP) = XLABEL(IESP,NPO)
              IF ( IXKABE .GT. 0 ) THEN
                XKABES(IESP) = XKABEL(IESP,NPO)
              ENDIF
            ENDDO
C
          ELSE
C
C         Interpolation
C
            DELT = TH(IT+1) - TH(IT)
            DO IESP = 1, NGAZG
C
C          Masse volumique de chaque constituant
C
              ALPRO = (RHOEL(IESP,IT+1)-RHOEL(IESP,IT))/DELT
              ROESP(IESP)  = RHOEL(IESP,IT) + ALPRO*(TP-TH(IT))
C
C          Viscosite de chaque constituant
C
              ALPVIS = (VISEL(IESP,IT+1)-VISEL(IESP,IT))/DELT
              VISESP(IESP) = VISEL(IESP,IT) + ALPVIS*(TP-TH(IT))
C
C          CP de chaque constituant
C
              ALPCP = (CPEL(IESP,IT+1)-CPEL(IESP,IT))/DELT
              CPESP(IESP) = CPEL(IESP,IT) + ALPCP*(TP-TH(IT))
C
C          Conductivite electrique (Sigma) de chaque constituant
C
              ALPSIG = (SIGEL(IESP,IT+1)-SIGEL(IESP,IT))/DELT
              SIGESP(IESP) = SIGEL(IESP,IT) + ALPSIG*(TP-TH(IT))
C
C          Conductivite thermique (Lambda) de chaque constituant
C
              ALPLAB = (XLABEL(IESP,IT+1)-XLABEL(IESP,IT))/DELT
              XLABES(IESP) = XLABEL(IESP,IT) + ALPLAB*(TP-TH(IT))
C
C          Emission nette radiative ou Terme source radiatif
C          de chaque constituant
C
              IF ( IXKABE .GT. 0 ) THEN
                ALPKAB = (XKABEL(IESP,IT+1)-XKABEL(IESP,IT))/DELT
                XKABES(IESP) = XKABEL(IESP,IT) + ALPKAB*(TP-TH(IT))
              ENDIF
C
            ENDDO
C
          ENDIF
C
C       Masse volumique du melange (sous relaxee eventuellement)
C       ==========================
C
          RHONP1 = 0.D0
          DO IESP = 1, NGAZG
            RHONP1 = RHONP1+YM(IESP)/ROESP(IESP)
          ENDDO
          RHONP1 = 1.D0/RHONP1
          IF(ISRROM.EQ.1) THEN
            PROPCE(IEL,IPCROM) =
     &           SRROM*PROPCE(IEL,IPCROM)+(1.D0-SRROM)*RHONP1
          ELSE
            PROPCE(IEL,IPCROM) = RHONP1
          ENDIF
C
C        Fraction volumique de chaque constituant
C
          DO IESP = 1, NGAZG
            YVOL(IESP) = YM(IESP)*ROESP(IESP)/PROPCE(IEL,IPCROM)
            IF ( YVOL(IESP) .LE. 0.D0 ) YVOL(IESP) = EPZERO**2
          ENDDO
C
C       Viscosite moleculaire dynamique en kg/(m s)
C       ==========================================
C
          DO IESP1 = 1, NGAZG
            DO IESP2 = 1, NGAZG
              COEF(IESP1,IESP2) = ( 1.D0
     &               + SQRT(VISESP(IESP1)/VISESP(IESP2))
     &                *SQRT(SQRT(ROESP(IESP2)/ROESP(IESP1))))**2.
     &                / SQRT( 1.D0 + ROESP(IESP1)/ROESP(IESP2) )
     &            / SQRT(8.D0)
            ENDDO
          ENDDO
C
          PROPCE(IEL,IPCVIS) = 0.D0
          DO IESP1=1,NGAZG
C
            SOMPHI = 0.D0
            DO IESP2=1,NGAZG
              IF ( IESP1 .NE. IESP2 ) THEN
                SOMPHI = SOMPHI
     &                  +COEF(IESP1,IESP2)*YVOL(IESP2)/YVOL(IESP1)
              ENDIF
            ENDDO
C
            PROPCE(IEL,IPCVIS) = PROPCE(IEL,IPCVIS)
     &                          +VISESP(IESP1)/(1.D0+SOMPHI)
C
          ENDDO
C
C       Chaleur specifique J/(kg degres)
C       ================================
C
          IF(ICP(IPHAS).GT.0) THEN
C
            PROPCE(IEL,IPCCP) = 0.D0
            DO IESP = 1, NGAZG
              PROPCE(IEL,IPCCP) = PROPCE(IEL,IPCCP )
     &                            +YM(IESP)*CPESP(IESP)
            ENDDO
C
          ENDIF
C
C       Lambda/Cp en kg/(m s)
C       ---------------------
C
          IF(IVISLS(ISCALT(IPHAS)).GT.0) THEN
C
            DO IESP1=1,NGAZG
              DO IESP2=1,NGAZG
                COEF(IESP1,IESP2) = ( 1.D0
     &                 + SQRT(XLABES(IESP1)/XLABES(IESP2))
     &                  *SQRT(SQRT(ROESP(IESP2)/ROESP(IESP1))))**2.D0
     &                  / SQRT( 1.D0 + ROESP(IESP1)/ROESP(IESP2) )
     &              / SQRT(8.D0)
              ENDDO
            ENDDO
C
C        On calcule d'abord juste Lambda
C
            PROPCE(IEL,IPCVSL) = 0.D0
            DO IESP1=1,NGAZG
C
              SOMPHI = 0.D0
              DO IESP2=1,NGAZG
                IF ( IESP1 .NE. IESP2 ) THEN
                  SOMPHI = SOMPHI
     &                    +COEF(IESP1,IESP2)*YVOL(IESP2)/YVOL(IESP1)
                ENDIF
              ENDDO
C
              PROPCE(IEL,IPCVSL) = PROPCE(IEL,IPCVSL)
     &                            +XLABES(IESP1)/(1.D0+1.065*SOMPHI)
C
            ENDDO
C
C        On divise par CP pour avoir Lambda/CP
C          On suppose Cp renseigne au prealable.
C
            IF(IPCCP.LE.0) THEN
C
C --- Si CP est uniforme, on utilise CP0(IPHAS)
C
              PROPCE(IEL,IPCVSL) = PROPCE(IEL,IPCVSL)/CP0(IPHAS)
C
            ELSE
C
C --- Si CP est non uniforme, on utilise le CP calcul au dessus
              PROPCE(IEL,IPCVSL) = PROPCE(IEL,IPCVSL)/PROPCE(IEL,IPCCP)
C
            ENDIF
          ENDIF
C
C       Conductivite electrique en S/m
C       ==============================
C
          IF ( IVISLS(IPOTR).GT.0 ) THEN
            PROPCE(IEL,IPCSIG) = 0.D0
            VAL = 0.D0
            DO IESP=1,NGAZG
              VAL = VAL + YVOL(IESP)/SIGESP(IESP)
            ENDDO
C
            PROPCE(IEL,IPCSIG) = 1.D0/VAL
          ENDIF
C
C       Emission nette radiative en W/m3
C       ================================
C
          IF ( IXKABE .GT. 0 ) THEN
            PROPCE(IEL,IPCRAY) = 0.D0
            VAL = 0.D0
            DO IESP=1,NGAZG
              VAL = VAL + YVOL(IESP)*XKABES(IESP)
            ENDDO
C
            PROPCE(IEL,IPCRAY) = VAL
          ENDIF
C
        ENDDO
C
C       Diffusivite variable a l'exclusion de l'enthalpie et de IPOTR
C       -------------------------------------------------------------
C         Il n'y a pas d'autres scalaires, et la boucle ne fait donc rien
C
C
        DO II = 1, NSCAPP
C
C --- Numero du scalaire
          ISCAL = ISCAPP(II)
C
C --- Si il s'agit de l'enthalpie son cas a deja ete traite plus haut
          ITH = 0
          IF (ISCAL.EQ.ISCALT(IPHAS)) ITH = 1
C
C --- Si il s'agit de Potentiel (IPOTR), son cas a deja ete traite
          IF (ISCAL.EQ.IPOTR) ITH = 1
C
C --- Si la variable est une fluctuation, sa diffusivite est
C       la meme que celle du scalaire auquel elle est rattachee :
C       il n'y a donc rien a faire ici : on passe directement
C       a la variable suivante sans renseigner PROPCE(IEL,IPCVSL).
C
          IF ( ITH.EQ.0 .AND. ISCAVR(ISCAL).LE.0) THEN
C
C --- On ne traite ici que les variables non thermiques
C                                        et pas le potentiel (sigma)
C                                   et qui ne sont pas des fluctuations
C
            IF(IVISLS(ISCAL).GT.0) THEN
C
C --- Rang de Lambda du scalaire
C     dans PROPCE, prop. physiques au centre des elements       : IPCVSL
C
              IPCVSL = IPPROC(IVISLS(ISCAL))
C
              DO IEL = 1, NCEL
                PROPCE(IEL,IPCVSL) = 1.D0
              ENDDO
C
            ENDIF
C
          ENDIF
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C 3 - CONDUCTION IONIQUE
C=======================================================================
C
C POUR LE MOMENT CETTE OPTION N'EST PAS ACTIVEE
C
      IF ( IPPMOD(IELION).GE.1  ) THEN
C
C       Masse volumique
C       ---------------
C
        IPCROM = IPPROC(IROM(IPHAS))
        DO IEL = 1, NCEL
          PROPCE(IEL,IPCROM) = 1.D0
        ENDDO
C
C       VISCOSITE
C       =========
C
        IPCVIS = IPPROC(IVISCL(IPHAS))
        DO IEL = 1, NCEL
          PROPCE(IEL,IPCVIS) = 1.D-2
        ENDDO
C
C       CHALEUR SPECIFIQUE VARIABLE J/(kg degres)
C       =========================================
C
        IF(ICP(IPHAS).GT.0) THEN
C
          IPCCP  = IPPROC(ICP   (IPHAS))
C
          DO IEL = 1, NCEL
            PROPCE(IEL,IPCCP ) = 1000.D0
          ENDDO
C
        ENDIF
C
C       Lambda/CP  VARIABLE en kg/(m s)
C       ===============================
C
        IF (IVISLS(ISCALT(IPHAS)).GT.0) THEN
C
          IPCVSL = IPPROC(IVISLS(ISCALT(IPHAS)))
C
          IF(IPCCP.LE.0) THEN
C
C --- Si CP est uniforme, on utilise CP0(IPHAS)
C
            DO IEL = 1, NCEL
              PROPCE(IEL,IPCVSL) = 1.D0/CP0(IPHAS)
            ENDDO
C
          ELSE
C
C --- Si CP est non uniforme, on utilise PROPCE ci dessus
            DO IEL = 1, NCEL
              PROPCE(IEL,IPCVSL) = 1.D0 /PROPCE(IEL,IPCCP)
            ENDDO
C
          ENDIF
C
        ENDIF
C
C       DIFFUSIVITE VARIABLE A L'EXCLUSION DE L'ENTHALPIE
C       ==================================================
C
        DO II = 1, NSCAPP
C
C --- Numero du scalaire
          ISCAL = ISCAPP(II)
C
C --- Si il s'agit de l'enthqlpie son cas a deja ete traite plus haut
          ITH = 0
          IF (ISCAL.EQ.ISCALT(IPHAS)) ITH = 1
C
C --- Si la variable est une fluctuation, sa diffusivite est
C       la meme que celle du scalaire auquel elle est rattachee :
C       il n'y a donc rien a faire ici : on passe directement
C       a la variable suivante sans renseigner PROPCE(IEL,IPCVSL).
C
          IF ( ITH.EQ.0 .AND. ISCAVR(ISCAL).LE.0) THEN
C
C --- On ne traite ici que les variables non thermiques
C                                   et qui ne sont pas des fluctuations
C
            IF(IVISLS(ISCAL).GT.0) THEN
C
C --- Rang de Lambda du scalaire
C     dans PROPCE, prop. physiques au centre des elements       : IPCVSL
C
              IPCVSL = IPPROC(IVISLS(ISCAL))
C
C --- Lambda en kg/(m s) au centre des cellules
C
C
              DO IEL = 1, NCEL
                PROPCE(IEL,IPCVSL) = 1.D0
              ENDDO
C
            ENDIF
C
          ENDIF
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C 4 - ON PASSE LA MAIN A L'UTILISATEUR (joule en particulier)
C=======================================================================
C
      MAXELT = MAX(NCELET, NFAC, NFABOR)
      ILS    = IDEBIA
      IFINIA = ILS + MAXELT
      CALL IASIZE('ELPHYV',IFINIA)
C
      CALL USELPH
C     ===========
     & ( IFINIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IBROM  , IZFPPP ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C
C
C
C La masse volumique au bord est traitee dans phyvar (recopie de la valeur
C     de la cellule de bord).
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(/,
     &' Module electrique: proprietes physiques lues sur fichier',/)
 9900 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR DANS ELPHYV (MODULE ELECTRIQUE)      ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@  Tabulation echoue avec une temperature TP = ', E14.5       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      RETURN
      END
c@z
