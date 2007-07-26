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
C     Copyright (C) 1998-2007 EDF S.A., France
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
                        SUBROUTINE USELPH
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
CFONC    REMPLISSAGE DES VARIABLES PHYSIQUES POUR LE MODULE ELECTRIQUE
CFONC
CFONC      ----> Effet Joule
CFONC      ----> Arc Electrique
CFONC      ----> Conduction Ionique
CFONC
CFONC       1) Masse Volumique
CFONC       2) Viscosite moleculaire
CFONC       3) Chaleur massique Cp
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
CFONC  Pour le module electrique, toutes les proprietes physiques sont
CFONC    supposees variables et contenues dans le tableau PROPCE
CFONC    (meme si elles sont physiquement constantes)
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
CFONC       - dans usiniv/useliv :
CFONC          . les variables de calcul  (initialisees a 0 par defaut
CFONC              ou a l
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
CFONC  Il est conseille de ne garder dans ce sous programme que
CFONC     le strict necessaire.
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
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IEL   , IPHAS
      INTEGER          IPCROM, IPCVIS, IPCCP , IPCVSL, IPCSIG
      INTEGER          MODE
C
      DOUBLE PRECISION TP
      DOUBLE PRECISION XKR   , XBR
      DOUBLE PRECISION ROM0  , TEMP0 , DILAR , AA    , BB    , CC
      DOUBLE PRECISION SRROM1, RHONP1
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C***********************************************************************
C
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
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
C       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
C       ET NON CELLE DE LA BIBLIOTHEQUE
C=======================================================================
C
C
C     En Joule, on s'arrete : il faut que l'utilisateur
C       donne les proprietes physiques
      IF ( IPPMOD(IELJOU).GE.1 ) THEN
C
        WRITE(NFECRA,9010)
        CALL CSEXIT (1)
C
C     En Arc on continue car on a un fichier de donnees
C       Un message indique que l'utilisateur n'a rien fourni
      ELSEIF(IPPMOD(IELARC).GE.1) THEN
C
        IF(IPASS.EQ.1) THEN
          WRITE(NFECRA,9011)
        ENDIF
C
        RETURN
C
      ENDIF
C
 9010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A LA DEFINITION DES PROP. PHYSIQUES   ',/,
     &'@    *********                                               ',/,
     &'@                      MODULE ELECTRIQUE                     ',/,
     &'@                                                            ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR uselph DOIT ETRE COMPLETE',/,
     &'@                                                            ',/,
     &'@     Ce sous-programme utilisateur permet de definir les    ',/,
     &'@       proprietes physiques. Il est indispensable.          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9011 FORMAT(/,
     &' Module arc electrique: pas d''intervention utilisateur pour ',/,
     &'                          le calcul des proprietes physiques.',/)
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C     Message au premier passage pour indiquer que l'utilisateur a
C       rapatrie le sous-programme.
      IF(IPASS.EQ.1) THEN
        WRITE(NFECRA,1000)
      ENDIF
C
C=======================================================================
C 1 - EFFET JOULE
C=======================================================================
C
      IF ( IPPMOD(IELJOU).GE.1 ) THEN
C
C
C     Attention, dans les modules electriques, la chaleur massique, la
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
C
C
C
C
C       Calcul de la temperature a partir de l'enthalpie
C       ------------------------------------------------
C
C       Ceci depend largement des choix utilisateur en
C         matiere de loi H-T (T en Kelvin)
C
C       On demande de fournir cette loi dans le sous programme usthht
C          (USERS/base/usthht.F)
C           usthht fournit en particulier un exemple d'interpolation
C            a partir d'une tabulation utilisateur
C           usthht en mode T->H sera utilise pour l'initialisation
C            de l'enthalpie dans useliv.
C
C       MODE = 1 : H=RTP(IEL,ISCA(IHM)) -> T=PROPCE(IEL,IPPROC(ITEMP))
        MODE = 1
C
        DO IEL = 1, NCEL
          CALL USTHHT (MODE,
     &         RTP(IEL,ISCA(IHM)),PROPCE(IEL,IPPROC(ITEMP)))
        ENDDO
C
C
C       Masse volumique au centre des cellules
C       --------------------------------------
C
C     ATTENTION :
C     =========
C       Dans le module electrique effet Joule, on fournira
C       OBLIGATOIREMENT la loi de variation de la masse volumique ici
C       en renseignant PROPCE(IEL,IPCROM)
C       (meme si elle est uniforme ou constante).
C
C
C     Masse Vol : RO = ROM0 / (1+DILAR*(T-T0)
C         (Choudhary) semblable a Plard (HE-25/94/017)
C
C          avec sous-relaxation (sauf au premier pas de temps)
C
        TEMP0  = 300.D0
        ROM0   = 2500.D0
        DILAR  = 7.5D-5
        IF(NTCABS.GT.1) THEN
          SRROM1 = SRROM
        ELSE
          SRROM1 = 0.D0
        ENDIF
C
        IPCROM = IPPROC(IROM(IPHAS))
        DO IEL = 1, NCEL
          RHONP1 = ROM0 /
     &            (1.D0+ DILAR * (PROPCE(IEL,IPPROC(ITEMP))-TEMP0) )
          PROPCE(IEL,IPCROM) =
     &         SRROM1*PROPCE(IEL,IPCROM)+(1.D0-SRROM1)*RHONP1
        ENDDO
C
C
C       Viscosite moleculaire dynamique en kg/(m s)
C        ------------------------------------------
C
C     ATTENTION :
C     =========
C       Dans le module electrique effet Joule, on fournira
C       OBLIGATOIREMENT la loi de variation de la viscosite ici
C       en renseignant PROPCE(IEL,IPCVIS)
C       (meme si elle est uniforme ou constante).
C
C
C     Viscosite : MU = EXP((AA/T-BB)-CC)
C          (Choudhary)
C      Plard (HE-25/94/017) ; limite a 1173K par C Delalondre
C
        IPCVIS = IPPROC(IVISCL(IPHAS))
        AA     = 10425.D0
        BB     =   500.D0
        CC     =-6.0917D0
C
        DO IEL = 1, NCEL
          IF ( PROPCE(IEL,IPPROC(ITEMP)) .GT. 1173.D0 ) THEN
            TP = PROPCE(IEL,IPPROC(ITEMP))
          ELSE
            TP= 1173.D0
          ENDIF
          PROPCE(IEL,IPCVIS) = EXP( (AA/(TP-BB))+CC )
        ENDDO
C
C
C       Chaleur specifique J/(kg degres)
C       --------------------------------
C
C     ATTENTION :
C     =========
C       Dans le module electrique effet Joule, on fournira
C       OBLIGATOIREMENT la loi de variation de la chaleur massique ici
C       en renseignant PROPCE(IEL,IPCPP)
C       (meme si elle est uniforme ou constante).
C
C
C        CP = 1381 (Choudhary)
C          coherent avec Plard (HE-25/94/017)
C
        IPCCP  = IPPROC(ICP(IPHAS))
        DO IEL = 1, NCEL
          PROPCE(IEL,IPCCP) = 1381.D0
        ENDDO
C
C
C       Lambda/Cp en kg/(m s)
C       ---------------------
C
C     ATTENTION :
C     =========
C       Dans le module electrique effet Joule, on fournira
C       OBLIGATOIREMENT la loi de variation de la conductivite ici
C       en renseignant PROPCE(IEL,IPCVSL)
C       (meme si elle est uniforme ou constante).
C
C
C         Lambda
C          On suppose Cp renseigne au prealable.
C
C          Plard (HE-25/94/017)
C
        IPCVSL = IPPROC(IVISLS(ISCALT(IPHAS)))
C
        DO IEL = 1, NCEL
          XBR = 85.25D0
     &         -5.93D-2*(PROPCE(IEL,IPPROC(ITEMP))-TKELVI)
     &         +2.39D-5*(PROPCE(IEL,IPPROC(ITEMP))-TKELVI)**2
          XKR = 16.D0*STEPHN*(1.4D0)**2*(PROPCE(IEL,IPPROC(ITEMP)))**3
     &         /(3.D0*XBR)
C
          PROPCE(IEL,IPCVSL) = 1.73D0 + XKR
        ENDDO
C
C --- On utilise CP calcule  dans PROPCE ci dessus
        DO IEL = 1, NCEL
          PROPCE(IEL,IPCVSL) = PROPCE(IEL,IPCVSL)/PROPCE(IEL,IPCCP)
        ENDDO
C
C
C       Conductivite electrique en S/m
C       ==============================
C
C     ATTENTION :
C     =========
C       Dans le module electrique effet Joule, on fournira
C       OBLIGATOIREMENT la loi de variation de la conductivite ici
C       en renseignant PROPCE(IEL,IPCSIG)
C       (meme si elle est uniforme ou constante).
C
C
C         SIGMA  (Plard HE-25/94/017)
C
        IPCSIG = IPPROC(IVISLS(IPOTR))
        DO IEL = 1, NCEL
          PROPCE(IEL,IPCSIG) =
     &         EXP(7.605D0-7200.D0/PROPCE(IEL,IPPROC(ITEMP)))
        ENDDO
C
C     La conductivite electrique pour le potentiel imaginaire est
C       toujours implicitement prise egale a la conductivite
C       utilisee pour le potentiel reel.
C       IL NE FAUT PAS la renseigner.
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C     Conductivite electrique imaginaire :
C     La conductivite reelle et imaginaire sont dans le meme tableau.
C       Ce choix est fait en dur dans varpos.
C       Les pointeurs pour les deux existent quand meme.
C     Sinon, on pourrait faire ceci :
        IF(1.EQ.0) THEN
          IF ( IPPMOD(IELJOU).EQ.2 .OR. IPPMOD(IELJOU).EQ.4 ) THEN
            DO IEL = 1, NCEL
              PROPCE(IEL,IPPROC(IVISLS(IPOTI))) =
     &             PROPCE(IEL,IPPROC(IVISLS(IPOTR)))
            ENDDO
          ENDIF
        ENDIF
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C
C       Diffusivite variable a l'exclusion de l'enthalpie et du potentiel
C       -----------------------------------------------------------------
C     Pour le moment, il n'y a pas d'autres scalaires et
C                                                  on ne fait donc rien
C
      ENDIF
C
C=======================================================================
C 2 - ARC ELECTRIQUE
C=======================================================================
C
C     Les proprietes physiques sont a priori fournies par fichier
C       de donnees. IL n'y a donc rien a faire ici.
C
C      IF ( IPPMOD(IELARC).GE.1 ) THEN
C      ENDIF
C
C
C=======================================================================
C 3 - CONDUCTION IONIQUE
C=======================================================================
C
C     CETTE OPTION N'EST PAS ACTIVABLE
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(/,
     &' Module electrique: intervention utilisateur pour        ',/,
     &'                      le calcul des proprietes physiques.',/)
C
C----
C FIN
C----
C
      RETURN
      END
c@z
