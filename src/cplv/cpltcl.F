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
                        SUBROUTINE CPLTCL
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , ITRIFB , ITYPFB , IZFPPP ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RCODCL ,
     &   W1     , W2     , W3     , W4     , W5     , W6     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
CFONC    --------------------------------------------------------------
CFONC
CFONC     ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE
CFONC
CFONC       COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
CFONC       TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON
CFONC
CFONC       CONDITIONS AUX LIMITES AUTOMATIQUES
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
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC            !
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !                                      !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR            !
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !                                      !
CARGU ! ICODCL       ! TE ! <-  ! CODE DE CONDITION LIMITES AUX FACES  !
CARGU !  (NFABOR,NVAR!    !     !  DE BORD                             !
CARGU !              !    !     ! = 1   -> DIRICHLET                   !
CARGU !              !    !     ! = 3   -> DENSITE DE FLUX             !
CARGU !              !    !     ! = 4   -> GLISSEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 5   -> FROTTEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 6   -> RUGOSITE ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 9   -> ENTREE/SORTIE LIBRE (VITESSE!
CARGU !              !    !     !  ENTRANTE EVENTUELLE     BLOQUEE     !
CARGU ! ITRIFB(NFABOR! TE !  -> ! INDIRECTION POUR TRI DES FACES DE BRD!
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS      )!    !     !                                      !
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
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS                !
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
CARGU ! RCODCL       ! TR ! <-  ! VALEUR DES CONDITIONS AUX LIMITES    !
CARGU !  (NFABOR,NVAR!    !     !  AUX FACES DE BORD                   !
CARGU !              !    !     ! RCODCL(1) = VALEUR DU DIRICHLET      !
CARGU !              !    !     ! RCODCL(2) = VALEUR DU COEF. D'ECHANGE!
CARGU !              !    !     !  EXT. (INFINIE SI PAS D'ECHANGE)     !
CARGU !              !    !     ! RCODCL(3) = VALEUR DE LA DENSITE DE  !
CARGU !              !    !     !  FLUX (NEGATIF SI GAIN) W/m2 OU      !
CARGU !              !    !     !  HAUTEUR DE RUGOSITE (m) si ICODCL=6 !
CARGU !              !    !     ! POUR LES VITESSES (VISTL+VISCT)*GRADU!
CARGU !              !    !     ! POUR LA PRESSION             DT*GRADP!
CARGU !              !    !     ! POUR LES SCALAIRES                   !
CARGU !              !    !     !        CP*(VISCLS+VISCT/SIGMAS)*GRADT!
CARGU ! W1,2,3,4,5,6 ! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU !  (NCELET     !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
CARGU ! COEFU        ! TR !  -  ! TAB DE TRAV                          !
CARGU !  (NFABOR,3)  !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
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
C ARGUMENTS
C
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
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
      INTEGER          ICODCL(NFABOR,NVAR)
      INTEGER          ITRIFB(NFABOR,NPHAS), ITYPFB(NFABOR,NPHAS)
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
      DOUBLE PRECISION RCODCL(NFABOR,NVAR,3)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION W4(NCELET),W5(NCELET),W6(NCELET)
      DOUBLE PRECISION COEFU(NFABOR,NDIM)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          II, IPHAS, IFAC, IZONE, MODE, IEL, IGE, IOK
      INTEGER          ICHA, IPBROM, ICKE, IPCVIS
      INTEGER          NBRVAL
      DOUBLE PRECISION QISQC, VISCLA, D2S3, UREF2, RHOMOY, DHY, XITURB
      DOUBLE PRECISION USTAR2, XKENT, XEENT, T1
      DOUBLE PRECISION H1    (NOZPPM)
      DOUBLE PRECISION QIMPC (NOZPPM) , QCALC(NOZPPM)
      DOUBLE PRECISION COEFE (NGAZEM)
      DOUBLE PRECISION F1MC  (NCHARM) , F2MC (NCHARM)
C
C***********************************************************************
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IPHAS = 1
      IPBROM = IPPROB(IROM  (IPHAS))
      IPCVIS = IPPROC(IVISCL(IPHAS))
C
      D2S3 = 2.D0/3.D0
C
C=======================================================================
C 1.  ECHANGES EN PARALLELE POUR LES DONNEES UTILISATEUR
C=======================================================================
C
C  En realite on pourrait eviter cet echange en modifiant uscpcl et en
C    demandant a l'utilisateur de donner les grandeurs dependant de la
C    zone hors de la boucle sur les faces de bord : les grandeurs
C    seraient ainsi disponibles sur tous les processeurs. Cependant,
C    ca rend le sous programme utilisateur un peu plus complique et
C    surtout, si l'utilisateur le modifie de travers, ca ne marche pas.
C  On suppose que toutes les grandeurs fournies sont positives, ce qui
C    permet d'utiliser un max pour que tous les procs les connaissent.
C    Si ce n'est pas le cas, c'est plus complique mais on peut s'en tirer
C    avec un max quand meme.
C
      IF(IRANGP.GE.0) THEN
        CALL PARIMX(NOZAPM,IQIMP )
C       ===========
        CALL PARIMX(NOZAPM,IENTAT)
C       ===========
        CALL PARIMX(NOZAPM,IENTCP)
C       ===========
        CALL PARRMX(NOZAPM,QIMPAT)
C       ===========
        CALL PARRMX(NOZAPM,TIMPAT)
C       ===========
        NBRVAL = NOZPPM*NCHARM
        CALL PARRMX(NBRVAL,QIMPCP)
C       ===========
        NBRVAL = NOZPPM*NCHARM
        CALL PARRMX(NBRVAL,TIMPCP)
C       ===========
        NBRVAL = NOZPPM*NCHARM*NCLCPM
        CALL PARRMX(NBRVAL,DISTCH)
C       ===========
      ENDIF
C
C
C=======================================================================
C 2.  CORRECTION DES VITESSES (EN NORME) POUR CONTROLER LES DEBITS
C     IMPOSES
C       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
C                     =========================
C=======================================================================
C
C --- Debit calcule
C
      DO IZONE = 1, NOZPPM
        QCALC(IZONE) = 0.D0
      ENDDO
      DO IFAC = 1, NFABOR
        IZONE = IZFPPP(IFAC)
        QCALC(IZONE) = QCALC(IZONE) - PROPFB(IFAC,IPBROM) *
     &      ( RCODCL(IFAC,IU(IPHAS),1)*SURFBO(1,IFAC) +
     &        RCODCL(IFAC,IV(IPHAS),1)*SURFBO(2,IFAC) +
     &        RCODCL(IFAC,IW(IPHAS),1)*SURFBO(3,IFAC) )
      ENDDO
C
      IF(IRANGP.GE.0) THEN
        CALL PARRSM(NOZAPM,QCALC )
      ENDIF
C
      DO IZONE = 1, NOZAPM
        IF ( IQIMP(IZONE).EQ.0 ) THEN
          QIMPC(IZONE) = QCALC(IZONE)
        ENDIF
      ENDDO
C
C --- Correction des vitesses en norme
C
      IOK = 0
      DO II = 1, NZFPPP
        IZONE = ILZPPP(II)
        IF ( IQIMP(IZONE).EQ.1 ) THEN
          IF(ABS(QCALC(IZONE)).LT.EPZERO) THEN
            WRITE(NFECRA,2001)IZONE,IQIMP(IZONE),QCALC(IZONE)
            IOK = IOK + 1
          ENDIF
        ENDIF
      ENDDO
      IF(IOK.NE.0) THEN
        CALL CSEXIT (1)
C       ===========
      ENDIF
      DO IFAC = 1, NFABOR
        IZONE = IZFPPP(IFAC)
        IF ( IQIMP(IZONE).EQ.1 ) THEN
          QIMPC(IZONE) = QIMPAT(IZONE)
          QISQC = QIMPC(IZONE)/QCALC(IZONE)
          RCODCL(IFAC,IU(IPHAS),1) = RCODCL(IFAC,IU(IPHAS),1)*QISQC
          RCODCL(IFAC,IV(IPHAS),1) = RCODCL(IFAC,IV(IPHAS),1)*QISQC
          RCODCL(IFAC,IW(IPHAS),1) = RCODCL(IFAC,IW(IPHAS),1)*QISQC
        ENDIF
C
      ENDDO
C
C
C
 2001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : MODULE PHYSIQUES PARTICULIERES              ',/,
     &'@    *********                                               ',/,
     &'@     COMBUSTION CHARBON PULVERISE COUPLE AU                 ',/,
     &'@     TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON :       ',/,
     &'@     PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,
     &'@                                                            ',/,
     &'@  Le debit est impose sur la zone IZONE =     ', I10         ,/,
     &'@    puisque                IQIMP(IZONE) =     ', I10         ,/,
     &'@  Or, sur cette zone, le produit RHO D S integre est nul :  ',/,
     &'@    il vaut                             = ',E14.5            ,/,
     &'@    (D est la direction selon laquelle est impose le debit).',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier uscpcl, et en particulier                        ',/,
     &'@    - que le vecteur  RCODCL(IFAC,IU(IPHAS),1),             ',/,
     &'@                      RCODCL(IFAC,IV(IPHAS),1),             ',/,
     &'@                      RCODCL(IFAC,IW(IPHAS),1) qui determine',/,
     &'@      la direction de la vitesse est non nul et n''est pas  ',/,
     &'@      uniformement perpendiculaire aux face d''entree       ',/,
     &'@    - que la surface de l''entree n''est pas nulle (ou que  ',/,
     &'@      le nombre de faces de bord dans la zone est non nul)  ',/,
     &'@    - que la masse volumique n''est pas nulle               ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C=======================================================================
C 3.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
C       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
C                     =========================
C         ON DETERMINE LA FAMILLE ET SES PROPRIETES
C           ON IMPOSE LES CONDITIONS AUX LIMITES
C           POUR LA TURBULENCE
C
C=======================================================================
C
      DO IFAC = 1, NFABOR
C
        IZONE = IZFPPP(IFAC)
C
C      ELEMENT ADJACENT A LA FACE DE BORD
C
        IF ( ITYPFB(IFAC,IPHAS).EQ.IENTRE ) THEN
C
C ----  Traitement automatique de la turbulence
C
          IF ( ICALKE(IZONE).NE.0 ) THEN
C
C       La turbulence est calculee par defaut si ICALKE different de 0
C          - soit a partir du diametre hydraulique, d'une vitesse
C            de reference adaptes a l'entree courante si ICALKE = 1
C          - soit a partir du diametre hydraulique, d'une vitesse
C            de reference et de l'intensite turvulente
C            adaptes a l'entree courante si ICALKE = 2
C
            UREF2 = RCODCL(IFAC,IU(IPHAS),1)**2
     &            + RCODCL(IFAC,IV(IPHAS),1)**2
     &            + RCODCL(IFAC,IW(IPHAS),1)**2
            UREF2 = MAX(UREF2,1.D-12)
            RHOMOY = PROPFB(IFAC,IPBROM)
            IEL    = IFABOR(IFAC)
            VISCLA = PROPCE(IEL,IPCVIS)
            ICKE   = ICALKE(IZONE)
            DHY    = DH(IZONE)
            XITURB = XINTUR(IZONE)
            USTAR2 = 0.D0
            XKENT = EPZERO
            XEENT = EPZERO
            IF (ICKE.EQ.1) THEN
              CALL KEENDB
C             ===========
     &        ( UREF2, DHY, RHOMOY, VISCLA, CMU, XKAPPA,
     &          USTAR2, XKENT, XEENT )
            ELSE IF (ICKE.EQ.2) THEN
              CALL KEENIN
C             ===========
     &        ( UREF2, XITURB, DHY, CMU, XKAPPA, XKENT, XEENT )
            ENDIF
C
            IF (ITYTUR(IPHAS).EQ.2) THEN
C
              RCODCL(IFAC,IK(IPHAS),1)  = XKENT
              RCODCL(IFAC,IEP(IPHAS),1) = XEENT
C
            ELSEIF (ITYTUR(IPHAS).EQ.3) THEN
C
              RCODCL(IFAC,IR11(IPHAS),1) = D2S3*XKENT
              RCODCL(IFAC,IR22(IPHAS),1) = D2S3*XKENT
              RCODCL(IFAC,IR33(IPHAS),1) = D2S3*XKENT
              RCODCL(IFAC,IR12(IPHAS),1) = 0.D0
              RCODCL(IFAC,IR13(IPHAS),1) = 0.D0
              RCODCL(IFAC,IR23(IPHAS),1) = 0.D0
              RCODCL(IFAC,IEP(IPHAS),1)  = XEENT
C
            ELSEIF (ITURB(IPHAS).EQ.50) THEN
C
              RCODCL(IFAC,IK(IPHAS),1)   = XKENT
              RCODCL(IFAC,IEP(IPHAS),1)  = XEENT
              RCODCL(IFAC,IPHI(IPHAS),1) = D2S3
              RCODCL(IFAC,IFB(IPHAS),1)  = 0.D0
C
            ELSEIF (ITURB(IPHAS).EQ.60) THEN
C
              RCODCL(IFAC,IK(IPHAS),1)   = XKENT
              RCODCL(IFAC,IOMG(IPHAS),1) = XEENT/CMU/XKENT
C
            ENDIF
C
          ENDIF
C
        ENDIF
C
      ENDDO
C
C=======================================================================
C 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
C       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
C                     =========================
C         ON DETERMINE LA FAMILLE ET SES PROPRIETES
C           ON IMPOSE LES CONDITIONS AUX LIMITES
C           POUR LES SCALAIRES
C=======================================================================
C
      DO II = 1, NZFPPP
C
        IZONE = ILZPPP(II)
C
C Une entree IENTRE est forcement du type
C            IENTAT = 1 ou IENTCP = 1
C
        IF ( IENTAT(IZONE).EQ.1 ) THEN
C
C ------ Calcul de H1(IZONE)
C
          DO IGE = 1, NGAZEM
            COEFE(IGE) = ZERO
          ENDDO
          COEFE(IO2) = WMOLE(IO2) / (WMOLE(IO2)+XSI*WMOLE(IN2))
          COEFE(IN2) = 1.D0 - COEFE(IO2)
          DO ICHA = 1, NCHARM
            F1MC(ICHA) = ZERO
            F2MC(ICHA) = ZERO
          ENDDO
          T1   = TIMPAT(IZONE)
          MODE = -1
          CALL CPTHP1
C         ===========
     &    ( MODE  , H1(IZONE) , COEFE  , F1MC   , F2MC   ,
     &      T1    )
C
          ENDIF
C
      ENDDO
C
      DO IFAC = 1, NFABOR
C
        IZONE = IZFPPP(IFAC)
C
C      ELEMENT ADJACENT A LA FACE DE BORD
C
        IF ( ITYPFB(IFAC,IPHAS).EQ.IENTRE ) THEN
C
C ----  Traitement automatique des scalaires physiques particulieres
C
          DO ICHA = 1, NCHARB
C
C ------ CL pour F1M et F2M du charbon ICHA
            RCODCL(IFAC,ISCA(IF1M(ICHA)),1) = ZERO
            RCODCL(IFAC,ISCA(IF2M(ICHA)),1) = ZERO
C
          ENDDO
C
C ------ CL pour F3M
          RCODCL(IFAC,ISCA(IF3M),1) = ZERO
C ------ CL pour FP4M
          RCODCL(IFAC,ISCA(IF4P2M),1)   = ZERO
C ------ CL pour HM
          DO IGE = 1, NGAZEM
            COEFE(IGE) = ZERO
          ENDDO
          COEFE(IO2) = WMOLE(IO2) / (WMOLE(IO2)+XSI*WMOLE(IN2))
          COEFE(IN2) = 1.D0 - COEFE(IO2)
          DO ICHA = 1, NCHARM
            F1MC(ICHA) = ZERO
            F2MC(ICHA) = ZERO
          ENDDO
          T1   = TIMPAT(IZONE)
          MODE = -1
          CALL CPTHP1
C         ===========
     &    ( MODE  , H1(IZONE) , COEFE  , F1MC   , F2MC   ,
     &      T1    )
          RCODCL(IFAC,ISCA(IHM),1) = H1(IZONE)
C
        ENDIF
C
      ENDDO
C
C
C----
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
