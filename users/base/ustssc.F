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
                        SUBROUTINE USTSSC
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCAL  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   CRVEXP , CRVIMP ,
     &   VISCF  , VISCB  , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     , W10    , W11    ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC ROUTINE UTILISATEUR
CFONC    ON PRECISE LES TERMES SOURCES UTILISATEURS
CFONC    POUR UN SCALAIRE SUR UN PAS DE TEMPS
CFONC
CFONC ON RESOUT RHO*VOLUME*D(VAR)/DT = CRVIMP*VAR + CRVEXP
CFONC
CFONC ON FOURNIT ICI CRVIMP ET CRVEXP (ILS CONTIENNENT RHO*VOLUME)
CFONC     CRVEXP en kg variable/s :
CFONC      ex : pour les temperatures      kg degres/s
CFONC           pour les enthalpies        Joules/s
CFONC     CRVIMP en kg /s :
CFONC
CFONC VEILLER A UTILISER UN CRVIMP NEGATIF
CFONC  (ON IMPLICITERA CRVIMP
CFONC   IE SUR LA DIAGONALE DE LA MATRICE, LE CODE AJOUTERA :
CFONC    MAX(-CRVIMP,0) EN SCHEMA STANDARD EN TEMPS
CFONC        -CRVIMP    SI LES TERMES SOURCES SONT A L'ORDRE 2
CFONC
CFONC CES TABLEAUX SONT INITIALISES A ZERO AVANT APPEL A CE SOUS
CFONC    PROGRAMME ET AJOUTES ENSUITE AUX TABLEAUX PRIS EN COMPTE
CFONC    POUR LA RESOLUTION
CFONC
CFONC  EN CAS D'ORDRE 2 DEMANDE SUR LES TERMES SOURCES, ON DOIT
CFONC    FOURNIR CRVEXP A L'INSTANT N     (IL SERA EXTRAPOLE) ET
CFONC            CRVIMP A L'INSTANT N+1/2 (IL EST  DANS LA MATRICE,
CFONC                                      ON LE SUPPOSE NEGATIF)
CFONC
CFONC L'identification des faces de bord concernees se fait grace
CFONC a la commande GETFBR.
CFONC
CFONC
CFONC IDENTIFICATION DES CELLULES/FACES DE BORD/FACES INTERNES
CFONC ========================================================
CFONC
CFONC   Les commandes GETCEL, GETFBR et GETFAC permettent d'identifier
CFONC   respectivement les cellules, faces ou faces de bord en fonction
CFONC   de differents criteres.
CFONC
CFONC   GETCEL(CHAINE,NLELT,LSTELT) :
CFONC   - CHAINE est une chaine de caractere fournie par l'utilisateur
CFONC     qui donne les criteres de selection
CFONC   - NLTELT est renvoye par la commande. C'est un entier qui
CFONC     correspond au nombre de cellules trouveees repondant au
CFONC     critere
CFONC   - LSTELT est renvoye par la commande. C'est un tableau d'entiers
CFONC     de taille NLTELT donnant la liste des cellules trouvees
CFONC     repondant au critere.
CFONC
CFONC   CHAINE peut etre constitue de :
CFONC   - references de couleurs (ex. : 1, 8, 26, ...
CFONC   - references de groupes (ex. : entrees, groupe1, ...)
CFONC   - criteres geometriques (ex. X<0.1, Y>=0.25, ...)
CFONC   Ces criteres peuvent etre combines par des operateurs logiques
CFONC   (AND et OR) et des parentheses
CFONC   ex. : '1 AND (groupe2 OR groupe3) AND Y<1' permettra de recuperer
CFONC   les cellules de couleur 1, appartenant aux groupes 'groupe2'
CFONC   ou 'groupe3' et de coordonnee Y inferieure a 1.
CFONC
CFONC   La syntaxe des commandes GETFBR et GETFAC est identique.
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
CARGU ! NPRFML       ! E  !  -> ! NOMBRE DE PROPRIETESE DES FAMILLES   !
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE SOMMETS                    !
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NCEPDP       ! E  !  -> ! NOMBRE DE CELLULES AVEC PDC          !
CARGU ! NCESMP       ! E  !  -> ! NOMBRE DE CELLULES A SOURCE DE MASSE !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! ISCAL        ! E  !  -> ! NUMERO DU SCALAIRE                   !
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
CARGU ! MAXELT       !  E !  -> ! NB MAX D'ELEMENTS (CELL,FAC,FBR)     !
CARGU ! LSTELT(MAXELT) TE !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! ICEPDC(NCELET! TE !  -> ! NUMERO DES NCEPDP CELLULES AVEC PDC  !
CARGU ! ICETSM(NCESMP! TE !  -> ! NUMERO DES CELLULES A SOURCE DE MASSE!
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
CARGU ! CKUPDC       ! TR !  -> ! TABLEAU DE TRAVAIL POUR PDC          !
CARGU !  (NCEPDP,6)  !    !     !                                      !
CARGU ! SMACEL       ! TR !  -> ! VALEUR DES VARIABLES ASSOCIEE A LA   !
CARGU ! (NCESMP,*   )!    !     !  SOURCE DE MASSE                     !
CARGU !              !    !     !  POUR IVAR=IPR, SMACEL=FLUX DE MASSE !
CARGU ! CRVEXP(NCELET! TR ! <-  ! TABLEAU DE TRAVAIL POUR PART EXPLICIT!
CARGU ! CRVIMP(NCELET! TR ! <-  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! VISCF(NFAC)  ! TR !  -  ! TABLEAU DE TRAVAIL    FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! TABLEAU DE TRAVAIL    FACES DE BORD  !
CARGU ! XAM(NFAC,2)  ! TR !  -  ! TABLEAU DE TRAVAIL    FACES DE BORD  !
CARGU ! W1..11(NCELET! TR !  -  ! TABLEAU DE TRAVAIL    CELLULES       !
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
      INCLUDE "paramx.h"
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
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
      INTEGER          NCEPDP , NCESMP
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          ISCAL
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          MAXELT, LSTELT(MAXELT)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICEPDC(NCEPDP)
      INTEGER          ICETSM(NCESMP), ITYPSM(NCESMP,NVAR)
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
      DOUBLE PRECISION CKUPDC(NCEPDP,6), SMACEL(NCESMP,NVAR)
      DOUBLE PRECISION CRVEXP(NCELET), CRVIMP(NCELET)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION XAM(NFAC,2)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION W10(NCELET), W11(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IVAR, IISCVR, IPCROM, IEL, IPHAS, IUTILE
      INTEGER          ILELT, NLELT
C
      DOUBLE PRECISION TAUF, PRODF, VOLF, PWATT
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C
      IF(1.EQ.1) RETURN
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C --- Numero du scalaire a traiter : ISCAL
C
C     ISCAL est une donnee d'entree de ce sous programme
C       qui indique quel scalaire on est en train de traiter
C       en effet ce sous programme est appele successivement pour
C       tous les scalaires du calcul.
C     ISCAL ne doit pas etre modifie par l'utilisateur
C
C     Si le calcul comporte plusieurs scalaires (de 1 a n par exemple)
C       et que l'on souhaite imposer des termes sources pour un
C       scalaire p donne, on doit utiliser ISCAL ici :
C       CRVIMP et CRVEXP ne seront renseignes que "IF(ISCAL.EQ.p)"
C
C
C --- Numero de la variable associee au scalaire a traiter ISCAL
      IVAR = ISCA(ISCAL)
C
C --- Nom de la variable associee au scalaire a traiter ISCAL
      CHAINE = NOMVAR(IPPRTP(IVAR))
C
C --- Indicateur de variance
C         Si ISCAVR = 0 :
C           le scalaire ISCAL n'est pas une variance
C         Si ISCAVR > 0 et ISCAVR < NSCAL + 1 :
C           le scalaire ISCAL est une variance associee
C           au scalaire ISCAVR
      IISCVR = ISCAVR(ISCAL)
C
C --- Numero de phase associee au scalaire ISCAL
      IPHAS = IPHSCA(ISCAL)
C
C --- Numero des grandeurs physiques (voir usclim)
      IPCROM = IPPROC(IROM(IPHAS))
C
      IF(IWARNI(IVAR).GE.1) THEN
        WRITE(NFECRA,1000) CHAINE(1:8)
      ENDIF
C
C
C=======================================================================
C 2. EXEMPLE DE TERME SOURCE ARBITRAIRE, POUR LA VARIABLE F :
C
C                             S = A * F + B
C
C            APPARAISSANT DANS LES EQUATIONS SOUS LA FORME :
C
C                       RHO VOLUME D(F)/Dt = VOLUME*S
C
C
C   CE TERME A UNE PARTIE QU'ON VEUT IMPLICITER         : A
C           ET UNE PARTIE QU'ON VA TRAITER EN EXPLICITE : B
C
C
C   ICI PAR EXEMPLE ARBITRAIREMENT :
C
C     A = - RHO / TAUF
C     B =   RHO * PRODF
C        AVEC
C     TAUF   = 10.D0  [secondes  ] (TEMPS DE DISSIPATION DE F)
C     PRODF  = 100.D0 [variable/s] (PRODUCTION DE F PAR UNITE DE TEMPS)
C
C   ON A ALORS
C     CRVIMP(IEL) = VOLUME(IEL)* A = - VOLUME(IEL) (RHO / TAUF )
C     CRVEXP(IEL) = VOLUME(IEL)* B =   VOLUME(IEL) (RHO * PRODF)
C
C   ON REMPLIT CI-DESSOUS CRVIMP ET CRVEXP CORRESPONDANTS.
C
C=======================================================================
C
C
C
C     ATTENTION, L'EXEMPLE EST COMPLETEMENT ARBITRAIRE
C     =========
C         ET DOIT ETRE REMPLACE PAR LES TERMES UTILISATEURS ADEQUATS
C
C
C ----------------------------------------------
C
C Il est assez courant que l'on oublie d'eliminer cet exemple
C   de ce sous-programme.
C On a donc prevu le test suivant pour eviter les mauvaises surprises
C
      IUTILE = 0
C
      IF(IUTILE.EQ.0) RETURN
C
C ----------------------------------------------
C
      TAUF  = 10.D0
      PRODF = 100.D0
C
      DO IEL = 1, NCEL
        CRVIMP(IEL) = - VOLUME(IEL)*PROPCE(IEL,IPCROM)/TAUF
      ENDDO
C
      DO IEL = 1, NCEL
        CRVEXP(IEL) =   VOLUME(IEL)*PROPCE(IEL,IPCROM)*PRODF
      ENDDO
C
C=======================================================================
C 3. EXEMPLE DE TERME SOURCE ARBITRAIRE, POUR LA VARIABLE F = ENTHALPIE :
C
C     IL S'AGIT D'UNE PARTICULARISATION DE L'EXEMPLE PRECEDENT AVEC
C
C                             S = B
C
C            APPARAISSANT DANS LES EQUATIONS SOUS LA FORME :
C
C                       RHO VOLUME D(F)/Dt = VOLUME*S
C
C
C   IL N'Y A RIEN A IMPLICITER, ON PEUT DONC IMPOSER
C
C     CRVIMP(IEL) = 0.D0
C
C
C   CE TERME A UNE PARTIE QU'ON VA TRAITER EN EXPLICITE : B
C
C     SI F EST UNE ENTHALPIE (J/kg), ON AURA ALORS B EN Watt/m3
C
C     SI ON CONNAIT LA PUISSANCE PWATT A INTRODUIRE UNIFORMEMENT DANS UN
C       VOLUME DONNE VOLF TEL QUE 0<X<1.2, 3.1<Y<4, ON POURRA ALORS
C       ECRIRE, POUR LES CELLULES DU MAILLAGE CONCERNEES :
C
C       CRVEXP(IEL) = VOLUME(IEL)* B =   VOLUME(IEL)*(PWATT/VOLF)
C
C
C=======================================================================
C
C
C
C     ATTENTION, L'EXEMPLE EST COMPLETEMENT ARBITRAIRE
C     =========
C         ET DOIT ETRE REMPLACE PAR LES TERMES UTILISATEURS ADEQUATS
C
C
C ----------------------------------------------
C
C Il est assez courant que l'on oublie d'eliminer cet exemple
C   de ce sous-programme.
C On a donc prevu le test suivant pour eviter les mauvaises surprises
C
      IUTILE = 0
C
      IF(IUTILE.EQ.0) RETURN
C
C ----------------------------------------------
C
C PUISSANCE A INTRODUIRE DANS VOLF
C   ATTENTION on suppose qu'on travaille en enthalpie,
C             si on travaille en temperature, il ne faut pas oublier de
C               diviser PWATT par Cp.
C
      PWATT = 100.D0
C
C CALCUL DE VOLF
C
      VOLF  = 0.D0
      CALL GETCEL('X > 0.0 and X < 1.2 and Y > 3.1 and'//
     &            'Y < 4.0',NLELT,LSTELT)
C
      DO ILELT = 1, NLELT
        IEL = LSTELT(ILELT)
        VOLF = VOLF + VOLUME(IEL)
      ENDDO
C
      DO ILELT = 1, NLELT
        IEL = LSTELT(ILELT)
C PAS DE TERME IMPLICITE
        CRVIMP(IEL) = 0.D0
C PUISSANCE EN EXPLICITE
        CRVEXP(IEL) = VOLUME(IEL)*PWATT/VOLF
      ENDDO
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(' TERMES SOURCES UTILISATEURS POUR LA VARIABLE ',A8,/)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
