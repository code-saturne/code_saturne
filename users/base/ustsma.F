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
                        SUBROUTINE USTSMA
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   NCESMP , IPHAS  , IAPPEL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC                     TERMES SOURCES DE MASSE
CFONC                      POUR LA PHASE IPHAS
CFONC
CFONC IAPPEL = 1 :
CFONC              CALCUL DU NOMBRE DE CELLULES AVEC SOURCE DE MASSE
CFONC               NCESMP
CFONC IAPPEL = 2 :
CFONC              REPERAGE DES CELLULES AVEC SOURCE DE MASSE
CFONC               ICETSM(NCESMP)
CFONC IAPPEL = 3 :
CFONC              CALCUL DES VALEURS DES COEFS DE SOURCE DE MASSE
CFONC
CFONC              L'equation de conservation de la masse devient :
CFONC
CFONC              d(rho)/dt + div(rho u) = GAMMA
CFONC                                  -
CFONC
CFONC              L'equation d'une variable f devient :
CFONC
CFONC              d(f)/dt = ..... + GAMMA*(f_i - f)
CFONC
CFONC              discretisee en :
CFONC
CFONC              RHO*(f^(n+1) - f^(n))/DT = .....
CFONC                                         + GAMMA*(f_i - f^(n+1))
CFONC
CFONC
CFONC              Deux possibilites pour chaque variable f :
CFONC                 - Flux de masse a la valeur de f ambiante
CFONC                            --> f_i = f^(n+1)
CFONC                    (l'equation de f n'est alors pas modifiee)
CFONC                 - Flux de masse avec une valeur donnee pour f
CFONC                            --> f_i specifie par l'utilisateur
CFONC
CFONC
CFONC
CFONC
CFONC             VARIABLES A REMPLIR PAR L'UTILISATEUR :
CFONC             =======================================
CFONC
CFONC              NCESMP : Nombre de cellules a source de masse
CFONC
CFONC              ICETSM(IELTSM) : Numero de la IELTSMieme cellule a
CFONC                                 source de masse (IELTSM<=NCESMP)
CFONC
CFONC              SMACEL(IEL,IPR(IPHAS)) : Valeur du flux de masse
CFONC                                  GAMMA (en kg/m^3/s)
CFONC                                  dans la IELieme cellule a source
CFONC                                  de masse
CFONC
CFONC              ITYPSM(IEL,IVAR) : type de flux associe a la variable
CFONC                                   IVAR dans la IELeme cellule a
CFONC                                   source de masse (pour toutes les
CFONC                                   variables sauf IVAR=IPR(IPHAS))
CFONC                  + ITYPSM = 0 --> injection a la valeur ambiante
CFONC                                   de IVAR
CFONC                  + ITYPSM = 1 --> injection a une valeur donnee
CFONC                                   de IVAR
CFONC
CFONC              SMACEL(IEL,IVAR) : Valeur de f_i pour la variable
CFONC                                 IVAR (pour toutes les variables
CFONC                                 sauf IVAR=IPR(IPHAS))
CFONC
CFONC
CFONC             REMARQUES
CFONC             =========
CFONC             * Si ITYPSM(IEL,IVAR)=0, SMACEL(IEL,IVAR)
CFONC                n'est pas utilise
CFONC             * Si SMACEL(IEL,IPR(IPHAS))<0,
CFONC                on enleve de la masse au systeme, Code_Saturne
CFONC                utilise donc automatiquement f_i=f^(n+1)
CFONC                QUELLES QUE SOIENT LES VALEURS DE ITYPSM(IEL)
CFONC                et SMACEL(IEL,IVAR)
CFONC
CFONC             * Si une variable n'est pas reliee a la phase iphas
CFONC                pour laquelle les informations precedentes ont ete
CFONC                completees, aucun terme source ne sera impose.
CFONC
CFONC             * Pour un scalaire qui ne repondrait pas a l'equation
CFONC                d(rho f)/dt + d(rho U f)/dx = ...
CFONC                (champ convecteur different par exemple)
CFONC                il est incorrect d'imposer le terme source de masse
CFONC                comme fait ici (sauf en cas d'injection a la valeur
CFONC                ambiante). Il faut introduire directement le terme
CFONC                source comme un terme source scalaire dans ustssc.
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
CARGU ! NCKPDP       ! E  !  -> ! NBR DE COEF DU TENSEUR DE PDC (3 OU 6!
CARGU ! NCESMP       ! E  ! <-> ! NOMBRE DE CELLULES AVEC TSM          !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE PHASE                      !
CARGU ! IAPPEL       ! E  !  -> ! INDIQUE LES DONNES A RENVOYER        !
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
CARGU ! ICETSM(NCESMP! TE ! <-> ! NUMERO DES CELLULES A SOURCE DE MASSE!
CARGU ! ITYPSM       ! TE ! <-> ! TYPE DE SOURCE DE MASSE POUR LES     !
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
CARGU ! RTPA         ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT            PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! CKUPDC(NCEPDP! TR !  -> ! TABLEAU DE TRAVAIL POUR PDC          !
CARGU !     , NCKPDP)!    !     !                                      !
CARGU ! SMACEL       ! TR ! <-> ! VALEUR DES VARIABLES ASSOCIEE A LA   !
CARGU ! (NCESMP,NVAR)!    !     !  SOURCE DE MASSE                     !
CARGU !              !    !     ! POUR IVAR=IPR, SMACEL=FLUX DE MASSE  !
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
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
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
      INTEGER          NCEPDP , NCKPDP , NCESMP
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IPHAS  , IAPPEL
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
      DOUBLE PRECISION DT(NCELET), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP)
      DOUBLE PRECISION SMACEL(NCESMP,NVAR)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IELTSM
      INTEGER          IFAC, IFML, ICOUL, IUTILE, II
      INTEGER          ILELT, NLELT
C
      DOUBLE PRECISION VENT, VENT2
      DOUBLE PRECISION DH, USTAR2
      DOUBLE PRECISION XKENT, XEENT
      DOUBLE PRECISION FLUCEL
      DOUBLE PRECISION VTOT  , GAMMA
C
C***********************************************************************
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IF(IAPPEL.EQ.1.OR.IAPPEL.EQ.2) THEN
C
C=======================================================================
C 1. POUR CHAQUE PHASE : UN OU DEUX APPELS
C
C    PREMIER APPEL :
C
C        IAPPEL = 1 : NCESMP : CALCUL DU NOMBRE DE CELLULES
C                                AVEC TERME SOURCE DE MASSE
C
C
C    DEUXIEME APPEL (POUR LES PHASES AVEC NCESMP > 0) :
C
C        IAPPEL = 2 : ICETSM : REPERAGE DU NUMERO DES CELLULES
C                                AVEC TERME SOURCE DE MASSE
C
C
C REMARQUES :
C
C        Ne pas utiliser SMACEL dans cette section
C          (il est rempli au troisieme appel, IAPPEL = 3)
C
C        Ne pas utiliser ICETSM dans cette section
C           au premier appel (IAPPEL = 1)
C
C        On passe ici a chaque pas de temps
C           (ATTENTION au cout calcul de vos developpements)
C
C=======================================================================
C
C ---------------------------------------------------------------------
C  1.1 A completer par l'utilisateur : selection des cellules
C  -----------------------------------------------------------
C
C Exemple 1 : Aucune source de masse (defaut)
        IELTSM = 0
C
C
C Exemple 2 : Sources de masse pour la phase 1
C               dans les cellules ayant une face de couleur 6
C            et dans les cellules dont le centre est a x < 2
C
C     Dans ce test en deux parties, il faut faire attention a
C       ne pas compter les cellules deux fois (une cellule dont une
C       face est couleur 6 peut aussi etre a x < 2, mais il ne
C       faut la compter qu'une seule fois). Il faut egalement noter
C       que lors du premier passage, le tableau ICETSM n'existe pas
C       encore (ne pas l'utiliser en dehors des tests IAPPEL.EQ.2).
C
        IUTILE = 0
        IF(IUTILE.EQ.1) THEN
C
          IF(IPHAS.EQ.1) THEN
C
            IELTSM = 0
C
C     Cellules dont le centre est a 250 < x < 500
C
            CALL GETCEL('X < 500.0 and X > 250.0',NLELT,LSTELT)
            DO ILELT = 1, NLELT
              II = LSTELT(ILELT)
              IELTSM = IELTSM + 1
              IF (IAPPEL.EQ.2) ICETSM(IELTSM) = II
            ENDDO
C
C
C     Cellules de bord dont une face est de couleur 3
C
            CALL GETFBR('3',NLELT,LSTELT)
            DO ILELT = 1, NLELT
              IFAC = LSTELT(ILELT)
              II   = IFABOR(IFAC)
C       On ne comptabilise cette cellule que si on ne l'a pas deja vue
C         pour cela on prend la negation du test precedent
              IF (.NOT.(XYZCEN(1,II).LT.500.D0.AND.
     &                  XYZCEN(1,II).GT.250.D0)    )THEN
                IELTSM = IELTSM + 1
                IF (IAPPEL.EQ.2) ICETSM(IELTSM) = II
              ENDIF
            ENDDO
C
          ELSE
            IELTSM = 0
          ENDIF
C
        ENDIF
C
C ---------------------------------------------------------------------
C  1.2 Sous section generique a ne pas modifier
C  ---------------------------------------------
C
C --- Pour IAPPEL = 1,
C      Renseigner NCESMP, nombre de cellules avec terme source de masse
C      Le bloc ci dessous est valable pourles 2 exemples ci dessus
C
        IF (IAPPEL.EQ.1) THEN
          NCESMP = IELTSM
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(IAPPEL.EQ.3) THEN
C
C=======================================================================
C
C 2. POUR CHAQUE PHASE AVEC NCESMP > 0 , TROISIEME APPEL
C
C      TROISIEME APPEL (POUR LES PHASES AVEC NCESMP > 0) :
C
C       IAPPEL = 3 : ITYPSM : TYPE DE SOURCE DE MASSE
C                    SMACEL : SOURCE DE MASSE
C
C
C   REMARQUE :
C
C      ATTENTION, Si on positionne ITYPSM(IEL,IVAR) A 1, IL FAUT
C      egalement renseigner SMACEL(IEL,IVAR)
C
C=======================================================================
C
C
C ---------------------------------------------------------------------
C  2.1 A completer par l'utilisateur ITYPSM et SMACEL
C  -----------------------------------------------------
C
C Exemple 1 : simulation d'une entree par des termes de source de masse
C           et affichage du flux de masse total pour la phase 1
C
        IF(IPHAS.EQ.1) THEN
C
          VENT = 0.1D0
          VENT2 = VENT**2
          DH     = 0.5D0
C         Calcul de la vitesse de frottement au carre (USTAR2)
C           et de k et epsilon en entree (XKENT et XEENT) a partir
C           de lois standards en conduite circulaire
C           (leur initialisation est inutile mais plus propre)
          USTAR2 = 0.D0
          XKENT  = EPZERO
          XEENT  = EPZERO
C
          CALL KEENDB
C         ===========
     &      ( VENT2, DH, RO0(IPHAS), VISCL0(IPHAS), CMU, XKAPPA,
     &        USTAR2, XKENT, XEENT )
C
          FLUCEL = 0.D0
          DO IELTSM = 1, NCESMP
            SMACEL(IELTSM,IPR(IPHAS)) = 30000.D0
            ITYPSM(IELTSM,IV(IPHAS)) = 1
            SMACEL(IELTSM,IV(IPHAS)) = VENT
            IF (ITYTUR(IPHAS).EQ.2) THEN
              ITYPSM(IELTSM,IK(IPHAS)) = 1
              SMACEL(IELTSM,IK(IPHAS)) = XKENT
              ITYPSM(IELTSM,IEP(IPHAS)) = 1
              SMACEL(IELTSM,IEP(IPHAS)) = XEENT
            ELSE IF (ITYTUR(IPHAS).EQ.3) THEN
              ITYPSM(IELTSM,IR11(IPHAS)) = 1
              ITYPSM(IELTSM,IR12(IPHAS)) = 1
              ITYPSM(IELTSM,IR13(IPHAS)) = 1
              ITYPSM(IELTSM,IR22(IPHAS)) = 1
              ITYPSM(IELTSM,IR23(IPHAS)) = 1
              ITYPSM(IELTSM,IR33(IPHAS)) = 1
              SMACEL(IELTSM,IR11(IPHAS)) = 2.D0/3.D0*XKENT
              SMACEL(IELTSM,IR12(IPHAS)) = 0.D0
              SMACEL(IELTSM,IR13(IPHAS)) = 0.D0
              SMACEL(IELTSM,IR22(IPHAS)) = 2.D0/3.D0*XKENT
              SMACEL(IELTSM,IR23(IPHAS)) = 0.D0
              SMACEL(IELTSM,IR33(IPHAS)) = 2.D0/3.D0*XKENT
              ITYPSM(IELTSM,IEP(IPHAS)) = 1
              SMACEL(IELTSM,IEP(IPHAS)) = XEENT
            ELSE IF (ITURB(IPHAS).EQ.50) THEN
              ITYPSM(IELTSM,IK(IPHAS)) = 1
              SMACEL(IELTSM,IK(IPHAS)) = XKENT
              ITYPSM(IELTSM,IEP(IPHAS)) = 1
              SMACEL(IELTSM,IEP(IPHAS)) = XEENT
              ITYPSM(IELTSM,IPHI(IPHAS)) = 1
              SMACEL(IELTSM,IPHI(IPHAS)) = 2.D0/3.D0
C     Il n'y a pas de terme source de masse dans l'equation de f_barre
            ELSE IF (ITURB(IPHAS).EQ.60) THEN
              ITYPSM(IELTSM,IK(IPHAS)) = 1
              SMACEL(IELTSM,IK(IPHAS)) = XKENT
              ITYPSM(IELTSM,IOMG(IPHAS))= 1
              SMACEL(IELTSM,IOMG(IPHAS))= XEENT/CMU/XKENT
            ENDIF
            IF(NSCAL.GT.0) THEN
              DO II = 1, NSCAL
                IF(IPHSCA(II).EQ.IPHAS) THEN
                  ITYPSM(IELTSM,ISCA(II)) = 1
                  SMACEL(IELTSM,ISCA(II)) = 1.D0
                ENDIF
              ENDDO
            ENDIF
            FLUCEL = FLUCEL+
     &                VOLUME(ICETSM(IELTSM))*SMACEL(IELTSM,IPR(IPHAS))
          ENDDO
C
          IF (IRANGP.GE.0) THEN
            CALL PARSOM (FLUCEL)
          ENDIF
C
          IF (IWARNI(IPR(IPHAS)).GE.1) THEN
            WRITE(NFECRA,1000) IPHAS, FLUCEL
          ENDIF
C
        ENDIF
C
C-----------------------------------------------------------------------
C
C Exemple 2 : simulation d'un soutirage (par une pompe par exemple)
C               total de 80 000 kg/s
C             On suppose que l'on souhaite repartir uniformement ce
C               puits de masse sur les NCESMP cellules selectionnees
C               plus haut.
C
C
C     Le test sur IUTILE permet de ne pas passer dans l'exemple si
C       on oublie de l'eliminer, lors de la mise en place d'un calcul
C       reel.
C
        IUTILE = 0
        IF(IUTILE.EQ.1) THEN
C
          IF(IPHAS.EQ.1) THEN
C
C     Calcul du volume total de la zone ou est impose le terme source
C       (le cas des calculs paralleles est prevu avec parsom)
            VTOT = 0.D0
            DO IELTSM = 1, NCESMP
              VTOT = VTOT + VOLUME(ICETSM(IELTSM))
            ENDDO
            IF (IRANGP.GE.0) THEN
              CALL PARSOM (VTOT)
            ENDIF
C
C     Le puits de masse est GAMMA = -80000/VTOT en kg/(m3 s)
C       (quel que soit le nombre de cellules NCESMP)
C     On l'impose ci-dessous (avec un test au cas ou VTOT=0)
C     On calcule au passage le terme puits total pour verification.
C
            IF (VTOT.GT.0.D0) THEN
              GAMMA = -80000.D0/VTOT
            ELSE
              WRITE(NFECRA,9000) IPHAS, VTOT
              CALL CSEXIT (1)
            ENDIF
C
            FLUCEL = 0.D0
            DO IELTSM = 1, NCESMP
              SMACEL(IELTSM,IPR(IPHAS)) = GAMMA
              FLUCEL = FLUCEL+
     &                VOLUME(ICETSM(IELTSM))*SMACEL(IELTSM,IPR(IPHAS))
            ENDDO
C
            IF (IRANGP.GE.0) THEN
              CALL PARSOM (FLUCEL)
            ENDIF
C
            IF (IWARNI(IPR(IPHAS)).GE.1) THEN
              WRITE(NFECRA,2000) IPHAS, FLUCEL, VTOT
            ENDIF
C
          ENDIF
C
        ENDIF
C
C-----------------------------------------------------------------------
C
      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(/,'PHASE ',I3,
     &         ' : FLUX DE MASSE GENERE DANS LE DOMAINE : ',E14.5,/)
C
 2000 FORMAT(/,'PHASE ',I3,
     &         ' : FLUX DE MASSE GENERE DANS LE DOMAINE : ',E14.5,/,
     &'                           REPARTI SUR UN VOLUME : ',E14.5)
C
 9000 FORMAT(/,'PHASE ',I3,
     &         ' : ERREUR DANS USTSMA ',/,
     &'   LE VOLUME DE LA ZONE AVEC PUITS DE MASSE VAUT = ',E14.5,/)
C
      RETURN
C
      END
c@z
