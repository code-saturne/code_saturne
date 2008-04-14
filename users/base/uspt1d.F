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
                        SUBROUTINE USPT1D
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NFPT1D , IPHAS  , IAPPEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IFPT1D , NPPT1D , ICLT1D ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO ,
     &   XYZNOD , VOLUME ,
     &   TPPT1D , RGPT1D , EPPT1D ,
     &   TEPT1D , HEPT1D , FEPT1D ,
     &   XLMT1D , RCPT1D , DTPT1D ,
     &   DT     , RTPA   ,
     &   PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RDEVEL , RTUSER , RA     )
C
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC      ENTREE DES DONNEES DU MODULE THERMIQUE EN 1D PAROI
CFONC
CFONC IAPPEL = 1 (un seul appel a l'initialisation) :
CFONC              CALCUL DU NOMBRE DE CELLULES OU L'ON IMPOSE UNE PAROI
CFONC
CFONC IAPPEL = 2 (un seul appel a l'initialisation) :
CFONC              REPERAGE DES CELLULES OU L'ON IMPOSE UNE PAROI
CFONC              DONNEES RELATIVE AU MAILLAGE
CFONC
CFONC IAPPEL = 3 (appel a chaque pas de temps) :
CFONC              VALEUR DES COEFFICIENTS PHYSIQUE DU CALCUL
CFONC              TYPE DE CONDITION LIMITE EN PAROI EXTERIEURE
CFONC              ICLT1D = 1 -> TEMPERATURE
CFONC              ICLT1D = 3 -> FLUX
CFONC              INITIALISATION DE LA TEMPERATURE
CFONC
CFONC
CFONC IDENTIFICATION DES FACES DE BORD
CFONC ================================
CFONC L'identification des faces de bord concernees se fait grace
CFONC a la commande GETFBR.
CFONC
CFONC   GETFBR(CHAINE,NLELT,LSTELT) :
CFONC   - CHAINE est une chaine de caractere fournie par l'utilisateur
CFONC     qui donne les criteres de selection
CFONC   - NLTELT est renvoye par la commande. C'est un entier qui
CFONC     correspond au nombre de faces de bord trouveees repondant au
CFONC     critere
CFONC   - LSTELT est renvoye par la commande. C'est un tableau d'entiers
CFONC     de taille NLTELT donnant la liste des faces de bord trouvees
CFONC     repondant au critere.
CFONC
CFONC   CHAINE peut etre constitue de :
CFONC   - references de couleurs (ex. : 1, 8, 26, ...
CFONC   - references de groupes (ex. : entrees, groupe1, ...)
CFONC   - criteres geometriques (ex. X<0.1, Y>=0.25, ...)
CFONC   Ces criteres peuvent etre combines par des operateurs logiques
CFONC   (AND et OR) et des parentheses
CFONC   ex. : '1 AND (groupe2 OR groupe3) AND Y<1' permettra de recuperer
CFONC   les faces de bord de couleur 1, appartenant aux groupes 'groupe2'
CFONC   ou 'groupe3' et de coordonnee Y inferieure a 1.
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
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NFML         ! E  !  -> ! NOMBRE DE FAMILLES D ENTITES         !
CARGU ! NPRFML       ! E  !  -> ! NOMBRE DE PROPRIETESE DES FAMILLES   !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NFPT1D       ! E  !  -> ! NOMBRE DE FACES AVEC MODULE THERM 1D !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IAPPEL       ! E  !  -> ! INDIQUE LES DONNES A RENVOYER        !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMFBR       ! TE !  -> ! NUMERO DE FAMILLE D'UNE FACE DE BORD !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IPRFML       ! TE !  -> ! PROPRIETES D'UNE FAMILLE             !
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! MAXELT       !  E !  -> ! NB MAX D'ELEMENTS (CELL,FAC,FBR)     !
CARGU ! LSTELT(MAXELT) TE !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! IFPT1D       ! TE !  -> ! NUMERO DE LA FACE EN TRAITEMENT      !
CARGU !              !    !     ! THERMIQUE EN PAROI                   !
CARGU ! NPPT1D       ! TE !  -> ! NOMBRE DE POINTS DE DISCRETISATION   !
CARGU !              !    !     ! DANS LA PAROI                        !
CARGU ! ICLT1D       ! TE !  -> ! TYPE DE CONDITION LIMITE             !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! EPPT1D       ! TR !  -> ! EPAISSEUR DE LA PAROI                !
CARGU ! RGPT1D       ! TR !  -> ! RAISON DU MAILLAGE                   !
CARGU ! TPPT1D       ! TR !  -> ! TEMPERATURE DE PAROI                 !
CARGU ! TEPT1D       ! TR !  -> ! TEMPERATURE EXTERIEURE               !
CARGU ! HEPT1D       ! TR !  -> ! COEFFICIENT D'ECHANGE EXTERIEUR      !
CARGU ! FEPT1D       ! TR !  -> ! FLUX EXTERIEUR                       !
CARGU ! XLMT1D       ! TR !  -> ! CONDUCTIVITE THERMIQUE DE LA PAROI   !
CARGU ! RCPT1D       ! TR !  -> ! ROCP DE LA PAROI                     !
CARGU ! DTPT1D       ! TR !  -> ! PAS DE TEMPS DE LA PAROI             !
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
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
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
C
C***********************************************************************
C
C Arguments
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS  , NFPT1D
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IPHAS  , IAPPEL
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          MAXELT, LSTELT(MAXELT)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IFPT1D(NFPT1D), NPPT1D(NFPT1D), ICLT1D(NFPT1D)
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
      DOUBLE PRECISION EPPT1D(NFPT1D) , RGPT1D(NFPT1D) , TPPT1D(NFPT1D)
      DOUBLE PRECISION TEPT1D(NFPT1D) , HEPT1D(NFPT1D) , FEPT1D(NFPT1D)
      DOUBLE PRECISION XLMT1D(NFPT1D) , RCPT1D(NFPT1D) , DTPT1D(NFPT1D)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C Variables locales
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IFBT1D , II , IFAC , IFML , ICOUL
      INTEGER          ILELT, NLELT
C
C***********************************************************************
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
C --- Relecture d'un fichier suite :
C     ISUIT1 = 0 ------> pas de relecture (reinitialisation du maillage et
C                                           de la temperature dans la paroi)
C     ISUIT1 = 1 ------> relecture du fichier suite de module thermique 1D
C     ISUIT1 = ISUITE -> relecture si le calcul fluide est une suite
C     L'initialisation de ISUIT1 dans uspt1d est obligatoire.
C
      ISUIT1 = ISUITE
C
      IFBT1D = 0
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C
      IF(1.EQ.1) RETURN
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
      IF(IAPPEL.EQ.1.OR.IAPPEL.EQ.2) THEN
C
C
C --- Determination des faces avec module thermique 1D
C     NFPT1D     : nb total de faces avec module thermique 1D
C     IFPT1D(II) : numero de la IIeme face avec module thermque 1D
C
C     Remarque : lors de la relecture d'un fichier suite, NFPT1D et IFPT1D
C                sont compares aux valeurs issues du fichier suite. Une
C                concordance totale est necessaire pour continuer le calcul.
C                En ce qui concerne le test sur IFPT1D, il necessite qui
C                le tableau soit range dans un ordre croissant
C                   ( IFPT1D(JJ) > IFPT1D(II) si JJ > II ).
C                Si ce n'est pas possible, contacter l'equipe de developpement
C                pour desactiver le test.
C
        CALL GETFBR('3',NLELT,LSTELT)
C       ===========
C
        DO ILELT = 1, NLELT
C
          IFAC = LSTELT(ILELT)
C
          IFBT1D =IFBT1D + 1
          IF (IAPPEL.EQ.2) IFPT1D(IFBT1D) = IFAC
C
        ENDDO
C
      ENDIF
C
      IF (IAPPEL.EQ.1) THEN
         NFPT1D = IFBT1D
      ENDIF
C
C --- Remplissage des parametres du maillage et d'initialisation
C     (un seul passage en debut de calcul)
C     NPPT1D(II) : nombre de points de discretisation associes a la
C                                   IIeme face avec module thermique 1D
C     EPPT1D(II) : epaisseur de la paroi associee a la IIeme face avec
C                                                   module thermique 1D
C     RGPT1D(II) : raison geometrique du raffinement du maillage associe
C                              a la IIeme face avec module thermique 1D
C                       ( RGPT1D(II) > 1 => petites mailles cote fluide )
C     TPPT1D(II) : temperature d'initialisation de la paroi associee a la
C                                   IIeme face avec module thermique 1D
C
C     Remarque : lors de la relecture d'un fichier suite de module thermique
C                1D, TPPT1D n'est pas utilise. NFPT1D, EPPT1D et RGPT1D sont
C                compares aux valeurs issues du fichier suite. Une
C                correspondance exacte est necessaire pour continuer le calcul.
      IF (IAPPEL.EQ.2) THEN
         IF(IPHAS.EQ.1) THEN
            DO II = 1, NFPT1D
              IFAC = IFPT1D(II)
              NPPT1D(II) = 8
              EPPT1D(II) = 0.01144D0
              RGPT1D(II) = 1.D0
              TPPT1D(II) = 25.D0
            ENDDO
         ENDIF
      ENDIF
C
C
C --- Remplissage des conditions aux limites en paroi externe
C     ICLT1D(II) : type de condition a la limite
C                  = 1 : condition de dirichlet, avec coefficient d'echange
C                  = 3 : condition de flux
C     TEPT1D(II) : temperature exterieure
C     HEPT1D(II) : coefficient d'echange exterieur
C     FEPT1D(II) : flux applique a l'exterieur ( flux<0 = flux entrant)
C     XLMT1D(II) : coefficient de conductivite lambda de la paroi (W/m/°C)
C     RCPT1D(II) : coefficient rho*Cp de la paroi (J/m3/°C)
C     DTPT1D(II) : pas de temps de resolution de l'equation thermique dans
C                  la IIeme face de bord avec module thermique 1D (s)
      IF (IAPPEL.EQ.3) THEN
         IF(IPHAS.EQ.1) THEN
            DO II = 1, NFPT1D
               ICLT1D(II) = 1
C     parametres physiques
               IFAC = IFPT1D(II)
               IF (CDGFBO(2,IFAC).LE.0.025D0) THEN
                 ICLT1D(II) = 3
                 FEPT1D(II) = -1.D4
               ELSE
                 ICLT1D(II) = 3
                 FEPT1D(II) =  1.D4
               ENDIF
               XLMT1D(II) = 31.5D0
               RCPT1D(II) = 3.5D6
               DTPT1D(II) = 0.3D0
            ENDDO
         ENDIF
      ENDIF
C
      RETURN
C
      END
c@z

