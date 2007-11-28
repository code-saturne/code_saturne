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
                        SUBROUTINE USVPST
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , IPART  ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NVLSTA ,
     &   NCELPS , NFACPS , NFBRPS ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ITYPPS , IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   LSTCEL , LSTFAC , LSTFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , STATIS ,
     &   TRACEL , TRAFAC , TRAFBR , RDEVEL , RTUSER , RA     )
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  ROUTINE UTILISATEUR POUR LA SORTIE DE VARIABLES SUR UN MAILLAGE
CFONC    DE POST TRAITEMENT DEJA DEFINI
CFONC
CFONC  PAR DEFAUT, DEUX MAILLAGES SONT DEFINIS AUTOMATIQUEMENT :
CFONC       - LE MAILLAGE VOLUMIQUE (IPART=-1) SELON LA VALEUR DE ICHRVL
CFONC       - LE MAILLAGE DE BORD   (IPART=-2) SELON LA VALEUR DE ICHRBO
CFONC  DES MAILLAGES SUPPLEMENTAIRES (CELLULES OU FACES INTERNES ET
CFONC    DE BORD) PEUVENT ETRE DEFINIS ET PARAMETRES A TRAVERS
CFONC    usdpst.F ET usmpst.F.
CFONC
CETTE  CETTE ROUTINE EST APPELEE UNE FOIS PAR MAILLAGE POST
CFONC    ET PAR PAS DE TEMPS AUQUEL CE MAILLAGE EST ACTIF
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
CARGU ! IPART        ! E  !  -> ! NUMERO DU MAILLAGE POST              !
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
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VARIABLES STAT. LAGRANGIEN !
CARGU ! NCELPS       ! E  !  -> ! NOMBRE DE CELLULES DU MAILLAGE POST  !
CARGU ! NFACPS       ! E  !  -> ! NOMBRE DE FACES INTERIEUR POST       !
CARGU ! NFBRPS       ! E  !  -> ! NOMBRE DE FACES DE BORD POST         !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! ITYPPS(3)    ! TE !  -> ! INDICATEUR DE PRESENCE (0 OU 1) DE   !
CARGU !              !    !     ! CELLULES (1), FACES (2), OU FACES DE !
CARGU !              !    !     ! DE BORD (3) DANS LE MAILLAGE POST    !
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
CARGU ! LSTCEL(NCELPS! TE !  -> ! LISTE DES CELLULES DU MAILLAGE POST  !
CARGU ! LSTFAC(NFACPS! TE !  -> ! LISTE DES FACES INTERIEURES POST     !
CARGU ! LSTFBR(NFBRPS! TE !  -> ! LISTE DES FACES DE BORD POST         !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TE !  -  ! MACRO TABLEAU ENTIER                 !
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
CARGU ! (NCELET)     !    !     !                                      !
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
CARGU ! STATIS       ! TR !  -> ! STATISTIQUES (LAGRANGIEN)            !
CARGU !NCELET,NVLSTA)!    !     !                                      !
CARGU ! TRACEL(*)    ! TR ! <-> ! TAB REEL VALEURS CELLULES POST       !
CARGU ! TRAFAC(*)    ! TR ! <-> ! TAB REEL VALEURS FACES INT. POST     !
CARGU ! TRAFBR(*)    ! TR ! <-> ! TAB REEL VALEURS FACES BORD POST     !
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
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "numvar.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          IPART
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS  , NVLSTA
      INTEGER          NCELPS , NFACPS , NFBRPS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          ITYPPS(3)
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          LSTCEL(NCELPS), LSTFAC(NFACPS), LSTFBR(NFBRPS)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTPA(NCELET,*), RTP(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA)
      DOUBLE PRECISION TRACEL(NCELPS*3)
      DOUBLE PRECISION TRAFAC(NFACPS*3), TRAFBR(NFBRPS*3)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*32     NAMEVR
C
      INTEGER          NTINDP
      INTEGER          IEL   , IFAC  , ILOC  , IPHAS, IVAR , ICLT
      INTEGER          IDIMT , II   , JJ
      INTEGER          IDIMTE, ITENSO, IENTLA, IVARPR
      INTEGER          IMOM1, IMOM2, IPCMO1, IPCMO2, IDTCM
      DOUBLE PRECISION POND
      DOUBLE PRECISION RBID(1)
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C
C***********************************************************************
C
C
C
C=======================================================================
C     1. TRAITEMENT DES VARIABLES A SORTIR
C         A RENSEIGNER PAR L'UTILISATEUR aux endroits indiques
C=======================================================================
C
C
C     Un maillage de posttraitement est une "part"
C       (au sens EnSight ; les équivalents MED et CGNS sont le maillage
C       et la base respectivement)
C     L'utilisateur aura defini ses maillages de posttraitement dans
C       usdpst (NBPART maillages de posttraitement)
C
C
C     La routine est appelee une fois pour chaque maillage IPART.
C     Pour chaque maillage et pour chacune des variables que l'on
C       souhaite posttraiter, on doit definir certains parametres et les
C       passer a la routine PSTEVA qui se charge de l'ecriture effective.
C     Ces parametres sont :
C       NAMEVR : nom de la variable
C       IDIMT  : dimension de la variable
C       IENTLA : dans le cas ou IDIMT est >1, IENTLA permet de specifier
C                si le tableau contenant la variable est range de maniere
C                entrelacee X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,... (IENTLA=1)
C                ou non entrelancee X1,X2,X3,...,Y1,Y2,Y3,...,Z1,Z2,Z3,...
C                                                               (IENTLA=0)
C       IVARPR : specifie si le tableau contenant la variable traitee est
C                defini sur le maillage "parent"/
C                En effet, meme si le maillage IPART considere contient
C                le meme nombre d'elements que le maillage complet "parent"
C                (NCELPS=NCEL), l'ordre de numerotation des elements n'est pas
C                forcement le meme. Le tableau TRACEL passe en argument de
C                PSTEVA est construit selon la numerotation du maillage IPART.
C                Pour posttraiter une variable contenue dans un tableau RTUSER
C                par exemple, il faut donc d'abord le reordonner dans le
C                tableau TRACEL :
C                  DO ILOC = 1, NCELPS
C                    IEL = LSTCEL(ILOC)
C                    TRACEL(ILOC) = RTUSER(IEL)
C                  ENDDO
C                Une alternative est cependant offerte, pour eviter des copies
C                inutiles. Si NCELPS=NCEL, on peut directement passer RTUSER
C                en argument de PSTEVA, en specifiant IVARPR=1, pour avertir
C                le code que la numerotation est celle du maillage "parent".
C                L'exemple ci-dessus concerne les cellules, mais la
C                problematique est la meme avec les faces internes ou faces de
C                bord.
C
C
C     Remarque : attention aux longueurs des noms de variables.
C
C                On autorise ici jusqu'à 32 caracteres , mais selon le
C                format utilise, les noms peuvent etre tronques :
C
C                  - a 19 caracteres pour une sortie EnSight
C                  - a 32 caracteres pour ue sortie MED 2.2
C
C                La longueur du nom n'est pas limitee en interne, et en
C                cas de deux variables aux noms tronques ne differant
C                qu'apres le 19ieme caractere, les lignes correspondantes
C                apparaitront normalement dans le fichier texte ".case"
C                EnSight, avec un meme champ de description ; il suffit
C                alors de renommer un de ces champs dans ce fichier
C                texte pour corriger le probleme.
C
C                Les caracteres blancs en debut ou fin de chaine sont
C                supprimes automatiquement. Selon le format utilise,
C                les caracteres interdits (sous EnSight, les caracteres
C                (  ) ] [ + - @ ! # * ^ $ / ainsi que les blancs et les
C                tabulations seront remplaces par le caractere _.
C                Ceci ne necessite aucune intervention utilisateur
C                (i.e. inutile d'utiliser ici les anciennes fonctions
C                Fortran VERLON et UNDSCR).
C
C
C     Exemples :
C               pour le maillage post 1, on sort
C                 la temperature interpolee a la face
C                   (et 0 sur les eventuelles faces interieures)
C
C               pour le maillage post 2, on sort
C                 la temperature
C
C
      IF (IPART.EQ.1 ) THEN
C
C
C       Initialisation
C         pas d'intervention utilisateur requise
        DO II = 1, 32
          NAMEVR (II:II) = ' '
        ENDDO
C
C       Nom de la variable
C         a renseigner par l'utilisateur
        NAMEVR = 'Temperature interpolee'
C
C       Dimension de la variable (3 = vecteur, 1=scalaire)
C         a renseigner par l'utilisateur
        IDIMT = 1
C
C       Valeurs entrelacées
        IENTLA = 0
C
C
C       Calcul des valeurs de la variable sur les faces internes
C         Pour simplifier l'exemple, on se contente ici d'une simple
C           interpolation lineaire.
C         Dans les calculs paralleles, si l'on utilise les voisins,
C           il faut faire un echange au prealable, comme d'habitude.
C         Dans les calculs avec periodicite, il faut egalement le faire
C         Pour les calculs ou periodicite et parallelisme coexistent,
C           l'appel a ces routines doit etre fait dans l'ordre
C           PARCOM puis PERCOM
C
C         a renseigner par l'utilisateur
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (RTP(1,ISCA(1)))
        ENDIF
C
        IF(IPERIO.EQ.1) THEN
          IVAR = ISCA(1)
          IDIMTE = 0
          ITENSO = 0
          CALL PERCOM
C         ===========
     &         ( IDIMTE , ITENSO ,
     &           RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &           RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &           RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR))
        ENDIF
C
        DO ILOC = 1, NFACPS
C
          IFAC = LSTFAC(ILOC)
          II = IFACEL(1, IFAC)
          JJ = IFACEL(2, IFAC)
          POND = RA(IPOND-1+IFAC)
C
          TRAFAC(ILOC)
     &         =         POND  * RTP(II, ISCA(1))
     &         + (1.D0 - POND) * RTP(JJ, ISCA(1))
C
        ENDDO
C
C
C       Ecriture effective des valeurs calculees
C
        IVARPR = 0
C
        CALL PSTEVA(IPART , NAMEVR, IDIMT, IENTLA, IVARPR,
C       ===========
     &              NTCABS, TTCABS, TRACEL, TRAFAC, TRAFBR)
C
C
C
C
      ELSE IF  (IPART.EQ.2) THEN
C
C
C       1.4.1 TRAITEMENT DE LA VITESSE
C       ------------------------------
C
C       Initialisation
C         pas d'intervention utilisateur requise
        DO II = 1, 32
          NAMEVR (II:II) = ' '
        ENDDO
C
C       Nom de la variable
C         a renseigner par l'utilisateur
        NAMEVR = 'Temperature'
C
C       Dimension de la variable (3 = vecteur, 1=scalaire)
C         a renseigner par l'utilisateur
        IDIMT = 1
C
C       Valeurs non entrelacées
        IENTLA = 0
C
        DO ILOC = 1, NCELPS
          IEL = LSTCEL(ILOC)
          TRACEL(ILOC) = RTP(IEL,ISCA(1))
        ENDDO
C
        IVARPR = 0
C
        CALL PSTEVA(IPART , NAMEVR, IDIMT, IENTLA, IVARPR,
C       ===========
     &              NTCABS, TTCABS, TRACEL, TRAFAC, TRAFBR)
C
      ENDIF
C     Fin du test sur le numero de maillage post.
C
C
      RETURN
C
      END
c@z
