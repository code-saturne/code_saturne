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
                        SUBROUTINE USCFPV
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC ROUTINE UTILISATEUR : REMPLISSAGE DES VARIABLES PHYSIQUES
CFONC     POUR LA PHYSIQUE PARTICULIERE : COMPRESSIBLE SANS CHOC
CFONC     PENDANT DE USPHYV.F
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
CFONC      automatiquement fait dans uscfth en fonction de la thermo
CFONC      choisie
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
CFONC      - de la viscosite moleculaire                VISCL  kg/(m s)
CFONC      - de la chaleur specifique associee          CP     J/(kg degres)
CFONC      - de la conductivite thermique associee      LAMBDA W/(m degres)
CFONC      - des "diffusivites" associees aux scalaires VISCLS kg/(m s)
CFONC
CFONC  La masse volumique ne doit pas être renseignée : en compressible,
CFONC    c'est une variable résolue, que l'on peut initialiser si
CFONC    nécessaire dans USCFXI (RTP).
CFONC
CFONC  On dispose des types de faces de bord au pas de temps
CFONC    precedent (sauf au premier pas de temps, ou les tableaux
CFONC    ITYPFB et ITRIFB n'ont pas ete renseignes)
CFONC
CFONC L'identification des cellules concernees peut s'appuyer sur la
CFONC commande GETCEL.
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
CARGU ! W1...3(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
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
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML), MAXELT, LSTELT(MAXELT)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
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
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IVART, IEL, IPHAS
      INTEGER          IPCVIS, IPCVSV, IPCCP
      INTEGER          IPCVSL, ITH, ISCAL, II, ICCFTH, IMODIF
      DOUBLE PRECISION VARAM, VARBM, VARCM, VARDM
      DOUBLE PRECISION VARAL, VARBL, VARCL, VARDL
      DOUBLE PRECISION VARAC, VARBC
      DOUBLE PRECISION XRTP
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C
      IF(1.EQ.1) THEN
        IUSCFP = 0
        RETURN
      ENDIF
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C 0. INITIALISATIONS A CONSERVER
C=======================================================================
C
C --- Initialisation memoire
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C
C
C   LES EXEMPLES FANTAISISTES SUIVANTS SONT A ADAPTER PAR L'UTILISATEUR
C   ====================================================================
C
C                ILS SONT TOUS SUIVIS DE "CALL CSEXIT (1)" PAR SECURITE
C                                (stop a enlever a l'utilisation)
C
C
C   Il est conseille de ne garder dans ce sous programme que
C     le strict necessaire.
C
C
C EXEMPLE 1 : VISCOSITE VARIABLE EN FONCTION DE LA TEMPERATURE
C
C EXEMPLE 2 : VISCOSITE EN VOLUME VARIABLE EN FONCTION DE LA TEMPERATURE
C
C EXEMPLE 3 : CHALEUR SPECIFIQUE VARIABLE EN FONCTION DE LA TEMPERATURE
C
C EXEMPLE 4 : CONDUCTIVITE THERMIQUE VARIABLE
C                  EN FONCTION DE LA TEMPERATURE
C
C EXEMPLE 5 : DIFFUSIVITE VARIABLE EN FONCTION DE LA TEMPERATURE
C                  POUR LES SCALAIRES
C=======================================================================
C
C
C
C
C
C
C
C
C=======================================================================
C  EXEMPLE 1 : VISCOSITE VARIABLE EN FONCTION DE LA TEMPERATURE
C ===========
C    Ci dessous on donne pour toutes les phases la meme loi pour
C       la viscosite
C    Les valeurs de cette propriete doivent etre fournies au centre des
C       cellules.
C  ===================================================================
C
C --- Boucle sur les phases : debut
      DO IPHAS = 1, NPHAS
C
C
C   Positions des variables, coefficients
C   -------------------------------------
C
C --- Numero de variable température pour la phase courante iphas
C       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
C          IVART = ISCA(2)
C
        IVART = ISCA(ITEMPK(IPHAS))
C
C --- Rang de la viscosite dynamique moleculaire de la phase IPHAS
C     dans PROPCE, prop. physiques au centre des elements       : IPCVIS
C
        IPCVIS = IPPROC(IVISCL(IPHAS))
C
C --- Coefficients des lois choisis et imposes par l'utilisateur
C       Les valeurs donnees ici sont fictives
C
        VARAM = -3.4016D-9
        VARBM =  6.2332D-7
        VARCM = -4.5577D-5
        VARDM =  1.6935D-3
C
C   Viscosite moleculaire dynamique en kg/(m s) au centre des cellules
C   ------------------------------------------------------------------
C       loi              MU        =
C                              T  *( T  *( AM  * T +  BM  )+ CM  )+ DM
C       soit    PROPCE(IEL,IPCVIS) =
C     &                       XRTP*(XRTP*(VARAM*XRTP+VARBM)+VARCM)+VARDM
C
        DO IEL = 1, NCEL
          XRTP = RTP(IEL,IVART)
          PROPCE(IEL,IPCVIS) =
     &                       XRTP*(XRTP*(VARAM*XRTP+VARBM)+VARCM)+VARDM
        ENDDO
C
C
      ENDDO
C --- Boucle sur les phases : fin
C
C
C
C --- A enlever a l'utilisation
      IF(1.EQ.1) THEN
        WRITE(NFECRA,9000)
        CALL CSEXIT (1)
      ENDIF
C
C
C
C
C
C=======================================================================
C  EXEMPLE 2: VISCOSITE EN VOLUME VARIABLE EN FONCTION DE LA TEMPERATURE
C ==========
C    Ci dessous on donne pour toutes les phases la meme loi pour
C       la viscosite
C    Les valeurs de cette propriete doivent etre fournies au centre des
C       cellules.
C  ===================================================================
C
C --- Boucle sur les phases : debut
      DO IPHAS = 1, NPHAS
C
C
C   Positions des variables, coefficients
C   -------------------------------------
C
C --- Numero de variable température pour la phase courante iphas
C       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
C          IVART = ISCA(2)
C
        IVART = ISCA(ITEMPK(IPHAS))
C
C --- Rang de la viscosite dynamique moleculaire en volume de la phase IPHAS
C     dans PROPCE, prop. physiques au centre des elements       : IPCVSV
C
        IF(IVISCV(IPHAS).GT.0) THEN
          IPCVSV = IPPROC(IVISCV(IPHAS))
        ELSE
          IPCVSV = 0
        ENDIF
C
C --- Stop si non variable
C
        IF(IPCVSV.LE.0) THEN
          WRITE(NFECRA,2000) IPHAS, IPHAS, IVISCV(IPHAS)
          CALL CSEXIT (1)
        ENDIF
C
C
C --- Coefficients des lois choisis et imposes par l'utilisateur
C       Les valeurs donnees ici sont fictives
C
        VARAM = -3.4016D-9
        VARBM =  6.2332D-7
        VARCM = -4.5577D-5
        VARDM =  1.6935D-3
C
C   Viscosite en volume en kg/(m s) au centre des cellules
C   ------------------------------------------------------
C       loi            KAPPA       =
C                              T  *( T  *( AM  * T +  BM  )+ CM  )+ DM
C       soit    PROPCE(IEL,IPCVSV) =
C     &                       XRTP*(XRTP*(VARAM*XRTP+VARBM)+VARCM)+VARDM
C
        DO IEL = 1, NCEL
          XRTP = RTP(IEL,IVART)
          PROPCE(IEL,IPCVSV) =
     &                       XRTP*(XRTP*(VARAM*XRTP+VARBM)+VARCM)+VARDM
        ENDDO
C
      ENDDO
C --- Boucle sur les phases : fin
C
C
C
C --- A enlever a l'utilisation
      IF(1.EQ.1) THEN
        WRITE(NFECRA,9000)
        CALL CSEXIT (1)
      ENDIF
C
C
C
C
C
C=======================================================================
C  EXEMPLE 3 : CHALEUR SPECIFIQUE VARIABLE EN FONCTION DE LA TEMPERATURE
C ===========
C
C    Ci dessous on donne pour toutes les phases la meme loi pour
C       la chaleur specifique
C    Les valeurs de cette propriete doivent etre fournies au centre des
C       cellules.
C
C***********************************************************************
C
C                             ATTENTION !
C
C     NE PAS RETIRER LA MISE A JOUR DE Cv A LA FIN DE L'EXEMPLE
C
C***********************************************************************
C
C=======================================================================
C
C --- Boucle sur les phases : debut
      DO IPHAS = 1, NPHAS
C
C
C   Positions des variables, coefficients
C   -------------------------------------
C
C --- Numero de variable température pour la phase courante iphas
C       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
C          IVART = ISCA(2)
C
        IVART = ISCA(ITEMPK(IPHAS))
C
C --- Rang de la chaleur specifique de la phase courante IPHAS
C     dans PROPCE, prop. physiques au centre des elements       : IPCCP
C
        IF(ICP(IPHAS).GT.0) THEN
          IPCCP  = IPPROC(ICP   (IPHAS))
        ELSE
          IPCCP  = 0
        ENDIF
C
C --- Stop si CP ou CV n'est pas variable
C
        IF(IPCCP.LE.0) THEN
          WRITE(NFECRA,1000) IPHAS, IPHAS, ICP(IPHAS)
          CALL CSEXIT (1)
        ENDIF
        IF(ICV(IPHAS).LE.0) THEN
          WRITE(NFECRA,1001) IPHAS, IPHAS, ICV(IPHAS)
          CALL CSEXIT (1)
        ENDIF
C
C
C --- Coefficients des lois choisis et imposes par l'utilisateur
C       Les valeurs donnees ici sont fictives
C
        VARAC = 0.00001D0
        VARBC = 1000.0D0
C
C
C
C   Chaleur specifique J/(kg degres) au centre des cellules
C   --------------------------------------------------------
C       loi              CP        =  AC  * T   +  BM
C       soit    PROPCE(IEL,IPCCP ) = VARAC*XRTP + VARBC
C
        DO IEL = 1, NCEL
          XRTP = RTP(IEL,IVART)
          PROPCE(IEL,IPCCP ) = VARAC*XRTP + VARBC
        ENDDO
C
C
C --- Mise a jour de Cv
C
        ICCFTH = 432
        IMODIF = 0
C
        CALL USCFTH
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   ICCFTH , IMODIF , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   PROPCE(1, IPPROC(ICV(IPHAS)) )    , W1     , W2     , W3     ,
C        ---------------------------------
     &   RDEVEL , RTUSER , RA     )
C
      ENDDO
C --- Boucle sur les phases : fin
C
C
C
C --- A enlever a l'utilisation
      IF(1.EQ.1) THEN
        WRITE(NFECRA,9000)
        CALL CSEXIT (1)
      ENDIF
C
C
C
C
C=======================================================================
C  EXEMPLE 4 : CONDUCTIVITE THERMIQUE VARIABLE EN FONCTION
C ===========    DE LA TEMPERATURE
C
C    Ci dessous on donne pour toutes les phases la meme loi pour lambda
C    Les valeurs de cette propriete doivent etre fournies au centre des
C       cellules.
C  ===================================================================
C
C --- Boucle sur les phases : debut
      DO IPHAS = 1, NPHAS
C
C
C   Positions des variables, coefficients
C   -------------------------------------
C
C --- Numero de variable température pour la phase courante iphas
C       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
C          IVART = ISCA(2)
C
        IVART = ISCA(ITEMPK(IPHAS))
C
C --- Rang de Lambda de la temperature de phase courante IPHAS
C     dans PROPCE, prop. physiques au centre des elements       : IPCVSL
C
        IF(IVISLS(ITEMPK(IPHAS)).GT.0) THEN
          IPCVSL = IPPROC(IVISLS(ITEMPK(IPHAS)))
        ELSE
          IPCVSL = 0
        ENDIF
C
C --- Stop si Lambda n'est pas variable
C
        IF(IPCVSL.LE.0) THEN
          WRITE(NFECRA,1010)
     &      ITEMPK(IPHAS), ITEMPK(IPHAS), IVISLS(ITEMPK(IPHAS))
          CALL CSEXIT (1)
        ENDIF
C
C --- Coefficients des lois choisis et imposes par l'utilisateur
C       Les valeurs donnees ici sont fictives
C
        VARAL = -3.3283D-7
        VARBL =  3.6021D-5
        VARCL =  1.2527D-4
        VARDL =  0.58923D0
C
C
C
C   Lambda en W/(m K) au centre des cellules
C   ----------------------------------------
C       loi       Lambda        =
C               {  T  *( T  *( AL  * T +  BL  )+ CL  )+ DL }
C       soit    PROPCE(IEL,IPCVSL) =
C     &         (XRTP*(XRTP*(VARAL*XRTP+VARBL)+VARCL)+VARDL)
C
        DO IEL = 1, NCEL
          XRTP = RTP(IEL,IVART)
          PROPCE(IEL,IPCVSL) =
     &         (XRTP*(XRTP*(VARAL*XRTP+VARBL)+VARCL)+VARDL)
        ENDDO
C
C
      ENDDO
C --- Boucle sur les phases : fin
C
C
C --- A enlever a l'utilisation
      IF(1.EQ.1) THEN
        WRITE(NFECRA,9000)
        CALL CSEXIT (1)
      ENDIF
C
C
C
C
C=======================================================================
C  EXEMPLE 5 : DIFFUSIVITE VARIABLE EN FONCTION DE LA TEMPERATURE
C ===========      POUR LES SCALAIRES UTILISATEURS
C     A l'exclusion de
C        temperature, enthalpie (traites plus haut)
C        variances de fluctuations (propriete egale a celle du
C                                                      scalaire associe)
C
C    Ci dessous on donne pour tous les scalaires (aux exclusions
C      ci-dessus pres) la meme loi pour la diffusivite
C    Les valeurs de cette propriete doivent etre fournies au centre des
C       cellules.
C  ===================================================================
C
C --- Boucle sur les phases : debut
      DO II = 1, NSCAUS
C
C --- Numero du scalaire utilisateur II dans la liste de tous les scalaires
        ISCAL = II
C
C
C --- S'il s'agit de la temperature,
C                                   son cas a deja ete traite plus haut
        ITH = 0
        DO IPHAS = 1, NPHAS
          IF (ISCAL.EQ.ITEMPK(IPHAS)) ITH = 1
        ENDDO
C
C --- Si la variable est une fluctuation, sa diffusivite est
C       la meme que celle du scalaire auquel elle est rattachee :
C       il n'y a donc rien a faire ici : on passe directement
C       a la variable suivante sans renseigner PROPCE(IEL,IPCVSL).
C
        IF (ITH.EQ.0.AND.ISCAVR(ISCAL).LE.0) THEN
C --- On ne traite ici que les variables non thermiques
C                                   et qui ne sont pas des fluctuations
C
C
C   Positions des variables, coefficients
C   -------------------------------------
C
C --- Numero de variable température pour la phase courante iphas
C       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
C          IVART = ISCA(2)
C
        IVART = ISCA(ITEMPK(IPHAS))
C
C
C --- Rang de Lambda du scalaire
C     dans PROPCE, prop. physiques au centre des elements       : IPCVSL
C
          IF(IVISLS(ISCAL).GT.0) THEN
            IPCVSL = IPPROC(IVISLS(ISCAL))
          ELSE
            IPCVSL = 0
          ENDIF
C
C --- Stop si Lambda n'est pas variable
C
          IF(IPCVSL.LE.0) THEN
            WRITE(NFECRA,1010) ISCAL, ISCAL, IVISLS(ISCAL)
            CALL CSEXIT (1)
          ENDIF
C
C --- Coefficients des lois choisis et imposes par l'utilisateur
C       Les valeurs donnees ici sont fictives
C
          VARAL = -3.3283D-7
          VARBL =  3.6021D-5
          VARCL =  1.2527D-4
          VARDL =  0.58923D0
C
C
C   Lambda en kg/(m s) au centre des cellules
C   ------------------------------------------
C       loi       Lambda           =
C                  T  *( T  *( AL  * T +  BL  )+ CL  )+ DL
C       soit    PROPCE(IEL,IPCVSL) =
C     &          XRTP*(XRTP*(VARAL*XRTP+VARBL)+VARCL)+VARDL
C
C
          DO IEL = 1, NCEL
            XRTP = RTP(IEL,IVART)
            PROPCE(IEL,IPCVSL) =
     &           (XRTP*(XRTP*(VARAL*XRTP+VARBL)+VARCL)+VARDL)
          ENDDO
C
C
        ENDIF
C --- Tests sur ITH et ISCAVR : fin
C
      ENDDO
C --- Boucle sur les scalaires : fin
C
C
C --- A enlever a l'utilisation
      IF(1.EQ.1) THEN
        WRITE(NFECRA,9000)
        CALL CSEXIT (1)
      ENDIF
C
C
C
C
C
C=======================================================================

C=======================================================================
C FORMATS
C----
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,
     &'@    *********                                               ',/,
     &'@    DONNEES DE CALCUL INCOHERENTES                          ',/,
     &'@                                                            ',/,
     &'@    Pour la phase ',I10                                      ,/,
     &'@      usini1 indique que la capacite calorifique a pression ',/,
     &'@      constante est uniforme ICP(',I10   ,') = ',I10   ,'   ',/,
     &'@      alors que l''on cherche a definir une capacite        ',/,
     &'@      calorifique variable dans uscfpv.                     ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Modifier usini1 ou uscfpv.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,
     &'@    *********                                               ',/,
     &'@    DONNEES DE CALCUL INCOHERENTES                          ',/,
     &'@                                                            ',/,
     &'@    Pour la phase ',I10                                      ,/,
     &'@      uscfth indique que la capacite calorifique a volume   ',/,
     &'@      constant  est uniforme ICV(',I10   ,') = ',I10   ,'   ',/,
     &'@      alors que l''on cherche a definir une capacite        ',/,
     &'@      calorifique variable dans uscfpv.                     ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Modifier usini1 ou uscfpv.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,
     &'@    *********                                               ',/,
     &'@    DONNEES DE CALCUL INCOHERENTES                          ',/,
     &'@                                                            ',/,
     &'@    Pour le scalaire ',I10                                   ,/,
     &'@      usini1 indique que la diffusivite est uniforme        ',/,
     &'@        IVISLS(',I10   ,') = ',I10   ,' alors que           ',/,
     &'@      uscfpv impose une diffusivite variable.               ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Modifier usini1 ou uscfpv.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,
     &'@    *********                                               ',/,
     &'@    DONNEES DE CALCUL INCOHERENTES                          ',/,
     &'@                                                            ',/,
     &'@    Pour la phase ',I10                                      ,/,
     &'@      uscfx2 indique que la viscosite en volume est uniforme',/,
     &'@        IVISCV(',I10   ,') = ',I10   ,' alors que           ',/,
     &'@      uscfpv impose une viscosite en volume variable.       ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Modifier uscfx2 ou uscfpv.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,
     &'@    *********                                               ',/,
     &'@    APPEL A csexit DANS LE SOUS PROGRAMME uscfpv            ',/,
     &'@                                                            ',/,
     &'@    Un appel a csexit (arret) a ete rencontre dans le sous  ',/,
     &'@      programme uscfpv. L''utilisateur est invite a verifier',/,
     &'@      que les exemples standard fournis par defaut ont bien ',/,
     &'@      ete elimines si besoin.                               ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Modifier uscfpv.                                        ',/,
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
