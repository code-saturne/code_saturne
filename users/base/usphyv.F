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
                        SUBROUTINE USPHYV
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IBROM  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   ,
     &   PROPCE , PROPFA , PROPFB ,
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
CFONC ROUTINE UTILISATEUR : REMPLISSAGE DES VARIABLES PHYSIQUES
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
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
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
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR), IBROM(NPHMX)
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
      INTEGER          IVART, ICLVAR, IEL, IPHAS
      INTEGER          IPCROM, IPBROM, IPCVIS, IPCCP
      INTEGER          IPCVSL, ITH, ISCAL, II
      INTEGER          IUTILE
      DOUBLE PRECISION VARA, VARB, VARC, VARAM, VARBM, VARCM, VARDM
      DOUBLE PRECISION                   VARAL, VARBL, VARCL, VARDL
      DOUBLE PRECISION                   VARAC, VARBC
      DOUBLE PRECISION XRTP
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
C 0. INITIALISATIONS A CONSERVER
C=======================================================================
C
C --- Initialisation memoire
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
C
C
C
C=======================================================================
C
C
C   LES EXEMPLES FANTAISISTES SUIVANTS SONT A ADAPTER PAR L'UTILISATEUR
C   ====================================================================
C
C   Chaque exemple est encadre par un test sur IUTILE, par securite.
C   Mettre IUTILE a 1 pour activer l'exemple.
C
C   Il est conseille de ne garder dans ce sous programme que
C     le strict necessaire.
C
C
C
C  EXEMPLE 1 : MASSE VOLUMIQUE VARIABLE EN FONCTION DE LA TEMPERATURE
C  EXEMPLE 2 : VISCOSITE       VARIABLE EN FONCTION DE LA TEMPERATURE
C  EXEMPLE 3 : CHALEUR SPECIFIQUE VARIABLE EN FONCTION DE LA TEMPERATURE
C  EXEMPLE 4 : Lambda/CP  VARIABLE EN FONCTION DE LA TEMPERATURE
C                  POUR LA TEMPERATURE OU L'ENTHALPIE
C  EXEMPLE 5 : DIFFUSIVITE VARIABLE EN FONCTION DE LA TEMPERATURE
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
C
C
C=======================================================================
C  EXEMPLE 1 : MASSE VOLUMIQUE VARIABLE EN FONCTION DE LA TEMPERATURE
C ===========
C    Ci dessous on donne pour toutes les phases la meme loi pour
C       la masse volumique
C    Les valeurs de cette propriete doivent etre fournies au centre des
C       cellules (et, de facon optionnelle, aux faces de bord).
C  ===================================================================
C
C     Le test sur IUTILE permet de desactiver les instructions (qui
C       ne sont fournies qu'a titre d'exemple a adapter)
C
      IUTILE = 0
      IF(IUTILE.EQ.1) THEN
C
C --- Boucle sur les phases : debut
        DO IPHAS = 1, NPHAS
C
C
C   Positions des variables, coefficients
C   -------------------------------------
C
C --- Numero de variable thermique pour la phase courante iphas
C       (et de ses conditions limites)
C       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
C          IVART = ISCA(2)
C
          IF (ISCALT(IPHAS).GT.0) THEN
            IVART = ISCA(ISCALT(IPHAS))
          ELSE
            WRITE(NFECRA,9010) ISCALT(IPHAS)
            CALL CSEXIT (1)
          ENDIF
C
C --- Position des conditions limites de la variable IVART
C
          ICLVAR = ICLRTP(IVART,ICOEF)
C
C --- Rang de la masse volumique de la phase courante IPHAS
C     dans PROPCE, prop. physiques au centre des elements       : IPCROM
C     dans PROPFB, prop. physiques au centre des faces de bord  : IPBROM
C
          IPCROM = IPPROC(IROM(IPHAS))
          IPBROM = IPPROB(IROM(IPHAS))
C
C --- Coefficients des lois choisis et imposes par l'utilisateur
C       Les valeurs donnees ici sont fictives
C
          VARA  = -4.0668D-3
          VARB  = -5.0754D-2
          VARC  =  1000.9D0
C
C
C
C   Masse volumique au centre des cellules
C   ---------------------------------------
C       loi              RHO       =   T  * (  A *  T +  B ) +   C
C       soit    PROPCE(IEL,IPCROM) = XRTP * (VARA*XRTP+VARB) + VARC
C
          DO IEL = 1, NCEL
            XRTP = RTP(IEL,IVART)
            PROPCE(IEL,IPCROM) = XRTP * (VARA*XRTP+VARB) + VARC
          ENDDO
C
C
C
C   Masse volumique aux faces de bord
C   ----------------------------------
C
C       Par defaut, la valeur de rho au bord est la valeur prise
C         au centre des elements voisins. C'est l'approche conseillee.
C       Pour etre dans ce cas il suffit de ne rien faire :
C         ne pas prescrire de valeur pour PROPFB(IFAC,IPBROM) et
C         ne pas modifier IBROM(IPHAS)
C         ---------------
C
C       Pour les utilisateurs qui ne souhaiteraient pas suivre ce
C         conseil, on precise que la temperature au bord peut etre
C         fictive, simplement destinee a conserver un flux (c'est
C         en particulier le cas en paroi). La valeur de rho calculee
C         au bord avec en introduisant cette temperature fictive
C         dans une loi physique peut donc etre totalement fausse
C         (negative par exemple).
C
C       Si malgre tout on souhaite imposer la loi :
C                        RHO       =   T  * (  A *  T +  B ) +   C
C       soit   PROPFB(IFAC,IPBROM) = XRTP * (VARA*XRTP+VARB) + VARC
C       T etant la temperature prise au centre des faces de bord,
C
C         on peut utiliser les lignes de code suivantes (volontairement
C         desactivees, car a manier avec precaution) :
C
C         Noter bien que dans le cas ou l'on impose la masse volumique
C         au bord, il faut le faire sur TOUTES les faces de bord.
C                                       ======
C
C          IBROM(IPHAS) = 1
C          DO IFAC = 1, NFABOR
C            IEL = IFABOR(IFAC)
C            XRTP = COEFA(IFAC,ICLVAR)+RTP(IEL,IVART)*COEFB(IFAC,ICLVAR)
C            PROPFB(IFAC,IPBROM) = XRTP * (VARA*XRTP+VARB) + VARC
C          ENDDO
C
C         IFABOR(IFAC) est l'element en regard de la face de bord
C
C         Attention IBROM(IPHAS) = 1 est indispensable pour que la loi
C           soit  prise en compte.       -------------
C
C
C
        ENDDO
C --- Boucle sur les phases : fin
      ENDIF
C --- Test sur IUTILE : fin
C
C
C
C
C
C
C=======================================================================
C  EXEMPLE 2 : VISCOSITE       VARIABLE EN FONCTION DE LA TEMPERATURE
C ===========
C    Ci dessous on donne pour toutes les phases la meme loi pour
C       la viscosite
C    Les valeurs de cette propriete doivent etre fournies au centre des
C       cellules.
C  ===================================================================
C
C     Le test sur IUTILE permet de desactiver les instructions (qui
C       ne sont fournies qu'a titre d'exemple a adapter)
C
      IUTILE = 0
      IF(IUTILE.EQ.1) THEN
C
C --- Boucle sur les phases : debut
        DO IPHAS = 1, NPHAS
C
C
C   Positions des variables, coefficients
C   -------------------------------------
C
C --- Numero de variable thermique pour la phase courante iphas
C       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
C          IVART = ISCA(2)
C
          IF (ISCALT(IPHAS).GT.0) THEN
            IVART = ISCA(ISCALT(IPHAS))
          ELSE
            WRITE(NFECRA,9010) ISCALT(IPHAS)
            CALL CSEXIT (1)
          ENDIF
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
C
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
     &           XRTP*(XRTP*(VARAM*XRTP+VARBM)+VARCM)+VARDM
          ENDDO
C
C
        ENDDO
C --- Boucle sur les phases : fin
      ENDIF
C --- Test sur IUTILE : fin
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
C  ===================================================================
C
C     Le test sur IUTILE permet de desactiver les instructions (qui
C       ne sont fournies qu'a titre d'exemple a adapter)
C
      IUTILE = 0
      IF(IUTILE.EQ.1) THEN
C
C --- Boucle sur les phases : debut
        DO IPHAS = 1, NPHAS
C
C
C   Positions des variables, coefficients
C   -------------------------------------
C
C --- Numero de variable thermique pour la phase courante iphas
C       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
C          IVART = ISCA(2)
C
          IF (ISCALT(IPHAS).GT.0) THEN
            IVART = ISCA(ISCALT(IPHAS))
          ELSE
            WRITE(NFECRA,9010) ISCALT(IPHAS)
            CALL CSEXIT (1)
          ENDIF
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
C --- Stop si CP n'est pas variable
C
          IF(IPCCP.LE.0) THEN
            WRITE(NFECRA,1000) IPHAS, IPHAS, ICP(IPHAS)
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
        ENDDO
C --- Boucle sur les phases : fin
      ENDIF
C --- Test sur IUTILE : fin
C
C
C
C
C
C
C=======================================================================
C  EXEMPLE 4 : Lambda/CP  VARIABLE EN FONCTION DE LA TEMPERATURE
C ===========      POUR LA TEMPERATURE OU L'ENTHALPIE
C
C    Ci dessous on donne pour toutes les phases la meme loi pour
C       le rapport lambda/Cp
C    Les valeurs de cette propriete doivent etre fournies au centre des
C       cellules.
C  ===================================================================
C
C     Le test sur IUTILE permet de desactiver les instructions (qui
C       ne sont fournies qu'a titre d'exemple a adapter)
C
      IUTILE = 0
      IF(IUTILE.EQ.1) THEN
C
C --- Boucle sur les phases : debut
        DO IPHAS = 1, NPHAS
C
C
C   Positions des variables, coefficients
C   -------------------------------------
C
C --- Numero de variable thermique pour la phase courante iphas
C       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
C          IVART = ISCA(2)
C
          IF (ISCALT(IPHAS).GT.0) THEN
            IVART = ISCA(ISCALT(IPHAS))
          ELSE
            WRITE(NFECRA,9010) ISCALT(IPHAS)
            CALL CSEXIT (1)
          ENDIF
C
C --- Rang de Lambda/CP de la variable thermique de phase courante IPHAS
C     dans PROPCE, prop. physiques au centre des elements       : IPCVSL
C
          IF(IVISLS(ISCALT(IPHAS)).GT.0) THEN
            IPCVSL = IPPROC(IVISLS(ISCALT(IPHAS)))
          ELSE
            IPCVSL = 0
          ENDIF
C
C --- Stop si Lambda/CP n'est pas variable
C
          IF(IPCVSL.LE.0) THEN
            WRITE(NFECRA,1010)
     &           ISCALT(IPHAS), ISCALT(IPHAS), IVISLS(ISCALT(IPHAS))
            CALL CSEXIT (1)
          ENDIF
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
C   Lambda/Cp en kg/(m s) au centre des cellules
C   ------------------------------------------
C       loi       Lambda/CP        =
C               {  T  *( T  *( AL  * T +  BL  )+ CL  )+ DL } / Cp
C       soit    PROPCE(IEL,IPCVSL) =
C     &         (XRTP*(XRTP*(VARAL*XRTP+VARBL)+VARCL)+VARDL)/CP0(IPHAS)
C
C       On suppose Cp renseigne au prealable.
C
          IF(IPCCP.LE.0) THEN
C
C --- Si CP est uniforme, on utilise CP0(IPHAS)
            DO IEL = 1, NCEL
              XRTP = RTP(IEL,IVART)
              PROPCE(IEL,IPCVSL) =
     &             (XRTP*(XRTP*(VARAL*XRTP+VARBL)+VARCL)+VARDL)
     &             /CP0(IPHAS)
            ENDDO
C
          ELSE
C
C --- Si CP est non uniforme, on utilise PROPCE ci dessus
            DO IEL = 1, NCEL
              XRTP = RTP(IEL,IVART)
              PROPCE(IEL,IPCVSL) =
     &             (XRTP*(XRTP*(VARAL*XRTP+VARBL)+VARCL)+VARDL)
     &             /PROPCE(IEL,IPCCP)
            ENDDO
C
          ENDIF
C
C
        ENDDO
C --- Boucle sur les phases : fin
      ENDIF
C --- Test sur IUTILE : fin
C
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
C     Le test sur IUTILE permet de desactiver les instructions (qui
C       ne sont fournies qu'a titre d'exemple a adapter)
C
      IUTILE = 0
      IF(IUTILE.EQ.1) THEN
C
C --- Boucle sur les scalaires : debut
        DO II = 1, NSCAUS
C
C --- Numero du scalaire utilisateur II dans la liste de tous les scalaires
          ISCAL = II
C
C
C --- S'il s'agit d'une variable thermique,
C                                   son cas a deja ete traite plus haut
          ITH = 0
          DO IPHAS = 1, NPHAS
            IF (ISCAL.EQ.ISCALT(IPHAS)) ITH = 1
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
C --- Numero de variable thermique pour la phase courante iphas
C       (Pour utiliser le scalaire utilisateur 2 a la place, ecrire
C          IVART = ISCA(2)
C
            IPHAS = IPHSCA(ISCAL)
            IF (ISCALT(IPHAS).GT.0) THEN
              IVART = ISCA(ISCALT(IPHAS))
            ELSE
              WRITE(NFECRA,9010) ISCALT(IPHAS)
              CALL CSEXIT (1)
            ENDIF
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
     &             (XRTP*(XRTP*(VARAL*XRTP+VARBL)+VARCL)+VARDL)
            ENDDO
C
C
          ENDIF
C --- Tests sur ITH et ISCAVR : fin
C
        ENDDO
C --- Boucle sur les scalaires : fin
      ENDIF
C --- Test sur IUTILE : fin
C
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
     &'@      usini1 indique que la chaleur specifique est uniforme ',/,
     &'@        ICP(',I10   ,') = ',I10   ,' alors que              ',/,
     &'@      usphyv impose une chaleur specifique variable.        ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Modifier usini1 ou usphyv.                              ',/,
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
     &'@      usphyv impose une diffusivite variable.               ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Modifier usini1 ou usphyv.                              ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,
     &'@    *********                                               ',/,
     &'@    APPEL A csexit DANS LE SOUS PROGRAMME usphyv            ',/,
     &'@                                                            ',/,
     &'@    La variable dont dependent les proprietes physiques ne  ',/,
     &'@      semble pas etre une variable de calcul.               ',/,
     &'@    En effet, on cherche a utiliser la temperature alors que',/,
     &'@      ISCALT(IPHAS) = ',I10                                  ,/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Verifier le codage de usphyv (et le test lors de la     ',/,
     &'@      definition de IVART).                                 ',/,
     &'@    Verifier la definition des variables de calcul dans     ',/,
     &'@      usini1. Si un scalaire doit jouer le role de la       ',/,
     &'@      temperature, verifier que ISCALT a ete renseigne.     ',/,
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
