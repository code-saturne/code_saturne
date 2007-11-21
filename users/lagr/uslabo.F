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
                        SUBROUTINE USLABO
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   KFACE  , NBPT   , ISUIVI ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , ITRIFB , IFRLAG , ITEPA  , INDEP  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   SURFBN , DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , ETTPA  , TEPA   , PARBOR , VITPAR , VITFLU , AUXL   ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC   SOUS-PROGRAMME UTILISATEUR (INTERVENTION NON OBLIGATOIRE)
CFONC
CFONC   SOUS-PROGRAMME UTILISATEUR DE GESTION DU COMPORTEMENT
CFONC     DES PARTICULES A UNE INTERACTION PARTICULE/FRONTIERE
CFONC     ET D'ENREGISTREMENT DES STATISTIQUES AUX FRONTIERES.
CFONC
CFONC   L'UTILISATEUR N'A PAS A MODIFIER LE PRESENT SOUS-PROGRAMME DANS
CFONC     LES CONDITIONS D'UTILISATION STANDARD. DANS LE CAS OU IL
CFONC     SOUHAITE AVOIR DES INTERACTIONS NON PREVU, IL DOIT INTERVENIR
CFONC     DANS LES RUBRIQUES NUMERO 8 et 10.
CFONC
CFONC   Gestion de l'interaction particule/face de frontiere
CFONC     en fonction des informations donnees par l'utilisateur
CFONC     (valeur de IUSCLB par zone) dans la routine USLAG2.
CFONC
CFONC   En fonction du nom stocke dans IUSCLB et affecte
CFONC     a la face de bord KFACE, on donne le type de comportement
CFONC     des particules.
CFONC
CFONC   En standard, on a :
CFONC
CFONC     IUSCLB conditions au bord pour les particules
CFONC     = IENTRL -> zone d'injection de particules
CFONC     = ISORTL -> sortie du domaine
CFONC     = IREBOL -> rebond des particules
CFONC     = IDEPO1 -> deposition definitive
CFONC     = IDEPO2 -> deposition definitive mais la particule reste en
CFONC                 memoire (utile si IENSI2 = 1 uniquement)
CFONC     = IDEPO3 -> deposition et remise en suspension possible
CFONC                 suivant les condition de l'ecoulement
CFONC     = IENCRL -> encrassement (Charbon uniquement IPHYLA = 2)
CFONC
CFONC
CFONC   De plus, si on souhaite un autre type d'interaction non prevue
CFONC     pour zone de faces de bord, dans USLAG2 on affecte dans
CFONC     IUSCLB(KZONE) un des noms suivants :
CFONC
CFONC                JBORD1, JBORD2, JBORD3, JBORD4, JBORD5
CFONC
CFONC     et dans cette routine on code le comportement
CFONC     des particules pour cette zone frontiere.
CFONC
CFONC   ATTENTION : en debut de routine la variable ISUIVI est
CFONC    initialisee a une valeur aberrante et DOIT etre modifiee
CFONC    avant la fin de la routine.
CFONC    Les vitesses de la particule et du fluide vu doivent etre
CFONC    modifiees en fonction de l'interaction via les tableaux
CFONC    VITPAR et VITFLU, elles NE DOIVENT PAS etre modifiees
CFONC    via ETTP et ETTPA dans ce sous-programme.
CFONC
CFONC   Regle de modification de ISUIVI :
CFONC   =================================
CFONC  1) Mettre ISUIVI = 0 si la particule ne doit pas etre
CFONC     suivi dans le maillage apres l'interaction de sa trajectoire
CFONC     et de la face de bord (ex : IENTRL, ISORTL, IDEPO1, IDEPO2).
CFONC
CFONC  2) mettre ISUIVI = 1 pour continuer a suivre la particule
CFONC     dans le maillage apres son interaction (ex : IREBOL, IDEPO3).
CFONC     On peut tres bien avoir ISUIVI = 0 ou ISUIVI = 1 pour un type
CFONC     d'interaction en fonction de comportement de la particule
CFONC     (ex : IENCRL).
CFONC
CFONC  REMARQUE : Lorsqu'il y a interaction, les calculs de la vitesse
CFONC  ^^^^^^^^   de la particule et de la vitesse du fluide vu sont
CFONC             forcement a l'ordre 1 (meme si on est a l'ordre 2
CFONC             par ailleurs).
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
CARGU ! LNDNOD       ! E  !  -> ! DIM. CONNECTIVITE CELLULES->FACES    !
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC           !
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR           !
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NTERSL       ! E  !  -> ! NBR TERMES SOURCES DE COUPLAGE RETOUR!
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VAR STATISTIQUES LAGRANGIEN!
CARGU ! NVISBR       ! E  !  -> ! NOMBRE DE STATISTIQUES AUX FRONTIERES!
CARGU ! KFACE        ! E  !  -> ! NUMERO DE LA FACE D'INTERACTION      !
CARGU ! NBPT         ! E  !  -> ! NUMERO DE LA PARTICULE TRAITEE       !
CARGU ! ISUIVI       ! E  ! <-  ! INDICATEUR DE SUIVI DE LA PARTICULE  !
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
CARGU !  (NFML,NPRFML!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC            !
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !                                      !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR            !
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !                                      !
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS)      !    !     !                                      !
CARGU ! ITRIFB(NFABOR! TE !  -> ! TAB D'INDIRECTION POUR TRI DES FACES !
CARGU !   NPHAS)     !    !     !                                      !
CARGU ! IFRLAG       ! TE !  -> ! NUMERO DE ZONE DE LA FACE DE BORD    !
CARGU !   (NFABOR)   !    !     !  POUR LE MODULE LAGRANGIEN           !
CARGU ! ITEPA        ! TE ! <-  ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
CARGU ! INDEP        ! TE ! <-  ! POUR CHAQUE PARTICULE :              !
CARGU !   (NBPMAX)   !    !     !   NUMERO DE LA CELLULE DE DEPART     !
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
CARGU ! VOLUME(NCELET! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! SURFBN(NFABOR! TR !  -> ! SURFACE DES FACES DE BORD            !
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT ET PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! ETTP         ! TR ! <-> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! PARBOR(NFABOR! TR ! <-> ! CUMUL DES STATISTIQUES AUX FRONTIERES!
CARGU !    NVISBR)   !    !     !                                      !
CARGU ! VITPAR       ! TR ! <-> ! VITESSE PARTICULE POUR LE TRAITEMENT !
CARGU !   (NBPMAX,3) !    !     !   INTERACTIONS PARTICULES/FRONTIERES !
CARGU ! VITFLU       ! TR ! <-> ! VITESSE FLUIDE VU POUR LE TRAITEMENT !
CARGU !   (NBPMAX,3) !    !     !   INTERACTIONS PARTICULES/FRONTIERES !
CARGU ! AUXL(NBPMAX,3! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "cstnum.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "cstphy.h"
      INCLUDE "pointe.h"
      INCLUDE "period.h"
      INCLUDE "parall.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "cpincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDNOD , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          KFACE  , NBPT   , ISUIVI
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1) , NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1) , NODFBR(LNDFBR)
      INTEGER          ITYPFB(NFABOR,NPHAS) , ITRIFB(NFABOR,NPHAS)
      INTEGER          IFRLAG(NFABOR) , ITEPA(NBPMAX,NIVEP)
      INTEGER          INDEP(NBPMAX)
      INTEGER          IDEVEL(NIDEVE) , ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD) , VOLUME(NCELET)
      DOUBLE PRECISION SURFBN(NFABOR)
      DOUBLE PRECISION DT(NCELET) , RTP(NCELET,*) , RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*) , COEFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION PARBOR(NFABOR,NVISBR) , AUXL(NBPMAX,3)
      DOUBLE PRECISION VITPAR(NBPMAX,3) , VITFLU(NBPMAX,3)
      DOUBLE PRECISION RDEVEL(NRDEVE) , RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IP , NFIN , KZONE , N1 , ICHA, IOK
C
      DOUBLE PRECISION AA
      DOUBLE PRECISION XP , YP , ZP
      DOUBLE PRECISION XQ , YQ , ZQ
      DOUBLE PRECISION XK , YK , ZK
      DOUBLE PRECISION XPQ , YPQ , ZPQ
      DOUBLE PRECISION XNN , YNN , ZNN
      DOUBLE PRECISION VNORL(1)  , ENC3 , VISCP , MASSE
      DOUBLE PRECISION DPINIT , DP03 , MP0 , TRAP , VNORM , ANG
      DOUBLE PRECISION ENERG , ENERGT
      DOUBLE PRECISION UXN   , VYN    , WZN
C
C***********************************************************************
C
C=======================================================================
C 0.  GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. Traitement en fonction du type de frontiere
C=======================================================================
C
      IOK = 0
C
C--> numero de la particule traitee
C
      IP = NBPT
C
C--> Zone de la face de bord a traitee
C
      KZONE = IFRLAG(KFACE)
C
C-->normale normee sortante de la face KFACE
C
      AA = 1.D0 / SURFBN(KFACE)
      XNN = SURFBO(1,KFACE) * AA
      YNN = SURFBO(2,KFACE) * AA
      ZNN = SURFBO(3,KFACE) * AA
C
C=======================================================================
C 2. Recherche du point d'intersection entre la face de bord et
C    le rayon. Les coordonnees sont stockees dans XK YK ZK
C=======================================================================
C
C Petit rappel de geometrie 3D :
C ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C            1)  equation d'un plan de normal (a,b,c) :
C                a*x + b*y + c*z + d = 0
C            2)  equation d'une droite qui passe par P et Q :
C                x = XP + (XQ-XP) * AA
C                y = YP + (YQ-YP) * AA
C                z = ZP + (ZQ-ZP) * AA
C                ou AA est un parametre qui varie dans l'ensemble des reels
C
C-->on determine le vecteur PQ
C
      XP = ETTPA(IP,JXP)
      YP = ETTPA(IP,JYP)
      ZP = ETTPA(IP,JZP)
C
      XQ = ETTP(IP,JXP)
      YQ = ETTP(IP,JYP)
      ZQ = ETTP(IP,JZP)
C
      XPQ = XQ - XP
      YPQ = YQ - YP
      ZPQ = ZQ - ZP
C
C-->si la particule n'a pas bougee ( si elle est deposee
C   sur la face de bord), elle n'est pas traitee.
C
      IF (XPQ.EQ.0.D0 .AND. YPQ.EQ.0.D0 .AND. ZPQ.EQ.0.D0) RETURN
C
C-->a partir de l'equation du plan de la face et de l'equation
C   parametre du rayon, on determine le parametre du point d'intersection
C
      AA = XPQ * SURFBO(1,KFACE)
     &   + YPQ * SURFBO(2,KFACE)
     &   + ZPQ * SURFBO(3,KFACE)
C
      IF ( AA.EQ.0.D0 ) THEN
        WRITE (NFECRA,9010) IP
        NBPERR = NBPERR + 1
        DNBPER = DNBPER + TEPA(IP,JRPOI)
        ISUIVI = 0
        ITEPA(IP,JISOR) = 0
        RETURN
      ENDIF
C
      AA =
     &     ( SURFBO(1,KFACE) * CDGFBO(1,KFACE)
     &     + SURFBO(2,KFACE) * CDGFBO(2,KFACE)
     &     + SURFBO(3,KFACE) * CDGFBO(3,KFACE)
     &     - SURFBO(1,KFACE) * XP
     &     - SURFBO(2,KFACE) * YP
     &     - SURFBO(3,KFACE) * ZP )
     &     / AA
C
C-->on reporte ce parametre dans l'equation de la droite du rayon pour
C   avoir le point d'intersection (XK YK ZK)
C
      XK = XP + XPQ * AA
      YK = YP + YPQ * AA
      ZK = ZP + ZPQ * AA
C
C=======================================================================
C 3. Depart de la particule du domaine de calcul,
C    ou deposition de la particule sur la frontiere
C=======================================================================
C
      IF ( IUSCLB(KZONE).EQ.ISORTL .OR.
     &     IUSCLB(KZONE).EQ.IENTRL .OR.
     &     IUSCLB(KZONE).EQ.IDEPO1      ) THEN
C
        ISUIVI = 0
        ITEPA(IP,JISOR) = 0
C
C      mise a jour du debit
C
        DEBLAG(KZONE) = DEBLAG(KZONE)-TEPA(IP,JRPOI)*ETTP(IP,JMP)
C
C--> La particule sort mais pour la visualisation ensight, on la place
C    correctement au point d'intersection.
C    Possibilite que le numero IP soit reutilise pour une autre
C    particule a entrer.
C
        ETTP(IP,JXP) = XK
        ETTP(IP,JYP) = YK
        ETTP(IP,JZP) = ZK
C
C=======================================================================
C 4. Deposition de la particule, celle-ci reste en memoire
C=======================================================================
C
      ELSE IF (IUSCLB(KZONE).EQ.IDEPO2) THEN
C
C--> La particule ne sort pas, elle n'est plus traitee,
C    mais toujours visualisee. Le numero IP n'est pas reutilisable.
C
          ISUIVI = 0
          ITEPA(IP,JISOR) = -ITEPA(IP,JISOR)
        ETTP(IP,JXP) = XK
        ETTP(IP,JYP) = YK
        ETTP(IP,JZP) = ZK
C
        DO N1 = 1,3
          VITPAR(IP,N1) = 0.D0
          VITFLU(IP,N1) = 0.D0
        ENDDO
C
C=======================================================================
C 5. Deposition de la particule, la remise en suspension est possible
C=======================================================================
C
      ELSE IF (IUSCLB(KZONE).EQ.IDEPO3) THEN
C
        ISUIVI = 0
        ITEPA(IP,JISOR) = IFABOR(KFACE)
        ETTP(IP,JXP) = XK
        ETTP(IP,JYP) = YK
        ETTP(IP,JZP) = ZK
C
        DO N1 = 1,3
          VITPAR(IP,N1) = 0.D0
          VITFLU(IP,N1) = 0.D0
        ENDDO
        DO N1 = JUP,JWF
          ETTPA(IP,N1) = 0.D0
        ENDDO
C
C=======================================================================
C 6. Deposition de la particule avec force d'attachement,
C    vitesse conservee et re-entrainement possible
C=======================================================================
C
      ELSE IF (IUSCLB(KZONE).EQ.IDEPFA) THEN
C
C
C Calcul du critere
C
        UXN = ETTP(IP,JUP)*XNN
        VYN = ETTP(IP,JVP)*YNN
        WZN = ETTP(IP,JWP)*ZNN
C
        ENERG = 0.5D0*ETTP(IP,JMP)*(UXN*UXN+VYN*VYN+WZN*WZN)
C
        ENERGT   = 3.34D-12*ETTP(IP,JDP)
C
        IF ( ENERG .GE. ENERGT )THEN
C
          ISUIVI = 0
          ITEPA(IP,JISOR) = IFABOR(KFACE)
          ETTP(IP,JXP) = XK
          ETTP(IP,JYP) = YK
          ETTP(IP,JZP) = ZK
C
          VITPAR(IP,1) = 0.D0
          VITPAR(IP,2) = 0.D0
          VITPAR(IP,3) = 0.D0
C
          VITFLU(IP,1) = 0.D0
          VITFLU(IP,2) = 0.D0
          VITFLU(IP,3) = 0.D0
C
        ELSE
C
          ISUIVI = 0
          ITEPA(IP,JISOR) = INDEP(IP)
          ETTP(IP,JXP) = ETTPA(IP,JXP)
          ETTP(IP,JYP) = ETTPA(IP,JYP)
          ETTP(IP,JZP) = ETTPA(IP,JZP)
C
C-->changement de la vitesse de la particule au point d'arrive
C   (comme pour une rebond elastique)
C
          AA = ABS(( VITPAR(IP,1)*XNN
     &              +VITPAR(IP,2)*YNN
     &              +VITPAR(IP,3)*ZNN) )*2.D0
C
          VITPAR(IP,1) = VITPAR(IP,1) - AA*XNN
          VITPAR(IP,2) = VITPAR(IP,2) - AA*YNN
          VITPAR(IP,3) = VITPAR(IP,3) - AA*ZNN
C
C-->Annule  la vitesse du fluide vu au point d'arrive
C   (comme pour une rebond elastique)
C
          AA = ABS( (VITFLU(IP,1)*XNN
     &             + VITFLU(IP,2)*YNN
     &             + VITFLU(IP,3)*ZNN) ) * 2.D0
C
          VITFLU(IP,1) = 0.D0
          VITFLU(IP,2) = 0.D0
          VITFLU(IP,3) = 0.D0
        ENDIF
C
C=======================================================================
C 7. Rebond elastique de la particule sur la frontiere
C=======================================================================
C
      ELSE IF (IUSCLB(KZONE).EQ.IREBOL) THEN
C
        ISUIVI = 1
        ITEPA(IP,JISOR) = IFABOR(KFACE)
C
C-->changement du point de depart
C
        ETTPA(IP,JXP) = XK
        ETTPA(IP,JYP) = YK
        ETTPA(IP,JZP) = ZK
C
          IF (IENSI1.EQ.1) THEN
            NFIN = 0
            CALL ENSLAG
C           ===========
     &       ( IDEBIA, IDEBRA  ,
     &         NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &         NFIN   , IP     ,
     &         ITEPA  ,
     &         ETTPA  , TEPA   , RA)
          ENDIF
C
C-->changement du point d'arrive
C   (la valeur absolue est pour eviter les p. scalaires
C    negatifs impossibles qui interviennent avec des erreurs d'arrondi
C    machines)
C
          AA = 2.D0 * ABS( (XQ-XK)*XNN + (YQ-YK)*YNN + (ZQ-ZK)*ZNN )
C
        ETTP(IP,JXP) = XQ - AA*XNN
        ETTP(IP,JYP) = YQ - AA*YNN
        ETTP(IP,JZP) = ZQ - AA*ZNN
C
C-->changement de la vitesse de la particule au point d'arrive
C
C
        AA = ABS( (VITPAR(IP,1)*XNN
     &          +  VITPAR(IP,2)*YNN
     &          +  VITPAR(IP,3)*ZNN) ) * 2.D0
C
        VITPAR(IP,1) = VITPAR(IP,1) - AA*XNN
        VITPAR(IP,2) = VITPAR(IP,2) - AA*YNN
        VITPAR(IP,3) = VITPAR(IP,3) - AA*ZNN
C
C-->changement de la vitesse du fluide vu au point d'arrive
C
        AA = ABS( (VITFLU(IP,1)*XNN
     &          +  VITFLU(IP,2)*YNN
     &          +  VITFLU(IP,3)*ZNN) ) * 2.D0
C
        VITFLU(IP,1) = VITFLU(IP,1) - AA*XNN
        VITFLU(IP,2) = VITFLU(IP,2) - AA*YNN
        VITFLU(IP,3) = VITFLU(IP,3) - AA*ZNN
C
C=======================================================================
C 8. Encrassement des grains de charbon
C=======================================================================
C
      ELSE IF (IUSCLB(KZONE).EQ.IENCRL) THEN
C
C--> Encrassement de la particules si ses proprietes le permettent
C    et en fonction d'une probabilite
C      ICI c'est si Tp     > TPENC
C                si VISCP  > VISCREF
C
        ICHA = ITEPA(IP,JINCH)
C
        IF ( ETTP(IP,JHP).GT.TPRENC(ICHA) ) THEN
C
          ENC3 = ( (1.D+7 * ENC1(ICHA))/((ETTP(IP,JHP)-150.D0)**2) )
     &           + ENC2(ICHA)
            VISCP = EXP( LOG(10.D0)*ENC3 ) * 0.1D0
C
            IF ( VISCP.GT.VISREF(ICHA) ) THEN
              N1 = 1
              CALL ZUFALL(N1,VNORL(1))
              TRAP = 1.D0-VISREF(ICHA) / VISCP
            ENDIF
C
C  Si VISCP <= VISREF ===> 100% de chance d'encrasser
C  Si VISCP  > VISREF ===> probabilte TRAP = 1-VISREF/VISCP
C                          d'encrasser
C                     ===> On encrasse si VNORL est compris
C                          entre TRAP et 1.
C
            IF ( VISCP.LE.VISREF(ICHA) .OR.
     &         (VISCP.GT.VISREF(ICHA) .AND. VNORL(1).GE.TRAP) ) THEN
C
C La calcul de la masse de grains de charbon encrasses est faite plus bas.
C
              NPENCR = NPENCR + 1
              ISUIVI = 0
              ITEPA(IP,JISOR)  =  0
            ETTP(IP,JXP) = XK
            ETTP(IP,JYP) = YK
            ETTP(IP,JZP) = ZK
C
            ENDIF
          ENDIF
C
C--> Si pas encrassement alors rebond elastique
C
          IF ( ITEPA(IP,JISOR).NE.0 ) THEN
C
            ISUIVI = 1
            ITEPA(IP,JISOR) = IFABOR(KFACE)
C
C-->changement du point de depart
C
          ETTPA(IP,JXP) = XK
          ETTPA(IP,JYP) = YK
          ETTPA(IP,JZP) = ZK
C
            IF (IENSI1.EQ.1) THEN
              NFIN = 0
              CALL ENSLAG
C             ===========
     &       ( IDEBIA, IDEBRA  ,
     &         NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &         NFIN   , IP     ,
     &         ITEPA  ,
     &         ETTPA  , TEPA   , RA)
            ENDIF
C
C-->changement du point d'arrive
C
            AA = 2.D0 * ABS((XQ-XK)*XNN + (YQ-YK)*YNN + (ZQ-ZK)*ZNN)
C
          ETTP(IP,JXP) = XQ - AA*XNN
          ETTP(IP,JYP) = YQ - AA*YNN
          ETTP(IP,JZP) = ZQ - AA*ZNN
C
          ENDIF
C
        IF (ITEPA(IP,JISOR).GT.0) THEN
C
C-->changement de la vitesse de la particule au point d'arrive
C
C
          AA = ABS( (VITPAR(IP,1)*XNN
     &            +  VITPAR(IP,2)*YNN
     &            +  VITPAR(IP,3)*ZNN) ) * 2.D0
C
          VITPAR(IP,1) = VITPAR(IP,1) - AA*XNN
          VITPAR(IP,2) = VITPAR(IP,2) - AA*YNN
          VITPAR(IP,3) = VITPAR(IP,3) - AA*ZNN
C
C-->changement de la vitesse du fluide vu au point d'arrive
C
          AA = ABS( (VITFLU(IP,1)*XNN
     &            +  VITFLU(IP,2)*YNN
     &            +  VITFLU(IP,3)*ZNN) ) * 2.D0
C
          VITFLU(IP,1) = VITFLU(IP,1) - AA*XNN
          VITFLU(IP,2) = VITFLU(IP,2) - AA*YNN
          VITFLU(IP,3) = VITFLU(IP,3) - AA*ZNN
C
          ENDIF
C
C=======================================================================
C 9. Interaction utilisateur numero 1 : JBORD1
C=======================================================================
C
C  ON PEUT FAIRE DE MEME AVEC JBORD2, JBORD3, JBORD4 et JBORD5
C    On ne donnne l'exemple que pour JBORD1
C
C     On regarde si on est dans la zone cherchee
C      ELSE IF (IUSCLB(KZONE).EQ.JBORD1) THEN
C
C       si on doit continuer a suivre la particule
C         ISUIVI = 0 OU 1
C
C       l'element de maillage concerne
C         ITEPA(IP,JISOR) =
C
C       on change la localisation de depart de la particule
C         ETTP(IP,JXP) =
C         ETTP(IP,JYP) =
C         ETTP(IP,JZP) =
C
C       la vitesse de la particule au point d'arrivee
C         VITPAR(IP,1) =
C         VITPAR(IP,2) =
C         VITPAR(IP,3) =
C
C       la vitesse du fluide vu au point d'arrivee
C         VITFLU(IP,1) =
C         VITFLU(IP,2) =
C         VITFLU(IP,3) =
C
C
      ELSE IF (IUSCLB(KZONE).EQ.JBORD1
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C    en standard, sans intervention de l'utilisateur,
C      on ne souhaite pas passer ici
C    mais on desire
C      que le test IUSCLB(KZONE).EQ.JBORD1 soit dans le us* exemple
C      que le source ci-dessous soit compile pour verifier
C         qu'il n y a pas d'erreur
     &         .AND.(0.EQ.1)
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
     &                                ) THEN
C
C     ----------------------------------------------------
C     EXEMPLE 1 : LA PARTICULE A UNE CHANCE SUR DEUX (TRAP)
C                 DE SE DEPOSER DEFINITIVEMENT A LA PAROI
C                 OU DE REBONDIR VERS L'ECOULEMENT.
C     ----------------------------------------------------
C
C
          N1 = 1
          CALL ZUFALL(N1,VNORL(1))
          TRAP = 0.5D0
C
          IF (VNORL(1).GE.TRAP) THEN
C
            ISUIVI = 0
            ITEPA(IP,JISOR)  =  0
            ETTP(IP,JXP) = XK
            ETTP(IP,JYP) = YK
            ETTP(IP,JZP) = ZK
C
          ELSE
C
            ISUIVI = 1
            ITEPA(IP,JISOR) = IFABOR(KFACE)
C
C-->changement du point de depart
C
            ETTPA(IP,JXP) = XK
            ETTPA(IP,JYP) = YK
            ETTPA(IP,JZP) = ZK
C
C-->changement du point d'arrive
C
            AA = 2.D0 * ABS((XQ-XK)*XNN + (YQ-YK)*YNN + (ZQ-ZK)*ZNN)
C
            ETTP(IP,JXP) = XQ - AA*XNN
            ETTP(IP,JYP) = YQ - AA*YNN
            ETTP(IP,JZP) = ZQ - AA*ZNN
C
          ENDIF
C
C-->Inutile de traiter les particule ITEPA(IP,JISOR)=0 car
C   Elle vont etre eliminees de la listes des particules.
C
          IF (ITEPA(IP,JISOR).GT.0) THEN
C
C-->changement de la vitesse de la particule au point d'arrive
C
C
            AA = ABS( (VITPAR(IP,1)*XNN
     &              +  VITPAR(IP,2)*YNN
     &              +  VITPAR(IP,3)*ZNN) ) * 2.D0
C
            VITPAR(IP,1) = VITPAR(IP,1) - AA*XNN
            VITPAR(IP,2) = VITPAR(IP,2) - AA*YNN
            VITPAR(IP,3) = VITPAR(IP,3) - AA*ZNN
C
C-->changement de la vitesse du fluide vu au point d'arrive
C
            AA = ABS( (VITFLU(IP,1)*XNN
     &              +  VITFLU(IP,2)*YNN
     &              +  VITFLU(IP,3)*ZNN) ) * 2.D0
C
            VITFLU(IP,1) = VITFLU(IP,1) - AA*XNN
            VITFLU(IP,2) = VITFLU(IP,2) - AA*YNN
            VITFLU(IP,3) = VITFLU(IP,3) - AA*ZNN
C
          ENDIF
C
C
C=======================================================================
C 10. Verification et sortie si erreur
C=======================================================================
C
      ELSE
        WRITE (NFECRA,9020) KZONE
        IOK = IOK + 1
      ENDIF
C
      IF (IOK.NE.0) THEN
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C=======================================================================
C 11. Enregistrement de l'interaction particule-frontiere s'il y a lieu
C=======================================================================
C
C     L'enregistrement des statistiques parietales debute des que
C     l'indicateur est a IENSI3 = 1. Cependant tant que le numero
C     de l'iteration Lagrangienne absolue est inferieur a NSTBOR,
C     ou que l'ecoulement est instationnaire (ISTTIO = 0), alors
C     le tableau PARBOR est remis a zero avant d'entrer dans ce
C     sous-programme.
C
C     NPSTF   : nombre d'iterations de calcul de stat aux frontieres
C               stationnaires
C
C     NPSTFT  : nombre d'iterations total des stats iaux frontieres
C               depuis le du calcul, partie instationnaire comprise
C               (a utiliser que pour les affichages listing)
C
C     TSTATP  : Temps physique d'enregistrement des statistiques
C               des interactions particules/frontiere stationnaires,
C               si instatinnaire alors contient DTP le dernier
C               pas de temps Lagrangien
C
C     Avant impression dans le listing, ou sortie pour le
C     post-processing, on applique aux statistiques sur les frontieres
C     une moyenne en fonction des infos fournies dans USLAG1 de la
C     maniere suivante :
C
C
C
C     CES LIGNES NE SONT DONNEES QU'A TITRE INDICATIF
C
C
C     DO IVAR = 1,NVISBR
C
C      IF (IMOYBR(IVAR).EQ.2) THEN
C
C         DO IFAC = 1,NFABOR
C           IF (PARBOR(IFAC,INBR).GT.SEUILF) THEN
C             PARBOR(IFAC,IVAR) = PARBOR(IFAC,IVAR) /PARBOR(IFAC,INBR)
C           ELSE
C             PARBOR(IFAC,IVAR) = 0.D0
C           ENDIF
C         ENDDO
C
C       ELSE IF (IMOYBR(IVAR).EQ.1) THEN
C
C         DO IFAC = 1,NFABOR
C           IF (PARBOR(IFAC,INBR).GT.SEUILF) THEN
C             PARBOR(IFAC,IVAR) = PARBOR(IFAC,IVAR) / TSTATP
C           ELSE
C             PARBOR(IFAC,IVAR) = 0.D0
C           ENDIF
C         ENDDO
C       ENDIF
C     ENDDO
C
C
C
C
      IF ( IENSI3.EQ.1 ) THEN
C
C--> exemple de types d'interaction sur lesquels on souhaite
C    enregistrer des informations
C
        IF ( IUSCLB(KZONE).EQ.IREBOL .OR.
     &       IUSCLB(KZONE).EQ.IDEPO1 .OR.
     &       IUSCLB(KZONE).EQ.IDEPO2 .OR.
     &       IUSCLB(KZONE).EQ.IDEPO3      ) THEN
C
          IF (INBRBD.EQ.1) THEN
            PARBOR(KFACE,INBR) = PARBOR(KFACE,INBR) + TEPA(IP,JRPOI)
          ENDIF
C
          IF (IFLMBD.EQ.1) THEN
              PARBOR(KFACE,IFLM) = PARBOR(KFACE,IFLM)
     &     + ( TEPA(IP,JRPOI) * ETTP(IP,JMP) /SURFBN(KFACE) )
          ENDIF
C
          IF (IANGBD.EQ.1) THEN
            VNORM = ETTP(IP,JUP) * ETTP(IP,JUP)
     &            + ETTP(IP,JVP) * ETTP(IP,JVP)
     &            + ETTP(IP,JWP) * ETTP(IP,JWP)
            VNORM = SQRT( VNORM )
            ANG =  ETTP(IP,JUP) * SURFBO(1,KFACE)
     &           + ETTP(IP,JVP) * SURFBO(2,KFACE)
     &           + ETTP(IP,JWP) * SURFBO(3,KFACE)
     &           / SURFBN(KFACE)
     &           / VNORM
            ANG = ACOS(ANG)
            PARBOR(KFACE,IANG) = PARBOR(KFACE,IANG) + ANG*TEPA(IP,JRPOI)
          ENDIF
C
          IF (IVITBD.EQ.1) THEN
            VNORM = ETTP(IP,JUP) * ETTP(IP,JUP)
     &            + ETTP(IP,JVP) * ETTP(IP,JVP)
     &            + ETTP(IP,JWP) * ETTP(IP,JWP)
            VNORM = SQRT( VNORM )
            PARBOR(KFACE,IVIT) =PARBOR(KFACE,IVIT) +VNORM*TEPA(IP,JRPOI)
          ENDIF
C
          IF (NUSBOR.GT.0) THEN
            DO N1 = 1,NUSBOR
              PARBOR(KFACE,IUSB(N1)) = 0.D0
            ENDDO
          ENDIF
C
C--> Cas particulier de la masse de grains de charbon encrasses
C
        ELSE IF ( IUSCLB(KZONE).EQ.IENCRL .AND. ISUIVI.EQ.0 ) THEN
C
          PARBOR(KFACE,INBR) = PARBOR(KFACE,INBR) + TEPA(IP,JRPOI)
C
          IF (IENCBD.EQ.1) THEN
C
            ICHA =  ITEPA(IP,JINCH)
            DPINIT = TEPA(IP,JRD0P)
C
            DP03 = DPINIT * DPINIT * DPINIT
            MP0  = PI * DP03 * RHO0CH(ICHA) / 6.D0
C
            MASSE = ETTP(IP,JMCH) + ETTP(IP,JMCK)
     &                                + XASHCH(ICHA) * MP0
C
            PARBOR(KFACE,IENC) = PARBOR(KFACE,IENC)
     &                         + TEPA(IP,JRPOI)*MASSE
C
          ENDIF
C
        ENDIF
C
      ENDIF
C
C
C=======================================================================
C Archives. Je laisse cette partie au cas ou...
C           Creation d'un repere local associe a la face frontiere
C=======================================================================
C
C Le repere local (T1,T2,N) est construit de facon a ce que N soit
C la normale normee de la face et que T1 et T2 appartiennent a la face
C
C-->1. Je connais les vecteurs N et PK, je definis T1 tel que
C      T1 = PK vectoriel N
C
c     XPK = XK - ETTPA(IP,JXP)
c     YPK = YK - ETTPA(IP,JYP)
c     ZPK = ZK - ETTPA(IP,JZP)
C
c     XT1 = YPK*ZNN - ZPK*YNN
c     YT1 = ZPK*XNN - XPK*ZNN
c     ZT1 = XPK*YNN - YPK*XNN
C
c     AA = SQRT(XT1*XT1 + YT1*YT1 + ZT1*ZT1)
c     XT1 = XT1 / AA
c     YT1 = YT1 / AA
c     ZT1 = ZT1 / AA
C
C-->2. Maintenant je peux construire T2 = - T1 vectoriel N
C
c     XT2 = YT1*ZNN - ZT1*YNN
c     YT2 = ZT1*XNN - XT1*ZNN
c     ZT2 = XT1*YNN - YT1*XNN
C
c     AA = SQRT(XT2*XT2 + YT2*YT2 + ZT2*ZT2)
c     XT2 = -XT2 / AA
c     YT2 = -YT2 / AA
c     ZT2 = -ZT2 / AA
C
C--------
C FORMATS
C--------
C
 9010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (USLABO)                                    ',/,
     &'@                                                            ',/,
     &'@  La normale d''une face de frontiere est perpendiculaire   ',/,
     &'@   a un rayon PQ : impossible                               ',/,
     &'@                                                            ',/,
     &'@  La particle ',I10,' est ELIMINEE                          ',/,
     &'@                                                            ',/,
     &'@  Contacter l''equipe de developpement.                     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (USLABO)                                    ',/,
     &'@                                                            ',/,
     &'@  Le type de condition aux limites IUSCLB                   ',/,
     &'@    n''est pas renseigne pour la frontiere NB = ',I10        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier USLAG2 et USLABO.                                ',/,
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
