c@a
c@versb
C-----------------------------------------------------------------------
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2009 EDF S.A., France
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
                        SUBROUTINE MEMTRI
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , IVERIF ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NCOFAB , NPROCE , NPROFA , NPROFB ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IISSTD , IFRCX  ,
     &   IDT    , IRTP   , IRTPA  , IPROPC , IPROPF , IPROPB ,
     &   ICOEFA , ICOEFB ,
     &   IFINIA , IFINRA )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC   GESTION MEMOIRE VARIABLES NON GEOMETRIQUES
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! IVERIF       ! E  !  -> ! INDICATEUR DES TESTS ELEMENTAIRES    !
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
CARGU ! NCOFAB       ! E  !  -> ! NOMBRE DE COUPLE DE CL A PREVOIR     !
CARGU ! NPROCE       ! E  !  -> ! NOMBRE DE PROP PHY AUX CENTRES       !
CARGU ! NPROFA       ! E  !  -> ! NOMBRE DE PROP PHY AUX FACES INTERNES!
CARGU ! NPROFB       ! E  !  -> ! NOMBRE DE PROP PHY AUX FACES DE BORD !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IISSTD       ! E  ! <-  ! "POINTEUR" SUR ISOSTD(REPERAGE SORTIE!
CARGU ! IDT          ! E  ! <-  ! "POINTEUR" SUR DT                    !
CARGU ! IRTP, IRTPA  ! E  ! <-  ! "POINTEUR" SUR RTP, RTPA             !
CARGU ! IPROPC       ! E  ! <-  ! "POINTEUR" SUR PROPCE                !
CARGU ! IPROPF       ! E  ! <-  ! "POINTEUR" SUR PROPFA                !
CARGU ! IPROPB       ! E  ! <-  ! "POINTEUR" SUR PROPFB                !
CARGU ! ICOEFA, B    ! E  ! <-  ! "POINTEUR" SUR COEFA, COEFB          !
CARGU ! IFRCX        ! E  ! <-  ! "POINTEUR" SUR FRCXT                 !
CARGU ! IFINIA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU ! IFINRA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU .______________.____._____.______________________________________.
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
      INCLUDE "dimfbr.h"
      INCLUDE "paramx.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "pointe.h"
      INCLUDE "albase.h"
      INCLUDE "period.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "cfpoin.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagdim.h"
      INCLUDE "lagran.h"
      INCLUDE "ihmpre.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          IVERIF
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NCOFAB , NPROCE , NPROFA , NPROFB
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IISSTD , IFRCX
      INTEGER          IDT
      INTEGER          IRTP   , IRTPA
      INTEGER          IPROPC , IPROPF , IPROPB
      INTEGER          ICOEFA , ICOEFB
      INTEGER          IFINIA , IFINRA
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IIS, IPHAS, IPPU, IPPV, IPPW, IVAR, IPROP
      INTEGER          IIPERO , IIIRIJ , IMOM, IDTNM
      INTEGER          IDPAR1 , IDPAR2 , IYPAR1, IIYPLB, IIFORB, IICOCI
C
C***********************************************************************
C
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
C
C
C
C=======================================================================
C 2. PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA
C=======================================================================
C
C --> Remarques :
C
C     IPUCOU = 1 ne depend pas de la phase
C
C     NCOFAB, NPROCE, NPROFA et NPROFB ont ete determines dans VARPOS
C         et ne servent en tant que dimensions que dans le present
C         sous programme. On pourrait les passer en common dans numvar.h
C
C     ITYPFB, ITRIFB et ISYMPA peuvent passer en entier dans certains
C         sous-pgm, il convient donc qu'ils soient en un seul bloc.
C         Le meilleur moyen de s'en assurer est de referencer ce bloc par un
C         pointeur unique independant de iphas.
C
C     Le tableau des zones frontieres des faces de bord pour les
C         physiques particulieres est de declare ci-dessous (voir PPCLIM)
C
C
C --> Preparations :
C
C     On regarde s'il existe au moins une periodicite de rotation
C                                 une phase avec Rij
C
      IIPERO = 0
      IF (IPEROT.GT.0) THEN
        IIPERO = 1
      ENDIF
      IIIRIJ = 0
      DO IPHAS = 1, NPHAS
        IF(ITYTUR(IPHAS).EQ.3) THEN
          IIIRIJ = 1
        ENDIF
      ENDDO
C
C
C     Distance a la paroi
C
C       On reserve ici idipar (distance a la paroi) : c'est oblige
C       On reserve aussi iyppar (yplus) : on pourrait s'en passer
C         et en faire un tableau local a reserver dans memdyp
C         mais ca facilite la visualisation et l'initialisation
C         de yplus pour son calcul a tous les pas de temps
C
C     Distance a la paroi, tableau selon le mode de calcul
      IDPAR1 = 0
      IDPAR2 = 0
C     Calcul par eq de diffusion
      IF(INEEDY.EQ.1.AND.ABS(ICDPAR).EQ.1) THEN
        IDPAR1 = 1
      ENDIF
C     Calcul direct
      IF(INEEDY.EQ.1.AND.ABS(ICDPAR).EQ.2) THEN
        IDPAR2 = 1
      ENDIF
C     Yplus associe (calcul par mode de diffusion et LES+VanDriest
C       ou lagrangien+IROULE=2)
      IYPAR1 = 0
      IF(INEEDY.EQ.1.AND.ABS(ICDPAR).EQ.1) THEN
        DO IPHAS = 1, NPHAS
          IF(ITYTUR(IPHAS).EQ.4) THEN
            IF(IDRIES(IPHAS).EQ.1) THEN
              IYPAR1 = 1
            ENDIF
          ENDIF
        ENDDO
        IF (IILAGR.GE.1 .AND. IROULE.EQ.2) IYPAR1 = 1
      ENDIF
C
C     Stockage supplémentaire si on initialise le gradient
C       par moindre carre
C
      IICOCI = 0
      IF(IMRGRA.EQ.4 .OR. IVERIF.EQ.1) THEN
        IICOCI = 1
      ENDIF
C
C
C     Post-traitements particuliers  (faces de bord)
C
C
      IIYPLB = 0
C     Yplus au bord
      IF(MOD(IPSTDV,IPSTYP).EQ.0) THEN
        IIYPLB = 1
      ENDIF
C     Efforts aux bords
      IIFORB = 0
      IF(INEEDF.EQ.1) THEN
        IIFORB = 1
      ENDIF
C
C
C --> Reservation de memoire entiere
C
C
      IITYPF = IDEBIA
      IITRIF = IITYPF + NFABOR *NPHAS
      IISYMP = IITRIF + NFABOR *NPHAS
      IFINIA = IISYMP + NFABOR *NPHAS
      DO IPHAS = 1, NPHAS
        IF(IDPAR2.EQ.1) THEN
          IIFAPA(IPHAS) = IFINIA
          IFINIA        = IIFAPA(IPHAS) + NCELET
        ELSE
C         cette valeur nulle est utilisee dans les tests
          IIFAPA(IPHAS) = 0
        ENDIF
      ENDDO
C
C  Zones de face de bord : on utilise provisoirement les zones des physiques
C    particulieres, meme sans physique particuliere
C    -> sera modifie lors de la restructuration des zones de bord
      IIZFPP = IFINIA
      IF(IPPMOD(IPHPAR).GE.1 .OR. IIHMPR.EQ.1) THEN
        IFINIA = IIZFPP + NFABOR
      ELSE
        IFINIA = IIZFPP
      ENDIF
C
      IISSTD = IFINIA
      IFINIA = IISSTD + (NFABOR+1)*NPHAS*IPHYDR
C
      IF(IPPMOD(ICOMPF).GE.0) THEN
        IIFBET = IFINIA
        IIFBRU = IIFBET + NFABOR*NPHAS
        IFINIA = IIFBRU + NFABOR*NPHAS
      ELSE
        IIFBET = 0
        IIFBRU = 0
      ENDIF
C
C --> Reservation de memoire reelle
C
      ICOEFA = IDEBRA
      ICOEFB = ICOEFA + NDIMFB *NCOFAB
      IRTP   = ICOEFB + NDIMFB *NCOFAB
      IRTPA  = IRTP   + NCELET *NVAR
      IPROPC = IRTPA  + NCELET *NVAR
      IPROPF = IPROPC + NCELET *NPROCE
      IPROPB = IPROPF + NFAC   *NPROFA
      IDT    = IPROPB + NDIMFB *NPROFB
      ICOCG  = IDT    + NCELET
      ICOCGB = ICOCG  + NCELET *9
      ICOCI  = ICOCGB + NCELBR *9
      ICOCIB = ICOCI  + NCELET *9 * IICOCI
      ITPUCO = ICOCIB + NCELBR *9 * IICOCI
      IDIPAR = ITPUCO + NCELET *NDIM*IPUCOU
      IYPPAR = IDIPAR + NCELET *IDPAR1
      IDUDXY = IYPPAR + NCELET *IYPAR1
      IWDUDX = IDUDXY + (NCELET-NCEL) * 3 * 3 * NPHAS * IIPERO
      IDRDXY = IWDUDX + (NCELET-NCEL) * 3 * 3 * NPHAS * IIPERO
      IWDRDX = IDRDXY + (NCELET-NCEL) * 6 * 3 * NPHAS * IIPERO*IIIRIJ
      IFRCX  = IWDRDX + (NCELET-NCEL) * 6 * 3 * NPHAS * IIPERO*IIIRIJ
      IYPLBR = IFRCX  + NCELET*NDIM*NPHAS*IPHYDR
      IFORBR = IYPLBR + NFABOR*NPHAS*IIYPLB
      IFINRA = IFORBR + NFABOR*NDIM*IIFORB
C
C     On rajoute des tableaux pour le k-omega SST si necessaire
C     En k-omega, on a besoin de calculer 2 Sij.Sij pour etre utilise
C     dans PHYVAR et dans TURBKW. On reserve un tableau pour divU en meme temps.
C     Les pointeurs IS2KW et IDVUKW sont fonction de IPHAS
C
      DO IPHAS = 1, NPHAS
        IS2KW(IPHAS)  = IFINRA
        IDVUKW(IPHAS) = IFINRA
        IF (ITURB(IPHAS).EQ.60) THEN
          IDVUKW(IPHAS) = IS2KW(IPHAS)  + NCELET
          IFINRA        = IDVUKW(IPHAS) + NCELET
        ENDIF
      ENDDO
C
C     En ALE, on reserve des tableaux supplementaires de position initiale
C       et de deplacement et le type de faces de bord
      IF (IALE.EQ.1) THEN
        IIMPAL = IFINIA
        IIALTY = IIMPAL + NNOD
        IFINIA = IIALTY + NFABOR
C
        IXYZN0 = IFINRA
        IDEPAL = IXYZN0 + NDIM*NNOD
        IFINRA = IDEPAL + NDIM*NNOD
      ELSE
        IIMPAL = 0
        IIALTY = 0
        IXYZN0 = 0
        IDEPAL = 0
      ENDIF
C
C --> Verification
C
      CALL IASIZE('MEMTRI',IFINIA)
C     ==========
C
      CALL RASIZE('MEMTRI',IFINRA)
C     ==========
C
C
C=======================================================================
C 3. CORRESPONDANCE POUR POST-TRAITEMENT
C=======================================================================
C
C --> Correspondance IPP2RA pour post-process
C       Variables de calcul et proprietes physiques
C
      DO IIS = 1 , NVPPMX
        IPP2RA(IIS) = 1
      ENDDO
C
C     IPPROC a ete complete au prealable dans VARPOS
C
      DO IPHAS = 1, NPHAS
C
        IVAR = IPR   (IPHAS)
        IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
        IVAR = IU    (IPHAS)
        IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
        IVAR = IV    (IPHAS)
        IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
        IVAR = IW    (IPHAS)
        IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
C
        IF    (ITYTUR(IPHAS).EQ.2) THEN
          IVAR = IK    (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
          IVAR = IEP   (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
          IVAR = IR11  (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
          IVAR = IR22  (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
          IVAR = IR33  (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
          IVAR = IR12  (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
          IVAR = IR13  (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
          IVAR = IR23  (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
          IVAR = IEP   (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
        ELSEIF(ITURB(IPHAS).EQ.50) THEN
          IVAR = IK    (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
          IVAR = IEP   (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
          IVAR = IPHI  (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
          IVAR = IFB   (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
        ELSEIF(ITURB(IPHAS).EQ.60) THEN
          IVAR = IK    (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
          IVAR = IOMG  (IPHAS)
          IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
        ENDIF
C
      ENDDO
C
      IF (IALE.EQ.1) THEN
        IVAR = IUMA
        IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
        IVAR = IVMA
        IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
        IVAR = IWMA
        IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
      ENDIF
C
C     Le choix fait dans VARPOS indique qu'on ne s'interessera
C       qu'aux proprietes au centre des cellules (pas au flux
C       de masse en particulier, ni a la masse volumique au bord)
C
      DO IPROP = 1, NPROCE
        IPP2RA(IPPPRO(IPROP)) = IPROPC+(IPROP-1)*NCELET
      ENDDO
C
C     Pour les moments, on repere dans IPPMOM le mode de division par le temps
C       = 0 : pas de division
C       > 0 : IPPMOM donne le pointeur dans RA sur le DT cumule
C                                                   (tableau NCEL dans PROPCE)
C       < 0 : IPPMOM donne le rang dans DTCMOM du DT cumule (uniforme)
      DO IPROP = 1, NVPPMX
        IPPMOM(IPROP) = 0
      ENDDO
      DO IMOM = 1, NBMOMT
C       Pointeur iprop des moments pour IPP2RA(IPPPRO(IPROP)) et IPPMOM(IPPPRO(IPROP))
        IPROP = IPPROC(ICMOME(IMOM))
C       Type de DT cumule et numero
        IDTNM = IDTMOM(IMOM)
        IF(IDTNM.GT.0) THEN
          IPPMOM(IPPPRO(IPROP)) =
     &         IPROPC+(IPPROC(ICDTMO(IDTNM))-1)*NCELET
        ELSEIF(IDTNM.LT.0) THEN
          IPPMOM(IPPPRO(IPROP)) = IDTNM
        ENDIF
      ENDDO
C
      DO IIS = 1 , NSCAL
        IVAR = ISCA  (IIS  )
        IPP2RA(IPPRTP(IVAR)) = IRTP  +(IVAR-1)*NCELET
      ENDDO
C
      IPP2RA(IPPDT ) = IDT

C     Couplage instationnaire vitesse/pression
      IF (IPUCOU.EQ.0) THEN
        IPP2RA(IPPTX)= 1
        IPP2RA(IPPTY)= 1
        IPP2RA(IPPTZ)= 1
      ELSE
        IPP2RA(IPPTX)= ITPUCO
        IPP2RA(IPPTY)= ITPUCO+NCELET
        IPP2RA(IPPTZ)= ITPUCO+2*NCELET
      ENDIF
C
C     Vecteur vitesse chrono
      DO IPHAS = 1, NPHAS
        IPPU = IPPRTP(IU(IPHAS))
        IPPV = IPPRTP(IV(IPHAS))
        IPPW = IPPRTP(IW(IPHAS))
        IF(ICHRVR(IPPU).EQ.1.AND.ICHRVR(IPPV).EQ.1.AND.
     &    ICHRVR(IPPW).EQ.1) THEN
          ICHRVR(IPPV) = 0
          ICHRVR(IPPW) = 0
          IPP2RA(IPPU) = - IPP2RA(IPPU)
        ENDIF
      ENDDO
C     Vecteur vitesse de maillage chrono
      IF (IALE.EQ.1) THEN
        IPPU = IPPRTP(IUMA)
        IPPV = IPPRTP(IVMA)
        IPPW = IPPRTP(IWMA)
        IF(ICHRVR(IPPU).EQ.1.AND.ICHRVR(IPPV).EQ.1.AND.
     &    ICHRVR(IPPW).EQ.1) THEN
          ICHRVR(IPPV) = 0
          ICHRVR(IPPW) = 0
          IPP2RA(IPPU) = - IPP2RA(IPPU)
        ENDIF
      ENDIF
C     Potentiel vecteur chrono
      IF(IPPMOD(IELARC).GE.2) THEN
        IPPU = IPPRTP(ISCA(IPOTVA(1)))
        IPPV = IPPRTP(ISCA(IPOTVA(2)))
        IPPW = IPPRTP(ISCA(IPOTVA(3)))
        IF(ICHRVR(IPPU).EQ.1.AND.ICHRVR(IPPV).EQ.1.AND.
     &                           ICHRVR(IPPW).EQ.1) THEN
          ICHRVR(IPPV) = 0
          ICHRVR(IPPW) = 0
          IPP2RA(IPPU) = - IPP2RA(IPPU)
        ENDIF
      ENDIF
C     Laplace vecteur chrono
      IF(IPPMOD(IELARC).GE.1) THEN
        IPPU = IPPPRO(IPPROC(ILAPLA(1)))
        IPPV = IPPPRO(IPPROC(ILAPLA(2)))
        IPPW = IPPPRO(IPPROC(ILAPLA(3)))
        IF(ICHRVR(IPPU).EQ.1.AND.ICHRVR(IPPV).EQ.1.AND.
     &                           ICHRVR(IPPW).EQ.1) THEN
          ICHRVR(IPPV) = 0
          ICHRVR(IPPW) = 0
          IPP2RA(IPPU) = - IPP2RA(IPPU)
        ENDIF
      ENDIF
C
      RETURN
      END
c@z
