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
                        SUBROUTINE EBUPHY
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
     &   YFUEGF , YOXYGF , YPROGF ,
     &   YFUEGB , YOXYGB , YPROGB ,
     &   TEMP   , MASMEL ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC ROUTINE PHYSIQUE PARTICULIERE : FLAMME DE PREMELAMGE MODELE EBU
CFONC  Calcul de RHO adiabatique ou permeatique (transport de H)
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
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
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
      DOUBLE PRECISION YFUEGF(NCELET),YOXYGF(NCELET),YPROGF(NCELET)
      DOUBLE PRECISION YFUEGB(NCELET),YOXYGB(NCELET),YPROGB(NCELET)
      DOUBLE PRECISION TEMP(NCELET),MASMEL(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IPCTEM, IPCFUE, IPCOXY, IPCPRO
      INTEGER          IGG, IPHAS, IEL, IPCROM
      INTEGER          IPCKAB, IPT4, IPT3
      INTEGER          IFAC, IZONE
      INTEGER          IPBROM, IPBYCG, IPCYCG, MODE
      DOUBLE PRECISION COEFG(NGAZGM), YGFM, YGBM, EPSI
      DOUBLE PRECISION NBMOL , TEMSMM , FMEL , CKABGF, CKABGB
      DOUBLE PRECISION MASMGB, HGB, TGB, MASMGF, MASMG
C
      INTEGER       IPASS
      DATA          IPASS /0/
      SAVE          IPASS
C
C***********************************************************************
C
C=======================================================================
C 0. ON COMPTE LES PASSAGES
C=======================================================================
C
      IPASS = IPASS + 1
C
C=======================================================================
C 1. INITIALISATIONS A CONSERVER
C=======================================================================
C
C --- Initialisation memoire
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C ---> Initialisation
C
      DO IGG = 1, NGAZGM
        COEFG(IGG) = ZERO
      ENDDO
C
C ---> Positions des variables, coefficients
C
      IPHAS = 1
      IPCROM = IPPROC(IROM(IPHAS))
      IPBROM = IPPROB(IROM(IPHAS))
      IPCTEM = IPPROC(ITEMP)
      IPCFUE = IPPROC(IYM(1))
      IPCOXY = IPPROC(IYM(2))
      IPCPRO = IPPROC(IYM(3))
      IF ( IRAYPP.GT.0 ) THEN
        IPCKAB = IPPROC(ICKABS)
        IPT4 = IPPROC(IT4M)
        IPT3 = IPPROC(IT3M)
      ENDIF
C
C
C=======================================================================
C 2. DETERMINATION DES GRANDEURS THERMOCHIMIQUES
C=======================================================================
C
C ---> Grandeurs GAZ FRAIS
C
C ----- Fournies par l'utilisateur
C       FMEL         --> Taux de melange
C                        constant pour les options 0 et 1
C                        variable sinon
C       TGF          --> Temperature gaz frais en K identique
C                        pour premelange frais et dilution
C ----- Deduites
C       YFUEGF(.)    --> Fraction massique fuel gaz frais
C       YOXYGF(.)    --> Fraction massique oxydant gaz frais
C       YPROGF(.)    --> Fraction massique produits gaz frais
C       HGF          --> Enthalpie massique gaz frais identique
C                        pour premelange frais et dilution
C       MASMGF       --> Masse molaire gaz frais
C       CKABGF       --> Coefficient d'absorption
C
C ---> Grandeurs GAZ BRULES
C
C ----- Deduites
C       TGB          --> Temperature gaz brules en K
C       YFUEGB(.)    --> Fraction massique fuel gaz brules
C       YOXYGB(.)    --> Fraction massique oxydant gaz brules
C       YPROGB(.)    --> Fraction massique produits gaz brules
C       MASMGB       --> Masse molaire gaz brules
C       CKABGB       --> Coefficient d'absorption
C
C ---> Grandeurs MELANGE
C
C       MASMEL           --> Masse molaire du melange
C       PROPCE(.,IPCTEM) --> Temperature du melange
C       PROPCE(.,IPCROM) --> Masse volumique du melange
C       PROPCE(,.F,O,P ) --> Fractions massiques en F, O, P
C       PROPCE(.,IPCKAB) --> Coefficient d'absorption
C       PROPCE(.,IPT4  ) --> terme T^4
C       PROPCE(.,IPT3  ) --> terme T^3
C
C
C ---> Fractions massiques des gaz frais et brules en F, O, P
C
      DO IEL = 1, NCEL
C
        IF ( IPPMOD(ICOEBU).EQ.0 .OR. IPPMOD(ICOEBU).EQ.1 ) THEN
          FMEL = FRMEL
        ELSE
          FMEL = RTP(IEL,ISCA(IFM))
        ENDIF
C
        YFUEGF(IEL) = FMEL
        YOXYGF(IEL) = 1.D0-FMEL
        YPROGF(IEL) = 0.D0
C
        YFUEGB(IEL) = MAX(ZERO,(FMEL-FS(1))/(1.D0-FS(1)))
        YPROGB(IEL) = (FMEL-YFUEGB(IEL))/FS(1)
        YOXYGB(IEL) = 1.D0 - YFUEGB(IEL) - YPROGB(IEL)
C
      ENDDO
C
      EPSI = 1.D-06
C
      DO IEL = 1, NCEL
C
C ---> Coefficients d'absorption des gaz frais et brules
C
        IF ( IRAYPP.GT.0 ) THEN
           CKABGF = YFUEGF(IEL)*CKABSG(1) + YOXYGF(IEL)*CKABSG(2)
     &            + YPROGF(IEL)*CKABSG(3)
           CKABGB = YFUEGB(IEL)*CKABSG(1) + YOXYGB(IEL)*CKABSG(2)
     &            + YPROGB(IEL)*CKABSG(3)
        ENDIF
C
C ---> Masse molaire des gaz frais
C
        COEFG(1) = YFUEGF(IEL)
        COEFG(2) = YOXYGF(IEL)
        COEFG(3) = YPROGF(IEL)
        NBMOL = 0.D0
        DO IGG = 1, NGAZG
          NBMOL = NBMOL + COEFG(IGG)/WMOLG(IGG)
        ENDDO
        MASMGF = 1.D0/NBMOL
C
C ---> Calcul de l'enthalpie des gaz frais
C
        MODE    = -1
        CALL COTHHT
C       ===========
     &  ( MODE   , NGAZG , NGAZGM  , COEFG  ,
     &    NPO    , NPOT   , TH     , EHGAZG ,
     &    HGF    , TGF    )
C
C ---> Masse molaire des gaz brules
C
        COEFG(1) = YFUEGB(IEL)
        COEFG(2) = YOXYGB(IEL)
        COEFG(3) = YPROGB(IEL)
        NBMOL = 0.D0
        DO IGG = 1, NGAZG
          NBMOL = NBMOL + COEFG(IGG)/WMOLG(IGG)
        ENDDO
        MASMGB = 1.D0/NBMOL
C
        YGFM = RTP(IEL,ISCA(IYGFM))
        YGBM = 1.D0 - YGFM
C
C ---> Masse molaire du melange
C
        MASMEL(IEL) = 1.D0 / ( YGFM/MASMGF + YGBM/MASMGB )
C
C ---> Calcul Temperature des gaz brules
C
        IF ( IPPMOD(ICOEBU).EQ.0 .OR. IPPMOD(ICOEBU).EQ.2 ) THEN
C ---- EBU Standard et modifie en conditions adiabatiques (sans H)
          HGB = HGF
        ELSE IF ( IPPMOD(ICOEBU).EQ.1 .OR. IPPMOD(ICOEBU).EQ.3 ) THEN
C ---- EBU Standard et modifie en conditions permeatiques (avec H)
          HGB = HGF
          IF ( YGBM.GT.EPSI ) THEN
            HGB = ( RTP(IEL,ISCA(IHM))-HGF*YGFM ) / YGBM
          ENDIF
        ENDIF
C
        MODE = 1
        CALL COTHHT
C       ===========
     &  ( MODE   , NGAZG , NGAZGM  , COEFG  ,
     &    NPO    , NPOT   , TH     , EHGAZG ,
     &    HGB    , TGB    )
C
        IF ( IPPMOD(ICOEBU).EQ.0 .OR. IPPMOD(ICOEBU).EQ.2 ) THEN
C ---- EBU Standard et modifie en conditions adiabatiques (sans H)
          TGBAD = TGB
        ENDIF
C
C ---> Temperature du melange
C      Rq PPl : Il serait plus judicieux de ponderer par les CP (GF et GB)
        PROPCE(IEL,IPCTEM) = YGFM*TGF + YGBM*TGB
C
C ---> Temperature / Masse molaire
C
        TEMSMM = YGFM*TGF/MASMGF + YGBM*TGB/MASMGB
C
C ---> Masse volumique du melange
C
        IF (IPASS.GT.1.OR.(ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1)) THEN
          PROPCE(IEL,IPCROM) = SRROM*PROPCE(IEL,IPCROM)
     &                       + (1.D0-SRROM)*
     &                         ( P0(IPHAS)/(RR*TEMSMM) )
        ENDIF
C
C ---> Fractions massiques des especes globales
C
        PROPCE(IEL,IPCFUE) = YFUEGF(IEL)*YGFM + YFUEGB(IEL)*YGBM
        PROPCE(IEL,IPCOXY) = YOXYGF(IEL)*YGFM + YOXYGB(IEL)*YGBM
        PROPCE(IEL,IPCPRO) = YPROGF(IEL)*YGFM + YPROGB(IEL)*YGBM
C
C ---> Grandeurs relatives au rayonnement
C
        IF ( IRAYPP.GT.0 ) THEN
          PROPCE(IEL,IPCKAB) = YGFM*CKABGF + YGBM*CKABGB
          PROPCE(IEL,IPT4)   = YGFM*TGF**4 + YGBM*TGB**4
          PROPCE(IEL,IPT3)   = YGFM*TGF**3 + YGBM*TGB**3
        ENDIF
C
      ENDDO
C
C
C=======================================================================
C 3. CALCUL DE RHO ET DES FRACTIONS MASSIQUES DES ESPECES GLOBALES
C    SUR LES BORDS
C=======================================================================
C
C --> Masse volumique au bord
C
      IPHAS = 1
      IBROM(IPHAS) = 1
C
C ---- Masse volumique au bord pour toutes les facettes
C      Les facettes d'entree seront recalculees apres
C      a partir des CL (si IPASS > 2).
C
C ---- Au premier passage sans suite ou si on n'a pas relu la
C      masse volumique dans le fichier suite, on n'a pas recalcule la
C      masse volumique ci-dessus, pas la peine de la reprojeter aux
C      faces.
C
      IF (IPASS.GT.1.OR.(ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1)) THEN
C
        DO IFAC = 1, NFABOR
          IEL = IFABOR(IFAC)
          PROPFB(IFAC,IPBROM) = PROPCE(IEL,IPCROM)
        ENDDO
C
      ENDIF
C
C
C ---- Masse volumique au bord pour les facettes d'entree UNIQUEMENT
C      Le test sur IZONE sert pour les reprises de calcul
C      On suppose implicitement que les elements ci-dessus ont ete relus
C      dans le fichier suite (i.e. pas de suite en combustion d'un calcul
C      a froid) -> sera pris en compte eventuellement dans les versions
C      suivantes
C
      IF ( IPASS.GT.1 .OR. ISUITE.EQ.1 ) THEN
        DO IFAC = 1, NFABOR
          IZONE = IZFPPP(IFAC)
          IF(IZONE.GT.0) THEN
            IF ( IENTGB(IZONE).EQ.1 .OR. IENTGF(IZONE).EQ.1 ) THEN
              COEFG(1) = FMENT(IZONE)
              COEFG(2) = 1.D0-FMENT(IZONE)
              COEFG(3) = ZERO
              IF ( IENTGB(IZONE).EQ.1 ) THEN
                COEFG(1) = MAX(ZERO,(FMENT(IZONE)-FS(1))/(1.D0-FS(1)))
                COEFG(3) = (FMENT(IZONE)-COEFG(1))/FS(1)
                COEFG(2) = 1.D0 - COEFG(1) - COEFG(3)
              ENDIF
              NBMOL = 0.D0
              DO IGG = 1, NGAZG
                NBMOL = NBMOL + COEFG(IGG)/WMOLG(IGG)
              ENDDO
             MASMG = 1.D0/NBMOL
             TEMSMM = TKENT(IZONE)/MASMG
             PROPFB(IFAC,IPBROM) = P0(IPHAS)/(RR*TEMSMM)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
C
C
C --> Fractions massiques des especes globales au bord
C     Uniquement si transport de H
C
      IF ( IRAYPP.GT.0 ) THEN
        DO IGG = 1, NGAZG
          IPBYCG = IPPROB(IYM(IGG))
          IPCYCG = IPPROC(IYM(IGG))
          DO IFAC = 1, NFABOR
            IEL = IFABOR(IFAC)
            PROPFB(IFAC,IPBYCG) = PROPCE(IEL,IPCYCG)
          ENDDO
        ENDDO
      ENDIF
C
C=======================================================================
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
