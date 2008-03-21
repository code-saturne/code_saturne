c@a
c@versb
C-----------------------------------------------------------------------
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
                        SUBROUTINE CPPHYV
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
     &   W1     , W2     , W3     , W4     ,
     &   W5     , W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC ROUTINE PHYSIQUE PARTICULIERE : COMBUSTION CHARBON PULVERISE
CFONC
CFONC  Calcul de RHO du melange
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
CFONC       - dans usppiv :
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
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET),W4(NCELET)
      DOUBLE PRECISION W5(NCELET),W6(NCELET),W7(NCELET),W8(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          NTBCPI, ICPWI, NTBCPR, ICPWR
      INTEGER          NTBMCI, IMCWI, NTBMCR, IMCWR
      INTEGER          NTBWOI, IWORI, NTBWOR, IWORR
      INTEGER          IFINIA, IFINRA
      INTEGER          IEL, ICHA, ICLA, IPHAS, IPCROM, IPCRO2
      INTEGER          IZONE, IFAC
      INTEGER          IPBROM, IPCX2C, IROMF
      DOUBLE PRECISION X1SRO1, X2SRO2, SRROM1, UNS1PW
      DOUBLE PRECISION X2TOT, WMOLME, UNSRO1
C
      INTEGER       IPASS
      DATA          IPASS /0/
      SAVE          IPASS
C
C***********************************************************************
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
c --- Initialisation des tableaux de travail
C
      DO IEL = 1, NCEL
        W1(IEL) = ZERO
        W2(IEL) = ZERO
        W3(IEL) = ZERO
        W4(IEL) = ZERO
        W5(IEL) = ZERO
        W6(IEL) = ZERO
        W7(IEL) = ZERO
        W8(IEL) = ZERO
      ENDDO
C
C     Pointeur sur masse volumique du gaz aux cellules
      IROMF = IPPROC(IROM1)
C
C=======================================================================
C 2. CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE DISPERSEE
C                    VALEURS CELLULES
C                    ----------------
C    FRACTION MASSIQUE DE SOLIDE
C    DIAMETRE
C    MASSE VOLUMIQUE
C=======================================================================
C
      CALL CPPHY2
C     ===========
     & ( NCELET , NCEL   ,
     &   RTP    , PROPCE )
C
C=======================================================================
C 3. CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE GAZEUSE
C                    VALEURS CELLULES
C                    ----------------
C    TEMPERATURE
C    MASSE VOLUMIQUE
C    CONCENTRATIONS DES ESPECES GAZEUSES
C=======================================================================
C
C --- Calcul de l'enthalpie du gaz     dans W8
C            de F1M                    dans W2
C            de F2M                    dans W3
C            de F3M                    dans W4
C            de F4M                    dans W5
C            de F5M                    dans W6
C            de F4P2M                  dans W7
C
C ---- W1 = - Somme des X2(i)
C
      DO ICLA = 1, NCLACP
        IPCX2C = IPPROC(IX2(ICLA))
        DO IEL = 1, NCEL
          W1(IEL) =  W1(IEL) - PROPCE(IEL,IPCX2C)
        ENDDO
      ENDDO
C
C
C ---- W2 = F1M = SOMME(F1M(ICHA))
C      W3 = F2M = SOMME(F2M(ICHA))
C      W4 = F3M
C      W6 = F5M
C      W5 = F4M = 1. - F1M - F2M - F3M - F5M
C      W7 = F4P2M
C
      DO ICHA = 1, NCHARB
        DO IEL = 1, NCEL
          W2(IEL) =  W2(IEL) + RTP(IEL,ISCA(IF1M(ICHA)))
          W3(IEL) =  W3(IEL) + RTP(IEL,ISCA(IF2M(ICHA)))
        ENDDO
      ENDDO
C
      DO IEL = 1, NCEL
        UNS1PW = 1.D0/(1.D0+W1(IEL))
        W2(IEL) =  W2(IEL)              *UNS1PW
        W3(IEL) =  W3(IEL)              *UNS1PW
        W4(IEL) =  RTP(IEL,ISCA(IF3M))  *UNS1PW
        IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
          W6(IEL) =  RTP(IEL,ISCA(IF5M))  *UNS1PW
        ELSE
          W6(IEL) = 0.D0
        ENDIF
        W5(IEL) = 1.D0 - W2(IEL) - W3(IEL) -W4(IEL)-W6(IEL)
        W7(IEL) =  RTP(IEL,ISCA(IF4P2M))*UNS1PW
      ENDDO
C
C ---- W8 = H1 (transport de H2)
C        Transport d'H2
C
      DO ICLA = 1, NCLACP
        DO IEL = 1, NCEL
          W8(IEL) =  W8(IEL) - RTP(IEL,ISCA(IH2(ICLA)))
        ENDDO
      ENDDO
      DO IEL = 1, NCEL
        W8(IEL) = (RTP(IEL,ISCA(IHM))+W8(IEL))/ ( 1.D0+W1(IEL) )
      ENDDO
C
C --- Gestion memoire
C     Autres tableaux
C
C ------ Macro tableau d'entiers TBCPI : NTBCPI
C        Macro tableau de reels  TBCPR : NTBCPR
C        Macro tableau d'entiers TBMCI : NTBMCI
C        Macro tableau de reels  TBMCR : NTBMCR
C        Macro tableau d'entiers TBWOI : NTBWOI
C        Macro tableau de reels  TBWOR : NTBWOR
C
      NTBCPI = 1
      NTBCPR = 15
      NTBMCI = 0
      NTBMCR = 2 + 2*NCHARB + 4
C  Ce sont en fait X1M, X2M,
C                  F1M(ICHA) et F2M(ICHA) pour chaque charbon
C                  ACHX1F1, ACHX2F2, ACOF1, ACOF2
      NTBWOI = 1
      NTBWOR = 4
C
      CALL MEMCP1
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NVAR   , NCELET , NCEL   , NFAC   , NFABOR ,
     &   NTBCPI , ICPWI  ,
     &   NTBCPR , ICPWR  ,
     &   NTBMCI , IMCWI  ,
     &   NTBMCR , IMCWR  ,
     &   NTBWOI , IWORI  ,
     &   NTBWOR , IWORR  ,
     &   IFINIA , IFINRA )
C
      CALL CPPHY1
C     ===========
     & ( IFINIA , IFINRA ,
     &   NCELET , NCEL   ,
     &   NTBCPI , NTBCPR , NTBMCI , NTBMCR , NTBWOI , NTBWOR ,
     &   W2     , W3     , W4     , W5     , W6    , W7     ,
C         F1M      F2M      F3M      F4M     F5M    F4P2M
     &   W8     ,
C         ENTH
     &   RTP    , PROPCE  , PROPCE(1,IROMF) ,
C                          ---------------- (masse vol. gaz)
     &   IA(ICPWI) , RA(ICPWR) ,
     &   IA(IMCWI) , RA(IMCWR) ,
     &   IA(IWORI) , RA(IWORR)  )
C
C
C=======================================================================
C 4. CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE DISPERSEE
C                    VALEURS CELLULES
C                    ----------------
C    TEMPERATURE
C=======================================================================
C
C --- Transport d'H2
C
      CALL CPTEH2
C     ===========
     & ( NCELET , NCEL   ,
     &   RTP    , PROPCE ,
     &   W3     , W4     )
C
C=======================================================================
C 5. CALCUL DES PROPRIETES PHYSIQUES DU MELANGE
C                    VALEURS CELLULES
C                    ----------------
C    MASSE VOLUMIQUE
C=======================================================================
C
C --- W2 = - Somme des X2(i)
C
      DO IEL = 1, NCEL
        W2(IEL) = ZERO
      ENDDO
C
      DO ICLA = 1, NCLACP
        IPCX2C = IPPROC(IX2(ICLA))
        DO IEL = 1, NCEL
          W2(IEL) =  W2(IEL) - PROPCE(IEL,IPCX2C)
        ENDDO
      ENDDO
C
C --- Calcul de Rho du melange : 1/Rho = X1/Rho1 + Somme(X2/Rho2)
C     On sous relaxe quand on a un rho n a disposition, ie
C       a partir du deuxieme passage ou
C       a partir du premier passage si on est en suite de calcul et
C         qu'on a relu la masse volumique dans le fichier suite.
C
      IPHAS = 1
      IPCROM = IPPROC(IROM(IPHAS))
C
      IF (IPASS.GT.1.OR.(ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1)) THEN
        SRROM1 = SRROM
      ELSE
        SRROM1 = 1.D0
      ENDIF
C
      DO IEL = 1, NCEL
        X2SRO2 = ZERO
        DO ICLA = 1, NCLACP
          IPCRO2 = IPPROC(IROM2(ICLA))
          IPCX2C = IPPROC(IX2(ICLA))
          X2SRO2 = X2SRO2 + PROPCE(IEL,IPCX2C) / PROPCE(IEL,IPCRO2)
        ENDDO
        X1SRO1 = (1.D0+W2(IEL)) / PROPCE(IEL,IROMF)
C ---- Sous relaxation eventuelle a donner dans ppini1.F
        PROPCE(IEL,IPCROM) = SRROM1*PROPCE(IEL,IPCROM)
     &                     + (1.D0-SRROM1)/(X1SRO1+X2SRO2)
      ENDDO
C
C
C=======================================================================
C 6. CALCUL DE RHO DU MELANGE
C
C                      VALEURS FACETTES
C                      ----------------
C=======================================================================
C
      IPHAS = 1
      IBROM(IPHAS) = 1
      IPBROM = IPPROB(IROM(IPHAS))
      IPCROM = IPPROC(IROM(IPHAS))
C
C ---> Masse volumique au bord pour toutes les facettes
C      Les facettes d'entree seront recalculees.
C
      DO IFAC = 1, NFABOR
        IEL = IFABOR(IFAC)
        PROPFB(IFAC,IPBROM) = PROPCE(IEL,IPCROM)
      ENDDO
C
C ---> Masse volumique au bord pour les facettes d'entree UNIQUEMENT
C     Le test sur IZONE sert pour les reprises de calcul
C
      IF ( IPASS.GT.1 .OR. ISUITE.EQ.1 ) THEN
        DO IFAC = 1, NFABOR
C
          IZONE = IZFPPP(IFAC)
          IF(IZONE.GT.0) THEN
            IF ( IENTAT(IZONE).EQ.1 .OR. IENTCP(IZONE).EQ.1 ) THEN
              X2SRO2 = ZERO
              X2TOT  = ZERO
              DO ICLA = 1, NCLACP
                X2SRO2 = X2SRO2 + X20(IZONE,ICLA)/RHO20(ICLA)
                X2TOT  = X2TOT  + X20(IZONE,ICLA)
              ENDDO
              WMOLME = (1.D0+XSI) / (WMOLE(IO2)+XSI*WMOLE(IN2))
              UNSRO1 = (WMOLME*RR*TIMPAT(IZONE)) / P0(IPHAS)
              X1SRO1 = (1.D0-X2TOT) * UNSRO1
              PROPFB(IFAC,IPBROM) = 1.D0 / (X1SRO1+X2SRO2)
            ENDIF
          ENDIF
C
        ENDDO
      ENDIF
C
C----
C FIN
C----
C
      RETURN
      END
c@z
