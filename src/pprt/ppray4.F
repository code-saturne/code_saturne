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
                        SUBROUTINE PPRAY4
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
C
     &   MODE   ,
C
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER ,
C
     &   TPAROP , HPAROP , TEMPK  ,
C
     &   RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC
CFONC    SOUS-PROGRAMME PHYSIQUES PARTICULIERES
CFONC    POUR LES CONVERSIONS TEMPERATURE <-> ENTHALPIE
CFONC       POUR LE MELANGE GAZEUX
CFONC
CFONC
CFONC
CFONC      A) L'ENTHALPIE DOIT ETRE CONVERTIE EN TEMPERATURE EN KELVIN
CFONC         IL FAUT ECRIRE UNE LOI DE CONVERSION :
CFONC         ENTHALPIE   -> TEMPERATURE (MODE =  1)
CFONC         POUR REMPLIR TEMPK
CFONC
CFONC
CFONC      B) SI DE PLUS ON CALCULE LES TEMPERATURES DE PAROI
CFONC         ALORS IL FAUT FOURNIR UNE LOI DE CONVERSION :
CFONC         TEMPERATURE -> ENTHALPIE   (MODE = -1)
CFONC         POUR REMPLIR HPAROP
CFONC
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
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE LA PHASE COURANTE          !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! MODE         ! E  !  -> ! TYPE DE CONVERSION ENTHAL<->TEMPK    !
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
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
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
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COFRUA,COFRUB! TR ! <-  ! CONDITIONS AUX LIMITES AUX           !
CARGU !(NFABOR)      !    !     !    FACES DE BORD POUR LA LUMINANCES  !
CARGU ! W1...6(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! TEMPK(NCELET)! TR ! <-  ! TEMPERATURE EN KELVIN                !
CARGU ! HPAROP(NFABOR! TR ! <-  ! ENTHALPIE MASSIQUE DE PAROI EN J/Kg  !
CARGU !              !    !     ! (EN DEGRES CELSIUS OU KELVIN)        !
CARGU ! TPAROP(NFABOR! TR !  -> ! TEMPERATURE DE PAROI EN KELVIN       !
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
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "radiat.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "fuincl.h"
      INCLUDE "ppincl.h"
C
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , IPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          MODE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML) , ITYPFB(NFABOR)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
C
      DOUBLE PRECISION TEMPK(NCELET)
      DOUBLE PRECISION TPAROP(NFABOR), HPAROP(NFABOR)
C
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IEL , IFAC , ICLA , ICHA , ISOL , IGE
      INTEGER          IPCX2C , IXCHCL, IXCKCL, IXNPCL, IGG, III
      INTEGER          IESP
      DOUBLE PRECISION COEFG(NGAZGM), COEFE(NGAZEM)
      DOUBLE PRECISION F1MC(NCHARM) , F2MC(NCHARM)
      DOUBLE PRECISION X2T , H2 , X2H2 , HF , XSOLID(NSOLIM),T1
      DOUBLE PRECISION YM    (NGAZGM)
      DOUBLE PRECISION DIAMGT,MASGUT,MKGOUT,MFGOUT,MKFINI,RHOFOL
C
C***********************************************************************
C
C=======================================================================
C 0 - GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1 - INITIALISATIONS GENERALES
C=======================================================================
C
C
C=======================================================================
C  2.1 - CALCUL DE LA TEMPERATURE EN KELVIN AUX CELLULES
C=======================================================================
C
C---> CONVERSION ENTHALPIE -> TEMPERATURE (MODE =  1)
C     AUX CELLULES FLUIDES
C     (UTILE SI NON TRANSPARENT)
C     -----------------------------------------------
C
C
      IF (MODE.EQ.1) THEN
C
C ---- Combustion gaz : Flamme de Premelange ou Flamme de Diffusion
C
        IF ( IPPMOD(ICOEBU).GE.0 .OR.
     &       IPPMOD(ICOD3P).GE.0      ) THEN
C
          DO IEL = 1,NCEL
            TEMPK(IEL) = PROPCE(IEL, IPPROC(ITEMP))
          ENDDO
C
C ---- Combustion charbon pulverise
C
        ELSE IF ( IPPMOD(ICP3PL).GE.0 ) THEN
C
          DO IEL = 1,NCEL
            TEMPK(IEL) = PROPCE(IEL, IPPROC(ITEMP1))
          ENDDO
C
C ---- Combustion fuel
C
        ELSE IF ( IPPMOD(ICFUEL).GE.0 ) THEN
C
          DO IEL = 1,NCEL
            TEMPK(IEL) = PROPCE(IEL, IPPROC(ITEMP1))
          ENDDO
C
C ---- Module Electrique
C
        ELSE IF ( IPPMOD(IELJOU).GE.1 .OR.
     &            IPPMOD(IELARC).GE.1 .OR.
     &            IPPMOD(IELION).GE.1       ) THEN
C
          DO IEL = 1,NCEL
            TEMPK(IEL) = PROPCE(IEL, IPPROC(ITEMP))
          ENDDO
C
        ENDIF
C
      ENDIF
C
C
C=======================================================================
C  2.2 - CALCUL DE L'ENTHALPIE AU FACES DE BORD
C=======================================================================
C
C---> CONVERSION TEMPERATURE -> ENTHALPIE (MODE = -1)
C     AUX FACES DE BORD DE PAROI
C     (INUTILE POUR LES PAROIS ISOTH(IFAC)=3 ET EPS(IFAC)=0)
C     ------------------------------------------------------
C
C
      IF (MODE.EQ.-1) THEN
C
        DO IFAC = 1,NFABOR
C
          IF (ITYPFB(IFAC).EQ.IPAROI .OR.
     &        ITYPFB(IFAC).EQ.IPARUG) THEN
C
C   Numero de la cellule en regard
C
            IEL = IFABOR(IFAC)
C
C ---- Combustion gaz : Flamme de Premelange ou Flamme de Diffusion
C
            IF ( IPPMOD(ICOEBU).GE.0 .OR.
     &           IPPMOD(ICOD3P).GE.0      ) THEN
C
              DO IGG = 1, NGAZGM
                COEFG(IGG) = ZERO
              ENDDO
              COEFG(1) = PROPFB(IFAC,IPPROB(IYM(1)))
              COEFG(2) = PROPFB(IFAC,IPPROB(IYM(2)))
              COEFG(3) = PROPFB(IFAC,IPPROB(IYM(3)))
              MODE     = -1
              CALL COTHHT
C             ===========
     &        ( MODE   , NGAZG , NGAZGM  , COEFG  ,
     &          NPO    , NPOT   , TH     , EHGAZG ,
     &          HPAROP(IFAC) , TPAROP(IFAC) )
C
C ---- Combustion charbon pulverise
C
            ELSE IF ( IPPMOD(ICP3PL).GE.0 ) THEN
C
              X2T  = ZERO
              X2H2 = ZERO
              DO ICLA = 1, NCLACP
                ICHA   = ICHCOR(ICLA)
                IXCHCL = ISCA(IXCH(ICLA))
                IXCKCL = ISCA(IXCK(ICLA))
                IXNPCL = ISCA(INP(ICLA ))
                IPCX2C = IPPROC(IX2(ICLA))
                X2T = X2T + PROPCE(IEL,IPCX2C)
                H2 = ZERO
                DO ISOL = 1, NSOLIM
                  XSOLID(ISOL) = ZERO
                ENDDO
                IF (PROPCE(IEL,IPCX2C).GT.EPSICP) THEN
                  XSOLID(ICH(ICHA) ) = RTP(IEL,IXCHCL)
     &                               / PROPCE(IEL,IPCX2C)
                  XSOLID(ICK(ICHA) ) = RTP(IEL,IXCKCL)
     &                               / PROPCE(IEL,IPCX2C)
                  XSOLID(IASH(ICHA)) = RTP(IEL,IXNPCL)*XMASH(ICLA)
     &                               / PROPCE(IEL,IPCX2C)
                  IF ( IPPMOD(ICP3PL).EQ.1 ) THEN
                    XSOLID(IWAT(ICHA)) = RTP(IEL,ISCA(IXWT(ICLA)))
     &                                  /PROPCE(IEL,IPCX2C)
                  ENDIF
                  III = ICLA
                  T1 = TPAROP(IFAC)
                  CALL CPTHP2
C                 ===========
     &            ( MODE   , III   , H2     , XSOLID ,
     &              TPAROP(IFAC) , T1 )
                ENDIF
                X2H2 = X2H2 + PROPCE(IEL,IPCX2C)*H2
              ENDDO
C
              DO IGE = 1, NGAZEM
                COEFE(IGE) = ZERO
              ENDDO
              DO ICHA = 1, NCHARB
                F1MC(ICHA) = RTP(IEL,ISCA(IF1M(ICHA)))
     &                     / (1.D0-X2T)
                F2MC(ICHA) =  RTP(IEL,ISCA(IF2M(ICHA)))
     &                     / (1.D0-X2T)
              ENDDO
              DO ICHA = (NCHARB+1), NCHARM
                F1MC(ICHA) = ZERO
                F2MC(ICHA) = ZERO
              ENDDO
              DO IGE = 1, (NGAZE-2*NCHARB)
                COEFE(IGE) = PROPCE(IEL,IPPROC(IYM1(IGE)))
              ENDDO
              CALL CPTHP1
C             ===========
     &        ( MODE   , HF     , COEFE       ,
     &          F1MC   , F2MC   ,TPAROP(IFAC) )
                HPAROP(IFAC) = (1.D0-X2T)*HF+X2H2
C
C ---- Combustion fuel
C
            ELSE IF ( IPPMOD(ICFUEL).GE.0 ) THEN
C
              X2T  = ZERO
              X2H2 = ZERO
C
              DO ICLA=1,NCLAFU
C
                X2T = X2T+RTPA(IEL,ISCA(IYFOL(ICLA)))
C
                MKFINI = RHO0FL*PI/6.D0*DINIKF(ICLA)**3
                RHOFOL = PROPCE(IEL,IPPROC(IROM3(ICLA)))
                DIAMGT = PROPCE(IEL,IPPROC(IDIAM3(ICLA)))
                MASGUT = RHOFOL*PI/6.D0*DIAMGT**3
                IF (DIAMGT.LE.DINIKF(ICLA)) THEN
                  MKGOUT = MASGUT
                ELSE
                  MKGOUT = MKFINI
                ENDIF
                MFGOUT = MASGUT - MKGOUT
                XSOLID(1) = 1.D0-FKC
                XSOLID(2) = FKC
                IF(MASGUT.GT.EPZERO) THEN
                  XSOLID(1) = MFGOUT / MASGUT
                  XSOLID(2) = MKGOUT / MASGUT
                ENDIF
                XSOLID(1) = MIN(1.D0,MAX(0.D0,XSOLID(1)))
                XSOLID(2) = MIN(1.D0,MAX(0.D0,XSOLID(2)))
C
                CALL FUTHP2
C               ===========
     &        ( MODE   , H2     , XSOLID , TPAROP(IFAC) )
C
                X2H2 =  X2H2 + RTPA(IEL,ISCA(IYFOL(ICLA)))*H2
C
              ENDDO
C
              DO IGE = 1,NGAZE
                COEFE(IGE) = PROPCE(IEL,IPPROC(IYM1(IGE)))
              ENDDO
              CALL FUTHP1
C             ===========
     &        ( MODE   , HF     , COEFE  , TPAROP(IFAC) )
                HPAROP(IFAC) = (1.D0-X2T)*HF+X2H2
C
C ---- Module Electrique
C
C
            ELSE IF ( IPPMOD(IELJOU).GE.1  ) THEN
C
              CALL USTHHT (MODE,HPAROP(IFAC),TPAROP(IFAC))
C
            ELSE IF ( IPPMOD(IELARC).GE.1  ) THEN
C
              IF ( NGAZG .EQ. 1 ) THEN
                YM(1) = 1.D0
                CALL ELTHHT(MODE,NGAZG,YM,HPAROP(IFAC),TPAROP(IFAC))
              ELSE
                YM(NGAZG) = 1.D0
                DO IESP = 1, NGAZG-1
                  YM(IESP) = RTP(IEL,ISCA(IYCOEL(IESP)))
                  YM(NGAZG) = YM(NGAZG) - YM(IESP)
                ENDDO
                CALL ELTHHT(MODE,NGAZG,YM,HPAROP(IFAC),TPAROP(IFAC))
              ENDIF
C
C     ELSE IF ( IPPMOD(IELION).GE.1 ) THEN
C     ... to be implemented ...
C
            ENDIF
C
          ENDIF
C
        ENDDO
C
      ENDIF
C
C----
C FIN
C----
C
      RETURN
      END
c@z
