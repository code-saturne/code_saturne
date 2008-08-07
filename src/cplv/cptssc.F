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
                        SUBROUTINE CPTSSC
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCAL  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IZFPPP , IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL ,
     &   SMBRS  , ROVSDT ,
     &   VISCF  , VISCB  , XAM    ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     , W10    , W11    ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC ROUTINE PHYSIQUE PARTICULIERE : FLAMME CHARBON PULVERISE
CFONC    ON PRECISE LES TERMES SOURCES POUR UN SCALAIRE PP
CFONC    SUR UN PAS DE TEMPS
CFONC
CFONC ATTENTION : LE TRAITEMENT DES TERMES SOURCES EST DIFFERENT
CFONC ---------   DE CELUI DE USTSSC.F
CFONC
CFONC ON RESOUT ROVSDT*D(VAR) = SMBRS
CFONC
CFONC ROVSDT ET SMBRS CONTIENNENT DEJA D'EVENTUELS TERMES SOURCES
CFONC   UTILISATEUR. IL FAUT DONC LES INCREMENTER ET PAS LES
CFONC   ECRASER
CFONC
CFONC POUR DES QUESTIONS DE STABILITE, ON NE RAJOUTE DANS ROVSDT
CFONC   QUE DES TERMES POSITIFS. IL N'Y A PAS DE CONTRAINTE POUR
CFONC   SMBRS
CFONC
CFONC DANS LE CAS D'UN TERME SOURCE EN CEXP + CIMP*VAR ON DOIT
CFONC ECRIRE :
CFONC           SMBRS  = SMBRS  + CEXP + CIMP*VAR
CFONC           ROVSDT = ROVSDT + MAX(-CIMP,ZERO)
CFONC
CFONC ON FOURNIT ICI ROVSDT ET SMBRS (ILS CONTIENNENT RHO*VOLUME)
CFONC     SMBRS en kg variable/s :
CFONC      ex : pour la vitesse            kg m/s2
CFONC           pour les temperatures      kg degres/s
CFONC           pour les enthalpies        Joules/s
CFONC     ROVSDT en kg /s
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
CARGU ! NCESMP       ! E  !  -> ! NOMBRE DE CELLULES A SOURCE DE MASSE !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! ISCAL        ! E  !  -> ! NUMERO DU SCALAIRE                   !
CARGU ! IFACEL       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMFBR       ! TE !  -> ! NUMERO DE FAMILLE D'UNE FACE DE BORD !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMCEL       ! TE !  -> ! NUMERO DE FAMILLE D'UNE CELLULE      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! IPRFML       ! TE !  -> ! PROPRIETES D'UNE FAMILLE             !
CARGU ! ITYPFB(NFABOR! TE ! <-  ! TYPE DES FACES DE BORD               !
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! ICEPDC(NCELET! TE !  -> ! NUMERO DES NCEPDP CELLULES AVEC PDC  !
CARGU ! ICETSM(NCESMP! TE !  -> ! NUMERO DES CELLULES A SOURCE DE MASSE!
CARGU ! ITYPSM       ! TE !  -> ! TYPE DE SOURCE DE MASSE POUR LES     !
CARGU ! (NCESMP,NVAR)!    !     !  VARIABLES (cf. USTSMA)              !
CARGU ! IZFPPP       ! TE ! <-  ! NUMERO DE ZONE DE LA FACE DE BORD    !
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
CARGU ! SMACEL       ! TR !  -> ! VALEUR DES VARIABLES ASSOCIEE A LA   !
CARGU ! (NCESMP,*   )!    !     !  SOURCE DE MASSE                     !
CARGU !              !    !     !  POUR IVAR=IPR, SMACEL=FLUX DE MASSE !
CARGU ! SMBRS(NCELET)! TR ! <-  ! SECOND MEMBRE EXPLICITE              !
CARGU ! ROVSDT(NCELET! TR ! <-  ! PARTIE DIAGONALE IMPLICITE           !
CARGU ! VISCF(NFAC)  ! TR !  -  ! TABLEAU DE TRAVAIL    FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! TABLEAU DE TRAVAIL    FACES DE BORD  !
CARGU ! XAM(NFAC,2)  ! TR !  -  ! TABLEAU DE TRAVAIL    FACES DE BORD  !
CARGU ! W1..11(NCELET! TR !  -  ! TABLEAU DE TRAVAIL    CELLULES       !
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
      INCLUDE "dimfbr.h"
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
      INCLUDE "ppcpfu.h"
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
      INTEGER          ISCAL
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML) , ITYPFB(NFABOR,NPHAS)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICEPDC(NCEPDP)
      INTEGER          ICETSM(NCESMP), ITYPSM(NCESMP,NVAR)
      INTEGER          IZFPPP(NFABOR)
      INTEGER          IDEVEL(NIDEVE)
      INTEGER          ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NDIMFB,*)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP), SMACEL(NCESMP,NVAR)
      DOUBLE PRECISION SMBRS(NCELET), ROVSDT(NCELET)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION XAM(NFAC,2)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
      DOUBLE PRECISION W7(NCELET), W8(NCELET), W9(NCELET)
      DOUBLE PRECISION W10(NCELET), W11(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IVAR , IPCROM, IEL, IPHAS
      INTEGER          NUMCLA , NUMCHA , ICLA , NUMTRA
      INTEGER          IPCGCH , IPCGD1 , IPCGD2 , IPCGHT , IPCSEC
      INTEGER          IXCHCL , IXCKCL , ISCALA
      INTEGER          IPCRO2 , IPCTE1 , IPCTE2 , IPCVSL , IPCCP
      INTEGER          IPCDIA , IPCVST
      INTEGER          MODE, IGE
      INTEGER          IFAC   , IFINRA , ICOEFA , ICOEFB
      INTEGER          IPCX2C , ICHA , IMODE , II
      INTEGER          ITERCH
      INTEGER          ITERMX,NBPAUV,NBRICH,NBEPAU,NBERIC,IPGHC2
C
      DOUBLE PRECISION EPSRGP , CLIMGP , EXTRAP , XNUSS
      DOUBLE PRECISION AUX, RHOVST
      DOUBLE PRECISION COEFE(NGAZEM)
      DOUBLE PRECISION T1, T2, HLV , HH2OV
      DOUBLE PRECISION F1MC(NCHARM), F2MC(NCHARM)
      DOUBLE PRECISION XHDEV1 , XHDEV2 , XHCO , XHO2 , GAMDV1 , GAMDV2
C
      DOUBLE PRECISION GAMHET
C
      DOUBLE PRECISION XXCO,XXO2,XXCO2,XXH2O,XCO2MX
      DOUBLE PRECISION XKP,XK0P,XKM,XK0M,WPLUS,WMOINS,T0P,T0M
      DOUBLE PRECISION AUXP,AUXM
C
      DOUBLE PRECISION ANMR,XCOT,XO2T,XCO2E,XO2E,XCOE,TAUCHI,TAUTUR
      DOUBLE PRECISION SQH2O , X2
C
      DOUBLE PRECISION ERR1MX,ERR2MX,ANMR0
C
      DOUBLE PRECISION ERRCH,TER1,DDELTA,XDEN
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C --- Numero du scalaire a traiter : ISCAL
C
C --- Numero de la variable associee au scalaire a traiter ISCAL
      IVAR = ISCA(ISCAL)
C
C --- Nom de la variable associee au scalaire a traiter ISCAL
      CHAINE = NOMVAR(IPPRTP(IVAR))
C
C --- Numero de phase associee au scalaire ISCAL
      IPHAS = IPHSCA(ISCAL)
C
C --- Numero des grandeurs physiques (voir usclim)
      IPCROM = IPPROC(IROM(IPHAS))
      IPCVST = IPPROC(IVISCT(IPHAS))
C
C=======================================================================
C 2. PRISE EN COMPTE DES TERMES SOURCES POUR LES VARIABLES RELATIVES
C    AUX CLASSES DE PARTICULES
C=======================================================================
C
C --> Terme source pour la fraction massique de charbon reactif
C
      IF ( IVAR.GE.ISCA(IXCH(1)) .AND. IVAR.LE.ISCA(IXCH(NCLACP)) ) THEN
C
        IF (IWARNI(IVAR).GE.1) THEN
          WRITE(NFECRA,1000) CHAINE(1:8)
        ENDIF
        NUMCLA = IVAR-ISCA(IXCH(1))+1
        IPCGCH = IPPROC(IGMDCH(NUMCLA))
C
        DO IEL = 1, NCEL
C
C ---- Calcul de  W1 = - rho.GMDCH > 0
C
          W1(IEL) = - PROPCE(IEL,IPCROM)*PROPCE(IEL,IPCGCH)
     &               *VOLUME(IEL)
C
C ---- Calcul des parties explicite et implicite du TS
C
          ROVSDT(IEL) = ROVSDT(IEL) + MAX(W1(IEL),ZERO)
          SMBRS(IEL) = SMBRS(IEL) - W1(IEL)*RTPA(IEL,IVAR)
C
        ENDDO
C
      ENDIF
C
C
C --> Terme source pour la fraction massique de coke
C
      IF ( IVAR.GE.ISCA(IXCK(1)) .AND. IVAR.LE.ISCA(IXCK(NCLACP)) ) THEN
C
        IF (IWARNI(IVAR).GE.1) THEN
          WRITE(NFECRA,1000) CHAINE(1:8)
        ENDIF
C
        NUMCLA = IVAR-ISCA(IXCK(1))+1
        IPCGCH = IPPROC(IGMDCH(NUMCLA))
        IPCGD1 = IPPROC(IGMDV1(NUMCLA))
        IPCGD2 = IPPROC(IGMDV2(NUMCLA))
        IXCHCL = ISCA(IXCH(NUMCLA))
        IXCKCL = ISCA(IXCK(NUMCLA))
        IPCGHT = IPPROC(IGMHET(NUMCLA))
        IPGHC2 = IPPROC(IGHCO2(NUMCLA))
C
        DO IEL = 1, NCEL
C
C ---- Calcul de W1 = - rho.Xch.GMDCH.Volume > 0
C
          W1(IEL) = - PROPCE(IEL,IPCROM)*RTP(IEL,IXCHCL)
     &                *PROPCE(IEL,IPCGCH)*VOLUME(IEL)
C
C AE : On prend RTP(IEL,IXCHCL) et pas RTPA(IEL,IXCHCL) afin
C      d'etre conservatif sur la masse
C
C ---- Calcul de W2 = rho.Xch.(GMDV1+GMDV2)Volume < 0
C
          W2(IEL) = PROPCE(IEL,IPCROM)*RTP(IEL,IXCHCL)
     &              *(PROPCE(IEL,IPCGD1)+PROPCE(IEL,IPCGD2))
     &              *VOLUME(IEL)
C
          IF ( RTPA(IEL,IXCKCL) .GT. EPSICP ) THEN
C
C Reaction C(s) + O2 ---> 0.5CO
C =============================
C
C ---- Calcul de la partie implicite  > 0 du TS relatif a GMHET

            W3(IEL) = - 2.D0/3.D0*PROPCE(IEL,IPCROM)*PROPCE(IEL,IPCGHT)
     &               /(RTPA(IEL,IXCKCL))**(1.D0/3.D0)*VOLUME(IEL)
C
C ---- Calcul de la partie explicite < 0 du TS relatif a GMHET
C
            W4(IEL) = PROPCE(IEL,IPCROM)*PROPCE(IEL,IPCGHT)
     &                * (RTPA(IEL,IXCKCL))**(2.D0/3.D0)*VOLUME(IEL)
C
          ELSE
            W3(IEL) = 0.D0
            W4(IEL) = 0.D0
          ENDIF
C
C ---- Calcul des parties explicite et implicite du TS
C
          ROVSDT(IEL) = ROVSDT(IEL) + MAX(W3(IEL),ZERO)
     &                              + MAX(W5(IEL),ZERO)
          SMBRS(IEL)  = SMBRS(IEL)  + W1(IEL) + W2(IEL) + W4(IEL)
     &                              + W6(IEL)
C
        ENDDO
C
        IF ( IHTCO2 .EQ. 1 ) THEN
C
          DO IEL = 1, NCEL
C
            IF ( RTPA(IEL,IXCKCL) .GT. EPSICP ) THEN
C
C Reaction C(s) + CO2 ---> 2CO
C =============================
C
C ---- Calcul de la partie implicite  > 0 du TS relatif a GMHET

              W5(IEL) = - 2.D0/3.D0*PROPCE(IEL,IPCROM)
     &                             *PROPCE(IEL,IPGHC2)
     &                 /(RTPA(IEL,IXCKCL))**(1.D0/3.D0)*VOLUME(IEL)
C
C ---- Calcul de la partie explicite < 0 du TS relatif a GMHET
C
              W6(IEL) = PROPCE(IEL,IPCROM)*PROPCE(IEL,IPGHC2)
     &                 *(RTPA(IEL,IXCKCL))**(2.D0/3.D0)*VOLUME(IEL)
C
            ELSE
              W5(IEL) = 0.D0
              W6(IEL) = 0.D0
            ENDIF
C
C ---- Calcul des parties explicite et implicite du TS
C
            ROVSDT(IEL) = ROVSDT(IEL) + MAX(W5(IEL),ZERO)
            SMBRS(IEL)  = SMBRS(IEL)  + W6(IEL)
C
          ENDDO
C
        ENDIF
C
      ENDIF
C
C
C --> Terme source pour la fraction massique de d'eau
C
      IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
C
        IF ( IVAR.GE.ISCA(IXWT(1)) .AND.
     &       IVAR.LE.ISCA(IXWT(NCLACP)) ) THEN
C
          IF (IWARNI(IVAR).GE.1) THEN
            WRITE(NFECRA,1000) CHAINE(1:8)
          ENDIF
          NUMCLA = IVAR-ISCA(IXWT(1))+1
          NUMCHA = ICHCOR(NUMCLA)
C
          IPCSEC = IPPROC(IGMSEC(NUMCLA))
          IPCX2C = IPPROC(IX2(NUMCLA))
C
          DO IEL = 1, NCEL
C
C ---- Calcul des parties explicite et implicite du TS
C
           IF ( RTPA(IEL,IVAR).GT. EPSICP .AND.
     &          XWATCH(NUMCHA).GT. EPSICP       ) THEN
             W1(IEL) = PROPCE(IEL,IPCROM)*PROPCE(IEL,IPCSEC)
     &                     *VOLUME(IEL)
     &                     *(1.D0/PROPCE(IEL,IPCX2C))
     &                     *(1.D0/XWATCH(NUMCHA))
C
              ROVSDT(IEL) = ROVSDT(IEL) + MAX(W1(IEL),ZERO)
              SMBRS(IEL)  = SMBRS(IEL)  - W1(IEL)*RTPA(IEL,IVAR)
           ENDIF
C
          ENDDO
C
        ENDIF
C
      ENDIF
C
C --> Terme source pour l'enthalpie du solide
C
      IF ( IVAR.GE.ISCA(IH2(1)) .AND. IVAR.LE.ISCA(IH2(NCLACP)) ) THEN
C
        IF (IWARNI(IVAR).GE.1) THEN
          WRITE(NFECRA,1000) CHAINE(1:8)
        ENDIF
C
        NUMCLA = IVAR-ISCA(IH2(1))+1
        NUMCHA = ICHCOR(NUMCLA)
        IXCHCL = ISCA(IXCH(NUMCLA))
        IXCKCL = ISCA(IXCK(NUMCLA))
        IPCX2C = IPPROC(IX2(NUMCLA))
        IPCRO2 = IPPROC(IROM2(NUMCLA ))
        IPCDIA = IPPROC(IDIAM2(NUMCLA))
        IPCTE2 = IPPROC(ITEMP2(NUMCLA))
        IPCTE1 = IPPROC(ITEMP1)
        IPCGHT = IPPROC(IGMHET(NUMCLA))
        IPGHC2 = IPPROC(IGMHET(NUMCLA))
C
        IPCGD1 = IPPROC(IGMDV1(NUMCLA))
        IPCGD2 = IPPROC(IGMDV2(NUMCLA))
        IPCGCH = IPPROC(IGMDCH(NUMCLA))
        IXCHCL = ISCA(IXCH(NUMCLA))
        IXCKCL = ISCA(IXCK(NUMCLA))
C
C ---- Calcul preliminaire : calcul de X2(NUMCLA) dans W11
C      ATTENTION tableau a conserver
C
        DO IEL = 1, NCEL
C Rq : on utilise PROPCE(IEL,IPCX2C) car on veut X2 a l'iteration n
C      (en RTPA)
          W11(IEL) = PROPCE(IEL,IPCX2C)
        ENDDO
C
C ---- Contribution aux bilans explicite et implicite
C        des echanges par diffusion moleculaire
C        6 Lambda Nu / diam**2 / Rho2 * Rho * (T1-T2)
C
C ------ Calcul de lambda dans W1
C
        XNUSS = 2.D0
        DO IEL = 1, NCEL
          IF ( IVISLS(IHM).GT.0 ) THEN
            IPCVSL = IPPROC(IVISLS(IHM))
            IF ( ICP(IPHAS).GT.0 ) THEN
              IPCCP   = IPPROC(ICP(IPHAS))
              W1(IEL) = PROPCE(IEL,IPCVSL) * PROPCE(IEL,IPCCP)
            ELSE
              W1(IEL) = PROPCE(IEL,IPCVSL) * CP0(IPHAS)
            ENDIF
          ELSE
            IF ( ICP(IPHAS).GT.0 ) THEN
              IPCCP   = IPPROC(ICP(IPHAS))
              W1(IEL) = VISLS0(IHM) * PROPCE(IEL,IPCCP)
            ELSE
              W1(IEL) = VISLS0(IHM) * CP0(IPHAS)
            ENDIF
          ENDIF
        ENDDO
C
C ------ Calcul du diametre des particules dans W2
C        On calcule le d20 = (A0.D0**2+(1-A0)*DCK**2)**0.5
C
        DO IEL = 1, NCEL
          W2(IEL) = ( XASHCH(NUMCHA)*DIAM20(NUMCLA)**2 +
     &                (1.D0-XASHCH(NUMCHA))*PROPCE(IEL,IPCDIA)**2
     &              )**0.5
        ENDDO
C
C ------ Contribution aux bilans explicite et implicite de
C        des echanges par diffusion moleculaire
C
        DO IEL = 1, NCEL
          IF ( W11(IEL).GT.EPSICP ) THEN
            AUX         = 6.D0 * W1(IEL) * XNUSS / W2(IEL)**2
     &                  / PROPCE(IEL,IPCRO2) * PROPCE(IEL,IPCROM)
     &                  * W11(IEL)
     &                  * VOLUME(IEL)
            RHOVST      = AUX / CP2CH(NUMCHA) /W11(IEL)
C
            SMBRS(IEL)  = SMBRS(IEL) - AUX *
     &                  ( PROPCE(IEL,IPCTE2) - PROPCE(IEL,IPCTE1) )
            ROVSDT(IEL) = ROVSDT(IEL) + MAX(ZERO,RHOVST)
          ENDIF
C
        ENDDO
C
C
C ---- Contribution aux bilans explicite et implicite
C        du terme echange d'energie entre les phases :
C        GAMA(dev1) H(mv1,T2)+GAMA(dev2) H(mv2,T2)
C
        DO IEL = 1, NCEL
C
C        Gama Dev1 et Gama Dev2
C
          GAMDV1 = PROPCE(IEL,IPCROM)*RTP(IEL,IXCHCL)
     &            *PROPCE(IEL,IPCGD1)
C
          GAMDV2 = PROPCE(IEL,IPCROM)*RTP(IEL,IXCHCL)
     &            *PROPCE(IEL,IPCGD2)
C
C        H(mv1,T2)
C
          DO IGE = 1, NGAZEM
            COEFE(IGE) = ZERO
          ENDDO
C
          COEFE(ICHX1) = A1(NUMCHA)*WMOLE(ICHX1C(NUMCHA))
     &  /( A1(NUMCHA)*WMOLE(ICHX1C(NUMCHA))
     &    +B1(NUMCHA)*WMOLE(ICO)
     &    +C1(NUMCHA)*WMOLE(IH2O) )
          COEFE(ICO  ) = B1(NUMCHA)*WMOLE(ICO)
     &   /( A1(NUMCHA)*WMOLE(ICHX1C(NUMCHA))
     &     +B1(NUMCHA)*WMOLE(ICO)
     &     +C1(NUMCHA)*WMOLE(IH2O) )
          COEFE(IH2O ) = C1(NUMCHA)*WMOLE(IH2O)
     &   /( A1(NUMCHA)*WMOLE(ICHX1C(NUMCHA))
     &     +B1(NUMCHA)*WMOLE(ICO)
     &     +C1(NUMCHA)*WMOLE(IH2O) )
C
          T2         = PROPCE(IEL,IPCTE2)
          DO ICHA = 1, NCHARM
            F1MC(ICHA) = ZERO
            F2MC(ICHA) = ZERO
          ENDDO
          F1MC(NUMCHA) = 1.D0
C
          MODE      = -1
          CALL CPTHP1
C         ===========
     &    ( MODE  , XHDEV1    , COEFE  , F1MC   , F2MC   ,
     &      T2    )
C
C        H(mv2,T2)
C
          DO IGE = 1, NGAZEM
            COEFE(IGE) = ZERO
          ENDDO
          COEFE(ICHX2) = A2(NUMCHA)*WMOLE(ICHX2C(NUMCHA))
     &   /( A2(NUMCHA)*WMOLE(ICHX2C(NUMCHA))
     &     +B2(NUMCHA)*WMOLE(ICO)
     &     +C2(NUMCHA)*WMOLE(IH2O) )
          COEFE(ICO  ) = B2(NUMCHA)*WMOLE(ICO)
     &   /( A2(NUMCHA)*WMOLE(ICHX2C(NUMCHA))
     &     +B2(NUMCHA)*WMOLE(ICO)
     &     +C2(NUMCHA)*WMOLE(IH2O) )
          COEFE(IH2O ) = C2(NUMCHA)*WMOLE(IH2O)
     &   /( A2(NUMCHA)*WMOLE(ICHX2C(NUMCHA))
     &     +B2(NUMCHA)*WMOLE(ICO)
     &     +C2(NUMCHA)*WMOLE(IH2O) )
C
          T2         = PROPCE(IEL,IPCTE2)
          DO ICHA = 1, NCHARM
            F1MC(ICHA) = ZERO
            F2MC(ICHA) = ZERO
          ENDDO
          F2MC(NUMCHA) = 1.D0
C
          MODE      = -1
          CALL CPTHP1
C         ===========
     &    ( MODE  , XHDEV2    , COEFE  , F1MC   , F2MC   ,
     &      T2    )
C
C         Contribution aux bilans explicite et implicite
C
          SMBRS(IEL) = SMBRS(IEL)
     &                 + (GAMDV1*XHDEV1+GAMDV2*XHDEV2)*VOLUME(IEL)
C
        ENDDO
C
C ------ combustion heterogene : C(s) + 02 ---> 0.5 C0
C        GamHET * (28/12 H(CO,T2)-16/12 H(O2,T1) )
C
        DO IEL = 1, NCEL
C
C        Calcul de HCO(T2)
C
          DO IGE = 1, NGAZEM
            COEFE(IGE) = ZERO
          ENDDO
          COEFE(ICO) = 1.D0
          DO ICHA = 1, NCHARM
            F1MC(ICHA) = ZERO
            F2MC(ICHA) = ZERO
          ENDDO
C
          T2        = PROPCE(IEL,IPCTE2)
          MODE      = -1
          CALL CPTHP1
C         ===========
     &    ( MODE  , XHCO    , COEFE  , F1MC   , F2MC   ,
     &      T2    )
C
C        Calcul de HO2(T1)
C
          DO IGE = 1, NGAZEM
            COEFE(IGE) = ZERO
          ENDDO
          COEFE(IO2) = 1.D0
          DO ICHA = 1, NCHARM
            F1MC(ICHA) = ZERO
            F2MC(ICHA) = ZERO
          ENDDO
C
          T1        = PROPCE(IEL,IPCTE1)
          MODE      = -1
          CALL CPTHP1
C         ===========
     &    ( MODE  , XHO2    , COEFE  , F1MC   , F2MC   ,
     &      T1    )
C
C         Contribution aux bilans explicite et implicite
C
          IF ( RTPA(IEL,IXCKCL) .GT. EPSICP ) THEN
C
            GAMHET = PROPCE(IEL,IPCROM)*PROPCE(IEL,IPCGHT)
     &               * ( (RTPA(IEL,IXCKCL))**(2.D0/3.D0) +
     &                2.D0/3.D0*(RTP(IEL,IXCKCL)-RTPA(IEL,IXCKCL))
     &                 /(RTPA(IEL,IXCKCL))**(1.D0/3.D0) )
C
          ELSE
            GAMHET = 0.D0
          ENDIF
C
         SMBRS(IEL) = SMBRS(IEL)
     &                 +GAMHET
     &                  *(28.D0/12.D0*XHCO-16.D0/12.D0*XHO2)
     &                  *VOLUME(IEL)
C
        ENDDO
C
C ------ combustion heterogene : C(s) + C02 ---> 2 C0
C        GamHET * (56/12 H(CO,T2)-44/12 H(O2,T1) )
C
        IF ( IHTCO2 .EQ. 1 ) THEN
          DO IEL = 1, NCEL
C
C        Calcul de HCO(T2)
C
            DO IGE = 1, NGAZEM
              COEFE(IGE) = ZERO
            ENDDO
            COEFE(ICO) = 1.D0
            DO ICHA = 1, NCHARM
              F1MC(ICHA) = ZERO
              F2MC(ICHA) = ZERO
            ENDDO
C
            T2        = PROPCE(IEL,IPCTE2)
            MODE      = -1
            CALL CPTHP1
C           ===========
     &      ( MODE  , XHCO    , COEFE  , F1MC   , F2MC   ,
     &        T2    )
C
C        Calcul de HCO2(T1)
C
            DO IGE = 1, NGAZEM
              COEFE(IGE) = ZERO
            ENDDO
            COEFE(ICO2) = 1.D0
            DO ICHA = 1, NCHARM
              F1MC(ICHA) = ZERO
              F2MC(ICHA) = ZERO
            ENDDO
C
            T1        = PROPCE(IEL,IPCTE1)
            MODE      = -1
            CALL CPTHP1
C           ===========
     &      ( MODE  , XHO2    , COEFE  , F1MC   , F2MC   ,
     &        T1    )
C
C         Contribution aux bilans explicite et implicite
C
            IF ( RTPA(IEL,IXCKCL) .GT. EPSICP ) THEN
C
              GAMHET = PROPCE(IEL,IPCROM)*PROPCE(IEL,IPGHC2)
     &                 * ( (RTPA(IEL,IXCKCL))**(2.D0/3.D0) +
     &                2.D0/3.D0*(RTP(IEL,IXCKCL)-RTPA(IEL,IXCKCL))
     &                 /(RTPA(IEL,IXCKCL))**(1.D0/3.D0) )
C
            ELSE
              GAMHET = 0.D0
            ENDIF
C
           SMBRS(IEL) = SMBRS(IEL)
     &                   +GAMHET
     &                    *(56.D0/12.D0*XHCO-44.D0/12.D0*XHO2)
     &                    *VOLUME(IEL)
C
          ENDDO
C
        ENDIF
C
C       --> Terme source sur H2 issu du sechage)
C
        IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
C
C ---- Contribution du TS interfacial aux bilans explicite et implicite
C
C
          IPCSEC = IPPROC(IGMSEC(NUMCLA))
          IPCTE2 = IPPROC(ITEMP2(NUMCLA))
C
          DO IEL = 1, NCEL
C
C          Calcul de H(H2O) a T2
C
            DO IGE = 1, NGAZEM
              COEFE(IGE) = ZERO
            ENDDO
            COEFE(IH2O) = 1.D0
            DO ICHA = 1, NCHARM
              F1MC(ICHA) = ZERO
              F2MC(ICHA) = ZERO
            ENDDO
C
            T2 = PROPCE(IEL,IPCTE2)
            IF ( T2 .GT. 100.D0 ) THEN
              T2 = 100.D0
            ENDIF
            MODE      = -1
            CALL CPTHP1
C           ===========
     &      ( MODE  , HH2OV    , COEFE  , F1MC   , F2MC   ,
     &        T2    )
C
            HLV = 2.263D+6
C
C         Contribution aux bilans explicite
C
            IF ( RTPA(IEL,ISCA(IXWT(NUMCLA))).GT. EPSICP .AND.
     &           XWATCH(NUMCHA) .GT. EPSICP       ) THEN
C
              AUX = -PROPCE(IEL,IPCROM)*PROPCE(IEL,IPCSEC)
     &       *(RTP(IEL,ISCA(IXWT(NUMCLA)))/PROPCE(IEL,IPCX2C))
     &       *(1.D0                    /XWATCH(NUMCHA))
     &             *HH2OV
C
            ELSE
              AUX = 0.D0
            ENDIF
C
            SMBRS(IEL) = SMBRS(IEL) + AUX*VOLUME(IEL)
C
          ENDDO
C
        ENDIF
C
      ENDIF
C
C=======================================================================
C 3. PRISE EN COMPTE DES TERMES SOURCES POUR LES VARIABLES RELATIVES
C    AU MELANGE GAZEUX
C=======================================================================
C
C --> Terme source pour les matieres volatiles legeres
C
      IF ( IVAR.GE.ISCA(IF1M(1)) .AND. IVAR.LE.ISCA(IF1M(NCHARB)) ) THEN
C
        IF (IWARNI(IVAR).GE.1) THEN
          WRITE(NFECRA,1000) CHAINE(1:8)
        ENDIF
C
C ---- Calcul de GMDEV1 = - SOMME (rho.XCH.GMDV1) > 0  --> W1
C
        NUMCHA = IVAR-ISCA(IF1M(1))+1
        DO IEL = 1, NCEL
          W1(IEL) = ZERO
        ENDDO
        DO ICLA = 1, NCLACP
          IPCGD1 = IPPROC(IGMDV1(ICLA))
          IXCHCL = ISCA(IXCH(ICLA))
          IF ( ICHCOR(ICLA).EQ.NUMCHA ) THEN
            DO IEL = 1, NCEL
              W1(IEL) = W1(IEL) - PROPCE(IEL,IPCROM)*RTP(IEL,IXCHCL)
     &                * PROPCE(IEL,IPCGD1)
            ENDDO
          ENDIF
        ENDDO
C
C ---- Contribution du TS interfacial aux bilans explicite et implicite
C
        DO IEL = 1, NCEL
          SMBRS(IEL)  = SMBRS(IEL)  + VOLUME(IEL) * W1(IEL)
        ENDDO
C
      ENDIF
C
C
C --> Terme source pour les matieres volatiles lourdes
C
      IF ( IVAR.GE.ISCA(IF2M(1)) .AND. IVAR.LE.ISCA(IF2M(NCHARB)) ) THEN
C
        IF (IWARNI(IVAR).GE.1) THEN
          WRITE(NFECRA,1000) CHAINE(1:8)
        ENDIF
C
C ---- Calcul de GMDEV2 = - SOMME (rho.XCH.GMDV2) >0 --> W1
C
        NUMCHA = IVAR-ISCA(IF2M(1))+1
        DO IEL = 1, NCEL
          W1(IEL) = ZERO
        ENDDO
        DO ICLA = 1, NCLACP
          IPCGD2 = IPPROC(IGMDV2(ICLA))
          IXCHCL = ISCA(IXCH(ICLA))
          IF ( ICHCOR(ICLA).EQ.NUMCHA ) THEN
            DO IEL = 1, NCEL
              W1(IEL) = W1(IEL) - PROPCE(IEL,IPCROM)*RTP(IEL,IXCHCL)
     &                * PROPCE(IEL,IPCGD2)
            ENDDO
          ENDIF
        ENDDO
C
C ---- Contribution du TS interfacial pour le bilan explicite
C
        DO IEL = 1, NCEL
          SMBRS(IEL)  = SMBRS(IEL)  + VOLUME(IEL) * W1(IEL)
        ENDDO
C
      ENDIF
C
C
C --> Terme source pour le traceur 3 (O2) (C de la comb. het.)
C
      IF ( IVAR.EQ.ISCA(IF3M) ) THEN
C
C RQ IMPORTANTE :  On prend les meme TS que pour Xck
C                  afin d'etre conservatif
C
        IF (IWARNI(IVAR).GE.1) THEN
          WRITE(NFECRA,1000) CHAINE(1:8)
        ENDIF
C
C ---- Calcul prelimimaire pour la partie explicite : W1
C                          pour la partie implicite : W2
C
        DO IEL = 1, NCEL
          W1(IEL) = ZERO
        ENDDO
C
        DO ICLA = 1, NCLACP
          IPCGHT = IPPROC(IGMHET(ICLA))
          IXCKCL = ISCA(IXCK(ICLA))
          DO IEL = 1, NCEL
            IF ( RTPA(IEL,IXCKCL) .GT. EPSICP ) THEN
              W1(IEL) = W1(IEL)
     &                 - PROPCE(IEL,IPCROM)*PROPCE(IEL,IPCGHT)
     &                 * ( (RTPA(IEL,IXCKCL))**(2.D0/3.D0) +
     &                    2.D0/3.D0*(RTP(IEL,IXCKCL)-RTPA(IEL,IXCKCL))
     &                    /(RTPA(IEL,IXCKCL))**(1.D0/3.D0) )
            ENDIF
          ENDDO
        ENDDO
C
C ---- Contribution du TS interfacial aux bilans explicite et implicite
C
        DO IEL = 1, NCEL
          SMBRS(IEL)  = SMBRS(IEL)  + VOLUME(IEL) * W1(IEL)
        ENDDO
C
      ENDIF
C
C
C
C --> Terme source pour le traceur 3 (CO2) (C de la comb. het.)
C
      IF ( IHTCO2 .EQ. 1 ) THEN
        IF ( IVAR.EQ.ISCA(IF3MC2) ) THEN
C
C RQ IMPORTANTE :  On prend les meme TS que pour Xck
C                  afin d'etre conservatif
C
          IF (IWARNI(IVAR).GE.1) THEN
            WRITE(NFECRA,1000) CHAINE(1:8)
          ENDIF
C
C ---- Calcul prelimimaire pour la partie explicite : W1
C                          pour la partie implicite : W2
C
          DO IEL = 1, NCEL
            W1(IEL) = ZERO
          ENDDO
C
          DO ICLA = 1, NCLACP
            IPCGHT = IPPROC(IGHCO2(ICLA))
            IXCKCL = ISCA(IXCK(ICLA))
            DO IEL = 1, NCEL
              IF ( RTPA(IEL,IXCKCL) .GT. EPSICP ) THEN
                W1(IEL) = W1(IEL)
     &                   - PROPCE(IEL,IPCROM)*PROPCE(IEL,IPCGHT)
     &                   * ( (RTPA(IEL,IXCKCL))**(2.D0/3.D0) +
     &                      2.D0/3.D0*(RTP(IEL,IXCKCL)-RTPA(IEL,IXCKCL))
     &                      /(RTPA(IEL,IXCKCL))**(1.D0/3.D0) )
              ENDIF
            ENDDO
          ENDDO
C
C ---- Contribution du TS interfacial aux bilans explicite et implicite
C
          DO IEL = 1, NCEL
            SMBRS(IEL)  = SMBRS(IEL)  + VOLUME(IEL) * W1(IEL)
          ENDDO
C
        ENDIF
C
      ENDIF
C
C --> Terme source pour la variance du traceur 4 (Air)
C
      IF ( IVAR.EQ.ISCA(IF4P2M) ) THEN
C
        IF (IWARNI(IVAR).GE.1) THEN
          WRITE(NFECRA,1000) CHAINE(1:8)
        ENDIF
C
C ---- Calcul des termes sources explicite et implicite
C      relatif aux echanges interfaciaux entre phases
C
       NUMTRA = 4
       CALL CPTSVI
C      ===========
     & ( NCELET , NCEL   , NUMTRA ,
     &   RTP    , PROPCE , VOLUME ,
     &   SMBRS  , ROVSDT ,
     &   W1     , W2     ,
     &   W3 )
C
C ---- Calcul des termes sources explicite et implicite
C      relatif aux termes de production et de dissipation
C
C      Pointeur relatif au scalaire associe
C      (0 si pas de scalaire associe)
       ISCALA = 0
C
       CALL CPTSVC
C      ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISCAL  , ISCALA ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   SMBRS  , ROVSDT ,
     &   VISCB  ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     ,
     &   RDEVEL , RTUSER , RA     )
C
      ENDIF
C
C --> Terme source pour le traceur 5 (Eau issue du séchage)
C
      IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
C
        IF ( IVAR.EQ.ISCA(IF5M) ) THEN
C
C
          IF (IWARNI(IVAR).GE.1) THEN
            WRITE(NFECRA,1000) CHAINE(1:8)
          ENDIF
C
C ---- Contribution du TS interfacial aux bilans explicite et implicite
C
          DO IEL = 1, NCEL
            W1(IEL) = ZERO
          ENDDO
C
          DO ICLA = 1, NCLACP
C
            IPCSEC = IPPROC(IGMSEC(ICLA))
            IPCX2C = IPPROC(IX2(ICLA))
            NUMCHA = ICHCOR(ICLA)
C
            DO IEL = 1, NCEL
C
              IF (  RTPA(IEL,ISCA(IXWT(ICLA))).GT. EPSICP
     &            .AND.
     &              XWATCH(NUMCHA) .GT. EPSICP       ) THEN
C
                W1(IEL) = W1(IEL)
     &      + PROPCE(IEL,IPCROM)*PROPCE(IEL,IPCSEC)
     &       *(RTP(IEL,ISCA(IXWT(ICLA)))/PROPCE(IEL,IPCX2C))
     &       *(1.D0                  /XWATCH(NUMCHA))
C
              ENDIF
C
            ENDDO
C
          ENDDO
C
          DO IEL = 1, NCEL
            SMBRS(IEL)  = SMBRS(IEL)  + VOLUME(IEL) * W1(IEL)
          ENDDO
C
        ENDIF
C
      ENDIF
C
C --> Terme source pour CO2
C
      IF ( IEQCO2 .GE. 1 ) THEN
C
        IF ( IVAR.EQ.ISCA(IYCO2) ) THEN
C
C
          IF (IWARNI(IVAR).GE.1) THEN
            WRITE(NFECRA,1000) CHAINE(1:8)
          ENDIF
C
C ---- Contribution du TS interfacial aux bilans explicite et implicite
C
C Oxydation du CO
C ===============
C
C  Dryer Glassman : XK0P en (moles/m3)**(-0.75) s-1
C
C          XK0P = 1.26D10
C
c           XK0P = 1.26D7 * (1.1)**(NTCABS)
c           IF ( XK0P .GT. 1.26D10 ) XK0P=1.26D10
C
C           T0P  = 4807.D0
C
C  Howard : XK0P en (moles/m3)**(-0.75) s-1
C
C             XK0P = 4.11D9
C             T0P  = 15090.D0
C
C  Westbrook & Dryer
C
          XK0P = 23.256D0
          T0P  = 20096.D0
C
C Dissociation du CO2 (Trinh Minh Chinh)
C ===================
C
C          XK0M = 5.D8
C          T0M  = 4807.D0
C
C          XK0M = 0.D0
C
C  Westbrook & Dryer
C
          XK0M = 20.03D0
          T0M  = 20096.D0
C
          ERR1MX = 0.D0
          ERR2MX = 0.D0
C
C Nombre d'iterations
C
          ITERMX = 500
C
C Nombre de points converges
C
         NBPAUV = 0
         NBEPAU = 0
         NBRICH = 0
         NBERIC = 0
C
C Precision pour la convergence
C
         ERRCH = 1.D-8
C
         DO IEL = 1, NCEL
C
           XXCO  = PROPCE(IEL,IPPROC(IYM1(ICO  )))/WMOLE(ICO)
     &            *PROPCE(IEL,IPPROC(IROM1))
           XXO2  = PROPCE(IEL,IPPROC(IYM1(IO2  )))/WMOLE(IO2)
     &             *PROPCE(IEL,IPPROC(IROM1))
           XXCO2 = PROPCE(IEL,IPPROC(IYM1(ICO2 )))/WMOLE(ICO2)
     &            *PROPCE(IEL,IPPROC(IROM1))
           XXH2O = PROPCE(IEL,IPPROC(IYM1(IH2O )))/WMOLE(IH2O)
     &            *PROPCE(IEL,IPPROC(IROM1))
C
           XXCO  = MAX(XXCO ,ZERO)
           XXO2  = MAX(XXO2 ,ZERO)
           XXCO2 = MAX(XXCO2,ZERO)
           XXH2O = MAX(XXH2O,ZERO)
C
           XCO2MX=XXCO2+MIN(XXCO,2.D0*XXO2)
C
C            XKP  = XK0P*EXP(-T0P/PROPCE(IEL,IPPROC(ITEMP1)))
C            XKM  = XK0M*EXP(-T0M/PROPCE(IEL,IPPROC(ITEMP1)))
C
           XKP = EXP(XK0P-T0P/PROPCE(IEL,IPPROC(ITEMP1)))
           XKM = EXP(XK0M-T0M/PROPCE(IEL,IPPROC(ITEMP1)))
C
C           Recherche de l'état d'équilibre
C
C           Recerche itérative sans controle de convergence
C            (pour conserver la parallelisation sur les mailles)
C           sur le nombre de moles de reaction séparant
C           l'etat avant réaction (tel que calculé par Cpcym)
C           de l'état d'équilibre
C
C           initialisation par l'état transporté
C
           ANMR  = XXCO2
           XCOT  = XXCO + XXCO2
           XO2T  = XXO2 + 0.5D0*XXCO2
           SQH2O = SQRT(XXH2O)
C
           ITERCH = 0
 10        CONTINUE
C
             XCO2E = ANMR
             XCOE  = XCOT-ANMR
             XO2E  = XO2T-0.5D0*ANMR
             ITERCH = ITERCH + 1
C
             TER1 = XKP*(XO2E**0.25D0)*SQH2O
             IF ( ABS(TER1) .GT. 1.D-15 ) THEN
               DDELTA = (TER1*XCOE-XKM*XCO2E)
     &                 /(TER1*(1.D0+XCOE/(8.D0*XO2E))+XKM)
             ELSE
               DDELTA = 0.D0
             ENDIF
C
             ANMR = ANMR + DDELTA
C
C  Clipping
C       Cote Pauvre
             IF ( XO2T .GT. 0.5D0*XCOT ) THEN
               ANMR = MAX(ZERO,MIN(XCOT,ANMR))
             ELSE
C       Cote Riche
               ANMR = MAX(ZERO,MIN(2.D0*XO2T,ANMR))
             ENDIF
C
C Test de convergence
C
             IF ( ABS(DDELTA) .GT. ERRCH
     &           .AND. ITERCH.LT.ITERMX       ) GOTO 10
C
C Calcul du nombre de points non convergé
C
           IF ( ITERCH.EQ.ITERMX
     &         .AND. ABS(DDELTA) .GT. ERRCH ) THEN
             NBERIC = NBERIC +1
           ENDIF
           ERR1MX = MAX(ERR1MX,ABS(DDELTA))
C
           XCO2E = ANMR
           XCOE  = XCOT-ANMR
           XO2E  = XO2T-0.5D0*ANMR
C
           XDEN = XKP*SQH2O*(0.5D0*(XO2T+XO2E))**0.25D0
           IF ( XDEN .NE. 0.D0 ) THEN
C
             TAUCHI = 1.D0/XDEN
             TAUTUR = RTPA(IEL,IK(IPHAS))/RTPA(IEL,IEP(IPHAS))
C
             X2 = 0.D0
             DO ICLA = 1, NCLACP
               X2 = X2 + PROPCE(IEL,IPPROC(IX2(ICLA)))
             ENDDO
C
             IF ( IEQCO2 .EQ. 1 ) THEN
C
C              On transporte CO2
C
               SMBRS(IEL)  = SMBRS(IEL)
     &                      +WMOLE(ICO2)/PROPCE(IEL,IPPROC(IROM1))
     &         * (XCO2E-XXCO2)/(TAUChI+TAUTUR)
     &         * (1.D0-X2)
     &         * VOLUME(IEL) * PROPCE(IEL,IPCROM)
C
             ELSE IF ( IEQCO2 .EQ. 2 ) THEN
C
C              On transporte CO
C
               SMBRS(IEL)  = SMBRS(IEL)
     &                      +WMOLE(ICO)/PROPCE(IEL,IPPROC(IROM1))
     &         * (XCOe-XXCO)/(TAUCHI+TAUTUR)
     &         * (1.D0-X2)
     &         * VOLUME(IEL) * PROPCE(IEL,IPCROM)
C
             ENDIF
C
             W1(IEL) = VOLUME(IEL)*PROPCE(IEL,IPCROM)/(TAUCHI+TAUTUR)
             ROVSDT(IEL) = ROVSDT(IEL) +   MAX(W1(IEL),ZERO)
C
           ELSE
             ROVSDT(IEL) = ROVSDT(IEL) + 0.D0
             SMBRS(IEL)  = SMBRS(IEL)  + 0.D0
           ENDIF
C
         ENDDO
C
         IF(IRANGP.GE.0) THEN
           CALL PARCPT(NBERIC)
           CALL PARMAX(ERR1MX)
         ENDIF
         write(NFECRA,*) ' Erreur max = ',ERR1MX
         write(NFECRA,*) ' Points non   ',NBERIC
C
C     Terme source : combustion heterogene par le CO2
C
         IF ( IHTCO2 .EQ. 1) THEN
C
           DO IEL = 1, NCEL
C
             AUX = 0.D0
             DO ICLA = 1,NCLACP
C
               IXCKCL = ISCA(IXCK(ICLA))
               IPGHC2 = IPPROC(IGHCO2(ICLA))
C
               AUX = AUX
     &              + PROPCE(IEL,IPCROM)*PROPCE(IEL,IPGHC2)
     &               *(RTPA(IEL,IXCKCL))**(2.D0/3.D0)*VOLUME(IEL)
C
             ENDDO
C
             ROVSDT(IEL) = ROVSDT(IEL) - AUX*(WMOLE(ICO2)/0.012)
C
           ENDDO
C
         ENDIF
C
        ENDIF
C
      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(' TERMES SOURCES PHYSIQUE PARTICULIERE POUR LA VARIABLE '
     &       ,A8,/)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
