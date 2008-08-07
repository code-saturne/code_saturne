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
                        SUBROUTINE FUTSSC
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
CFONC ROUTINE PHYSIQUE PARTICULIERE : FLAMME FUEL
CFONC    ON PRECISE LES TERMES SOURCES POUR UN SCALAIRE PP
CFONC    SUR UN PAS DE TEMPS
CFONC
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
      INCLUDE "fuincl.h"
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
      INTEGER          NCEPDP , NCKPDP , NCESMP, NUMTRA
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
      INTEGER          IVAR , IPCROM, IEL, IPHAS , ICLA , NUMCLA
      INTEGER          IEXP1 , IEXP2 , IEXP3 , IFAC
      INTEGER          ISCALA
      INTEGER          IPCRO2 , IPCTE1 , IPCTE2
      INTEGER          IPCDIA , IPCVST
      INTEGER          IMODE  , IESP
      INTEGER          IPCGEV , IPCGHT , IPCHGL
      INTEGER          ITERMX,NBPAUV,NBRICH,NBEPAU,NBERIC
C
      DOUBLE PRECISION AUX, RHOVST
      DOUBLE PRECISION T1, H1, T2, H2
      DOUBLE PRECISION XNG,RHOFOL , ROM
      DOUBLE PRECISION GAMEVA,FDEV,FSD,FTRAC,HFOV
      DOUBLE PRECISION HO2,HCO,XESP(NGAZEM),XCOKE,T2MT1
      DOUBLE PRECISION GMECH,GMVAP,GMHET
C
      DOUBLE PRECISION XXCO,XXO2,XXCO2,XXH2O,XCO2MX
      DOUBLE PRECISION XKP,XK0P,XKM,XK0M,WPLUS,WMOINS,T0P,T0M
      DOUBLE PRECISION AUXP,AUXM, AUX1 , AUX2 , AUX3
C
      DOUBLE PRECISION XEQ,ANMR,XCOT,XO2T,XCO2E,XO2E,XCOE,TAUCHI,TAUTUR
      DOUBLE PRECISION SQH2O , X2 , WMHCN , WMNO
C
      INTEGER          ITERCH
      DOUBLE PRECISION ERR1MX,ERR2MX,ANMR0
C
      DOUBLE PRECISION ERRCH,TER1,DDELTA,FN,QPR
      DOUBLE PRECISION AUXMAX,AUXMIN
      DOUBLE PRECISION YMOY,VOLM,VOLMP,DMP
C
C***********************************************************************
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
C --- Temperature phase gaz
C
      IPCTE1 = IPPROC(ITEMP1)
C
C=======================================================================
C 2. PRISE EN COMPTE DES TERMES SOURCES POUR LES VARIABLES RELATIVES
C    AUX CLASSES DE PARTICULES
C=======================================================================
C
C --> Terme source pour l'enthalpie du liquide
C
      IF ( IVAR .GE. ISCA(IHLF(1))     .AND.
     &     IVAR .LE. ISCA(IHLF(NCLAFU))      ) THEN
C
        IF (IWARNI(IVAR).GE.1) THEN
          WRITE(NFECRA,1000) CHAINE(1:8)
        ENDIF
C
        NUMCLA = IVAR-ISCA(IHLF(1))+1
C
        IPCRO2 = IPPROC(IROM3 (NUMCLA))
        IPCDIA = IPPROC(IDIAM3(NUMCLA))
        IPCTE2 = IPPROC(ITEMP3(NUMCLA))
        IPCGEV = IPPROC(IGMEVA(NUMCLA))
        IPCGHT = IPPROC(IGMHTF(NUMCLA))
        IPCHGL = IPPROC(IH1HLF(NUMCLA))
C
C       La variable est l'enthalpie du liquide ramenée à
C       la masse de mélange
C       Les flux interfaciaux contribuent à faire varier l'enthalpie
C       du liquide
C       La vapeur emporte son enthalpie
C       flux = PROPCE(IEL,IPPROC(IGMEVA))
C       enthalpie massique reconstituée à partir de EHGAZE(IFOV )
C       à la température de la goutte
C       L'oxydation héterogène comporte un flux "entrant" d'O2
C       un flux sortant de CO
C       le flux net est celui du carbone
C       fluxIN  = 16/12 * PROPCE(IEL,IPPROC(IGMHTF))
C       fluxOUT = 28/12 * PROPCE(IEL,IPPROC(IGMHTF))
C       Enthalpie entrante reconstituée à partir de EHGAZE(IO2 )
C       à la température du gaz environnant
C       Enthalpie sortante reconstituée à partir de EHGAZE(ICO )
C       à la température du grain
C
        IMODE = -1
        DO IEL = 1, NCEL
c
          ROM = PROPCE(IEL,IPCROM)
c
          DO IESP = 1, NGAZEM
            XESP(IESP) = ZERO
          ENDDO
C
          XESP(IFOV) = 1.D0
          CALL FUTHP1(IMODE,HFOV,XESP,PROPCE(IEL,IPCTE2))
C
          XESP(IFOV) = ZERO
          XESP(IO2)  = 1.D0
          CALL FUTHP1(IMODE,HO2 ,XESP,PROPCE(IEL,IPCTE1))
C
          XESP(IO2)  = ZERO
          XESP(ICO)  = 1.D0
          CALL FUTHP1(IMODE,HCO ,XESP,PROPCE(IEL,IPCTE2))
C
          T2MT1 = PROPCE(IEL,IPCTE2)-PROPCE(IEL,IPCTE1)
C
          GMECH = -PROPCE(IEL,IPCHGL)*T2MT1
          GMVAP = PROPCE(IEL,IPCGEV)*HFOV*T2MT1
          GMHET = 16.D0/12.D0*PROPCE(IEL,IPCGHT)*HO2
     &           -28.D0/12.D0*PROPCE(IEL,IPCGHT)*HCO
C
          SMBRS(IEL) = SMBRS(IEL) +
     &         ( GMECH+GMVAP+GMHET )*ROM*VOLUME(IEL)
C
          RHOVST = ( PROPCE(IEL,IPCHGL)
     &              -PROPCE(IEL,IPCGEV)*HFOV )/CP2FOL
     &            *ROM*VOLUME(IEL)

          ROVSDT(IEL) = ROVSDT(IEL) +  MAX(ZERO,RHOVST)
C
        ENDDO
C
C --> T.S. pour la masse de liquide
C
      ELSEIF ( IVAR .GE. ISCA(IYFOL(1))     .AND.
     &         IVAR .LE. ISCA(IYFOL(NCLAFU))        ) THEN
C
        IF (IWARNI(IVAR).GE.1) THEN
          WRITE(NFECRA,1000) CHAINE(1:8)
        ENDIF
C
        NUMCLA = IVAR-ISCA(IYFOL(1))+1
C
        IPCRO2 = IPPROC(IROM3 (NUMCLA))
        IPCDIA = IPPROC(IDIAM3(NUMCLA))
        IPCTE2 = IPPROC(ITEMP3(NUMCLA))
        IPCGEV = IPPROC(IGMEVA(NUMCLA))
        IPCGHT = IPPROC(IGMHTF(NUMCLA))
        IPCHGL = IPPROC(IH1HLF(NUMCLA))
C
        DO IEL = 1, NCEL
C
          T2MT1 =  PROPCE(IEL,IPCTE2)-PROPCE(IEL,IPCTE1)
          GMVAP = -PROPCE(IEL,IPCGEV)*T2MT1
          GMHET = -PROPCE(IEL,IPCGHT)
C
          SMBRS(IEL) = SMBRS(IEL)
     &         - PROPCE(IEL,IPCROM)*VOLUME(IEL)*(GMVAP+GMHET)
C          ROVSDT(IEL) = ROVSDT(IEL) + ZERO
C
        ENDDO
C
C --> T.S. pour le traceur de la vapeur
C
      ELSEIF ( IVAR .EQ. ISCA(IFVAP) ) THEN
C
        IF (IWARNI(IVAR).GE.1) THEN
          WRITE(NFECRA,1000) CHAINE(1:8)
        ENDIF
C
        DO ICLA = 1, NCLAFU
C
          IPCTE2 = IPPROC(ITEMP3(ICLA))
          IPCGEV = IPPROC(IGMEVA(ICLA))
C
          DO IEL = 1, NCEL
C
            T2MT1 = PROPCE(IEL,IPCTE2)-PROPCE(IEL,IPCTE1)
            GMVAP = -PROPCE(IEL,IPCGEV)*T2MT1
C
            SMBRS(IEL) = SMBRS(IEL)
     &                 + GMVAP*PROPCE(IEL,IPCROM)*VOLUME(IEL)
C           ROVSDT(IEL) = ROVSDT(IEL) + ZERO
C
          ENDDO
C
        ENDDO
C
C --> T.S. pour le traceur du C ex réaction heterogene
C
      ELSEIF ( IVAR .EQ. ISCA(IFHTF) ) THEN
C
        IF (IWARNI(IVAR).GE.1) THEN
          WRITE(NFECRA,1000) CHAINE(1:8)
        ENDIF
C
        DO ICLA = 1, NCLAFU
C
          IPCGHT = IPPROC(IGMHTF(ICLA))
C
          DO IEL = 1, NCEL
C
            SMBRS(IEL) = SMBRS(IEL)
     &           -PROPCE(IEL,IPCROM)*PROPCE(IEL,IPCGHT)*VOLUME(IEL)
C           ROVSDT(IEL) = ROVSDT(IEL) + ZERO
C
          ENDDO
C
        ENDDO
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
       CALL FUTSVI
C      ===========
     & ( NCELET , NCEL   , NUMTRA ,
     &   RTP    , PROPCE , VOLUME ,
     &   SMBRS  , ROVSDT , W1      )
C
C ---- Calcul des termes sources explicite et implicite
C      relatif aux termes de production et de dissipation
C
C      Pointeur relatif au scalaire associe
C      (0 si pas de scalaire associe)
C
        ISCALA = 0
C
        CALL FUTSVC
C       ===========
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
          XK0P = 23.256
          T0P  = 20096.
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
          XK0M = 20.03
          T0M  = 20096.
C
          ERR1MX = 0.D0
          ERR2MX = 0.D0
C
C Nombre d'iterations
C
          ITERMX = 5000
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
     &            *PROPCE(IEL,IPCROM)
           XXO2  = PROPCE(IEL,IPPROC(IYM1(IO2  )))/WMOLE(IO2)
     &            *PROPCE(IEL,IPCROM)
           XXCO2 = PROPCE(IEL,IPPROC(IYM1(ICO2 )))/WMOLE(ICO2)
     &            *PROPCE(IEL,IPCROM)
           XXH2O = PROPCE(IEL,IPPROC(IYM1(IH2O )))/WMOLE(IH2O)
     &            *PROPCE(IEL,IPCROM)
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
C           if ( iel .eq. 36 ) then
           XKP = EXP(XK0P-T0P/PROPCE(IEL,IPPROC(ITEMP1)))
           XKM = EXP(XK0M-T0M/PROPCE(IEL,IPPROC(ITEMP1)))
           XEQ = EXP(XK0M-XK0P)
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
           XO2E = XO2T-0.5D0*ANMR
C
C    Si presence de O2 et H2O
C
           IF ( XO2E .GT. 0.D0 .AND. SQH2O .GT. 0.D0 ) THEN
C
               IF ( IEL.EQ.36 .OR. IEL.EQ.15151) THEN
                 write(*,*) ' TEST  = ',SQH2O,SQH2O .GT. 0.D0
               ENDIF
C
             ITERCH = 0
 10          CONTINUE
C
               XCO2E = ANMR
               XCOE = XCOT-ANMR
               XO2E = XO2T-0.5D0*ANMR
               ITERCH = ITERCH + 1
C
               IF ( IEL.EQ. 36 .OR. IEL.EQ.15151) THEN
                 write(*,*) ' ITERCH  = ',ITERCH
               ENDIF
C
               TER1 = XKP*SQH2O
C
               DDELTA = (TER1*XCOE*XO2E-XKM*XCO2E*(XO2E**0.75D0))
     &                 /(TER1*(XO2E+XCOE/8.D0)+XKM*(XO2E**0.75D0))
C
C
               ANMR = ANMR + DDELTA
C
C  Clipping
C       Cote Pauvre
              IF ( XO2T .GT. 0.5D0*XCOT ) THEN
                ANMR = MAX( ZERO, MIN(XCOT,ANMR))
              ELSE
C       Cote Riche
                ANMR = MAX( ZERO, MIN(2.D0*XO2T,ANMR))
              ENDIF
C
C Test de convergence
C
              IF ( ABS(DDELTA) .GT. ERRCH
     &             .AND. ITERCH.LT.ITERMX       ) GOTO 10
C
C Calcul du nombre de points non convergé
C
              IF ( ITERCH.EQ.ITERMX
     &            .AND. ABS(DDELTA) .GT. ERRCH ) THEN
C
                NBERIC = NBERIC +1
C
              ENDIF
 20          CONTINUE
             ERR1MX = MAX(ERR1MX,ABS(DDELTA))
C
             XCO2E = ANMR
             XCOE  = XCOT-ANMR
             XO2E  = XO2T-0.5D0*ANMR
C
             TAUChI = 1.D0/(XKP*SQH2O*(0.5D0*(XO2T+XO2E))**0.25D0)
             TAUTUR = RTPA(IEL,IK(IPHAS))/RTPA(IEL,IEP(IPHAS))
C
             X2 = 0.D0
             DO ICLA = 1,NCLAFU
               X2 = X2 + RTP(IEL,ISCA(IYFOL(ICLA)))
             ENDDO
C
C            On transporte CO2
C
             SMBRS(IEL)  = SMBRS(IEL)
     &                   +WMOLE(ICO2)/PROPCE(IEL,IPCROM)
     &        * (XCO2E-XXCO2)/(TAUChI+TAUTUR) * VOLUME(IEL)
     &        * (1.D0-X2)
C
             W1(IEL) = VOLUME(IEL)/(TAUCHI+TAUTUR)
             ROVSDT(IEL) = ROVSDT(IEL) +   MAX(W1(IEL),ZERO)
C
           ELSE
C
C     Pas de CO
C
             SMBRS(IEL)  = SMBRS(IEL) + 0.D0
             ROVSDT(IEL) = ROVSDT(IEL)+ 0.D0
C
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
        ENDIF
C
      ENDIF
C
C
C --> Terme source pour HCN et NO : uniquement a partir de la 2eme
C                                   iter
C
      IF ( IEQNOX .EQ. 1 .AND. NTCABS .GT. 1) THEN
C
        IF ( IVAR.EQ.ISCA(IYHCN) .OR. IVAR.EQ.ISCA(IYNO) ) THEN
C
          IEXP1  = IPPROC(IGHCN1)
          IEXP2  = IPPROC(IGHCN2)
          IEXP3  = IPPROC(IGNOTH)
C
C QPR= %N libéré pendant l'evaporation/taux de matieres volatiles
C          moyen
C
          QPR = 0.5D0
C
C YMOY = % vapeur en sorties
C
          YMOY  = 0.D0
          VOLM  = 0.D0
          DMP   = 0.D0
          VOLMP = 0.D0
          DO IFAC=1,NFABOR
C
            IF ( ITYPFB(IFAC,IPHAS).EQ. ISOLIB ) THEN
C
              IEL = IFABOR(IFAC)
C
              YMOY = YMOY + RTPA(IEL,ISCA(IFVAP))*VOLUME(IEL)
              VOLM = VOLM + VOLUME(IEL)
              DO ICLA=1,NCLAFU
                DMP  = DMP  + RTPA(IEL,ISCA(ING(ICLA)))
     &                       *VOLUME(IEL)*DINIFL(ICLA)
                VOLMP = VOLMP + VOLUME(IEL)
              ENDDO
C
            ENDIF
          ENDDO
C
          IF ( IRANGP .GE. 0 ) THEN
            CALL PARSOM(YMOY)
            CALL PARSOM(VOLM)
            CALL PARSOM(DMP)
            CALL PARSOM(VOLMP)
          ENDIF
C
          IF ( YMOY .GT. 0.D0 .AND. VOLM .GT. 0.D0 ) THEN
            YMOY = YMOY/VOLM
            IF ( DMP .GT. 0.D0 .AND. VOLMP .GT.0.D0 ) THEN
              YMOY = YMOY/(DMP/VOLMP)
              WRITE(NFECRA,*) ' Taux YMOY = ',YMOY
            ELSE
              WRITE(NFECRA,2000)
              YMOY = 0.5D0
            ENDIF
          ELSE
            WRITE(NFECRA,2000)
            YMOY = 0.5D0
          ENDIF
C
C Azote dans le fuel
C
          FN = 0.015
C
C Masse molaire
C
          WMHCN = 0.027D0
          WMNO  = 0.030D0
C
          IF ( IVAR.EQ.ISCA(IYHCN) ) THEN
C
C        Terme source HCN
C
            IF (IWARNI(IVAR).GE.1) THEN
              WRITE(NFECRA,1000) CHAINE(1:8)
            ENDIF
C
            AUXMIN = 1.D+20
            AUXMAX =-1.D+20
C
            DO IEL=1,NCEL
C
              AUX = VOLUME(IEL)*PROPCE(IEL,IPCROM)
     &             *(PROPCE(IEL,IEXP1)+PROPCE(IEL,IEXP2))
     &             *RTPA(IEL,ISCA(IYNO))
     &             *PROPCE(IEL,IPPROC(IMMEL))
     &             /WMOLE(IO2)
C
              SMBRS(IEL)  = SMBRS(IEL)
     &                     -AUX*RTPA(IEL,IVAR)
              ROVSDT(IEL) = ROVSDT(IEL) + AUX
C
              GMVAP = 0.D0
              GMHET = 0.D0
              DO ICLA=1,NCLAFU
C
                IPCGEV = IPPROC(IGMEVA(ICLA))
                IPCGHT = IPPROC(IGMHTF(ICLA))
                IPCTE2 = IPPROC(ITEMP3(ICLA))
                IPCTE1 = IPPROC(ITEMP1)
C
                GMVAP = GMVAP
     &                 + PROPCE(IEL,IPCROM)*VOLUME(IEL)
     &                  *PROPCE(IEL,IPCGEV)
     &                  *(PROPCE(IEL,IPCTE2)-PROPCE(IEL,IPCTE1))
C
                GMHET = GMHET
     &                 +PROPCE(IEL,IPCROM)*VOLUME(IEL)
     &                                    *PROPCE(IEL,IPCGHT)
C
              ENDDO
C
              AUX = -FN*WMHCN/(WMOLE(IN2)/2.D0)
     &                 *(QPR*GMVAP+(1.D0-QPR*YMOY)/YMOY*GMHET)
              SMBRS(IEL)  = SMBRS(IEL) + AUX
              AUXMIN = MIN(AUXMIN,AUX)
              AUXMAX = MAX(AUXMAX,AUX)
C
            ENDDO
            write(*,*) ' AUX MIN MAX = ',AUXMIN,AUXMAX
C
          ENDIF
C
          IF ( IVAR.EQ.ISCA(IYNO) ) THEN
C
C        Terme source NO
C
            IF (IWARNI(IVAR).GE.1) THEN
              WRITE(NFECRA,1000) CHAINE(1:8)
            ENDIF
C
            DO IEL=1,NCEL
C
              AUX1 = VOLUME(IEL)*PROPCE(IEL,IPCROM)
     &              *PROPCE(IEL,IEXP1)*RTPA(IEL,ISCA(IYHCN))
     &              *PROPCE(IEL,IPPROC(IMMEL))/WMHCN
              AUX2 = VOLUME(IEL)*PROPCE(IEL,IPCROM)
     &              *PROPCE(IEL,IEXP2)*RTPA(IEL,ISCA(IYHCN))
     &              *WMNO/WMHCN
              AUX3 = VOLUME(IEL)*PROPCE(IEL,IPCROM)**1.5D0
     &              *PROPCE(IEL,IEXP3)
     &              *(PROPCE(IEL,IPPROC(IYM1(IO2))))*0.5D0
     &              * PROPCE(IEL,IPPROC(IYM1(IN2)))
C
              SMBRS(IEL)  = SMBRS(IEL) - AUX1*RTPA(IEL,IVAR)
     &                               + AUX2 + AUX3
              ROVSDT(IEL) = ROVSDT(IEL) + AUX1
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
 2000 FORMAT(' ATTENTION : LE TAUX DE VAPEUR EN SORTIE',
     &       ' EST NULLE',/,
     &       ' ON PREND QPR = 0.5')
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
