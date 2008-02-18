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
      INTEGER          IVAR , IPCROM, IEL, IPHAS
      INTEGER          IPCGHT
      INTEGER          ISCALA
      INTEGER          IPCRO2 , IPCTE1 , IPCTE2
      INTEGER          IPCDIA , IPCVST
      INTEGER          IESP
      INTEGER          IMODE
C
      DOUBLE PRECISION RHOVST
      DOUBLE PRECISION ROM
      DOUBLE PRECISION HFOV
      DOUBLE PRECISION HO2,HCO,XESP(NGAZEM),T2MT1
      DOUBLE PRECISION GMECH,GMVAP,GMHET
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
C=======================================================================
C 2. PRISE EN COMPTE DES TERMES SOURCES POUR LES VARIABLES RELATIVES
C    AUX CLASSES DE PARTICULES
C=======================================================================
C
      IPCRO2 = IPPROC(IROM3)
      IPCDIA = IPPROC(IDIAM3)
      IPCTE2 = IPPROC(ITEMP3)
      IPCTE1 = IPPROC(ITEMP1)
      IPCGHT = IPPROC(IGMHTF)
C
C
C --> Terme source pour l'enthalpie du liquide
C
      IF ( IVAR .EQ. ISCA(IHLF) ) THEN
C
        IF (IWARNI(IVAR).GE.1) THEN
          WRITE(NFECRA,1000) CHAINE(1:8)
        ENDIF
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
          GMECH = -PROPCE(IEL,IPPROC(IH1HLF))*T2MT1
          GMVAP = PROPCE(IEL,IPPROC(IGMEVA))*HFOV*T2MT1
          GMHET = 16.D0/12.D0*PROPCE(IEL,IPPROC(IGMHTF))*HO2
     &           -28.D0/12.D0*PROPCE(IEL,IPPROC(IGMHTF))*HCO
C
          SMBRS(IEL) = SMBRS(IEL) +
     &         ( GMECH+GMVAP+GMHET )*ROM*VOLUME(IEL)
C
          RHOVST = ( PROPCE(IEL,IPPROC(IH1HLF))
     &              -PROPCE(IEL,IPPROC(IGMEVA))*HFOV )/CP2FOL
     &            *ROM*VOLUME(IEL)

          ROVSDT(IEL) = ROVSDT(IEL) +  MAX(ZERO,RHOVST)
C
        ENDDO
C
C --> T.S. pour la masse de liquide
C
      ELSEIF ( IVAR .EQ. ISCA(IYFOL) ) THEN
C
        IF (IWARNI(IVAR).GE.1) THEN
          WRITE(NFECRA,1000) CHAINE(1:8)
        ENDIF
C
        DO IEL = 1, NCEL
C
          T2MT1 =  PROPCE(IEL,IPCTE2)-PROPCE(IEL,IPCTE1)
          GMVAP = -PROPCE(IEL,IPPROC(IGMEVA))*T2MT1
          GMHET = -PROPCE(IEL,IPPROC(IGMHTF))
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
        DO IEL = 1, NCEL
C
          T2MT1 = PROPCE(IEL,IPCTE2)-PROPCE(IEL,IPCTE1)
          GMVAP = -PROPCE(IEL,IPPROC(IGMEVA))*T2MT1
C
          SMBRS(IEL) = SMBRS(IEL)
     &               + GMVAP*PROPCE(IEL,IPCROM)*VOLUME(IEL)
C          ROVSDT(IEL) = ROVSDT(IEL) + ZERO
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
        DO IEL = 1, NCEL
C
          SMBRS(IEL) = SMBRS(IEL)
     &         -PROPCE(IEL,IPCROM)*PROPCE(IEL,IPCGHT)*VOLUME(IEL)
C          ROVSDT(IEL) = ROVSDT(IEL) + ZERO
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
      IF ( IEQCO2 .EQ. 1 ) THEN
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
          DO IEL = 1, NCEL
C
C            SMBRS (IEL) = SMBRS (IEL) + ZERO
C            ROVSDT(IEL) = ROVSDT(IEL) + ZERO
C
          ENDDO
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
