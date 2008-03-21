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
                        SUBROUTINE EBUINI
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB , COEFA  , COEFB  ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  INITIALISATION DES VARIABLES DE CALCUL
CFONC     POUR LA PHYSIQUE PARTICULIERE : COMBUSTION GAZ MODELE EBU
CFONC     PENDANT DE USINIV.F
CFONC
CFONC  Cette routine est appelee en debut de calcul (suite ou non)
CFONC      avant le debut de la boucle en temps
CFONC
CFONC  Elle permet d'INITIALISER ou de MODIFIER (pour les calculs suite)
CFONC      les variables de calcul,
CFONC      les valeurs du pas de temps
CFONC
CFONC
CFONC  On dispose ici de ROM et VISCL initialises par RO0 et VISCL0
CFONC      ou relues d'un fichier suite
CFONC  On ne dispose des variables VISCLS, CP (quand elles sont
CFONC      definies) que si elles ont pu etre relues dans un fichier
CFONC      suite de calcul
CFONC
CFONC  Les proprietes physiques sont accessibles dans le tableau
CFONC      PROPCE (prop au centre), PROPFA (aux faces internes),
CFONC      PROPFB (prop aux faces de bord)
CFONC      Ainsi,
CFONC       PROPCE(IEL,IPPROC(IROM  (IPHAS))) designe ROM   (IEL ,IPHAS)
CFONC       PROPCE(IEL,IPPROC(IVISCL(IPHAS))) designe VISCL (IEL ,IPHAS)
CFONC       PROPCE(IEL,IPPROC(ICP   (IPHAS))) designe CP    (IEL ,IPHAS)
CFONC       PROPCE(IEL,IPPROC(IVISLS(ISCAL))) designe VISLS (IEL ,ISCAL)
CFONC
CFONC       PROPFA(IFAC,IPPROF(IFLUMA(IVAR ))) designe FLUMAS(IFAC,IVAR)
CFONC
CFONC       PROPFB(IFAC,IPPROB(IROM  (IPHAS))) designe ROMB  (IFAC,IPHAS)
CFONC       PROPFB(IFAC,IPPROB(IFLUMA(IVAR ))) designe FLUMAB(IFAC,IVAR)
CFONC
CFONC  LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
CFONC      SE FERA EN STANDARD DANS LE SOUS PROGRAMME PPPHYV
CFONC      ET PAS ICI
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
CARGU ! DT(NCELET)   ! TR ! <-> ! VALEUR DU PAS DE TEMPS               !
CARGU ! RTP          ! TR ! <-> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA COEFB  ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
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
      INCLUDE "period.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
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
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IEL, MODE, IGG, IPHAS, IZONE
      INTEGER          ISCAL, IVAR, II, IDIMTE, ITENSO
      DOUBLE PRECISION HINIT, COEFG(NGAZGM), HAIR, TINITK
      DOUBLE PRECISION SOMMQF, SOMMQT, SOMMQ, TENTM, FMELM
      DOUBLE PRECISION VALMAX, VALMIN, XKENT, XEENT, D2S3
C
C NOMBRE DE PASSAGES DANS LA ROUTINE
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATION VARIABLES LOCALES
C=======================================================================
C
      IPASS = IPASS + 1
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      DO IGG = 1, NGAZGM
        COEFG(IGG) = ZERO
      ENDDO
C
      IPHAS  = 1
C
      D2S3 = 2.D0/3.D0
C
C=======================================================================
C 2. INITIALISATION DES INCONNUES :
C      UNIQUEMENT SI ON NE FAIT PAS UNE SUITE
C=======================================================================
C
      IF ( ISUITE.EQ.0 ) THEN
C
C ---> Initialisation au 1er passage avec de l'air a TINITK
C                                    ======================
C
        IF ( IPASS.EQ.1 ) THEN
C
C ----- Temperature du melange : air a TINITK
          TINITK = T0(IPHAS)
C
C ----- Enthalpie de l'air a TINITK
          IF ( IPPMOD(ICOEBU).EQ.1 .OR. IPPMOD(ICOEBU).EQ.3 ) THEN
            COEFG(1) = ZERO
            COEFG(2) = 1.D0
            COEFG(3) = ZERO
            MODE     = -1
            CALL COTHHT
C           ===========
     &        ( MODE   , NGAZG , NGAZGM  , COEFG  ,
     &          NPO    , NPOT   , TH     , EHGAZG ,
     &          HAIR   , TINITK )
          ENDIF
C
C ----- On en profite pour initialiser FRMEL et TGF
C       CAR on n'a pas encore vu usebuc.F
C
          FRMEL = ZERO
          TGF   = 300.D0
C
C ---- Initialisation de k et epsilon
C
          XKENT = 1.D-10
          XEENT = 1.D-10
C
          DO IEL = 1, NCEL
C
C ---- TURBULENCE
C
            IF (ITYTUR(IPHAS).EQ.2) THEN
C
              RTP(IEL,IK(IPHAS))  = XKENT
              RTP(IEL,IEP(IPHAS)) = XEENT
C
            ELSEIF (ITYTUR(IPHAS).EQ.3) THEN
C
              RTP(IEL,IR11(IPHAS)) = D2S3*XKENT
              RTP(IEL,IR22(IPHAS)) = D2S3*XKENT
              RTP(IEL,IR33(IPHAS)) = D2S3*XKENT
              RTP(IEL,IR12(IPHAS)) = 0.D0
              RTP(IEL,IR13(IPHAS)) = 0.D0
              RTP(IEL,IR23(IPHAS)) = 0.D0
              RTP(IEL,IEP(IPHAS))  = XEENT
C
            ELSEIF (ITURB(IPHAS).EQ.50) THEN
C
              RTP(IEL,IK(IPHAS))   = XKENT
              RTP(IEL,IEP(IPHAS))  = XEENT
              RTP(IEL,IPHI(IPHAS)) = D2S3
              RTP(IEL,IFB(IPHAS))  = 0.D0
C
            ELSEIF (ITURB(IPHAS).EQ.60) THEN
C
              RTP(IEL,IK(IPHAS))   = XKENT
              RTP(IEL,IOMG(IPHAS)) = XEENT/CMU/XKENT
C
            ENDIF
C
C ----- Fraction massique de gaz frais
C
            RTP(IEL,ISCA(IYGFM)) = 1.D0
C
C ----- Fraction de melange
C
            IF ( IPPMOD(ICOEBU).EQ.2 .OR. IPPMOD(ICOEBU).EQ.3 ) THEN
              RTP(IEL,ISCA(IFM)) = ZERO
            ENDIF
C
C ----- Enthalpie du melange
C
            IF ( IPPMOD(ICOEBU).EQ.1 .OR. IPPMOD(ICOEBU).EQ.3 ) THEN
              RTP(IEL,ISCA(IHM)) = HAIR
            ENDIF
C
          ENDDO
C
C ---> Initialisation au 2eme passage
C
        ELSE IF ( IPASS.EQ.2 ) THEN
C
C ----- Calculs preliminaires : Fraction de melange, T, H
C     (la valeur NOZAPM est utilisee pour inclure les aspects parall)
          SOMMQF = ZERO
          SOMMQ  = ZERO
          SOMMQT = ZERO
          DO IZONE = 1, NOZAPM
            SOMMQF = SOMMQF + QIMP(IZONE)*FMENT(IZONE)
            SOMMQT = SOMMQT + QIMP(IZONE)*TKENT(IZONE)
            SOMMQ  = SOMMQ  + QIMP(IZONE)
          ENDDO
C
          IF(ABS(SOMMQ).GT.EPZERO) THEN
            FMELM = SOMMQF / SOMMQ
            TENTM = SOMMQT / SOMMQ
          ELSE
            FMELM = ZERO
            TENTM = T0(IPHAS)
          ENDIF
C
C ----- Enthalpie du melange HINIT
          IF ( IPPMOD(ICOEBU).EQ.1 .OR. IPPMOD(ICOEBU).EQ.3 ) THEN
            COEFG(1) = FMELM
            COEFG(2) = (1.D0-FMELM)
            COEFG(3) = ZERO
            MODE     = -1
            CALL COTHHT
C           ===========
     &        ( MODE   , NGAZG , NGAZGM  , COEFG  ,
     &          NPO    , NPOT   , TH     , EHGAZG ,
     &          HINIT  , TENTM )
          ENDIF
C
C
          DO IEL = 1, NCEL
C
C ----- Fraction massique de gaz frais
C
            RTP(IEL,ISCA(IYGFM)) = 5.D-1
C
C ----- Fraction de melange
C
            IF ( IPPMOD(ICOEBU).EQ.2 .OR. IPPMOD(ICOEBU).EQ.3 ) THEN
              RTP(IEL,ISCA(IFM)) = FMELM
            ENDIF
C
C ----- Enthalpie du melange
C
            IF ( IPPMOD(ICOEBU).EQ.1 .OR. IPPMOD(ICOEBU).EQ.3 ) THEN
              RTP(IEL,ISCA(IHM)) = HINIT
            ENDIF
C
          ENDDO
C
C ----- On donne la main a l'utilisateur
C
          CALL USEBUI
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB , COEFA  , COEFB  ,
     &   RDEVEL , RTUSER , RA     )
C
C ----- En periodique et en parallele,
C       il faut echanger ces initialisations (qui sont en fait dans RTPA)
C
C     Parallele
          IF(IRANGP.GE.0) THEN
            CALL PARCOM(RTP(1,ISCA(IYGFM)))
C           ===========
            IF ( IPPMOD(ICOEBU).EQ.2 .OR. IPPMOD(ICOEBU).EQ.3 ) THEN
              CALL PARCOM(RTP(1,ISCA(IFM  )))
C             ===========
            ENDIF
            IF ( IPPMOD(ICOEBU).EQ.1 .OR. IPPMOD(ICOEBU).EQ.3 ) THEN
              CALL PARCOM(RTP(1,ISCA(IHM  )))
C             ===========
            ENDIF
          ENDIF
C
C     Periodique
          IF(IPERIO.EQ.1) THEN
            IDIMTE = 0
            ITENSO = 0
            IVAR   = ISCA(IYGFM)
            CALL PERCOM
C           ===========
     &      ( IDIMTE , ITENSO ,
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR))
            IF ( IPPMOD(ICOEBU).EQ.2 .OR. IPPMOD(ICOEBU).EQ.3 ) THEN
              IDIMTE = 0
              ITENSO = 0
              IVAR   = ISCA(IFM  )
              CALL PERCOM
C             ===========
     &      ( IDIMTE , ITENSO ,
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR))
            ENDIF
            IF ( IPPMOD(ICOEBU).EQ.1 .OR. IPPMOD(ICOEBU).EQ.3 ) THEN
              IDIMTE = 0
              ITENSO = 0
              IVAR   = ISCA(IHM  )
              CALL PERCOM
C             ===========
     &      ( IDIMTE , ITENSO ,
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR),
     &        RTP(1,IVAR), RTP(1,IVAR), RTP(1,IVAR))
            ENDIF
          ENDIF
C
C
C      Impressions de controle
C
          WRITE(NFECRA,2000)
C
          DO II  = 1, NSCAPP
            ISCAL = ISCAPP(II)
            IVAR  = ISCA(ISCAL)
            VALMAX = -GRAND
            VALMIN =  GRAND
            DO IEL = 1, NCEL
              VALMAX = MAX(VALMAX,RTP(IEL,IVAR))
              VALMIN = MIN(VALMIN,RTP(IEL,IVAR))
            ENDDO
            CHAINE = NOMVAR(IPPRTP(IVAR))
            IF (IRANGP.GE.0) THEN
              CALL PARMIN(VALMIN)
C             ===========
              CALL PARMAX(VALMAX)
C             ===========
            ENDIF
            WRITE(NFECRA,2010)CHAINE(1:8),VALMIN,VALMAX
          ENDDO
C
          WRITE(NFECRA,2020)
C
        ENDIF
C
      ENDIF
C
C----
C FORMATS
C----
C
C
 2000 FORMAT(
     &'                                                             ',/,
     &' ----------------------------------------------------------- ',/,
     &'                                                             ',/,
     &'                                                             ',/,
     &' ** INITIALISATION DES VARIABLES PROPRES AU GAZ (FL PRE EBU) ',/,
     &'    -------------------------------------------------------- ',/,
     &'           2eme PASSAGE                                      ',/,
     &' ---------------------------------                           ',/,
     &'  Variable  Valeur min  Valeur max                           ',/,
     &' ---------------------------------                           '  )
C
 2010 FORMAT(
     & 2X,     A8,      E12.4,      E12.4                              )
C
 2020 FORMAT(
     &' ---------------------------------                           ',/)
C
C----
C FIN
C----
C
      RETURN
      END
c@z
