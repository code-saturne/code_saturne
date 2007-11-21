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
                        SUBROUTINE CFPHYV
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
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC ROUTINE PHYSIQUE PARTICULIERE : COMPRESSIBLE SANS CHOC
CFONC
CFONC  Calcul des proprietes physiques variables
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
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! W1...3(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA, IFINIA
      INTEGER          IPHAS , IEL
      INTEGER          IFAC
      INTEGER          IIROM , IIROMB
      INTEGER          MAXELT, ILS
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C***********************************************************************
C=======================================================================
C 1. INITIALISATIONS A CONSERVER
C=======================================================================
C
C --- Initialisation memoire
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
C=======================================================================
C 2. ON DONNE LA MAIN A L'UTILISATEUR
C=======================================================================
C
      MAXELT = MAX(NCELET, NFAC, NFABOR)
      ILS    = IDEBIA
      IFINIA = ILS + MAXELT
      CALL IASIZE('CFPHYV',IFINIA)
C
      IUSCFP = 1
      CALL USCFPV
C     ===========
     & ( IFINIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , NPHMX  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER , RA     )
C
C     Si IUSCFP = 0, l'utilisateur n'a pas inclus le ss pgm uscfpv dans
C       ses sources. C'est une erreur si Cp, Cv ou Lambda est variable.
C     On se contente de faire le test au premier passage.
      IF(IPASS.EQ.0) THEN
        IPASS = IPASS + 1
        DO IPHAS = 1, NPHAS
          IF((IVISLS(ITEMPK(IPHAS)).GT.0.OR.
     &        ICP(IPHAS).GT.0.OR.ICV(IPHAS).GT.0).AND.IUSCFP.EQ.0) THEN
            WRITE(NFECRA,1000)
     &           IVISLS(ITEMPK(IPHAS)),ICP(IPHAS),ICV(IPHAS)
            CALL CSEXIT (1)
C           ===========
          ENDIF
        ENDDO
      ENDIF
C
C=======================================================================
C 3. MISE A JOUR DE LAMBDA/CV
C=======================================================================
C
C On a vérifié auparavant que CV0 était non nul.
C Si CV variable est nul, c'est une erreur utilisateur. On fait
C     un test à tous les passages (pas optimal), sachant que pour
C     le moment, on est en gaz parfait avec CV constant : si quelqu'un
C     essaye du CV variable, ce serait dommage que cela lui explose à la
C     figure pour de mauvaises raisons.
C Si IVISLS(IENERG(IPHAS)).EQ.0, on a forcement IVISLS(ITEMPK(IPHAS)).EQ.0
C     et ICV(IPHAS).EQ.0, par construction de IVISLS(IENERG(IPHAS)) dans
C     le sous-programme cfvarp
C
      DO IPHAS = 1, NPHAS
C
        IF(IVISLS(IENERG(IPHAS)).GT.0) THEN
C
          IF(IVISLS(ITEMPK(IPHAS)).GT.0) THEN
C
            DO IEL = 1, NCEL
              PROPCE(IEL,IPPROC(IVISLS(IENERG(IPHAS)))) =
     &             PROPCE(IEL,IPPROC(IVISLS(ITEMPK(IPHAS))))
            ENDDO
C
          ELSE
            DO IEL = 1, NCEL
              PROPCE(IEL,IPPROC(IVISLS(IENERG(IPHAS)))) =
     &             VISLS0(ITEMPK(IPHAS))
            ENDDO
C
          ENDIF
C
          IF(ICV(IPHAS).GT.0) THEN
C
            DO IEL = 1, NCEL
              IF(PROPCE(IEL,IPPROC(ICV(IPHAS))).LE.0.D0) THEN
                WRITE(NFECRA,2000)IEL,PROPCE(IEL,IPPROC(ICV(IPHAS)))
                CALL CSEXIT (1)
C               ===========
              ENDIF
            ENDDO
C
            DO IEL = 1, NCEL
              PROPCE(IEL,IPPROC(IVISLS(IENERG(IPHAS)))) =
     &             PROPCE(IEL,IPPROC(IVISLS(IENERG(IPHAS))))
     &             / PROPCE(IEL,IPPROC(ICV(IPHAS)))
            ENDDO
C
          ELSE
C
            DO IEL = 1, NCEL
              PROPCE(IEL,IPPROC(IVISLS(IENERG(IPHAS)))) =
     &             PROPCE(IEL,IPPROC(IVISLS(IENERG(IPHAS))))
     &             / CV0(IPHAS)
            ENDDO
C
          ENDIF
C
        ELSE
C
          VISLS0(IENERG(IPHAS)) = VISLS0(ITEMPK(IPHAS))/CV0(IPHAS)
C
        ENDIF
C
C
      ENDDO
C
C
C=======================================================================
C 3. MISE A JOUR DE ROM et ROMB :
C     On ne s'en sert a priori pas, mais la variable existe
C     On a ici des valeurs issues du pas de temps précédent (y compris
C       pour les conditions aux limites) ou issues de valeurs initiales
C     L'échange pério/parall sera fait dans phyvar.
C=======================================================================
C
      DO IPHAS = 1, NPHAS
C
        IIROM  = IPPROC(IROM  (IPHAS))
        IIROMB = IPPROB(IROM  (IPHAS))
C
        DO IEL = 1, NCEL
          PROPCE(IEL,IIROM)  = RTPA(IEL,ISCA(IRHO(IPHAS)))
        ENDDO
C
        DO IFAC = 1, NFABOR
          IEL = IFABOR(IFAC)
          PROPFB(IFAC,IIROMB) =
     &         COEFA(IFAC,ICLRTP(ISCA(IRHO(IPHAS)),ICOEF))
     &         + COEFB(IFAC,ICLRTP(ISCA(IRHO(IPHAS)),ICOEF))
     &         * RTPA(IEL,ISCA(IRHO(IPHAS)))
        ENDDO
C
      ENDDO
C
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION (MODULE COMPRESSIBLE)  ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@  Une ou plusieurs des propriétés suivantes a été déclarée  ',/,
     &'@    variable (repérée ci-dessous par un indicateur non nul) ',/,
     &'@    et une loi doit être fournie dans uscfpv.               ',/,
     &'@         propriété                               indicateur ',/,
     &'@     - conductivité thermique                    ',I10       ,/,
     &'@     - capacité calorifique à pression constante ',I10       ,/,
     &'@     - capacité calorifique à volume constant    ',I10       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Renseigner uscfpv ou déclarer les propriétés constantes et',/,
     &'@    uniformes (uscfx2 pour la conductivité thermique,       ',/,
     &'@    uscfth pour les capacités calorifiques).                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION (MODULE COMPRESSIBLE)  ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@  La capacité calorifique à volume constant présente (au    ',/,
     &'@    moins) une valeur négative ou nulle :                   ',/,
     &'@    cellule ',I10,   '  Cv = ',E18.9                         ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier uscfpv.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
