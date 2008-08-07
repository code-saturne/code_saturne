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
                        SUBROUTINE PPCABS
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER ,
     &   CK   ,   RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC
CFONC    SOUS-PROGRAMME PHYSIQUES PARTICULIERES
CFONC
CFONC   DONNE LA VALEUR DU COEFFICIENT D'ABSORPTION POUR
CFONC     LE MELANGE GAZEUX ET LES PARTICULES POUR LE CP.
CFONC
CFONC   NOTER QUE IPHAS VAUT 1.
CFONC
CFONC   DANS LE CAS DU MODELE P-1 ON VERIFIE QUE LA LONGUEUR OPTIQUE
CFONC     DU MILIEU EST AU MINIMUM DE L'ORDRE DE L'UNITE
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
CARGU !  NPHAS      )!    !     !                                      !
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
CARGU ! W1...3(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! CK (NCELET,*)! TR ! <-  ! COEFFICIENT D'ABSORPTION DU MILIEU   !
CARGU !              !    !     ! (NUL SI TRANSPARENT) ET DES CHARBONS !
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
      INCLUDE "parall.h"
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
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML) , ITYPFB(NFABOR)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
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
      DOUBLE PRECISION CK(NCELET,*)
C
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA, IEL, IFAC, ICLA, IPCK, ICHA, IOK
      DOUBLE PRECISION XM, D2, VV, SF, XLC, XKMIN, PP
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
C  1 - COEFFICIENT D'ABSORPTION DU MELANGE GAZEUX (m-1)
C=======================================================================
C
      IF ( IPPMOD(ICOD3P).GE.0 .OR. IPPMOD(ICOEBU).GE.0 ) THEN
C
C ----> Combustion gaz : Flamme de diffusion
C                        Flamme de premelange (Modele EBU)
C
        IF (IMODAK.EQ.1) THEN
C
          DO IEL = 1, NCEL
            XM = 1.D0/ (  PROPCE(IEL,(IYM(1)))/WMOLG(1)
     &                  + PROPCE(IEL,(IYM(2)))/WMOLG(2)
     &                  + PROPCE(IEL,(IYM(3)))/WMOLG(3) )
            W1(IEL) = PROPCE(IEL,(IYM(3)))*XM/WMOLG(3)*XCO2
            W2(IEL) = PROPCE(IEL,(IYM(3)))*XM/WMOLG(3)*XH2O
            W3(IEL) = 0.D0
          ENDDO
          CALL RAYDAK(NCEL,NCELET,
C         ===========
     &                CK,W1,W2,W3,PROPCE(1,IPPROC(ITEMP)))

          write(NFECRA,*) ' a verifier '
          write(NFECRA,*) ' a finir   : raydak '
          write(NFECRA,*) ' Le codage est a terminer par le groupe I81'
          write(NFECRA,*) '                     13-10-03 22:38:03      '
          CALL CSEXIT(1)
C
        ELSE
          DO IEL = 1, NCEL
            CK(IEL,1) = PROPCE(IEL,IPPROC(ICKABS))
          ENDDO
        ENDIF
C
      ELSE IF ( IPPMOD(ICP3PL).GE.0 ) THEN
C
C ---->  Charbon
C
        IF (IMODAK.EQ.1) THEN
C
          DO IEL = 1,NCEL
C concentration volumique en CO2
            W1(IEL) = PROPCE(IEL,IPPROC(IMMEL))/WMOLE(ICO2)
     &               *PROPCE(IEL,IPPROC(IYM1(ICO2)))
C concentration volumique en H20
            W2(IEL) = PROPCE(IEL,IPPROC(IMMEL))/WMOLE(IH2O)
     &               *PROPCE(IEL,IPPROC(IYM1(IH2O)))
C fraction volumique de suies
            W3(IEL) = 0.D0
C
          ENDDO
C
          CALL RAYDAK(NCEL,NCELET,
C         ===========
     &                CK,W1,W2,W3,PROPCE(1,IPPROC(ITEMP1)))
C
        ELSE
          DO IEL = 1, NCEL
            CK(IEL,1) = CKABS1
          ENDDO
        ENDIF
C
      ELSE IF ( IPPMOD(ICFUEL).GE.0 ) THEN
C
C ---->  Fuel
C
        IF (IMODAK.EQ.1) THEN
C
          DO IEL = 1,NCEL
C concentration volumique en CO2
            W1(IEL) = PROPCE(IEL,IPPROC(IMMEL))/WMOLE(ICO2)
     &               *PROPCE(IEL,IPPROC(IYM1(ICO2)))
C concentration volumique en H20
            W2(IEL) = PROPCE(IEL,IPPROC(IMMEL))/WMOLE(IH2O)
     &               *PROPCE(IEL,IPPROC(IYM1(IH2O)))
C fraction volumique de suies
            W3(IEL) = 0.D0
C
          ENDDO
C
          CALL RAYDAK(NCEL,NCELET,
C         ===========
     &                CK,W1,W2,W3,PROPCE(1,IPPROC(ITEMP1)))
C
        ELSE
          DO IEL = 1, NCEL
            CK(IEL,1) = CKABS1
          ENDDO
        ENDIF
C
      ENDIF
C
C
C=======================================================================
C  2 - COEFFICIENT D'ABSORPTION DES PARTICULES PAR CLASSE K2/X2 (m-1)
C=======================================================================
C
C ---->  Charbon
C
      IF ( IPPMOD(ICP3PL).GE.0 ) THEN
C
        DO ICLA = 1, NCLACP
C
          IPCK = 1 + ICLA
          ICHA = ICHCOR(ICLA)
C
          DO IEL = 1, NCEL
C
C ---> Calcul du diametre des particules
C
            D2 = ( XASHCH(ICHA)*DIAM20(ICLA)**2 +
     &           ( 1.D0-XASHCH(ICHA))
     &             *PROPCE(IEL,IPPROC(IDIAM2(ICLA)))**2 )**0.5D0
C
C ---> Calcul du coeficient d'absorption des particules K2/X2
C         3./2. ROM/(ROM2*D2)
C
            CK(IEL,IPCK) = 1.5D0*PROPCE(IEL,IPPROC(IROM(IPHAS)))
     &                       / ( PROPCE(IEL,IPPROC(IROM2(ICLA)))*D2)
C
          ENDDO
C
        ENDDO
C
      ENDIF
C
C ---->  Fuel
C
      IF ( IPPMOD(ICFUEL).GE.0 ) THEN
C
        DO ICLA = 1, NCLAFU
C
          IPCK = 1 + ICLA
C
          DO IEL = 1, NCEL
C
C ---> Calcul du coeficient d'absorption des particules K2/X2
C         3./2. ROM/(ROM2*D2)
C
            CK(IEL,IPCK) = 1.5D0*PROPCE(IEL,IPPROC(IROM(IPHAS)))
     &                 / ( PROPCE(IEL,IPPROC(IROM3(ICLA)))
     &                    *PROPCE(IEL,IPPROC(IDIAM3(ICLA))) )
C
          ENDDO
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C  3 - COEFFICIENT D'ABSORPTION GAZ (Arc Electrique)
C=======================================================================
C
C
      IF ( IPPMOD(IELARC).GE.1 ) THEN
C
        DO IEL = 1, NCEL
C
C ---> Directement donne par le fichier dp_elec
C
            CK(IEL,1) = PROPCE(IEL,IPPROC(IDRAD))
C
        ENDDO
C
      ENDIF
C
C=======================================================================
C  4 - CLIPPING DU COEFFICIENT D'ABSORPTION DANS LA CAS DE L'APPROX P-1
C=======================================================================
C
C
C--> MODELE P-1 : Controle standard des valeurs du coefficient
C                 d'absorption. Ce coefficient doit assurer une
C                 longueur optique au minimum de l'ordre de l'unite.
C
C
        IF (IRAYON(1).EQ.2) THEN
C
C         Coefficient d'absorption du melange gaz-particules de charbon
C
           DO IEL = 1, NCEL
             W3(IEL) =  CK(IEL,1)
           ENDDO
C
           IF ( IPPMOD(ICP3PL).GE.0 ) THEN
             DO ICLA = 1,NCLACP
               IPCK = 1+ICLA
               DO IEL = 1,NCEL
                 W3(IEL) = W3(IEL)
     &                    + PROPCE(IEL,IPPROC(IX2(ICLA)))
     &                     *CK(IEL,IPCK)
               ENDDO
             ENDDO
           ELSEIF ( IPPMOD(ICFUEL).GE.0 ) THEN
             DO ICLA = 1,NCLAFU
               IPCK = 1+ICLA
               DO IEL = 1,NCEL
                 W3(IEL) = W3(IEL)
     &                    + PROPCE(IEL,IPPROC(IYFOL(ICLA)))
     &                     *CK(IEL,IPCK)
               ENDDO
             ENDDO
           ENDIF
C
C         Calcul de la longueur caractéristique XLC du domaine de calcul
C
          SF = 0.D0
          VV = 0.D0
C
          DO IFAC = 1,NFABOR
             SF = SF + SQRT(
     &               SURFBO(1,IFAC)**2 +
     &               SURFBO(2,IFAC)**2 +
     &               SURFBO(3,IFAC)**2 )
          ENDDO
          IF (IRANGP.GE.0) THEN
            CALL PARSOM(SF)
C           ===========
          ENDIF
C
          DO IEL = 1,NCEL
             VV = VV + VOLUME(IEL)
          ENDDO
          IF (IRANGP.GE.0) THEN
            CALL PARSOM(VV)
C           ===========
          ENDIF
C
          XLC = 3.6D0 * VV / SF
C
C         Clipping de CK a XKMIN
C
          XKMIN = 1.D0 / XLC

          IOK = 0.D0
          DO IEL = 1,NCEL
            IF (W3(IEL).LT.XKMIN) THEN
              IOK = IOK +1
            ENDIF
          ENDDO
C
C     Arret en fin de pas de temps si epaisseur optique trop grande
          PP = XNP1MX/100.0D0
          IF (DBLE(IOK).GT.PP*DBLE(NCEL)) THEN
             WRITE(NFECRA,1000) XKMIN, DBLE(IOK)/DBLE(NCEL)*100.D0,
     &                          XNP1MX
             ISTPP1 = 1
C             CALL CSEXIT (1)
C            ===========
          ENDIF
C
        ENDIF
C
C -------
C FORMAT
C -------
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT APPROXIMATION P-1 (PPCABS)      ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    LA LONGUEUR OPTIQUE DU MILIEU SEMI-TRANSPARENT          ',/,
     &'@      DOIT AU MOINS ETRE DE L''ORDRE DE L''UNITE POUR ETRE  ',/,
     &'@      DANS LE DOMAINE D''APPLICATION DE L''APPROXIMATION P-1',/,
     &'@    CELA NE SEMBLE PAS ETRE LE CAS ICI.                     ',/,
     &'@                                                            ',/,
     &'@    LE COEFFICIENT D''ABSORPTION MINIMUM POUR ASSURER CETTE ',/,
     &'@      LONGUEUR OPTIQUE EST XKMIN = ',E10.4                   ,/,
     &'@    CETTE VALEUR N''EST PAS ATTEINTE POUR ', E10.4,'%       ',/,
     &'@      DES CELLULES DU MAILLAGE.                             ',/,
     &'@    LE POURCENTAGE DE CELLULES DU MAILLAGE POUR LESQUELLES  ',/,
     &'@      ON ADMET QUE CETTE CONDITION SOIT VIOLEE EST IMPOSE   ',/,
     &'@      PAR DEFAUT OU DANS USINI1 A XNP1MX = ', E10.4,'%      ',/,
     &'@                                                            ',/,
     &'@    Le calcul est interrompu.                               ',/,
     &'@                                                            ',/,
     &'@    Verifier les valeurs du coefficient d''absorption CK    ',/,
     &'@      dans PPCABS, USRAY3 ou Fichier thermochimie.          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      RETURN
C
      END
c@z
