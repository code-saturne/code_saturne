c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS
C
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
                        SUBROUTINE USLAEN
C                       *****************
C     ---------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NVLSTA ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVARL  , IVARL1 , IVARLM , IFLU   , ILPD1  , ICLA   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , STATIS , STATIV , TRACEL ,
     &   RDEVEL , RTUSER , RA     )
C     ---------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC     SOUS-PROGRAMME UTILISATEUR (INTERVENTION NON OBLIGATOIRE)
CFONC
CFONC     POUR LES ECRITURES LISTING ET LE POST-PROCESSING :
CFONC     MOYENNE DES VARIABLES STATISTIQUES VOLUMIQUES LAGRANGIENNES
CFONC     INTERVENTIONS UTILISATEUR POSSIBLES
CFONC
c@fonce
C-----------------------------------------------------------------------
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
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC           !
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR           !
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VARIABLES STAT. LAGRANGIEN !
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
CARGU !  (NFML,NPRFML!    !     !                                      !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC            !
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !                                      !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR            !
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !                                      !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IVARL        !  E !  -> ! NUM DE LA STAT (ENTRE 1 ET NVLSTA)   !
CARGU ! IVARL1       !  E !  -> ! NUM DE LA STAT GLOBALE+GROUPE        !
CARGU !              !    !     !   (MOYENNE OU VARIANCE)              !
CARGU ! IVARLM       !  E !  -> ! NUM DE LA MOYENNE STAT GLOBALE+GROUPE!
CARGU ! IFLU         !  E !  -> ! 0 : MOYENNE DE LA STAT IVARL/IVARL1  !
CARGU !              !    !     ! 1 : VARIANCE DE LA STAT IVARL/IVARL1 !
CARGU ! ILPD1        !  E !  -> ! "POINTEUR" SUR LE POIDS STAT GLOBALE !
CARGU !              !    !     !    OU SUR LE POIDS STAT D'UN GROUPE  !
CARGU ! ICLA         !  E !  -> !   0 : STATISTIQUE GLOBALE            !
CARGU !              !    !  -> ! !=0 : STATISTIQUE POUR LE GROUPE ICLA!
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
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS                !
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
CARGU ! STATIS(NCELET! TR !  -> ! CUMUL DES STATISTIQUES VOLUMIQUES    !
CARGU !   NVLSTA)    !    !     !                                      !
CARGU ! STATIV       ! TR !  -> ! CUMUL POUR LES VARIANCES DES         !
CARGU !(NCELET,      !    !     !    STATISTIQUES VOLUMIQUES           !
CARGU !   NVLSTA-1)  !    !     !                                      !
CARGU ! TRACEL(NCELET! TR ! <-> ! TAB REEL VALEURS CELLULES POST       !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TABLEAU DES STATISTQUES SUR LES      !
CARGU !              !    !     ! PARTICULES                           !
CARGU ! RA(*)        ! TR !  -  ! MACRO TABLEAU REEL                   !
CARGU !______________!____!_____!______________________________________!

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
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS  , NVLSTA
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IVARL , IVARL1 , IVARLM , IFLU , ILPD1 , ICLA
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1) , NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1) , NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE) , ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD) , VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET) , RTP(NCELET,*) , RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*) , COEFB(NFABOR,*)
      DOUBLE PRECISION TRACEL(NCELET)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA)
      DOUBLE PRECISION STATIV(NCELET,NVLSTA-1)
      DOUBLE PRECISION RDEVEL(NRDEVE) , RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA , IEL
      DOUBLE PRECISION AA
C
C***********************************************************************
C=======================================================================
C 0.  PAR DEFAUT, ON CONSIDERE QUE LE SOUS PROGRAMME CI-DESSOUS CONVIENT
C       A L'UTILISATEUR, C'EST-A-DIRE QUE LA MISE EN OEUVRE DU LAGRANGIEN
C       DECLENCHE LA PRODUCTION DE STATISTIQUES STANDARD.
C
C     L'UTILISATEUR N'A PAS A MODIFIER LE PRESENT SOUS-PROGRAMME DANS
C       LES CONDITIONS D'UTILISATION STANDARD.
C     DANS LE CAS OU IL SOUHAITE SORTIR DES STATISTIQUES NON PREVUES,
C       IL DOIT INTERVENIR DANS LA RUBRIQUE NUMERO 2.
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1 . ZONE DE STATISTIQUES STANDARD
C=======================================================================
C
C--> Cas general :
C      Vitesse X des particules : IVARL=ILVX
C      Vitesse Y des particules : IVARL=ILVY
C      Vitesse Z des particules : IVARL=ILVZ
C      Temperature des particules : IVARL=ILTP
C      Diametre des particules : IVARL=ILDP
C      Masse des particules : IVARL=ILMP
C      Temperature des grains de charbon : IVARL=ILHP
C      Masse de Charbon reactif des grains de charbon : IVARL=ILMCH
C      Masse de Coke des grains de charbon : IVARL=ILMCK
C      Diametre du coeur retrecissant des grains de charbon: IVARL=ILMCK
C    sauf fraction volumique (IVARL = ILFV) et
C    somme des poids statistiques (IVARL = ILPD)
C
C
      IF (IVARL.NE.ILFV .AND. IVARL.NE.ILPD) THEN
C
C
C-----> Moyenne
C
C
        IF (IFLU.EQ.0) THEN
C
          DO IEL = 1, NCEL
            IF (STATIS(IEL,ILPD1).GT.SEUIL ) THEN
              TRACEL(IEL) = STATIS(IEL,IVARL1) / STATIS(IEL,ILPD1)
            ELSE
              TRACEL(IEL) = ZERO
            ENDIF
          ENDDO
C
C-----> Variance
C
        ELSE
C
          DO IEL = 1, NCEL
            IF ( STATIS(IEL,ILPD1).GT.SEUIL ) THEN
              AA = STATIS(IEL,IVARLM)/STATIS(IEL,ILPD1)
              TRACEL(IEL) =  STATIV(IEL,IVARL1)/STATIS(IEL,ILPD1)
     &                    -( AA * AA )
              TRACEL(IEL) = SQRT( MAX(ZERO,TRACEL(IEL)))
            ELSE
              TRACEL(IEL) = ZERO
            ENDIF
          ENDDO
C
        ENDIF
C
C--> Fraction volumique (ILFV)
C
      ELSE IF (IVARL.EQ.ILFV) THEN
C
C-----> Moyenne
C
        IF (IFLU.EQ.0) THEN

          DO IEL = 1, NCEL
            IF (STATIS(IEL,ILPD1).GT.SEUIL .AND. NPST.GT.0) THEN
              TRACEL(IEL) = STATIS(IEL,ILFV)
     &                      / (DBLE(NPST) * VOLUME(IEL))
            ELSE IF (STATIS(IEL,ILPD1).GT.SEUIL .AND.
     &                IPLAS.GE.IDSTNT                  ) THEN
              TRACEL(IEL) = STATIS(IEL,ILFV) / VOLUME(IEL)
            ELSE
              TRACEL(IEL) = ZERO
            ENDIF
          ENDDO
C
        ELSE
C
C-----> Variance
C
          DO IEL = 1, NCEL
C
            IF (STATIS(IEL,ILPD1).GT.SEUIL .AND. NPST.GT.0) THEN
C
              AA = STATIS(IEL,IVARLM) / (DBLE(NPST) * VOLUME(IEL))
              TRACEL(IEL) = STATIV(IEL,IVARL1)
     &                 / ( DBLE(NPST) * VOLUME(IEL) * VOLUME(IEL))
     &                 - AA*AA
              TRACEL(IEL) = SQRT( MAX(ZERO,TRACEL(IEL)) )
C
            ELSE IF ( STATIS(IEL,ILPD1).GT.SEUIL .AND.
     &                IPLAS.GE.IDSTNT                  ) THEN
C
              AA =  STATIS(IEL,IVARLM) / VOLUME(IEL)
              TRACEL(IEL) = STATIV(IEL,IVARL1) / VOLUME(IEL)
     &                         - AA*AA
              TRACEL(IEL) = SQRT( MAX(ZERO,TRACEL(IEL)))
C
            ELSE
              TRACEL(IEL) = ZERO
            ENDIF
C
          ENDDO
        ENDIF
C
C--> Somme des poids statistiques
C
      ELSE IF (IVARL.EQ.ILPD) THEN
C
        IF (IFLU .EQ.0) THEN
          DO IEL = 1, NCEL
            TRACEL(IEL) = STATIS(IEL,IVARL1)
          ENDDO
        ELSE
          WRITE(NFECRA,9000) IFLU
          DO IEL = 1, NCEL
            TRACEL(IEL) = ZERO
          ENDDO
        ENDIF
C
      ENDIF
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR DANS LE MODULE LAGRANGIEN (USLAEN)   ',/,
     &'@    *********                                               ',/,
     &'@  IL N''EST PAS POSSIBLE DE CALCULER LA VARIANCE DU POIDS   ',/,
     &'@     STATISTIQUE                                            ',/,
     &'@                                                            ',/,
     &'@  La variance du poids statistique a ete demande            ',/,
     &'@    dans USLAEN (IVARL=',   I10,' et IFLU=',  I10,').       ',/,
     &'@                                                            ',/,
     &'@  L''appel au sous-programme USLAEN doit etre verifie.      ',/,
     &'@                                                            ',/,
     &'@  Le calcul continue.                                       ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C=======================================================================
C 2. ZONE D'INTERVENTION UTILISATEUR
C=======================================================================
C
C    --------------------------------------------------
C    EXEMPLE 1 : STATISTIQUE CALCULEE DANS USLAST.F ET
C                STOCKEE DANS LE TABLEAU STATIS
C    --------------------------------------------------
C
        IF (NVLSTS.GT.0) THEN
C
          IF (IVARL.EQ.ILVU(1)) THEN
C
C-----> Moyenne pour la concentration massique
C
            IF (IFLU.EQ.0) THEN
C
              DO IEL = 1, NCEL
                IF ( STATIS(IEL,ILPD1).GT.SEUIL .AND. NPST.GT.0 ) THEN
                  TRACEL(IEL) = STATIS(IEL,IVARL1)
     &                        / ( DBLE(NPST) *RO0(ILPHAS) *VOLUME(IEL) )
                ELSE IF ( STATIS(IEL,ILPD1).GT.SEUIL .AND.
     &                  IPLAS.GE.IDSTNT                  ) THEN
                  TRACEL(IEL) = STATIS(IEL,IVARL1)
     &                        / ( RO0(ILPHAS) *VOLUME(IEL) )
                ELSE
                  TRACEL(IEL) = ZERO
                ENDIF
              ENDDO
C
            ELSE
C
C-----> Variance de la concentration massique
C
              DO IEL = 1, NCEL
                IF (STATIS(IEL,ILPD1).GT.SEUIL) THEN
                  AA = STATIS(IEL,IVARLM)/STATIS(IEL,ILPD1)
                  TRACEL(IEL) = STATIV(IEL,IVARL1)/STATIS(IEL,ILPD1)
     &                        -( AA * AA)
                  TRACEL(IEL) = SQRT( MAX(ZERO,TRACEL(IEL)))
                ELSE
                  TRACEL(IEL) = ZERO
                ENDIF
C
              ENDDO
C
            ENDIF
C
          ENDIF
C
        ENDIF
C
C=======================================================================
C
C
      END
c@z
