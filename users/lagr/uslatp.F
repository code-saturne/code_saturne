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
                        SUBROUTINE USLATP
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   ,
     &   NPRFML , NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   NUMPT  , ITEPA  , IDEVEL , ITUSER , IA     ,
     &   REP    , UVWR   , ROMF   , ROMP   , XNUL   , TAUP   ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , PROPCE , PROPFA , PROPFB ,
     &   ETTP   , ETTPA  , TEPA   ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC     SOUS-PROGRAMME UTILISATEUR (INTERVENTION NON OBLIGATOIRE)
CFONC
CFONC     MODIFICATION DU CALCUL DU TEMPS DE RELAXATION DES PARTICULES
CFONC       EN FONCTION DE LA FORMULATION CHOISIE POUR LE COEFFICIENT
CFONC       DE TRAINEE
CFONC
CFONC
CFONC             rho             4 d
CFONC                p               p
CFONC       Tau = ---- --------------------------------
CFONC          p
CFONC             rho   3 C     | U [X (t),t] - V (t) |
CFONC                f     drag    f  p          p
CFONC
CFONC      Tau  : TEMPS DE RELAXATION PARTICULAIRE (VALEUR A CALCULER)
CFONC         p
CFONC
CFONC      rho  : MASSE VOLUMIQUE DU FLUIDE
CFONC         p
CFONC
CFONC      rho  : MASSE VOLUMIQUE DE LA PARTICULE
CFONC         f
CFONC
CFONC      C     : COEFFICIENT DE TRAINEE
CFONC       drag
CFONC
CFONC      d    : DIAMETRE DE LA PARTICULE
CFONC       p
CFONC
CFONC      U [X (t),t] : VITESSE INSTANTANEE DU FLUIDE VU
CFONC       f  p
CFONC
CFONC      V (t) : VITESSE DE LA PARTICULE
CFONC       p
CFONC
CFONC     CE SOUS PROGRAMME EST APPELE DANS UNE BOUCLE SUR
CFONC       LES PARTICULES : ATTENTION DE NE PAS TROP LE CHARGER
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
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! NUMPT        ! E  !  -> ! NUMERO DE LA PARTICULE COURANTE      !
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TE !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! REP          ! R  !  -> ! NOMBRE DE REYNOLDS PARTICULAIRE      !
CARGU !              !    !     ! REP = UVWR * ETTP(NUMPT,JDP) / XNUL  !
CARGU ! UVWR         ! R  !  -> ! VITESSE RELATIVE DE LA PARTICULE     !
CARGU !              !    !     ! UVWR= |VIT FLUIDE VU - VIT PARTICULE|!
CARGU ! ROMF         ! R  !  -> ! MASSE VOLUMIQUE DU FLUIDE A LA       !
CARGU !              !    !     ! POSITION DE LA PARTICULE             !
CARGU ! ROMP         ! R  !  -> ! MASSE VOLUMIQUE DE LA PARTICULE      !
CARGU ! XNUL         ! R  !  -> ! VISCOSITE CINEMATIQUE DU FLUIDE A LA !
CARGU !              !    !     ! POSITION DE LA PARTICULE             !
CARGU ! TAUP         ! R  ! <-  ! TEMPS DE RELAXATION PARTICULAIRE     !
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
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
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
C
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
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "cpincl.h"
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
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          NUMPT
      INTEGER          ITEPA(NBPMAX,NIVEP)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION REP    , UVWR   , ROMF   , ROMP   , XNUL  , TAUP
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD) , VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET) , RTP(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IP
      DOUBLE PRECISION FDR
C
C VARIABLES LOCALES UTILISATEUR
C
      DOUBLE PRECISION CD1 , CD2 , D2
      DOUBLE PRECISION REC1, REC2, REC3, REC4
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C
      IF(1.EQ.1) RETURN
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C 0.  GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. INITIALISATIONS
C=======================================================================
C
      IP = NUMPT
C
C=======================================================================
C 2. TEMPS DE RELAXATION AVEC LE COEFFICIENT DE TRAINEE STANDARD
C=======================================================================
C
C   Cet exemple est desactive, il donne le temps de relaxation
C   standard a titre indicatif.
C
      IF (1.EQ.0) THEN
C
        CD1  = 0.15D0
        CD2  = 0.687D0
C
        IF (REP.LE.1000) THEN
            D2 = ETTP(IP,JDP) * ETTP(IP,JDP)
            FDR = 18.D0 * XNUL * (1.D0 + CD1 * REP**CD2) / D2
        ELSE
            FDR = (0.44D0 * 3.D0 / 4.D0) * UVWR / ETTP(IP,JDP)
        ENDIF
C
        TAUP = ROMP / ROMF / FDR
C
      ENDIF
C
C=======================================================================
C 3. CALCUL DU TEMPS DE RELAXATION AVEC LE COEFFICIENT DE TRAINEE DE
C    MORSI S.A. et ALEXANDER A.J.,
C    Journal of Fluid Mechanics Vol.55, pp 193-208 (1972)
C=======================================================================
C
      REC1 =  0.1D0
      REC2 =  1.0D0
      REC3 =  10.D0
      REC4 = 200.D0
C
      D2 = ETTP(IP,JDP) * ETTP(IP,JDP)
C
      IF ( REP.LE.REC1 ) THEN
        FDR = 18.D0 * XNUL / D2
C
      ELSE IF ( REP.LE.REC2 ) THEN
        FDR = 3.D0/4.D0 * XNUL / D2
     &      * (22.73D0 + 0.0903D0/REP + 3.69D0*REP )
C
      ELSE IF ( REP.LE.REC3 ) THEN
        FDR = 3.D0/4.D0 * XNUL / D2
     &      * (29.1667D0 - 3.8889D0/REP + 1.222D0*REP)
C
      ELSE IF ( REP.LE.REC4 ) THEN
          FDR = 18.D0*XNUL/D2 *(1.D0 + 0.15D0*REP**0.687D0)
C
      ELSE
         FDR = (0.44D0 * 3.D0 / 4.D0) * UVWR / ETTP(IP,JDP)
      ENDIF
C
      TAUP = ROMP / ROMF / FDR
C
C
C**********************************************************************
C
C--------
C FORMATS
C--------
C
C
C----
C FIN
C----
C
      END
c@z
