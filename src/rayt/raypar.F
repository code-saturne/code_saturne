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
                        SUBROUTINE RAYPAR
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , ISOTHP , IZFRAP ,
     &   IDEVEL , ITUSER , IA     ,
     &   TMIN   , TMAX   , TX     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB , RCODCL ,
     &   COEFA  , COEFB  ,
     &   RDEVEL , RTUSER ,
     &   TPAROP , QINCIP , TEXTP  , TINTP  ,
     &   XLAMP  , EPAP   , EPSP   ,
     &   HFCONP , FLCONP , TEMPKP ,
C
     &   RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
CFONC    -----------------------------------------
CFONC
CFONC    CALCUL DES TEMPERATURES DE PAROIS AVEC UN BILAN DE FLUX
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
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE LA PHASE                   !
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
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! ICODCL       ! TE ! <-  ! CODE DE CONDITION LIMITES AUX FACES  !
CARGU !  (NFABOR,NVAR!    !     !  DE BORD                             !
CARGU !              !    !     ! = 1   -> DIRICHLET                   !
CARGU !              !    !     ! = 3   -> DENSITE DE FLUX             !
CARGU !              !    !     ! = 4   -> GLISSEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 5   -> FROTTEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 6   -> RUGOSITE ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 9   -> ENTREE/SORTIE LIBRE (VITESSE!
CARGU !              !    !     !  ENTRANTE EVENTUELLE     BLOQUEE     !
CARGU ! ISOTHP(NFABOR! TE !  -> ! LISTE DES FRONTIERES ISOTHERMES      !
CARGU ! IZFRAP(NFABOR! TE !  -> ! NUMERO DE ZONE DES FACES DE BORD     !
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
CARGU ! SURFBO       ! TR !  -> ! VECTEUR SURFACE DES FACES DE BORD    !
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
CARGU ! RCODCL       ! TR ! <-  ! VALEUR DES CONDITIONS AUX LIMITES    !
CARGU !  (NFABOR,NVAR!    !     !  AUX FACES DE BORD                   !
CARGU !              !    !     ! RCODCL(1) = VALEUR DU DIRICHLET      !
CARGU !              !    !     ! RCODCL(2) = VALEUR DU COEF. D'ECHANGE!
CARGU !              !    !     !  EXT. (INFINIE SI PAS D'ECHANGE)     !
CARGU !              !    !     ! RCODCL(3) = VALEUR DE LA DENSITE DE  !
CARGU !              !    !     !  FLUX (NEGATIF SI GAIN) W/m2 OU      !
CARGU !              !    !     !  HAUTEUR DE RUGOSITE (m) si ICODCL=6 !
CARGU !              !    !     ! POUR LES VITESSES (VISTL+VISCT)*GRADU!
CARGU !              !    !     ! POUR LA PRESSION             DT*GRADP!
CARGU !              !    !     ! POUR LES SCALAIRES                   !
CARGU !              !    !     !        CP*(VISCLS+VISCT/SIGMAS)*GRADT!
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! TPAROP(NFABOR! TR ! <-  ! TEMPERATURE DE PAROI EN KELVIN       !
CARGU ! QINCIP(NFABOR! TR !  -> ! DENSITE DE FLUX RADIATIF AUX BORDS   !
CARGU ! TEXTP(NFABOR)! TR !  -> ! TEMPERATURE DE BORD EXTERNE          !
CARGU !              !    !     ! EN DEGRES CELCIUS                    !
CARGU ! TINTP(NFABOR)! TR !  -> ! TEMPERATURE DE BORD INTERNE          !
CARGU !              !    !     ! EN DEGRES CELCIUS                    !
CARGU ! XLAMP(NFABOR)! TR !  -> ! COEFFICIENT DE CONDUCTIVITE THERMIQUE!
CARGU !              !    !     ! DES FACETTES DE PAROI (W/m/K)        !
CARGU ! EPAP(NFABOR) ! TR !  -> ! EPAISSEUR DES FACETTES DE PAROI (m)  !
CARGU ! EPSP(NFABOR) ! TR !  -> ! EMISSIVITE DES FACETTES DE BORD      !
CARGU ! HFCONP(NFABOR! TR !  -> ! COEFFICIENT D'ECHANGE FLUIDE AUX     !
CARGU !              !    !     ! FACES DE BORD                        !
CARGU ! FLCONP(NFABOR! TR !  -> ! DENSITE DE FLUX CONVECTIF AUX FACES  !
CARGU ! TEMPKP(NCELET! TR !  -> ! TEMPERATURE EN KELVIN                !
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
      INCLUDE "parall.h"
      INCLUDE "radiat.h"
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
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML) , ITYPFB(NFABOR)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      INTEGER          ISOTHP(NFABOR), IZFRAP(NFABOR)
      INTEGER          ICODCL(NFABOR,NVAR)
C
      DOUBLE PRECISION TMIN , TMAX , TX
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION RCODCL(NFABOR,NVAR,3)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
C
      DOUBLE PRECISION TPAROP(NFABOR), QINCIP(NFABOR)
      DOUBLE PRECISION TEXTP(NFABOR), TINTP(NFABOR)
      DOUBLE PRECISION XLAMP(NFABOR), EPAP(NFABOR), EPSP(NFABOR)
      DOUBLE PRECISION HFCONP(NFABOR) , FLCONP(NFABOR)
      DOUBLE PRECISION TEMPKP(NCELET)
C
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
C
      INTEGER          IVART, IFAC, IEL, IZONE
      INTEGER          IFACMX, IFACMN
      INTEGER          N1MIN,N1MAX,NRELAX,NMOINS,NPLUS
      INTEGER          IITPIM,IIPGRN,IIPREF,IIFGRN,IIFREF
      INTEGER          INDTP(NOZRDM)
      INTEGER          NBRVAL, INDTPM
      INTEGER          NB1INT  ,NB2INT  ,NBRRDP
      PARAMETER       (NB1INT=5,NB2INT=5,NBRRDP=5)
      INTEGER          INTTM1(NB1INT),INTTM2(NB2INT)
C
      DOUBLE PRECISION ESL,EPP,SIGT3,SIGT4
      DOUBLE PRECISION TPMIN,TPMAX,RAPP,RAPMAX,ABRAPP
      DOUBLE PRECISION QCMAX,QRMAX
      DOUBLE PRECISION QCMIN,QRMIN
      DOUBLE PRECISION QRAYT,QCONV,QINCI,DETEP
      DOUBLE PRECISION SURFBN,UND0,TP4
      DOUBLE PRECISION XTPMAX,YTPMAX,ZTPMAX,XTPMIN,YTPMIN,ZTPMIN
      DOUBLE PRECISION TZOMAX(NOZRDM),TZOMIN(NOZRDM),TZOMOY(NOZRDM)
      DOUBLE PRECISION FLUNET(NOZRDM),RADIOS(NOZRDM),SURFT(NOZRDM)
      DOUBLE PRECISION RDPTMP(NBRRDP)
C
C
C***********************************************************************
C
C=======================================================================
C 0. GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
C
      UND0 = 1.D0
C
      IVART = ISCA(ISCALT(IPHAS))
C
C
      TPMAX = -GRAND
      TPMIN =  GRAND
C
      RAPMAX = 0.D0
      N1MIN  = 0
      N1MAX  = 0
      NRELAX = 0
      NMOINS = 0
      NPLUS  = 0
      IITPIM = 0
      IIPGRN = 0
      IIPREF = 0
      IIFGRN = 0
      IIFREF = 0
C
      IFACMX = 0
      IFACMN = 0
C
C=======================================================================
C 2. Decompte du nombre de couleurs differentes de facettes de bord
C    On suppose qu'on en a au moins une
C=======================================================================
C
C
      DO IZONE = 1, NOZRDM
        INDTP (IZONE) =      0
        TZOMAX(IZONE) = -GRAND
        TZOMIN(IZONE) =  GRAND
      ENDDO
C
C=======================================================================
C 3. Calcul des temperatures de paroi
C=======================================================================
C
C     3.1) parois isothermes
C     ----------------------
C
C
      DO IFAC = 1, NFABOR
        IF (ISOTHP(IFAC).EQ.ITPIMP) THEN
          IZONE = IZFRAP(IFAC)
C
C      Calcul
          TPAROP(IFAC)  = TINTP(IFAC)
C
C      Max-Min
          QCONV  = FLCONP(IFAC)
          QINCI  = QINCIP(IFAC)
          SIGT4  = STEPHN*TPAROP(IFAC)**4
          EPP    = EPSP(IFAC)
          QRAYT  = EPP *( QINCI - SIGT4 )
C
          IF (TPMAX.LE.TPAROP(IFAC)) THEN
            IFACMX = IFAC
            TPMAX = TPAROP(IFAC)
            QCMAX = QCONV
            QRMAX = QRAYT
          ENDIF
C
          IF (TPMIN.GE.TPAROP(IFAC)) THEN
            IFACMN = IFAC
            TPMIN = TPAROP(IFAC)
            QCMIN = QCONV
            QRMIN = QRAYT
          ENDIF
C
          TZOMAX(IZONE) = MAX(TZOMAX(IZONE),TPAROP(IFAC))
          TZOMIN(IZONE) = MIN(TZOMIN(IZONE),TPAROP(IFAC))
C
C      Reperage pour impression
          IITPIM = 1
          INDTP(IZONE) = ITPIMP
C
C
        ENDIF
      ENDDO
C
C
C     3.2) parois grises ou noires (non reflechissantes)
C     --------------------------------------------------
C
      DO IFAC = 1, NFABOR
        IF (ISOTHP(IFAC).EQ.IPGRNO) THEN
          IZONE = IZFRAP(IFAC)
C
C       Calcul
          ESL    = EPAP(IFAC) /XLAMP(IFAC)
          QCONV  = FLCONP(IFAC)
          QINCI  = QINCIP(IFAC)
          EPP    = EPSP(IFAC)
C
          SIGT3  = STEPHN*TPAROP(IFAC)**3
          QRAYT  = EPP *( QINCI -SIGT3*TPAROP(IFAC))
C
          DETEP  = ( ESL*(QCONV+QRAYT) - (TPAROP(IFAC)-TEXTP(IFAC)) )
     &           / ( 1.D0+ 4.D0*ESL*EPP*SIGT3 + ESL*HFCONP(IFAC) )
C
          RAPP   = DETEP /TPAROP(IFAC)
          ABRAPP = ABS(RAPP)
C
C       Relaxation
          IF (ABRAPP.GE.TX) THEN
            NRELAX = NRELAX +1
            TPAROP(IFAC) = TPAROP(IFAC) * (1.D0+ TX*SIGN(UND0,RAPP))
          ELSE
            TPAROP(IFAC) = TPAROP(IFAC) + DETEP
          ENDIF
C
          RAPMAX = MAX(RAPMAX,ABRAPP)
C
          IF (RAPP.LE.0.D0) THEN
            NMOINS = NMOINS+1
          ELSE
            NPLUS  = NPLUS +1
          ENDIF
C
C       Clipping
          IF (TPAROP(IFAC).LT.TMIN) THEN
            N1MIN = N1MIN +1
            TPAROP(IFAC) = TMIN
          ENDIF
          IF (TPAROP(IFAC).GT.TMAX) THEN
            N1MAX = N1MAX +1
            TPAROP(IFAC) = TMAX
          ENDIF
C
C       Max-Min
          IF (TPMAX.LE.TPAROP(IFAC)) THEN
            IFACMX = IFAC
            TPMAX = TPAROP(IFAC)
            QCMAX = QCONV
            QRMAX = QRAYT
          ENDIF
C
          IF (TPMIN.GE.TPAROP(IFAC)) THEN
            IFACMN = IFAC
            TPMIN = TPAROP(IFAC)
            QCMIN = QCONV
            QRMIN = QRAYT
          ENDIF
C
          TZOMAX(IZONE) = MAX(TZOMAX(IZONE),TPAROP(IFAC))
          TZOMIN(IZONE) = MIN(TZOMIN(IZONE),TPAROP(IFAC))
C
C       Reperage pour impression
          IIPGRN = 1
          INDTP(IZONE) = IPGRNO
C
        ENDIF
      ENDDO
C
C
C     3.3) parois reflechissantes
C     ---------------------------
C
C
      DO IFAC = 1, NFABOR
        IF (ISOTHP(IFAC).EQ.IPREFL) THEN
          IZONE = IZFRAP(IFAC)
C
C       Calcul
          ESL    = EPAP(IFAC) /XLAMP(IFAC)
          QCONV  = FLCONP(IFAC)
C
          DETEP  = ( ESL*QCONV - (TPAROP(IFAC)-TEXTP(IFAC)) )
     &           / ( 1.D0 + ESL*HFCONP(IFAC) )
C
          RAPP = DETEP /TPAROP (IFAC)
          ABRAPP = ABS(RAPP)
C
C       Relaxation
          IF (ABRAPP.GE.TX) THEN
            NRELAX = NRELAX +1
            TPAROP(IFAC) = TPAROP(IFAC) * (1.D0+ TX*SIGN(UND0,RAPP))
          ELSE
            TPAROP(IFAC) = TPAROP(IFAC) + DETEP
          ENDIF
C
          RAPMAX = MAX(RAPMAX,ABRAPP)
C
          IF (RAPP.LE.0.D0) THEN
            NMOINS = NMOINS+1
          ELSE
            NPLUS  = NPLUS +1
          ENDIF
C
C       Clipping
          IF (TPAROP(IFAC).LT.TMIN) THEN
            N1MIN = N1MIN +1
            TPAROP(IFAC) = TMIN
          ENDIF
          IF (TPAROP(IFAC).GT.TMAX) THEN
            N1MAX = N1MAX +1
            TPAROP(IFAC) = TMAX
          ENDIF
C
C       Max-Min
          IF (TPMAX.LE.TPAROP(IFAC)) THEN
            IFACMX = IFAC
            TPMAX = TPAROP(IFAC)
            QCMAX = QCONV
            QRMAX = 0.D0
          ENDIF
C
          IF (TPMIN.GE.TPAROP(IFAC)) THEN
            IFACMN = IFAC
            TPMIN = TPAROP(IFAC)
            QCMIN = QCONV
            QRMIN = 0.D0
          ENDIF
C
          TZOMAX(IZONE) = MAX(TZOMAX(IZONE),TPAROP(IFAC))
          TZOMIN(IZONE) = MIN(TZOMIN(IZONE),TPAROP(IFAC))
C
C       Reperage pour impression
          IIPREF = 1
          INDTP(IZONE) = IPREFL
C
        ENDIF
      ENDDO
C
C
C
C     3.4) parois a flux de conduction impose non reflechissante
C          si le flux impose est nul, la paroi est adiabatique
C          (la transmition de chaleur par rayonnement est
C          equilibree avec celle de la convection)
C     ----------------------------------------------------------
C
      DO IFAC = 1, NFABOR
        IF (ISOTHP(IFAC).EQ.IFGRNO) THEN
          IZONE = IZFRAP(IFAC)
C
C       Calcul
          QCONV  = FLCONP(IFAC)
          QINCI  = QINCIP(IFAC)
          EPP    = EPSP(IFAC)
C
          SIGT3  = STEPHN*TPAROP(IFAC)**3
          QRAYT  = EPP *( QINCI -SIGT3*TPAROP(IFAC))
C
          DETEP  = ( QCONV + QRAYT - RCODCL(IFAC,IVART,3) )
     &           / ( 4.D0*EPP*SIGT3 + HFCONP(IFAC) )
C
          RAPP   = DETEP /TPAROP (IFAC)
          ABRAPP = ABS(RAPP)
C
C       Relaxation
          IF (ABRAPP.GE.TX) THEN
            NRELAX = NRELAX +1
            TPAROP(IFAC) = TPAROP(IFAC) * (1.D0+ TX*SIGN(UND0,RAPP))
          ELSE
            TPAROP(IFAC) = TPAROP(IFAC) + DETEP
          ENDIF
C
          RAPMAX = MAX(RAPMAX,ABRAPP)
C
          IF (RAPP.LE.0.D0) THEN
            NMOINS = NMOINS+1
          ELSE
            NPLUS  = NPLUS +1
          ENDIF
C
C       Clipping
          IF (TPAROP(IFAC).LT.TMIN) THEN
            N1MIN = N1MIN +1
            TPAROP(IFAC) = TMIN
          ENDIF
          IF (TPAROP(IFAC).GT.TMAX) THEN
            N1MAX = N1MAX +1
            TPAROP(IFAC) = TMAX
          ENDIF
C
C       Max-Min
          IF (TPMAX.LE.TPAROP(IFAC)) THEN
            IFACMX = IFAC
            TPMAX = TPAROP(IFAC)
            QCMAX = QCONV
            QRMAX = QRAYT
          ENDIF
C
          IF (TPMIN.GE.TPAROP(IFAC)) THEN
            IFACMN = IFAC
            TPMIN = TPAROP(IFAC)
            QCMIN = QCONV
            QRMIN = QRAYT
          ENDIF
C
          TZOMAX(IZONE) = MAX(TZOMAX(IZONE),TPAROP(IFAC))
          TZOMIN(IZONE) = MIN(TZOMIN(IZONE),TPAROP(IFAC))
C
C       Reperage pour impression
          IIFGRN = 1
          INDTP(IZONE) = IFGRNO
C
        ENDIF
      ENDDO
C
C
C     3.5) parois a flux de conduction imposee et reflechissante
C          equivalent a impose un flux total au fluide
C     -----------------------------------------------------------
C
C
      DO IFAC = 1, NFABOR
        IF (ISOTHP(IFAC).EQ.IFREFL) THEN
          IZONE = IZFRAP(IFAC)
C
C       Calcul
          IEL = IFABOR(IFAC)
C
          TPAROP(IFAC)= HFCONP(IFAC)*TEMPKP(IEL)-RCODCL(IFAC,IVART,3)
          TPAROP(IFAC)= TPAROP(IFAC)/MAX(HFCONP(IFAC),EPZERO)
C
C       Clipping
          IF (TPAROP(IFAC).LT.TMIN) THEN
            N1MIN = N1MIN +1
            TPAROP(IFAC) = TMIN
          ENDIF
          IF (TPAROP(IFAC).GT.TMAX) THEN
            N1MAX = N1MAX +1
            TPAROP(IFAC) = TMAX
          ENDIF
C
C       Max-Min
          QCONV  = FLCONP(IFAC)
          QRAYT  = 0.D0
C
          IF (TPMAX.LE.TPAROP(IFAC)) THEN
            IFACMX = IFAC
            TPMAX = TPAROP(IFAC)
            QCMAX = QCONV
            QRMAX = QRAYT
          ENDIF
C
          IF (TPMIN.GE.TPAROP(IFAC)) THEN
            IFACMN = IFAC
            TPMIN = TPAROP(IFAC)
            QCMIN = QCONV
            QRMIN = QRAYT
          ENDIF
C
          TZOMAX(IZONE) = MAX(TZOMAX(IZONE),TPAROP(IFAC))
          TZOMIN(IZONE) = MIN(TZOMIN(IZONE),TPAROP(IFAC))
C
C       Reperage pour impression
          IIFREF = 1
          INDTP(IZONE) = IFREFL
C
        ENDIF
      ENDDO
C
C
C=======================================================================
C 4. IMPRESSIONS
C=======================================================================
C
C
C
C Test pour savoir s'il y a des zones de paroi
C  (donc si des impressions sont necessaires)
C
      IF (IRANGP.GE.0) THEN
        CALL PARIMX(NOZARM(IPHAS),INDTP )
C       ===========
      ENDIF
C
      INDTPM = 0
      DO IZONE = 1, NOZRDM
        IF (INDTP(IZONE).NE.0) THEN
          INDTPM = 1
        ENDIF
      ENDDO
C
C Si il y a des parois
C
      IF(INDTPM.GT.0) THEN
C
C   Si on veut imprimer en niveau 1 au moins
C
        IF (IIMPAR.GE.1) THEN
C
C      Calcul de TZOMOY FLUNET RADIOS SURFT
C
          DO IZONE = 1, NOZRDM
            TZOMOY(IZONE) = ZERO
            FLUNET(IZONE) = ZERO
            RADIOS(IZONE) = ZERO
            SURFT (IZONE) = ZERO
          ENDDO
          DO IFAC = 1, NFABOR
            SURFBN = RA(ISRFBN-1+IFAC)
            IZONE = IZFRAP(IFAC)
            IF (INDTP(IZONE).NE.0) THEN
              TP4 = TPAROP(IFAC)**4
              TZOMOY(IZONE)= TZOMOY(IZONE) + TPAROP(IFAC)*SURFBN
              FLUNET(IZONE)= FLUNET(IZONE)
     &             + EPSP(IFAC) *(QINCIP(IFAC) -STEPHN*TP4 )*SURFBN
              RADIOS(IZONE)= RADIOS(IZONE)
     &             - ( EPSP(IFAC)       *STEPHN*TP4
     &             +(1.D0-EPSP(IFAC))*QINCIP(IFAC)      )*SURFBN
              SURFT (IZONE) = SURFT(IZONE) + SURFBN
            ENDIF
          ENDDO
C
          IF(IRANGP.GE.0) THEN
            CALL PARRSM(NOZARM(IPHAS),TZOMOY)
C           ===========
            CALL PARRSM(NOZARM(IPHAS),FLUNET)
C           ===========
            CALL PARRSM(NOZARM(IPHAS),RADIOS)
C           ===========
            CALL PARRSM(NOZARM(IPHAS),SURFT )
C           ===========
          ENDIF
C
          DO IZONE = 1, NOZRDM
            IF (INDTP(IZONE).NE.0) THEN
              TZOMOY(IZONE) = TZOMOY(IZONE)/SURFT(IZONE)
              RADIOS(IZONE) = RADIOS(IZONE)/SURFT(IZONE)
            ENDIF
          ENDDO
C
C
C      Determination de la TPMAX TPMIN et des grandeurs associees
C
          IF(IFACMX.GT.0) THEN
            XTPMAX = XYZCEN(1,IFABOR(IFACMX))
            YTPMAX = XYZCEN(2,IFABOR(IFACMX))
            ZTPMAX = XYZCEN(3,IFABOR(IFACMX))
          ELSE
            XTPMAX = 0.D0
            YTPMAX = 0.D0
            ZTPMAX = 0.D0
          ENDIF
          IF(IFACMN.GT.0) THEN
            XTPMIN = XYZCEN(1,IFABOR(IFACMN))
            YTPMIN = XYZCEN(2,IFABOR(IFACMN))
            ZTPMIN = XYZCEN(3,IFABOR(IFACMN))
          ELSE
            XTPMIN = 0.D0
            YTPMIN = 0.D0
            ZTPMIN = 0.D0
          ENDIF
C
          IF(IRANGP.GE.0) THEN
C
            RDPTMP(1) = XTPMAX
            RDPTMP(2) = YTPMAX
            RDPTMP(3) = ZTPMAX
            RDPTMP(4) = QCMAX
            RDPTMP(5) = QRMAX
            NBRVAL = NBRRDP
            CALL PARMXL(NBRVAL,TPMAX,RDPTMP)
C           ===========
            XTPMAX =  RDPTMP(1)
            YTPMAX =  RDPTMP(2)
            ZTPMAX =  RDPTMP(3)
            QCMAX  =  RDPTMP(4)
            QRMAX  =  RDPTMP(5)
C
            RDPTMP(1) = XTPMIN
            RDPTMP(2) = YTPMIN
            RDPTMP(3) = ZTPMIN
            RDPTMP(4) = QCMIN
            RDPTMP(5) = QRMIN
            NBRVAL = NBRRDP
            CALL PARMNL(NBRVAL,TPMIN,RDPTMP)
C           ===========
            XTPMIN =  RDPTMP(1)
            YTPMIN =  RDPTMP(2)
            ZTPMIN =  RDPTMP(3)
            QCMIN  =  RDPTMP(4)
            QRMIN  =  RDPTMP(5)
C
          ENDIF
C
C
C      Determination des compteurs et autres
C
          IF(IRANGP.GE.0) THEN
C
            CALL PARMAX(RAPMAX)
C           ===========
C
            INTTM2(1) = NMOINS
            INTTM2(2) = NPLUS
            INTTM2(3) = N1MIN
            INTTM2(4) = N1MAX
            INTTM2(5) = NRELAX
            NBRVAL = NB2INT
            CALL PARISM(NBRVAL,INTTM2)
C           ===========
            NMOINS = INTTM2(1)
            NPLUS  = INTTM2(2)
            N1MIN  = INTTM2(3)
            N1MAX  = INTTM2(4)
            NRELAX = INTTM2(5)
C
            CALL PARRMX(NOZARM(IPHAS),TZOMAX)
C           ===========
            CALL PARRMN(NOZARM(IPHAS),TZOMIN)
C           ===========
C
            INTTM1(1) = IITPIM
            INTTM1(2) = IIPGRN
            INTTM1(3) = IIPREF
            INTTM1(4) = IIFGRN
            INTTM1(5) = IIFREF
            NBRVAL    = NB1INT
            CALL PARIMX(NBRVAL,INTTM1)
C           ===========
            IITPIM = INTTM1(1)
            IIPGRN = INTTM1(2)
            IIPREF = INTTM1(3)
            IIFGRN = INTTM1(4)
            IIFREF = INTTM1(5)
C
          ENDIF
C
C
C      Impressions
C
          WRITE(NFECRA,1000)
          WRITE(NFECRA,1010)
C
          IF (NRELAX.GT.0) THEN
            WRITE(NFECRA,2010) TX*100.D0 ,NRELAX
            WRITE(NFECRA,1010)
          ENDIF
C
          IF ( N1MIN.GT.0 .OR. N1MAX.GT.0 ) THEN
            WRITE(NFECRA,2020)
            WRITE(NFECRA,2030) N1MIN
            WRITE(NFECRA,2040) N1MAX
            WRITE(NFECRA,1010)
          ENDIF
C
          IF (RAPMAX.GT.0 .OR. NMOINS.GT.0 .OR. NPLUS.GT.0 ) THEN
            WRITE(NFECRA,2050) RAPMAX * 100.D0
            WRITE(NFECRA,2060) NMOINS
            WRITE(NFECRA,2070) NPLUS
            WRITE(NFECRA,1010)
          ENDIF
C
          IF (IITPIM.EQ.1) THEN
            WRITE(NFECRA,3020)
            DO IZONE = 1, NOZRDM
              IF (INDTP(IZONE).EQ.ITPIMP) THEN
                WRITE(NFECRA,3000) IZONE,
     &               TZOMAX(IZONE)-TKELVI,
     &               TZOMIN(IZONE)-TKELVI,
     &               TZOMOY(IZONE)-TKELVI,
     &               FLUNET(IZONE)
              ENDIF
            ENDDO
            WRITE(NFECRA,1010)
          ENDIF
C
          IF (IIPGRN.EQ.1) THEN
            WRITE(NFECRA,3010)
            DO IZONE = 1, NOZRDM
              IF (INDTP(IZONE).EQ.IPGRNO) THEN
                WRITE(NFECRA,3000) IZONE,
     &               TZOMAX(IZONE)-TKELVI,
     &               TZOMIN(IZONE)-TKELVI,
     &               TZOMOY(IZONE)-TKELVI,
     &               FLUNET(IZONE)
              ENDIF
            ENDDO
            WRITE(NFECRA,1010)
          ENDIF
C
          IF (IIPREF.EQ.1) THEN
            WRITE(NFECRA,3030)
            DO IZONE = 1, NOZRDM
              IF (INDTP(IZONE).EQ.IPREFL) THEN
                WRITE(NFECRA,3000) IZONE,
     &               TZOMAX(IZONE)-TKELVI,
     &               TZOMIN(IZONE)-TKELVI,
     &               TZOMOY(IZONE)-TKELVI,
     &               FLUNET(IZONE)
              ENDIF
            ENDDO
            WRITE(NFECRA,1010)
          ENDIF
C
          IF (IIFGRN.EQ.1) THEN
            WRITE(NFECRA,3040)
            DO IZONE = 1, NOZRDM
              IF (INDTP(IZONE).EQ.IFGRNO) THEN
                WRITE(NFECRA,3000) IZONE,
     &               TZOMAX(IZONE)-TKELVI,
     &               TZOMIN(IZONE)-TKELVI,
     &               TZOMOY(IZONE)-TKELVI,
     &               FLUNET(IZONE)
              ENDIF
            ENDDO
            WRITE(NFECRA,1010)
          ENDIF
C
          IF (IIFREF.EQ.1) THEN
            WRITE(NFECRA,3050)
            DO IZONE = 1, NOZRDM
              IF (INDTP(IZONE).EQ.IFREFL) THEN
                WRITE(NFECRA,3000) IZONE,
     &               TZOMAX(IZONE)-TKELVI,
     &               TZOMIN(IZONE)-TKELVI,
     &               TZOMOY(IZONE)-TKELVI,
     &               FLUNET(IZONE)
              ENDIF
            ENDDO
            WRITE(NFECRA,1010)
          ENDIF
C
        ENDIF
C
C
C   Si on veut imprimer EN PLUS en niveau 2
C                       =======
C
        IF (IIMPAR.GE.2) THEN
C
          WRITE(NFECRA,5010) TPMAX-TKELVI
          WRITE(NFECRA,5020) XTPMAX, YTPMAX, ZTPMAX
          WRITE(NFECRA,5030) QCMAX
          WRITE(NFECRA,5040) QRMAX
          WRITE(NFECRA,5050) TPMIN-TKELVI
          WRITE(NFECRA,5060) XTPMIN, YTPMIN, ZTPMIN
          WRITE(NFECRA,5070) QCMIN
          WRITE(NFECRA,5080) QRMIN
C
        ENDIF
C
      ENDIF
C
C -------
C FORMATS
C -------
C
 1000 FORMAT (/, 3X,'** INFORMATIONS SUR LA TEMPERATURE DES PAROIS',/,
     &           3X,'   ------------------------------------------')
C
 1010 FORMAT('------------------------------------'
     &         ,'-----------------------------------')
C
 2010 FORMAT('ATTENTION, temperature de paroi relaxee a ',G7.2,'%'
     &      ,' en (',I8,' points)')
C
 2020 FORMAT('ATTENTION, temperature de paroi CLIPPE AU MIN-MAX :')
C
 2030 FORMAT('NOMBRE DE POINTS CLIPPE AU MINIMUM : ',I8)
C
 2040 FORMAT('NOMBRE DE POINTS CLIPPE AU MAXIMUM : ',I8)
C
 2050 FORMAT('Variation maximale : ',G9.4,'%')
C
 2060 FORMAT('Temperature de paroi diminuant  : ',I8,' faces de bord')
C
 2070 FORMAT('Temperature de paroi augmentant : ',I8,' faces de bord')
C
 3000 FORMAT(I10,8X,E10.4,4X,E10.4,4X,E10.4,4X,E10.4)
C
 3010 FORMAT('Grises ou noires Temp max (C)  '
     &      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')
C
 3020 FORMAT('Profils imposes  Temp max (C)  '
     &      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')
C
 3030 FORMAT('Parois a EPS=0   Temp max (C)  '
     &      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')
C
 3040 FORMAT('Flux imp EPS!=0  Temp max (C)  '
     &      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')
C
 3050 FORMAT('Flux imp EPS=0   Temp max (C)  '
     &      ,'Temp min (C)  Temp moy (C)  Flux Net (W)')
C
 5010 FORMAT(/,12X,'TEMPERATURE PAROI MAXIMALE (Degre celsius) = ',
     &         G15.7)
 5020 FORMAT(15X,' AU POINT X Y Z =  ',E10.4,2X,E10.4,2X,E10.4)
C
 5030 FORMAT(15X,' FLUX CONVECTIF  = ',G15.7)
C
 5040 FORMAT(15X,' FLUX RADIATIF  = ',G15.7)
C
 5050 FORMAT(/,12X,'TEMPERATURE PAROI MINIMALE (Degre celsius) = ',
     &         G15.7)
C
 5060 FORMAT(15X,' AU POINT X Y Z =  ',E10.4,2X,E10.4,2X,E10.4)
C
 5070 FORMAT(15X,' FLUX CONVECTIF  = ',G15.7)
C
 5080 FORMAT(15X,' FLUX RADIATIF  = ',G15.7,/)
C
C
C ---
C FIN
C ---
C
      RETURN
C
      END
c@z
C
C   ---------------------------------------------------------------------------------
C   Variation maximale : .1484E-13%
C   Temperature de paroi diminuant  :       33
C   Temperature de paroi augmentant :       27
C   ---------------------------------------------------------------------------------
C   Zones de parois  Temp max (C)    Temp min (C)    Temp moy (C)    Flux Net (W)
C           15        0.7941E+03      0.7410E+03      0.7688E+03      0.7688E+03
C   ---------------------------------------------------------------------------------
C   Profils imposes  Temp max (C)    Temp min (C)    Temp moy (C)    Flux Net (W)
C           15        0.7941E+03      0.7410E+03      0.7688E+03      0.7688E+03
C   ---------------------------------------------------------------------------------
C
