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
                        SUBROUTINE MTPROJ
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVEP   , NIVEP  , NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ITEPA  ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , ETTPA  , TEPA   , STATIS , TSLAGR , PARBOR ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC MODIFICATION UTILISATEUR EN FIN DE PAS DE TEMPS POUR MATISSE
CFONC
CFONC   COPIE ET SPECIALISATION DE USPROJ
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
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULES AUTORISE    !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NTERSL       ! E  !  -> ! NBR TERMES SOURCES DE COUPLAGE RETOUR!
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VAR STATISTIQUES LAGRANGIEN!
CARGU ! NVISBR       ! E  !  -> ! NOMBRE DE STATISTIQUES AUX FRONTIERES!
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
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (NFAC+1)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (LNDFAC)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !  (NFABOR+1)  !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (LNDFBR  ) !    !     !  (OPTIONNEL)                         !
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
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
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! ETTPA        ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE PRECEDENTE    !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! STATIS       ! TR !  -> ! MOYENNES STATISTIQUES                !
CARGU !(NCELET,NVLSTA!    !     !                                      !
CARGU ! TSLAGR       ! TR !  -> ! TERME DE COUPLAGE RETOUR DU          !
CARGU !(NCELET,NTERSL!    !     !   LAGRANGIEN SUR LA PHASE PORTEUSE   !
CARGU ! PARBOR       ! TR !  -> ! INFOS SUR INTERACTION DES PARTICULES !
CARGU !(NFABOR,NVISBR!    !     !   AUX FACES DE BORD                  !
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
      INCLUDE "dimfbr.h"
      INCLUDE "paramx.h"
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "matiss.h"
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
      INTEGER          NBPMAX , NVP    , NVEP   , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ITEPA(NBPMAX,NIVEP)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NDIMFB,*)
      DOUBLE PRECISION COEFA(NDIMFB,*), COEFB(NDIMFB,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , ETTPA(NBPMAX,NVP)
      DOUBLE PRECISION TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION STATIS(NCELET,NVLSTA) , TSLAGR(NCELET,NTERSL)
      DOUBLE PRECISION PARBOR(NFABOR,NVISBR)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IEL    , IFAC   , IFML   , ICOUL
      INTEGER          IPHAS  , IUIPH  , IFLMAB
      DOUBLE PRECISION TS0    , VS0    , TAAMAX , TPCMAX , TPPMAX
      DOUBLE PRECISION FLMASS , BILENT , POTFLO
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
C --- Gestion memoire
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C --- Une seule phase
C
      IPHAS = 1
C
C=======================================================================
C 1. AFFICHAGES
C=======================================================================
C
C
C --- Affichages au dernier pas de temps seulement
C
      IF (NTMABS .EQ. NTCABS) THEN
C
C --- Nombre de Richardson : calcul et affichage par
C       . mtimpi en convection forcee
C       . mttsns en convection naturelle
C
C
C --- Puissance totale
        IF (IRANGP.LE.0) THEN
          WRITE(IMPMAT,1001) PUITOT
        ENDIF
C
C
C --- Debit enthalpique
        IF (IRANGP.LE.0) THEN
          WRITE(IMPMAT,1002) DEBCON
        ENDIF
C
C
C --- Coeff d'echange
        IF (IRANGP.LE.0) THEN
          WRITE(IMPMAT,1011) CFECCA
          WRITE(IMPMAT,1012) CFECMA
        ENDIF
C
C --- Bilan masse en sortie
C       (la correction par FRDTRA est la correction correspondant
C        au rapport d'echelle transverse entre le reel et le modele)
        IUIPH  = IU(IPHAS)
        IFLMAB = IPPROB(IFLUMA(IUIPH))
        FLMASS = 0.D0
        DO IFAC = 1, NFABOR
          IFML  = IFMFBR(IFAC)
          ICOUL = IPRFML(IFML,1)
          IF (ICOUL.EQ.ICMTFO) THEN
            FLMASS = FLMASS + PROPFB(IFAC,IFLMAB)
          ENDIF
        ENDDO
        FLMASS = FLMASS * FRDTRA
        IF (IRANGP.GE.0) CALL PARSOM(FLMASS)
C
        IF (IRANGP.LE.0) THEN
          WRITE(IMPMAT,1021) FLMASS
        ENDIF
C
C --- Temperature moyenne dans la cheminee d'evacuation
C     (scalaire 1)
        TS0 = 0.D0
        VS0 = 0.D0
        DO IEL = 1, NCEL
          IFML  = IFMCEL(IEL   )
          ICOUL = IPRFML(IFML,1)
          IF(ICOUL.EQ.ICMTCO) THEN
            TS0 = TS0 + VOLUME(IEL)*RTP(IEL,ISCA(ITAAMT))
            VS0 = VS0 + VOLUME(IEL)
          ENDIF
        ENDDO
        TS0 = TS0/MAX(VS0,EPZERO)
        IF (IRANGP.LE.0) THEN
          WRITE(IMPMAT,1022) TS0
        ENDIF
C
C
C --- Calcul du bilan enthalpique en pourcentage
C       On calcule le rapport de rhoUS * Cp * Delta T (Watt)
C         a PUITOT (multiplication par 1.D3 car PUITOT est en kW,
C         division par 100 pour obtenir une donnee en %)
C       L'ecart de temperature est pris comme l'ecart entre la
C         temperature de sortie TS0 et la temperature d'entree TINIT
        BILENT =
     &       CP0(IPHAS)*FLMASS*(TS0-TINIT)/(PUITOT*1.D3)*100.D0
        IF (IRANGP.LE.0) THEN
          WRITE(IMPMAT,1031) BILENT
        ENDIF
C
C
C --- Calcul du potentiel de flottabilite
C       Calcul de delta_rho * g * delta_h en Pascal avec, a pression
C         constante : delta_rho = rho_ref/T_ref * delta_T
C         D'ou POTFLO = rho_ref/T_ref * delta_T * g * Delta_h
C       L'ecart de temperature est pris comme l'écart entre la
C         temperature de sortie TS0 et la temperature d'entree TINIT
C       L'ecart de hauteur est pris entre le haut de la cheminee de
C         sortie et la mi hauteur de la zone de stockage (z=0 au sol)
        POTFLO = RRFMAT/(TRFMAT+TKELVI)*(TS0-TINIT)
     &       * SQRT(GX**2+GY**2+GZ**2)*(HCHEVA-0.5D0*EPCHEL*NCHEST)
        IF (IRANGP.LE.0) THEN
          WRITE(IMPMAT,1032) POTFLO
        ENDIF
C
C
C --- Calcul des max des scalaires
C     . TAA* : Temperature Air Ambiant (scalaire ITAAMT)
C     . TPC* : Temperature Peau Colis  (scalaire ITPCMT)
C     . TPP* : Temperature Peau Paroi  (scalaire ITPPMT)
C
        TAAMAX = 0.D0
        TPCMAX = 0.D0
        TPPMAX = 0.D0
C
        DO IEL = 1, NCEL
          TAAMAX = MAX(TAAMAX , RTP(IEL,ISCA(ITAAMT)))
          TPCMAX = MAX(TPCMAX , RTP(IEL,ISCA(ITPCMT)))
          TPPMAX = MAX(TPPMAX , RTP(IEL,ISCA(ITPPMT)))
        ENDDO
C
        IF (IRANGP.GE.0) THEN
          CALL PARMAX(TAAMAX)
          CALL PARMAX(TPCMAX)
          CALL PARMAX(TPPMAX)
        ENDIF
C
        IF (IRANGP.LE.0) THEN
          WRITE(IMPMAT,1041) TAAMAX
          WRITE(IMPMAT,1042) TPCMAX
          WRITE(IMPMAT,1043) TPPMAX
        ENDIF
C
        CLOSE(IMPMAT)
C
C       Fin du test sur NTCABS
      ENDIF
C
C
C--------
C FORMATS
C--------
C
C
 1001 FORMAT(' Puissance totale de l''installation                   ',
     &'  :',E12.5, ' kW')
 1002 FORMAT(' Debit enthalpique vers le ciel de l''entrepot         ',
     &'  :',E12.5, ' kW')
C
 1011 FORMAT(' Coefficient d''echange moyen conteneur/air            ',
     &'  :',E12.5,' W/m2/C')
 1012 FORMAT(' Coefficient d''echange moyen mur/air                  ',
     &'  :',E12.5,' W/m2/C')
C
 1021 FORMAT(' Debit masse de circulation d''air                     ',
     &'  :',E12.5, ' kg/s')
 1022 FORMAT(' Temperature d''air en sortie                          ',
     &'  :',E12.5, ' °C')
C
 1031 FORMAT(' Bilan enthalpique                                    ',
     &'  :',E12.5, ' %')
 1032 FORMAT(' Potentiel de flottabilite                            ',
     &'  :',E12.5, ' Pa')
C
 1041 FORMAT(' Temperature maximale d''air ambiant                   ',
     &'  :',E12.5, ' °C')
 1042 FORMAT(' Temperature maximale des conteneurs                  ',
     &'  :',E12.5, ' °C')
 1043 FORMAT(' Temperature maximale des murs                        ',
     &'  :',E12.5, ' °C')
C
      RETURN
      END
c@z
