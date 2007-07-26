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
                        SUBROUTINE SCHTMP
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IAPPEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ISOSTD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  TRAITEMENT DU FLUX DE MASSE, DE LA VISCOSITE, DE LA MASSE
CFONC  VOLUMIQUE, DE LA CHALEUR SPECIFIQUE ET DU TABLEAU TSNSA
CFONC  DANS LE CAS D'UN THETA SCHEMA
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
CARGU ! IAPPEL       ! E  !  -> ! NUMERO DE L'APPEL (AVANT OU APRES    !
CARGU !              !    !     ! PHYVAR                               !
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
CARGU !   (NFABOR+1) !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (LNDFBR)   !    !     !  (OPTIONNEL)                         !
CARGU ! ISOSTD       ! TE !  -> ! INDICATEUR DE SORTIE STANDARD        !
CARGU !    (NFABOR+1)!    !     !  +NUMERO DE LA FACE DE REFERENCE     !
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
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "cstphy.h"
      INCLUDE "pointe.h"
      INCLUDE "period.h"
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
      INTEGER          NVAR   , NSCAL  , NPHAS  , IAPPEL
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ISOSTD(NFABOR+1,NPHAS)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IEL    , IFAC   , ISCAL
      INTEGER          IPHAS  , IUIPH
      INTEGER          IPCROM , IPCROA
      INTEGER          IPBROM , IPBROA
      INTEGER          IFLMAS , IFLMAB , IFLMBA, IFLMSA
      INTEGER          IPCVIS , IPCVST
      INTEGER          IPCVSA , IPCVTA , IPCVSL
      INTEGER          IICP   , IICPA
      DOUBLE PRECISION FLUX   , THETA  , AA, BB, VISCOS, XMASVO, VARCP
C
C***********************************************************************
C
C
C=======================================================================
C 0. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1.  AU TOUT DEBUT DE LA BOUCLE EN TEMPS
C=======================================================================
C
      IF(IAPPEL.EQ.1) THEN
C
C --- Application du schema en temps sur le flux de masse
C     Soit F le flux de masse
C       - si ISTMPF = 2 (schema non standard, theta>0 en fait = 0.5)
C         PROPFA(1,IFLMAS) contient F_(n-2+theta) et on y met F(n-1+theta)
C         PROPFA(1,IFLMSA) contient F_(n-1+theta) et
C                                    on y met une extrapolation en n+theta
C       - si ISTMPF = 1 (schema non standard, theta=0)
C         PROPFA(1,IFLMAS) contient deja F_n et
C         PROPFA(1,IFLMSA) n'est pas utilise : on ne fait rien
C       - sinon : ISTMPF = 0 (schema standard, theta= -999)
C         PROPFA(1,IFLMAS) et PROPFA(1,IFLMSA) contiennent tous les deux F(n)
C                                            : on ne fait rien
C
C
C     Ordre 2 en temps pour le flux (theta = 0.5) a entrer dans navsto
C       Flux convectif = 2F(n-1+theta)-F(n-2+theta)
C       Flux conserve  =  F(n-1+theta)
C     Au premier pas de temps, l'ancien a ete initialise dans inivar (a 0)
C       en suite de calcul, les deux ont ete relus.
C
        DO IPHAS = 1, NPHAS
          IF(ISTMPF(IPHAS).EQ.2) THEN
            IUIPH  = IU(IPHAS)
            IFLMAS = IPPROF(IFLUMA(IUIPH))
            IFLMAB = IPPROB(IFLUMA(IUIPH))
            IFLMSA = IPPROF(IFLUAA(IUIPH))
            IFLMBA = IPPROB(IFLUAA(IUIPH))
            DO IFAC = 1 , NFAC
              FLUX                =      PROPFA(IFAC,IFLMAS)
              PROPFA(IFAC,IFLMAS) = 2.D0*PROPFA(IFAC,IFLMAS)
     &                                 - PROPFA(IFAC,IFLMSA)
              PROPFA(IFAC,IFLMSA) = FLUX
            ENDDO
            DO IFAC = 1 , NFABOR
              FLUX                =      PROPFB(IFAC,IFLMAB)
              PROPFB(IFAC,IFLMAB) = 2.D0*PROPFB(IFAC,IFLMAB)
     &                                 - PROPFB(IFAC,IFLMBA)
              PROPFB(IFAC,IFLMBA) = FLUX
            ENDDO
          ENDIF
        ENDDO
C
C     Les valeurs courantes ecrasent les valeurs anterieures
C       en cas d'extrapolation en temps (theta > 0)
C       Pour RHO, on le fait en double si ICALHY = 1 (et sur NCELET)
C     Au debut du calcul les flux nouveau et ancien ont ete initialises inivar
        DO IPHAS = 1, NPHAS
          IF(IROEXT(IPHAS).GT.0) THEN
            IPCROM = IPPROC(IROM  (IPHAS))
            IPCROA = IPPROC(IROMA (IPHAS))
            DO IEL = 1, NCELET
              PROPCE(IEL,IPCROA) = PROPCE(IEL,IPCROM)
            ENDDO
            IPBROM = IPPROB(IROM  (IPHAS))
            IPBROA = IPPROB(IROMA (IPHAS))
            DO IFAC = 1, NFABOR
              PROPFB(IFAC,IPBROA) = PROPFB(IFAC,IPBROM)
            ENDDO
          ENDIF
        ENDDO
        DO IPHAS = 1, NPHAS
          IF(IVIEXT(IPHAS).GT.0) THEN
            IPCVIS = IPPROC(IVISCL(IPHAS))
            IPCVST = IPPROC(IVISCT(IPHAS))
            IPCVSA = IPPROC(IVISLA(IPHAS))
            IPCVTA = IPPROC(IVISTA(IPHAS))
            DO IEL = 1, NCEL
              PROPCE(IEL,IPCVSA) = PROPCE(IEL,IPCVIS)
              PROPCE(IEL,IPCVTA) = PROPCE(IEL,IPCVST)
            ENDDO
          ENDIF
        ENDDO
        DO IPHAS = 1, NPHAS
          IF(ICPEXT(IPHAS).GT.0) THEN
            IF(ICP   (IPHAS).GT.0) THEN
              IICP   = IPPROC(ICP   (IPHAS))
              IICPA  = IPPROC(ICPA  (IPHAS))
              DO IEL = 1, NCEL
                PROPCE(IEL,IICPA ) = PROPCE(IEL,IICP  )
              ENDDO
            ENDIF
          ENDIF
        ENDDO
C
C     Remarque : si on faisant cette operation pour tous les
C       scalaires, on la ferait plusieurs fois pour les scalaires
C       ayant une variance
        IF (NSCAL.GE.1) THEN
          DO ISCAL = 1, NSCAL
            IF(IVISLS(ISCAL).GT.0.AND.ISCAVR(ISCAL).LE.0) THEN
              IF(IVSEXT(ISCAL).GT.0) THEN
                IPCVSL = IPPROC(IVISLS(ISCAL))
                IPCVSA = IPPROC(IVISSA(ISCAL))
                DO IEL = 1, NCEL
                  PROPCE(IEL,IPCVSA) = PROPCE(IEL,IPCVSL)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF
C
        RETURN
C
C
C=======================================================================
C 2.  JUSTE APRES PHYVAR (ET DONC AVANT NAVSTO)
C=======================================================================
C
C
      ELSEIF(IAPPEL.EQ.2) THEN
C
C 2.1 MISE A JOUR DES VALEURS ANCIENNES
C =====================================
C
C     On passe ici dans le cas ou l'on suspecte que la valeur portee
C       par les variables "anciennes" n'est pas satisfaisante pour
C       extrapoler.
C     On passe ici au premier pas de temps et lorsque le fichier suite
C       ne comportait pas grandeur requise.
C
        DO IPHAS = 1, NPHAS
          IF(INITRO(IPHAS).NE.1) THEN
            INITRO(IPHAS) = 1
            IF(IROEXT(IPHAS).GT.0) THEN
              IPCROM = IPPROC(IROM  (IPHAS))
              IPCROA = IPPROC(IROMA (IPHAS))
              DO IEL = 1, NCELET
                PROPCE(IEL,IPCROA) = PROPCE(IEL,IPCROM)
              ENDDO
              IPBROM = IPPROB(IROM  (IPHAS))
              IPBROA = IPPROB(IROMA (IPHAS))
              DO IFAC = 1, NFABOR
                PROPFB(IFAC,IPBROA) = PROPFB(IFAC,IPBROM)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
        DO IPHAS = 1, NPHAS
          IF(INITVI(IPHAS).NE.1) THEN
            INITVI(IPHAS) = 1
            IF(IVIEXT(IPHAS).GT.0) THEN
              IPCVIS = IPPROC(IVISCL(IPHAS))
              IPCVST = IPPROC(IVISCT(IPHAS))
              IPCVSA = IPPROC(IVISLA(IPHAS))
              IPCVTA = IPPROC(IVISTA(IPHAS))
              DO IEL = 1, NCEL
                PROPCE(IEL,IPCVSA) = PROPCE(IEL,IPCVIS)
                PROPCE(IEL,IPCVTA) = PROPCE(IEL,IPCVST)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
        DO IPHAS = 1, NPHAS
          IF(INITCP(IPHAS).NE.1) THEN
            INITCP(IPHAS) = 1
            IF(ICPEXT(IPHAS).GT.0) THEN
              IF(ICP   (IPHAS).GT.0) THEN
                IICP   = IPPROC(ICP   (IPHAS))
                IICPA  = IPPROC(ICPA  (IPHAS))
                DO IEL = 1, NCEL
                  PROPCE(IEL,IICPA ) = PROPCE(IEL,IICP  )
                ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDDO
C
C     Remarque : si on faisant cette operation pour tous les
C       scalaires, on la ferait plusieurs fois pour les scalaires
C       ayant une variance
        IF (NSCAL.GE.1) THEN
          DO ISCAL = 1, NSCAL
            IF(INITVS(ISCAL).NE.1) THEN
              INITVS(ISCAL) = 1
              IF(IVISLS(ISCAL).GT.0.AND.ISCAVR(ISCAL).LE.0) THEN
                IF(IVSEXT(ISCAL).GT.0) THEN
                  IPCVSL = IPPROC(IVISLS(ISCAL))
                  IPCVSA = IPPROC(IVISSA(ISCAL))
                  DO IEL = 1, NCEL
                    PROPCE(IEL,IPCVSA) = PROPCE(IEL,IPCVSL)
                  ENDDO
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF
C
C
C 2.2 EXTRAPOLATION DES NOUVELLES VALEURS
C =======================================
C
C --- Extrapolation de la viscosite et de la masse volumique dans le cas d'un
C     theta schema
C     A partir de Fn-1 et Fn on calcule Fn+theta
C     On conserve les nouvelles valeurs dans l'ancien tableau pour
C     retablir en fin de pas de temps
C
C     Le calcul pour Rho est fait sur NCELET afin d'economiser un echange.
C
        DO IPHAS = 1, NPHAS
          IF(IROEXT(IPHAS).GT.0) THEN
            IPCROM = IPPROC(IROM  (IPHAS))
            IPCROA = IPPROC(IROMA (IPHAS))
            THETA  = THETRO(IPHAS)
            DO IEL = 1, NCELET
              XMASVO = PROPCE(IEL,IPCROM)
              PROPCE(IEL,IPCROM) = (1.D0+THETA) * PROPCE(IEL,IPCROM)
     &                           -       THETA  * PROPCE(IEL,IPCROA)
              PROPCE(IEL,IPCROA) = XMASVO
            ENDDO
            IPBROM = IPPROB(IROM  (IPHAS))
            IPBROA = IPPROB(IROMA (IPHAS))
            DO IFAC = 1, NFABOR
              XMASVO = PROPFB(IFAC,IPBROM)
              PROPFB(IFAC,IPBROM) = (1.D0+THETA) * PROPFB(IFAC,IPBROM)
     &                            -       THETA  * PROPFB(IFAC,IPBROA)
              PROPFB(IFAC,IPBROA) = XMASVO
            ENDDO
          ENDIF
          IF(IVIEXT(IPHAS).GT.0) THEN
            IPCVIS = IPPROC(IVISCL(IPHAS))
            IPCVST = IPPROC(IVISCT(IPHAS))
            IPCVSA = IPPROC(IVISLA(IPHAS))
            IPCVTA = IPPROC(IVISTA(IPHAS))
            THETA  = THETVI(IPHAS)
            DO IEL = 1, NCEL
              VISCOS = PROPCE(IEL,IPCVIS)
              PROPCE(IEL,IPCVIS) = (1.D0+THETA) * PROPCE(IEL,IPCVIS)
     &                           -       THETA  * PROPCE(IEL,IPCVSA)
              PROPCE(IEL,IPCVSA) = VISCOS
              VISCOS = PROPCE(IEL,IPCVST)
              PROPCE(IEL,IPCVST) = (1.D0+THETA) * PROPCE(IEL,IPCVST)
     &                           -       THETA  * PROPCE(IEL,IPCVTA)
              PROPCE(IEL,IPCVTA) = VISCOS
            ENDDO
          ENDIF
          IF(ICPEXT(IPHAS).GT.0) THEN
            IF(ICP(IPHAS).GT.0) THEN
              IICP   = IPPROC(ICP   (IPHAS))
              IICPA  = IPPROC(ICPA  (IPHAS))
              THETA  = THETCP(IPHAS)
              DO IEL = 1, NCEL
                VARCP  = PROPCE(IEL,IICP  )
                PROPCE(IEL,IICP  ) = (1.D0+THETA) * PROPCE(IEL,IICP  )
     &                             -       THETA  * PROPCE(IEL,IICPA )
                PROPCE(IEL,IPCVSA) = VARCP
              ENDDO
            ENDIF
          ENDIF
        ENDDO
C
C     Remarque : si on faisant cette operation pour tous les
C       scalaires, on la ferait plusieurs fois pour les scalaires
C       ayant une variance ET CE SERAIT FAUX
        IF (NSCAL.GE.1) THEN
          DO ISCAL = 1, NSCAL
            IF(IVISLS(ISCAL).GT.0.AND.ISCAVR(ISCAL).LE.0) THEN
              IF(IVSEXT(ISCAL).GT.0) THEN
                THETA  = THETVS(ISCAL)
                IPCVSL = IPPROC(IVISLS(ISCAL))
                IPCVSA = IPPROC(IVISSA(ISCAL))
                IF(IPCVSL.GT.0) THEN
                  DO IEL = 1, NCEL
                    VISCOS = PROPCE(IEL,IPCVSA)
                    PROPCE(IEL,IPCVSL) = (1.D0+THETA)*PROPCE(IEL,IPCVSL)
     &                                 -       THETA *PROPCE(IEL,IPCVSA)
                    PROPCE(IEL,IPCVSA) = VISCOS
                  ENDDO
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF
C
        RETURN
C
C
C
C=======================================================================
C 3.  JUSTE APRES NAVSTO, DANS LES BOUCLES U/P ET ALE
C=======================================================================
C
      ELSEIF(IAPPEL.EQ.3) THEN
C
C     On traite ici le flux de masse uniquement
C        On suppose qu'il n'y en a qu'un seul.
C
C     Si ISTMPF = 1 : standard : on ne fait rien
C     Si ISTMPF = 2 : ordre 2 (THETFL > 0 : = 0.5)
C       on calcule F(n+theta) par interpolation a partir
C       de F_(n-1+theta) et F(n+1) et on le met dans PROPFA(1,IFLMAS)
C     Si ISTMPF = 0 : explicite (THETFL = 0) : on remet F(n) dans
C       PROPFA(1,IFLMAS) sauf a la derniere iteration (un traitement
C       complementaire sera fait en IAPPEL=4)
C
C     Dans le cas ou on itere sur navsto, on passe ici
C     - a toutes les iterations si ISTMPF.NE.0
C     - a toutes les iterations sauf la derniere si ISTMPF.EQ.0
C
        DO IPHAS = 1, NPHAS
          IUIPH  = IU(IPHAS)
          IFLMAS = IPPROF(IFLUMA(IUIPH))
          IFLMAB = IPPROB(IFLUMA(IUIPH))
          IF(ISTMPF(IPHAS).EQ.2) THEN
            IFLMSA = IPPROF(IFLUAA(IUIPH))
            IFLMBA = IPPROB(IFLUAA(IUIPH))
            THETA  = THETFL(IPHAS)
            AA = 1.D0/(2.D0-THETA)
            BB = (1.D0-THETA)/(2.D0-THETA)
            DO IFAC = 1 , NFAC
              PROPFA(IFAC,IFLMAS) = AA * PROPFA(IFAC,IFLMAS)
     &                            + BB * PROPFA(IFAC,IFLMSA)
            ENDDO
            DO IFAC = 1 , NFABOR
              PROPFB(IFAC,IFLMAB) = AA * PROPFB(IFAC,IFLMAB)
     &                            + BB * PROPFB(IFAC,IFLMBA)
            ENDDO
          ELSEIF(ISTMPF(IPHAS).EQ.0) THEN
            IFLMSA = IPPROF(IFLUAA(IUIPH))
            IFLMBA = IPPROB(IFLUAA(IUIPH))
            DO IFAC = 1 , NFAC
              PROPFA(IFAC,IFLMAS) = PROPFA(IFAC,IFLMSA)
            ENDDO
            DO IFAC = 1 , NFABOR
              PROPFB(IFAC,IFLMAB) = PROPFB(IFAC,IFLMBA)
            ENDDO
          ENDIF
        ENDDO
C
        RETURN
C
C=======================================================================
C 3.  JUSTE APRES NAVSTO (ET STRDEP), HORS LES BOUCLES U/P ET ALE
C=======================================================================
C
      ELSEIF(IAPPEL.EQ.4) THEN
C
C     On traite ici le flux de masse uniquement
C        On suppose qu'il n'y en a qu'un seul.
C
C     Si ISTMPF = 1 : standard : on ne fait rien
C     Si ISTMPF = 2 : ordre 2 (THETFL > 0 : = 0.5)
C       on calcule F(n+theta) par interpolation a partir
C       de F_(n-1+theta) et F(n+1) et on le met dans PROPFA(1,IFLMAS)
C     Si ISTMPF = 0 : explicite (THETFL = 0)
C       On sauvegarde F_(n+1) dans PROPFA(1,IFLMSA),mais on continue
C       les calculs avec F_(n) mis dans PROPFA(1,IFLMAS)
C
C     On retablira au dernier appel de schtmp pour ISTMPF = 0
C
C     Dans le cas ou on itere sur navsto, on passe ici
C       - a toutes les iterations            si ISTMPF.NE.0
C       - uniquement a la derniere iteration si ISTMPF.EQ.0
C         (ce faisant, a partir de la deuxieme sous-iteration,
C          le calcul sera fait avec F(n+1) et plus F(n), mais on
C          suppose que l'utilisateur a choisi de faire des sous-iter
C          aussi pour impliciter le flux de masse)
C
        DO IPHAS = 1, NPHAS
          IUIPH  = IU(IPHAS)
          IFLMAS = IPPROF(IFLUMA(IUIPH))
          IFLMAB = IPPROB(IFLUMA(IUIPH))
          IF(ISTMPF(IPHAS).EQ.2) THEN
            IFLMSA = IPPROF(IFLUAA(IUIPH))
            IFLMBA = IPPROB(IFLUAA(IUIPH))
            THETA  = THETFL(IPHAS)
            AA = 1.D0/(2.D0-THETA)
            BB = (1.D0-THETA)/(2.D0-THETA)
            DO IFAC = 1 , NFAC
              PROPFA(IFAC,IFLMAS) = AA * PROPFA(IFAC,IFLMAS)
     &                            + BB * PROPFA(IFAC,IFLMSA)
            ENDDO
            DO IFAC = 1 , NFABOR
              PROPFB(IFAC,IFLMAB) = AA * PROPFB(IFAC,IFLMAB)
     &                            + BB * PROPFB(IFAC,IFLMBA)
            ENDDO
          ELSEIF(ISTMPF(IPHAS).EQ.0) THEN
            IFLMSA = IPPROF(IFLUAA(IUIPH))
            IFLMBA = IPPROB(IFLUAA(IUIPH))
            DO IFAC = 1 , NFAC
              FLUX = PROPFA(IFAC,IFLMAS)
              PROPFA(IFAC,IFLMAS) = PROPFA(IFAC,IFLMSA)
              PROPFA(IFAC,IFLMSA) = FLUX
            ENDDO
            DO IFAC = 1 , NFABOR
              FLUX = PROPFB(IFAC,IFLMAB)
              PROPFB(IFAC,IFLMAB) = PROPFB(IFAC,IFLMBA)
              PROPFB(IFAC,IFLMBA) = FLUX
            ENDDO
          ENDIF
        ENDDO
C
        RETURN
C
C=======================================================================
C 3.  JUSTE APRES SCALAI
C=======================================================================
C
      ELSEIF(IAPPEL.EQ.5) THEN
C
C 3.1 RETABLISSEMENT POUR LE FLUX DE MASSE
C ========================================
C
C     On corrige les manipulations sur le flux de masse faites dans
C       l'appel precedent afin d'etre pret pour le pas de temps suivant.
C
C     Si ISTMPF = 1 : standard : on ne fait rien
C     Si ISTMPF = 2 : ordre 2 (THETFL > 0 : = 0.5) : on ne fait rien
C     Si ISTMPF = 0 : explicite (THETFL = 0)
C       on remet F_(n+1) (stocke dans PROPFA(1,IFLMSA)) dans PROPFA(1,IFLMAS)
C       de sorte que les deux flux contiennent la meme chose
C
        DO IPHAS = 1, NPHAS
          IF(ISTMPF(IPHAS).EQ.0) THEN
            IUIPH  = IU(IPHAS)
            IFLMAS = IPPROF(IFLUMA(IUIPH))
            IFLMAB = IPPROB(IFLUMA(IUIPH))
            IFLMSA = IPPROF(IFLUAA(IUIPH))
            IFLMBA = IPPROB(IFLUAA(IUIPH))
            DO IFAC = 1 , NFAC
              PROPFA(IFAC,IFLMAS) = PROPFA(IFAC,IFLMSA)
            ENDDO
            DO IFAC = 1 , NFABOR
              PROPFB(IFAC,IFLMAB) = PROPFB(IFAC,IFLMBA)
            ENDDO
          ENDIF
        ENDDO
C
C 3.1 RETABLISSEMENT POUR LES PROPRIETES PHYSIQUES
C ================================================
C
C     Le calcul pour Rho est fait sur NCELET afin d'economiser un echange.
C
        DO IPHAS = 1, NPHAS
          IF(IROEXT(IPHAS).GT.0) THEN
            IPCROM = IPPROC(IROM  (IPHAS))
            IPCROA = IPPROC(IROMA (IPHAS))
            DO IEL = 1, NCELET
              PROPCE(IEL,IPCROM) = PROPCE(IEL,IPCROA)
            ENDDO
            IPBROM = IPPROB(IROM  (IPHAS))
            IPBROA = IPPROB(IROMA (IPHAS))
            DO IFAC = 1, NFABOR
              PROPFB(IFAC,IPBROM) = PROPFB(IFAC,IPBROA)
            ENDDO
          ENDIF
          IF(IVIEXT(IPHAS).GT.0) THEN
            IPCVIS = IPPROC(IVISCL(IPHAS))
            IPCVST = IPPROC(IVISCT(IPHAS))
            IPCVSA = IPPROC(IVISLA(IPHAS))
            IPCVTA = IPPROC(IVISTA(IPHAS))
            DO IEL = 1, NCEL
              PROPCE(IEL,IPCVIS) = PROPCE(IEL,IPCVSA)
              PROPCE(IEL,IPCVST) = PROPCE(IEL,IPCVTA)
            ENDDO
          ENDIF
          IF(ICPEXT(IPHAS).GT.0) THEN
            IF(ICP(IPHAS).GT.0) THEN
              IICP   = IPPROC(ICP   (IPHAS))
              IICPA  = IPPROC(ICPA  (IPHAS))
              DO IEL = 1, NCEL
                PROPCE(IEL,IICP  ) = PROPCE(IEL,IICPA )
              ENDDO
            ENDIF
          ENDIF
        ENDDO
C
C     Remarque : si on faisant cette operation pour tous les
C       scalaires, on la ferait plusieurs fois pour les scalaires
C       ayant une variance
        IF (NSCAL.GE.1) THEN
          DO ISCAL = 1, NSCAL
            IF(IVISLS(ISCAL).GT.0.AND.ISCAVR(ISCAL).LE.0) THEN
              IF(IVSEXT(ISCAL).GT.0) THEN
                IPCVSL = IPPROC(IVISLS(ISCAL))
                IPCVSA = IPPROC(IVISSA(ISCAL))
                DO IEL = 1, NCEL
                  PROPCE(IEL,IPCVSL) = PROPCE(IEL,IPCVSA)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF
C
        RETURN
C
      ENDIF
C
C=======================================================================
C
C--------
C FORMATS
C--------
C
C----
C FIN
C----
C
      END
c@z
