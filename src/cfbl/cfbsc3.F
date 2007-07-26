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
                        SUBROUTINE CFBSC3
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , ICONVP , IDIFFP , NSWRGP , IMLIGP , IRCFLP ,
     &   ISCHCP , ISSTPP , INC    , IMRGRA , ICCOCG ,
     &   IPP    , IWARNP ,
     &   BLENCP , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   PVAR   , COEFAP , COEFBP , COFAFP , COFBFP ,
     &   FLUMAS , FLUMAB , VISCF  , VISCB  ,
     &   FLVARF , FLVARB ,
     &   DPDX   , DPDY   , DPDZ   , DPDXA  , DPDYA  , DPDZA  ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ---------
c@foncb
CFONC
CFONC CALCUL DU FLUX DE CONVECTION-DIFFUSION D'UNE VARIABLE AUX FACES
CFONC
CFONC                     .                    ----->        -->
CFONC  FLVARF (FACEij) =  m   PVAR  - Visc   ( grad PVAR )  . n
CFONC                      ij     ij      ij             ij    ij
CFONC
CFONC                   .                 ----->       -->
CFONC  FLVARB (FABi) =  m  PVAR - Visc  ( grad PVAR ) . n
CFONC                    i     i      i             i    i
CFONC
CFONC
CFONC CALCUL EN UPWIND
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                 ARGUMENTS
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
CARGU ! IVAR         ! E  !  -> ! NUMERO DE LA VARIABLE                !
CARGU ! ICONVP       ! E  !  -> ! INDICATEUR = 1 CONVECTION, 0 SINON   !
CARGU ! IDIFFP       ! E  !  -> ! INDICATEUR = 1 DIFFUSION , 0 SINON   !
CARGU ! NSWRGP       ! E  !  -> ! NOMBRE DE SWEEP POUR RECONSTRUCTION  !
CARGU !              !    !     !             DES GRADIENTS            !
CARGU ! IMLIGP       ! E  !  -> ! METHODE DE LIMITATION DU GRADIENT    !
CARGU !              !    !     !  < 0 PAS DE LIMITATION               !
CARGU !              !    !     !  = 0 A PARTIR DES GRADIENTS VOISINS  !
CARGU !              !    !     !  = 1 A PARTIR DU GRADIENT MOYEN      !
CARGU ! IRCFLP       ! E  !  -> ! INDICATEUR = 1 REC FLUX, 0 SINON     !
CARGU ! ISCHCP       ! E  !  -> ! INDICATEUR = 1 CENTRE , 0 2ND ORDER  !
CARGU ! ISSTPP       ! E  !  -> ! INDICATEUR = 1 SANS TEST DE PENTE    !
CARGU !              !    !     !            = 0 AVEC TEST DE PENTE    !
CARGU ! INC          ! E  !  -> ! INDICATEUR = 0 RESOL SUR INCREMENT   !
CARGU !              !    !     !              1 SINON                 !
CARGU ! IMRGRA       ! E  !  -> ! INDICATEUR = 0 GRADRC 97             !
CARGU !              ! E  !  -> !            = 1 GRADMC 99             !
CARGU ! ICCOCG       ! E  !  -> ! INDICATEUR = 1 POUR RECALCUL DE COCG !
CARGU !              !    !     !              0 SINON                 !
CARGU ! IPP          ! E  !  -> ! NUMERO DE VARIABLE POUR POST         !
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! BLENCP       ! R  !  -> ! 1 - PROPORTION D'UPWIND              !
CARGU ! EPSRGP       ! R  !  -> ! PRECISION RELATIVE POUR LA           !
CARGU !              !    !     !  RECONSTRUCTION DES GRADIENTS 97     !
CARGU ! CLIMGP       ! R  !  -> ! COEF GRADIENT*DISTANCE/ECART         !
CARGU ! EXTRAP       ! R  !  -> ! COEF EXTRAP GRADIENT                 !
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
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TE !  -  ! MACRO TABLEAU ENTIER                 !
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
CARGU ! PVAR (NCELET ! TR !  -> ! VARIABLE RESOLUE (INSTANT PRECEDENT) !
CARGU ! COEFAP, B    ! TR !  -> ! TABLEAUX DES COND LIM POUR P         !
CARGU !   (NFABOR)   !    !     !  SUR LA NORMALE A LA FACE DE BORD    !
CARGU ! COFAFP, B    ! TR !  -> ! TABLEAUX DES COND LIM POUR LE FLUX DE!
CARGU !   (NFABOR)   !    !     !  DIFFUSION DE P                      !
CARGU ! FLUMAS(NFAC) ! TR !  -> ! FLUX DE MASSE AUX FACES INTERNES     !
CARGU ! FLUMAB(NFABOR! TR !  -> ! FLUX DE MASSE AUX FACES DE BORD      !
CARGU ! VISCF (NFAC) ! TR !  -> ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU !              !    !     !  POUR SECOND MEMBRE                  !
CARGU ! VISCB (NFABOR! TR !  -> ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU !              !    !     !  POUR SECOND MEMBRE                  !
CARGU ! FLVARF(NFAC) ! TR ! <-  ! FLUX DE CONVECTION-DIFFUSION         !
CARGU !              !    !     !  AUX FACES INTERNES                  !
CARGU ! FLVARB(NFABOR! TR ! <-  ! FLUX DE CONVECTION-DIFFUSION         !
CARGU !              !    !     !  AUX FACES DE BORD                   !
CARGU ! DPDX,Y,Z     ! TR !  -  ! TABLEAU DE TRAVAIL POUR LE GRAD DE P !
CARGU !    (NCELET)  !    !     !                                      !
CARGU ! DPDXA,YA,ZA  ! TR !  -  ! TABLEAU DE TRAVAIL POUR LE GRAD DE P !
CARGU !    (NCELET)  !    !     !  AVEC DECENTREMENT AMONT             !
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
      INCLUDE "pointe.h"
      INCLUDE "vector.h"
      INCLUDE "entsor.h"
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
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IVAR   , ICONVP , IDIFFP , NSWRGP , IMLIGP
      INTEGER          IRCFLP , ISCHCP , ISSTPP
      INTEGER          INC    , IMRGRA , ICCOCG
      INTEGER          IWARNP , IPP
      DOUBLE PRECISION BLENCP , EPSRGP , CLIMGP, EXTRAP
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION PVAR (NCELET), COEFAP(NFABOR), COEFBP(NFABOR)
      DOUBLE PRECISION                COFAFP(NFABOR), COFBFP(NFABOR)
      DOUBLE PRECISION FLUMAS(NFAC), FLUMAB(NFABOR)
      DOUBLE PRECISION VISCF (NFAC), VISCB (NFABOR)
      DOUBLE PRECISION FLVARF(NFAC), FLVARB(NFABOR)
      DOUBLE PRECISION DPDX (NCELET),DPDY (NCELET),DPDZ (NCELET)
      DOUBLE PRECISION DPDXA(NCELET),DPDYA(NCELET),DPDZA(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      CHARACTER*80     CHAINE
      CHARACTER*8      CNOM
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC,II,JJ,INFAC,IEL, IIJ, III
      INTEGER          IPHYDP
      DOUBLE PRECISION PFAC,PFACD,PIP,PJP,FLUI,FLUJ,FLUX
      DOUBLE PRECISION PIF,PJF
      DOUBLE PRECISION DPXF,DPYF,DPZF
      DOUBLE PRECISION DIJPFX, DIJPFY, DIJPFZ
      DOUBLE PRECISION DIIPFX, DIIPFY, DIIPFZ
      DOUBLE PRECISION DJJPFX, DJJPFY, DJJPFZ
      DOUBLE PRECISION DIIPBX, DIIPBY, DIIPBZ
      DOUBLE PRECISION POND
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATION
c=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      CHAINE = NOMVAR(IPP)
      CNOM   = CHAINE(1:8)
C
C
C=======================================================================
C 2.  CALCUL DU BILAN AVEC TECHNIQUE DE RECONSTRUCTION
C=======================================================================
C
C ======================================================================
C ---> CALCUL DU GRADIENT DE PVAR
C ======================================================================
C    DPDX sert pour la reconstruction des flux de diffusion
C       (convection en upwind)
C    On doit donc le calculer uniquement si on a de la diffusion
C       et qu'on reconstruit les flux
C
      IF( IDIFFP.NE.0 .AND. IRCFLP.EQ.1 ) THEN
C
        IPHYDP = 0
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP ,  IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DPDXA  , DPDXA  , DPDXA  ,
     &   PVAR   , COEFAP , COEFBP ,
     &   DPDX   , DPDY   , DPDZ   ,
C        ------   ------   ------
     &   DPDXA  , DPDYA  , DPDZA  ,
     &   RDEVEL , RTUSER , RA     )
C
      ELSE
        DO IEL = 1, NCELET
          DPDX(IEL) = 0.D0
          DPDY(IEL) = 0.D0
          DPDZ(IEL) = 0.D0
        ENDDO
      ENDIF
C
C
C ======================================================================
C ---> ASSEMBLAGE A PARTIR DES FACETTES FLUIDES
C ======================================================================
C
      INFAC = 0
C
      DO IFAC = 1, NFAC
        FLVARF(IFAC) = 0.D0
      ENDDO
C
      DO IFAC = 1, NFABOR
        FLVARB(IFAC) = 0.D0
      ENDDO
C
C
C  --> FLUX UPWIND PUR
C  =====================
C
      IF (IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
        DO IFAC = 1, NFAC
C
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
C
          IIJ = IDIJPF-1+3*(IFAC-1)
          DIJPFX = RA(IIJ+1)
          DIJPFY = RA(IIJ+2)
          DIJPFZ = RA(IIJ+3)
C
          POND   = RA(IPOND-1+IFAC)
C
C ON RECALCULE A CE NIVEAU II' ET JJ'
C
          DIIPFX = CDGFAC(1,IFAC) - (XYZCEN(1,II)+
     &               (1.D0-POND) * DIJPFX)
          DIIPFY = CDGFAC(2,IFAC) - (XYZCEN(2,II)+
     &               (1.D0-POND) * DIJPFY)
          DIIPFZ = CDGFAC(3,IFAC) - (XYZCEN(3,II)+
     &               (1.D0-POND) * DIJPFZ)
          DJJPFX = CDGFAC(1,IFAC) -  XYZCEN(1,JJ)+
     &                   POND  * DIJPFX
          DJJPFY = CDGFAC(2,IFAC) -  XYZCEN(2,JJ)+
     &                   POND  * DIJPFY
          DJJPFZ = CDGFAC(3,IFAC) -  XYZCEN(3,JJ)+
     &                   POND  * DIJPFZ
C
          DPXF = 0.5D0*(DPDX(II) + DPDX(JJ))
          DPYF = 0.5D0*(DPDY(II) + DPDY(JJ))
          DPZF = 0.5D0*(DPDZ(II) + DPDZ(JJ))
C
C     reconstruction uniquement si IRCFLP = 1
          PIP = PVAR(II)
     &           + IRCFLP*(DPXF*DIIPFX+DPYF*DIIPFY+DPZF*DIIPFZ)
          PJP = PVAR(JJ)
     &           + IRCFLP*(DPXF*DJJPFX+DPYF*DJJPFY+DPZF*DJJPFZ)
C
          FLUI = 0.5D0*( FLUMAS(IFAC) +ABS(FLUMAS(IFAC)) )
          FLUJ = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
C
          PIF = PVAR(II)
          PJF = PVAR(JJ)
          INFAC = INFAC+1
C
          FLUX = ICONVP*( FLUI*PIF +FLUJ*PJF )
     &           + IDIFFP*VISCF(IFAC)*( PIP -PJP )
C
C --- FLVARF(IFAC) : flux de convection-diffusion de la variable
C                    a la face ij
C
          FLVARF(IFAC) = FLUX
C
        ENDDO
C
      ELSE
C
C VECTORISATION NON FORCEE
        DO IFAC = 1, NFAC
C
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
C
          IIJ = IDIJPF-1+3*(IFAC-1)
          DIJPFX = RA(IIJ+1)
          DIJPFY = RA(IIJ+2)
          DIJPFZ = RA(IIJ+3)
C
          POND   = RA(IPOND-1+IFAC)
C
C ON RECALCULE A CE NIVEAU II' ET JJ'
C
          DIIPFX = CDGFAC(1,IFAC) - (XYZCEN(1,II)+
     &               (1.D0-POND) * DIJPFX)
          DIIPFY = CDGFAC(2,IFAC) - (XYZCEN(2,II)+
     &               (1.D0-POND) * DIJPFY)
          DIIPFZ = CDGFAC(3,IFAC) - (XYZCEN(3,II)+
     &               (1.D0-POND) * DIJPFZ)
          DJJPFX = CDGFAC(1,IFAC) -  XYZCEN(1,JJ)+
     &                   POND  * DIJPFX
          DJJPFY = CDGFAC(2,IFAC) -  XYZCEN(2,JJ)+
     &                   POND  * DIJPFY
          DJJPFZ = CDGFAC(3,IFAC) -  XYZCEN(3,JJ)+
     &                   POND  * DIJPFZ
C
          DPXF = 0.5D0*(DPDX(II) + DPDX(JJ))
          DPYF = 0.5D0*(DPDY(II) + DPDY(JJ))
          DPZF = 0.5D0*(DPDZ(II) + DPDZ(JJ))
C
          PIP = PVAR(II)
     &           + IRCFLP*(DPXF*DIIPFX+DPYF*DIIPFY+DPZF*DIIPFZ)
          PJP = PVAR(JJ)
     &           + IRCFLP*(DPXF*DJJPFX+DPYF*DJJPFY+DPZF*DJJPFZ)
C
          FLUI = 0.5D0*( FLUMAS(IFAC) +ABS(FLUMAS(IFAC)) )
          FLUJ = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
C
          PIF = PVAR(II)
          PJF = PVAR(JJ)
          INFAC = INFAC+1
C
          FLUX = ICONVP*( FLUI*PIF +FLUJ*PJF )
     &           + IDIFFP*VISCF(IFAC)*( PIP -PJP )
C
C --- FLVARF(IFAC) : flux de convection-diffusion de la variable
C                    a la face ij
C
          FLVARF(IFAC) = FLUX
C
        ENDDO
C
      ENDIF
C
C
C ======================================================================
C ---> ASSEMBLAGE A PARTIR DES FACETTES DE BORD
C ======================================================================
C
      IF (IVECTB.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
        DO IFAC = 1, NFABOR
C
          II = IFABOR(IFAC)
C
          III = IDIIPB-1+3*(IFAC-1)
          DIIPBX = RA(III+1)
          DIIPBY = RA(III+2)
          DIIPBZ = RA(III+3)
C
          FLUI = 0.5D0*( FLUMAB(IFAC) +ABS(FLUMAB(IFAC)) )
          FLUJ = 0.5D0*( FLUMAB(IFAC) -ABS(FLUMAB(IFAC)) )
C
          PIP = PVAR(II)
     &       +IRCFLP*(DPDX(II)*DIIPBX+DPDY(II)*DIIPBY+DPDZ(II)*DIIPBZ)
C
          PFAC  = INC*COEFAP(IFAC) +COEFBP(IFAC)*PIP
          PFACD = INC*COFAFP(IFAC) +COFBFP(IFAC)*PIP
C
          FLUX = ICONVP*( FLUI*PVAR(II) +FLUJ*PFAC )
     &         + IDIFFP*VISCB(IFAC)*( PIP -PFACD )
C
C --- FLVARB(IFAC) : flux de convection-diffusion de la variable
C                    a la face de bord i
C
          FLVARB(IFAC) = FLUX
C
        ENDDO
C
      ELSE
C
        DO IFAC = 1, NFABOR
C
          II = IFABOR(IFAC)
C
          III = IDIIPB-1+3*(IFAC-1)
          DIIPBX = RA(III+1)
          DIIPBY = RA(III+2)
          DIIPBZ = RA(III+3)
C
          FLUI = 0.5D0*( FLUMAB(IFAC) +ABS(FLUMAB(IFAC)) )
          FLUJ = 0.5D0*( FLUMAB(IFAC) -ABS(FLUMAB(IFAC)) )
C
          PIP = PVAR(II)
     &       +IRCFLP*(DPDX(II)*DIIPBX+DPDY(II)*DIIPBY+DPDZ(II)*DIIPBZ)
C
          PFAC  = INC*COEFAP(IFAC) +COEFBP(IFAC)*PIP
          PFACD = INC*COFAFP(IFAC) +COFBFP(IFAC)*PIP
C
          FLUX = ICONVP*( FLUI*PVAR(II) +FLUJ*PFAC )
     &         + IDIFFP*VISCB(IFAC)*( PIP -PFACD )
C
C --- FLVARB(IFAC) : flux de convection-diffusion de la variable
C                    a la face de bord i
C
          FLVARB(IFAC) = FLUX
C
        ENDDO
C
      ENDIF
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
      RETURN
C
      END
c@z
