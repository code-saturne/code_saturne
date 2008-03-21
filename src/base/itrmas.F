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
                        SUBROUTINE ITRMAS
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   INIT   , INC    , IMRGRA , ICCOCG , NSWRGP , IMLIGP ,
     &   IPHYDP , IWARNP , NFECRA ,
     &   EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   FEXTX  , FEXTY  , FEXTZ  ,
     &   PVAR   , COEFAP , COEFBP , VISCF  , VISCB  ,
     &   VISELX , VISELY , VISELZ ,
     &   FLUMAS , FLUMAB ,
     &   DPDX   , DPDY   , DPDZ   , DPDXA  , DPDYA  , DPDZA  ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC INCREMENTATION DU FLUX DE MASSE A PARTIR DE GRAD(P)
CFONC  P = PRESSION, INCREMENT DE PRESSION, DOUBLE INCREMENT DE PRESSION
CFONC  grad(P) = GRADIENT FACETTE POUR L'INSTANT
CFONC
CFONC   .      .               --->     -->
CFONC   m   =  m     -  Visc   grad(P) . n
CFONC    ij     ij          ij       ij   ij
CFONC
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
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! INIT         ! E  !  -> ! > 0 : INITIALISATION DU FLUX DE MASSE!
CARGU ! INC          ! E  !  -> ! INDICATEUR = 0 RESOL SUR INCREMENT   !
CARGU !              !    !     !              1 SINON                 !
CARGU ! IMRGRA       ! E  !  -> ! INDICATEUR = 0 GRADRC 97             !
CARGU !              ! E  !  -> !            = 1 GRADMC 99             !
CARGU ! ICCOCG       ! E  !  -> ! INDICATEUR = 1 POUR RECALCUL DE COCG !
CARGU !              !    !     !              0 SINON                 !
CARGU ! NSWRGP       ! E  !  -> ! NOMBRE DE SWEEP POUR RECONSTRUCTION  !
CARGU !              !    !     !             DES GRADIENTS            !
CARGU ! IMLIGP       ! E  !  -> ! METHODE DE LIMITATION DU GRADIENT    !
CARGU !              !    !     !  < 0 PAS DE LIMITATION               !
CARGU !              !    !     !  = 0 A PARTIR DES GRADIENTS VOISINS  !
CARGU !              !    !     !  = 1 A PARTIR DU GRADIENT MOYEN      !
CARGU ! IWARNP       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! IPHYDP       ! E  !  -> ! INDICATEUR DE PRISE EN COMPTE DE LA  !
CARGU !              !    !     ! PRESSION HYDROSTATIQUE               !
CARGU ! NFECRA       ! E  !  -> ! UNITE DU FICHIER SORTIE STD          !
CARGU ! EPSRGP       ! R  !  -> ! PRECISION RELATIVE POUR LA           !
CARGU !              !    !     !  RECONSTRUCTION DES GRADIENTS 97     !
CARGU ! CLIMGP       ! R  !  -> ! COEF GRADIENT*DISTANCE/ECART         !
CARGU ! EXTRAP       ! R  !  -> ! COEF EXTRAP GRADIENT                 !
CARGU ! PVAR  (NCELET! TR !  -> ! VARIABLE (PRESSION)                  !
CARGU ! COEFAP, B    ! TR !  -> ! TABLEAUX DES COND LIM POUR PVAR      !
CARGU !   (NFABOR)   !    !     !  SUR LA NORMALE A LA FACE DE BORD    !
CARGU ! VISCF (NFAC) ! TR !  -> ! "VISCOSITE" FACE INTERNE(DT*SURF/DIST!
CARGU ! VISCB (NFABOR! TR !  -> ! "VISCOSITE" FACE DE BORD(DT*SURF/DIST!
CARGU ! VISELX(NCELET! TR !  -> ! "VISCOSITE" PAR CELLULE  DIR X       !
CARGU ! VISELY(NCELET! TR !  -> ! "VISCOSITE" PAR CELLULE  DIR Y       !
CARGU ! VISELZ(NCELET! TR !  -> ! "VISCOSITE" PAR CELLULE  DIR Z       !
CARGU ! FLUMAS(NFAC) ! TR ! <-> ! FLUX DE MASSE AUX FACES INTERNES     !
CARGU ! FLUMAB(NFABOR! TR ! <-> ! FLUX DE MASSE AUX FACES DE BORD      !
CARGU ! FEXTX,Y,Z    ! TR !  -> ! FORCE EXTERIEURE GENERANT LA PRESSION!
CARGU !   (NCELET)   !    !     !  HYDROSTATIQUE                       !
CARGU ! DPD.(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR LE GRAD DE P !
CARGU ! DPDXA,Y,Z    ! TR !  -  ! TABLEAU DE TRAVAIL POUR LE GRAD DE P !
CARGU !    (NCELET   !    !     !  AVEC DECENTREMENT AMONT             !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
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
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
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
      INTEGER          INIT   , INC    , IMRGRA , ICCOCG
      INTEGER          NSWRGP , IMLIGP
      INTEGER          IWARNP , IPHYDP , NFECRA
      DOUBLE PRECISION EPSRGP , CLIMGP , EXTRAP
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
      DOUBLE PRECISION PVAR(NCELET), COEFAP(NFABOR), COEFBP(NFABOR)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION VISELX(NCELET), VISELY(NCELET), VISELZ(NCELET)
      DOUBLE PRECISION FLUMAS(NFAC), FLUMAB(NFABOR)
      DOUBLE PRECISION FEXTX(NCELET),FEXTY(NCELET),FEXTZ(NCELET)
      DOUBLE PRECISION DPDX (NCELET),DPDY (NCELET),DPDZ (NCELET)
      DOUBLE PRECISION DPDXA(NCELET),DPDYA(NCELET),DPDZA(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC, II, JJ, IIJ, III
      INTEGER          IDIMTE, ITENSO
      DOUBLE PRECISION PFAC,PIP
      DOUBLE PRECISION DPXF  , DPYF  , DPZF  ,SURFAN, DIST
      DOUBLE PRECISION DIJPFX, DIJPFY, DIJPFZ
      DOUBLE PRECISION DIIPBX, DIIPBY, DIIPBZ
      DOUBLE PRECISION DIJX  , DIJY  , DIJZ
C
C***********************************************************************
C
C
C     C'est penible de passer viscf et viscel qui portent une info
C                similaire
C
C
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IF( INIT.GE.1 ) THEN
        DO IFAC = 1, NFAC
          FLUMAS(IFAC) = 0.D0
        ENDDO
        DO IFAC = 1, NFABOR
          FLUMAB(IFAC) = 0.D0
        ENDDO
      ELSEIF(INIT.NE.0) THEN
        WRITE(NFECRA,1000) INIT
        CALL CSEXIT (1)
      ENDIF
C
C ---> TRAITEMENT DU PARALLELISME
C
      IF(IRANGP.GE.0) CALL PARCOM (PVAR)
C                     ===========
C
C ---> TRAITEMENT DE LA PERIODICITE
C
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     &  ( IDIMTE , ITENSO ,
     &    PVAR   , PVAR   , PVAR  ,
     &    PVAR   , PVAR   , PVAR  ,
     &    PVAR   , PVAR   , PVAR  )
      ENDIF
C
C=======================================================================
C 2.  INCREMENT DU FLUX DE MASSE SS TECHNIQUE DE RECONSTRUCTION
C=======================================================================
C
      IF( NSWRGP.LE.1 ) THEN
C
C     FLUX DE MASSE SUR LES FACETTES FLUIDES
C
        DO IFAC = 1, NFAC
C
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
C
          FLUMAS(IFAC) = FLUMAS(IFAC) + VISCF(IFAC)*(PVAR(II) -PVAR(JJ))
C
        ENDDO
C
C     FLUX DE MASSE SUR LES FACETTES DE BORD
C
        DO IFAC = 1, NFABOR
C
          II = IFABOR(IFAC)
          PFAC = INC*COEFAP(IFAC) +COEFBP(IFAC)*PVAR(II)
C
          FLUMAB(IFAC) = FLUMAB(IFAC) +VISCB(IFAC)*( PVAR(II) -PFAC )
C
        ENDDO
C
      ENDIF
C
C
C=======================================================================
C 3.  INCREMENTATION DU FLUX DE MASSE AVEC TECHNIQUE DE
C         RECONSTRUCTION SI LE MAILLAGE EST NON ORTHOGONAL
C=======================================================================
C
      IF( NSWRGP.GT.1 ) THEN
C
C     CALCUL DU GRADIENT
C
        CALL GRDCEL
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPR(1) , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   FEXTX  , FEXTY  , FEXTZ  ,
     &   PVAR   , COEFAP , COEFBP ,
     &   DPDX   , DPDY   , DPDZ   ,
C        ------   ------   ------
     &   DPDXA  , DPDYA  , DPDZA  ,
     &   RDEVEL , RTUSER , RA     )
C
C ---> TRAITEMENT DU PARALLELISME
C
      IF(IRANGP.GE.0) THEN
        CALL PARCOM (VISELX)
C       ===========
        CALL PARCOM (VISELY)
C       ===========
        CALL PARCOM (VISELZ)
C       ===========
      ENDIF
C
C ---> TRAITEMENT DE LA PERIODICITE
C
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 1
        ITENSO = 0
        CALL PERCOM
C       ===========
     &  ( IDIMTE , ITENSO ,
     &    VISELX , VISELX , VISELX ,
     &    VISELY , VISELY , VISELY ,
     &    VISELZ , VISELZ , VISELZ )
      ENDIF
C
C     FLUX DE MASSE SUR LES FACETTES FLUIDES
C
        DO IFAC = 1, NFAC
C
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
C
          DPXF = 0.5D0*(VISELX(II)*DPDX(II) + VISELX(JJ)*DPDX(JJ))
          DPYF = 0.5D0*(VISELY(II)*DPDY(II) + VISELY(JJ)*DPDY(JJ))
          DPZF = 0.5D0*(VISELZ(II)*DPDZ(II) + VISELZ(JJ)*DPDZ(JJ))
C
          IIJ = IDIJPF-1+3*(IFAC-1)
          DIJPFX = RA(IIJ+1)
          DIJPFY = RA(IIJ+2)
          DIJPFZ = RA(IIJ+3)
C
C---> DIJ = IJ - (IJ.N) N
          DIJX = (XYZCEN(1,JJ)-XYZCEN(1,II))-DIJPFX
          DIJY = (XYZCEN(2,JJ)-XYZCEN(2,II))-DIJPFY
          DIJZ = (XYZCEN(3,JJ)-XYZCEN(3,II))-DIJPFZ
C
          DIST   = RA(IDIST-1+IFAC)
          SURFAN = RA(ISRFAN-1+IFAC)
C
          FLUMAS(IFAC) = FLUMAS(IFAC)
     &     + VISCF(IFAC)*( PVAR(II) -PVAR(JJ) )
     &     + ( DPXF * DIJX
     &     +   DPYF * DIJY
     &     +   DPZF * DIJZ )*SURFAN/DIST
C
        ENDDO
C
C     FLUX DE MASSE SUR LES FACETTES DE BORD
C
        DO IFAC = 1, NFABOR
C
          II = IFABOR(IFAC)
          III = IDIIPB-1+3*(IFAC-1)
          DIIPBX = RA(III+1)
          DIIPBY = RA(III+2)
          DIIPBZ = RA(III+3)
C
          PIP = PVAR(II)
     &        + DPDX(II)*DIIPBX
     &        + DPDY(II)*DIIPBY + DPDZ(II)*DIIPBZ
          PFAC = INC*COEFAP(IFAC) +COEFBP(IFAC)*PIP
C
          FLUMAB(IFAC) = FLUMAB(IFAC) +VISCB(IFAC)*( PIP -PFAC )
C
        ENDDO
C
      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT('ITRMAS APPELE AVEC INIT = ',I10)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
