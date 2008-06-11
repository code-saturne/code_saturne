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
                        SUBROUTINE DISTPR
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ITYPFB ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DISTPA ,
     &   VISCF  , VISCB  , DAM    , XAM    , SMBDP  , ROVSDP ,
     &   RTPDP  , COEFAD , COEFBD ,
     &   W1     , W2     , W3     , W4     , W5     , W6     , W7     ,
     &   W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC CALCUL DE LA DISTANCE A LA PAROI PAR INVERSION DE LA SOLUTION 3D
CFONC    DE L'EQUATION DE DIFFUSION D'UN SCALAIRE
CFONC
CFONC   On resout
CFONC     div[grad(T)] = -1
CFONC       avec :
CFONC       T(bord)   = 1 en paroi
CFONC       grad(T).n = 0 ailleurs
CFONC
CFONC   La distance a la paroi vaut alors :
CFONC
CFONC    d ~ -|grad(T)|+[grad(T).grad(T)+2.T]^(1/2)
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
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS      )!    !     !                                      !
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
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! DISTPA(NCELET! TR ! <-  ! TAB DES DISTANCES A LA PAROI         !
CARGU ! VISCF(NFAC)  ! TR !  -  ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! SMBRP(NCELET)! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! ROVSDP(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! RTPDP(NCELET)! TR !  -  ! VAR DE TRAVAIL DU SCLAIRE DIFFUSE    !
CARGU ! DRTP(NCELET) ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! COEFAD,COEFBD! TR !  -  ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR)    !    !     ! FACES DE BORD DU SCALAIRE DIFFUSE    !
CARGU ! W1...9(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
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
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "coincl.h"
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
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ITYPFB(NFABOR,NPHAS)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DISTPA(NCELET), VISCF (NFAC)  , VISCB (NFABOR)
      DOUBLE PRECISION DAM   (NCELET), XAM   (NFAC,2)
      DOUBLE PRECISION SMBDP (NCELET), ROVSDP(NCELET)
      DOUBLE PRECISION RTPDP (NCELET)
      DOUBLE PRECISION COEFAD(NFABOR), COEFBD(NFABOR)
      DOUBLE PRECISION W1    (NCELET), W2    (NCELET), W3    (NCELET)
      DOUBLE PRECISION W4    (NCELET), W5    (NCELET), W6    (NCELET)
      DOUBLE PRECISION W7    (NCELET), W8    (NCELET), W9    (NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          NDIRCP, ICONVP, IDIFFP, ISYM
      INTEGER          IPOL  , IRESLP, IPP
      INTEGER          NITERF, ICYCLE, NCYMXP, NITMFP
      INTEGER          IINVPE, IDIMTE, ITENSO
      INTEGER          ISQRT , IEL   , IFAC  , IPHAS
      INTEGER          INC   , ICCOCG, IPHYDP, IVAR
      INTEGER          ISWEEP, NITTOT, IDTVA0
C
      DOUBLE PRECISION RELAXP, THETAP, RNORM, RESIDU, RNOINI
      DOUBLE PRECISION DISMAX, DISMIN
C***********************************************************************
C
C
C
C=======================================================================
C 1. INITIALISATIONS
C=======================================================================
C     Memoire
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C     On travaille sur la phase 1
      IPHAS = 1
C     Nombre d'iteration totale pour l'inversion
      NITTOT = 0
C
C     La distance a la paroi est initialisee a 0 pour la reconstruction
C
      DO IEL = 1, NCELET
        DISTPA(IEL) = 0.D0
        RTPDP(IEL)  = 0.D0
      ENDDO
C
C=======================================================================
C 2.CONDITIONS LIMITES
C=======================================================================
C
C     Conditions aux limites pour le scalaire resolu T
C       Dirichlet a 0 en paroi
C       Neumann nul ailleurs
C     On test aussi la pressence d'un Dirichlet

      NDIRCP = 0
C
      DO IFAC = 1, NFABOR
        IF(ITYPFB(IFAC,IPHAS).EQ.IPAROI .OR.
     &     ITYPFB(IFAC,IPHAS).EQ.IPARUG) THEN
          COEFAD(IFAC) = 0.0D0
          COEFBD(IFAC) = 0.0D0
          NDIRCP = 1
        ELSE
          COEFAD(IFAC) = 0.0D0
          COEFBD(IFAC) = 1.0D0
        ENDIF
      ENDDO
C
C=======================================================================
C 3. PREPARATION DU SYSTEME A RESOUDRE
C=======================================================================
C
C -- Diagonale
C
      DO IEL = 1, NCEL
        ROVSDP(IEL) = 0.D0
      ENDDO
C
C -- Diffusion aux faces
C
      DO IEL = 1, NCEL
        W1(IEL) = 1.D0
      ENDDO
C
      CALL VISCFA
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IMVISF ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     ,
     &   VISCF  , VISCB  ,
     &   RDEVEL , RTUSER , RA     )
C
      ICONVP = 0
      IDIFFP = 1
      ISYM   = 1
      THETAP = 1.D0
C
      CALL MATRIX
C     ===========
     & ( NCELET , NCEL   , NFAC   , NFABOR ,
     &   ICONVP , IDIFFP , NDIRCP ,
     &   ISYM   , NFECRA ,
     &   THETAP ,
     &   IFACEL , IFABOR ,
     &   COEFBD , ROVSDP ,
     &   VISCF  , VISCB  , VISCF  , VISCB  ,
     &   DAM    , XAM    )
C
C -- Second membre
C
      DO IEL = 1, NCEL
        SMBDP(IEL)  = VOLUME(IEL)
      ENDDO
C=======================================================================
C 4. RESOLUTION DU SYSTEME
C=======================================================================
      IPP = 1
      NOMVAR(IPP) = 'DisParoi'
      IPOL   = 0
      IRESLP = 0
C Pas de multigrille (NCYMXP,NITMFP sont arbitraies)
      NCYMXP = 100
      NITMFP = 10
C Periodicite
      IINVPE = 0
      IF(IPERIO.EQ.1) IINVPE = 1
      ISQRT = 1
C
C -- Boucle de reconstruction
C
C Si NSWRSY = 1, on doit faire 2 inversions
      DO ISWEEP = 0, NSWRSY
C
        CALL PRODSC(NCELET,NCEL,ISQRT,SMBDP,SMBDP,RNORM)
        IF (IWARNY.GE.2) THEN
           WRITE(NFECRA,5000) NOMVAR(IPP),ISWEEP,RNORM
        ENDIF
        IF (ISWEEP.LE.1) RNOINI = RNORM
C Test de convergence
        IF(RNORM.LE.10.D0*EPSILY*RNOINI) GOTO 100
C
        DO IEL = 1, NCELET
          RTPDP(IEL) = 0.D0
        ENDDO
C
        CALL INVERS
C       ===========
     & ( NOMVAR(IPP)     , IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   ISYM   , IPOL   , IRESLP , NITMAY , IMGRPY ,
     &   NCYMXP , NITMFP ,
     &   IWARNY , NFECRA , NITERF , ICYCLE , IINVPE ,
     &   EPSILY , RNORM  , RESIDU ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DAM    , XAM    , SMBDP  , RTPDP  ,
     &   W3     , W4     , W5     , W6     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
        NITTOT = NITTOT + NITERF
        DO IEL = 1, NCEL
          DISTPA(IEL) = DISTPA(IEL) + RTPDP(IEL)
        ENDDO
C
C    - Echange pour le parallelisme
C
        IF(IRANGP.GE.0) THEN
          CALL PARCOM (RTPDP)
C         ===========
        ENDIF
C
C    - Echange pour la periodicite
        IF(IPERIO.EQ.1) THEN
          IDIMTE = 0
          ITENSO = 0
          CALL PERCOM
C         ===========
     &      ( IDIMTE , ITENSO ,
     &        RTPDP  , RTPDP  , RTPDP  ,
     &        RTPDP  , RTPDP  , RTPDP  ,
     &        RTPDP  , RTPDP  , RTPDP  )
        ENDIF
C
        IF(ISWEEP.LT.NSWRSY) THEN
          INC    = 0
          ICCOCG = 1
          IVAR = 0
          IDTVA0 = 0
          RELAXP = 1.D0
C
          CALL BILSC2
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IDTVA0 , IVAR   , ICONVP , IDIFFP , NSWRGY , IMLIGY , IRCFLY ,
     &   ISCHCY , ISSTPY , INC    , IMRGRA , ICCOCG ,
     &   IPP    , IWARNY ,
     &   BLENCY , EPSRGY , CLIMGY , EXTRAY , RELAXP , THETAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   RTPDP  , RTPDP  , COEFAD , COEFBD , COEFAD , COEFBD ,
     &   VISCF  , VISCB  , VISCF  , VISCB  ,
     &   SMBDP  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
        ENDIF
      ENDDO
C
 100  CONTINUE
C
C On travail ensuite sur RTPDP pour calculer DISTPA
      DO IEL=1,NCEL
        RTPDP(IEL)  = DISTPA(IEL)
      ENDDO
C
C=======================================================================
C 5. CALCUL DE LA DISTANCE A LA PAROI
C=======================================================================v
C
C    - Echange pour le parallelisme
C
      IF(IRANGP.GE.0) THEN
        CALL PARCOM (RTPDP)
C       ===========
      ENDIF
C
C    - Echange pour la periodicite
C
      IF(IPERIO.EQ.1) THEN
        IDIMTE = 0
        ITENSO = 0
        CALL PERCOM
C       ===========
     &      ( IDIMTE , ITENSO ,
     &        RTPDP  , RTPDP  , RTPDP  ,
     &        RTPDP  , RTPDP  , RTPDP  ,
     &        RTPDP  , RTPDP  , RTPDP  )
C
      ENDIF
C
C    - Calcul du gradient
C
      INC    = 1
      ICCOCG = 1
      IPHYDP = 0
      IVAR   = 0
C
      CALL GRDCEL
C     ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR   , IMRGRA , INC    , ICCOCG , NSWRGY , IMLIGY , IPHYDP ,
     &   IWARNY , NFECRA , EPSRGY , CLIMGY , EXTRAY ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W2     , W3     ,
     &   RTPDP  , COEFAD , COEFBD ,
     &   W4     , W5     , W6     ,
     &   W7     , W8     , W9     ,
     &   RDEVEL , RTUSER , RA     )
C
      DO IEL = 1, NCEL
        W1(IEL) = W4(IEL)**2.D0+W5(IEL)**2.D0+W6(IEL)**2.D0
        IF(W1(IEL)+2.D0*RTPDP(IEL).GT.0.D0) THEN
          DISTPA(IEL) = - SQRT(W1(IEL))
     &                  + SQRT(W1(IEL)+2.D0*RTPDP(IEL))
        ELSE
          WRITE(NFECRA,8000)IEL, XYZCEN(1,IEL)
     &                     ,XYZCEN(2,IEL),XYZCEN(3,IEL)
        ENDIF
      ENDDO
C
C=======================================================================
C 6. CALCUL DES BORNES ET IMPRESSIONS
C=======================================================================
        DISMAX = -GRAND
        DISMIN =  GRAND
C
        DO IEL = 1, NCEL
            DISMIN = MIN(DISTPA(IEL),DISMIN)
            DISMAX = MAX(DISTPA(IEL),DISMAX)
        ENDDO
C
        IF (IRANGP.GE.0) THEN
          CALL PARMIN(DISMIN)
          CALL PARMAX(DISMAX)
        ENDIF
C
C     Impressions
C
        WRITE(NFECRA,1000)DISMIN, DISMAX, NITTOT
C
C
C=======================================================================
C 7. FORMATS
C=======================================================================
C
 1000 FORMAT(
     &'                                                             ',/,
     &' ** DISTANCE A LA PAROI                                      ',/,
     &'    -------------------                                      ',/,
     &'                                                             ',/,
     &'   Distance min = ',E14.5    ,'  Distance max = ',E14.5       ,/,
     &'                                                             ',/,
     &'     (Calcul de la distance realise en ',I10   ,' iterations)',/)
C
 5000 FORMAT(1X,A8,' : SWEEP = ',I5,' NORME SECOND MEMBRE = ',E14.6)
C
 8000   FORMAT(
     &'@                                                            ',/,
     &'@ @@ ATTENTION : Calcul de la distance a la paroi            ',/,
     &'@    *********                                               ',/,
     &'@  La variable associee ne converge pas a la cellule ',I10    ,/,
     &'@       Coord X      Coord Y      Coord Z                    ',/,
     &'@ ',3E13.5                                                    ,/)
C
      RETURN
      END
c@z
