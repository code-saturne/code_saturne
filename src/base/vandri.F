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
                        SUBROUTINE VANDRI
C                       *****************
C     --------------------------------------------------------------
     & (  NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NPHAS ,
     &    NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &    ITYPFB , IFABOR , IFAPAT , IDEVEL , ITUSER , IA    ,
     &    XYZCEN , CDGFBO , UETBOR , VISVDR , YPLUSC , PROPCE ,
     &    RDEVEL , RTUSER , RA    )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC IMPOSITION D'UN AMORTISSEMENT DE TYPE VAN DRIEST POUR LA LES
CFONC nut est amortie par (1-exp(-y+/d+))**2 ou d+ est mis par defaut a 26
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE                !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE LA PHASE TRAITEE           !
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS)      !    !     !                                      !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFAPAT       ! TE !  -> ! NO DE FACE DE BRD CODE 5 LA + PROCHE !
CARGU ! (NCELET)     !    !     !    (Rij et echo de paroi      )      !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! XYZCEN       ! TR !  -> ! POINT ASSOCIES AUX VOLUMES DE CONTROL!
CARGU ! (NDIM,NCELET !    !     !                                      !
CARGU ! CDGFBO       ! TR !  -> ! CENTRE DE GRAVITE DES FACES DE BORD  !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! UETBOR       ! TR !  -> ! VITESSE DE FROTTEMENT AU BORD        !
CARGU ! (NFABOR,NPHAS!    !     !  POUR VAN DRIEST EN LES              !
CARGU ! VISVDR(NPHAS)! TR !  -> ! VISCOSITE DYNAMIQUE DS LES CELLULES  !
CARGU ! (NCELET,NPHAS!    !     !  DE BORD APRES AMORTISST DE V DRIEST !
CARGU ! YPLUSC       ! TR !  -> ! VALEUR DE YPLUS AUX CELLULES         !
CARGU ! (NCELET  )   !    !     !    DANS LE CAS ABS(ICDPAR).EQ.1      !
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
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
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NDIM, NCELET , NCEL   , NFAC   , NFABOR, NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IPHAS
      INTEGER          ITYPFB(NFABOR,NPHAS),IFABOR(NFABOR)
      INTEGER          IFAPAT(NCELET)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET),CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION UETBOR(NFABOR,NPHAS), VISVDR(NCELET,NPHAS)
      DOUBLE PRECISION YPLUSC(NCELET)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IEL   , IFAC  , IPCVIS, IPCVST, IPCROM
      DOUBLE PRECISION YPLUS , YMINPA, VISCOS
C
C***********************************************************************
C
      IPCROM = IPPROC(IROM  (IPHAS))
      IPCVIS = IPPROC(IVISCL(IPHAS))
      IPCVST = IPPROC(IVISCT(IPHAS))
C
C     Calcul direct de la distance a la paroi (non compatible parall/perio)
      IF(ABS(ICDPAR).EQ.2) THEN
C
C     En sequentiel, RAS
        IF(IRANGP.LT.0) THEN
          DO IEL = 1, NCEL
            IFAC = IFAPAT(IEL)
            VISCOS = PROPCE(IEL,IPCVIS)/PROPCE(IEL,IPCROM)
            YMINPA = SQRT((CDGFBO(1,IFAC)-XYZCEN(1,IEL))**2
     &           +        (CDGFBO(2,IFAC)-XYZCEN(2,IEL))**2
     &           +        (CDGFBO(3,IFAC)-XYZCEN(3,IEL))**2)
            YPLUS = UETBOR(IFAC,IPHAS) * YMINPA/ VISCOS
            PROPCE(IEL,IPCVST) = PROPCE(IEL,IPCVST)*
     &           (1.0D0-EXP(-YPLUS/CDRIES(IPHAS)))**2
          ENDDO
C     En parallele, on n'amortit que la premiere maille de paroi :
C     dangereux mais a priori inutile (car l'utilisation de
C     ICDPAR=+/-2 en parallele est bloque dans verini)
        ELSE
          WRITE(NFECRA,1000)
          DO IFAC = 1, NFABOR
            IF(ITYPFB(IFAC,IPHAS).EQ.IPAROI .OR.
     &         ITYPFB(IFAC,IPHAS).EQ.IPARUG ) THEN
              IEL = IFABOR(IFAC)
              VISCOS = PROPCE(IEL,IPCVIS)/PROPCE(IEL,IPCROM)
              YMINPA = SQRT((CDGFBO(1,IFAC)-XYZCEN(1,IEL))**2
     &             +        (CDGFBO(2,IFAC)-XYZCEN(2,IEL))**2
     &             +        (CDGFBO(3,IFAC)-XYZCEN(3,IEL))**2)
              YPLUS = UETBOR(IFAC,IPHAS) * YMINPA/ VISCOS
              PROPCE(IEL,IPCVST) = PROPCE(IEL,IPCVST)*
     &             (1.0D0-EXP(-YPLUS/CDRIES(IPHAS)))**2
            ENDIF
          ENDDO
      ENDIF
C
C     Nouveau mode de calcul : c'est plus simple
      ELSEIF(ABS(ICDPAR).EQ.1) THEN
        DO IEL = 1, NCEL
          YPLUS = YPLUSC(IEL)
          PROPCE(IEL,IPCVST) = PROPCE(IEL,IPCVST)*
     &           (1.0D0-EXP(-YPLUS/CDRIES(IPHAS)))**2
        ENDDO
      ENDIF
C
C     Pour les cellules de paroi on remet la viscosite turbulente
C     qui avait ete amortie dans clptur et qui a servi a calculer
C     les conditions aux limites
      DO IEL = 1, NCEL
        IF (VISVDR(IEL,IPHAS).GT.-900.D0)
     &       PROPCE(IEL,IPCVST) = VISVDR(IEL,IPHAS)
      ENDDO
C
C--------
C FORMATS
C--------
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : DANS LE CAS DE LA LES AVEC AMORTISSEMENT    ',/,
     &'@    *********                                               ',/,
     &'@    L''AMORTISSEMENT DE VAN DRIEST N''EST FAIT QUE SUR LA   ',/,
     &'@    PREMIERE CELLULE A LA PAROI EN CAS DE PARALLELISME      ',/,
     &'@                                                            ',/,
     &'@  Le calcul se poursuit.                                    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C----
C FIN
C----
C
C
      RETURN
      END
c@z
