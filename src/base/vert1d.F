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
                        SUBROUTINE VERT1D
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NFABOR , NFPT1D , IAPPEL ,
     &   IFPT1D , NPPT1D , ICLT1D , IA     ,
     &   RGPT1D , EPPT1D ,
     &   XLMBT1 , RCPT1D , DTPT1D , RA     )
C
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC   VERIFICATION DES DONNEES UTILISATEUR DU MODULE THERMIQUE 1D EN PAROI
CFONC
CFONC IAPPEL = 1 (un seul appel a l'initialisation) :
CFONC              VERIFICATION DU NOMBRE DE CELLULES OU L'ON IMPOSE UNE PAROI
CFONC              VERIFICATION DE ISUIT1
CFONC
CFONC IAPPEL = 2 (un seul appel a l'initialisation) :
CFONC              VERIFICATION DU REPERAGE DES CELLULES OU L'ON IMPOSE
CFONC                                                                UNE PAROI
CFONC              VERIFICATION DES DONNEES RELATIVE AU MAILLAGE
CFONC
CFONC IAPPEL = 3 (appel a chaque pas de temps) :
CFONC              VERIFICATION DES TYPES DE C.L. EN PAROI EXTERIEURE
CFONC              VERIFICATION DES VALEURS DES COEFFICIENTS PHYSIQUE DU CALCUL
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NFML         ! E  !  -> ! NOMBRE DE FAMILLES D ENTITES         !
CARGU ! NPRFML       ! E  !  -> ! NOMBRE DE PROPRIETESE DES FAMILLES   !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NFPT1D       ! E  !  -> ! NOMBRE DE FACES AVEC MODULE THERM 1D !
CARGU ! IFPT1D       ! E  !  -> ! NUMERO DE LA FACE EN TRAITEMENT      !
CARGU !              !    !     ! THERMIQUE EN PAROI                   !
CARGU ! NPPT1D       ! E  !  -> ! NOMBRE DE POINTS DE DISCRETISATION   !
CARGU !              !    !     ! DANS LA PAROI                        !
CARGU ! EPPT1D       ! R  !  -> ! EPAISSEUR DE LA PAROI                !
CARGU ! RGPT1D       ! R  !  -> ! RAISON DU MAILLAGE                   !
CARGU ! ICLT1D       ! E  !  -> ! TYPE DE CONDITION LIMITE             !
CARGU ! XLMBT1       ! R  !  -> ! DIFFUSIVITE THERMIQUE                !
CARGU ! RCPT1D       ! R  !  -> ! ROCP                                 !
CARGU ! DTPT1D       ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IAPPEL       ! E  !  -> ! INDIQUE LES DONNES A RENVOYER        !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMFBR       ! TE !  -> ! NUMERO DE FAMILLE D'UNE FACE DE BORD !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IPRFML       ! TE !  -> ! PROPRIETES D'UNE FAMILLE             !
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
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
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
C
C***********************************************************************
C
C Arguments
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NFABOR , NFPT1D
      INTEGER          IAPPEL
C
      INTEGER          IFPT1D(NFPT1D) , NPPT1D(NFPT1D) , ICLT1D(NFPT1D)
      INTEGER          IA(*)
C
      DOUBLE PRECISION EPPT1D(NFPT1D) , RGPT1D(NFPT1D)
      DOUBLE PRECISION XLMBT1(NFPT1D) , RCPT1D(NFPT1D) , DTPT1D(NFPT1D)
      DOUBLE PRECISION RA(*)
C
C Variables locales
      INTEGER          IDEBIA , IDEBRA
      INTEGER          II, IFAC
C
C***********************************************************************
C
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
      IF(IAPPEL.EQ.1) THEN
C
        IF ((NFPT1D.LT.0).OR.(NFPT1D.GT.NFABOR)) THEN
          WRITE(NFECRA,1000)NFABOR,NFPT1D
          CALL CSEXIT(1)
        ENDIF
        IF (ISUIT1.NE.0 .AND. ISUIT1.NE.1) THEN
          WRITE(NFECRA,1010)ISUIT1
          CALL CSEXIT(1)
        ENDIF
C
      ELSE IF (IAPPEL.EQ.2) THEN
C
        DO II = 1, NFPT1D
          IF (IFPT1D(II).LT.0 .OR. IFPT1D(II).GT.NFABOR) THEN
            WRITE(NFECRA,2000)NFABOR,II,IFPT1D(II)
            CALL CSEXIT(1)
          ENDIF
        ENDDO
C
        DO II = 1, NFPT1D
          IFAC = IFPT1D(II)
          IF (NPPT1D(II).LE.0) THEN
            WRITE(NFECRA,2010)II,IFAC,NPPT1D(II)
            CALL CSEXIT(1)
          ENDIF
          IF (EPPT1D(II).LE.0.D0) THEN
            WRITE(NFECRA,2020)'EPPT1D','EPPT1D',II,EPPT1D(II),IFAC
            CALL CSEXIT(1)
          ENDIF
          IF (RGPT1D(II).LE.0.D0) THEN
            WRITE(NFECRA,2020)'RGPT1D','RGPT1D',II,RGPT1D(II),IFAC
            CALL CSEXIT(1)
          ENDIF
        ENDDO
C
      ELSE IF (IAPPEL.EQ.3) THEN
C
        DO II = 1, NFPT1D
          IFAC = IFPT1D(II)
          IF (ICLT1D(II).NE.1 .AND. ICLT1D(II).NE.3) THEN
            WRITE(NFECRA,3000)'ICLT1D','ICLT1D',II,ICLT1D(II),IFAC
            CALL CSEXIT(1)
          ENDIF
          IF (XLMBT1(II).LE.0.D0) THEN
            WRITE(NFECRA,2020)'XLMBT1','XLMBT1',II,XLMBT1(II),IFAC
            CALL CSEXIT(1)
          ENDIF
          IF (RCPT1D(II).LE.0.D0) THEN
            WRITE(NFECRA,2020)'RCPT1D','RCPT1D',II,RCPT1D(II),IFAC
            CALL CSEXIT(1)
          ENDIF
          IF (DTPT1D(II).LE.0.D0) THEN
            WRITE(NFECRA,2020)'DTPT1D','DTPT1D',II,DTPT1D(II),IFAC
            CALL CSEXIT(1)
          ENDIF
        ENDDO
C
      ENDIF
C
C----
C FORMATS
C----
C
 1000   FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    MODULE THERMIQUE 1D EN PAROI                            ',/,
     &'@                                                            ',/,
     &'@    NFPT1D DOIT ETRE POSITIF ET INFERIEUR A NFABOR          ',/,
     &'@    ON A ICI                                                ',/,
     &'@       NFABOR=',I10                                          ,/,
     &'@       NFPT1D=',I10                                          ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier uspt1d.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1010   FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    MODULE THERMIQUE 1D EN PAROI                            ',/,
     &'@                                                            ',/,
     &'@    ISUIT1 DOIT ETRE UN ENTIER EGAL A 0 OU 1                ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier uspt1d.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2000   FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    MODULE THERMIQUE 1D EN PAROI                            ',/,
     &'@                                                            ',/,
     &'@    LE TABLEAU IFPT1D(II) DOIT RENVOYER A UN NUMERO DE FACE ',/,
     &'@                                                     DE BORD',/,
     &'@    ON A ICI                                                ',/,
     &'@       NFABOR            =',I10                              ,/,
     &'@       IFPT1D(',I10,   ')=',I10                              ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier uspt1d.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2010   FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    MODULE THERMIQUE 1D EN PAROI                            ',/,
     &'@                                                            ',/,
     &'@    LE TABLEAU NPPT1D(II) DOIT RENVOYER A UN ENTIER POSITIF ',/,
     &'@    ON A ICI                                                ',/,
     &'@       NPPT1D(',I10,   ')=',I10                              ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier uspt1d.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2020   FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    MODULE THERMIQUE 1D EN PAROI                            ',/,
     &'@                                                            ',/,
     &'@    LE TABLEAU ',A6,' DOIT RENVOYER A UN REEL > 0           ',/,
     &'@    ON A ICI                                                ',/,
     &'@       ',A6,'(',I10,   ')=',E14.5                            ,/,
     &'@       (FACE DE BORD NUMERO ',I10   ,')                     ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier uspt1d.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3000   FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    MODULE THERMIQUE 1D EN PAROI                            ',/,
     &'@                                                            ',/,
     &'@    LE TABLEAU ',A6,' NE PEUT CONTENIR QUE LES VALEURS      ',/,
     &'@                                                     1 OU 3',/,
     &'@    ON A ICI                                                ',/,
     &'@       ',A6,'(',I10,   ')=',I10                              ,/,
     &'@       (FACE DE BORD NUMERO ',I10   ,')                     ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier uspt1d.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      RETURN
C
      END
c@z

