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
                        SUBROUTINE CPPDFR
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   ,
     &   INDPDF , SC     , SD     , SP2M   ,
     &   DSC    , DSD    , SDEB   , SFIN   ,
     &   HREC   )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  CALCUL DES PARAMETRES DE LA PDF
CFONC  PDF RECTANGLE - PICS DE DIRAC CENTREE PPl - AE
CFONC
CFONC  LE RESULTAT EST :
CFONC  ---------------
CFONC     CALCUL DES PARAMETRES ASSOCIES AUX FONCTIONS RECTANGLE - DIRAC
CFONC
CFONC          DSC contient   le Dirac en SC
CFONC          DSD  - - - - - le Dirac en SD
CFONC          SDEB - - - - - l'abcisse de debut du rectangle
CFONC          SFIN - - - - - - - - - - -  fin - - - - - - -
CFONC          HREC - - - - - la hauteur du rectangle
CFONC          INDPDF indique le passage ou non par la pdf
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! SP2M         ! TR !  -> ! VARIANCE DE S                        !
CARGU ! SC           ! TR !  -> ! BORNE MIN DE S                       !
CARGU ! SD           ! TR !  -> ! BORNE MAX DE S                       !
CARGU ! DSC          ! TR !  <- ! DIRAC EN C                           !
CARGU ! DSD          ! TR !  <- ! DIRAC EN D                           !
CARGU ! SDEB         ! TR !  <- ! ABSCISSE DEBUT RECTANGLE             !
CARGU ! SFIN         ! TR !  <- ! ABSCISSE FIN RECTANGLE               !
CARGU ! HREC         ! TR !  <- ! HAUTEUR RECTANGLE                    !
CARGU ! INDPDF       ! TR !  <- ! INDICATEUR PASSAGE OU NON PAR PDF    !
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
C**********************************************************************
C     DONNEES EN COMMON
C**********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "pointe.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET , NCEL
      INTEGER          INDPDF(NCELET)
C
      DOUBLE PRECISION SP2M(NCELET)
      DOUBLE PRECISION SC(NCELET)   , SD(NCELET)
      DOUBLE PRECISION DSC(NCELET)  , DSD(NCELET)
      DOUBLE PRECISION SDEB(NCELET) , SFIN(NCELET)
      DOUBLE PRECISION HREC(NCELET)

C
C VARIABLES LOCALES
C
      INTEGER          IEL
      INTEGER          N1     , N2     , N3     , N4     , N5
#ifdef DEBUG
      DOUBLE PRECISION MOM0   , MOM1   , MOM2
#endif
      DOUBLE PRECISION T1     , T2     , T3     , EPSI
C
C
C***********************************************************************
C
C=======================================================================
C 1.  CALCULS PRELIMINAIRES
C=======================================================================
C
      EPSI = 1.D-06
C Parametre relatif a la variance
      T1 = 1.D-05
C Parametre relatif a la moyenne
      T2 = 5.D-04
C Parametre supplementaire
      T3   = SQRT(3.D0*T1)
C
      DO IEL = 1, NCEL
        IF ( INDPDF(IEL).GT.0 ) THEN
          IF ( (SP2M(IEL).LT.T1) .OR. (ABS(SD(IEL)+SC(IEL)).LT.T2) )
     &      INDPDF(IEL) = 0
        ENDIF
      ENDDO
C
C
C=======================================================================
C 2.  CALCUL DES PARAMETRES DE LA FONCTION DENSITE DE PROBABILITE
C
C     DONNEES POUR LES FONCTIONS RECTANGLES
C       SC   contient la borne inferieure de S
C       SD   - - - -  la borne superieure de S
C       SP2M - - - -  la variance de S
C       SM   - - - -  la moyenne de S = 0
C
C     SC -- DEB ---- 0 ---- FIN ---------- SD
C
C     PARAMETRES ASSOCIES AUX FONCTIONS RECTANGLES
C       DSC    contient le pic de Dirac en C
C       DSD    - - - -  le pic de Dirac en D
C       SDEB   - - - -  l'abcisse de debut du rectangle
C       SFIN   - - - - - - - - - - - fin   - - - - - -
C       HREC   - - - -  la hauteur du rectangle
C       INDPDF - - - -  le type de PDF : 2,3,11,12,13
C                   (= 0 si pas de passage par les pdf)
C
C=======================================================================
C
      DO IEL = 1, NCEL
C
        IF (INDPDF(IEL).GT.0) THEN
C
          IF ( ( SD(IEL)  .GE.(-SC(IEL)) .AND.
     &         SP2M(IEL).LE.(SC(IEL)**2)/3.D0  ) .OR.
     &         ( SD(IEL)  .LT.(-SC(IEL)) .AND.
     &         SP2M(IEL).LE.(SD(IEL)**2)/3.D0  )       ) THEN
C
C --> Rectangle seul
C
            HREC(IEL) = SQRT(3.D0*SP2M(IEL))
C     ATTENTION HREC n'est pas la hauteur du rectangle
C                    mais un intermediaire de calcul
            DSC(IEL) = ZERO
            DSD(IEL) = ZERO
            SDEB(IEL) = - HREC(IEL)
            SFIN(IEL) = + HREC(IEL)
          ELSEIF ( SD(IEL)  .GE.(-SC(IEL)) .AND.
     &           SP2M(IEL).LE.(-SC(IEL)/3.D0*(SC(IEL)+2.D0*SD(IEL)))
     &           ) THEN
C
C --> Rectangle et un Dirac en SC
C
            SDEB(IEL) = SC(IEL)
            SFIN(IEL) = (3.D0*SP2M(IEL)+SDEB(IEL)**2)/(-2.D0*SDEB(IEL))
            DSC(IEL)= (SFIN(IEL)+SDEB(IEL))/(SFIN(IEL)-SDEB(IEL))
            DSD(IEL) = ZERO
          ELSEIF ( SD(IEL)  .LT.(-SC(IEL)) .AND.
     &           SP2M(IEL).LT. (-SD(IEL)/3.D0*(2.D0*SC(IEL)+SD(IEL)))
     &           ) THEN
C
C --> Rectangle et un Dirac en SD
C
            SFIN(IEL) = SD(IEL)
            SDEB(IEL) = (3.D0*SP2M(IEL)+SFIN(IEL)**2)/(-2.D0*SFIN(IEL))
            DSD(IEL) = (-SFIN(IEL)-SDEB(IEL))/(SFIN(IEL)-SDEB(IEL))
            DSC(IEL)= ZERO
          ELSE
C
C --> Rectangle et deux Diracs en SC et SD
C
            SDEB(IEL) = SC(IEL)
            SFIN(IEL) = SD(IEL)
            DSD(IEL) = (3.D0*SP2M(IEL)+SDEB(IEL)**2
     &           +2.D0*SDEB(IEL)*SFIN(IEL))
     &           / ((SFIN(IEL)-SDEB(IEL))**2)
            DSC(IEL)= DSD(IEL)+(SFIN(IEL)+SDEB(IEL))
     &           /(SFIN(IEL)-SDEB(IEL))
C
          ENDIF
C
          IF ( ABS(SFIN(IEL)-SDEB(IEL)).GT.EPSI ) THEN
            HREC(IEL) = (1.D0-DSC(IEL)-DSD(IEL))/(SFIN(IEL)-SDEB(IEL))
          ELSE
            SDEB(IEL) = MIN(SD(IEL),MAX(SC(IEL),-T3))
            SFIN(IEL) = MIN(SD(IEL),MAX(SC(IEL),+T3))
            HREC(IEL) = (1.D0-DSC(IEL)-DSD(IEL))/(SFIN(IEL)-SDEB(IEL))
          ENDIF
C
        ENDIF
C
      ENDDO
C
C
C=======================================================================
C 3. IMPRESSION
C=======================================================================
C
      N1 = 0
      N2 = 0
      N3 = 0
      N4 = 0
      N5 = 0
C
      DO IEL = 1, NCEL
        IF ( INDPDF(IEL).NE.0  ) N1 = N1+1
        IF ( INDPDF(IEL).EQ.3  ) N2 = N2+1
        IF ( INDPDF(IEL).EQ.12 ) N3 = N3+1
        IF ( INDPDF(IEL).EQ.13 ) N4 = N4+1
        IF ( INDPDF(IEL).EQ.11 ) N5 = N5+1
C
C ---- Test
#ifdef DEBUG
        IF ( INDPDF(IEL).NE.0  ) THEN
          MOM0 = HREC(IEL)*(SFIN(IEL)-SDEB(IEL))
     &         + DSC(IEL)+DSD(IEL)
          IF ( ABS(MOM0-1.D0).GT.EPSI ) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'MOM0 VC VA', MOM0, 1.D0
          ENDIF
          MOM1 = HREC(IEL)*(SFIN(IEL)**2-SDEB(IEL)**2)/2.D0
     &         + DSC(IEL)*SC(IEL)+DSD(IEL)*SD(IEL)
          IF ( MOM1.GT.EPSI ) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'MOM1 VC VA', MOM1, 0.D0
          ENDIF
          MOM2 = HREC(IEL)*(SFIN(IEL)**3-SDEB(IEL)**3)/3.D0
     &         + DSC(IEL)*SC(IEL)**2+DSD(IEL)*SD(IEL)**2
          IF ( ABS(MOM2-SP2M(IEL)).GT.EPSI ) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'MOM2 VC VA',  MOM2, SP2M(IEL)
          ENDIF
        ENDIF
#endif
C ---- Fin Test
C
      ENDDO
C
      IF (IRANGP.GE.0) THEN
        CALL PARCPT (N1)
C       ===========
        CALL PARCPT (N2)
C       ===========
        CALL PARCPT (N3)
C       ===========
        CALL PARCPT (N4)
C       ===========
        CALL PARCPT (N5)
C       ===========
      ENDIF
C
      WRITE(NFECRA,1000) NCEL, N1, N2, N3, N4, N5
C
C
C----
C FORMATS
C----
C
 1000 FORMAT (/,
     *'CONTROLE DES PARAMETRES DANS CPPDFR.F',/,
     *'======================================',/,
     *' Nb de points de calculs                                  = ',
     *   I6,/,
     *' Nb de points turbulents (passage par les PDF)            = ',
     *   I6,/,
C ..v.7..1....v....2....v....3....v....4....v....5....v....6....v....7.I
     *' Nb de points turbulents pour lesquels support PDF = I4M  = ',
     *   I6,/,
     *' Nb de points turbulents pour lesquels C app. [I4,L3]     = ',
     *   I6,/,
     *' - - - - - - - - - - - - pour lesquels C app. [I4,L5]     = ',
     *   I6,/,
     *' - - - - - - - - - - - - pour lesquels C app. [L5,I3max]  = ',
     *   I6)
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
