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
                        SUBROUTINE CLIPKE
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NVAR   , NPHAS  ,
     &   IPHAS  , ICLIP  , IWARNK ,
     &   PROPCE , RTP    )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC CLIPPING DE K ET EPSILON
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NVAR         ! E  !  -> ! NOMBRE DE VARIABLES                  !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE LA PHASE A TRAITER         !
CARGU ! ICLIP        ! E  !  -> ! INDICATEUR = 0 ON UTILISE VISCL0     !
CARGU !              !    !     !            SINON ON UTILISE VISCL    !
CARGU ! IWARNK       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! PROPCE       ! TR ! <-> ! TABLEAUX DES VARIABLES AU PDT COURANT!
CARGU !(NCELET,*     !    !     !                                      !
CARGU ! RTP          ! TR ! <-> ! TABLEAUX DES VARIABLES AU PDT COURANT!
CARGU ! (NCELET     )!    !     !                                      !
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
      INCLUDE "numvar.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NVAR, NCELET, NCEL, NPHAS
      INTEGER          IPHAS, ICLIP, IWARNK
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION RTP(NCELET,NVAR)
C
C VARIABLES LOCALES
C
      INTEGER          ICLPKE,IEL,ICLPK2,ICLPE2
      INTEGER          IVAR,IPP,II,IIKIPH,IIEIPH,IIVISC,IIROMC
      DOUBLE PRECISION XEPMIN,XEPM,XE,XKMIN,XKM,XK, VMIN, VMAX, VAR
      DOUBLE PRECISION EPZ2
C
C***********************************************************************
C
      IIKIPH = IK(IPHAS)
      IIEIPH = IEP(IPHAS)
      IIVISC = IPPROC(IVISCL(IPHAS))
      IIROMC = IPPROC(IROM(IPHAS))
C
C
C Une petite valeur pour eviter des valeurs exactement nulles.
C
      EPZ2 = EPZERO**2
C
C=======================================================================
C ---> Stockage Min et Max pour listing
C=======================================================================
C
      DO II = 1, 2
        IF(II.EQ.1) THEN
          IVAR = IIKIPH
        ELSEIF(II.EQ.2) THEN
          IVAR = IIEIPH
        ENDIF
        IPP  = IPPRTP(IVAR)
C
        VMIN =  GRAND
        VMAX = -GRAND
        DO IEL = 1, NCEL
          VAR = RTP(IEL,IVAR)
          VMIN = MIN(VMIN,VAR)
          VMAX = MAX(VMAX,VAR)
        ENDDO
        IF (IRANGP.GE.0) THEN
          CALL PARMAX (VMAX)
C         ===========
          CALL PARMIN (VMIN)
C         ===========
        ENDIF
        VARMNA(IPP) = VMIN
        VARMXA(IPP) = VMAX
C
      ENDDO
C
C=======================================================================
C ---> Detection des valeurs hors norme "physiques"
C       uniquement pour avertissement
C       ou dans le cas ICLKEP = 1
C=======================================================================
C
      IF (IWARNK.GE.2.OR.ICLKEP(IPHAS).EQ.1) THEN
C
        IF(ICLIP.EQ.1) THEN
C
          XKM = 1296.D0*SQRT(CMU)/ALMAX(IPHAS)**2
          XEPM = 46656.D0*CMU/ALMAX(IPHAS)**4
          ICLPKE = 0
          DO IEL=1,NCEL
            XK = RTP(IEL,IIKIPH)
            XE = RTP(IEL,IIEIPH)
            XKMIN = XKM*(PROPCE(IEL,IIVISC)/PROPCE(IEL,IIROMC))**2
            XEPMIN = XEPM*(PROPCE(IEL,IIVISC)/PROPCE(IEL,IIROMC))**3
            IF(XK.LE.XKMIN.OR.XE.LE.XEPMIN) THEN
              IF(ICLKEP(IPHAS).EQ.1) THEN
                RTP(IEL,IIKIPH)  = XKMIN
                RTP(IEL,IIEIPH) = XEPMIN
              ENDIF
              ICLPKE = ICLPKE + 1
            ENDIF
          ENDDO
C
        ELSEIF(ICLIP.EQ.0) THEN
C
          XKMIN = 1296.D0*SQRT(CMU)/ALMAX(IPHAS)**2*
     &            (VISCL0(IPHAS)/RO0(IPHAS))**2
          XEPMIN = 46656.D0*CMU/ALMAX(IPHAS)**4*
     &            (VISCL0(IPHAS)/RO0(IPHAS))**3
          ICLPKE = 0
          DO IEL=1,NCEL
            XK = RTP(IEL,IIKIPH)
            XE = RTP(IEL,IIEIPH)
            IF(XK.LE.XKMIN.OR.XE.LE.XEPMIN) THEN
              IF(ICLKEP(IPHAS).EQ.1) THEN
                RTP(IEL,IIKIPH)  = XKMIN
                RTP(IEL,IIEIPH) = XEPMIN
              ENDIF
              ICLPKE = ICLPKE + 1
            ENDIF
          ENDDO
C
        ELSE
C
          WRITE(NFECRA,1000)ICLIP
          CALL CSEXIT (1)
C
        ENDIF
C
        IF (IRANGP.GE.0) CALL PARCPT (ICLPKE)
C                        ===========
C
C ---  Impression eventuelle
C
        IF(IWARNK.GE.2) THEN
C
          WRITE(NFECRA,1010)ICLPKE
C
        ENDIF
C
C ---  Stockage nb de clippings pour listing
C
        IF(ICLKEP(IPHAS).EQ.1) THEN
          ICLPMN(IPPRTP(IIKIPH)) = ICLPKE
          ICLPMN(IPPRTP(IIEIPH)) = ICLPKE
        ENDIF
C
      ENDIF
C
C=======================================================================
C ---> Clipping "standard" ICLKEP = 0
C=======================================================================
C
      IF(ICLKEP(IPHAS).EQ.0) THEN
C
        ICLPK2 = 0
        ICLPE2 = 0
        DO IEL = 1, NCEL
          XK = RTP(IEL,IIKIPH)
          XE = RTP(IEL,IIEIPH)
          IF (ABS(XK).LE.EPZ2) THEN
            ICLPK2 = ICLPK2 + 1
            RTP(IEL,IIKIPH) = MAX(RTP(IEL,IIKIPH),EPZ2)
          ELSEIF(XK.LE.0.D0) THEN
            ICLPK2 = ICLPK2 + 1
            RTP(IEL,IIKIPH) = -XK
          ENDIF
          IF (ABS(XE).LE.EPZ2) THEN
            ICLPE2 = ICLPE2 + 1
            RTP(IEL,IIEIPH) = MAX(RTP(IEL,IIEIPH),EPZ2)
          ELSEIF(XE.LE.0.D0) THEN
            ICLPE2 = ICLPE2 + 1
            RTP(IEL,IIEIPH) = -XE
          ENDIF
        ENDDO
C
        IF (IRANGP.GE.0) THEN
          CALL PARCPT (ICLPK2)
C         ===========
          CALL PARCPT (ICLPE2)
C         ===========
        ENDIF
C
C ---  Stockage nb de clippings pour listing
C
        ICLPMN(IPPRTP(IIKIPH)) = ICLPK2
        ICLPMN(IPPRTP(IIEIPH)) = ICLPE2
C
      ENDIF
C
C
C=======================================================================
C ---> Formats
C=======================================================================
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS clipke                           ',/,
     &'@    *********                                               ',/,
     &'@     APPEL DE clipke              AVEC OPTION = ',I10        ,/,
     &'@     Phase : ',I10                                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1010 FORMAT(
     & I10,' VALEURS DU K-EPS AU DELA DES ECHELLES BASEES SUR ALMAX')
C
      RETURN
C
      END
c@z
