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
                        SUBROUTINE MODPAR
C                       *****************
C     -------------------------------------------------------------
     & ( NTCABS , NTMABS )
C     -------------------------------------------------------------
C***********************************************************************
c@foncb
CFONC
CFONC FONCTION :
CFONC --------
CFONC           MODIFICATION DE NTMABS AU COURS DU CALCUL
CFONC           POUR ARRET INTERACTIF
CFONC
c@fonce
C-----------------------------------------------------------------------
C                       ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NTCABS       ! E  !  -> ! NUMERO ABSOLU DU PAS DE TEMPS COURANT!
CARGU ! NTMABS       ! E  !  -> ! NUMERO ABSOLU DU PAS DE TEMPS FINAL  !
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
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NTCABS , NTMABS
C
C VARIABLES LOCALES
C
      INTEGER          IRANGS, LNG, ITMP(1)
      LOGICAL EXSTP
C
C***********************************************************************
C
C On passe ici en sequentiel, ou en parallele avec le proc 0.
      IF (IRANGP.LE.0) THEN
C---> ARRET D'URGENCE
C
        INQUIRE (FILE=FICSTP,EXIST=EXSTP)
C
C       SI UN FICHIER FICSTP EXISTE
C
        IF (EXSTP) THEN
C
C         LIRE LE NOMBRE MAX D'ITERATIONS (ABSOLU)
C
          OPEN (FILE=FICSTP,UNIT=IMPSTP)
          READ (IMPSTP,*,ERR=5200,END=5200)
 5200     READ (IMPSTP,*,ERR=5100,END=5100)NTMABS
 5100     CONTINUE
          CLOSE (IMPSTP,STATUS='DELETE')
C
C         COMPARER LE TEMPS ECOULE ET LE TEMPS MAX
C          MODIFIER FICSTP SI BESOIN
C
          IF(NTCABS.GT.NTMABS)THEN
            NTMABS = NTCABS
          ENDIF
C
C         SORTIES
C
          WRITE (NFECRA,1000) NTCABS,NTMABS
C
          OPEN (FILE=FICSTP//'.mod',UNIT=IMPSTP)
          WRITE (IMPSTP,1000) NTCABS,NTMABS
          CLOSE (IMPSTP)
        ENDIF
C
      ENDIF
C
C En parallele, bcast.
      IF(IRANGP.GE.0) THEN
        IRANGS  = 0
        LNG     = 1
        ITMP(1) = NTMABS
        CALL PARBCI(IRANGS,LNG,ITMP)
        NTMABS = ITMP(1)
      ENDIF
C
 1000 FORMAT(/,
     &'*************************************************************',/,
     &'            NTCABS COURANT  = ',I10                           ,/,
     &'            NTMABS RESET TO = ',I10                           ,/,
     &'*************************************************************',/,
     &                                                                /)
C
      END
c@z
