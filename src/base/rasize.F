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
                        SUBROUTINE RASIZE
C                       *****************
C     -------------------------------------------------------------
     & ( SOUSPG , MEMRDP )
C     -------------------------------------------------------------
C***********************************************************************
C
C  FONCTION :
C  --------
c@foncb
CFONC
CFONC   TEST POUR VERIFIER SI LA RESERVATION EN MEMOIRE
CFONC    CENTRALE DANS IA FAITE PAR SOUSPG NE DEPASSE PAS MTOTRA
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! SOUSPG       ! A  !  -> ! SOUS PROGRAMME APPELANT              !
CARGU ! MEMRDP       ! E  !  -> ! INDICE DE LA DERNIERE CASE UTILISEE  !
CARGU !              !    !     !  DANS RA                             !
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
      INCLUDE "memoir.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
      CHARACTER SOUSPG*6
      INTEGER   MEMRDP
C
C---> REMANENCE ASSUREE POUR SPGMAX ET MEMIMX
C
      CHARACTER     SPGMAX*6
      DATA          SPGMAX /'DEFAUT'/
      INTEGER       MEMRMX
      DATA          MEMRMX /0/
      SAVE          SPGMAX, MEMRMX
C
C***********************************************************************
C---
C 1.  SI SOUSPG = 'FINFIN'  IMPRESSION DE MEMRMX UTILISE
C---
      IF ( SOUSPG.EQ.'FINFIN' ) THEN
        IF (IRANGP.GE.0) THEN
           CALL PARMEM(MEMRMX, SPGMAX)
        ENDIF
        WRITE(NFECRA,2000) MEMRMX,SPGMAX
        RETURN
      ENDIF
C
C---
C 3.  MISE A JOUR DE MEMRMX ET SPGMAX
C---
      IF ( MEMRDP.GT.MEMRMX ) THEN
        MEMRMX  = MEMRDP
        SPGMAX = SOUSPG
      ENDIF
C
C---
C 4.  TEST SUR LE DEPASSEMEMNT OU NON DEPASSEMENT DE MTOTRA
C---
      IF ( MEMRDP.LE.MTOTRA ) THEN
C
C--->      SI TAILLE SUFFISANTE: PAS D'IMPRESSION  ET RETURN
C
        RETURN
C
      ELSE
C
C--->      SI TAILLE INSUFFISANTE: MESSAGE ET STOP
C
        WRITE (NFECRA,1000) SOUSPG,MEMRDP,MTOTRA,MEMRDP-MTOTRA
        CALL CSEXIT (1)
C
      ENDIF
C
C---> FORMATS
C
 1000 FORMAT (//,
     & ' SOUS PROGRAMME APPELANT RASIZE            :    ',A6,/,
     & ' MEMOIRE NECESSAIRE RA (EN NBRE DE REELS)  = '   ,I12,/,
     & '         DISPONIBLE                        = '   ,I12,/,
     & ' ----> AUGMENTER LONGRA DANS LE LANCEUR DE = '   ,I12,' REELS'/)
 2000 FORMAT(I12,' REELS   NECESSAIRES : MAXIMUM ATTEINT DANS ', A6)
C
      END
c@z
