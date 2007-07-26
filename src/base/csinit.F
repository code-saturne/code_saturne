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
                        SUBROUTINE CSINIT
C                       *****************
C     -------------------------------------------------------------
     & ( ARGIFO , IRGPAR , NRGPAR , ILISR0, ILISRP )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  INIT DU LISTING ET DE PARAMETRES ASSOCIES AU PREPROCESSEUR
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! ARGIFO       ! E  !  -> ! VALEUR DE IFOENV ;                   !
CARGU !              !    !     ! FORMAT DE COMMUNICATION AVEC         !
CARGU !              !    !     ! LE PREPROCESSEUR                     !
CARGU !              !    !     !   0 : PAS DE COMMUNICATIONS          !
CARGU !              !    !     !   1 : COMMUNICATION PAR FICHIERS     !
CARGU ! IRGPAR       ! E  !  -> ! RANG SI PARALLELE ; -1 SI SEQUENTIEL !
CARGU ! NRGPAR       ! E  !  -> ! NOMBRE DE PROCESSUS ; 1 SI SEQUENTIEL!
CARGU ! ILISR0       ! E  !  -> ! OPTION DE SORTIE DU LISTING :        !
CARGU !              !    !     !   0 : RANG 0 NON REDIRIGE            !
CARGU !              !    !     !   1 : RANG 0 DANS FICHIER listing,   !
CARGU ! ILISRP       ! E  !  -> ! OPTION DE SORTIE DU LISTING :        !
CARGU !              !    !     !   0 : RANGS > 0 NON REDIRIGES        !
CARGU !              !    !     !       (POUR DEBUGGER PAR EXEMPLE)    !
CARGU !              !    !     !   1 : RANGS > 0 REDIRIGES DANS       !
CARGU !              !    !     !       FICHIERS listing_n*            !
CARGU !              !    !     !   2 : RANGS > 0 REDIRIGES DANS       !
CARGU !              !    !     !       /dev/null (SUPPRESSION)        !
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
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
      INTEGER          ARGIFO, IRGPAR, NRGPAR, ILISR0, ILISRP
C
      CHARACTER        NAME*300
C
C***********************************************************************
C
C=======================================================================
C Initialisation du common IPARAL
C=======================================================================
C
      IRANGP = IRGPAR
      NRANGP = NRGPAR
C
C=======================================================================
C Initialisation des paramètres de lecture des données Préprocesseur
C=======================================================================
C
      IFOENV = ARGIFO
C
C-----------------------------------
C --- stdout :
C     NFECRA = 6 par défaut
C   Les autres fichiers listing sont fermes dans csclli a la fin de cs_main.
C-----------------------------------
C
      NFECRA = 6
C
      IF (IRANGP.LE.0) THEN
        IF (ILISR0.EQ.1) THEN
          NFECRA = 9
          NAME = 'listing'
        ENDIF
      ELSE
        IF (ILISRP.EQ.1) THEN
          NFECRA = 9
          WRITE (NAME,'(A9,I4.4)') 'listing_n', IRANGP + 1
        ELSE IF (ILISRP.EQ.2) THEN
          NFECRA = 9
          NAME = '/dev/null'
        ENDIF
      ENDIF
C
      IF (NFECRA.EQ.9) THEN
         OPEN (FILE=NAME, UNIT=NFECRA,
     &         FORM='FORMATTED', STATUS='UNKNOWN', ERR=900)
      ENDIF
C
      GOTO 950
C
 900  WRITE (0, 999) NAME
      CALL CSEXIT (1)
C
 950  CONTINUE
C
 999  FORMAT(/,
     &'Code_Saturne : Erreur d''initialisation :',/,
     &'Impossible d''ouvrir le fichier : ',A,/)
C
      RETURN
      END
c@z
