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
                        SUBROUTINE RAYSCA
C                       *****************
C     -------------------------------------------------------------
     & ( IISCA  , IPHAS  ,
     &   NCELET , NCEL   ,
     &   SMBRS  , ROVSDT , VOLUME ,
     &   RAYE   , RAYI   )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
CFONC    -----------------------------------------
CFONC
CFONC        PRISE EN COMPTE DES TERMES SOURCES RADIATIFS
CFONC        IMPLICITE ET EXPLICITE
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IISCA        ! E  !  -> ! NUM SCALAIRE TEMPERATURE OU ENTHALPIE!
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE LA PHASE COURANTE          !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! SMBRS(NCELET)! TR ! <-> ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! ROVSDT(NCELET! TR ! <-> ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! VOLUME(NCELET! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! RAYE(NCELET) ! TR !  -> ! TERME SOURCE RAYONNEMENT EXPLICITE   !
CARGU ! RAYI(NCELET) ! TR !  -> ! TERME SOURCE RAYONNEMENT IMPLICITE   !
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
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "radiat.h"
      INCLUDE "entsor.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IISCA , NCELET , NCEL, IPHAS
C
      DOUBLE PRECISION VOLUME(NCELET)
      DOUBLE PRECISION SMBRS(NCELET)
      DOUBLE PRECISION ROVSDT(NCELET)
      DOUBLE PRECISION RAYE(NCELET,NPHAST)
      DOUBLE PRECISION RAYI(NCELET,NPHAST)
C
C VARIABLES LOCALES
C
      INTEGER          IEL , IPH, IIPH
C
C***********************************************************************
C=======================================================================
C 1. RECHERCHE DE LA ZONE MEMOIRE (IPH) EN FONCTION DU NUMERO DE PHASE
C    COURANT IPHAS POUR TROUVER LES BONS TERMES SOURCES
C=======================================================================
C
C     On determine le numero de phase reduit pour le rayonnement
      IPH = 0
      DO IIPH = 1, NPHAST
        IF (IRAPHA(IIPH).EQ.IPHAS) THEN
          IPH = IIPH
        ENDIF
      ENDDO
      IF(IPH.EQ.0) THEN
        WRITE(NFECRA,9001)
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C=======================================================================
C 2. PRISE EN COMPTE DES TERMES SOURCES RADIATIFS
C=======================================================================
C
C
      IF (ABS(ISCSTH(IISCA)).EQ.1 .OR. ISCSTH(IISCA).EQ.2) THEN
C
        DO IEL = 1,NCEL
          RAYI(IEL,IPH) = MAX(-RAYI(IEL,IPH),ZERO)
        ENDDO
C
        DO IEL = 1,NCEL
C
C--> PARTIE EXPLICITE
C
          SMBRS(IEL)  = SMBRS(IEL) +  RAYE(IEL,IPH)*VOLUME(IEL)
C
C--> PARTIE IMPLICITE
C
          ROVSDT(IEL) = ROVSDT(IEL) + RAYI(IEL,IPH)*VOLUME(IEL)
C
        ENDDO
C
      ENDIF
C
C--------
C FORMATS
C--------
C
 9001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT                                 ',/,
     &'@    *********                                               ',/,
     &'@      ERREUR SUR LE NOMBRE DE PHASES DANS RAYSCA            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
