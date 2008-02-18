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
                        SUBROUTINE CPRAYS
C                       *****************
C     -------------------------------------------------------------
     &  ( IVAR   , NCELET , NCEL   ,
     &    VOLUME , PROPCE , RAYE   , RAYI   , SMBRS  , ROVSDT )
C     -------------------------------------------------------------
C
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC ROUTINE PHYSIQUE PARTICULIERE : FLAMME CHARBON PULVERISE
CFONC    PRISE EN COMPTE DES TERMES SOURCES RADIATIFS
CFONC    IMPLICITE ET EXPLICITE DANS L'EQUATION DES PARTICULES
CFONC    DE LA CLASSE ICLA
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IVAR         ! E  !  -> ! NUMERO DE LA VARIABLE SCALAIRE       !
CARGU !              !    !     !   ENERGIE (ENTHALPIE H2) POUR LE     !
CARGU !              !    !     !   CHARBON                            !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! VOLUME(NCELET! TR !  -> ! VOLUME DES CELLULES                  !
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! RAYE(NCELET  ! TR !  -> ! TERME SOURCE RADIATIF EXPLICITE      !
CARGU ! RAYI(NCELET  ! TR ! <-> ! TERME SOURCE RADIATIF IMPLICITE      !
CARGU ! SMBRS(NCELET ! TR ! <-> ! SECOND MEMBRE DU SYSTEME             !
CARGU ! ROVSDT(NCELET! TR ! <-> ! DIAGONALE DU SYSTEME                 !
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
      INCLUDE "numvar.h"
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
      INTEGER          IVAR , NCELET, NCEL
C
      DOUBLE PRECISION VOLUME(NCELET)
      DOUBLE PRECISION SMBRS(NCELET)
      DOUBLE PRECISION ROVSDT(NCELET)
      DOUBLE PRECISION RAYE(NCELET,NPHASC)
      DOUBLE PRECISION RAYI(NCELET,NPHASC)
      DOUBLE PRECISION PROPCE(NCELET,*)
C
C VARIABLES LOCALES
C
      INTEGER          IEL , NUMCLA , IPCL
C
C***********************************************************************
C
C=======================================================================
C 1. RECHERCHE DE LA ZONE MEMOIRE (IPH) EN FONCTION DU NUMERO DE PHASE
C    COURANT IPHAS POUR TROUVER LES BONS TERMES SOURCES
C=======================================================================
C
      NUMCLA = IVAR-ISCA(IH2(1))+1
      IPCL   = 1+NUMCLA
C
C=======================================================================
C 2. PRISE EN COMPTE DES TERMES SOURCES RADIATIFS
C=======================================================================
C
C
      DO IEL = 1,NCEL
        RAYI(IEL,IPCL) = MAX(-RAYI(IEL,IPCL),ZERO)
      ENDDO
C
      DO IEL = 1,NCEL
        IF ( PROPCE(IEL,IPPROC(IX2(NUMCLA))) .GT. EPZERO ) THEN
C
C--> PARTIE EXPLICITE
C
          SMBRS(IEL)  = SMBRS(IEL) +  RAYE(IEL,IPCL)*VOLUME(IEL)
     &                               *PROPCE(IEL,IPPROC(IX2(NUMCLA)))
C
C--> PARTIE IMPLICITE
C
          ROVSDT(IEL) = ROVSDT(IEL) + RAYI(IEL,IPCL)*VOLUME(IEL)
        ENDIF
C
      ENDDO
C
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
