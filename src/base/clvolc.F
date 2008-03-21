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
                        SUBROUTINE CLVOLC
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   ,
     &   VOLTOT , VOLUME )
C     -------------------------------------------------------------
C***********************************************************************
C
C  FONCTION  :
C  --------
c@foncb
CFONC
CFONC      CALCUL DU VOLUME GEOMETRIQUE DES ELEMENTS
CFONC      FORMULE DE GREEN : 3*VOLUME = SOMME DIV(R) AVEC R = t(X,Y,Z)
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! VOLTOT       ! R  ! <-  ! VOLUME TOTAL DU DOMAINE              !
CARGU ! VOLUME       ! TR ! <-  ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU .______________.____._____.______________________________________.
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
      INCLUDE "period.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER NCELET, NCEL
      DOUBLE PRECISION VOLTOT
      DOUBLE PRECISION VOLUME(NCELET)
C
      INTEGER IEL, IDIMTE, ITENSO
      DOUBLE PRECISION VOLMIN
C
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C=======================================================================
C 1. INITIALISATION EFFECTUEE DANS cs_maillage_grd.c
C=======================================================================
C
      IPASS = IPASS + 1
C
C=======================================================================
C 2. ON PREND LES MIN ET MAX
C=======================================================================
C
      VOLMIN = 1.D+12
      VOLTOT = 0.D0
C
      DO IEL  = 1, NCEL
C
         VOLMIN = MIN(VOLMIN, VOLUME(IEL))
         VOLTOT = VOLTOT + VOLUME(IEL)
C
      ENDDO
C
C       On communique le volume au halo pour les filtrages des modeles
C       dynamiques de L.E.S.
C       Dans le cas d'un voisinage etendu traite separement, l'appel a
C       PARCVE est fait directement dans la routine CFLITR (ce qui d'ailleurs
C       pourrait etre optimise)
      IF (IRANGP.GE.0) THEN
        CALL PARCOM(VOLUME)
C       ===========
        CALL PARMIN (VOLMIN)
C       ===========
        CALL PARSOM (VOLTOT)
C       ===========
      ENDIF
      IF(IPERIO.EQ.1) THEN
C
        IDIMTE = 0
        ITENSO = 0
C
        CALL PERCOM
C       ===========
     &       ( IDIMTE , ITENSO ,
     &         VOLUME , VOLUME , VOLUME ,
     &         VOLUME , VOLUME , VOLUME ,
     &         VOLUME , VOLUME , VOLUME )
C
      ENDIF
C
C     En ALE, on passe plusieurs fois ici.
C     Au premier passage (avant calculs) on ecrit, on teste et on s'arrete
C       si pb.
C     Aux passages suivants, on n'ecrit pas, on teste et on finit le pas
C       de temps si pb.
      IF (IPASS.EQ.1) THEN
        WRITE(NFECRA,1000) VOLMIN, VOLTOT
        IF (VOLMIN.LE.0.D0) THEN
          WRITE(NFECRA,1002)
          CALL CSEXIT (1)
        ENDIF
      ELSE
        IF (VOLMIN.LE.0.D0) THEN
          WRITE(NFECRA,1001) VOLMIN, VOLTOT
          WRITE(NFECRA,1002)
          NTMABS = NTCABS
        ENDIF
      ENDIF
C=======================================================================
C 5. FIN
C=======================================================================
C
 1000 FORMAT(
     &' --- Information sur les volumes                             ',/,
     &'       Volume de controle minimal = ',4X,E18.9                ,/,
     &'       Volume total du domaine    = ',4X,E18.9                 )
 1001 FORMAT(/,' CLVOLC : VOLUME DE CONTROLE MINIMAL     = ',E18.9,/,
     &         '          VOLUME TOTAL DU DOMAINE        = ',E18.9,/,/)
 1002 FORMAT(/,' CLVOLC : ARRET SUITE A LA DETECTION D''UN',/,
     &         '          VOLUME NEGATIF',/)
C
      RETURN
      END
c@z
