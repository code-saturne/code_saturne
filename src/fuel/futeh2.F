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
                        SUBROUTINE FUTEH2
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NRTUSE ,
     &   RTP    , PROPCE , RTUSER )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC CALCUL DE LA TEMPERATURE DES PARTICULES
CFONC   EN FONCTION DE L'ENTHALPIE DU FOL ET DES CONCENTRATIONS
CFONC
c@fonce
C                             ARGUMENTS
c@argub
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT)        !
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! EH0          ! TR ! <-> ! TABLEAU REEL DE TRAVAIL              !
CARGU ! EH1          ! TR ! <-> ! TABLEAU REEL DE TRAVAIL              !
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
C     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
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
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "fuincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET, NCEL , NRTUSE
      DOUBLE PRECISION RTP(NCELET,*), PROPCE(NCELET,*)
      DOUBLE PRECISION RTUSER(NRTUSE)
C
C VARIABLES LOCALES
C
      INTEGER          ICEL
      INTEGER          IPCTE1 , IPCTE2
      INTEGER          MODE
      DOUBLE PRECISION EH2, XSOLID(2)
      DOUBLE PRECISION MKFINI,DIAMGT
      DOUBLE PRECISION MASGUT,MFGOUT,MKGOUT,RHOFOL
C
C***********************************************************************
C RQ IMPORTANTE : On suppose pour l'instant que H2 = H02 + CP2(T2-TREF)
C
C
C
C=======================================================================
C 1. CALCULS PRELIMINAIRES
C=======================================================================
C
C --- Initialisation des tableaux
C
      IPCTE1 = IPPROC(ITEMP1)
      IPCTE2 = IPPROC(ITEMP3)
C
      DO ICEL = 1, NCEL
C
C --- Initialisation de T2 a T20
C
        PROPCE(ICEL,IPCTE2) = 373.D0
C 20/09/05 Pour l'instant TinFol en dur, un jour on recuperera dans USFUCL
      ENDDO
C
C=======================================================================
C 2. CALCUL DE LA TEMPERATURE DES PARTICULES
C=======================================================================
C
      MODE = 1
C
      MKFINI = RHO0FL*PI/6.D0*DINIKF**3
C
      DO ICEL = 1, NCEL
C
        RHOFOL = PROPCE(ICEL,IPPROC(IROM3))
        DIAMGT = PROPCE(ICEL,IPPROC(IDIAM3))
        MASGUT = RHO0FL*PI/6.D0*DIAMGT**3
        IF (DIAMGT.LE.DINIKF) THEN
          MKGOUT = MASGUT
        ELSE
          MKGOUT = MKFINI
        ENDIF
        MFGOUT = MASGUT - MKGOUT
        XSOLID(1) = 1.D0-FKC
        XSOLID(2) = FKC
        IF(MASGUT.GT.ZERO) THEN
          XSOLID(1) = MFGOUT / MASGUT
          XSOLID(2) = MKGOUT / MASGUT
        ENDIF
        XSOLID(1) = MIN(1.D0,MAX(0.D0,XSOLID(1)))
        XSOLID(2) = Min(1.D0,Max(0.D0,XSOLID(2)))
C
        IF ( RTP(ICEL,ISCA(IYFOL)) .GT. (3.D3*EPSIFL) ) THEN
          EH2 =  RTP(ICEL,ISCA(IHLF))/RTP(ICEL,ISCA(IYFOL))
          CALL FUTHP2
C         ============
     &    ( MODE , EH2 , XSOLID , PROPCE(ICEL,IPCTE2) )
        ELSE
          PROPCE(ICEL,IPCTE2) = PROPCE(ICEL,IPCTE1)
        ENDIF
C
      ENDDO
C
C----
C FIN
C----
C
      RETURN
      END
c@z
