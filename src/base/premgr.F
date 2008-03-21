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
                        SUBROUTINE PREMGR
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NCEL   , NFAC   , IFACEL ,
     &   IP     , INDIC  , IW     , ICELFA  )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC        MULTIGRILLE ALGEBRIQUE :
CFONC   VERIFICATIONRS PRELIMINAIRES
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS MAILLAGE FIN       !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES MAILL. FIN  !
CARGU ! IFACEL       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFAC)    !    !     ! SUR MAIL FIN                         !
CARGU ! IP(NCEL)     ! TE !  -  ! POINTEURS SUR VOISINS POUR CONNECTI- !
CARGU !              !    !     ! VITE INVERSE                         !
CARGU ! INDIC(NFAC)  ! TE !  -  ! TABLEAU DE TRAVAIL
CARGU ! IW(NCEL)     ! TE !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! ICELFA       ! TE !  -  ! CONNCETIVITE CELLULES->FACES MAILLA- !
CARGU !  (2*NFAC)    !    !     ! GE FIN                               !
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
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NCEL, NFAC
      INTEGER          IFACEL(2,NFAC), IP(NCEL), INDIC(NFAC), IW(NCEL)
      INTEGER          ICELFA(2*NFAC)
C
C VARIABLES LOCALES
C
      INTEGER          IEL , IFAC, II, JJ, IV, JV, NI, NJ
      INTEGER          IDEBI, IFINI, IDEBJ, IFINJ
      INTEGER          IDECOU
C
C***********************************************************************
C
C=======================================================================
C 1.  CONSTRUCTION DE LA CONNECTIVITE        INVERSE : IP , ICELFA
C=======================================================================
C
C INITIALISATION
C
      DO IEL = 1, NCEL
        IW(IEL) = 0
      ENDDO
C
C CREATION DU TABLEAU IP DE POINTEURS SUR LES VOISINS DE CELLULES
C IP(1)                   : NOMBRE DE VOISINS DE LA CELLULE 1
C IP(K) - IP(K-1), K >= 2 : NOMBRE DE VOISINS DE LA CELLULE K
C IP(NCEL)                : = 2*NFAC
      DO IFAC = 1, NFAC
        II = IFACEL(1,IFAC)
        JJ = IFACEL(2,IFAC)
        IW(II) = IW(II) + 1
        IW(JJ) = IW(JJ) + 1
      ENDDO
      DO IEL = 2, NCEL
        IW(IEL) = IW(IEL) + IW(IEL-1)
      ENDDO
C RQ : CI-DESSOUS, IP EST DECREMENTE DU NOMBRE DE VOISINS
C      CAR IL SERA REINCREMENTE DANS LA BOUCLE SUIVANTE
      IP(1) = 0
      DO IEL = 2, NCEL
        IP(IEL) = IW(IEL-1)
      ENDDO
C
      DO IFAC = 1, NFAC
        II = IFACEL(1,IFAC)
        JJ = IFACEL(2,IFAC)
        IP(II) = IP(II) + 1
        IP(JJ) = IP(JJ) + 1
        ICELFA(IP(II)) = IFAC
        ICELFA(IP(JJ)) = IFAC
      ENDDO
C
C=======================================================================
C 2.  DETECTION DES FACES REDECOUPEES
C     Pour l'instant, on ne fait que verifier qu'il n'y a pas
C     de faces redecoupees sinon on stoppe
C=======================================================================
C
      IDECOU = 0
      DO IFAC = 1, NFAC
        INDIC(IFAC) = 0
      ENDDO
C
      DO IFAC = 1, NFAC
        II = IFACEL(1,IFAC)
        JJ = IFACEL(2,IFAC)
C
        IF (INDIC(IFAC).EQ.0) THEN
          IF (II.EQ.1) THEN
            IDEBI = 1
            IFINI = IP(1)
          ELSE
            IDEBI = IP(II-1) + 1
            IFINI = IP(II)
          ENDIF
          IF (JJ.EQ.1) THEN
            IDEBJ = 1
            IFINJ = IP(1)
          ELSE
            IDEBJ = IP(JJ-1) + 1
            IFINJ = IP(JJ)
          ENDIF
C
          DO NI = IDEBI, IFINI
            DO NJ = IDEBJ, IFINJ
              IV = ICELFA(NI)
              JV = ICELFA(NJ)
              IF (IV.EQ.JV.AND.IV.NE.IFAC) THEN
                INDIC(IFAC) = 1
                INDIC(IV)   = 1
                IDECOU = IDECOU + 1
                GOTO 111
              ENDIF
            ENDDO
          ENDDO
 111      CONTINUE
        ENDIF
      ENDDO
C
      IF (IDECOU.NE.0) THEN
        WRITE(NFECRA,1000) IDECOU
        CALL CSEXIT (1)
      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS premgr (multigrille)             ',/,
     &'@    *********                                               ',/,
     &'@    DETECTION DE ',I10   ,' FACES REDECOUPEES               ',/,
     &'@                                                            ',/,
     &'@  Le maillage utilise contient des faces redecoupees, i.e.  ',/,
     &'@    deux faces distinctes peuvent avoir le meme couple de   ',/,
     &'@    cellules voisines.                                      ',/,
     &'@  Cette configuration n''est pas prise en compte dans la    ',/,
     &'@    version actuelle lorsque le multigrille est active.     ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Desactiver l''option multigrille (IMGR = 1) dans usini1.  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      END
c@z
