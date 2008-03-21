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
                        SUBROUTINE CATSMA
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NCESMP , ITERNS , ISNEXP ,
     &   THETV  ,
     &   ICETSM , ITPSMP ,
     &   VOLUME , PVARA  , SMCELP , GAMMA  ,
     &   TSEXP  , TSIMP  , GAPINJ )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  CALCUL DES TERMES SOURCES IMPLICITE ET EXPLICITE
CFONC   VENANT DES SOURCES DE MASSE
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE DE CELLULES                   !
CARGU ! NCESMP       ! E  !  -> ! NOMBRE DE CELLULES A SOURCE DE MASSE !
CARGU ! ITERNS       ! E  !  -> ! NUMERO D'ITERATION SUR NAVSTO        !
CARGU ! ISNEXP       ! E  !  -> ! INDICATEUR POUR EXTRAPOLATION DES    !
CARGU !              !    !     !  TERMES SOURCES DE LA PHASE TRAITEE  !
CARGU ! THETV        ! R  !  -> ! THETA SCH. POUR LA VARIABLE          !
CARGU !              !    !     !    THETV  V(N+1) + (1-THETV) V(N)    !
CARGU ! ICETSM(NCESMP! TE !  -> ! NUMERO DES CELLULES A SOURCE DE MASSE!
CARGU ! ITPSMP(NCESMP! TE !  -> ! TYPE DE SOURCE DE MASSE POUR LA      !
CARGU !              !    !     !  VARIABLE TRAITEE (cf. USTSMA)       !
CARGU ! VOLUME(NCEL) ! TR !  -> ! VOLUME DES CELLULES                  !
CARGU ! PVARA(NCEL)  ! TR !  -> ! VALEUR DE LA VARIABLE EN DEBUT DE PAS!
CARGU !              !    !     !  DE TEMPS                            !
CARGU ! SMCELP(NCESMP! TR !  -> ! VALEUR DE LA VARIABLE ASSOCIEE A LA  !
CARGU !              !    !     !  SOURCE DE MASSE                     !
CARGU ! GAMMA(NCESMP)! TR !  -> ! VALEUR DU FLUX DE MASSE              !
CARGU ! TSEXP(NCEL)  ! TR ! <-> ! PART DU TERME SOURCE EXPLICITE       !
CARGU !              !    !     !  LINEAIRE PAR RAPPORT A LA VARIABLE  !
CARGU ! TSIMP(NCEL)  ! TR ! <-> ! VALEUR ASSOCIEE A TSEXP QUI IRA DANS !
CARGU !              !    !     !  LA MATRICE                          !
CARGU ! GAPINJ(NCEL) ! TR ! <-  ! PART DU TERME SOURCE EXPLICITE       !
CARGU !              !    !     !  INDEPENDANT DE LA VARIABLE          !
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
      INTEGER          NCELET, NCEL  , NCESMP, ITERNS, ISNEXP
      INTEGER          ICETSM(NCESMP), ITPSMP(NCESMP)
      DOUBLE PRECISION THETV
      DOUBLE PRECISION VOLUME(NCELET)
      DOUBLE PRECISION PVARA (NCELET)
      DOUBLE PRECISION SMCELP(NCESMP), GAMMA (NCESMP)
      DOUBLE PRECISION TSEXP (NCELET), TSIMP (NCELET), GAPINJ(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER II, IEL
C
C***********************************************************************
C
C Explication des tests GAMMA(II).GT.0.D0 .AND. ITPSMP(II).EQ.1 :
C     Si on enleve de la matiere ou si on entre a la valeur de la
C       cellule alors l'equation de IVAR n'est pas modifiee
C     Sinon, on ajoute le terme source GAMMA*(f_i-f^(n+1))
C
C     Dans TSIMP, on ajoute le terme qui ira sur la diagonale,
C       soit Gamma (*theta en cas d'ordre 2)
C     Dans TSEXP on ajoute le terme correspondant du second membre
C       cad Gamma * Pvar (avec Pvar)
C     Dans GAPINJ on place le terme Gamma Pinj qui ira au second membre
C
C     La distinction entre TSEXP et W1 (qui vont finalement tous les
C       deux au second membre) sert pour l'ordre 2 en temps.
C
      IF(ITERNS.EQ.1) THEN
        DO IEL = 1, NCEL
          GAPINJ(IEL) = 0.D0
        ENDDO
        DO II = 1, NCESMP
          IEL = ICETSM(II)
          IF (GAMMA(II).GT.0.D0 .AND. ITPSMP(II).EQ.1) THEN
            TSEXP(IEL) = TSEXP(IEL) - VOLUME(IEL)*GAMMA(II) * PVARA(IEL)
            GAPINJ(IEL) = VOLUME(IEL)*GAMMA(II) * SMCELP(II)
          ENDIF
        ENDDO
      ENDIF
C
C     Sur la diagonale
      IF(ISNEXP.GT.0) THEN
        DO II = 1, NCESMP
          IEL = ICETSM(II)
          IF (GAMMA(II).GT.0.D0 .AND. ITPSMP(II).EQ.1) THEN
            TSIMP(IEL) = TSIMP(IEL) + VOLUME(IEL)*GAMMA(II)*THETV
          ENDIF
        ENDDO
      ELSE
        DO II = 1, NCESMP
          IEL = ICETSM(II)
          IF (GAMMA(II).GT.0.D0 .AND. ITPSMP(II).EQ.1) THEN
            TSIMP(IEL) = TSIMP(IEL) + VOLUME(IEL)*GAMMA(II)
          ENDIF
        ENDDO
      ENDIF
C
      RETURN
      END
c@z
