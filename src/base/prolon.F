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
                        SUBROUTINE PROLON
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NCELG  , IPROL  , NFECRA ,
     &   IRESPR , RF     , RG     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC PROLONGATION DE RG DANS RF SI IPROL = 1 ( ALORS RG,RF = VARIABLES)
CFONC RESTRICTION  DE RF DANS RG SI IPROL = 2 ( ALORS RG,RF = RESIDUS  )
CFONC
c@fonce
C-----------------------------------------------------------------------
C                       ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NCELG        ! E  !  -> ! NOMBRE DE MACRO CELLULES             !
CARGU ! IPROL        ! E  !  -> ! INDICATEUR = 1 PROLONGEMENT          !
CARGU !              !    !     !            = 2 RESTRICTION           !
CARGU ! IRESPR(NCELET! TE !  -> ! RESTR. DES RESIDUS ET PROL. DES SOL  !
CARGU ! RF   (NCELET ! TR ! <-> ! GRANDEUR SUR MAILLAGE FIN            !
CARGU ! RG   (NCELG) ! TR ! <-> ! GRANDEUR SUR MAILLAGE GROSSIER       !
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
C***********************************************************************
C
      INTEGER          NCELET , NCEL   , NCELG  , IPROL, NFECRA
C
      INTEGER          IRESPR(NCELET)
      DOUBLE PRECISION RF(NCELET), RG(NCELG)
C
C VARIABLES LOCALES
C
      INTEGER II, IEL
C
C***********************************************************************
C
C     PROLONGATION DE RG DANS RF
C
      IF( IPROL.EQ.1 ) THEN
        DO IEL=1,NCEL
          II = IRESPR(IEL)
          RF(IEL) = RG(II)
        ENDDO
C
C
C     RESTRICTION DES RF DANS RG
C
      ELSEIF( IPROL.EQ.2 ) THEN
        DO IEL=1,NCELG
          RG(IEL) = 0.D0
        ENDDO
        DO IEL=1,NCEL
          II = IRESPR(IEL)
          RG(II) = RG(II) +RF(IEL)
        ENDDO
C
C
      ELSE
        WRITE(NFECRA,1000)
        CALL CSEXIT (1)
C
      ENDIF
C
 1000 FORMAT(' PROLON APPELE AVEC IPROL =', I10)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
