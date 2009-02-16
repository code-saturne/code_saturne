c@a
c@versb
C-----------------------------------------------------------------------
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2009 EDF S.A., France
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
                        SUBROUTINE TSTVEC
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NFAC   , NFABOR ,
     &   IFACEL , IFABOR ,
     &   IWORKF , ISMBS  , ISMBV  ,
     &   RWORKF , RSMBS  , RSMBV )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ---------
c@foncb
CFONC
CFONC VERIFICATION DE LA VECTORISATION APRES RENUMEROTATION
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                 ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC  /NFABOR! E  !  -> ! NOMBRE TOTAL DE FACES INTERNES/DE BRD!
CARGU ! IFACEL       ! TE ! <-> ! No DES ELTS VOISINS D'UNE FACE INTERN!
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE ! <-> ! No DE L'ELT VOISIN D'UNE FACE DE BORD!
CARGU ! NFABOR  )    !    !     !                                      !
CARGU ! IWORKF(*     ! TE !  -  ! TAB DE TRAV DE DIM MAX(NFAC,NFABOR)  !
CARGU ! ISMBS (NCELET! TE !  -  ! TAB DE TRAV POUR ASSEMBLAGE SCALAIRE !
CARGU ! ISMBV (NCELET! TE !  -  ! TAB DE TRAV POUR ASSEMBLAGE VECTORIEL!
CARGU ! IPNFAW       ! TE !  -  ! TAB DE TRAV POUR IPNFAC              !
CARGU !   (NFAC+1)   !    !     !                                      !
CARGU ! NODFAW       ! TE !  -  ! TAB DE TRAV POUR NODFAC              !
CARGU !   (LNDFAC)   !    !     !                                      !
CARGU ! IPNFBW       ! TE !  -  ! TAB DE TRAV POUR IPNFBR              !
CARGU !   (NFABOR+1) !    !     !                                      !
CARGU ! NODFBW       ! TE !  -  ! TAB DE TRAV POUR NODFBR              !
CARGU !   (LNDFBR)   !    !     !                                      !
CARGU ! RWORKF(*     ! TR !  -  ! TAB DE TRAV DE DIM MAX(NFAC,NFABOR)  !
CARGU ! RSMBS (NCELET! TR !  -  ! TAB DE TRAV POUR ASSEMBLAGE SCALAIRE !
CARGU ! RSMBV (NCELET! TR !  -  ! TAB DE TRAV POUR ASSEMBLAGE VECTORIEL!
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
      INCLUDE "vector.h"
      INCLUDE "entsor.h"
C
C ARGUMENTS
C
      INTEGER          NCELET, NCEL, NFAC, NFABOR
      INTEGER          IRVECI, IRVECB
      INTEGER          IFACEL(2,NFAC),IFABOR(NFABOR)
      INTEGER          IWORKF(*), ISMBS(NCELET), ISMBV(NCELET)
      DOUBLE PRECISION RWORKF(*), RSMBS(NCELET), RSMBV(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          II, JJ, IOK, IFAC
      INTEGER          IEL, ISTOP
C
C
C***********************************************************************
C
C -----> Test d'assemblage sur des entiers
C
      ISTOP = 0
C
C --- Faces internes
C
      IF(IVECTI.EQ.1) THEN
C
        DO IFAC = 1, NFAC
          IWORKF(IFAC) = 1
        ENDDO
        DO IEL = 1, NCELET
          ISMBV(IEL) = 0
          ISMBS(IEL) = 0
        ENDDO
C
!OCL NOVREC,VRL(16)
        DO IFAC = 1, NFAC
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
          ISMBV(II) = ISMBV(II) + IWORKF(IFAC)
          ISMBV(JJ) = ISMBV(JJ) + IWORKF(IFAC)
        ENDDO
C
C VECTORISATION NON FORCEE
        DO IFAC = 1, NFAC
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
          ISMBS(II) = ISMBS(II) + IWORKF(IFAC)
          ISMBS(JJ) = ISMBS(JJ) + IWORKF(IFAC)
        ENDDO
C
        IOK = 0
        DO IEL = 1, NCEL
          IF(ISMBS(IEL).NE.ISMBV(IEL)) THEN
            IOK = IOK - 1
          ENDIF
        ENDDO
C
        IF(IOK.NE.0) THEN
          WRITE(NFECRA,3101)IOK
          ISTOP = 1
        ENDIF
C
      ENDIF
C
C --- Faces de bord
C
      IF(IVECTB.EQ.1) THEN
C
        DO IFAC = 1, NFABOR
          IWORKF(IFAC) = 1
        ENDDO
        DO IEL = 1, NCEL
          ISMBV(IEL) = 0
          ISMBS(IEL) = 0
        ENDDO
C
!OCL NOVREC,VRL(16)
        DO IFAC = 1, NFABOR
          II = IFABOR(IFAC)
          ISMBV(II) = ISMBV(II) + IWORKF(IFAC)
        ENDDO
C
C VECTORISATION NON FORCEE
        DO IFAC = 1, NFABOR
          II = IFABOR(IFAC)
          ISMBS(II) = ISMBS(II) + IWORKF(IFAC)
        ENDDO
C
        IOK = 0
        DO IEL = 1, NCEL
          IF(ISMBS(IEL).NE.ISMBV(IEL)) THEN
            IOK = IOK - 1
          ENDIF
        ENDDO
C
        IF(IOK.NE.0) THEN
          WRITE(NFECRA,3201)IOK
          ISTOP = 1
        ENDIF
C
      ENDIF
C
C
C -----> Test d'assemblage sur des reels
C
C --- Faces internes
C
      IF(IVECTI.EQ.1) THEN
C
        DO IFAC = 1, NFAC
          RWORKF(IFAC) = 1.D0
        ENDDO
        DO IEL = 1, NCELET
          RSMBV(IEL) = 0.D0
          RSMBS(IEL) = 0.D0
        ENDDO
C
!OCL NOVREC,VRL(16)
        DO IFAC = 1, NFAC
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
          RSMBV(II) = RSMBV(II) + RWORKF(IFAC)
          RSMBV(JJ) = RSMBV(JJ) + RWORKF(IFAC)
        ENDDO
C
C VECTORISATION NON FORCEE
        DO IFAC = 1, NFAC
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
          RSMBS(II) = RSMBS(II) + RWORKF(IFAC)
          RSMBS(JJ) = RSMBS(JJ) + RWORKF(IFAC)
        ENDDO
C
        IOK = 0
        DO IEL = 1, NCEL
          IF(ABS(RSMBS(IEL)-RSMBV(IEL)).GT.1.D-20) THEN
            IOK = IOK - 1
          ENDIF
        ENDDO
C
        IF(IOK.NE.0) THEN
          WRITE(NFECRA,3102)IOK
          ISTOP  = 1
        ENDIF
C
      ENDIF
C
C --- Faces de bord
C
      IF(IVECTB.EQ.1) THEN
C
        DO IFAC = 1, NFABOR
          RWORKF(IFAC) = 1.D0
        ENDDO
        DO IEL = 1, NCEL
          RSMBV(IEL) = 0.D0
          RSMBS(IEL) = 0.D0
        ENDDO
C
!OCL NOVREC,VRL(16)
        DO IFAC = 1, NFABOR
          II = IFABOR(IFAC)
          RSMBV(II) = RSMBV(II) + RWORKF(IFAC)
        ENDDO
C
C VECTORISATION NON FORCEE
        DO IFAC = 1, NFABOR
          II = IFABOR(IFAC)
          RSMBS(II) = RSMBS(II) + RWORKF(IFAC)
        ENDDO
C
        IOK = 0
        DO IEL = 1, NCEL
          IF(ABS(RSMBS(IEL)-RSMBV(IEL)).GT.1.D-20) THEN
            IOK = IOK - 1
          ENDIF
        ENDDO
C
        IF(IOK.NE.0) THEN
          WRITE(NFECRA,3202)IOK
          ISTOP  = 1
        ENDIF
C
      ENDIF
C
      IF(ISTOP.NE.0) THEN
        WRITE(NFECRA,9000)
        CALL CSEXIT (1)
      ENDIF
C
C=======================================================================
C 6. FORMATS
C=======================================================================
C
#if defined(_CS_LANG_FR)
C
 3101 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS tstvec                           ',/,
     &'@    *********                                               ',/,
     &'@    PROBLEME D''ASSEMBLAGE D''ENTIERS AUX FACES INTERNES    ',/,
     &'@    ',I10   ,' OCCURRENCES DU PROBLEME                      ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Verifier le maillage.                                     ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3201 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS tstvec                           ',/,
     &'@    *********                                               ',/,
     &'@    PROBLEME D''ASSEMBLAGE D''ENTIERS AUX FACES DE BORD     ',/,
     &'@    ',I10   ,' OCCURRENCES DU PROBLEME                      ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Verifier le maillage.                                     ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3102 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS tstvec                           ',/,
     &'@    *********                                               ',/,
     &'@    PROBLEME D''ASSEMBLAGE DE REELS   AUX FACES INTERNES    ',/,
     &'@    ',I10   ,' OCCURRENCES DU PROBLEME                      ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Verifier le maillage.                                     ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3202 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS tstvec                           ',/,
     &'@    *********                                               ',/,
     &'@    PROBLEME D''ASSEMBLAGE DE REELS   AUX FACES DE BORD     ',/,
     &'@    ',I10   ,' OCCURRENCES DU PROBLEME                      ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Verifier le maillage.                                     ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS tstvec                           ',/,
     &'@    *********                                               ',/,
     &'@     CONFIGURATION IMPREVUE A LA RENUMEROTATION DES FACES   ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Verifier le maillage.                                     ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#else
C
 3101 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN tstvec                                ',/,
     &'@    ********                                                ',/,
     &'@    PROBLEM IN ASSEMBLING INTEGERS AT THE INTERIOR FACES    ',/,
     &'@    ',I10   ,' OCCURRENCES OF THE PROBLEM                   ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify the mesh.                                          ',/,
     &'@  Contact the support.                                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3201 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN tstvec                                ',/,
     &'@    ********                                                ',/,
     &'@    PROBLEM IN ASSEMBLING INTEGERS AT THE BOUNDARY FACES    ',/,
     &'@    ',I10   ,' OCCURRENCES OF THE PROBLEM                   ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify the mesh.                                          ',/,
     &'@  Contact the support.                                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3102 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN tstvec                                ',/,
     &'@    ********                                                ',/,
     &'@    PROBLEM IN ASSEMBLING REALS AT THE INTERIOR FACES       ',/,
     &'@    ',I10   ,' OCCURRENCES OF THE PROBLEM                   ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify the mesh.                                          ',/,
     &'@  Contact the support.                                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3202 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN tstvec                                ',/,
     &'@    ********                                                ',/,
     &'@    PROBLEM IN ASSEMBLING REALS AT THE BOUNDARY FACES       ',/,
     &'@    ',I10   ,' OCCURRENCES OF THE PROBLEM                   ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify the mesh.                                          ',/,
     &'@  Contact the support.                                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ WARNING: ABORT IN tstvec                                ',/,
     &'@    ********                                                ',/,
     &'@     UNEXPECTED CONFIGURATION DURING THE FACES RENUMBERING  ',/,
     &'@                                                            ',/,
     &'@  The calculation will not be run.                          ',/,
     &'@                                                            ',/,
     &'@  Verify the mesh.                                          ',/,
     &'@  Contact the support.                                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
#endif
C
      RETURN
      END
c@z
